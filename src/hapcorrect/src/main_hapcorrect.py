#!/usr/bin/env python3
import statistics
import sys
import shutil

import numpy
import pysam
import argparse
import logging
import os
import pandas as pd

from multiprocessing import Pool
from collections import defaultdict

logger = logging.getLogger()

from src.hapcorrect.src.process_bam import (get_all_reads_parallel, update_coverage_hist, get_segments_coverage, haplotype_update_all_bins_parallel,
                                            get_snps_frequencies, tumor_bam_haplotag)
from src.hapcorrect.src.process_vcf import (vcf_parse_to_csv_for_het_phased_snps_phasesets, get_snp_frequencies_segments, snps_frequencies_chrom_mean,
                                            get_snps_frquncies_coverage, vcf_parse_to_csv_for_snps, index_vcf, rephase_vcf, get_phasingblocks,
                                            snps_frequencies_chrom_mean_phasesets, snps_frequencies_homo_chrom_mean)
from src.hapcorrect.src.phase_correction import (generate_phasesets_bins, phase_flips_cis_trans, switch_inter_phaseblocks_bins,
                                                 phaseblock_flipping_simple_heuristics, check_missing_phasesets_original, reintroduce_broken_phasesets,
                                                 update_remaining_phasesets, subtract_intervals, without_phasesets_bins_correction,
                                                 remove_centromere_phaseblocks)
from src.hapcorrect.src.utils import (get_chromosomes_bins, write_segments_coverage, csv_df_chromosomes_sorter, get_snps_frquncies_coverage_from_bam,
                                      adjust_loh_cent_phaseblocks, extract_centromere_regions, infer_missing_phaseblocks, df_chromosomes_sorter,
                                      is_phasesets_check_simple_heuristics, write_df_csv, loh_regions_events, snps_frequencies_chrom_genes,
                                      genes_segments_coverage, genes_segments_list, extend_snps_ratios_df, get_chromosomes_regions,
                                      add_breakpoints, update_hp_assignment_loh_segments)
from src.hapcorrect.src.extras import get_contigs_list
from src.hapcorrect.src.plots import plot_coverage_data, change_point_detection, plot_coverage_data_after_correction, loh_plots_genome
from src.hapcorrect.src.cpd import cpd_positions_means
from src.hapcorrect.src.loh import detect_loh_centromere_regions, plot_snps

MIN_SV_SIZE = 50


def safe_rmtree(path):
        if os.path.islink(path):
            os.unlink(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)
        elif os.path.exists(path):
            os.remove(path)


def main_process(args):
    centromere_regions = extract_centromere_regions(args)
    if args.control_bam is None:
        args.control_bam = []
    all_bams = args.target_bam + args.control_bam
    target_genomes = set(os.path.basename(b) for b in args.target_bam)
    control_genomes = set(os.path.basename(b) for b in args.control_bam)

    # TODO: check that all bams have the same reference
    first_bam = all_bams[0]
    ref_lengths = None
    with pysam.AlignmentFile(first_bam, "rb") as a:
        ref_lengths = dict(zip(a.references, a.lengths))

    if not os.path.isdir(args.out_dir_plots):
        os.mkdir(args.out_dir_plots)

    if os.path.exists(args.out_dir_plots+'/data_phasing'):
        safe_rmtree(args.out_dir_plots+'/data_phasing')
        os.mkdir(args.out_dir_plots+'/data_phasing')
    else:
        os.mkdir(args.out_dir_plots+'/data_phasing')

    if os.path.exists(args.out_dir_plots+'/coverage_data'):
        safe_rmtree(args.out_dir_plots+'/coverage_data')
        os.mkdir(args.out_dir_plots+'/coverage_data')
    else:
        os.mkdir(args.out_dir_plots+'/coverage_data')

    thread_pool = Pool(args.threads)

    if args.histogram_coverage or args.tumor_phased_vcf:
        segments_by_read = defaultdict(list)
        genome_ids = []
        for bam_file in all_bams:
            genome_id = os.path.basename(bam_file)
            genome_ids.append(genome_id)
            logger.info("Parsing reads from %s", args.target_bam[0])
            segments_by_read_bam = get_all_reads_parallel(bam_file, thread_pool, ref_lengths,
                                                          args.min_mapping_quality, genome_id, MIN_SV_SIZE)
            segments_by_read.update(segments_by_read_bam)
            logger.info("Parsed %s segments", len(segments_by_read_bam))

        logger.info('Computing coverage histogram')
        coverage_histograms = update_coverage_hist(genome_ids, ref_lengths, segments_by_read, args.min_mapping_quality,
                                                   args.max_read_error, args)
        del segments_by_read

    if args.tumor_phased_vcf and not args.quick_start:
        segments = get_chromosomes_bins(args.target_bam[0], args.bin_size, args)
        segments_coverage = get_segments_coverage(segments, coverage_histograms)
        logger.info('Writing tumor coverage for bins')
        write_segments_coverage(segments_coverage, 'coverage_tumor.csv', args)

    #cancer_genes_df_all = []
    #logger.info('Computing coverage for genes')
    #genes_segments = genes_segments_list(args.target_bam[0], args)
    #genes_coverage = get_segments_coverage(genes_segments, coverage_histograms)
    #ancer_genes_df_all = genes_segments_coverage(genes_coverage, args)

    chroms = get_contigs_list(args.contigs)

    logger.info('Parsing phaseblocks information')
    if args.normal_phased_vcf:
        output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(args.normal_phased_vcf, args)
    else:
        output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(args.tumor_phased_vcf, args)
    phasesets_segments = generate_phasesets_bins(args.target_bam[0], output_phasesets_file_path, args.bin_size, args)  # TODO update for multiple bam files

    logger.info('Computing coverage for bins')
    segments = get_chromosomes_bins(args.target_bam[0], args.bin_size, args)

    if args.quick_start:
        csv_df_phasesets = csv_df_chromosomes_sorter(args.quick_start_coverage_path + '/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
        csv_df_coverage = csv_df_chromosomes_sorter(args.quick_start_coverage_path + '/coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
    elif args.histogram_coverage:
        segments_coverage = get_segments_coverage(segments, coverage_histograms)
        logger.info('Writing coverage for bins')
        write_segments_coverage(segments_coverage, 'coverage.csv', args)

        #if args.breakpoints:
        #    phasesets_segments = add_breakpoints(args, phasesets_segments, breakpoints_additional)
        logger.info('Computing coverage for phaseblocks')
        phasesets_coverage = get_segments_coverage(phasesets_segments, coverage_histograms)

        logger.info('Writing coverage for phaseblocks')
        write_segments_coverage(phasesets_coverage, 'coverage_ps.csv', args)

        logger.info('Loading coverage (bins) and coverage (phaseblocks) files...')
        csv_df_phasesets = csv_df_chromosomes_sorter(args.out_dir_plots+'/coverage_data/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
        csv_df_coverage = csv_df_chromosomes_sorter(args.out_dir_plots+'/coverage_data/coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

        #Missing phaseblocks and coverage added
        #phasesets_coverage_missing = update_phasesets_coverage_with_missing_phasesets(chroms, csv_df_phasesets, args.target_bam[0], coverage_histograms)
        #write_segments_coverage(phasesets_coverage_missing, 'coverage_ps_missing.csv')
        #csv_df_phasesets_missing = csv_df_chromosomes_sorter('data/coverage_ps_missing.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

        del coverage_histograms
    else:
        csv_df_coverage = pd.DataFrame([sublist[1:4] for sublist in segments], columns=["chr", "start", "end"])
        csv_df_phasesets = pd.DataFrame([sublist[1:4] for sublist in phasesets_segments], columns=["chr", "start", "end"])

    csv_df_phasesets = csv_df_phasesets[(csv_df_phasesets["end"] - csv_df_phasesets["start"]) >= args.min_phaseblock]

    get_snp_frequencies_segments(args, args.target_bam[0], thread_pool)
    df_snps_frequencies = csv_df_chromosomes_sorter(args.out_dir_plots+'/data_phasing/snps_frequencies.csv', ['chr', 'pos', 'freq_value_a', 'hp_a', 'freq_value_b', 'hp_b'])

    cancer_genes_df_all = snps_frequencies_chrom_genes(df_snps_frequencies, args)

    if args.tumor_phased_vcf:
        output_phasesets_file_path = vcf_parse_to_csv_for_snps(args.tumor_phased_vcf, args)
        df_snps_in_csv_loh = csv_df_chromosomes_sorter(output_phasesets_file_path, ['chr', 'pos', 'qual', 'gt', 'dp', 'vaf'])
        #df_snps_frequencies_homo = csv_df_chromosomes_sorter(args.out_dir_plots + '/data_phasing/snps_frequencies_homo.csv', ['chr', 'pos', 'freq_value_a', 'hp_a', 'freq_value_b', 'hp_b'])

        if args.quick_start:
            csv_df_coverage_tumor = csv_df_chromosomes_sorter(args.quick_start_coverage_path + '/' + 'coverage_tumor.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
        else:
            csv_df_coverage_tumor = csv_df_chromosomes_sorter(args.out_dir_plots + '/coverage_data/' + 'coverage_tumor.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

    if not os.path.isdir(args.out_dir_plots + '/phasing_output'):
        os.mkdir(args.out_dir_plots + '/phasing_output')

    filename = f"{os.path.join(args.out_dir_plots, 'phasing_output', 'PHASE_CORRECTION_INDEX.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")

    start_values_phasesets_contiguous_all = []
    loh_regions_events_all = []
    df_updated_coverage = []
    df_coverage = []
    df_phasesets = []
    df_snps_ratios = []
    offset = 0
    for index, chrom in enumerate(chroms):
        if chrom in chroms:
            regions = get_chromosomes_regions(args)
            logger.info('Loading coverage (bins) and coverage (phaseblocks) datasets for ' + chrom)
            csv_df_phaseset = csv_df_phasesets[csv_df_phasesets['chr'] == chrom]
            ref_start_values_phasesets = csv_df_phaseset.start.values.tolist()
            ref_end_values_phasesets = csv_df_phaseset.end.values.tolist()
            if args.histogram_coverage:
                haplotype_1_values_phasesets = csv_df_phaseset.hp1.values.tolist()
                haplotype_2_values_phasesets = csv_df_phaseset.hp2.values.tolist()

            centromere_region = centromere_regions[centromere_regions['chr'] == chrom]

            csv_df_coverage_chrom = csv_df_coverage[csv_df_coverage['chr'] == chrom]
            ref_start_values = csv_df_coverage_chrom.start.values.tolist()
            ref_end_values = csv_df_coverage_chrom.end.values.tolist()
            loh_chrom = pd.DataFrame({'chr': [], 'start': [], 'end': []})

            if args.histogram_coverage:
                unphased_reads_values = csv_df_coverage_chrom.hp3.values.tolist()
                haplotype_1_values = csv_df_coverage_chrom.hp1.values.tolist()
                haplotype_2_values = csv_df_coverage_chrom.hp2.values.tolist()
            else:
                unphased_reads_values = [0 for i in range(len(ref_start_values))]

            #ref_start_values_phasesets, ref_end_values_phasesets = add_breakpoints(ref_start_values_phasesets, ref_end_values_phasesets, breakpoints)
            if args.histogram_coverage:
                snps_haplotype1_mean = haplotype_1_values
                snps_haplotype2_mean = haplotype_2_values
            else:
                snps_haplotype1_mean, snps_haplotype2_mean  = snps_frequencies_chrom_mean(df_snps_frequencies, ref_start_values, chrom, args)
                haplotype_1_values_phasesets, haplotype_2_values_phasesets = \
                        snps_frequencies_chrom_mean_phasesets(df_snps_frequencies, ref_start_values_phasesets, ref_end_values_phasesets, chrom, args)

            #TODO create coverage and phaseset-coverage CSVs
            if not args.histogram_coverage:
                ##################################
                df_coverage.append(pd.DataFrame(list(zip([chrom for ch in range(len(ref_start_values))], 
                                                ref_start_values, ref_end_values, snps_haplotype1_mean, 
                                                snps_haplotype2_mean, unphased_reads_values)), columns=['chr', 'start', 'end', 'hp1', 'hp2', 'hp3']))
                df_phasesets.append(pd.DataFrame(list(zip([chrom for ch in range(len(ref_start_values_phasesets))], 
                                                 ref_start_values_phasesets, ref_end_values_phasesets, 
                                                 haplotype_1_values_phasesets, haplotype_2_values_phasesets)), columns=['chr', 'start', 'end', 'hp1', 'hp2']))
                ##################################
            plot_coverage_data(html_graphs, args, chrom, ref_start_values, ref_end_values, 
                               snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values, 
                               haplotype_1_values_phasesets, haplotype_2_values_phasesets, 
                               ref_start_values_phasesets, ref_end_values_phasesets, "without_phase_correction")
            ##################################
            if args.tumor_phased_vcf:
                (ref_start_values_updated, snps_het_counts, snps_homo_counts,
                centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends) = \
                    get_snps_frquncies_coverage(df_snps_in_csv_loh, chrom, ref_start_values, args.bin_size_snps,
                                                args.hets_ratio, args.hets_smooth_window, args)
                df_snps_ratios_chrom = extend_snps_ratios_df(chrom, offset, ref_start_values_updated, snps_het_counts, snps_homo_counts)
                df_snps_ratios.append(df_snps_ratios_chrom)
                offset += regions[index]

                if args.without_phasing:
                    detect_loh_centromere_regions(None, chrom, args, centromere_region_starts, centromere_region_ends,
                                                  loh_region_starts, loh_region_ends, ref_start_values, 
                                                  ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean, 
                                                  unphased_reads_values, haplotype_1_values_phasesets, 
                                                  haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
                else:
                    csv_df_coverage_tumor_chrom = csv_df_coverage_tumor[csv_df_coverage_tumor['chr'] == chrom]
                    (snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values,
                     haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                     ref_start_values_phasesets, ref_end_values_phasesets, loh_region_starts, loh_region_ends, hp) = \
                             detect_loh_centromere_regions(csv_df_coverage_tumor_chrom, chrom, args, centromere_region_starts,
                                                           centromere_region_ends, loh_region_starts, loh_region_ends,
                                                           ref_start_values, ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean,
                                                           unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                                                           ref_start_values_phasesets, ref_end_values_phasesets)
                    if len(haplotype_1_values_phasesets):
                        (haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                         ref_start_values_phasesets, ref_end_values_phasesets) = \
                                 adjust_loh_cent_phaseblocks(args, loh_region_starts, loh_region_ends, centromere_region,
                                                             snps_haplotype1_mean, snps_haplotype2_mean, haplotype_1_values_phasesets,
                                                             haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
                    if len(loh_region_starts):
                        chr_list = [chrom for ch in range(len(loh_region_starts))]
                        loh_chrom = pd.DataFrame(list(zip(chr_list, loh_region_starts, loh_region_ends, hp)), columns=['chr', 'start', 'end', 'hp'])
                        loh_regions_events_all.append(loh_chrom)
            else:
               loh_region_starts = []
               loh_region_ends = []

            updated_df_phasesets = pd.DataFrame(list(zip(ref_start_values_phasesets, ref_end_values_phasesets,
                                                         haplotype_1_values_phasesets, haplotype_2_values_phasesets)),
                                                columns=['start', 'end', 'hp1', 'hp2'])

            #phaseblocks flipping main algo
            is_simple_heuristics = True
            #is_phasesets_check_simple_heuristics(ref_start_values_phasesets, ref_end_values_phasesets)
            if len(ref_start_values_phasesets) >= 1:
                is_simple_heuristics = False
            if args.enable_simple_heuristics:
                is_simple_heuristics = True

            if is_simple_heuristics:
                # #plot resultant
                (snps_haplotype1_mean, snps_haplotype2_mean, haplotype_1_values_phasesets, 
                 haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets) = \
                    phaseblock_flipping_simple_heuristics(chrom, args, is_simple_heuristics, snps_haplotype1_mean, 
                                                          snps_haplotype2_mean, ref_start_values, ref_start_values, 
                                                          haplotype_1_values_phasesets, haplotype_2_values_phasesets, 
                                                          ref_start_values_phasesets, ref_end_values_phasesets)
                plot_coverage_data(html_graphs, args, chrom, ref_start_values, ref_end_values, 
                                   snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values,
                                   haplotype_1_values_phasesets, haplotype_2_values_phasesets, 
                                   ref_start_values_phasesets, ref_end_values_phasesets, "phase_correction_1")
            else:
                # #detect centromeres
                #ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = detect_centromeres(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, snps_haplotype1_mean, snps_haplotype2_mean, ref_start_values, args.bin_size'])
                #infer missing phaseblocks
                #ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = infer_missing_phaseblocks(ref_start_values, ref_end_values, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, snps_haplotype1_mean, snps_haplotype2_mean, args.bin_size'])

                haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = \
                        remove_centromere_phaseblocks(haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                                                      ref_start_values_phasesets, ref_end_values_phasesets, centromere_region)

                breakpoints_chrom = []
                (broken_phasesets, mean_cis_trans_ps, snps_haplotype1_mean,
                 snps_haplotype2_mean, haplotype_1_values_phasesets,
                 haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets) = \
                         phase_flips_cis_trans(chrom, args, breakpoints_chrom, snps_haplotype1_mean, snps_haplotype2_mean,
                                               ref_start_values, ref_end_values, haplotype_1_values_phasesets,
                                               haplotype_2_values_phasesets, ref_start_values_phasesets,
                                               ref_end_values_phasesets, True, False, False, False, False)
                                               #internal_ps=False, bins_adjust=False, merge=False, merge_final=False
                plot_coverage_data(html_graphs, args, chrom, ref_start_values, ref_end_values,
                                   snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values,
                                   haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                                   ref_start_values_phasesets, ref_end_values_phasesets, "phase_correction_0.1")
                # switch inter phaseblocks bins
                snps_haplotype1_mean, snps_haplotype2_mean = \
                        switch_inter_phaseblocks_bins(chrom, args, ref_start_values, snps_haplotype1_mean, snps_haplotype2_mean,
                                                      ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets,
                                                      haplotype_2_values_phasesets)
                plot_coverage_data(html_graphs, args, chrom, ref_start_values, ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean,
                                   unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                                   ref_start_values_phasesets, ref_end_values_phasesets, "phase_correction_0.2")

                ###########################################################
                #merge ps
                from src.hapcorrect.src.phase_correction import merge_contiguous_indices, find_indices_to_be_merged
                indices_merge = find_indices_to_be_merged(mean_cis_trans_ps, ref_start_values_phasesets,
                                                          ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets)
                ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = \
                        merge_contiguous_indices(indices_merge, snps_haplotype1_mean, snps_haplotype2_mean, ref_start_values,
                                                 ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets,
                                                 haplotype_2_values_phasesets, False)
                plot_coverage_data(html_graphs, args, chrom, ref_start_values, ref_end_values, snps_haplotype1_mean,
                                   snps_haplotype2_mean, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                                   ref_start_values_phasesets, ref_end_values_phasesets, "phase_correction_0.3")
                ###########################################################
                #snps_haplotype1_mean, snps_haplotype2_mean = switch_inter_phaseblocks_bins(chrom, args, ref_start_values, snps_haplotype1_mean, snps_haplotype2_mean, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets)

                #ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = update_remaining_phasesets(indices_merge, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets_merged, ref_end_values_phasesets_merged, haplotype_1_values_phasesets_merged, haplotype_2_values_phasesets_merged)
                (_, mean_cis_trans_ps, snps_haplotype1_mean, snps_haplotype2_mean,
                 haplotype_1_values_phasesets_merged, haplotype_2_values_phasesets_merged,
                 ref_start_values_phasesets_merged, ref_end_values_phasesets_merged) = \
                         phase_flips_cis_trans(chrom, args, breakpoints_chrom, snps_haplotype1_mean, snps_haplotype2_mean,
                                               ref_start_values, ref_end_values, haplotype_1_values_phasesets,
                                               haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets,
                                               False, False, False, True, True) #internal_ps=False, bins_adjust=False, merge=False, merge_final=False
                if broken_phasesets:
                    (haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                     ref_start_values_phasesets, ref_end_values_phasesets) = \
                             reintroduce_broken_phasesets(broken_phasesets, ref_start_values, snps_haplotype1_mean, snps_haplotype2_mean,
                                                          updated_df_phasesets.start.values.tolist(), updated_df_phasesets.end.values.tolist(),
                                                          ref_start_values_phasesets_merged, ref_end_values_phasesets_merged,
                                                          haplotype_1_values_phasesets_merged, haplotype_2_values_phasesets_merged)
            snps_haplotype1_mean, snps_haplotype2_mean = \
                    without_phasesets_bins_correction(loh_chrom, args, chrom, ref_start_values, snps_haplotype1_mean, snps_haplotype2_mean,
                                                      ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets,
                                                      haplotype_2_values_phasesets)
            plot_coverage_data(html_graphs, args, chrom, ref_start_values, ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean,
                               unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                               ref_start_values_phasesets, ref_end_values_phasesets, "phase_correction_1")
            from src.smoothing import smoothing
            snps_haplotype1_mean_smoothed, snps_haplotype2_mean_smoothed, unphased_reads_values_smoothed = \
                    smoothing(snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values, conv_window_size=5)
            plot_coverage_data(html_graphs, args, chrom, ref_start_values, ref_end_values, snps_haplotype1_mean_smoothed,
                               snps_haplotype2_mean_smoothed, unphased_reads_values_smoothed, haplotype_1_values_phasesets,
                               haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, "phase_correction_2")
            ##################################
            chr_list = [chrom for ch in range(len(snps_haplotype1_mean))]
            df_updated_coverage.append(pd.DataFrame(list(zip(chr_list, ref_start_values, ref_end_values, snps_haplotype1_mean,
                                                             snps_haplotype2_mean, unphased_reads_values)),
                                                    columns=['chr', 'start', 'end', 'hp1', 'hp2', 'hp3']))
            ##################################
            chr_list = [chrom for ch in range(len(ref_start_values_phasesets))]

            updated_ps = pd.DataFrame(list(zip(chr_list, ref_start_values_phasesets, ref_end_values_phasesets)), columns=['chr', 'start', 'end'])
            if args.tumor_phased_vcf and not loh_chrom.empty and not updated_ps.empty:
                updated_ps = subtract_intervals(updated_ps, loh_chrom)
            updated_ps = check_missing_phasesets_original(updated_ps, csv_df_phaseset)
            updated_ps = updated_ps.drop_duplicates()
            updated_ps_haplotype_1_values_phasesets, updated_ps_haplotype_2_values_phasesets = \
                    snps_frequencies_chrom_mean_phasesets(df_snps_frequencies, updated_ps.start.values.tolist(), updated_ps.end.values.tolist(), chrom, args)
            updated_ps['hp1'] = updated_ps_haplotype_1_values_phasesets
            updated_ps['hp2'] = updated_ps_haplotype_2_values_phasesets
            start_values_phasesets_contiguous_all.append(updated_ps)

    html_graphs.write("</body></html>")

    updated_ps = pd.concat(start_values_phasesets_contiguous_all)
    write_df_csv(updated_ps[['chr', 'start', 'end']], args.out_dir_plots + '/data_phasing/' + args.genome_name + '_phasesets.csv')
    write_df_csv(updated_ps, args.out_dir_plots + '/coverage_data/' + 'coverage_ps.csv')

    write_df_csv(pd.concat(df_updated_coverage), args.out_dir_plots+'/coverage_data/phase_corrected_coverage.csv')
    if not args.histogram_coverage:
        write_df_csv(pd.concat(df_coverage), args.out_dir_plots + '/coverage_data/coverage.csv')
        write_df_csv(pd.concat(df_phasesets), args.out_dir_plots + '/coverage_data/coverage_ps.csv')
    if not cancer_genes_df_all.empty:
        write_df_csv(cancer_genes_df_all, args.out_dir_plots + '/coverage_data/cancer_genes_coverage.csv')
    else:
        write_df_csv(pd.DataFrame(columns=['chr', 'start', 'end', 'gene', 'hp1', 'hp2']), args.out_dir_plots + '/coverage_data/cancer_genes_coverage.csv')

    if args.tumor_phased_vcf and len(loh_regions_events_all):
        loh_df_final = pd.concat(loh_regions_events_all)
        loh_df_final_filtered = loh_df_final[loh_df_final['end'] - loh_df_final['start'] >= 20000]
        #update HP in LOH regions for flipped segments
        loh_df = update_hp_assignment_loh_segments(args, loh_df_final_filtered, pd.concat(df_updated_coverage))
        write_df_csv(loh_df, args.out_dir_plots+'/coverage_data/'+args.genome_name+'_loh_segments.csv')
        csv_df_loh_regions = csv_df_chromosomes_sorter(args.out_dir_plots + '/coverage_data/'+args.genome_name + '_loh_segments.csv', ['chr', 'start', 'end', 'hp'])
    else:
        csv_df_loh_regions = pd.DataFrame(columns=['chr', 'start', 'end', 'hp'])

    if os.path.isfile(args.out_dir_plots+'/data_phasing/' + args.genome_name + '_phase_change_segments.csv'):
        csv_df_phase_change_segments = csv_df_chromosomes_sorter(args.out_dir_plots+'/data_phasing/' + args.genome_name + '_phase_change_segments.csv', ['chr', 'start', 'end'])
        csv_df_phasesets_segments = csv_df_chromosomes_sorter(args.out_dir_plots+'/data_phasing/' + args.genome_name + '_phasesets.csv', ['chr', 'start', 'end'])

        if args.normal_phased_vcf:
            get_phasingblocks(args.normal_phased_vcf)
            logger.info('VCF edit for phase change segments')
            out_vcf = os.path.join(args.out_dir_plots, 'phasing_output', 'rephased.vcf.gz')
            rephase_vcf(csv_df_phase_change_segments, csv_df_phasesets_segments,
                        csv_df_loh_regions, args.normal_phased_vcf, out_vcf, args)
            index_vcf(out_vcf)
            get_phasingblocks(out_vcf)

        elif not args.normal_phased_vcf and args.tumor_phased_vcf:
            get_phasingblocks(args.tumor_phased_vcf)
            logger.info('VCF edit for phase change segments')
            out_vcf = os.path.join(args.out_dir_plots, 'phasing_output', 'rephased.vcf.gz')
            rephase_vcf(csv_df_phase_change_segments, csv_df_phasesets_segments,
                        csv_df_loh_regions, args.tumor_phased_vcf, out_vcf, args)
            index_vcf(out_vcf)
            get_phasingblocks(out_vcf)

        if args.rehaplotag_tumor_bam:
            logger.info('Rehaplotagging tumor BAM')
            tumor_bam_haplotag(args, out_vcf)
    else:
        logger.info('Note: No phase-change segments detected for phase correction heauristics, in output bins output, it could be handled as simple heuristics, high coverage as HP-1 and low coverage bins as HP-2.')
        logger.info('Note: Also, no updated phased VCF generated.')

    if args.tumor_phased_vcf:
        loh_plots_genome(pd.concat(df_snps_ratios),  args, csv_df_loh_regions)
        plot_snps(args, df_snps_in_csv_loh)
