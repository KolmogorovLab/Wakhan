#!/usr/bin/env python3

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

from hapcorrect.src.process_bam import get_all_reads_parallel, update_coverage_hist, get_segments_coverage, haplotype_update_all_bins_parallel, get_snps_frequencies, tumor_bam_haplotag
from hapcorrect.src.process_vcf import vcf_parse_to_csv_for_het_phased_snps_phasesets, get_snp_frequencies_segments, snps_frequencies_chrom_mean, get_snps_frquncies_coverage, vcf_parse_to_csv_for_snps, index_vcf, rephase_vcf, get_phasingblocks
from hapcorrect.src.phase_correction import generate_phasesets_bins, phaseblock_flipping, phase_correction_centers, contiguous_phaseblocks, detect_centromeres, flip_phaseblocks_contigous, remove_overlaping_contiguous, switch_inter_phaseblocks_bins
from hapcorrect.src.utils import get_chromosomes_bins, write_segments_coverage, csv_df_chromosomes_sorter, get_snps_frquncies_coverage_from_bam, \
                    infer_missing_phaseblocks, df_chromosomes_sorter, is_phasesets_check_simple_heuristics, write_df_csv, loh_regions_events
from hapcorrect.src.extras import get_contigs_list
from hapcorrect.src.plots import plot_coverage_data, change_point_detection, plot_coverage_data_after_correction
from hapcorrect.src.cpd import cpd_positions_means
from hapcorrect.src.loh import detect_loh_centromere_regions, plot_snps

def main_process():
    # default tunable parameters
    MAX_READ_ERROR = 0.1
    MIN_MAPQ = 10
    MIN_SV_SIZE = 50
    BIN_SIZE = 50000
    BIN_SIZE_SNPS = 50000
    MAX_CUT_THRESHOLD = 100
    MIN_ALIGNED_LENGTH = 5000
    MAX_CUT_THRESHOLD_SNPS_COUNTS = 50

    SAMTOOLS_BIN = "samtools"
    BCFTOOLS_BIN = "bcftools"

    DEFAULT_CONTIGS = 'chr1-22' #('chr1-22' '1-22') ('chr1-22,chrX' '1-22,X')

    parser = argparse.ArgumentParser \
        (description="Plot coverage and copy number profiles from a bam and phased VCF files")

    parser.add_argument("--target-bam", dest="target_bam",
                        metavar="path", required=True, default=None, nargs="+",
                        help="path to one or 1/4 multiple target haplotagged bam files (must be indexed)")
    parser.add_argument("--control-bam", dest="control_bam",
                        metavar="path", required=False, default=None, nargs="+",
                        help="path to one or multiple control haplotagged bam files (must be indexed)")
    parser.add_argument("--reference", dest="reference",
                        metavar="path", required=False, default=None,
                        help="path to reference")

    parser.add_argument("--out-dir-plots", dest="out_dir_plots",
                        default=None, required=True,
                        metavar="path", help="Output directory")

    parser.add_argument("--normal-phased-vcf", dest="normal_phased_vcf",
                        metavar="path", required=False, default=None,
                        help="Path to normal phased vcf")
    parser.add_argument("--tumor-vcf", dest="tumor_vcf",
                        metavar="path", required=False, default=None,
                        help="Path to tumor VCF for LOH detection")

    parser.add_argument("--breakpoints", dest="breakpoints", metavar="path", required=True, default=None, help="Path to breakpoints/SVs VCF file")
    parser.add_argument("--cpd-internal-segments", dest="cpd_internal_segments", metavar="path", required=False, default=None,
                        help="change point detection algo for more precise segments after breakpoint segments")

    parser.add_argument("--genome-name", dest="genome_name",
                        required=True, default=None,
                        help="Genome sample/cellline name to be displayed on plots")
    parser.add_argument("--contigs", dest="contigs",
                        required=False, default=DEFAULT_CONTIGS,
                        help="List of contigs (choromosomes) to be included in the plots [e.g., chr1-22,X,Y]")

    parser.add_argument("--bin-size", "--bin_size", dest="bin_size",
                        default=BIN_SIZE, metavar="int", type=int, help="coverage (readdepth) bin size [50k]")
    parser.add_argument("--bin-size-snps", "--bin_size_snps", dest="bin_size_snps",
                        default=BIN_SIZE_SNPS, metavar="int", type=int, help="SNPs bin size [50k]")

    parser.add_argument("--cut-threshold", "--cut_threshold", dest="cut_threshold",
                        default=MAX_CUT_THRESHOLD, metavar="int", type=int, help="Maximum cut threshold for coverage (readdepth) [100]")
    parser.add_argument("--cut-threshold-snps-counts", "--cut_threshold_snps_counts", dest="cut_threshold_snps_counts",
                        default=MAX_CUT_THRESHOLD_SNPS_COUNTS, metavar="int", type=int,
                        help="Maximum cut threshold for SNPs counts [50]")

    parser.add_argument("--no-of-clusters", dest="no_of_clusters",
                        required=False, default=None, metavar="int", type=int,
                        help="Number of clusters for bins clustering")

    parser.add_argument("--min-aligned-length", "--min_aligned_length", dest="min_aligned_length",
                        default=MIN_ALIGNED_LENGTH, metavar="int", type=int, help="Minimum aligned reads length [5000]")

    parser.add_argument('--pdf-enable',  dest="pdf_enable", required=False,
                        default=False, help="Enabling PDF output coverage plots")
    parser.add_argument('--html-enable',  dest="html_enable", required=False,
                        default=True, help="Enabling HTML output coverage plots")

    parser.add_argument('--unphased-reads-coverage-enable',  dest="unphased_reads_coverage_enable", required=False,
                        default=True, help="Enabling unphased reads coverage output in plots")
    parser.add_argument('--without-phasing', dest="without_phasing", required=False,
                        default=False, help="Enabling coverage and copynumbers without phasing in plots")

    parser.add_argument('--phaseblock-flipping-enable',  dest="phaseblock_flipping_enable", required=False,
                        default=True, help="Enabling phaseblock flipping in coverage plots")
    parser.add_argument('--smoothing-enable',  dest="smoothing_enable", required=False,
                        default=False, help="Enabling smoothing in coverage plots")
    parser.add_argument('--phaseblocks-enable',  dest="phaseblocks_enable", required=False,
                        default=False, help="Enabling phaseblocks display in coverage plots")

    parser.add_argument('--het-phased-snps-freq-enable',  dest="het_phased_snps_freq_enable", required=False,
                        default=False, help="Enabling hetrozygous phased snps frequencies in coverage plots")
    parser.add_argument('--snps-freq-vcf-enable',  dest="snps_freq_vcf_enable", required=False,
                        default=False, help="Enabling snps frequencies in coverage plots")
    parser.add_argument('--breakpoints-enable',  dest="breakpoints_enable", required=False,
                        default=False, help="Enabling breakpoints in coverage plots")
    parser.add_argument('--copynumbers-enable', dest="copynumbers_enable", required=False,
                        default=True, help="Enabling copy number in coverage plots")
    parser.add_argument('--copynumbers-subclonal-enable', dest="copynumbers_subclonal_enable", required=False,
                        default=False, help="Enabling subclonal copy number in coverage plots")

    parser.add_argument('--enable-debug', dest="enable_debug", required=False,
                        default=False, help="Enabling debug")

    parser.add_argument('--variable-size-bins', dest="variable_size_bins", required=False,
                        default=False, help="enable variable size bins to use breakpoints")
    parser.add_argument('--enable-simple-heuristics', dest="enable_simple_heuristics", required=False,
                        default=False, help="enable simple heuristics")

    parser.add_argument('--rephase-normal-vcf', dest="rephase_normal_vcf", required=False,
                        default=False, help="enable rephase normal vcf")
    parser.add_argument('--rephase-tumor-vcf', dest="rephase_tumor_vcf", required=False,
                        default=False, help="enable rephase tumor vcf")
    parser.add_argument('--rehaplotag-tumor-bam', dest="rehaplotag_tumor_bam", required=False,
                        default=False, help="enable rehaplotag tumor bam")

    parser.add_argument("--bins-cluster-means", dest="bins_cluster_means",
                        default=None, required=False, type=lambda s: [int(item) for item in s.split(',')],
                        help="bins cluster means")

    parser.add_argument("-t", "--threads", dest="threads",
                        default=1, metavar="int", type=int, help="number of parallel threads [8]")
    parser.add_argument('--dryrun', dest="dryrun", required=False,
                        default=False, help="Enabling dryrun")
    parser.add_argument("--dryrun-path", dest="dryrun_path",
                        default=None, required=False,
                        metavar="path", help="dryrun data directory")

    parser.add_argument("--max-read-error", dest="max_read_error",
                        default=MAX_READ_ERROR, metavar="float", type=float,
                        help=f"maximum base alignment error [{MAX_READ_ERROR}]")
    parser.add_argument("--min-mapq", dest="min_mapping_quality",
                        default=MIN_MAPQ, metavar="int", type=int,
                        help=f"minimum mapping quality for aligned segment [{MIN_MAPQ}]")

    args = parser.parse_args()

    arguments = {
        "target_bam": args.target_bam,
        "reference": args.reference,
        "phaseblock_flipping_enable": args.phaseblock_flipping_enable,
        "unphased_reads_coverage_enable": args.unphased_reads_coverage_enable,
        "without_phasing": args.without_phasing,
        "smoothing_enable": args.smoothing_enable,
        "phaseblocks_enable": args.phaseblocks_enable,
        "copynumbers_enable": args.copynumbers_enable,
        "copynumbers_subclonal_enable": args.copynumbers_subclonal_enable,
        "het_phased_snps_freq_enable": args.het_phased_snps_freq_enable,
        "snps_freq_vcf_enable": args.snps_freq_vcf_enable,
        "out_dir_plots": args.out_dir_plots,
        "normal_phased_vcf": args.normal_phased_vcf,
        "tumor_vcf": args.tumor_vcf,
        "breakpoints": args.breakpoints,
        "cpd_internal_segments": args.cpd_internal_segments,
        "breakpoints_enable": args.breakpoints_enable,
        "genome_name": args.genome_name,
        "bin_size": args.bin_size,
        "bins_cluster_means": args.bins_cluster_means,
        "bin_size_snps": args.bin_size_snps,
        "pdf_enable": args.pdf_enable,
        "cut_threshold": args.cut_threshold,
        "cut_threshold_snps_counts": args.cut_threshold_snps_counts,
        "min_aligned_length": args.min_aligned_length,
        "no_of_clusters": args.no_of_clusters,
        "contigs": args.contigs,
        "threads": args.threads,
        "dryrun": args.dryrun,
        "dryrun_path": args.dryrun_path,
        "enable_debug": args.enable_debug,
        "variable_size_bins": args.variable_size_bins,
        "enable_simple_heuristics": args.enable_simple_heuristics,
        "rephase_normal_vcf": args.rephase_normal_vcf,
        "rephase_tumor_vcf": args.rephase_tumor_vcf,
        "rehaplotag_tumor_bam": args.rehaplotag_tumor_bam,

    }
    logging.basicConfig(level=logging.DEBUG)

    if args.control_bam is None:
        args.control_bam = []
    all_bams = args.target_bam + args.control_bam
    target_genomes = set(os.path.basename(b) for b in args.target_bam)
    control_genomes = set(os.path.basename(b) for b in args.control_bam)

    if not shutil.which(SAMTOOLS_BIN):
        print("samtools not found", file=sys.stderr)
        return 1

    # TODO: check that all bams have the same reference
    first_bam = all_bams[0]
    ref_lengths = None
    with pysam.AlignmentFile(first_bam, "rb") as a:
        ref_lengths = dict(zip(a.references, a.lengths))

    if not os.path.isdir(args.out_dir_plots):
        os.mkdir(args.out_dir_plots)
    #if not os.path.isdir('data'):
    if os.path.exists(arguments['out_dir_plots']+'/data_phasing'):
        shutil.rmtree(arguments['out_dir_plots']+'/data_phasing')
        os.mkdir(arguments['out_dir_plots']+'/data_phasing')
    else:
        os.mkdir(arguments['out_dir_plots']+'/data_phasing')

    thread_pool = Pool(args.threads)

    segments_by_read = defaultdict(list)
    genome_ids = []
    for bam_file in all_bams:
        genome_id = os.path.basename(bam_file)
        genome_ids.append(genome_id)
        print("Parsing reads from", genome_id, file=sys.stderr)
        segments_by_read_bam = get_all_reads_parallel(bam_file, thread_pool, ref_lengths,
                                                      args.min_mapping_quality, genome_id, MIN_SV_SIZE)
        segments_by_read.update(segments_by_read_bam)
        print("Parsed {0} segments".format(len(segments_by_read_bam)), file=sys.stderr)

    logging.info('Computing coverage histogram')
    coverage_histograms = update_coverage_hist(genome_ids, ref_lengths, segments_by_read, args.min_mapping_quality,
                                               args.max_read_error, arguments)
    del segments_by_read
    chroms = get_contigs_list(arguments['contigs'])

    if arguments['dryrun']:
        csv_df_phasesets = csv_df_chromosomes_sorter(arguments['dryrun_path'] + arguments['genome_name'] + '/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
        #csv_df_phasesets_missing = csv_df_chromosomes_sorter(arguments['dryrun_path'] + arguments['genome_name'] + '/coverage_ps_missing.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
        #csv_df_phasesets = df_chromosomes_sorter(pd.concat([csv_df_phasesets, csv_df_phasesets_missing], ignore_index=True), ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

        csv_df_coverage = csv_df_chromosomes_sorter(arguments['dryrun_path'] + arguments['genome_name'] + '/coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
    else:
        logging.info('Computing coverage for bins')
        segments = get_chromosomes_bins(args.target_bam[0], arguments['bin_size'], arguments)
        segments_coverage = get_segments_coverage(segments, coverage_histograms)
        logging.info('Writing coverage for bins')
        write_segments_coverage(segments_coverage, 'coverage.csv', arguments)

        logging.info('Parsing phaseblocks information')
        if arguments['normal_phased_vcf']:
            output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(arguments['normal_phased_vcf'], arguments)
        else:
            output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(arguments['tumor_vcf'], arguments)
        phasesets_segments = generate_phasesets_bins(args.target_bam[0], output_phasesets_file_path, arguments['bin_size'], arguments) #TODO update for multiple bam files
        logging.info('Computing coverage for phaseblocks')
        phasesets_coverage = get_segments_coverage(phasesets_segments, coverage_histograms)

        logging.info('Writing coverage for phaseblocks')
        write_segments_coverage(phasesets_coverage, 'coverage_ps.csv', arguments)

        logging.info('Loading coverage (bins) and coverage (phaseblocks) files...')
        csv_df_phasesets = csv_df_chromosomes_sorter(arguments['out_dir_plots']+'/data_phasing/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
        csv_df_coverage = csv_df_chromosomes_sorter(arguments['out_dir_plots']+'/data_phasing/coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

        #Missing phaseblocks and coverage added
        #phasesets_coverage_missing = update_phasesets_coverage_with_missing_phasesets(chroms, csv_df_phasesets, args.target_bam[0], coverage_histograms)
        #write_segments_coverage(phasesets_coverage_missing, 'coverage_ps_missing.csv')
        #csv_df_phasesets_missing = csv_df_chromosomes_sorter('data/coverage_ps_missing.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

        del coverage_histograms

    if arguments['normal_phased_vcf']:
        get_snp_frequencies_segments(arguments, arguments['target_bam'][0], thread_pool)
        df_snps_frequencies = csv_df_chromosomes_sorter(arguments['out_dir_plots']+'/data_phasing/snps_frequencies.csv', ['chr', 'pos', 'freq_value_a', 'hp_a', 'freq_value_b', 'hp_b'])
        df_snps_frequencies = df_snps_frequencies.drop(df_snps_frequencies[(df_snps_frequencies.chr == "chrY")].index)

    if arguments['tumor_vcf']:
        output_phasesets_file_path = vcf_parse_to_csv_for_snps(arguments['tumor_vcf'], arguments)
        df_snps_in_csv = csv_df_chromosomes_sorter(output_phasesets_file_path, ['chr', 'pos', 'qual', 'gt', 'dp', 'vaf'])
        # plot all SNPs frequencies, ratios and counts
        #if arguments['enable_debug']:
        #    plot_snps(arguments, df_snps_in_csv)

    if not os.path.isdir(arguments['out_dir_plots'] + '/phasing_output'):
        os.mkdir(arguments['out_dir_plots'] + '/phasing_output')

    filename = f"{os.path.join(arguments['out_dir_plots'], 'phasing_output', 'PHASE_CORRECTION_INDEX.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")
    start_values_phasesets_contiguous_all = []
    loh_regions_events_all = []
    df_updated_coverage = []
    for index, chrom in enumerate(chroms):
        if chrom in chroms:# and (chrom == '1' or chrom == '18'):# and (chrom == 'chr5' or chrom == 'chr16'):
            logging.info('Loading coverage (bins) and coverage (phaseblocks) datasets for ' + chrom)
            csv_df_phaseset = csv_df_phasesets[csv_df_phasesets['chr'] == chrom]
            haplotype_1_values_phasesets = csv_df_phaseset.hp1.values.tolist()
            haplotype_2_values_phasesets = csv_df_phaseset.hp2.values.tolist()
            ref_start_values_phasesets = csv_df_phaseset.start.values.tolist()
            ref_end_values_phasesets = csv_df_phaseset.end.values.tolist()

            csv_df_coverage_chrom = csv_df_coverage[csv_df_coverage['chr'] == chrom]
            unphased_reads_values = csv_df_coverage_chrom.hp3.values.tolist()
            haplotype_1_values = csv_df_coverage_chrom.hp1.values.tolist()
            haplotype_2_values = csv_df_coverage_chrom.hp2.values.tolist()
            ref_start_values = csv_df_coverage_chrom.start.values.tolist()
            ref_end_values = csv_df_coverage_chrom.end.values.tolist()

            if arguments['normal_phased_vcf']:
                snps_haplotype1_mean, snps_haplotype2_mean  = snps_frequencies_chrom_mean(df_snps_frequencies, ref_start_values, chrom, arguments)
            else:
                snps_haplotype1_mean = haplotype_1_values
                snps_haplotype2_mean = haplotype_2_values
            plot_coverage_data(html_graphs, arguments, chrom, ref_start_values, ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, "without_phase_correction")
            ##################################
            if arguments['tumor_vcf']:
                ref_start_values_updated, snps_het_counts, snps_homo_counts, centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends = get_snps_frquncies_coverage(df_snps_in_csv, chrom, ref_start_values, arguments['bin_size'])
                if arguments['without_phasing']:
                    detect_loh_centromere_regions(chrom, arguments, centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends, ref_start_values, ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
                else:
                    snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, loh_region_starts, loh_region_ends = detect_loh_centromere_regions(chrom, arguments, centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends, ref_start_values, ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
                    loh_regions_events_all.extend(loh_regions_events(chrom, loh_region_starts, loh_region_ends))
            else:
                loh_region_starts =[]
                loh_region_ends = []
            ##################################
            # change_point_detection(snps_haplotype1_mean, ref_start_values, ref_end_values, arguments, chrom,
            #                        html_graphs, 1, color='#6A5ACD')
            # change_point_detection(snps_haplotype2_mean, ref_start_values, ref_end_values, arguments, chrom,
            #                        html_graphs, 2, color='#2E8B57')

            # cpd_haplotype1_mean, cpd_haplotype1_start, cpd_haplotype1_end,  cpd_haplotype2_mean, cpd_haplotype2_start, cpd_haplotype2_end = cpd_positions_means(snps_haplotype1_mean, snps_haplotype2_mean, arguments)

            # change_point_detection(snps_haplotype1_mean, ref_start_values, ref_end_values, arguments, chrom,
            #                        html_graphs, 1, color='#6A5ACD')
            # change_point_detection(snps_haplotype2_mean, ref_start_values, ref_end_values, arguments, chrom,
            #                        html_graphs, 2, color='#2E8B57')

            # phase_correction_centers(arguments, cpd_haplotype1_mean, cpd_haplotype1_start, cpd_haplotype1_end,  cpd_haplotype2_mean, cpd_haplotype2_start, cpd_haplotype2_end, snps_haplotype1_mean, snps_haplotype2_mean)
            # print(cpd_haplotype1_mean)
            # print([(cpd_haplotype1_end - cpd_haplotype1_start) //50000 for cpd_haplotype1_start, cpd_haplotype1_end in zip(cpd_haplotype1_start, cpd_haplotype1_end)])
            # print(cpd_haplotype2_mean)
            # print([(cpd_haplotype2_end - cpd_haplotype2_start) // 50000 for cpd_haplotype2_start, cpd_haplotype2_end in
            #       zip(cpd_haplotype2_start, cpd_haplotype2_end)])

            # print(cpd_haplotype2_mean,cpd_haplotype2_end-cpd_haplotype2_start)

            ##################################
            #phaseblocks flipping main algo
            is_simple_heuristics = True
            #is_phasesets_check_simple_heuristics(ref_start_values_phasesets, ref_end_values_phasesets)
            if len(ref_start_values_phasesets) >= 1:
                is_simple_heuristics = False
            if arguments['enable_simple_heuristics']:
                is_simple_heuristics = True

            if is_simple_heuristics:
                # #plot resultant
                snps_haplotype1_mean, snps_haplotype2_mean, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = \
                    phaseblock_flipping(chrom, arguments, is_simple_heuristics, snps_haplotype1_mean, snps_haplotype2_mean, ref_start_values, ref_start_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
                plot_coverage_data(html_graphs, arguments, chrom, ref_start_values, ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, "phase_correction_1")
            else:
                # #detect centromeres
                #ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = detect_centromeres(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, snps_haplotype1_mean, snps_haplotype2_mean, ref_start_values, arguments['bin_size'])
                #infer missing phaseblocks
                #ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = infer_missing_phaseblocks(ref_start_values, ref_end_values, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, snps_haplotype1_mean, snps_haplotype2_mean, arguments['bin_size'])

                snps_haplotype1_mean, snps_haplotype2_mean, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = \
                    phaseblock_flipping(chrom, arguments, is_simple_heuristics, snps_haplotype1_mean, snps_haplotype2_mean, ref_start_values, ref_start_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

                #switch inter phaseblocks bins
                snps_haplotype1_mean, snps_haplotype2_mean = switch_inter_phaseblocks_bins(chrom, arguments, ref_start_values, snps_haplotype1_mean, snps_haplotype2_mean, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets)

                plot_coverage_data(html_graphs, arguments, chrom, ref_start_values, ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values,
                                   haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, "phase_correction_0")

                ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = contiguous_phaseblocks(haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, loh_region_starts, loh_region_ends)

            #flip phaseblocks based on more contiguous phaseblocks
            #haplotype_1_values_phasesets, haplotype_2_values_phasesets, snps_haplotype1_mean, snps_haplotype2_mean = flip_phaseblocks_contigous(chrom, arguments, ref_start_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, snps_haplotype1_mean, snps_haplotype2_mean)
            plot_coverage_data(html_graphs, arguments, chrom, ref_start_values, ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, "phase_correction_1")

            ##################################
            chr_list = [chrom for ch in range(len(snps_haplotype1_mean))]
            df_updated_coverage.append(pd.DataFrame(list(zip(chr_list, ref_start_values, ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean, unphased_reads_values)),
                                             columns=['chr', 'start', 'end', 'hp1', 'hp2', 'hp3']))

            ##################################
            # start_values_phasesets_contiguous_df_all = remove_overlaping_contiguous(chrom, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2)
            chr_list = [chrom for ch in range(len(ref_start_values_phasesets))]
            start_values_phasesets_contiguous_all.append(pd.DataFrame(list(zip(chr_list, ref_start_values_phasesets)), columns=['chr', 'start']))

    html_graphs.write("</body></html>")

    write_segments_coverage(loh_regions_events_all, arguments['genome_name'] + '_loh_segments.csv', arguments)
    write_df_csv(pd.concat(start_values_phasesets_contiguous_all), arguments['out_dir_plots']+'/data_phasing/'+arguments['genome_name']+'_phasesets.csv')
    write_df_csv(pd.concat(df_updated_coverage), arguments['out_dir_plots']+'/data_phasing/'+arguments['genome_name']+'_coverage.csv')

    csv_df_phase_change_segments = csv_df_chromosomes_sorter(arguments['out_dir_plots']+'/data_phasing/' + arguments['genome_name'] + '_phase_change_segments.csv', ['chr', 'start', 'end'])
    csv_df_phasesets_segments = csv_df_chromosomes_sorter(arguments['out_dir_plots']+'/data_phasing/' + arguments['genome_name'] + '_phasesets.csv', ['chr', 'start'])
    csv_df_loh_regions = csv_df_chromosomes_sorter(arguments['out_dir_plots']+'/data_phasing/' + arguments['genome_name'] + '_loh_segments.csv', ['chr', 'start', 'end'])

    if arguments["normal_phased_vcf"]:
        get_phasingblocks(arguments["normal_phased_vcf"])
        logging.info('VCF edit for phase change segments')
        out_vcf = os.path.join(arguments['out_dir_plots'], 'phasing_output', arguments['genome_name']+'.rephased.vcf.gz')
        rephase_vcf(csv_df_phase_change_segments, csv_df_phasesets_segments, arguments["normal_phased_vcf"], out_vcf)
        index_vcf(out_vcf)
        get_phasingblocks(out_vcf)

    elif not arguments["normal_phased_vcf"] and arguments["tumor_vcf"]:
        get_phasingblocks(arguments["tumor_vcf"])
        logging.info('VCF edit for phase change segments')
        out_vcf = os.path.join(arguments['out_dir_plots'], 'phasing_output', arguments['genome_name']+'.rephased.vcf.gz')
        rephase_vcf(csv_df_phase_change_segments, csv_df_phasesets_segments, arguments["tumor_vcf"], out_vcf)
        index_vcf(out_vcf)
        get_phasingblocks(out_vcf)

    if arguments["rehaplotag_tumor_bam"]:
        logging.info('Rehaplotagging tumor BAM')
        tumor_bam_haplotag(arguments, out_vcf)

    return 0


#Tumor-normal (tumor and normal VCFs)
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam  --tumor-vcf /home/rezkuh/gits/data/HG008_HiFi/HG008_HiFi.vcf.gz  --normal-phased-vcf /home/rezkuh/gits/data/HG008_HiFi/HG008BL_HiFi.vcf.gz --genome-name HG008_HiFi --out-dir-plots HG008_HiFi --cut-threshold 150 --rephase-normal-vcf True
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam  --tumor-vcf /home/rezkuh/gits/data/1395/1395.vcf.gz  --normal-phased-vcf /home/rezkuh/gits/data/1395/1395BL.vcf.gz --genome-name 1395 --out-dir-plots 1395 --cut-threshold 150 --rephase-normal-vcf True

#Tumor-normal (normal VCF)
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam   --normal-phased-vcf /home/rezkuh/gits/data/1437/1437BL.vcf.gz --genome-name 1437 --out-dir-plots 1437 --cut-threshold 150

#Tumor only (tumor VCF)
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam   --tumor-vcf /home/rezkuh/gits/data/colo357_R10/colo357.vcf.gz --genome-name colo357 --out-dir-plots colo357 --cut-threshold 150