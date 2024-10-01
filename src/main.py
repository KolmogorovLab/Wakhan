#!/usr/bin/env python3

import os
import statistics
import sys
import pysam
import shutil
import logging
import argparse
import pandas as pd
import numpy as np

from multiprocessing import Pool
from collections import defaultdict

from hapcorrect.src.main_hapcorrect import main_process

from bam_processing import get_all_reads_parallel, update_coverage_hist, get_segments_coverage, haplotype_update_all_bins_parallel, get_snps_frequencies
from utils import get_chromosomes_bins, write_segments_coverage, write_segments_coverage_dict, csv_df_chromosomes_sorter,\
    seperate_dfs_coverage, flatten_smooth, get_contigs_list, write_copynumber_segments_csv, integer_fractional_cluster_means, \
    adjust_diversified_segments, get_chromosomes_bins_bam, normal_genome_proportion, update_subclonal_means_states, adjust_first_copy_mean, \
    merge_adjacent_regions_cn, parse_sv_vcf, weigted_means_ploidy, average_p_value_genome, check_adjust_last_cn_states
from plots import coverage_plots_chromosomes, copy_number_plots_genome, plots_genome_coverage, copy_number_plots_chromosomes, copy_number_plots_genome_breakpoints_unphased, \
    copy_number_plots_genome_breakpoints_unphased_test, copy_number_plots_genome_breakpoints, copy_number_plots_genome_breakpoints_subclonal, copy_number_plots_genome_details
from vcf_processing import vcf_parse_to_csv_for_het_phased_snps_phasesets
from snps_loh import plot_snps_frequencies, plot_snps_ratios_genome
from phasing_correction import generate_phasesets_bins
from optimization import peak_detection_optimization
from extras import sv_vcf_bps_cn_check

def main():
    # default tunable parameters
    MAX_READ_ERROR = 0.1
    MIN_MAPQ = 10
    MIN_SV_SIZE = 50
    BIN_SIZE = 50000
    BIN_SIZE_SNPS = 50000
    MAX_CUT_THRESHOLD = 100
    MIN_ALIGNED_LENGTH = 5000
    MAX_CUT_THRESHOLD_SNPS_COUNTS = 50
    HETS_RATIO_LOH = 0.3
    HETS_SMOOTH_WINDOW = 45
    HETS_LOH_SEG_SIZE = 2000000


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

    parser.add_argument("--breakpoints", dest="breakpoints",
                        metavar="path", required=True, default=None,
                        help="Path to breakpoints/SVs VCF file")
    parser.add_argument("--cpd-internal-segments", dest="cpd_internal_segments",
                        metavar="path", required=False, default=None,
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
    parser.add_argument("--hets-ratio", "--hets_ratio", dest="hets_ratio",
                        default=HETS_RATIO_LOH, metavar="float", type=float, help="Hetrozygous SNPs ratio threshold for LOH detection [0.3]")
    parser.add_argument("--hets-smooth-window", "--hets_smooth_window", dest="hets_smooth_window",
                        default=HETS_SMOOTH_WINDOW, metavar="int", type=int, help="Hetrozygous SNPs ratio smoothing window size for LOH detection [45]")
    parser.add_argument("--hets-loh-seg-size", "--hets_loh_seg_size", dest="hets_loh_seg_size",
                        default=HETS_LOH_SEG_SIZE, metavar="int", type=int, help="LOH detection minimum segment size where Het SNPs ratio is dropped [2M]")

    parser.add_argument("--cut-threshold", "--cut_threshold", dest="cut_threshold",
                        default=MAX_CUT_THRESHOLD, metavar="int", type=int, help="Maximum cut threshold for coverage (readdepth) [100]")
    parser.add_argument("--cut-threshold-snps-counts", "--cut_threshold_snps_counts", dest="cut_threshold_snps_counts",
                        default=MAX_CUT_THRESHOLD_SNPS_COUNTS, metavar="int", type=int,
                        help="Maximum cut threshold for SNPs counts [50]")

    parser.add_argument("--min-aligned-length", "--min_aligned_length", dest="min_aligned_length",
                        default=MIN_ALIGNED_LENGTH, metavar="int", type=int, help="Minimum aligned reads length [5000]")

    parser.add_argument('--pdf-enable', action="store_true",  dest="pdf_enable", required=False,
                        default=False, help="Enabling PDF output coverage plots")

    parser.add_argument('--unphased-reads-coverage-disable', action="store_true",  dest="unphased_reads_coverage_disable", required=False,
                        default=False, help="Disabling unphased reads coverage output in plots")
    parser.add_argument('--without-phasing', action="store_true", dest="without_phasing", required=False,
                        default=False, help="Enabling coverage and copynumbers without phasing in plots")

    parser.add_argument('--phaseblock-flipping-disable', action="store_true",  dest="phaseblock_flipping_disable", required=False,
                        default=False, help="Disabling phaseblock flipping in coverage plots")
    parser.add_argument('--smoothing-enable', action="store_true",  dest="smoothing_enable", required=False,
                        default=False, help="Enabling smoothing in coverage plots")
    parser.add_argument('--phaseblocks-enable', action="store_true",  dest="phaseblocks_enable", required=False,
                        default=False, help="Enabling phaseblocks display in coverage plots")

    parser.add_argument('--copynumbers-disable', action="store_true", dest="copynumbers_disable", required=False,
                        default=False, help="Disabling copy number in coverage plots")
    parser.add_argument('--copynumbers-subclonal-enable', action="store_true", dest="copynumbers_subclonal_enable", required=False,
                        default=False, help="Enabling subclonal copy number in coverage plots")

    parser.add_argument('--enable-debug', action="store_true", dest="enable_debug", required=False,
                        default=False, help="Enabling debug")

    parser.add_argument('--rephase-normal-vcf', action="store_true", dest="rephase_normal_vcf", required=False,
                        default=False, help="enable rephase normal vcf")
    parser.add_argument('--rephase-tumor-vcf', action="store_true", dest="rephase_tumor_vcf", required=False,
                        default=False, help="enable rephase tumor vcf")
    parser.add_argument('--rehaplotag-tumor-bam', action="store_true", dest="rehaplotag_tumor_bam", required=False,
                        default=False, help="enable rehaplotag tumor bam")

    parser.add_argument('--variable-size-bins', action="store_true", dest="variable_size_bins", required=False,
                        default=False, help="enable variable size bins to use breakpoints")

    parser.add_argument('--enable-simple-heuristics', action="store_true", dest="enable_simple_heuristics", required=False,
                        default=False, help="enable simple heuristics")

    parser.add_argument("--bins-cluster-means", dest="bins_cluster_means",
                        default=None, required=False, type=lambda s: [int(item) for item in s.split(',')],
                        help="bins cluster means")

    parser.add_argument("--tumor-purity", dest="tumor_purity", default=0.0, metavar="float", type=float, help="user input tumor purity")
    parser.add_argument("--tumor-ploidy", dest="tumor_ploidy", default=0.0, metavar="float", type=float, help="user input tumor ploidy")
    parser.add_argument("--confidence-subclonal-score", dest="confidence_subclonal_score", default=0.6, metavar="float", type=float, help="user input p-value to detect if a segment is subclonal/off to integer copynumber")

    parser.add_argument("-t", "--threads", dest="threads",
                        default=1, metavar="int", type=int, help="number of parallel threads [8]")
    parser.add_argument('--dryrun', action="store_true", dest="dryrun", required=False,
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

    logging.basicConfig(level=logging.INFO)

    pbs = parse_sv_vcf(args.breakpoints)

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
    if os.path.exists(args.out_dir_plots+'/data'):
        shutil.rmtree(args.out_dir_plots+'/data')
        os.mkdir(args.out_dir_plots+'/data')
    else:
        os.mkdir(args.out_dir_plots+'/data')

    thread_pool = Pool(args.threads)

    if args.phaseblock_flipping_disable:
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
        coverage_histograms = update_coverage_hist(genome_ids, ref_lengths, segments_by_read, args.min_mapping_quality, args.max_read_error, args)
        del segments_by_read

    if not args.phaseblock_flipping_disable:
       main_process() #hapcorrect
       if args.dryrun:
           if args.without_phasing:
               df = pd.read_csv(args.dryrun_path + args.genome_name + '/coverage_hps.csv', sep='\t', names=['chr', 'start', 'end', 'hp1', 'hp2', 'unphased'])
               df['hp1'] = df['hp1'] + df['hp2'] + df['un']
               df.to_csv(args.dryrun_path + args.genome_name + '/coverage.csv', sep='\t', columns=['chr', 'start', 'end', 'hp1'], index=False, header=False)

               csv_df_coverage = csv_df_chromosomes_sorter(args.dryrun_path + args.genome_name + '/coverage.csv', ['chr', 'start', 'end', 'coverage'])
               csv_df_phasesets = csv_df_chromosomes_sorter(args.dryrun_path + args.genome_name + '/coverage_ps.csv', ['chr', 'start', 'end', 'coverage'])
           else:
               csv_df_phasesets = csv_df_chromosomes_sorter(args.dryrun_path + args.genome_name + '/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
               csv_df_coverage = csv_df_chromosomes_sorter(args.out_dir_plots+'/data_phasing/'+args.genome_name+'_coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
       else:
           csv_df_phasesets = csv_df_chromosomes_sorter(args.out_dir_plots+'/data_phasing/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
           csv_df_coverage = csv_df_chromosomes_sorter(args.out_dir_plots+'/data_phasing/'+args.genome_name+'_coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
    else:
        if args.dryrun and not args.without_phasing:
            csv_df_coverage = csv_df_chromosomes_sorter(args.dryrun_path + args.genome_name + '/coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
            csv_df_phasesets = csv_df_chromosomes_sorter(args.dryrun_path + args.genome_name + '/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
        else:
            logging.info('Computing coverage for bins')
            segments = get_chromosomes_bins_bam(args.target_bam[0], args.bin_size, args)
            segments_coverage = get_segments_coverage(segments, coverage_histograms)
            logging.info('Writing coverage for bins')
            write_segments_coverage_dict(segments_coverage, 'coverage.csv', args)

            logging.info('Parsing phaseblocks information')
            if args.normal_phased_vcf:
                output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(args.normal_phased_vcf, args)
            else:
                output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(args.tumor_vcf, args)
            phasesets_segments = generate_phasesets_bins(args.target_bam[0], output_phasesets_file_path, args.bin_size, args) #TODO update for multiple bam files
            logging.info('Computing coverage for phaseblocks')
            phasesets_coverage = get_segments_coverage(phasesets_segments, coverage_histograms)
            logging.info('Writing coverage for phaseblocks')
            write_segments_coverage_dict(phasesets_coverage, 'coverage_ps.csv', args)
            del coverage_histograms

            logging.info('Loading coverage (bins) and coverage (phaseblocks) files...')
            csv_df_phasesets = csv_df_chromosomes_sorter(args.out_dir_plots+'/data/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
            csv_df_coverage = csv_df_chromosomes_sorter(args.out_dir_plots+'/data/coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

    #TODO add chrX,chrY support later on
    #csv_df_coverage = csv_df_coverage.drop(csv_df_coverage[(csv_df_coverage.chr == "chrX") | (csv_df_coverage.chr == "chrY")].index)
    #csv_df_phasesets = csv_df_phasesets.drop(csv_df_phasesets[(csv_df_phasesets.chr == "chrX") | (csv_df_phasesets.chr == "chrY")].index)

    logging.info('Generating coverage plots chromosomes-wise')
    haplotype_1_values_updated, haplotype_2_values_updated, unphased, csv_df_snps_mean, snps_cpd_means, snps_cpd_points_weights, snps_cpd_means_df = \
        coverage_plots_chromosomes(csv_df_coverage, csv_df_phasesets, args, thread_pool)

    if args.without_phasing == False:
        df_hp1, df_hp2, df_unphased = seperate_dfs_coverage(args, csv_df_snps_mean, csv_df_snps_mean.hp1.tolist(), csv_df_snps_mean.hp2.tolist(), csv_df_snps_mean.hp3.tolist())

    logging.info('Generating optimal clusters plots for bins')
    df_segs_hp1 = snps_cpd_means_df[0]
    df_segs_hp2 = snps_cpd_means_df[1]

    print(snps_cpd_means)
    print(snps_cpd_points_weights)

    centers, subclonals, x_axis, observed_hist, single_copy_cov = peak_detection_optimization(snps_cpd_means, snps_cpd_points_weights)
    #centers, subclonals, x_axis, observed_hist = peak_detection_optimization(csv_df_snps_mean.hp1.tolist()+csv_df_snps_mean.hp2.tolist(), [1 for i in range(2*len(csv_df_snps_mean.hp2.tolist()))])
    tumor_cov = statistics.mean([sum(x) for x in zip(haplotype_1_values_updated, haplotype_2_values_updated, unphased)])
    # if args.tumor_purity and args.tumor_ploidy:
    #     #print(normal_genome_proportion(0.45, 8, 100))
    #     #tumor_cov = 100
    #     normal_cov = 60
    #     tumor_cov = statistics.mean([sum(x) for x in zip(haplotype_1_values_updated, haplotype_2_values_updated, unphased)])
    #     purity = (tumor_cov / args.tumor_ploidy) / ((normal_cov / 2) + (tumor_cov / args.tumor_ploidy))
    #     print("tumor_coverage:", tumor_cov)
    #     _, _, _, normal_fraction = normal_genome_proportion(args.tumor_purity, args.tumor_ploidy, tumor_cov)
    #     print("purity, normal fraction:", purity, normal_fraction)
    #     centers = [normal_fraction] + [normal_fraction + (i * centers[1]) for i in range(1, len(centers))]

    if args.without_phasing:
        logging.info('Generating coverage/copy numbers plots genome wide')
        integer_fractional_means = sorted([i for i in range(0, len(centers))])
        df_segs_hp1, df_segs_hp2 = adjust_diversified_segments(centers, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)

        write_copynumber_segments_csv(df_segs_hp1, args, centers, integer_fractional_means, None)
        copy_number_plots_genome(centers, integer_fractional_means, df_hp1, df_segs_hp1, df_hp2, df_segs_hp2, df_hp1, args, x_axis, observed_hist)
        copy_number_plots_genome_breakpoints_unphased_test(centers, integer_fractional_means, df_hp1, df_segs_hp1, df_hp2, df_segs_hp2, df_hp1, args)
    else:
        logging.info('Generating coverage/copy numbers plots genome wide')
        ################################
        data = []
        average_p_value = []
        for j in 0, single_copy_cov -1, single_copy_cov, single_copy_cov + 1, single_copy_cov + 2, single_copy_cov + 3:
            cen_out = [j] + [j + (i * single_copy_cov) for i in range(1, len(centers))]
            df_segs_hp1_updated, df_segs_hp2_updated = adjust_diversified_segments(cen_out, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)
            average_p_value.append(average_p_value_genome(args, cen_out, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2))

            overall_ploidy = weigted_means_ploidy(df_segs_hp1_updated, df_segs_hp2_updated, cen_out, sorted([i for i in range(0, len(cen_out))]))
            tumor_purity = (tumor_cov / overall_ploidy) / (((j*2) / 2) + (tumor_cov / overall_ploidy))
            data.append([overall_ploidy, tumor_purity, cen_out])
            print("overall_ploidy: ", overall_ploidy, "tumor_purity:", tumor_purity, "average_p_value:", average_p_value, "for i:", j, "centers: ", cen_out[0:4])

        for j in 0, np.argmax(average_p_value[1:]) + 1:
            if j == 0:
                cen_out = data[0][2]
                integer_fractional_means = sorted([i for i in range(0, len(cen_out))])
                args.tumor_ploidy = round(data[0][0], 2)
                args.tumor_purity = round(data[0][1], 2)
            else:
                cen_out = data[np.argmax(average_p_value[1:])+1][2]
                integer_fractional_means = sorted([i for i in range(0, len(cen_out))])
                args.tumor_ploidy = round(data[np.argmax(average_p_value[1:])+1][0], 2)
                args.tumor_purity = round(data[np.argmax(average_p_value[1:])+1][1], 2)

            df_segs_hp1_updated, df_segs_hp2_updated = adjust_diversified_segments(cen_out, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)
            write_copynumber_segments_csv(df_segs_hp1_updated, args, cen_out, integer_fractional_means, 1, '_'+ str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_copynumbers_segments.bed')
            write_copynumber_segments_csv(df_segs_hp2_updated, args, cen_out, integer_fractional_means, 2, '_'+ str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_copynumbers_segments.bed')

            copy_number_plots_genome(cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_unphased, args, x_axis, observed_hist)
            copy_number_plots_genome_breakpoints(cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args)

            if args.copynumbers_subclonal_enable:
                df_segs_hp1_updated, df_segs_hp2_updated = adjust_diversified_segments(cen_out, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)
                df_segs_hp1_updated, df_segs_hp2_updated = update_subclonal_means_states(cen_out, subclonals, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2, args)
                df_segs_hp1_updated = merge_adjacent_regions_cn(df_segs_hp1_updated, args)
                df_segs_hp2_updated = merge_adjacent_regions_cn(df_segs_hp2_updated, args)
                cen_out = [int(i) for i in cen_out]
                copy_number_plots_genome_breakpoints_subclonal(cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args)
                write_copynumber_segments_csv(df_segs_hp1_updated, args, cen_out, integer_fractional_means, 1, '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_copynumbers_subclonal_segments.bed')
                write_copynumber_segments_csv(df_segs_hp2_updated, args, cen_out, integer_fractional_means, 2, '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_copynumbers_subclonal_segments.bed')

        # plot_snps_frequencies(args, csv_df_snps_mean, df_segs_hp1_updated, df_segs_hp2_updated, centers, integer_fractional_means, thread_pool)
        ################################

        ################################
        # for j in 0, single_copy_cov -1, single_copy_cov, single_copy_cov + 1, single_copy_cov + 2:
        #     cen_out = [j] + [j + (i * single_copy_cov) for i in range(1, len(centers))]
        #     df_segs_hp1_updated, df_segs_hp2_updated = adjust_diversified_segments(cen_out, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)
        #     average_p_value = average_p_value_genome(args, cen_out, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2)
        #
        #     overall_ploidy = weigted_means_ploidy(df_segs_hp1_updated, df_segs_hp2_updated, cen_out, sorted([i for i in range(0, len(cen_out))]))
        #     tumor_purity = (tumor_cov / overall_ploidy) / (((j*2) / 2) + (tumor_cov / overall_ploidy))
        #     print("overall_ploidy: ", overall_ploidy, "tumor_purity:", tumor_purity, "average_p_value:", average_p_value, "for i:", j, "centers: ", cen_out[0:4])
        ################################

        # df_segs_hp1_updated, df_segs_hp2_updated = adjust_diversified_segments(centers, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)
        # centers_ = adjust_first_copy_mean(args, centers, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2)
        # integer_fractional_means = sorted([i for i in range(0, len(centers))])
        # df_segs_hp1_updated, df_segs_hp2_updated = adjust_diversified_segments(centers, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)
        #
        # tumor_cov = statistics.mean([sum(x) for x in zip(haplotype_1_values_updated, haplotype_2_values_updated, unphased)])
        # average_p_value = average_p_value_genome(args, centers, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2)
        # print("average_p_value:", average_p_value)
        # overall_ploidy = weigted_means_ploidy(df_segs_hp1_updated, df_segs_hp2_updated, centers, integer_fractional_means)
        # print("overall_ploidy: ", overall_ploidy)
        # #tumor_purity = (tumor_cov / overall_ploidy) / ((optimized_normal_coverage / 2) + (tumor_cov / overall_ploidy))
        # #print("overall_ploidy: ", overall_ploidy, "tumor_purity:", tumor_purity)
        #
        # print(centers)
        # print(integer_fractional_means)
        # write_copynumber_segments_csv(df_segs_hp1_updated, args, centers, integer_fractional_means, 1, '_copynumbers_segments.bed')
        # write_copynumber_segments_csv(df_segs_hp2_updated, args, centers, integer_fractional_means, 2, '_copynumbers_segments.bed')
        #
        # copy_number_plots_genome(centers, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_unphased, args, x_axis, observed_hist)
        # copy_number_plots_genome_breakpoints(centers, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args)
        # plot_snps_frequencies(args, csv_df_snps_mean, df_segs_hp1_updated, df_segs_hp2_updated, centers, integer_fractional_means, thread_pool)

        # if args.copynumbers_subclonal_enable:
        #     df_segs_hp1_updated, df_segs_hp2_updated = adjust_diversified_segments(centers, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)
        #     df_segs_hp1_updated, df_segs_hp2_updated = update_subclonal_means_states(centers, subclonals, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2, args)
        #     df_segs_hp1_updated = merge_adjacent_regions_cn(df_segs_hp1_updated, args)
        #     df_segs_hp2_updated = merge_adjacent_regions_cn(df_segs_hp2_updated, args)
        #     centers  =  [int(i) for i in centers]
        #     copy_number_plots_genome_breakpoints_subclonal(centers, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args)
        #     write_copynumber_segments_csv(df_segs_hp1_updated, args, centers, integer_fractional_means, 1, '_'+ str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) + '_copynumbers_subclonal_segments.bed')
        #     write_copynumber_segments_csv(df_segs_hp2_updated, args, centers, integer_fractional_means, 2, '_'+ str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_copynumbers_subclonal_segments.bed')

    #SNPs LOH and plots
    if args.tumor_vcf:
        plot_snps_ratios_genome(args)

    if os.path.exists(args.out_dir_plots+'/data'): #
        shutil.rmtree(args.out_dir_plots+'/data')
    if os.path.exists(args.out_dir_plots+'/data_phasing'):
        shutil.rmtree(args.out_dir_plots+'/data_phasing')

    return 0

if __name__ == "__main__":
    main()

#UCSC tumor/normal celllines
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots 1954  --normal-phased-vcf /home/rezkuh/gits/data/1954/1954BL.vcf.gz --copynumbers-enable True  --unphased-reads-coverage-enable True  --phaseblock-flipping-enable True --phaseblocks-enable True   --genome-name 1954  --cut-threshold 150

#pancreatic_organoid data
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots  --phased-vcf /home/rezkuh/gits/data/pancreatic_organoid/pancBL.vcf.gz  --copynumbers-enable True  --unphased-reads-coverage-enable True --snps-freq-vcf-enable True --phaseblock-flipping-enable True --phaseblocks-enable True  --genome-name pancreatic_organoid  --cut-threshold 150

#Tumor only (HPV)
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots --snps-freq-vcf-enable True --cut-threshold 50 --copynumbers-enable True --phaseblock-flipping-enable True   --snps-freq-vcf-enable True  --phased-vcf /home/rezkuh/gits/data/R10/HT3/HT3.vcf.gz  --genome-name R10/HT3
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots --snps-freq-vcf-enable True --cut-threshold 50 --copynumbers-enable True --phaseblock-flipping-enable True   --snps-freq-vcf-enable True  --phased-vcf /home/rezkuh/gits/data/R10/CaSki/CaSki.vcf.gz  --genome-name R10/CaSki

#NIST GIAB HG008
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots HG008_HiFi  --normal-phased-vcf /home/rezkuh/gits/data/HG008_HiFi/HG008BL_HiFi.vcf.gz --tumor-vcf /home/rezkuh/gits/data/HG008_HiFi/HG008_HiFi.vcf.gz --copynumbers-enable True  --unphased-reads-coverage-enable True  --phaseblock-flipping-enable True --phaseblocks-enable True   --genome-name HG008_HiFi  --cut-threshold 150

#colo829
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots colo829  --normal-phased-vcf /home/rezkuh/gits/data/colo829/colo829_pepper_normal.phased.vcf.gz  --copynumbers-enable True   --phaseblock-flipping-enable True --phaseblocks-enable True   --genome-name colo829  --cut-threshold 150

#colo829-porec
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots colo829-porec  --normal-phased-vcf /home/rezkuh/gits/data/colo829-porec/colo829.vcf.gz  --copynumbers-enable True    --phaseblocks-enable True   --genome-name colo829-porec  --cut-threshold 150

#Mouse unphased data
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam    --copynumbers-enable True    --unphased-reads-coverage-enable True    --cut-threshold 75  --without-phasing True --tumor-vcf /home/rezkuh/gits/data/mouse/somatic_calls/C1_somatic_calls_pass_snp.vcf.gz --out-dir-plots C1 --contigs 1-19 --genome-name C1 --bin-size-snps 1000000
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam    --copynumbers-enable True    --unphased-reads-coverage-enable True    --cut-threshold 75  --without-phasing True --tumor-vcf /home/rezkuh/gits/data/mouse/somatic_calls/C23_somatic_calls_pass_snp.vcf.gz --out-dir-plots C23 --contigs 1-19,X,Y --genome-name C23 --bin-size-snps 1000000

#Dog data
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam   --normal-phased-vcf /home/rezkuh/gits/data/OT4/ON2.vcf.gz --tumor-vcf /home/rezkuh/gits/data/OT2/OT2.vcf.gz     --copynumbers-enable True    --genome-name OT2 --out-dir-plots OT2  --cut-threshold 60 --phaseblock-flipping-enable True --phaseblocks-enable True --contigs chr1-38

#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots 1437  --normal-phased-vcf /home/rezkuh/gits/data/1437_re/1437BL.vcf.gz --copynumbers-enable True  --unphased-reads-coverage-enable True  --phaseblocks-enable True   --genome-name 1437  --cut-threshold 150

#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots 1937  --normal-phased-vcf /home/rezkuh/gits/data/1937_re/1937BL.vcf.gz --copynumbers-enable True  --unphased-reads-coverage-enable True  --phaseblocks-enable True   --genome-name 1937 --contigs chr1-22,X --cut-threshold 150 --breakpoints /home/rezkuh/gits/data/1937_re/ --cpd-internal-segments True
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots colo357_R10  --tumor-vcf /home/rezkuh/gits/data/colo357_R10/colo357.vcf.gz --copynumbers-enable True  --unphased-reads-coverage-enable True  --phaseblocks-enable True   --genome-name colo357_R10 --contigs chr1-22,X --cut-threshold 150 --breakpoints /home/rezkuh/gits/data/colo357_R10/ --cpd-internal-segments True
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots R10/colo357_R10  --tumor-vcf /home/rezkuh/gits/data/R10/colo357_R10/colo357.vcf.gz --copynumbers-enable True  --unphased-reads-coverage-enable True  --phaseblocks-enable True   --genome-name R10/colo357_R10 --contigs chr1-22,X --cut-threshold 150 --breakpoints /home/rezkuh/gits/data/R10/colo357_R10/

#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots colo829 --genome-name colo829  --cut-threshold 200   --normal-phased-vcf /home/rezkuh/gits/data/colo829/colo829BL.vcf.gz --breakpoints /home/rezkuh/gits/data/colo829/severus_somatic.vcf
#--dryrun True --dryrun-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots colo357_R10  --copynumbers-enable True  --unphased-reads-coverage-enable True    --genome-name colo357_R10  --cut-threshold 150  --tumor-vcf /home/rezkuh/gits/data/R10/colo357_R10/colo357.vcf.gz  --phaseblock-flipping-enable True --breakpoints /home/rezkuh/gits/data/R10/colo357_R10/
