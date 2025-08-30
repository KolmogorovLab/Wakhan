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
from scipy.signal import find_peaks

from src.hapcorrect.src.main_hapcorrect import main_process
from src.__version__ import __version__

from src.bam_processing import get_all_reads_parallel, update_coverage_hist, get_segments_coverage, haplotype_update_all_bins_parallel, get_snps_frequencies
from src.utils import get_chromosomes_bins, write_segments_coverage, write_segments_coverage_dict, csv_df_chromosomes_sorter,\
    seperate_dfs_coverage, flatten_smooth, get_contigs_list, write_copynumber_segments_csv, integer_fractional_cluster_means, \
    adjust_diversified_segments, get_chromosomes_bins_bam, normal_genome_proportion, update_subclonal_means_states, adjust_first_copy_mean, \
    merge_adjacent_regions_cn, merge_adjacent_regions_cn_unphased, parse_sv_vcf, weigted_means_ploidy, average_p_value_genome, collect_loh_centromere_regions, centromere_regions_blacklist, \
    extract_centromere_regions, update_genes_phase_corrected_coverage, weighted_means, extract_breakpoints_additional, write_df_csv, adjust_bps_cn_segments_boundries, dna_purity_to_cell_purity, move_100pct_purity_sol, find_p_values_peaks, add_confidence_score_cn_segemnts, \
    centromere_regions_blacklist_bins
from src.plots import coverage_plots_chromosomes, copy_number_plots_genome_details, copy_number_plots_genome, plots_genome_coverage, copy_number_plots_chromosomes, breakpoints_segments_means, \
    copy_number_plots_genome_breakpoints, copy_number_plots_genome_breakpoints_subclonal, copy_number_plots_genome_subclonal, genes_copy_number_plots_genome, genes_plots_genome, heatmap_copy_number_plots_genome, plot_ploidy_purity_p_values
from src.vcf_processing import vcf_parse_to_csv_for_het_phased_snps_phasesets
from src.snps_loh import plot_snps_frequencies_without_phasing, plot_snps_frequencies, plot_snps_ratios_genome, snps_df_loh, variation_plots, write_loh_regions
from src.phasing_correction import generate_phasesets_bins, fix_inter_cn_phase_switch_errors, bins_correction_phaseblocks
from src.optimization import peak_detection_optimization
from src.generate_vcf import read_cn_segments_process_vcf
from src.extras import sv_vcf_bps_cn_check

logger = logging.getLogger()

solutions_df = pd.DataFrame(columns=['repository_name', 'dna_purity', 'cell_purity', 'ploidy', 'confidence'])

def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: " "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: " "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.INFO)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)

def _version():
    return __version__

def find_peaks_with_min_distance(signal, min_distance=2):
    peaks, _ = find_peaks(signal, distance=min_distance)
    return peaks.tolist()

def copy_numbers_assignment_haplotypes(args, tumor_cov, max_limit, single_copy_cov, centers, subclonals, df_hp1, df_hp2, df_segs_hp1, df_segs_hp2, snps_cpd_means_df, csv_df_snps_mean, df_snps_in_csv, df_unphased, x_axis, observed_hist, is_half):
    data = []
    average_p_value = []
    #max_limit = 3/2 #debug
    for normal_coverage in np.arange(0, max_limit, 0.1):
        cen_out = [normal_coverage] + [normal_coverage + (i * single_copy_cov) for i in range(1, len(centers))]
        df_segs_hp1_updated, df_segs_hp2_updated = adjust_diversified_segments(cen_out, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)
        p_value = average_p_value_genome(args, cen_out, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2)

        overall_ploidy = weigted_means_ploidy(args, df_segs_hp1_updated, df_segs_hp2_updated, cen_out, sorted([i for i in range(0, len(cen_out))]))

        if overall_ploidy < 0.1:
            continue
        normal_coverage = normal_coverage * 2 #make it diploid
        dna_tumor_purity = (tumor_cov - normal_coverage) / tumor_cov
        cellular_tumor_purity = ((tumor_cov - normal_coverage) / overall_ploidy) / ((normal_coverage / 2) + ((tumor_cov - normal_coverage) / overall_ploidy))

        if args.dna_purity:
            tumor_purity = dna_tumor_purity
        else:  # cell purity
            tumor_purity = cellular_tumor_purity
        # tumor_purity = dna_purity_to_cell_purity(tumor_purity, overall_ploidy)
        _, _, _, normal_fraction = normal_genome_proportion(tumor_purity, overall_ploidy, tumor_cov)
        if normal_coverage == 0 and p_value > 0:
            average_p_value.append(p_value)
            data.append([overall_ploidy, tumor_purity, cen_out, p_value, dna_tumor_purity])
            logger.info("overall_ploidy: %s, dna_tumor_purity: %s, cell_tumor_purity: %s, average_p_value: %s, for i: %s,  centers: %s, norm frac: %s",
                overall_ploidy, dna_tumor_purity, cellular_tumor_purity, p_value, normal_coverage, cen_out[0:4], normal_fraction)
            continue

        if (float(args.purity_range.split('-')[0]) <= tumor_purity <= float(args.purity_range.split('-')[1])) and (float(args.ploidy_range.split('-')[0]) <= overall_ploidy <= float(args.ploidy_range.split('-')[1])):
            average_p_value.append(p_value)
            data.append([overall_ploidy, tumor_purity, cen_out, p_value, dna_tumor_purity])
            logger.info("overall_ploidy: %s, dna_tumor_purity: %s, cell_tumor_purity: %s, average_p_value: %s, for i: %s,  centers: %s, norm frac: %s",
                overall_ploidy, dna_tumor_purity, cellular_tumor_purity, p_value, normal_coverage, cen_out[0:4], normal_fraction)

    plot_ploidy_purity_p_values(args, [data[n][0] for n in range(len(data))], [data[n][1] for n in range(len(data))], [data[n][3] for n in range(len(data))])
    if average_p_value:
        optimized_normal = find_p_values_peaks(average_p_value)
        for j in optimized_normal:
            args.tumor_ploidy = round(data[j][0], 2)
            args.tumor_purity = round(data[j][1], 2)
            dna_tumor_fraction = round(data[j][4], 2)
            cen_out = data[j][2]
            logger.info('Normal optimized clusters means: %s', cen_out)
            integer_fractional_means = sorted([i for i in range(0, len(cen_out))])
            p_value_confidence = round(data[j][3], 2)
            logger.info('Generating coverage/copy numbers plots genome wide for solution with tumor cellular fraction {0}, ploidy {1} and tumor dna fraction {2}'.format(args.tumor_purity, args.tumor_ploidy, dna_tumor_fraction))
            solutions_df.loc[j] = [str(args.tumor_ploidy)+'_'+str(args.tumor_purity)+'_'+str(p_value_confidence), dna_tumor_fraction, args.tumor_purity, args.tumor_ploidy, p_value_confidence]

            #if args.breakpoints:
            #    df_segs_hp1, df_segs_hp2 = check_dv_support_in_breakpoints(cen_out[1], df_hp1, df_hp2, df_segs_hp1, df_segs_hp2)
            df_segs_hp1_updated, df_segs_hp2_updated = adjust_diversified_segments(cen_out, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)
            loh_regions = collect_loh_centromere_regions(df_segs_hp1_updated, df_segs_hp2_updated, cen_out, integer_fractional_means, args)

            #df_segs_hp1_updated = merge_adjacent_regions_cn(df_segs_hp1_updated, args)
            #df_segs_hp2_updated = merge_adjacent_regions_cn(df_segs_hp2_updated, args)

            #df_segs_hp1_updated = merge_bps_regions_cn(df_segs_hp1_updated, args)

            write_copynumber_segments_csv(df_hp1, df_hp2, df_segs_hp1_updated, args, cen_out, integer_fractional_means, 1, '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value_confidence) + '_copynumbers_segments_HP_1.bed', p_value_confidence, is_half)
            write_copynumber_segments_csv(df_hp1, df_hp2, df_segs_hp2_updated, args, cen_out, integer_fractional_means, 2, '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value_confidence) + '_copynumbers_segments_HP_2.bed', p_value_confidence, is_half)
            write_loh_regions(loh_regions, 'loh_regions.bed', args, p_value_confidence, is_half)
            read_cn_segments_process_vcf(args, str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value_confidence), 'integers')

            df_segs_hp1_updated, df_segs_hp2_updated = centromere_regions_blacklist(args, df_segs_hp1_updated, df_segs_hp2_updated)
            df_hp1, df_hp2 = centromere_regions_blacklist_bins(args, df_hp1, df_hp2, df_segs_hp1_updated, df_segs_hp2_updated)

            copy_number_plots_genome_details(cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_unphased, args, x_axis, observed_hist, p_value_confidence, is_half)

            # df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2 = fix_inter_cn_phase_switch_errors(args, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2)

            df_segs_hp1_updated, df_segs_hp2_updated = add_confidence_score_cn_segemnts(centers, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2, args)

            if args.breakpoints:
                copy_number_plots_genome_breakpoints(cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args, p_value_confidence, loh_regions, df_snps_in_csv, is_half)
            else:
                copy_number_plots_genome(cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args, p_value_confidence, loh_regions, df_snps_in_csv, is_half)
            df_segs_hp1_updated = df_segs_hp1_updated.drop('confidence_value', axis=1)
            df_segs_hp2_updated = df_segs_hp2_updated.drop('confidence_value', axis=1)
            variation_plots(args, csv_df_snps_mean, df_segs_hp1_updated, df_segs_hp2_updated, cen_out, integer_fractional_means, df_snps_in_csv, loh_regions, p_value_confidence, is_half)
            #debug uncomment
            if not args.phaseblock_flipping_disable:
               df_genes = update_genes_phase_corrected_coverage(args, df_segs_hp1_updated, df_segs_hp2_updated, p_value_confidence, cen_out, integer_fractional_means, is_half)
               genes_plots_genome(df_genes, cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args, p_value_confidence, loh_regions, df_snps_in_csv, is_half)

                # genes_copy_number_plots_genome(df_genes, cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args, p_value_confidence, loh_regions, df_snps_in_csv)
                #heatmap_copy_number_plots_genome(df_genes, cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args, p_value_confidence, loh_regions, df_snps_in_csv, is_half)

            if args.copynumbers_subclonal_enable:
                df_segs_hp1_updated, df_segs_hp2_updated = adjust_diversified_segments(cen_out, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)
                df_segs_hp1_updated, df_segs_hp2_updated = update_subclonal_means_states(cen_out, subclonals, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2, args)
                loh_regions = collect_loh_centromere_regions(df_segs_hp1_updated, df_segs_hp2_updated, cen_out, integer_fractional_means, args)
                write_loh_regions(loh_regions, 'loh_regions_subclonal.bed', args, p_value_confidence, is_half)
                # df_segs_hp1_updated = merge_adjacent_regions_cn(df_segs_hp1_updated, args)
                # df_segs_hp2_updated = merge_adjacent_regions_cn(df_segs_hp2_updated, args)

                cen_out = [int(i) for i in cen_out]
                if args.breakpoints:
                    copy_number_plots_genome_breakpoints_subclonal(cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args, p_value_confidence, loh_regions, df_snps_in_csv, is_half)
                else:
                    copy_number_plots_genome_subclonal(cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args, p_value_confidence, loh_regions, df_snps_in_csv, is_half)
                write_copynumber_segments_csv(df_hp1, df_hp2, df_segs_hp1_updated, args, cen_out, integer_fractional_means, 1, '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value_confidence) + '_copynumbers_subclonal_segments_HP_1.bed', p_value_confidence, is_half)
                write_copynumber_segments_csv(df_hp1, df_hp2, df_segs_hp2_updated, args, cen_out, integer_fractional_means, 2, '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value_confidence) + '_copynumbers_subclonal_segments_HP_2.bed', p_value_confidence, is_half)
                read_cn_segments_process_vcf(args, str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value_confidence), 'subclonal')
    else:
        logger.info(
            'No estimated purity [%s] and ploidy [%s] value detected inside given ranges or overall ploidy is less than 0.1', args.purity_range, args.ploidy_range)

def parse_aruments_(args):
    parser = argparse.ArgumentParser(description="CNA/BNA tool with default arguments")

    subparsers = parser.add_subparsers(dest="command")  # store which subcommand is used

    # Default arguments (when no command given)
    parser.add_argument("--input", required=True, help="Default input file")
    parser.add_argument("--output", help="Default output file")

    # CNA subcommand
    cna_parser = subparsers.add_parser("cna", help="Run CNA analysis")
    cna_parser.add_argument("--cna-file", required=False, help="CNA input file")
    cna_parser.add_argument("--input", required=True, help="CNA input file")
    cna_parser.add_argument("--cna-param", type=int, default=10, help="CNA parameter")

    # BNA subcommand
    bna_parser = subparsers.add_parser("bna", help="Run BNA analysis")
    bna_parser.add_argument("--bna-file", required=True, help="BNA input file")
    bna_parser.add_argument("--bna-threshold", type=float, default=0.05, help="BNA threshold")

    args = parser.parse_args()

    # If no subcommand, use default group
    if args.command is None:
        print("Running default mode...")
        print(f"Input: {args.input}")
        print(f"Output: {args.output}")
    elif args.command == "cna":
        print("Running CNA mode...")
        print(f"CNA file: {args.cna_file}")
        print(f"CNA param: {args.cna_param}")
    elif args.command == "bna":
        print("Running BNA mode...")
        print(f"BNA file: {args.bna_file}")
        print(f"BNA threshold: {args.bna_threshold}")
def build_parser():
    # default tunable parameters
    MAX_READ_ERROR = 0.1
    MIN_MAPQ = 10
    BIN_SIZE = 50000
    BIN_SIZE_SNPS = 200000
    MAX_CUT_THRESHOLD = 100
    MIN_ALIGNED_LENGTH = 5000
    MAX_CUT_THRESHOLD_SNPS_COUNTS = 50
    HETS_RATIO_LOH = 0.15
    HETS_SMOOTH_WINDOW = 1
    HETS_LOH_SEG_SIZE = 200000 #2000000
    BP_MIN_LENGTH = 10000

    DEFAULT_CONTIGS = 'chr1-22,chrX' #('chr1-22' '1-22') ('chr1-22,chrX' '1-22,X')
    DEFAULT_PURITY = '0.25-1.0'
    DEFAULT_PLOIDY = '1.0-5.5'
    DEFAULT_REF = 'grch38'

    global_parser = argparse.ArgumentParser(prog="wakhan", description="Wakhan plots coverage and copy number profiles from a bam and phased VCF files")
    ###################################################################################
    global_parser.add_argument("command", nargs="?", choices=["cna", "hapcorrect"], help="Optional cna or hapcorrect modes. If omitted, runs default mode with both modes." )
    ###################################################################################
    global_parser.add_argument("--target-bam", dest="target_bam", required=True, metavar="path", default=None, nargs="+", help="path to one or 1/4 multiple target haplotagged bam files (must be indexed)")
    global_parser.add_argument("--control-bam", dest="control_bam", metavar="path", default=None, nargs="+", help="path to one or multiple control haplotagged bam files (must be indexed)")
    global_parser.add_argument("--reference", dest="reference", required=True, metavar="path", default=None, help="path to reference")

    global_parser.add_argument("--out-dir-plots", dest="out_dir_plots", required=True, default=None, metavar="path", help="Output directory")

    global_parser.add_argument("--normal-phased-vcf", dest="normal_phased_vcf", metavar="path", default=None, help="Path to normal phased vcf")
    global_parser.add_argument("--tumor-phased-vcf", dest="tumor_phased_vcf", metavar="path", default=None, help="Path to tumor phased VCF for LOH detection")

    global_parser.add_argument("--centromere", dest="centromere", metavar="path", default='annotations/grch38.cen_coord.curated.bed', help="Path to centromere annotations BED file")
    global_parser.add_argument("--cancer-genes", dest="cancer_genes", metavar="path", default='annotations/COSMIC_cancer_genes.tsv', help="Path to default COSMIC Cancer Genes TSV file")
    global_parser.add_argument("--user-input-genes", dest="user_input_genes", metavar="path", help="Path to user input genes names in *.bed file [each name in a single line (annotations/user_input_genes_example_1.bed), or tab delimated genes etries (annotations/user_input_genes_example_2.bed)], these genes will be used in plots instead of default COSMIC cancer genes")
    global_parser.add_argument("--reference-name", dest="reference_name", default=DEFAULT_REF, help="Default reference name: grch38 [grch38, chm13]")
    global_parser.add_argument("--histogram-coverage", dest="histogram_coverage", action="store_true", help="use histogram coverage instead of SNPs pileup")

    global_parser.add_argument("--genome-name", dest="genome_name", default='Sample', help="Genome sample/cell line name to be displayed on plots")
    global_parser.add_argument("--contigs", dest="contigs", default=DEFAULT_CONTIGS, help="List of contigs (choromosomes) to be included in the plots, default chr1-22,chrX [e.g., chr1-22,X,Y], Note: Please use 1-22,X [e.g., 1-22,X,Y] in case REF, BAM, and VCFs entries don't contain `chr` name/notion")

    global_parser.add_argument("--bin-size", "--bin_size", dest="bin_size", default=BIN_SIZE, metavar="int", type=int, help="coverage (readdepth) bin size [50k]")
    global_parser.add_argument("--bin-size-snps", "--bin_size_snps", dest="bin_size_snps", default=BIN_SIZE_SNPS, metavar="int", type=int, help="SNPs bin size [50k]")
    global_parser.add_argument("--hets-ratio", "--hets_ratio", dest="hets_ratio", default=HETS_RATIO_LOH, metavar="float", type=float, help="Hetrozygous SNPs ratio threshold for LOH detection [0.25]")
    global_parser.add_argument("--hets-smooth-window", "--hets_smooth_window", dest="hets_smooth_window", default=HETS_SMOOTH_WINDOW, metavar="int", type=int, help="Hetrozygous SNPs ratio smoothing window size for LOH detection [45]")
    global_parser.add_argument("--hets-loh-seg-size", "--hets_loh_seg_size", dest="hets_loh_seg_size", default=HETS_LOH_SEG_SIZE, metavar="int", type=int, help="LOH detection minimum segment size where Het SNPs ratio is dropped [2M]")
    global_parser.add_argument('--loh-enable', action="store_true",  dest="loh_enable", default=True, help="Enabling LOH regions in CN plots")

    global_parser.add_argument("--cut-threshold", "--cut_threshold", dest="cut_threshold", default=MAX_CUT_THRESHOLD, metavar="int", type=int, help="Maximum cut threshold for coverage (readdepth) [100]")
    global_parser.add_argument("--cut-threshold-snps-counts", "--cut_threshold_snps_counts", dest="cut_threshold_snps_counts", default=MAX_CUT_THRESHOLD_SNPS_COUNTS, metavar="int", type=int, help="Maximum cut threshold for SNPs counts [50]")

    global_parser.add_argument("--min-aligned-length", "--min_aligned_length", dest="min_aligned_length", default=MIN_ALIGNED_LENGTH, metavar="int", type=int, help="Minimum aligned reads length [5000]")

    global_parser.add_argument('--pdf-enable', action="store_true",  dest="pdf_enable", default=False, help="Enabling PDF output coverage plots")

    global_parser.add_argument('--unphased-reads-coverage-disable', action="store_true",  dest="unphased_reads_coverage_disable", default=False, help="Disabling unphased reads coverage output in plots")
    global_parser.add_argument('--without-phasing', action="store_true", dest="without_phasing", default=False, help="Enabling coverage and copynumbers without phasing in plots")

    global_parser.add_argument('--phaseblock-flipping-disable', action="store_true",  dest="phaseblock_flipping_disable", default=False, help="Disabling phaseblock flipping in coverage plots")
    global_parser.add_argument('--smoothing-enable', action="store_true",  dest="smoothing_enable", default=False, help="Enabling smoothing in coverage plots")
    global_parser.add_argument('--phaseblocks-enable', action="store_true",  dest="phaseblocks_enable", default=False, help="Enabling phaseblocks display in coverage plots")

    global_parser.add_argument('--enable-debug', action="store_true", dest="enable_debug", default=False, help="Enabling debug")

    global_parser.add_argument('--rephase-normal-vcf', action="store_true", dest="rephase_normal_vcf", default=False, help="enable rephase normal vcf")
    global_parser.add_argument('--rephase-tumor-vcf', action="store_true", dest="rephase_tumor_vcf", default=False, help="enable rephase tumor vcf")
    global_parser.add_argument('--rehaplotag-tumor-bam', action="store_true", dest="rehaplotag_tumor_bam", default=False, help="enable rehaplotag tumor bam")

    global_parser.add_argument('--variable-size-bins', action="store_true", dest="variable_size_bins", default=False, help="enable variable size bins to use breakpoints")

    global_parser.add_argument('--enable-simple-heuristics', action="store_true", dest="enable_simple_heuristics", default=False, help="enable simple heuristics")

    global_parser.add_argument("-t", "--threads", dest="threads", default=1, metavar="int", type=int, help="number of parallel threads [8]")
    global_parser.add_argument('--quick-start', action="store_true", dest="quick_start", default=False, help="Enabling quick_start")
    global_parser.add_argument("--quick-start-coverage-path", dest="quick_start_coverage_path", default=None, metavar="path", help="quick start coverage data directory")

    global_parser.add_argument("--max-read-error", dest="max_read_error", default=MAX_READ_ERROR, metavar="float", type=float, help=f"maximum base alignment error [{MAX_READ_ERROR}]")
    global_parser.add_argument("--min-mapq", dest="min_mapping_quality", default=MIN_MAPQ, metavar="int", type=int, help=f"minimum mapping quality for aligned segment [{MIN_MAPQ}]")

    global_parser.add_argument("--cpd-internal-segments", dest="cpd_internal_segments", action="store_true", help="change point detection algo for more precise segments after breakpoint/cpd segments")
    global_parser.add_argument("--change-point-detection-for-cna", dest="change_point_detection_for_cna", action="store_true", help="use change point detection algo for more cna segmentation instead of breakpoints")
    global_parser.add_argument('--consider-wgd', action="store_true", dest="consider_wgd", default=False, help="Consider half peak in first copy estimation and optimization for WGD")
    global_parser.add_argument('--copynumbers-disable', action="store_true", dest="copynumbers_disable", default=False, help="Disabling copy number in coverage plots")
    global_parser.add_argument('--copynumbers-subclonal-enable', action="store_true", dest="copynumbers_subclonal_enable", default=True, help="Enabling subclonal copy number in coverage plots")

    global_parser.add_argument("--purity-range", dest="purity_range", default=DEFAULT_PURITY, help="Estimated tumor purity range (fraction) between [default: 0.5-1.0]")
    global_parser.add_argument("--ploidy-range", dest="ploidy_range", default=DEFAULT_PLOIDY, help="Estimated tumor ploidy range between [default: 1.0-5.5]")
    global_parser.add_argument('--dna-purity', action="store_true", dest="dna_purity", default=False, help="Enabling DNA purity instead of default cell purity estimation")
    global_parser.add_argument("--tumor-purity", dest="tumor_purity", default=0.0, metavar="float", type=float, help="user input tumor purity")
    global_parser.add_argument("--tumor-ploidy", dest="tumor_ploidy", default=0.0, metavar="float", type=float, help="user input tumor ploidy")
    global_parser.add_argument("--confidence-subclonal-score", dest="confidence_subclonal_score", default=0.6, metavar="float", type=float, help="user input p-value to detect if a segment is subclonal/off to integer copynumber")

    global_parser.add_argument("--first-copy", dest="first_copy", default=0, metavar="int", type=int, help="simulation only first copy")
    global_parser.add_argument("--first-copy-breakpoints-filter", dest="first_copy_breakpoints_filter", default=0, metavar="int", type=int, help="Breakpoints filter only first copy")

    global_parser.add_argument("--breakpoints", dest="breakpoints", metavar="path", required=False, default=None, help="Path to breakpoints/SVs VCF file")
    global_parser.add_argument("--breakpoints-min-length", dest="breakpoints_min_length", default=BP_MIN_LENGTH, metavar="int", type=int, help="breakpoints minimum length to include [10k]")
    global_parser.add_argument("--use-sv-haplotypes", dest="use_sv_haplotypes", action="store_true", default=False, required=False, help="Enable using phased SVs/breakpoints")

    return global_parser

def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    missing = []
    if not args.tumor_phased_vcf and not args.normal_phased_vcf:
        missing.append("--normal-phased-vcf/--tumor-phased-vcf")
    if missing:
        logger.info(f"Missing required at least one normal/tumor phased VCF argument: {', '.join(missing)}")
        return 0
    if not args.command == 'hapcorrect' and not args.change_point_detection_for_cna and not args.breakpoints:
        logger.info('At least one parameter --breakpoints <SV VCF path> which is highly recommneded or --change-point-detection-for-cna should be used for copy number segmentation model.')
        return 0

    if not os.path.isdir(args.out_dir_plots):
        os.mkdir(args.out_dir_plots)

    log_file = os.path.join(args.out_dir_plots, "wakhan.log")
    _enable_logging(log_file, debug=False, overwrite=False)

    # Disable propagation for this logger
    logger.propagate = False

    logger.info("Starting Wakhan " + _version())
    logger.info("Cmd: %s", " ".join(sys.argv))
    logger.info("Python version: " + sys.version)

    fileDir = os.path.dirname(__file__)
    args.cancer_genes = os.path.join(fileDir, args.cancer_genes)
    args.centromere = os.path.join(fileDir, args.centromere)

    if os.path.exists(args.out_dir_plots+'/data'):
        shutil.rmtree(args.out_dir_plots+'/data')
        os.mkdir(args.out_dir_plots+'/data')
    else:
        os.mkdir(args.out_dir_plots+'/data')

    if args.command == 'cna':
        logger.info('Starting cna() module...')
        cna_process(args) # cna()
        logger.info('cna() module finished successfully.')
        return 0
    elif args.command == 'hapcorrect':
        logger.info('Starting hapcorrect() module...')
        main_process(args)  # hapcorrect()
        logger.info('hapcorrect() module finished successfully.')
        return 0
    elif args.command is None:
        wakhan_all(args) #hapcorrect + cna

def wakhan_all(args):
    logger.info('Starting hapcorrect() module...')
    main_process(args)  # hapcorrect
    logger.info('hapcorrect() module finished successfully.')

    logger.info('Starting cna() module...')
    cna_process(args) #cna
    logger.info('cna() module finished successfully.')

    return 0

def cna_process(args):
    thread_pool = Pool(args.threads)
    SAMTOOLS_BIN = "samtools"
    BCFTOOLS_BIN = "bcftools"
    MIN_SV_SIZE = 50

    if not args.quick_start:# and not args.without_phasing and not args.phaseblock_flipping_disable and not args.histogram_coverage:
        args.quick_start = True
        args.quick_start_coverage_path = args.out_dir_plots+'/coverage_data'

    if args.control_bam is None:
        args.control_bam = []
    all_bams = args.target_bam + args.control_bam
    target_genomes = set(os.path.basename(b) for b in args.target_bam)
    control_genomes = set(os.path.basename(b) for b in args.control_bam)

    if not shutil.which(SAMTOOLS_BIN):
        logger.info("samtools not found")
        return 1

    # TODO: check that all bams have the same reference
    first_bam = all_bams[0]
    ref_lengths = None
    with pysam.AlignmentFile(first_bam, "rb") as a:
        ref_lengths = dict(zip(a.references, a.lengths))

    if args.phaseblock_flipping_disable and args.histogram_coverage:
        segments_by_read = defaultdict(list)
        genome_ids = []
        for bam_file in all_bams:
            genome_id = os.path.basename(bam_file)
            genome_ids.append(genome_id)
            logger.info("Parsing reads from %s", args.target_bam)
            segments_by_read_bam = get_all_reads_parallel(bam_file, thread_pool, ref_lengths, args.min_mapping_quality, genome_id, MIN_SV_SIZE)
            segments_by_read.update(segments_by_read_bam)
            logger.info("Parsed %s segments", len(segments_by_read_bam))

        logger.info('Computing coverage histogram')
        coverage_histograms = update_coverage_hist(genome_ids, ref_lengths, segments_by_read, args.min_mapping_quality, args.max_read_error, args)
        del segments_by_read

    breakpoints_additional = extract_breakpoints_additional(args)
    if not args.phaseblock_flipping_disable:
       if args.enable_debug or (args.quick_start and os.path.exists(args.quick_start_coverage_path + '/phase_corrected_coverage.csv') and os.path.exists(args.quick_start_coverage_path + '/coverage_ps.csv') and os.path.exists(args.quick_start_coverage_path + '/coverage.csv')):
            logger.info('Using existing phase corrected coverage data')
       else:
           logger.info('Please use hapcorrect() module first to generate phase corrected coverage data')
           #main_process(args, breakpoints_additional, centromere_regions) #hapcorrect
       if args.quick_start:
           if args.without_phasing:
               df = pd.read_csv(args.quick_start_coverage_path + '/coverage_hps.csv', sep='\t', names=['chr', 'start', 'end', 'hp1', 'hp2', 'unphased'])
               df['hp1'] = df['hp1'] + df['hp2'] + df['un']
               df.to_csv(args.quick_start_coverage_path + '/coverage.csv', sep='\t', columns=['chr', 'start', 'end', 'hp1'], index=False, header=False)

               csv_df_coverage = csv_df_chromosomes_sorter(args.quick_start_coverage_path + '/coverage.csv', ['chr', 'start', 'end', 'coverage'])
               csv_df_phasesets = csv_df_chromosomes_sorter(args.quick_start_coverage_path + '/coverage_ps.csv', ['chr', 'start', 'end', 'coverage'])
           else:
               csv_df_phasesets = csv_df_chromosomes_sorter(args.quick_start_coverage_path + '/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
               csv_df_coverage = csv_df_chromosomes_sorter(args.quick_start_coverage_path+'/phase_corrected_coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
       else:
           csv_df_phasesets = csv_df_chromosomes_sorter(args.out_dir_plots+'/coverage_data/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
           csv_df_coverage = csv_df_chromosomes_sorter(args.out_dir_plots+'/coverage_data/phase_corrected_coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
    else:
        if args.quick_start and not args.without_phasing:
            csv_df_coverage = csv_df_chromosomes_sorter(args.quick_start_coverage_path + '/coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
            csv_df_phasesets = csv_df_chromosomes_sorter(args.quick_start_coverage_path + '/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
        #elif args.quick_start and args.without_phasing:
            #df = pd.read_csv(args.quick_start_coverage_path + '/coverage_hps.csv', sep='\t', names=['chr', 'start', 'end', 'hp1', 'hp2', 'unphased'])
            #df['unphased'] = df['hp1'] + df['hp2'] + df['unphased']
            #df.to_csv(args.quick_start_coverage_path + '/coverage.csv', sep='\t', columns=['chr', 'start', 'end', 'unphased'], index=False, header=False)

            #csv_df_coverage = csv_df_chromosomes_sorter(args.quick_start_coverage_path + '/coverage.csv', ['chr', 'start', 'end', 'coverage'])
            #csv_df_phasesets = csv_df_chromosomes_sorter(args.quick_start_coverage_path + '/coverage_ps.csv', ['chr', 'start', 'end', 'coverage'])
        elif args.histogram_coverage:
            if os.path.exists(args.out_dir_plots + '/coverage_data'):
                shutil.rmtree(args.out_dir_plots + '/coverage_data')
                os.mkdir(args.out_dir_plots + '/coverage_data')
            else:
                os.mkdir(args.out_dir_plots + '/coverage_data')
            logger.info('Computing coverage for bins')
            segments = get_chromosomes_bins_bam(args.target_bam[0], args.bin_size, args)
            segments_coverage = get_segments_coverage(segments, coverage_histograms)
            logger.info('Writing coverage for bins')
            write_segments_coverage_dict(segments_coverage, 'coverage.csv', args)
            logger.info('Parsing phaseblocks information')
            if args.normal_phased_vcf:
                output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(args.normal_phased_vcf, args)
            else:
                output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(args.tumor_phased_vcf, args)
            if not args.without_phasing:
                phasesets_segments = generate_phasesets_bins(args.target_bam[0], output_phasesets_file_path, args.bin_size, args) #TODO update for multiple bam files
                logger.info('Computing coverage for phaseblocks')
                phasesets_coverage = get_segments_coverage(phasesets_segments, coverage_histograms)
                logger.info('Writing coverage for phaseblocks')
                write_segments_coverage_dict(phasesets_coverage, 'coverage_ps.csv', args)
            del coverage_histograms
            logger.info('Loading coverage (bins) and coverage (phaseblocks) files...')
            if args.without_phasing:
                df = pd.read_csv(args.out_dir_plots+'/coverage_data/coverage.csv', sep='\t', names=['chr', 'start', 'end', 'hp1', 'hp2', 'unphased'])
                df['unphased'] = df['hp1'] + df['hp2'] + df['unphased']
                df.to_csv(args.out_dir_plots+'/coverage_data/coverage.csv', sep='\t', columns=['chr', 'start', 'end', 'unphased'], index=False, header=False)

                csv_df_coverage = csv_df_chromosomes_sorter(args.out_dir_plots+'/coverage_data/coverage.csv', ['chr', 'start', 'end', 'coverage'])
                csv_df_phasesets = pd.DataFrame(columns=['chr', 'start', 'end', 'coverage'])#csv_df_chromosomes_sorter(args.out_dir_plots+'/coverage_data/coverage_ps.csv', ['chr', 'start', 'end', 'coverage'])
            else:
                csv_df_phasesets = csv_df_chromosomes_sorter(args.out_dir_plots+'/coverage_data/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
                csv_df_coverage = csv_df_chromosomes_sorter(args.out_dir_plots+'/coverage_data/coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

    #TODO add chrX,chrY support later on
    #csv_df_coverage = csv_df_coverage.drop(csv_df_coverage[(csv_df_coverage.chr == "chrX") | (csv_df_coverage.chr == "chrY")].index)
    #csv_df_phasesets = csv_df_phasesets.drop(csv_df_phasesets[(csv_df_phasesets.chr == "chrX") | (csv_df_phasesets.chr == "chrY")].index)

    logger.info('Generating coverage plots chromosomes-wise')
    haplotype_1_values_updated, haplotype_2_values_updated, unphased, csv_df_snps_mean, snps_cpd_means, snps_cpd_points_weights, snps_cpd_means_df = \
        coverage_plots_chromosomes(csv_df_coverage, csv_df_phasesets, args, thread_pool)

    if not args.without_phasing:
        df_hp1, df_hp2, df_unphased = seperate_dfs_coverage(args, csv_df_snps_mean, csv_df_snps_mean.hp1.tolist(), csv_df_snps_mean.hp2.tolist(), csv_df_snps_mean.hp3.tolist())

    logger.info('Generating optimal clusters plots for bins')
    if args.without_phasing:
        df_hp1 = csv_df_snps_mean
        df_hp2 = csv_df_snps_mean

    # seg_min = []
    # seg_min_weights = []
    # for i in range(len(snps_cpd_means)):
    #     if snps_cpd_points_weights[i] * args.bin_size > 20000000:
    #         seg_min.append(snps_cpd_means[i])
    #         seg_min_weights.append(snps_cpd_points_weights[i])
    # if len(seg_min) > 1:
    #     max_limit = int(weighted_means(seg_min, seg_min_weights))
    #     logger.info('Max limit for normal optimization through 20Mb segments: %s', max_limit)
    # else:
    max_limit = int(weighted_means(snps_cpd_means, snps_cpd_points_weights))
    logger.info('Max limit for normal optimization: %s', max_limit)

    if args.histogram_coverage:
        tumor_cov = statistics.mean([sum(x) for x in zip(haplotype_1_values_updated, haplotype_2_values_updated, unphased)])
    else:
        tumor_cov = statistics.mean([sum(x) for x in zip(haplotype_1_values_updated, haplotype_2_values_updated)])
    logger.info('Tumor coverage: %s', tumor_cov)

    centers, is_half_peak, centers_half, subclonals, x_axis, observed_hist, single_copy_cov, single_copy_cov_half = peak_detection_optimization(args, snps_cpd_means, snps_cpd_points_weights, tumor_cov)
    args.first_copy_breakpoints_filter = single_copy_cov

    logger.info('Initial detected clusters means: %s', centers)
    args.cut_threshold = single_copy_cov * len(centers)

    #df_segs_hp1, df_segs_hp2 = breakpoints_segments_means(args, csv_df_snps_mean)
    haplotype_1_values_updated, haplotype_2_values_updated, unphased, csv_df_snps_mean, snps_cpd_means, snps_cpd_points_weights, snps_cpd_means_df = breakpoints_segments_means(csv_df_coverage, csv_df_phasesets, args, thread_pool)
    if args.without_phasing:
        df_segs_hp1 = snps_cpd_means_df
        df_segs_hp2 = snps_cpd_means_df
        df_hp1 = csv_df_snps_mean
        df_hp2 = csv_df_snps_mean
    else:
        df_segs_hp1 = snps_cpd_means_df[0]
        df_segs_hp2 = snps_cpd_means_df[1]

    #TODO Debug
    #df_snps_in_csv = pd.DataFrame(columns=['chr', 'pos', 'vaf'])

    if os.path.exists(args.quick_start_coverage_path + '/baf.csv'):
        df_snps_in_csv = csv_df_chromosomes_sorter(args.quick_start_coverage_path + '/baf.csv', ['chr', 'pos', 'vaf'], ',')
    else:
        df_snps_in_csv = snps_df_loh(args, thread_pool, df_hp1)
        df_snps_in_csv.to_csv(args.out_dir_plots + '/coverage_data/'+'baf.csv', index=False, header=None)

    if args.without_phasing:
        logger.info('Generating coverage/copy numbers plots genome wide')
        integer_fractional_means = sorted([i for i in range(0, len(centers))])
        df_segs_hp1_updated, df_segs_hp2_updated = adjust_diversified_segments(centers, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)
        df_segs_hp1_updated = merge_adjacent_regions_cn_unphased(df_segs_hp1_updated, args)
        df_segs_hp2_updated = merge_adjacent_regions_cn_unphased(df_segs_hp2_updated, args)
        if args.quick_start:
            loh_regions = pd.DataFrame(columns=['chr', 'start', 'end'])
        else:
            #TODO remove it, temp for debug
            #loh_regions = pd.DataFrame(columns=['chr', 'start', 'end'])
            loh_regions = plot_snps_frequencies_without_phasing(args, csv_df_snps_mean, df_segs_hp1_updated, df_segs_hp2_updated, centers, integer_fractional_means)
        write_copynumber_segments_csv(df_hp1, df_hp1, df_segs_hp1_updated, args, centers, integer_fractional_means, None, '_copynumbers_segments.bed', None, False)

        ###############coverage data save#################
        # if not 'state' in df_hp1.columns:
        #     states_coverage = []
        #     coverage_only = df_hp1.coverage.values.tolist()
        #     for i, (start, end) in enumerate(zip(df_hp1.start.values.tolist(), df_hp1.end.values.tolist())):
        #         states_coverage.append(min(centers, key=lambda x: abs(x - coverage_only[i])))
        #     df_hp1['state'] = states_coverage
        #     df = df_hp1.copy()
        # write_copynumber_segments_csv(df, args, centers, integer_fractional_means, None, '_coverage_segments.bed', None)
        ##################################################
        copy_number_plots_genome_details(centers, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args, x_axis, observed_hist, None, False)
        if args.breakpoints:
            copy_number_plots_genome_breakpoints(centers, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args, None, loh_regions, df_snps_in_csv, False)
        else:
            copy_number_plots_genome(centers, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args, None, loh_regions, df_snps_in_csv, False)

        if args.copynumbers_subclonal_enable:
            df_segs_hp1_updated, df_segs_hp2_updated = adjust_diversified_segments(centers, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args)
            df_segs_hp1_updated, df_segs_hp2_updated = update_subclonal_means_states(centers, subclonals, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2, args)
            #df_segs_hp1_updated = merge_adjacent_regions_cn(df_segs_hp1_updated, args)
            #df_segs_hp2_updated = merge_adjacent_regions_cn(df_segs_hp2_updated, args)
            write_copynumber_segments_csv(df_hp1, df_hp1, df_segs_hp1_updated, args, centers, integer_fractional_means, None, '_copynumbers_subclonal_segments.bed', None, False)
            cen_out = [int(i) for i in centers]
            if args.breakpoints:
                copy_number_plots_genome_breakpoints_subclonal(cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args, None, loh_regions, df_snps_in_csv, False)
            else:
                copy_number_plots_genome_subclonal(cen_out, integer_fractional_means, df_hp1, df_segs_hp1_updated, df_hp2, df_segs_hp2_updated, df_hp1, args, None, loh_regions, df_snps_in_csv, False)
    else:
        ################################
        copy_numbers_assignment_haplotypes(args, tumor_cov, max_limit, single_copy_cov, centers, subclonals, df_hp1, df_hp2, df_segs_hp1, df_segs_hp2, snps_cpd_means_df, csv_df_snps_mean, df_snps_in_csv, df_unphased, x_axis, observed_hist, False)
        if is_half_peak:
            copy_numbers_assignment_haplotypes(args, tumor_cov, max_limit, single_copy_cov_half, centers_half, subclonals, df_hp1, df_hp2, df_segs_hp1, df_segs_hp2, snps_cpd_means_df, csv_df_snps_mean, df_snps_in_csv, df_unphased, x_axis, observed_hist, False)

        if len(solutions_df):
            solutions_df['solution_rank'] = solutions_df['confidence'].rank(method='first', ascending=False).astype(int)
            solutions_df.sort_values(by='solution_rank', inplace=True)
            solutions_df.to_csv(args.out_dir_plots+'/solutions_ranks.tsv', sep='\t', header=True, index=False, mode='w')

            for index, link_info in solutions_df.iterrows():
                target_folder = os.path.abspath(args.out_dir_plots) + '/' + link_info['repository_name']
                link_path = os.path.abspath(args.out_dir_plots) + '/' + 'solution_' + str(link_info['solution_rank'])
                try:
                    os.symlink(target_folder, link_path, target_is_directory=True)
                    logger.info(f"Solution symbolic link created: {link_path} -> {target_folder}")
                except OSError as e:
                    logger.info(f"Error creating symbolic link: {e}")

    #if average_p_value:
    #    #SNPs ratios and LOH and plots
    #    plot_snps_ratios_genome(args, df_snps_in_csv, loh_regions)
    ################################
    if os.path.exists(args.out_dir_plots+'/data'): #
        shutil.rmtree(args.out_dir_plots+'/data')
    if os.path.exists(args.out_dir_plots+'/data_phasing'):
        shutil.rmtree(args.out_dir_plots+'/data_phasing')

    #TODO, not required yet
    #move_100pct_purity_sol(args)

    return 0

#UCSC tumor/normal celllines
#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots 1954  --normal-phased-vcf /home/rezkuh/gits/data/1954/1954BL.vcf.gz --copynumbers-enable True  --unphased-reads-coverage-enable True  --phaseblock-flipping-enable True --phaseblocks-enable True   --genome-name 1954  --cut-threshold 150

#pancreatic_organoid data
#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots  --phased-vcf /home/rezkuh/gits/data/pancreatic_organoid/pancBL.vcf.gz  --copynumbers-enable True  --unphased-reads-coverage-enable True --snps-freq-vcf-enable True --phaseblock-flipping-enable True --phaseblocks-enable True  --genome-name pancreatic_organoid  --cut-threshold 150

#Tumor only (HPV)
#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots --snps-freq-vcf-enable True --cut-threshold 50 --copynumbers-enable True --phaseblock-flipping-enable True   --snps-freq-vcf-enable True  --phased-vcf /home/rezkuh/gits/data/R10/HT3/HT3.vcf.gz  --genome-name R10/HT3
#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots --snps-freq-vcf-enable True --cut-threshold 50 --copynumbers-enable True --phaseblock-flipping-enable True   --snps-freq-vcf-enable True  --phased-vcf /home/rezkuh/gits/data/R10/CaSki/CaSki.vcf.gz  --genome-name R10/CaSki

#NIST GIAB HG008
#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots HG008_HiFi  --normal-phased-vcf /home/rezkuh/gits/data/HG008_HiFi/HG008BL_HiFi.vcf.gz --tumor-vcf /home/rezkuh/gits/data/HG008_HiFi/HG008_HiFi.vcf.gz --copynumbers-enable True  --unphased-reads-coverage-enable True  --phaseblock-flipping-enable True --phaseblocks-enable True   --genome-name HG008_HiFi  --cut-threshold 150

#colo829
#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots colo829  --normal-phased-vcf /home/rezkuh/gits/data/colo829/colo829_pepper_normal.phased.vcf.gz  --copynumbers-enable True   --phaseblock-flipping-enable True --phaseblocks-enable True   --genome-name colo829  --cut-threshold 150

#colo829-porec
#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots colo829-porec  --normal-phased-vcf /home/rezkuh/gits/data/colo829-porec/colo829.vcf.gz  --copynumbers-enable True    --phaseblocks-enable True   --genome-name colo829-porec  --cut-threshold 150

#Mouse unphased data
#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam    --copynumbers-enable True    --unphased-reads-coverage-enable True    --cut-threshold 75  --without-phasing True --tumor-vcf /home/rezkuh/gits/data/mouse/somatic_calls/C1_somatic_calls_pass_snp.vcf.gz --out-dir-plots C1 --contigs 1-19 --genome-name C1 --bin-size-snps 1000000
#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam    --copynumbers-enable True    --unphased-reads-coverage-enable True    --cut-threshold 75  --without-phasing True --tumor-vcf /home/rezkuh/gits/data/mouse/somatic_calls/C23_somatic_calls_pass_snp.vcf.gz --out-dir-plots C23 --contigs 1-19,X,Y --genome-name C23 --bin-size-snps 1000000

#Dog data
#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam   --normal-phased-vcf /home/rezkuh/gits/data/OT4/ON2.vcf.gz --tumor-vcf /home/rezkuh/gits/data/OT2/OT2.vcf.gz     --copynumbers-enable True    --genome-name OT2 --out-dir-plots OT2  --cut-threshold 60 --phaseblock-flipping-enable True --phaseblocks-enable True --contigs chr1-38

#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots R10/colo357_R10  --tumor-vcf /home/rezkuh/gits/data/R10/colo357_R10/colo357.vcf.gz --copynumbers-enable True  --unphased-reads-coverage-enable True  --phaseblocks-enable True   --genome-name R10/colo357_R10 --contigs chr1-22,X --cut-threshold 150 --breakpoints /home/rezkuh/gits/data/R10/colo357_R10/

#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots colo829 --genome-name colo829  --cut-threshold 200   --normal-phased-vcf /home/rezkuh/gits/data/colo829/colo829BL.vcf.gz --breakpoints /home/rezkuh/gits/data/colo829/severus_somatic.vcf
#--quick-start True --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots colo357_R10  --copynumbers-enable True  --unphased-reads-coverage-enable True    --genome-name colo357_R10  --cut-threshold 150  --tumor-vcf /home/rezkuh/gits/data/R10/colo357_R10/colo357.vcf.gz  --phaseblock-flipping-enable True --breakpoints /home/rezkuh/gits/data/R10/colo357_R10/

#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --cut-threshold 150  --normal-phased-vcf /home/rezkuh/gits/data/1437/1437BL.vcf.gz   --breakpoints /home/rezkuh/gits/data/1437/severus_somatic.vcf  --out-dir-plots 1437_merged_merged --genome-name 1437_merged_merged --copynumbers-subclonal-enable --loh-enable --purity-range 0.5-1.0 --ploidy-range 2-8
#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --cut-threshold 150  --normal-phased-vcf /home/rezkuh/gits/data/2009/2009BL.vcf.gz   --breakpoints /home/rezkuh/gits/data/2009/severus_somatic.vcf  --out-dir-plots 2009_merged --genome-name 2009_merged --copynumbers-subclonal-enable --loh-enable --purity-range 0.5-1.0 --ploidy-range 2-8
#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/C15_1_30000000.bam --cut-threshold 75  --without-phasing --phaseblock-flipping-disable --tumor-vcf /home/rezkuh/gits/data/mouse/somatic_calls/C23_somatic_calls_pass_snp.vcf.gz --breakpoints /home/rezkuh/gits/data/mouse/severus_somatic.vcf  --out-dir-plots C23 --contigs 1-19,X --genome-name C23 --bin-size-snps 1000000 --centromere annotations/mouse.bed --hets-smooth-window 10 --copynumbers-subclonal-enable --loh-enable

#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/colo357_R10 --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam  --out-dir-plots colo357_R10  --tumor-vcf /home/rezkuh/gits/data/colo357_R10/colo357_phased.vcf.gz  --phaseblocks-enable  --genome-name colo357_R10  --cut-threshold 75  --breakpoints /home/rezkuh/gits/data/colo357_R10/severus_somatic.vcf --contigs chr1-22 --copynumbers-subclonal-enable --loh-enable  --purity-range 0.5-1.0 --ploidy-range 1-4

#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/ --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --cut-threshold 50  --tumor-vcf /home/rezkuh/gits/data/dbg/longphase.vcf.gz   --breakpoints /home/rezkuh/gits/data/dbg/severus_somatic.vcf  --out-dir-plots dbg --genome-name dbg --copynumbers-subclonal-enable --loh-enable --phaseblock-flipping-disable --contigs chr1-22,X,Y

#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/2009_merged --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --cut-threshold 150  --normal-phased-vcf /home/rezkuh/gits/data/2009/2009BL.vcf.gz    --out-dir-plots 2009_merged --genome-name 2009_merged --copynumbers-subclonal-enable --loh-enable --purity-range 0.5-1.0 --ploidy-range 2-8 --phaseblocks-enable --change-point-detection-for-cna
#--quick-start  --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/C15_1_30000000.bam --cut-threshold 75  --without-phasing --phaseblock-flipping-disable    --contigs 1-19,X  --bin-size-snps 1000000 --centromere annotations/mouse.bed --hets-smooth-window 10 --copynumbers-subclonal-enable --loh-enable --tumor-vcf  /home/rezkuh/gits/data/mouse/somatic_calls/C23_somatic_calls_pass_snp.vcf.gz --breakpoints /home/rezkuh/gits/data/mouse/C23.haplotagged.vcf  --quick-start-coverage-path /home/rezkuh/gits/data/C23 --genome-name C23 --cpd-internal-segments  --out-dir-plots C23

#1437 BPs
#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/1437_merged --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --cut-threshold 150  --normal-phased-vcf /home/rezkuh/gits/data/1437/1437BL.vcf.gz    --out-dir-plots 1437_merged --genome-name 1437_merged --copynumbers-subclonal-enable --loh-enable --purity-range 0.5-1.0 --ploidy-range 2-8 --phaseblocks-enable --breakpoints /home/rezkuh/gits/data/1437/severus_somatic.vcf
#1437 BPs and input genes
#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/1437_merged --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --cut-threshold 150  --normal-phased-vcf /home/rezkuh/gits/data/1437/1437BL.vcf.gz    --out-dir-plots 1437_merged --genome-name 1437_merged --copynumbers-subclonal-enable --loh-enable --purity-range 0.5-1.0 --ploidy-range 2-8 --phaseblocks-enable --breakpoints /home/rezkuh/gits/data/1437/severus_somatic.vcf --user-input-genes /home/rezkuh/gits/backup/872024/Wakhan/src/annotations/user_input_genes.txt

#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/rouf3 --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --cut-threshold 50  --normal-phased-vcf /home/rezkuh/gits/data/rouf3/rouf3BL.vcf.gz    --out-dir-plots rouf3 --genome-name rouf3 --copynumbers-subclonal-enable --loh-enable --purity-range 0.5-1.0 --ploidy-range 0.5-8 --phaseblocks-enable --breakpoints /home/rezkuh/gits/data/rouf3/severus_somatic.vcf --breakpoints-min-length 100

#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/HG008_ONT_10k --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --cut-threshold 150  --normal-phased-vcf /home/rezkuh/gits/data/HG008/hic.phased.vcf.gz    --out-dir-plots /home/rezkuh/gits/backup/HG008_ONT_10k  --genome-name HG008_ONT_10k  --copynumbers-subclonal-enable --loh-enable  --phaseblocks-enable --breakpoints /home/rezkuh/gits/data/HG008_ONT_10k/severus_somatic.vcf --bin-size 10000

#Wakhan TODOs
#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/1437/80pct.60x --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam  --out-dir-plots 1437_80pct.60x_unphased_bam  --normal-phased-vcf /home/rezkuh/gits/data/1437/80pct.60x/BL1437.ONT.30x.longphase.vcf.gz  --phaseblocks-enable  --genome-name 1437_80pct_60x  --cut-threshold 100  --breakpoints /home/rezkuh/gits/data/1437/80pct.60x/80pct.60x_severus_somatic.vcf --contigs chr1-22 --copynumbers-subclonal-enable --loh-enable  --purity-range 0.5-1.0 --ploidy-range 1-4
#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/1437/80pct.60x_to --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam  --out-dir-plots 1437_80pct.60x_to_unphased_bam  --tumor-vcf /home/rezkuh/gits/data/1437/80pct.60x_to/80pct.60x_longphase.vcf.gz  --phaseblocks-enable  --genome-name 1437_80pct_60x  --cut-threshold 100  --breakpoints /home/rezkuh/gits/data/1437/80pct.60x_to/80pct.60x_severus_somatic.vcf --contigs chr1-22 --copynumbers-subclonal-enable --loh-enable  --purity-range 0.5-1.0 --ploidy-range 1-4
#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/1954/60pct.20x_to --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam  --out-dir-plots 1954_60pct.20x_to  --tumor-vcf /home/rezkuh/gits/data/1954/60pct.20x_to/longphase.vcf.gz  --phaseblocks-enable  --genome-name 60pct.20x_to  --cut-threshold 100  --breakpoints /home/rezkuh/gits/data/1954/60pct.20x_to/severus_somatic.vcf --contigs chr1-22 --copynumbers-subclonal-enable --loh-enable --bin-size-snps 100000

#--quick-start --quick-start-coverage-path /home/rezkuh/gits/data/BT474_err --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam  --out-dir-plots BT474_err  --tumor-vcf /home/rezkuh/gits/data/BT474_err/longphase.vcf.gz  --phaseblocks-enable  --genome-name BT474_err  --cut-threshold 100  --breakpoints /home/rezkuh/gits/data/BT474_err/severus_somatic.vcf --contigs chr1-22 --copynumbers-subclonal-enable --loh-enable --bin-size-snps 100000 --bin-size 10000  --hets-ratio 0.25 --ploidy-range 1-6 --cpd-internal-segments
#cna --quick-start --quick-start-coverage-path /home/rezkuh/gits/data/SKBR3 --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots SKBR3 --tumor-phased-vcf /home/rezkuh/gits/data/SKBR3/longphase.vcf.gz  --genome-name SKBR3 --cut-threshold 70  --loh-enable  --bin-size 10000 --phaseblocks-enable --breakpoints /home/rezkuh/gits/data/SKBR3/severus_somatic.vcf