#!/usr/bin/env python3

import sys
import shutil

import numpy
import pysam
import argparse
import logging
import os

from multiprocessing import Pool
from collections import defaultdict
from bam_processing import get_all_reads_parallel, update_coverage_hist, get_segments_coverage, haplotype_update_all_bins_parallel, get_snps_frequencies
from cnvlib import descriptives
from cnvlib import cluster
from utils import get_chromosomes_bins, write_segments_coverage, csv_df_chromosomes_sorter, generate_phasesets_bins, csv_df_chromosomes_sorter_snps,\
    apply_copynumbers, csv_df_chromosomes_sorter_copyratios, seperate_dfs_coverage, flatten_smooth, get_contigs_list, \
    csv_df_chromosomes_sorter_snps_from_bam
from plots import coverage_plots_chromosomes, copy_number_plots_genome, plots_genome_coverage, copy_number_plots_chromosomes
from vcf_processing import vcf_parse_to_csv_for_het_phased_snps_phasesets, get_snp_segments
from bam_processing import process_bam_for_snps_freqs

#remove
from utils import get_chromosomes_bins_replica

def main():
    # default tunable parameters
    MAX_READ_ERROR = 0.1
    MIN_MAPQ = 10
    MIN_SV_SIZE = 50
    BIN_SIZE = 50000
    MAX_CUT_THRESHOLD = 100
    MIN_ALIGNED_LENGTH = 5000

    SAMTOOLS_BIN = "samtools"
    BCFTOOLS_BIN = "bcftools"

    DEFAULT_CONTIGS = 'chr1-22'

    parser = argparse.ArgumentParser \
        (description="Find breakpoints and build breakpoint graph from a bam file")

    parser.add_argument("--target-bam", dest="target_bam",
                        metavar="path", required=True, default=None, nargs="+",
                        help="path to one or multiple target haplotagged bam files (must be indexed)")
    parser.add_argument("--control-bam", dest="control_bam",
                        metavar="path", required=False, default=None, nargs="+",
                        help="path to one or multiple control haplotagged bam files (must be indexed)")
    parser.add_argument("--reference", dest="reference",
                        metavar="path", required=False, default=None,
                        help="path to reference")

    parser.add_argument("--out-dir-plots", dest="out_dir_plots",
                        default=None, required=True,
                        metavar="path", help="Output plots directory")

    parser.add_argument("--phased-vcf", dest="phased_vcf",
                        metavar="path", required=True, default=None,
                        help="Path to phased vcf")
    parser.add_argument("--phased-vcf-snps-freqs", dest="phased_vcf_snps_freqs",
                        metavar="path", required=False, default=None,
                        help="Path to phased vcf to plot snps frequencies coverage")

    parser.add_argument("--breakpoints-file", dest="breakpoints_file",
                        metavar="path", required=False, default=None,
                        help="Path to breakpoints file to plot inconjunction with coverage plots")

    parser.add_argument("--genome-name", dest="genome_name",
                        required=True, default=None,
                        help="Genome sample/cellline name to be displayed on plots")
    parser.add_argument("--contigs", dest="contigs",
                        required=False, default=DEFAULT_CONTIGS,
                        help="List of contigs (choromosomes) to be included in the plots [e.g., chr1-22,X,Y]")

    parser.add_argument("--bin-size", "--bin_size", dest="bin_size",
                        default=BIN_SIZE, metavar="int", type=int, help="coverage (readdepth) bin size [50k]")
    parser.add_argument("--cut-threshold", "--cut_threshold", dest="cut_threshold",
                        default=MAX_CUT_THRESHOLD, metavar="int", type=int, help="Maximum cut threshold for coverage (readdepth) [100]")
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
                        default=False, help="Enabling unphased reads coverage output in plots")

    parser.add_argument('--phaseblock-flipping-enable',  dest="phaseblock_flipping_enable", required=False,
                        default=False, help="Enabling phaseblock flipping in coverage plots")
    parser.add_argument('--smoothing-enable',  dest="smoothing_enable", required=False,
                        default=False, help="Enabling smoothing in coverage plots")
    parser.add_argument('--phaseblocks-enable',  dest="phaseblocks_enable", required=False,
                        default=False, help="Enabling phaseblocks display in coverage plots")
    parser.add_argument('--het-phased-snps-freq-enable',  dest="het_phased_snps_freq_enable", required=False,
                        default=False, help="Enabling hetrozygous phased snps frequencies in coverage plots")
    parser.add_argument('--breakpoints-enable',  dest="breakpoints_enable", required=False,
                        default=False, help="Enabling breakpoints in coverage plots")
    parser.add_argument('--copynumbers-enable', dest="copynumbers_enable", required=False,
                        default=False, help="Enabling copy number in coverage plots")
    parser.add_argument('--copynumbers-subclonal-enable', dest="copynumbers_subclonal_enable", required=False,
                        default=False, help="Enabling subclonal copy number in coverage plots")

    parser.add_argument("-t", "--threads", dest="threads",
                        default=1, metavar="int", type=int, help="number of parallel threads [8]")

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
        "smoothing_enable": args.smoothing_enable,
        "phaseblocks_enable": args.phaseblocks_enable,
        "copynumbers_enable": args.copynumbers_enable,
        "copynumbers_subclonal_enable": args.copynumbers_subclonal_enable,
        "het_phased_snps_freq_enable": args.het_phased_snps_freq_enable,
        "out_dir_plots": args.out_dir_plots,
        "phased_vcf": args.phased_vcf,
        "phased_vcf_snps_freqs": args.phased_vcf_snps_freqs,
        "breakpoints_file": args.breakpoints_file,
        "breakpoints_enable": args.breakpoints_enable,
        "genome_name": args.genome_name,
        "bin_size": args.bin_size,
        "pdf_enable": args.pdf_enable,
        "cut_threshold": args.cut_threshold,
        "min_aligned_length": args.min_aligned_length,
        "no_of_clusters": args.no_of_clusters,
        "contigs": args.contigs,
        "threads": args.threads,
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
    if not os.path.isdir('data'):
        os.mkdir('data')

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

    if arguments['het_phased_snps_freq_enable']:
        get_snp_segments(arguments, args.target_bam[0], thread_pool)
        csv_df_snps = csv_df_chromosomes_sorter_snps_from_bam('data/snps_frequencies.csv')

    logging.info('Computing coverage histogram')
    coverage_histograms = update_coverage_hist(genome_ids, ref_lengths, segments_by_read, args.min_mapping_quality,
                                               args.max_read_error, arguments)

    logging.info('Computing coverage for bins')
    segments = get_chromosomes_bins(args.target_bam[0], arguments['bin_size'], arguments)
    #segments.append(('colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam', 'chr7', 78318498, 78486891))
    segments_coverage = get_segments_coverage(segments, coverage_histograms)
    logging.info('Writing coverage for bins')
    write_segments_coverage(segments_coverage, 'coverage.csv')

    logging.info('Parsing phaseblocks information')
    output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(arguments['phased_vcf'])
    phasesets_segments = generate_phasesets_bins(args.target_bam[0], output_phasesets_file_path, arguments['bin_size'], arguments) #TODO update for multiple bam files
    logging.info('Computing coverage for phaseblocks')
    phasesets_coverage = get_segments_coverage(phasesets_segments, coverage_histograms)
    logging.info('Writing coverage for phaseblocks')
    write_segments_coverage(phasesets_coverage, 'coverage_ps.csv')

    logging.info('Loading coverage (bins) and coverage (phaseblocks) files...')
    csv_df_phasesets = csv_df_chromosomes_sorter('data/coverage_ps.csv')
    csv_df_coverage = csv_df_chromosomes_sorter('data/coverage.csv')
    #csv_df_phasesets = csv_df_chromosomes_sorter('/home/rezkuh/gits/data/'+arguments['genome_name']+'/coverage_ps.csv')
    #csv_df_coverage = csv_df_chromosomes_sorter('/home/rezkuh/gits/data/'+arguments['genome_name']+'/coverage.csv')

    logging.info('Generating coverage plots chromosomes-wise')
    haplotype_1_values_updated, haplotype_2_values_updated, unphased, csv_df_snps_mean = \
        coverage_plots_chromosomes(csv_df_coverage, csv_df_phasesets, csv_df_snps, arguments)

    csv_df_coverage = csv_df_coverage.drop(csv_df_coverage[(csv_df_coverage.chr == "chrX") | (csv_df_coverage.chr == "chrY")].index)
    df_hp1, df_hp2, df_unphased = seperate_dfs_coverage(csv_df_coverage, haplotype_1_values_updated, haplotype_2_values_updated, unphased)

    #cluster.plot_optimal_clusters(haplotype_1_values_updated, haplotype_1_values_updated, unphased,  arguments, "coverage")

    logging.info('Generating coverage plots genome wide')
    plots_genome_coverage(df_hp1, df_hp2, df_unphased, arguments)

    logging.info('Generating optimal clusters plots for bins')

    #haplotype_1_values_updated, haplotype_2_values_updated, unphased = flatten_smooth(haplotype_1_values_updated, haplotype_2_values_updated, unphased)
    #cluster.plot_optimal_clusters(haplotype_1_values_updated, haplotype_2_values_updated, unphased,  arguments)

    #TODO Covergae CopyNumbers
    #df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2 = apply_copynumbers(csv_df_coverage, haplotype_1_values_updated, haplotype_2_values_updated, arguments)
    #csv_df_snps_mean_clean = csv_df_snps_mean[(csv_df_snps_mean['hp1'] > 3) & (csv_df_snps_mean['hp2'] > 3)]

    cluster.plot_optimal_clusters(csv_df_snps_mean.hp1.values, csv_df_snps_mean.hp2.values, unphased,  arguments, None)
    #cluster.plot_optimal_clusters(csv_df_snps_mean.hp1.clip(upper=arguments['cut_threshold']).values, csv_df_snps_mean.hp2.clip(upper=arguments['cut_threshold']).values, unphased,  arguments, None)

    ############################################
    #csv_df_snps_mean.to_csv('data/'+arguments['genome_name']+'_snps.csv', sep='\t')
    ############################################
    #TODO SNPs CopyNumbers
    df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2 = apply_copynumbers(csv_df_snps_mean, csv_df_snps_mean.hp1.values.tolist(), csv_df_snps_mean.hp2.values.tolist(), arguments)

    logging.info('Generating copy number log2 ratios plots chromosomes-wise')
    #copy_number_plots_chromosomes(df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, arguments)

    logging.info('Generating coverage/copy numbers plots genome wide')
    copy_number_plots_genome(df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, df_unphased, arguments)

if __name__ == "__main__":
    main()

#--phaseblock-flipping-enable True
# --phaseblocks-enable True --unphased-reads-coverage-enable True

#--smoothing-enable True
#--pdf-enable True
#--copynumbers-enable True

