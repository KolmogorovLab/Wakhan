#!/usr/bin/env python3

import sys
import shutil
import pysam
import argparse
import logging
import os
import pandas as pd

from multiprocessing import Pool
from collections import defaultdict
from bam_processing import get_all_reads_parallel, update_coverage_hist, get_segments_coverage, haplotype_update_all_bins_parallel, get_snps_frequencies
from cnvlib import descriptives
from cnvlib import cluster
from utils import get_chromosomes_bins, write_segments_coverage, write_segments_coverage_dict, csv_df_chromosomes_sorter,\
    apply_copynumbers, seperate_dfs_coverage, flatten_smooth, get_contigs_list, write_copynumber_segments_csv, integer_fractional_cluster_means, \
    adjust_diversified_segments
from plots import coverage_plots_chromosomes, copy_number_plots_genome, plots_genome_coverage, copy_number_plots_chromosomes
from vcf_processing import vcf_parse_to_csv_for_het_phased_snps_phasesets
from snps_loh import plot_snps_frequencies
from phasing_correction import generate_phasesets_bins

#remove
from utils import get_chromosomes_bins_replica

def main():
    # default tunable parameters
    MAX_READ_ERROR = 0.1
    MIN_MAPQ = 10
    MIN_SV_SIZE = 50
    BIN_SIZE = 50000
    BIN_SIZE_SNPS = 50000
    MAX_CUT_THRESHOLD = 100
    MIN_ALIGNED_LENGTH = 5000

    SAMTOOLS_BIN = "samtools"
    BCFTOOLS_BIN = "bcftools"

    DEFAULT_CONTIGS = 'chr1-22' #('chr1-22' '1-22') ('chr1-22,chrX' '1-22,X')

    parser = argparse.ArgumentParser \
        (description="Plot coverage and copy number profiles from a bam and phased VCF files")

    parser.add_argument("--target-bam", dest="target_bam",
                        metavar="path", required=True, default=None, nargs="+",
                        help="path to one or1/4 multiple target haplotagged bam files (must be indexed)")
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
                        metavar="path", required=True, default=None,
                        help="Path to normal phased vcf")
    parser.add_argument("--tumor-vcf", dest="tumor_vcf",
                        metavar="path", required=False, default=None,
                        help="Path to tumor VCF for LOH detection")

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
    parser.add_argument("--bin-size-snps", "--bin_size_snps", dest="bin_size_snps",
                        default=BIN_SIZE_SNPS, metavar="int", type=int, help="SNPs bin size [50k]")

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
    parser.add_argument('--without-phasing', dest="without_phasing", required=False,
                        default=False, help="Enabling coverage and copynumbers without phasing in plots")

    parser.add_argument('--phaseblock-flipping-enable',  dest="phaseblock_flipping_enable", required=False,
                        default=False, help="Enabling phaseblock flipping in coverage plots")
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
                        default=False, help="Enabling copy number in coverage plots")
    parser.add_argument('--copynumbers-subclonal-enable', dest="copynumbers_subclonal_enable", required=False,
                        default=False, help="Enabling subclonal copy number in coverage plots")

    parser.add_argument('--enable-debug', dest="enable_debug", required=False,
                        default=False, help="Enabling debug")

    parser.add_argument("-t", "--threads", dest="threads",
                        default=1, metavar="int", type=int, help="number of parallel threads [8]")
    parser.add_argument('--dryrun', dest="dryrun", required=False,
                        default=False, help="Enabling dryrun")

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
        "breakpoints_file": args.breakpoints_file,
        "breakpoints_enable": args.breakpoints_enable,
        "genome_name": args.genome_name,
        "bin_size": args.bin_size,
        "bin_size_snps": args.bin_size_snps,
        "pdf_enable": args.pdf_enable,
        "cut_threshold": args.cut_threshold,
        "min_aligned_length": args.min_aligned_length,
        "no_of_clusters": args.no_of_clusters,
        "contigs": args.contigs,
        "threads": args.threads,
        "dryrun": args.dryrun,
        "enable_debug": args.enable_debug,
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
    if os.path.exists('data'):
        shutil.rmtree('data')
        os.mkdir('data')
    else:
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

    logging.info('Computing coverage histogram')
    coverage_histograms = update_coverage_hist(genome_ids, ref_lengths, segments_by_read, args.min_mapping_quality,
                                               args.max_read_error, arguments)
    del segments_by_read

    if arguments['dryrun']:
        if arguments['without_phasing']:
            df = pd.read_csv('/home/rezkuh/gits/data/' + arguments['genome_name'] + '/coverage_hps.csv', sep='\t', names=['chr', 'start', 'end', 'hp1', 'hp2', 'un'])
            df['hp1'] = df['hp1'] + df['hp2'] + df['un']
            df.to_csv('/home/rezkuh/gits/data/' + arguments['genome_name'] + '/coverage.csv', sep='\t', columns=['chr', 'start', 'end', 'hp1'], index=False, header=False)

            csv_df_coverage = csv_df_chromosomes_sorter('/home/rezkuh/gits/data/' + arguments['genome_name'] + '/coverage.csv', ['chr', 'start', 'end', 'coverage'])
            csv_df_phasesets = csv_df_chromosomes_sorter('/home/rezkuh/gits/data/' + arguments['genome_name'] + '/coverage_ps.csv', ['chr', 'start', 'end', 'coverage'])

        else:
            csv_df_phasesets = csv_df_chromosomes_sorter('/home/rezkuh/gits/data/' + arguments['genome_name'] + '/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
            csv_df_coverage = csv_df_chromosomes_sorter('/home/rezkuh/gits/data/' + arguments['genome_name'] + '/coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
    else:
        logging.info('Computing coverage for bins')
        segments = get_chromosomes_bins(args.target_bam[0], arguments['bin_size'], arguments)
        segments_coverage = get_segments_coverage(segments, coverage_histograms)
        logging.info('Writing coverage for bins')
        write_segments_coverage_dict(segments_coverage, 'coverage.csv')

        logging.info('Parsing phaseblocks information')
        if arguments['normal_phased_vcf']:
            output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(arguments['normal_phased_vcf'])
        else:
            output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(arguments['tumor_vcf'])
        phasesets_segments = generate_phasesets_bins(args.target_bam[0], output_phasesets_file_path, arguments['bin_size'], arguments) #TODO update for multiple bam files
        logging.info('Computing coverage for phaseblocks')
        phasesets_coverage = get_segments_coverage(phasesets_segments, coverage_histograms)
        logging.info('Writing coverage for phaseblocks')
        write_segments_coverage_dict(phasesets_coverage, 'coverage_ps.csv')
        del coverage_histograms

        logging.info('Loading coverage (bins) and coverage (phaseblocks) files...')
        csv_df_phasesets = csv_df_chromosomes_sorter('data/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
        csv_df_coverage = csv_df_chromosomes_sorter('data/coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

    #TODO add chrX,chrY support later on
    csv_df_coverage = csv_df_coverage.drop(csv_df_coverage[(csv_df_coverage.chr == "chrX") | (csv_df_coverage.chr == "chrY")].index)
    csv_df_phasesets = csv_df_phasesets.drop(csv_df_phasesets[(csv_df_phasesets.chr == "chrX") | (csv_df_phasesets.chr == "chrY")].index)

    logging.info('Generating coverage plots chromosomes-wise')
    haplotype_1_values_updated, haplotype_2_values_updated, unphased, csv_df_snps_mean, snps_cpd_means, snps_cpd_means_df = \
        coverage_plots_chromosomes(csv_df_coverage, csv_df_phasesets, arguments, thread_pool)

    if arguments['without_phasing'] == False:
        df_hp1, df_hp2, df_unphased = seperate_dfs_coverage(arguments, csv_df_coverage, haplotype_1_values_updated, haplotype_2_values_updated, unphased)

    logging.info('Generating optimal clusters plots for bins')

    #haplotype_1_values_updated, haplotype_2_values_updated, unphased = flatten_smooth(haplotype_1_values_updated, haplotype_2_values_updated, unphased)
    #cluster.plot_optimal_clusters(haplotype_1_values_updated, haplotype_2_values_updated, unphased,  arguments)

    #TODO Covergae CopyNumbers
    #df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2 = apply_copynumbers(csv_df_coverage, haplotype_1_values_updated, haplotype_2_values_updated, arguments)
    #csv_df_snps_mean_clean = csv_df_snps_mean[(csv_df_snps_mean['hp1'] > 3) & (csv_df_snps_mean['hp2'] > 3)]

    #cluster.plot_optimal_clusters(csv_df_snps_mean.hp1.values, csv_df_snps_mean.hp2.values, unphased,  arguments, None)
    #cluster.plot_optimal_clusters(csv_df_snps_mean.hp1.clip(upper=arguments['cut_threshold']).values, csv_df_snps_mean.hp2.clip(upper=arguments['cut_threshold']).values, unphased,  arguments, None)

    ############################################
    #csv_df_snps_mean.to_csv('data/'+arguments['genome_name']+'_snps.csv', sep='\t')
    ############################################

    if arguments['without_phasing']:
        df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, states, centers, stdev = apply_copynumbers(csv_df_coverage, csv_df_coverage.coverage.values.tolist(), csv_df_coverage.coverage.values.tolist(), arguments, snps_cpd_means, snps_cpd_means)
        integer_fractional_means = integer_fractional_cluster_means(arguments, df_segs_hp1, df_segs_hp2, centers)
        logging.info('Correcting diversified segments')
        df_segs_hp1, df_segs_hp2 = adjust_diversified_segments(centers, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, arguments)
        logging.info('Generating copy number plots chromosomes-wise')
        logging.info('Generating coverage/copy numbers plots genome wide')
        write_copynumber_segments_csv(df_segs_hp1, arguments, centers, integer_fractional_means, None)
        centers, integer_fractional_centers = copy_number_plots_genome(centers, integer_fractional_means, df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, df_cnr_hp1, arguments)
    else:
        df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, states, centers, stdev = apply_copynumbers(csv_df_snps_mean, csv_df_snps_mean.hp1.values.tolist(), csv_df_snps_mean.hp2.values.tolist(), arguments, snps_cpd_means, snps_cpd_means)
        logging.info('Generating coverage/copy numbers plots genome wide')
        integer_fractional_means = integer_fractional_cluster_means(arguments, df_segs_hp1, df_segs_hp2, centers)
        df_segs_hp1, df_segs_hp2 = adjust_diversified_segments(centers, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, arguments)
        write_copynumber_segments_csv(df_segs_hp1, arguments, centers, integer_fractional_means, 1)
        write_copynumber_segments_csv(df_segs_hp2, arguments, centers, integer_fractional_means, 2)
        copy_number_plots_genome(centers, integer_fractional_means, df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, df_unphased, arguments)

    #SNPs LOH and plots
    plot_snps_frequencies(arguments, csv_df_snps_mean, df_segs_hp1, df_segs_hp2, centers, integer_fractional_means)

    return 0

if __name__ == "__main__":
    main()

#--phaseblock-flipping-enable True
# --phaseblocks-enable True --unphased-reads-coverage-enable True

#--smoothing-enable True
#--pdf-enable True
#--copynumbers-enable True

#--threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots  --phased-vcf /home/rezkuh/gits/data/C1/C1.vcf.gz    --copynumbers-enable True  --without-phasing True   --genome-name C1  --cut-threshold 60 --contigs 1-19

#--snps-freq-vcf-enable True
#--dryrun True --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots  --phased-vcf /home/rezkuh/gits/data/colo357/colo357.vcf.gz    --unphased-reads-coverage-enable True --smoothing-enable True --copynumbers-enable True  --unphased-reads-coverage-enable True --snps-freq-vcf-enable True --phaseblock-flipping-enable True   --genome-name colo357  --cut-threshold 150

#--het-phased-snps-freq-enable
#--dryrun True --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots  --phased-vcf /home/rezkuh/gits/data/1937/1937BL.vcf.gz --smoothing-enable True --copynumbers-enable True  --unphased-reads-coverage-enable True --het-phased-snps-freq-enable True --phaseblock-flipping-enable True  --genome-name 1937 --cut-threshold 250
#--dryrun True --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots  --phased-vcf /home/rezkuh/gits/data/1395/1395BL.vcf.gz --copynumbers-enable True  --unphased-reads-coverage-enable True --het-phased-snps-freq-enable True --phaseblock-flipping-enable True  --genome-name 1395 --cut-threshold 150
#--dryrun True --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots  --phased-vcf /home/rezkuh/gits/data/OT4/ON2.vcf.gz    --copynumbers-enable True --het-phased-snps-freq-enable True   --genome-name OT4  --cut-threshold 60 --phaseblock-flipping-enable True --contigs chr1-38

#--dryrun True --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots  --phased-vcf /home/rezkuh/gits/data/CBL/CBL.vcf.gz --smoothing-enable True --copynumbers-enable True  --unphased-reads-coverage-enable True --snps-freq-vcf-enable True  --genome-name C13 --cut-threshold 1 --without-phasing True --contigs 1-19

#--dryrun True --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots  --phased-vcf /home/rezkuh/gits/data/pancreatic_organoid/pancBL.vcf.gz  --copynumbers-enable True  --unphased-reads-coverage-enable True --snps-freq-vcf-enable True --phaseblock-flipping-enable True --phaseblocks-enable True  --genome-name pancreatic_organoid  --cut-threshold 150

#--dryrun True --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots  --phased-vcf /home/rezkuh/gits/data/mouse/somatic_calls/C14_somatic_calls_pass.vcf.gz    --unphased-reads-coverage-enable True --snps-freq-vcf-enable True   --cut-threshold 0 --genome-name C14 --without-phasing True  --contigs 1-19,X,Y


#--dryrun True --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots --snps-freq-vcf-enable True --cut-threshold 70 --copynumbers-enable True --without-phasing True  --contigs 1-19,X,Y --phased-vcf /home/rezkuh/gits/data/mouse/somatic_calls/C13_somatic_calls_pass_snp.vcf.gz  --genome-name C13
#--dryrun True --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots --snps-freq-vcf-enable True --cut-threshold 50 --copynumbers-enable True --phaseblock-flipping-enable True   --snps-freq-vcf-enable True  --phased-vcf /home/rezkuh/gits/data/R10/HT3/HT3.vcf.gz  --genome-name R10/HT3
#--dryrun True --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots coverage_plots --snps-freq-vcf-enable True --cut-threshold 50 --copynumbers-enable True --phaseblock-flipping-enable True   --snps-freq-vcf-enable True  --phased-vcf /home/rezkuh/gits/data/R10/CaSki/CaSki.vcf.gz  --genome-name R10/CaSki

#--dryrun True --threads 1 --reference /home/rezkuh/GenData/reference/GRCh38_no_alt_analysis_set.fasta  --target-bam /home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam --out-dir-plots HG008_HiFi  --normal-phased-vcf /home/rezkuh/gits/data/HG008_HiFi/HG008BL_HiFi.vcf.gz --tumor-vcf /home/rezkuh/gits/data/HG008_HiFi/HG008_HiFi.vcf.gz --copynumbers-enable True  --unphased-reads-coverage-enable True  --phaseblock-flipping-enable True --phaseblocks-enable True   --genome-name HG008_HiFi  --cut-threshold 150