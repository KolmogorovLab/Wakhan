import pathlib
import subprocess
import statistics
import os
import gzip
import pandas as pd
import logging

logger = logging.getLogger()

from src.utils import write_segments_coverage_dict, csv_df_chromosomes_sorter, smoothing
from src.hapcorrect.src.process_bam import process_bam_for_snps_freqs
# Import common VCF utilities from hapcorrect
from src.hapcorrect.src.process_vcf import count_bases, chunks

import csv
import multiprocessing
from collections import Counter
from typing import List, Tuple


def get_snps_counts(snps_df_sorted, chrom, ref_start_values, bin_size):  # TODO This module needs better implementation, currently slow
    snps_df = snps_df_sorted[snps_df_sorted['chr'] == chrom]
    if 'gt' in snps_df.columns:
        snps_df['gt'].astype(str)
    #snps_df = snps_df[(snps_df['qual'] > 15)]

    if snps_df.vaf.dtype == object:
        snps_df_vaf = [eval(i) for i in snps_df.vaf.str.split(',').str[0].values.tolist()]
    else:
        snps_df_vaf = snps_df.vaf.values.tolist()
    snps_df_gt = snps_df['gt'].tolist()

    snps_df_pos = snps_df.pos.values.tolist()
    snps_het = []
    snps_homo = []
    snps_het_pos = []
    snps_homo_pos = []
    for index, (gt,vaf) in enumerate(zip(snps_df_gt,snps_df_vaf)):
        if gt == '1/1' or  gt == '1|1':
            snps_homo.append(vaf)
            snps_homo_pos.append(snps_df_pos[index])
        else:
            snps_het.append(vaf)
            snps_het_pos.append(snps_df_pos[index])

    snps_het_counts = []
    snps_het_pos_bins = []
    for index, pos in enumerate(ref_start_values):
        l2 = [i for i in snps_het_pos if i > pos and i < pos+bin_size]
        snps_het_counts.append(len(l2))
        snps_het_pos_bins.append(pos)

    snps_homo_counts = []
    snps_homo_pos_bins = []
    for index, pos in enumerate(ref_start_values):
        l2 = [i for i in snps_homo_pos if i > pos and i < pos+bin_size]
        snps_homo_counts.append(len(l2))
        snps_homo_pos_bins.append(pos)

    return snps_het_counts, snps_homo_counts, snps_het_pos_bins, snps_homo_pos_bins

def get_snps_counts_cn_regions(snps_df_sorted, chrom, ref_start_values, ref_end_values):  # TODO This module needs better implementation, currently slow
    snps_df = snps_df_sorted[snps_df_sorted['chr'] == chrom]
    if 'gt' in snps_df.columns:
        snps_df['gt'].astype(str)    #snps_df = snps_df[(snps_df['qual'] > 15)]

    if snps_df.vaf.dtype == object:
        snps_df_vaf = [eval(i) for i in snps_df.vaf.str.split(',').str[0].values.tolist()]
    else:
        snps_df_vaf = snps_df.vaf.values.tolist()
    snps_df_gt = snps_df['gt'].tolist()

    snps_df_pos = snps_df.pos.values.tolist()
    snps_het = []
    snps_homo = []
    snps_het_pos = []
    snps_homo_pos = []
    for index, (gt,vaf) in enumerate(zip(snps_df_gt,snps_df_vaf)):
        if gt == '1/1' or  gt == '1|1':
            snps_homo.append(vaf)
            snps_homo_pos.append(snps_df_pos[index])
        else:
            snps_het.append(vaf)
            snps_het_pos.append(snps_df_pos[index])

    snps_het_counts = []
    snps_het_pos_bins = []
    for index, (start,end) in enumerate(zip(ref_start_values, ref_end_values)):
        l2 = [i for i in snps_het_pos if i > start and i < end]
        snps_het_counts.append(len(l2))

    snps_homo_counts = []
    snps_homo_pos_bins = []
    for index, (start,end) in enumerate(zip(ref_start_values,ref_end_values)):
        l2 = [i for i in snps_homo_pos if i > start and i < end]
        snps_homo_counts.append(len(l2))

    return snps_het_counts, snps_homo_counts, ref_start_values, ref_end_values


def get_snps_frquncies_genome(snps_df):  # TODO This module needs better implementation, currently slow
    if 'gt' in snps_df.columns:
        snps_df['gt'].astype(str)
    #snps_df = snps_df[(snps_df['qual'] > 15)]

    if snps_df.vaf.dtype == object:
        snps_df_vaf = [eval(i) for i in snps_df.vaf.str.split(',').str[0].values.tolist()]
    else:
        snps_df_vaf = snps_df.vaf.values.tolist()
    snps_df_gt = snps_df['gt'].tolist()

    snps_df_pos = snps_df.pos.values.tolist()
    snps_het = []
    snps_homo = []
    snps_het_pos = []
    snps_homo_pos = []
    for index, (gt,vaf) in enumerate(zip(snps_df_gt,snps_df_vaf)):
        if gt == '1/1' or  gt == '1|1':
            snps_homo.append(vaf)
            snps_homo_pos.append(snps_df_pos[index])
        else:
            snps_het.append(vaf)
            snps_het_pos.append(snps_df_pos[index])

    return snps_het, snps_homo, snps_het_pos, snps_homo_pos


def het_snps_means_df(args, chrom, snps_df_sorted, ref_start_values):
    snps_df = snps_df_sorted[snps_df_sorted['chr'] == chrom]
    snps_df['gt'].astype(str)
    snps_df = snps_df[(snps_df['qual'] > 15)]

    snps_df_A_allele = snps_df[(snps_df['gt'] == '0|1') | (snps_df['gt'] == '1|0')]
    snps_df_A_allele.reindex(snps_df_A_allele)
    if snps_df_A_allele.vaf.dtype == object:
        snps_df_A_allele_vaf = [eval(i) for i in snps_df_A_allele.vaf.str.split(',').str[0].values.tolist()]
    else:
        snps_df_A_allele_vaf = [i for i in snps_df_A_allele.vaf.values.tolist()]
    snps_df_A_allele_pos = [i for i in snps_df_A_allele.pos.values.tolist()]
    snps_df_A_allele_vaf = [1 - x if x > 0.5 else x for x in snps_df_A_allele_vaf]

    # snps_df_B_allele = snps_df[(snps_df['gt'] == '1|0') | (snps_df['gt'] == '1/0')]
    # snps_df_B_allele.reindex(snps_df_B_allele)
    # if snps_df_B_allele.vaf.dtype == object:
    #     snps_df_B_allele_vaf = [eval(i) for i in snps_df_B_allele.vaf.str.split(',').str[0].values.tolist() if eval(i) > 0.5]
    # else:
    #     snps_df_B_allele_vaf = [i for i in snps_df.vaf.values.tolist() if i > 0.5]

    snps_A_allele_mean = []
    total = 0
    for index, i in enumerate(ref_start_values):
        len_cov = len(snps_df_A_allele[(snps_df_A_allele.pos >= i) & (snps_df_A_allele.pos < i + args.bin_size)])
        if len_cov == 0:
            snps_A_allele_mean.append(0)
        else:
            sub_list = snps_df_A_allele_vaf[total:(total + len_cov)]
            if sub_list:
                snps_A_allele_mean.append(statistics.mean(sub_list))
                #snps_A_allele_mean.append(min(sub_list)/sum(sub_list))
            else:
                snps_A_allele_mean.append(0)
        total += len_cov

    return pd.DataFrame(list(zip([chrom for ch in range(len(ref_start_values))], ref_start_values, snps_A_allele_mean)), columns=['chr', 'start', 'vaf_mean'])

def get_snps_frquncies_coverage(snps_df_sorted, chrom, ref_start_values, args): #TODO This module needs better implementation, currently slow
    snps_df = snps_df_sorted[snps_df_sorted['chr'] == chrom]
    if 'gt' in snps_df.columns:
        snps_df['gt'].astype(str)
    #snps_df = snps_df[(snps_df['qual'] > 15)]

    # snps_df_haplotype1 = snps_df[(snps_df['gt'] == '0|1') | (snps_df['gt'] == '1|0')]
    # snps_df_haplotype1.reindex(snps_df_haplotype1)
    # #snps_df_haplotype1_list = snps_df_haplotype1['vaf'].str.split(',').str[0].values.tolist()
    # snps_df_haplotype1_vaf = snps_df_haplotype1.vaf.values.tolist() #[eval(i) for i in snps_df_haplotype1_list]
    # snps_df_haplotype1_dp = snps_df_haplotype1.dp.values.tolist()
    # snps_df_haplotype1_pos = snps_df_haplotype1.pos.values.tolist()
    # snps_haplotype1 = []
    # snps_haplotype1_pos = []
    # for i in range(0, len(snps_df_haplotype1_vaf)):
    #     snps_haplotype1.append(snps_df_haplotype1_vaf[i] * snps_df_haplotype1_dp[i])
    #     snps_haplotype1_pos.append(snps_df_haplotype1_pos[i])
    #
    # snps_df_haplotype2 = snps_df[(snps_df['gt'] == '0/0') | (snps_df['gt'] == '1/1') ]
    # snps_df_haplotype2.reindex(snps_df_haplotype2)
    # #snps_df_haplotype2_list = snps_df_haplotype2['vaf'].str.split(',').str[0].values.tolist()
    # snps_df_haplotype2_vaf = snps_df_haplotype2.vaf.values.tolist() #[eval(i) for i in snps_df_haplotype2_list]
    # snps_df_haplotype2_dp = snps_df_haplotype2.dp.values.tolist()
    # snps_df_haplotype2_pos = snps_df_haplotype2.pos.values.tolist()
    # snps_haplotype2 = []
    # snps_haplotype2_pos = []
    # for i in range(0, len(snps_df_haplotype2_vaf)):
    #     snps_haplotype2.append(snps_df_haplotype2_vaf[i] * snps_df_haplotype2_dp[i])
    #     snps_haplotype2_pos.append(snps_df_haplotype2_pos[i])
    #
    # snps_haplotype1_mean = []
    # total = 0
    # for index, i in enumerate(ref_start_values):
    #     len_cov = len(snps_df_haplotype1[(snps_df_haplotype1.pos >= i) & (snps_df_haplotype1.pos < i + bin_size)])
    #     if len_cov == 0:
    #         snps_haplotype1_mean.append(0)
    #     else:
    #         sub_list = snps_haplotype1[total:(total + len_cov)]
    #         snps_haplotype1_mean.append(statistics.median(sub_list))
    #     total += len_cov
    #
    # snps_haplotype2_mean = []
    # total = 0
    # for index, i in enumerate(ref_start_values):
    #     len_cov = len(snps_df_haplotype2[(snps_df_haplotype2.pos >= i) & (snps_df_haplotype2.pos < i + bin_size)])
    #     if len_cov == 0:
    #         snps_haplotype2_mean.append(0)
    #     else:
    #         sub_list = snps_haplotype2[total:(total + len_cov)]
    #         snps_haplotype2_mean.append(statistics.median(sub_list))
    #     total += len_cov
    #
    # snps_haplotype1_counts = []
    # total = 0
    # for index, i in enumerate(ref_start_values):
    #     len_cov = len(snps_df_haplotype1[(snps_df_haplotype1.pos >= i) & (snps_df_haplotype1.pos < i + bin_size)])
    #     if len_cov == 0:
    #         snps_haplotype1_counts.append(0)
    #     else:
    #         sub_list = snps_haplotype1[total:(total + len_cov)]
    #         snps_haplotype1_counts.append(len(sub_list))
    #     total += len_cov
    #
    # snps_haplotype2_counts = []
    # total = 0
    # for index, i in enumerate(ref_start_values):
    #     len_cov = len(snps_df_haplotype2[(snps_df_haplotype2.pos >= i) & (snps_df_haplotype2.pos < i + bin_size)])
    #     if len_cov == 0:
    #         snps_haplotype2_counts.append(0)
    #     else:
    #         sub_list = snps_haplotype2[total:(total + len_cov)]
    #         snps_haplotype2_counts.append(len(sub_list))
    #     total += len_cov

    #snps_df['vaf'] = snps_df['vaf'].astype(float)

    if snps_df.vaf.dtype == object:
        snps_df_vaf = [eval(i) for i in snps_df.vaf.str.split(',').str[0].values.tolist()]
    else:
        snps_df_vaf = snps_df.vaf.values.tolist()
    snps_df_gt = snps_df['gt'].tolist()

    snps_df_pos = snps_df.pos.values.tolist()
    snps_het = []
    snps_homo = []
    snps_het_pos = []
    snps_homo_pos = []
    for index, (gt,vaf) in enumerate(zip(snps_df_gt,snps_df_vaf)):
        if gt == '1/1' or  gt == '1|1':
        #if vaf > 0.9 or vaf < 0.1:
            snps_homo.append(vaf)
            snps_homo_pos.append(snps_df_pos[index])
        else:
            snps_het.append(vaf)
            snps_het_pos.append(snps_df_pos[index])

    snps_het_counts = []
    for index, pos in enumerate(ref_start_values):
        l2 = [i for i in snps_het_pos if i > pos and i < pos + args.bin_size_snps]
        #snps_het_counts.append(len(l2)/len(snps_het_pos)*bin_size_snps)#((len(l2)/bin_size_snps)*100)
        #snps_het_counts.append(len(l2))
        if l2:
            snps_het_counts.append(snps_het[snps_het_pos.index(max(l2))])
        else:
            snps_het_counts.append(0)

    snps_homo_counts = []
    for index, pos in enumerate(ref_start_values):
        l2 = [i for i in snps_homo_pos if i > pos and i < pos + args.bin_size_snps]
        #snps_homo_counts.append(len(l2)/len(snps_homo_pos)*bin_size_snps)
        if l2:
            snps_homo_counts.append(snps_homo[snps_homo_pos.index(max(l2))])
        else:
            snps_homo_counts.append(0)

    snps_homo_counts = []
    snps_het_counts = []
    loh_regions = []
    centromere_region = []
    for index, pos in enumerate(ref_start_values):
        l1 = [i for i in snps_het_pos if i > pos and i < pos + args.bin_size_snps]
        l2 = [i for i in snps_homo_pos if i > pos and i < pos + args.bin_size_snps]
        # snps_homo_counts.append(len(l2)/len(snps_homo_pos)*bin_size_snps)
        if l1:
            het_ratio = len(l1) / (len(l1) + len(l2))
            snps_het_counts.append(het_ratio)
        else:
            snps_het_counts.append(0)
            het_ratio = 0
        if l2:
            homo_ratio = len(l2) / (len(l1) + len(l2))
            snps_homo_counts.append(homo_ratio)
        else:
            snps_homo_counts.append(0)
            homo_ratio = 0

        if len(l1) == 0 and len(l2) == 0:
            centromere_region.append(pos)
            centromere_region.append(pos + args.bin_size_snps)

    snps_het_counts_updated = []
    snps_homo_counts_updated = []
    ref_start_values_updated = []
    for i, (het_counts, homo_counts, start_values) in enumerate(zip(snps_het_counts, snps_homo_counts, ref_start_values)):
        if not (het_counts == 0 and homo_counts == 0):
            snps_het_counts_updated.append(het_counts)
            snps_homo_counts_updated.append(homo_counts)
            ref_start_values_updated.append(start_values)

    if snps_het_counts:
        snps_het_counts_updated, snps_homo_counts_updated, _ = smoothing(snps_het_counts_updated, snps_homo_counts_updated, snps_homo_counts_updated, conv_window_size=args.hets_smooth_window) #40, 15

    for index, pos in enumerate(ref_start_values_updated):
        #if snps_het_counts[index] < 0.7 and snps_homo_counts[index] > 0.14:
        if snps_het_counts_updated[index] < args.hets_ratio: #and snps_homo_counts_updated[index] > 1 - args.hets_ratio:
            loh_regions.append(pos)
            loh_regions.append(pos + args.bin_size_snps)

    centromere_region_starts, centromere_region_ends = squash_regions(centromere_region, args.bin_size_snps)
    loh_region_starts, loh_region_ends = squash_regions(loh_regions, args.bin_size_snps)

    return ref_start_values_updated, snps_het_counts_updated, snps_homo_counts_updated, centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends

def squash_regions(region, bin_size):
    region = pd.Series(region, dtype=object).drop_duplicates().tolist()
    region_starts = [v for i, v in enumerate(region) if i == 0 or region[i] > region[i - 1] + bin_size]
    region_ends = []
    for i, val in enumerate(region_starts):
        region_ends.append(region[region.index(val) - 1] + bin_size - 1)
    region_ends = sorted(region_ends)

    return region_starts, region_ends

def snps_mean(df_snps, ref_start_values, ref_end_values, chrom, args):
    df = df_snps[df_snps['chr'] == chrom]

    #df = dict(tuple(df_snps.groupby('hp')))
    haplotype_1_position = df.pos.values.tolist()
    haplotype_1_coverage = df.freq_value_b.values.tolist()
    haplotype_2_position = df.pos.values.tolist()
    haplotype_2_coverage = df.freq_value_a.values.tolist()

    snps_haplotype1_mean=[]
    total=0
    for index, (i,j) in enumerate(zip(ref_start_values,ref_end_values)):
        len_cov= len(df[(df.pos >= i) & (df.pos < i + (j-i))])
        if len_cov ==0:
            snps_haplotype1_mean.append(0)
        else:
            sub_list=haplotype_1_coverage[total:(total + len_cov)]
            if sub_list:
                snps_haplotype1_mean.append(statistics.mean(sub_list))
            else:
                snps_haplotype1_mean.append(0)
        total += len_cov

    snps_haplotype2_mean = []
    total = 0
    for index, (i,j) in enumerate(zip(ref_start_values, ref_end_values)):
        len_cov= len(df[(df.pos >= i) & (df.pos < i + (j-i))])
        if len_cov ==0:
            snps_haplotype2_mean.append(0)
        else:
            sub_list=haplotype_2_coverage[total:(total + len_cov)]
            if sub_list:
                snps_haplotype2_mean.append(statistics.mean(sub_list))
            else:
                snps_haplotype2_mean.append(0)
        total += len_cov
    return snps_haplotype1_mean, snps_haplotype2_mean




def vcf_parse_to_csv_for_het_phased_snps(input_vcf, args):
    #pathlib.Path(input_vcf).suffix #extension
    # TODO add output check conditions with all these processes
    basefile = pathlib.Path(input_vcf).stem #filename without extension
    output_vcf = basefile + '_het_phased_snps_bafs.vcf.gz'
    output_vcf = f"{os.path.join(args.out_dir_plots, 'data', output_vcf)}"

    output_csv = basefile + '_bafs.csv'
    output_csv = f"{os.path.join(args.out_dir_plots, 'data', output_csv)}"

    # logger.info('bcftools -> Filtering out hetrozygous and phased SNPs and generating a new VCF')
    # # Filter out het, phased SNPs
    # cmd = ['bcftools', 'view', '--threads', '$(nproc)',  '-g', 'het', '--types', 'snps', input_vcf, '-Oz', '-o', output_vcf]
    # process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    # process.wait()

    logger.info('bcftools -> Query for phasesets and GT, DP, VAF feilds by creating a CSV file')
    # bcftools query for phasesets and GT,DP,VAF
    cmd = ['bcftools', 'query', '-f',  '%CHROM\t%POS\t[%PS]\n', '-i PS>1', input_vcf, '-o', output_csv] #
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    return output_csv

def get_snp_segments(args, target_bam, thread_pool):
    if not args.phaseblock_flipping_disable:
        normal_vcf = os.path.join(args.out_dir_plots, 'phasing_output', args.genome_name+'.rephased.vcf.gz')
    elif args.normal_phased_vcf:
        normal_vcf = args.normal_phased_vcf
    else:
        normal_vcf = args.tumor_phased_vcf

    basefile = pathlib.Path(normal_vcf).stem
    output_csv = basefile + '_het_snps.csv'
    output_csv = f"{os.path.join(args.out_dir_plots, 'data', output_csv)}"

    output_bed = basefile + '_het_snps.bed'
    output_bed = f"{os.path.join(args.out_dir_plots, 'data', output_bed)}"

    output_acgts = basefile + '_het_snps_freqs.csv'
    output_acgts = f"{os.path.join(args.out_dir_plots, 'data', output_acgts)}"

    # logger.info('bcftools -> Filtering out hetrozygous and phased SNPs and generating a new VCF')
    # # Filter out het, phased SNPs
    # cmd = ['bcftools', 'view', '--threads', '$(nproc)',  '--phased', '-g', 'het', '--types', 'snps', args.phased_vcf'], '-Oz', '-o', args.phased_vcf']]
    # process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    # process.wait()

    logger.info('bcftools -> Query for het SNPs and creating a %s CSV file', output_csv)
    # bcftools query for phasesets and GT,DP,VAF
    cmd = ['bcftools', 'query', '-i', 'GT="het"', '-f',  '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n', normal_vcf, '-o', output_csv] #
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    cmd = ['bcftools', 'query', '-i', 'GT="het"', '-f',  '%CHROM\t%POS\n', normal_vcf, '-o', output_bed] #
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    logger.info('SNPs frequency -> CSV to dataframe conversion')
    dataframe_snps = csv_df_chromosomes_sorter(output_csv, ['chr', 'start', 'ref', 'alt', 'gt'])
    logger.info('SNPs frequency -> Comuting het SNPs frequency from tumor BAM')

    if not args.phaseblock_flipping_disable and not args.quick_start:
        output_pileups = args.out_dir_plots + '/coverage_data/' +  'pileup_SNPs.csv'
    elif args.quick_start:
            output_pileups = args.quick_start_coverage_path + '/' + 'pileup_SNPs.csv'
    else:
        output_pileups = process_bam_for_snps_freqs(args, thread_pool)  # TODO Updated

    compute_acgt_frequency(output_pileups, output_acgts, args)
    dataframe_acgt_frequency = csv_df_chromosomes_sorter(output_acgts, ['chr', 'start', 'a', 'c', 'g', 't'], ',')
    dataframe_acgt_frequency = pd.merge(dataframe_snps, dataframe_acgt_frequency, on=['chr', 'start'])
    snp_segments_frequencies = get_snp_segments_frequencies_final(dataframe_acgt_frequency)
    write_segments_coverage_dict(snp_segments_frequencies, 'snps_frequencies.csv', args)

def get_snp_segments_frequencies_final(dataframe_acgt_frequency):
    snp_segments = dataframe_acgt_frequency.values.tolist()
    snp_segments_final = []
    for i in range(len(snp_segments)):
        contig, pos, ref, alt, gt, a, c, g, t = snp_segments[i]

        acgts = {}
        acgts['A'] = a
        acgts['C'] = c
        acgts['G'] = g
        acgts['T'] = t

        freq_value_a = 0
        freq_value_b = 0

        if gt == '1|0' or gt == '1/0' :  # ref freqs
            freq_value_a = acgts.get(ref)
            freq_value_b = acgts.get(alt)

        elif gt == '0|1' or gt == '0/1':  # alt freqs
            freq_value_a = acgts.get(alt)
            freq_value_b = acgts.get(ref)

        if freq_value_a == None:
            freq_value_a = 0
        if freq_value_b == None:
            freq_value_b = 0

        hp_a = 1
        hp_b = 2

        snp_segments_final.append((contig + '\t' + str(pos) + '\t' + str(freq_value_a) + '\t'+ str(hp_a) + '\t' + str(freq_value_b) + '\t' + str(hp_b)))

    return snp_segments_final



def compute_acgt_frequency(pileup, snps_frequency, args): #https://www.biostars.org/p/95700/
    with open(pileup, 'r') as f:
        input_data = f.readlines()
    base_counts = []
    positions = []
    for line in input_data:
        elements = line.strip().split('\t')
        positions.append((elements[0], int(elements[1]), elements[2], elements[4]))
    with multiprocessing.Pool(processes=args.threads) as pool:
        base_counts = list(pool.imap(count_bases, chunks(positions, 1000)))
    base_counts = [item for sublist in base_counts for item in sublist]
    with open(snps_frequency, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(base_counts)



