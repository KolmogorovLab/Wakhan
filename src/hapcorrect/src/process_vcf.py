import pathlib
import subprocess
import statistics
import logging
import os
import gzip
import pysam
import bisect
import pandas as pd
import numpy as np
from collections import  defaultdict, Counter

logger = logging.getLogger()

from src.hapcorrect.src.utils import write_segments_coverage_snps, csv_df_chromosomes_sorter
from src.hapcorrect.src.process_bam import process_bam_for_snps_freqs
from src.hapcorrect.src.smoothing import smoothing

import csv
import multiprocessing
from collections import Counter
from typing import List, Tuple

def af_field_selection(vcf_path):
    vcf_file = open(vcf_path, "r") if "gz" not in vcf_path else gzip.open(vcf_path, "rt")
    af_field = 'NULL'
    for line in vcf_file:
        line = line.rstrip("\n")
        # skip header
        if not line.startswith("#"):
            [_,_, _, _, _, _,_,_, format_, sample_] = line.split("\t")
            # searching in INFO for AF/VAF
            if format_.__contains__("AF"):
                try:
                    format_.split(":").index("AF")
                    af_index = format_.split(":").index("AF")
                    #af = float(sample_.split(":")[af_index])
                    af_field = 'AF'
                except ValueError:
                    if format_.__contains__("VAF"):
                        af_index = format_.split(":").index("VAF")
                        #af = float(sample_.split(":")[af_index])
                        af_field = 'VAF'
            break
    vcf_file.close()
    return af_field

def get_snps_counts(snps_df_sorted, chrom, ref_start_values, bin_size):  # TODO This module needs better implementation, currently slow
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
    for index, pos in enumerate(ref_start_values):
        l2 = [i for i in snps_het_pos if i > pos and i < pos+bin_size]
        if l2:
            snps_het_counts.append(len(l2))
            snps_het_pos_bins.append(pos)

    snps_homo_counts = []
    snps_homo_pos_bins = []
    for index, pos in enumerate(ref_start_values):
        l2 = [i for i in snps_homo_pos if i > pos and i < pos+bin_size]
        if l2:
            snps_homo_counts.append(len(l2))
            snps_homo_pos_bins.append(pos)

    return snps_het_counts, snps_homo_counts, snps_het_pos_bins, snps_homo_pos_bins
def get_snps_frquncies(snps_df_sorted, chrom):  # TODO This module needs better implementation, currently slow
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
    for index, (gt, vaf) in enumerate(zip(snps_df_gt, snps_df_vaf)):
        if gt == '1/1' or gt == '1|1':
            snps_homo.append(vaf)
            snps_homo_pos.append(snps_df_pos[index])
        else:
            snps_het.append(vaf)
            snps_het_pos.append(snps_df_pos[index])

    return snps_het, snps_homo, snps_het_pos, snps_homo_pos

def het_homo_snps_gts(snps_df_sorted, chrom, ref_start_values,
                                bin_size):
    snps_df = snps_df_sorted[snps_df_sorted['chr'] == chrom]
    if 'gt' in snps_df.columns:
        snps_df['gt'].astype(str)
    snps_df = snps_df[(snps_df['qual'] > 15)]

    snps_df_haplotype1 = snps_df[(snps_df['gt'] == '0|1') | (snps_df['gt'] == '1|0')]
    snps_df_haplotype1.reindex(snps_df_haplotype1)
    if snps_df.vaf.dtype == object:
        snps_df_haplotype1_vaf = [eval(i) for i in snps_df.vaf.str.split(',').str[0].values.tolist()]
    else:
        snps_df_haplotype1_vaf = snps_df.vaf.values.tolist()
    snps_df_haplotype1_pos = snps_df_haplotype1.pos.values.tolist()

    snps_df_haplotype2 = snps_df[(snps_df['gt'] == '0/0') | (snps_df['gt'] == '1/1') ]
    snps_df_haplotype2.reindex(snps_df_haplotype2)
    if snps_df.vaf.dtype == object:
        snps_df_haplotype2_vaf = [eval(i) for i in snps_df.vaf.str.split(',').str[0].values.tolist()]
    else:
        snps_df_haplotype2_vaf = snps_df.vaf.values.tolist()
    snps_df_haplotype2_pos = snps_df_haplotype2.pos.values.tolist()

    return snps_df_haplotype1_vaf, snps_df_haplotype2_vaf, snps_df_haplotype1_pos, snps_df_haplotype2_pos

def get_snps_frquncies_coverage(snps_df_sorted, chrom, ref_start_values, bin_size, hets_ratio, smooth_loh_conv_window, args): #TODO This module needs better implementation, currently slow

    snps_df = snps_df_sorted[snps_df_sorted['chr'] == chrom]
    if snps_df.empty:
        return [], [], [], [], [], [], []

    ref_start_values = [i for i in range(1, ref_start_values[-1], bin_size)]
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
        #if vaf > 0.9 or vaf < 0.1:
            snps_homo.append(vaf)
            snps_homo_pos.append(snps_df_pos[index])
        else:
            snps_het.append(vaf)
            snps_het_pos.append(snps_df_pos[index])

    snps_het_counts = []
    for index, pos in enumerate(ref_start_values):
        l2 = [i for i in snps_het_pos if i > pos and i < pos+bin_size]
        #snps_het_counts.append(len(l2)/len(snps_het_pos)*bin_size)#((len(l2)/bin_size)*100)
        #snps_het_counts.append(len(l2))
        if l2:
            snps_het_counts.append(snps_het[snps_het_pos.index(max(l2))])
        else:
            snps_het_counts.append(0)

    snps_homo_counts = []
    for index, pos in enumerate(ref_start_values):
        l2 = [i for i in snps_homo_pos if i > pos and i < pos+bin_size]
        #snps_homo_counts.append(len(l2)/len(snps_homo_pos)*bin_size)
        if l2:
            snps_homo_counts.append(snps_homo[snps_homo_pos.index(max(l2))])
        else:
            snps_homo_counts.append(0)

    snps_homo_counts = []
    snps_het_counts = []
    loh_regions = []
    centromere_region = []
    for index, pos in enumerate(ref_start_values):
        l1 = [i for i in snps_het_pos if i > pos and i < pos + bin_size]
        l2 = [i for i in snps_homo_pos if i > pos and i < pos + bin_size]
        # snps_homo_counts.append(len(l2)/len(snps_homo_pos)*bin_size)
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
            centromere_region.append(pos + bin_size)

    snps_het_counts_updated = []
    snps_homo_counts_updated = []
    ref_start_values_updated = []
    for i, (het_counts, homo_counts, start_values) in enumerate(zip(snps_het_counts, snps_homo_counts, ref_start_values)):
        if not (het_counts == 0 and homo_counts == 0):
            snps_het_counts_updated.append(het_counts)
            snps_homo_counts_updated.append(homo_counts)
            ref_start_values_updated.append(start_values)

    if snps_het_counts:
        snps_het_counts_updated, snps_homo_counts_updated, _ = smoothing(snps_het_counts_updated, snps_homo_counts_updated, snps_homo_counts_updated, conv_window_size=smooth_loh_conv_window)

    for index, pos in enumerate(ref_start_values_updated):
        #if snps_het_counts[index] < 0.7 and snps_homo_counts[index] > 0.14:
        if snps_het_counts_updated[index] < hets_ratio : #and snps_homo_counts_updated[index] > 1 - hets_ratio:
            loh_regions.append(pos)
            loh_regions.append(pos + bin_size)

    centromere_region_starts, centromere_region_ends = squash_regions(centromere_region, bin_size)
    loh_region_starts, loh_region_ends = squash_regions(loh_regions, bin_size)

    return ref_start_values_updated, snps_het_counts_updated, snps_homo_counts_updated, centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends

def squash_regions(region, bin_size):
    series = pd.Series(region, dtype="int")
    region = series.drop_duplicates().tolist()
    region_starts = [v for i, v in enumerate(region) if i == 0 or region[i] > region[i - 1] + bin_size]
    region_ends = [v + bin_size - 1 for i, v in enumerate(region) if i == len(region) - 1 or region[i + 1] - region[i] > bin_size]

    # region_ends = []
    # for i, val in enumerate(region_starts):
    #     region_ends.append(region[region.index(val) - 1] + bin_size - 1)
    # region_ends = sorted(region_ends)

    return region_starts, region_ends

def get_vafs_from_normal_phased_vcf(df_snps, chroms):
    df_final = []
    for index, chrom in enumerate(chroms):
        df = df_snps[df_snps['chr'] == chrom]
        vaf = []
        # df = dict(tuple(df_snps.groupby('hp')))
        haplotype_1_position = df.pos.values.tolist()
        haplotype_1_coverage = df.freq_value_b.values.tolist()
        haplotype_2_position = df.pos.values.tolist()
        haplotype_2_coverage = df.freq_value_a.values.tolist()
        # vaf = a/a+b
        vaf = [round(i / (i + j + 0.0000000001), 3) for i, j in zip(haplotype_1_coverage, haplotype_2_coverage)]
        df_final.append(pd.DataFrame(list(zip([df['chr'].iloc[0] for ch in range(len(haplotype_1_position))], haplotype_1_position, vaf)), columns=['chr', 'pos', 'vaf']))

    return pd.concat(df_final)

def snps_frequencies_chrom_mean_phasesets(df_snps, ref_start_values, ref_end_values, chrom, args):
    df = df_snps[df_snps['chr'] == chrom]

    #df = dict(tuple(df_snps.groupby('hp')))
    haplotype_1_position = df.pos.values.tolist()
    haplotype_1_coverage = df.freq_value_b.values.tolist()
    haplotype_2_position = df.pos.values.tolist()
    haplotype_2_coverage = df.freq_value_a.values.tolist()
    #vaf = a/a+b
    #from normal as well
    snps_haplotype1_mean = []
    for index, (i,j) in enumerate(zip(ref_start_values, ref_end_values)):
        sub_list = haplotype_1_coverage[haplotype_1_position.index(min(haplotype_1_position, key=lambda x:abs(x-i))):haplotype_1_position.index(min(haplotype_1_position, key=lambda x:abs(x-j)))]
        sub_list = [x for x in sub_list if x != 0]
        if sub_list:

            #snps_haplotype1_mean.append(statistics.mean(sub_list))
            snps_haplotype1_mean.append(statistics.median(sub_list))
        else:
            snps_haplotype1_mean.append(0)

    snps_haplotype2_mean = []
    for index, (i, j) in enumerate(zip(ref_start_values, ref_end_values)):
        sub_list = haplotype_2_coverage[haplotype_2_position.index(min(haplotype_2_position, key=lambda x:abs(x-i))):haplotype_2_position.index(min(haplotype_2_position, key=lambda x:abs(x-j)))]
        sub_list = [x for x in sub_list if x != 0]
        if sub_list:
            #snps_haplotype2_mean.append(statistics.mean(sub_list))
            snps_haplotype2_mean.append(statistics.median(sub_list))
        else:
            snps_haplotype2_mean.append(0)

    return snps_haplotype1_mean, snps_haplotype2_mean

def remove_outliers_iqr(data, threshold=5):
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    lower_bound = q1 - threshold * iqr
    upper_bound = q3 + threshold * iqr
    outliers_removed = data[(data >= lower_bound) & (data <= upper_bound)]
    return outliers_removed
def snps_frequencies_chrom_mean(df_snps, ref_start_values, chrom, args):
    df = df_snps[df_snps['chr'] == chrom]

    #df = dict(tuple(df_snps.groupby('hp')))
    haplotype_1_position = df.pos.values.tolist()
    haplotype_1_coverage = df.freq_value_b.values.tolist()
    haplotype_2_position = df.pos.values.tolist()
    haplotype_2_coverage = df.freq_value_a.values.tolist()

    snps_haplotype1_mean=[]
    total=0
    for index, i in enumerate(ref_start_values):
        # if chrom == 'chr14' and i == 106260001:
        #     print('here')
        len_cov= len(df[(df.pos >= i) & (df.pos < i + args.bin_size)])
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
    for index, i in enumerate(ref_start_values):
        len_cov= len(df[(df.pos >= i) & (df.pos < i + args.bin_size)])
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

def snps_frequencies_homo_chrom_mean(file, chrom, args, loh_region_starts, loh_region_ends, haplotype_1_values, haplotype_2_values):

    if args.quick_start:
        csv_df_coverage = csv_df_chromosomes_sorter(args.quick_start_coverage_path+'/'+file, ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
    else:
        csv_df_coverage = csv_df_chromosomes_sorter(args.out_dir_plots+'/coverage_data/'+file, ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

    df = csv_df_coverage[csv_df_coverage['chr'] == chrom]

    #df = dict(tuple(df_snps.groupby('hp')))
    haplotype_1_position = df.pos.values.tolist()
    haplotype_1_coverage = df.freq_value_b.values.tolist()
    haplotype_2_position = df.pos.values.tolist()
    haplotype_2_coverage = df.freq_value_a.values.tolist()

    for j, (starts, ends) in enumerate(zip(loh_region_starts, loh_region_ends)):
        if ends - starts > 2000000:
            total = starts // args.bin_size
            for i in range(starts//args.bin_size,ends//args.bin_size):
                #len_cov = len(df[(df.pos >= i) & (df.pos < i + args.bin_size)])
                len_cov = len(df[(df['pos'] >= i*args.bin_size) & (df['pos'] <= (i+1)*args.bin_size)])
                if len_cov == 0:
                    haplotype_1_sub_value = 0
                else:
                    sub_list = haplotype_1_coverage[total:(total + len_cov)]
                    if sub_list:
                        haplotype_1_sub_value = statistics.mean(sub_list)
                    else:
                        haplotype_1_sub_value = 0

                if len_cov == 0:
                    haplotype_2_sub_value = 0
                else:
                    sub_list = haplotype_2_coverage[total:(total + len_cov)]
                    if sub_list:
                        haplotype_2_sub_value = statistics.mean(sub_list)
                    else:
                        haplotype_2_sub_value = 0

                total += len_cov

                haplotype_1_values[i] = haplotype_1_sub_value + haplotype_2_sub_value
                haplotype_2_values[i] = 0

    return haplotype_1_values, haplotype_2_values

def cpd_mean(haplotype1_means, haplotype2_means, ref_values, chrom, args):
    import ruptures as rpt
    import numpy as np

    data = np.array(haplotype1_means, dtype='int') #numpy.clip(haplotype1_means, a_min=1, a_max=300)
    algo = rpt.Pelt(model="rbf").fit(data)
    result = algo.predict(pen=10)
    change_points = [i for i in result if i <= len(data)]

    snps_haplotype1_mean = []
    snps_haplotype1_pos = []
    start = 0
    snps_haplotype1_pos.append(0)
    for index, point in enumerate(change_points):
        sub_list=haplotype1_means[start:point]
        if sub_list:
            snps_haplotype1_mean.append(statistics.mean(sub_list))
        else:
            snps_haplotype1_mean.append(0)
        start = point + 1
        snps_haplotype1_pos.append(point * args.bin_size)
    snps_haplotype1_pos.append(ref_values[-1] * args.bin_size)
    ############################################################
    data = np.array(haplotype2_means, dtype='int') #numpy.clip(haplotype2_means, a_min=1, a_max=300)
    algo = rpt.Pelt(model="rbf").fit(data)
    result = algo.predict(pen=10)
    change_points = [i for i in result if i <= len(data)]

    snps_haplotype2_mean = []
    snps_haplotype2_pos = []
    start = 0
    snps_haplotype2_pos.append(0)
    for index, point in enumerate(change_points):
        sub_list=haplotype2_means[start:point]
        if sub_list:
            snps_haplotype2_mean.append(statistics.mean(sub_list))
        else:
            snps_haplotype2_mean.append(0)
        start = point + 1
        snps_haplotype2_pos.append(point * args.bin_size)
    snps_haplotype2_pos.append(ref_values[-1] * args.bin_size)

    chr = range(len(snps_haplotype1_pos))
    df_cpd_hp1 = pd.DataFrame(list(zip([chrom for ch in chr], snps_haplotype1_pos, snps_haplotype1_mean)), columns=['chr', 'start', 'hp1'])
    df_cpd_hp2 = pd.DataFrame(list(zip([chrom for ch in chr], snps_haplotype2_pos, snps_haplotype2_mean)), columns=['chr', 'start', 'hp2'])

    slices = slice_when(lambda x, y: y - x > 10, sorted(snps_haplotype1_mean + snps_haplotype2_mean))
    return df_cpd_hp1, df_cpd_hp2

def slice_when(predicate, iterable):
  i, x, size = 0, 0, len(iterable)
  while i < size-1:
    if predicate(iterable[i], iterable[i+1]):
      yield iterable[x:i+1]
      x = i + 1
    i += 1
  yield iterable[x:size]


def vcf_parse_to_csv_for_snps(input_vcf, args):
    # pathlib.Path(input_vcf).suffix #extension
    # TODO add output check conditions with all these processes
    basefile = pathlib.Path(input_vcf).stem  # filename without extension
    output_vcf = basefile + '_snps.vcf.gz'
    output_vcf = f"{os.path.join(args.out_dir_plots, 'data_phasing', output_vcf)}"

    output_csv = basefile + '_snps.csv'
    output_csv = f"{os.path.join(args.out_dir_plots, 'data_phasing', output_csv)}"

    # logger.info('bcftools -> Filtering out hetrozygous and phased SNPs and generating a new VCF')
    # # Filter out het, phased SNPs
    # #cmd = ['bcftools', 'view', '--threads', '$(nproc)',  '--phased', '-g', 'het', '--types', 'snps', input_vcf, '-Oz', '-o', output_vcf]
    # cmd = ['bcftools', 'view', '--threads', '$(nproc)',  '-g', 'het', '--types', 'snps', input_vcf, '-Oz', '-o', output_vcf]
    # process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    # process.wait()

    af_field = af_field_selection(input_vcf)
    logger.info('bcftools -> Query for phasesets and GT, DP, VAF feilds by creating a CSV file')
    # bcftools query for phasesets and GT,DP,VAF
    query = '%CHROM\t%POS\t%QUAL\t[%GT]\t[%DP]\t[%'+af_field+']\n'
    cmd = ['bcftools', 'query', '-f', query, '-i', 'FILTER="PASS"', input_vcf, '-o', output_csv]  #
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    return output_csv
def vcf_parse_to_csv_for_het_phased_snps_phasesets(input_vcf, args):
    #pathlib.Path(input_vcf).suffix #extension
    # TODO add output check conditions with all these processes
    basefile = pathlib.Path(input_vcf).stem #filename without extension
    output_vcf = basefile + '_het_phased_snps.vcf.gz'
    output_vcf = f"{os.path.join(args.out_dir_plots, 'data_phasing', output_vcf)}"

    output_csv = basefile + '_phasesets.csv'
    output_csv = f"{os.path.join(args.out_dir_plots, 'data_phasing', output_csv)}"

    # logger.info('bcftools -> Filtering out hetrozygous and phased SNPs and generating a new VCF')
    # # Filter out het, phased SNPs
    # cmd = ['bcftools', 'view', '--threads', '$(nproc)',  '--phased', '-g', 'het', '--types', 'snps', input_vcf, '-Oz', '-o', output_vcf]
    # process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    # process.wait()

    logger.info('bcftools -> Query for phasesets and GT, DP, VAF feilds by creating a CSV file')
    # bcftools query for phasesets and GT,DP,VAF
    cmd = ['bcftools', 'query', '-f',  '%CHROM\t%POS\t[%PS]\n', '-i PS>1', input_vcf, '-o', output_csv] #
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    return output_csv

def get_snp_frequencies_segments(args, target_bam, thread_pool):
    if args.normal_phased_vcf:
        vcf_input = args.normal_phased_vcf
    else:
        vcf_input = args.tumor_phased_vcf

    basefile = pathlib.Path(vcf_input).stem
    output_csv = basefile + '_het_snps.csv'
    output_csv = f"{os.path.join(args.out_dir_plots, 'data_phasing', output_csv)}"

    output_csv_homo = basefile + '_homo_snps.csv'
    output_csv_homo = f"{os.path.join(args.out_dir_plots, 'data_phasing', output_csv_homo)}"

    output_acgts = basefile + '_het_snps_freqs.csv'
    output_acgts = f"{os.path.join(args.out_dir_plots, 'data_phasing', output_acgts)}"

    output_acgts_homo = basefile + '_homo_snps_freqs.csv'
    output_acgts_homo = f"{os.path.join(args.out_dir_plots, 'data_phasing', output_acgts_homo)}"

    logger.info('bcftools -> Query for het SNPs and creating a %s CSV file', output_csv)
    # bcftools query for phasesets and GT,DP,VAF
    cmd = ['bcftools', 'query', '-i', 'GT="het"', '-f',  '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n', vcf_input, '-o', output_csv] #
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    logger.info('SNPs frequency -> CSV to dataframe conversion for heterozygous SNPs')
    dataframe_snps = csv_df_chromosomes_sorter(output_csv, ['chr', 'start', 'ref', 'alt', 'gt'])

    logger.info('SNPs frequency -> Computing SNPs frequency from tumor BAM')

    if args.quick_start:
        output_pileups = args.quick_start_coverage_path + '/'+ 'pileup_SNPs.csv'
    else:
        output_pileups = process_bam_for_snps_freqs(args, thread_pool)

    logger.info('SNPs frequency -> Computing ACGTs frequencies for heterozygous SNPs')
    compute_acgt_frequency(output_pileups, output_acgts, args)
    dataframe_acgt_frequency = csv_df_chromosomes_sorter(output_acgts, ['chr', 'start', 'a', 'c', 'g', 't'], ',')
    dataframe_acgt_frequency = pd.merge(dataframe_snps, dataframe_acgt_frequency, on=['chr', 'start'])
    snp_segments_frequencies = get_snp_segments_frequencies_final(dataframe_acgt_frequency)
    write_segments_coverage_snps(snp_segments_frequencies, 'snps_frequencies.csv', args)

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

        hp_a = 1
        hp_b = 2

        freq_value_a = 0
        freq_value_b = 0

        if gt == '1|0':  # ref freqs
            hp = 2
            freq_value_a = acgts.get(ref)
            freq_value_b = acgts.get(alt)

        elif gt == '0|1':  # alt freqs
            hp = 1
            freq_value_a = acgts.get(alt)
            freq_value_b = acgts.get(ref)

        if freq_value_a == None:
            freq_value_a = 0
        if freq_value_b == None:
            freq_value_b = 0

        #snp_segments_final.append((contig+'\t'+str(pos)+'\t'+ref+'\t'+alt+'\t'+str(ref_value_new)+'\t'+str(alt_value_new)+'\t'+str(hp)))
        snp_segments_final.append((contig + '\t' + str(pos) + '\t' + str(freq_value_a) + '\t'+ str(hp_a) + '\t' + str(freq_value_b) + '\t' + str(hp_b)))

    return snp_segments_final
def bam_pileups_snps(snps_list, target_bam, args):
    basefile = pathlib.Path(target_bam).stem
    output_csv = basefile + '_snps_pileup.csv'
    output_csv = f"{os.path.join(args.out_dir_plots, 'data_phasing', output_csv)}"

    #target_bam = '/home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7_haplotagged.bam'

    cmd = ['samtools', 'mpileup', '-l',  snps_list, '-f', args.reference, target_bam, '-o', output_csv] #
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    return output_csv

def count_bases(positions: List[Tuple[str, int, str, str]]) -> List[Tuple[str, int, int, int, int, int]]:
    base_counts = []
    for seq_name, pos, ref_base, reads in positions:
        counts = Counter(reads.replace(",", ref_base).replace(".", ref_base))
        base_counts.append((seq_name, pos, counts['A'] + counts['a'], counts['C'] + counts['c'], counts['G'] + counts['g'], counts['T'] + counts['t']))
    return base_counts

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

from collections import Counter

def parse_pileup_line(line):
    chrom, pos, ref, depth, bases, qual = line.strip().split("\t")[:6]
    pos = int(pos)
    depth = int(depth)

    # Remove symbols that are not base letters
    i = 0
    clean_bases = ""
    while i < len(bases):
        c = bases[i]
        if c == '^':  # start of a read segment, skip the next character (mapping quality)
            i += 2
        elif c == '$':  # end of a read segment
            i += 1
        elif c in '+-':  # insertion or deletion
            i += 1
            indel_len = ''
            while i < len(bases) and bases[i].isdigit():
                indel_len += bases[i]
                i += 1
            i += int(indel_len)  # skip indel bases
        else:
            clean_bases += c
            i += 1

    # Normalize bases (case-insensitive) and convert to ACGT
    base_map = {'.': ref.upper(), ',': ref.upper(), 'A': 'A', 'a': 'A',
                'C': 'C', 'c': 'C', 'G': 'G', 'g': 'G', 'T': 'T', 't': 'T'}

    converted = [base_map.get(b, 'N') for b in clean_bases if base_map.get(b, 'N') in 'ACGT']
    base_counts = Counter(converted)
    return {
        'chrom': chrom,
        'pos': pos,
        #'ref': ref,
        #'depth': depth,
        'A': base_counts['A'],
        'C': base_counts['C'],
        'G': base_counts['G'],
        'T': base_counts['T']
    }

def compute_acgt_frequency(pileup, snps_frequency, args): #https://www.biostars.org/p/95700/
    with open(pileup) as infile, open(snps_frequency, "w", newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=["chrom", "pos", "A", "C", "G", "T"])
        #writer.writeheader()
        for line in infile:
            row = parse_pileup_line(line)
            writer.writerow(row)

    ##############################################################
    # with open(pileup, 'r') as f:
    #     input_data = f.readlines()
    # base_counts = []
    # positions = []
    # for line in input_data:
    #     elements = line.strip().split('\t')
    #     positions.append((elements[0], int(elements[1]), elements[2], elements[4]))
    # with multiprocessing.Pool(processes=args.threads) as pool:
    #     base_counts = list(pool.imap(count_bases, chunks(positions, 1000)))
    # base_counts = [item for sublist in base_counts for item in sublist]
    # with open(snps_frequency, 'w') as f:
    #     writer = csv.writer(f)
    #     writer.writerows(base_counts)
    ###################################################################
    # cmd = ['sequenza-utils', 'pileup2acgt', '-p', pileup,  '-o', "pileup2acgt.tsv"]
    # process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    # process.wait()
    #
    # # Read the TSV file
    # df = pd.read_csv("pileup2acgt.tsv", sep="\t", header=None, skiprows=1) #names=['chr',	'n_base',	'ref_base',	'read.depth',	'A',	'C',	'G',	'T',	'strand'])
    # # Select only columns 0, 1, 4, 5, 6, 7
    # df_filtered = df.iloc[:, [0, 1, 4, 5, 6, 7]]
    # # Save as a CSV file (comma-separated)
    # df_filtered.to_csv(snps_frequency, index=False, header=False)

def rephase_vcf(df, id_df, loh_df, vcf_in, out_vcf):
    chr_list = list(set(df['chr']))
    chr_list_loh = []
    if len(loh_df):
        chr_list_loh = list(set(loh_df['chr']))
    start_pos = defaultdict(list)
    end_pos = defaultdict(list)
    idls = defaultdict(list)
    start_pos_loh = defaultdict(list)
    end_pos_loh = defaultdict(list)
    hp = defaultdict(list)
    for seq in chr_list:
        start_pos[seq] = sorted([key for key, val in Counter(df.loc[df['chr'] == seq, 'start']).items() if val%2 == 1])
        end_pos[seq] = sorted([key for key, val in Counter(df.loc[df['chr'] == seq, 'end']).items() if val%2 == 1])
        idls[seq] = sorted(id_df.loc[id_df['chr'] == seq, 'start'])
        if seq in chr_list_loh:
            start_pos_loh[seq] = list(loh_df.loc[loh_df['chr'] == seq, 'start'])
            end_pos_loh[seq] = list(loh_df.loc[loh_df['chr'] == seq, 'end'])
            hp[seq] = list(loh_df.loc[loh_df['chr'] == seq, 'hp'])
    vcf_in=pysam.VariantFile(vcf_in,"r")
    vcf_out = pysam.VariantFile(out_vcf, 'w', header=vcf_in.header)
    for var in vcf_in:
        sample = var.samples.keys()[0]
        if var.samples[sample].phased:
            ind = bisect.bisect_right(idls[var.chrom], var.pos)
            if ind > 0:
                var.samples[sample]['PS'] = idls[var.chrom][ind-1]
            strt = bisect.bisect_right(start_pos[var.chrom], var.pos)
            end = bisect.bisect_right(end_pos[var.chrom], var.pos)
            if strt == end + 1:
                (a,b) = var.samples[sample]['GT']
                new_gt = (abs(a-1), abs(b-1))
                var.samples[sample]['GT'] = new_gt
                var.samples[sample].phased = True
        elif var.samples[sample]['GT'] == (1,1):
            strt = bisect.bisect_right(start_pos_loh[var.chrom], var.pos)
            end = bisect.bisect_right(end_pos_loh[var.chrom], var.pos)
            if strt == end + 1:
                var.samples[sample]['GT'] = (1,0) if hp[var.chrom][end] == 1 else (0,1)
                var.samples[sample].phased = True
                var.samples[sample]['PS'] = int(start_pos_loh[var.chrom][end])
        vcf_out.write(var)
    vcf_out.close()

def index_vcf(out_vcf):
    bcf_cmd = ['bcftools', 'index', out_vcf]
    bcf_1 = subprocess.Popen(bcf_cmd, stdout=subprocess.PIPE)
    bcf_1.wait()
    if bcf_1.returncode != 0:
        raise ValueError('bcftols index subprocess returned nonzero value: {}'.format(bcf_1.returncode))
def _calc_nx(lengths, norm_len, rate):
    n50 = 0
    sum_len = 0
    l50 = 0
    for l in sorted(lengths, reverse=True):
        sum_len += l
        l50 += 1
        if sum_len > rate * norm_len:
            n50 = l
            break
    return l50, n50

def get_phasingblocks(hb_vcf):
    MIN_BLOCK_LEN = 10000
    MIN_SNP = 10
    vcf = pysam.VariantFile(hb_vcf)
    haplotype_blocks = defaultdict(list)
    startpoint_list = defaultdict(list)
    endpoint_list = defaultdict(list)
    switch_points = defaultdict(list)
    id_list = defaultdict(list)
    for var in vcf:
        if 'PS' in var.samples.items()[0][1].items()[-1] and var.samples.items()[0][1]['PS']:
            haplotype_blocks[(var.chrom, var.samples.items()[0][1]['PS'])].append(var.pos)
    phased_lengths = []
    for (chr_id, block_name), coords in haplotype_blocks.items():
        if max(coords) - min(coords) > MIN_BLOCK_LEN and len(coords) >= MIN_SNP:
            startpoint_list[chr_id].append(min(coords))
            endpoint_list[chr_id].append(max(coords))
            phased_lengths.append(max(coords) - min(coords))
            id_list[chr_id].append(block_name)
    total_phased = sum(phased_lengths)
    _l50, n50 = _calc_nx(phased_lengths, total_phased, 0.50)
    logger.info("Total phased length: %s", total_phased)
    logger.info("Phase blocks N50: %s", n50)
    return (phased_lengths, haplotype_blocks)