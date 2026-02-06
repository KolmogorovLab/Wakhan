import pathlib
import subprocess
import statistics
import logging
import os
import gzip
import pysam
import pandas as pd
import numpy as np
from intervaltree import IntervalTree, Interval
from collections import  defaultdict, Counter

logger = logging.getLogger()

from src.utils.chromosome import csv_df_chromosomes_sorter
from src.output.writers import write_segments_coverage_snps
from src.file_tools.process_bam import process_bam_for_snps_freqs
from src.coverage.smoothing import smoothing

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
    #snps_het = []
    #snps_homo = []
    snps_het_pos = []
    snps_homo_pos = []
    for index, (gt,vaf) in enumerate(zip(snps_df_gt,snps_df_vaf)):
        if gt == '1/1' or  gt == '1|1':
        #if vaf > 0.9 or vaf < 0.1:
            #snps_homo.append(vaf)
            snps_homo_pos.append(snps_df_pos[index])
        else:
            #snps_het.append(vaf)
            snps_het_pos.append(snps_df_pos[index])

    """
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
    """

    snps_homo_counts = []
    snps_het_counts = []
    loh_regions = []
    centromere_region = []
    #num_hets = []
    #num_homs = []
    for index, pos in enumerate(ref_start_values):
        l1 = [i for i in snps_het_pos if i > pos and i < pos + bin_size]
        l2 = [i for i in snps_homo_pos if i > pos and i < pos + bin_size]
        #num_hets.append(len(l1))
        #num_homs.append(len(l2))
        # snps_homo_counts.append(len(l2)/len(snps_homo_pos)*bin_size)
        if l1:
            het_ratio = len(l1) / (len(l1) + len(l2))
            snps_het_counts.append(het_ratio)
        else:
            snps_het_counts.append(0)
            #het_ratio = 0
        if l2:
            homo_ratio = len(l2) / (len(l1) + len(l2))
            snps_homo_counts.append(homo_ratio)
        else:
            snps_homo_counts.append(0)
            #homo_ratio = 0

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

    #Smoothing was shifting values by 1 to the left, now disabled
    #if snps_het_counts:
    #    snps_het_counts_updated, snps_homo_counts_updated, _ = \
    #            smoothing(snps_het_counts_updated, snps_homo_counts_updated, snps_homo_counts_updated, conv_window_size=smooth_loh_conv_window)

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

    #Note: removing 1st and last bin from LOH since they only partially contain LOH and may contain phased regions with breakpoints
    region_starts = [v + bin_size for i, v in enumerate(region) if i == 0 or region[i] > region[i - 1] + bin_size]
    region_ends = [v - 1 for i, v in enumerate(region) if i == len(region) - 1 or region[i + 1] - region[i] > bin_size]

    # region_ends = []
    # for i, val in enumerate(region_starts):
    #     region_ends.append(region[region.index(val) - 1] + bin_size - 1)
    # region_ends = sorted(region_ends)

    return region_starts, region_ends

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
        sub_list = []
        try:
            sub_list = haplotype_1_coverage[haplotype_1_position.index(min(haplotype_1_position, key=lambda x:abs(x-i))):haplotype_1_position.index(min(haplotype_1_position, key=lambda x:abs(x-j)))]
        except ValueError:
            logger.info('No Hets pileup found!')
        sub_list = [x for x in sub_list if x != 0]
        if sub_list:

            #snps_haplotype1_mean.append(statistics.mean(sub_list))
            snps_haplotype1_mean.append(statistics.median(sub_list))
        else:
            snps_haplotype1_mean.append(0)

    snps_haplotype2_mean = []
    for index, (i, j) in enumerate(zip(ref_start_values, ref_end_values)):
        sub_list = []
        try:
            sub_list = haplotype_2_coverage[haplotype_2_position.index(min(haplotype_2_position, key=lambda x:abs(x-i))):haplotype_2_position.index(min(haplotype_2_position, key=lambda x:abs(x-j)))]
        except ValueError:
            logger.info('No Hets pileup found!')
        sub_list = [x for x in sub_list if x != 0]
        if sub_list:
            #snps_haplotype2_mean.append(statistics.mean(sub_list))
            snps_haplotype2_mean.append(statistics.median(sub_list))
        else:
            snps_haplotype2_mean.append(0)

    return snps_haplotype1_mean, snps_haplotype2_mean

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

def vcf_parse_to_csv_for_snps(input_vcf, args, output_subdir):
    # pathlib.Path(input_vcf).suffix #extension
    # TODO add output check conditions with all these processes
    basefile = pathlib.Path(input_vcf).stem  # filename without extension
    output_vcf = basefile + '_snps.vcf.gz'
    output_vcf = f"{os.path.join(args.out_dir_plots, output_subdir, output_vcf)}"

    output_csv = basefile + '_snps.csv'
    output_csv = f"{os.path.join(args.out_dir_plots, output_subdir, output_csv)}"

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


def vcf_parse_to_csv_for_het_phased_snps_phasesets(input_vcf, args, output_subdir):
    #pathlib.Path(input_vcf).suffix #extension
    # TODO add output check conditions with all these processes
    basefile = pathlib.Path(input_vcf).stem #filename without extension
    output_vcf = basefile + '_het_phased_snps.vcf.gz'
    output_vcf = f"{os.path.join(args.out_dir_plots, output_subdir, output_vcf)}"

    output_csv = basefile + '_phasesets.csv'
    output_csv = f"{os.path.join(args.out_dir_plots, output_subdir, output_csv)}"

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
    subprocess.check_call(cmd)

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

def count_bases(positions: List[Tuple[str, int, str, str]]) -> List[Tuple[str, int, int, int, int, int]]:
    base_counts = []
    for seq_name, pos, ref_base, reads in positions:
        counts = Counter(reads.replace(",", ref_base).replace(".", ref_base))
        base_counts.append((seq_name, pos, counts['A'] + counts['a'], counts['C'] + counts['c'], counts['G'] + counts['g'], counts['T'] + counts['t']))
    return base_counts

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


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


def rephase_vcf(flip_bins_df, phasesets_df, loh_df, vcf_in, out_vcf, args):
    chr_list = list(set(flip_bins_df['chr']))
    chr_list_loh = []
    if len(loh_df):
        chr_list_loh = list(set(loh_df['chr']))
    flip_bins = defaultdict(IntervalTree)
    new_phasesets = defaultdict(IntervalTree)
    loh_regions = defaultdict(IntervalTree)

    for seq in chr_list:
        flip_counter = defaultdict(int)
        for index, row in flip_bins_df[flip_bins_df['chr'] == seq].iterrows():
            flip_counter[(row['start'], row['end'])] += 1
        for (start, end), counter in flip_counter.items():
            if counter % 2 == 1:    #can flip multiple times, only include bins flipped odd number of times
                flip_bins[seq].add(Interval(start, end))
        for index, row in phasesets_df[phasesets_df['chr'] == seq].iterrows():
            new_phasesets[seq].add(Interval(row['start'], row['end']))

        for index, row in loh_df[loh_df['chr'] == seq].iterrows():
            loh_regions[seq].add(Interval(row['start'], row['end'], row['hp']))

    pysam_verbosity = pysam.set_verbosity(0)
    with pysam.VariantFile(vcf_in, 'r') as vcf_reader:
        original_ps = defaultdict(IntervalTree)
        blocks_pos = defaultdict(list)
        for var in vcf_reader:
            sample = var.samples.keys()[0]
            if var.samples[sample].phased:
                blocks_pos[(var.chrom, var.samples[sample]['PS'])].append(var.pos)
        for (chrom, ps), positions in blocks_pos.items():
            original_ps[chrom].add(Interval(min(positions), max(positions) + 1, ps))

    vcf_reader = pysam.VariantFile(vcf_in,"r")
    vcf_out = pysam.VariantFile(out_vcf, 'w', header=vcf_reader.header)

    for var in vcf_reader:
        sample = var.samples.keys()[0]
        old_ps = list(original_ps[var.chrom][var.pos])
        new_ps = list(new_phasesets[var.chrom][var.pos])

        if var.samples[sample].phased:
            if len(new_ps) > 0:
                var.samples[sample]['PS'] = new_ps[0][0]

            #if original block is small - it has not been phased by Wakhan.
            #it should be then set to unphased
            if old_ps[0][1] - old_ps[0][0] < args.min_phaseblock:
                var.samples[sample]['PS'] = None
                var.samples[sample]['GT'] = (0, 1)
                var.samples[sample].phased = False
                vcf_out.write(var)
                continue

            flip_ovlps = flip_bins[var.chrom][var.pos]
            flip_ovlps_next = flip_bins[var.chrom][var.pos + args.bin_size]
            flip_ovlps_prev = flip_bins[var.chrom][var.pos - args.bin_size]

            #block start and end that does not fully cover a bin, follows the next/previous full bin
            need_flip = False
            if var.pos - old_ps[0][0] < args.bin_size:
                need_flip = len(flip_ovlps_next) > 0
            elif old_ps[0][1] - var.pos > args.bin_size:
                need_flip = len(flip_ovlps) > 0
            else:
                need_flip = len(flip_ovlps_prev) > 0

            if need_flip > 0:
                (a,b) = var.samples[sample]['GT']
                new_gt = (abs(a-1), abs(b-1))
                var.samples[sample]['GT'] = new_gt
                var.samples[sample].phased = True

        #unphased
        elif var.samples[sample]['GT'] == (1,1):
            #check that variant is within LOH, minus flanking regions
            within_loh = False
            loh_ovlps = list(loh_regions[var.chrom][var.pos])
            if len(loh_ovlps) > 0:
                #if min(var.pos - loh_ovlps[0][0], loh_ovlps[0][1] - var.pos) > args.bin_size_snps:
                within_loh = True

            #within LOH and does not belong to original phse blocks:
            if within_loh and len(old_ps) == 0:
                var.samples[sample]['GT'] = (1,0) if loh_ovlps[0][2] == 1 else (0,1)
                var.samples[sample].phased = True
                var.samples[sample]['PS'] = loh_ovlps[0][0]

        vcf_out.write(var)

    vcf_out.close()
    pysam.set_verbosity(pysam_verbosity)


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
