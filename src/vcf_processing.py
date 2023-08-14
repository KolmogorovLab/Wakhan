import pathlib
import subprocess
import statistics
import logging
import os
import pandas as pd
from bam_processing import get_snps_frequencies, process_bam_for_snps_freqs
from utils import write_segments_coverage, csv_df_chromosomes_sorter_snps_frequency

import csv
import multiprocessing
from collections import Counter
from typing import List, Tuple

def get_snps_frquncies_coverage(snps_df_sorted, chrom, ref_start_values, bin_size): #TODO This module needs better implementation, currently slow
    snps_df = snps_df_sorted[snps_df_sorted['chr'] == chrom]
    snps_df['gt'].astype(str)

    snps_df_haplotype1 = snps_df[(snps_df['gt'] == '0|1') | (snps_df['gt'] == '0|2')]
    snps_df_haplotype1.reindex(snps_df_haplotype1)
    snps_df_haplotype1_list = snps_df_haplotype1['vaf'].str.split(',').str[0].values.tolist()
    snps_df_haplotype1_vaf = [eval(i) for i in snps_df_haplotype1_list]
    snps_df_haplotype1_dp = snps_df_haplotype1.dp.values.tolist()
    snps_haplotype1 = []
    for i in range(0, len(snps_df_haplotype1_vaf)):
        snps_haplotype1.append(snps_df_haplotype1_vaf[i] * snps_df_haplotype1_dp[i])

    snps_df_haplotype2 = snps_df[(snps_df['gt'] == '1|0') | (snps_df['gt'] == '2|0') ]
    snps_df_haplotype2.reindex(snps_df_haplotype2)
    snps_df_haplotype2_list = snps_df_haplotype2['vaf'].str.split(',').str[0].values.tolist()
    snps_df_haplotype2_vaf = [eval(i) for i in snps_df_haplotype2_list]
    snps_df_haplotype2_dp = snps_df_haplotype2.dp.values.tolist()
    snps_haplotype2 = []
    for i in range(0, len(snps_df_haplotype2_vaf)):
        snps_haplotype2.append(snps_df_haplotype2_vaf[i] * snps_df_haplotype2_dp[i])

    snps_haplotype1_mean=[]
    total=0
    for index, i in enumerate(ref_start_values):
        len_cov= len(snps_df_haplotype1[(snps_df_haplotype1.pos >= i) & (snps_df_haplotype1.pos < i + bin_size)])
        if len_cov ==0:
            snps_haplotype1_mean.append(0)
        else:
            sub_list=snps_haplotype1[total:(total + len_cov)]
            snps_haplotype1_mean.append(statistics.median(sub_list))
        total += len_cov

    snps_haplotype2_mean = []
    total = 0
    for index, i in enumerate(ref_start_values):
        len_cov= len(snps_df_haplotype2[(snps_df_haplotype2.pos >= i) & (snps_df_haplotype2.pos < i + bin_size)])
        if len_cov ==0:
            snps_haplotype2_mean.append(0)
        else:
            sub_list=snps_haplotype2[total:(total + len_cov)]
            snps_haplotype2_mean.append(statistics.median(sub_list))
        total += len_cov
    return snps_haplotype1_mean, snps_haplotype2_mean

def snps_mean(df_snps, ref_start_values, chrom):
    df = df_snps[df_snps['chr'] == chrom]

    #df = dict(tuple(df_snps.groupby('hp')))
    haplotype_1_position = df.pos.values.tolist()
    haplotype_1_coverage = df.freq_value_b.values.tolist()
    haplotype_2_position = df.pos.values.tolist()
    haplotype_2_coverage = df.freq_value_a.values.tolist()

    snps_haplotype1_mean=[]
    total=0
    for index, i in enumerate(ref_start_values):
        len_cov= len(df[(df.pos >= i) & (df.pos < i + 50000)])
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
        len_cov= len(df[(df.pos >= i) & (df.pos < i + 50000)])
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

def vcf_parse_to_csv_for_het_phased_snps_phasesets(input_vcf):
    #pathlib.Path(input_vcf).suffix #extension
    # TODO add output check conditions with all these processes
    basefile = pathlib.Path(input_vcf).stem #filename without extension
    output_vcf = basefile + '_het_phased_snps.vcf.gz'
    output_vcf = f"{os.path.join('data', output_vcf)}"

    output_csv = basefile + '_phasesets.csv'
    output_csv = f"{os.path.join('data', output_csv)}"

    logging.info('bcftools -> Filtering out hetrozygous and phased SNPs and generating a new VCF')
    # Filter out het, phased SNPs
    cmd = ['bcftools', 'view', '--threads', '$(nproc)',  '--phased', '-g', 'het', '--types', 'snps', input_vcf, '-Oz', '-o', output_vcf]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    logging.info('bcftools -> Query for phasesets and GT, DP, VAF feilds by creating a CSV file')
    # bcftools query for phasesets and GT,DP,VAF
    cmd = ['bcftools', 'query', '-f',  '%CHROM\t%POS\t%QUAL\t%FILTER\t[%PS]\t[%GT]\t[%DP]\t[%VAF]\n', '-i PS>1', output_vcf, '-o', output_csv] #
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    return output_csv

def get_snp_segments(arguments, target_bam, thread_pool):
    basefile = pathlib.Path(arguments['phased_vcf']).stem
    output_csv = basefile + '_het_snps.csv'
    output_csv = f"{os.path.join('data', output_csv)}"

    output_bed = basefile + '_het_snps.bed'
    output_bed = f"{os.path.join('data', output_bed)}"

    output_acgts = basefile + '_het_snps_freqs.csv'
    output_acgts = f"{os.path.join('data', output_acgts)}"

    # logging.info('bcftools -> Filtering out hetrozygous and phased SNPs and generating a new VCF')
    # # Filter out het, phased SNPs
    # cmd = ['bcftools', 'view', '--threads', '$(nproc)',  '--phased', '-g', 'het', '--types', 'snps', arguments['phased_vcf'], '-Oz', '-o', arguments['phased_vcf']]
    # process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    # process.wait()

    logging.info('bcftools -> Query for het SNPs and creating a %s CSV file', output_csv)
    # bcftools query for phasesets and GT,DP,VAF
    cmd = ['bcftools', 'query', '-i', 'GT="het"', '-f',  '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n', arguments['phased_vcf'], '-o', output_csv] #
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    cmd = ['bcftools', 'query', '-i', 'GT="het"', '-f',  '%CHROM\t%POS\n', arguments['phased_vcf'], '-o', output_bed] #
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    #output_csv = '/home/rezkuh/GenData/COLO829/COLO829_chr7_details.csv'
    #output_bed = '/home/rezkuh/GenData/COLO829/colo829_chr7.csv'

    logging.info('SNPs frequency -> CSV to dataframe conversion')
    dataframe_snps = pd.read_csv(output_csv, sep='\t', names=['chr', 'start', 'ref', 'alt', 'gt'])

    logging.info('SNPs frequency -> Comuting het SNPs frequency from tumor BAM')

    #output_pileups = bam_pileups_snps(output_bed, target_bam, arguments)
    output_pileups = process_bam_for_snps_freqs(arguments, thread_pool) #TODO Updated
    #output_pileups = '/home/rezkuh/gits/data/'+arguments['genome_name']+'/'+arguments['genome_name']+'_SNPs.csv'

    compute_acgt_frequency(output_pileups, output_acgts)
    dataframe_acgt_frequency = csv_df_chromosomes_sorter_snps_frequency(output_acgts)
    dataframe_acgt_frequency = pd.merge(dataframe_snps, dataframe_acgt_frequency, on=['chr', 'start'])
    snp_segments_frequencies = get_snp_segments_frequencies_final(dataframe_acgt_frequency)
    write_segments_coverage(snp_segments_frequencies, 'snps_frequencies.csv')

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

        # hp = 0
        # if gt == '1|0':  # ref freqs
        #     hp = 1
        # elif gt == '0|1':  # alt freqs
        #     hp = 2
        #
        # freq_value = acgts.get(ref)
        ######################################
        # if (hp == 1):
        #     freq_value = acgts.get(alt)
        # elif (hp == 2):
        #     freq_value = acgts.get(ref)
        ######################################
        hp = 0

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
def bam_pileups_snps(snps_list, target_bam, arguments):
    basefile = pathlib.Path(target_bam).stem
    output_csv = basefile + '_snps_pileup.csv'
    output_csv = f"{os.path.join('data', output_csv)}"

    #target_bam = '/home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7_haplotagged.bam'

    cmd = ['samtools', 'mpileup', '-l',  snps_list, '-f', arguments['reference'], target_bam, '-o', output_csv] #
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

def compute_acgt_frequency(pileup, snps_frequency): #https://www.biostars.org/p/95700/
    with open(pileup, 'r') as f:
        input_data = f.readlines()
    base_counts = []
    positions = []
    for line in input_data:
        elements = line.strip().split('\t')
        positions.append((elements[0], int(elements[1]), elements[2], elements[4]))
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        base_counts = list(pool.imap(count_bases, chunks(positions, 1000)))
    base_counts = [item for sublist in base_counts for item in sublist]
    with open(snps_frequency, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(base_counts)



