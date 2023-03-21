import pathlib
import subprocess
import statistics
import logging
import os
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
    cmd = ['bcftools', 'view',  '--phased', '-g', 'het', '--types', 'snps', input_vcf, '-Oz', '-o', output_vcf]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    logging.info('tabix -> Index generation for new VCF')
    # Tabix VCF
    cmd = ['tabix', '-p',  'vcf', output_vcf]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    logging.info('bcftools -> Query for phasesets and GT,DP,VAF feilds by creating a CSV file')
    # bcftools query for phasesets and GT,DP,VAF
    cmd = ['bcftools', 'query', '-f',  '%CHROM\t%POS\t%QUAL\t%FILTER\t[%PS]\t[%GT]\t[%DP]\t[%VAF]\n', '-i PS>1', output_vcf, '-o', output_csv] #
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    return output_csv


