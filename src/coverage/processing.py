import pandas as pd
import statistics
import logging

logger = logging.getLogger()

from src.coverage.smoothing import smoothing


"""
def get_snps_frquncies_coverage_from_bam(df, chrom):
    df = df[df['chr'] == chrom]
    #df = dict(tuple(df.groupby('hp')))
    haplotype_1_position = df.pos.values.tolist()
    haplotype_1_coverage = df.freq_value_b.values.tolist()
    haplotype_2_position = df.pos.values.tolist()
    haplotype_2_coverage = df.freq_value_a.values.tolist()

    return haplotype_1_position, haplotype_1_coverage, haplotype_2_position, haplotype_2_coverage
"""


def flatten(values):
    return [item for sublist in values for item in sublist]


def flatten_smooth(hp1, hp2, unphased):
    hp1 = flatten(hp1)
    hp2 = flatten(hp2)
    unphased = flatten(unphased)
    unphased, hp1, hp2 = smoothing(unphased, hp1, hp2, conv_window_size=15)

    return hp1, hp2, unphased


def seperate_dfs_coverage(args, df, haplotype_1_values_updated, haplotype_2_values_updated, unphased):
    if args.without_phasing:
        return df[['chr', 'start', 'end', 'coverage']].copy()
    else:
        df_hp1 = df[['chr', 'start','end', 'hp1']].copy()
        df_hp2 = df[['chr', 'start','end', 'hp2']].copy()
        df_unphased = df[['chr', 'start','end', 'hp3']].copy()
        df_hp1['hp1'] = haplotype_1_values_updated
        df_hp2['hp2'] = haplotype_2_values_updated
        df_unphased['hp3'] = unphased
        return df_hp1, df_hp2, df_unphased


def get_vafs_from_normal_phased_vcf(df_snps, df_coverages, chroms, args):
    df_final = []
    logger.info('Generating BAFs')
    for index, chrom in enumerate(chroms):
        df = df_snps[df_snps['chr'] == chrom]
        df_coverage = df_coverages[df_coverages['chr'] == chrom]
        starts_pos = df_coverage.start.values.tolist()
        vaf = []
        # df = dict(tuple(df_snps.groupby('hp')))
        haplotype_1_position = df.pos.values.tolist()
        haplotype_1_coverage = df.freq_value_b.values.tolist()
        haplotype_2_position = df.pos.values.tolist()
        haplotype_2_coverage = df.freq_value_a.values.tolist()
        # vaf = a/a+b
        vaf = [round(min(i,j) / (i + j + 0.0000001), 3) for i, j in zip(haplotype_1_coverage, haplotype_2_coverage)]
        #vaf = list(map(lambda x: 1 - x if x > 0.5 else x, vaf))

        snps_het_vaf = []
        for index, pos in enumerate(starts_pos):
            l2 = [i for i in haplotype_1_position if i > pos and i < pos + args.bin_size]
            if l2:
                snps_het_vaf.append(vaf[haplotype_1_position.index(max(l2))])
            else:
                snps_het_vaf.append(0)

        df_final.append(pd.DataFrame(list(zip([df['chr'].iloc[0] for ch in range(len(starts_pos))], starts_pos, snps_het_vaf)), columns=['chr', 'pos', 'vaf']))

    return pd.concat(df_final)


def get_vafs_from_tumor_phased_vcf(df_snps, df_coverages, chroms, args):
    df_final = []

    for index, chrom in enumerate(chroms):
        df = df_snps[df_snps['chr'] == chrom]
        logger.info('Generating BAFs for %s', chrom)
        df_coverage = df_coverages[df_coverages['chr'] == chrom]
        if df.empty or df_coverage.empty:
            continue
        starts_pos = df_coverage.start.values.tolist()
        vaf = []
        # df = dict(tuple(df_snps.groupby('hp')))

        snps_df_haplotype1 = df[(df['gt'] == '0|1') | (df['gt'] == '0/1') | (df['gt'] == '1|0') | (df['gt'] == '1/0')]
        snps_df_haplotype1.reindex(snps_df_haplotype1)
        if df.vaf.dtype == object:
            haplotype_1_coverage = [eval(i) for i in df.vaf.str.split(',').str[0].values.tolist()]
        else:
            haplotype_1_coverage = df.vaf.values.tolist()
        haplotype_1_position = snps_df_haplotype1.pos.values.tolist()

        # snps_df_haplotype1 = df[(df['gt'] == '1|0') | (df['gt'] == '1/0')]
        # snps_df_haplotype1.reindex(snps_df_haplotype1)
        # if df.vaf.dtype == object:
        #     haplotype_2_coverage = [eval(i) for i in df.vaf.str.split(',').str[0].values.tolist()]
        # else:
        #     haplotype_2_coverage = df.vaf.values.tolist()
        # haplotype_2_position = snps_df_haplotype1.pos.values.tolist()

        # vaf = a/a+b
        #vaf = [round(min(i,j) / (i + j + 0.0000001), 3) for i, j in zip(haplotype_1_coverage, haplotype_2_coverage)]

        haplotype_1_coverage = list(map(lambda x: 1 - x if x > 0.5 else x, haplotype_1_coverage))

        snps_het_vaf = []
        for index, pos in enumerate(starts_pos):
            l2 = [i for i in haplotype_1_position if i > pos and i < pos + args.bin_size]
            #print(haplotype_1_position.index(min(l2)), ':', haplotype_1_position.index(max(l2)))
            #data = haplotype_1_coverage[haplotype_1_position.index(min(l2)):haplotype_1_position.index(max(l2))]
            if len(l2) == 1:
                snps_het_vaf.append(haplotype_1_coverage[haplotype_1_position.index(l2[0])])
            elif len(l2) > 1:
                snps_het_vaf.append(statistics.mean(haplotype_1_coverage[haplotype_1_position.index(min(l2)):haplotype_1_position.index(max(l2))]))
            else:
                snps_het_vaf.append(0)

        df_final.append(pd.DataFrame(list(zip([df['chr'].iloc[0] for ch in range(len(starts_pos))], starts_pos, snps_het_vaf)), columns=['chr', 'pos', 'vaf']))

    return pd.concat(df_final)


def extend_snps_ratios_df(chrom, offset, ref_start_values_updated, snps_het_counts, snps_homo_counts):
    chr_list = [chrom for ch in range(len(ref_start_values_updated))]

    df_snps_ratios_chrom = pd.DataFrame(
        list(zip(chr_list, ref_start_values_updated, snps_het_counts, snps_homo_counts)),
        columns=['chr', 'start', 'hets_ratios', 'homos_ratios'])

    df_snps_ratios_chrom['start_overall'] = df_snps_ratios_chrom['start'].apply(lambda x: x + offset)

    return df_snps_ratios_chrom
