import numpy as np
import pandas as pd
import scipy.stats as stats
import os
import shutil
import statistics
import logging

logger = logging.getLogger()

from src.breakpoint.breakpoints import sv_vcf_bps_cn_check
from src.utils_tmp.chromosome import csv_df_chromosomes_sorter, df_chromosomes_sorter, get_contigs_list
from src.utils_tmp.statistics import remove_outliers_iqr, weighted_means

def normal_genome_proportion(p, l, C):
    #C — average coverage in tumor genome
    #p - tumor purity [1-p — normal admixture]
    #l  - tumor ploidy (average number of chromosomes)
    #2 - normal ploidy

    average_contribution_by_normal_genome = 2 * (1 - p) / ((2 * (1 - p) + l * p))
    average_contribution_by_tumor_genome = l * p / ((2 * (1 - p) + l * p))
    combined_normal_coverage_fraction = C * 2 * (1 - p) / ((2 * (1 - p) + l * p))
    per_haplotype = C * 2 * (1 - p) / ((2 * (1 - p) + l * p) * 2)

    return average_contribution_by_normal_genome, average_contribution_by_tumor_genome, combined_normal_coverage_fraction, per_haplotype


def weigted_means_ploidy(args, df_hp_1, df_hp_2, centers, integer_fractional_means):
    df_1 = df_hp_1.copy()
    df_2 = df_hp_2.copy()

    fileDir = os.path.dirname(__file__)  # os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)

    for i in range(len(integer_fractional_means)):
        df_1['state'].mask(df_1['state'] == centers[i], integer_fractional_means[i], inplace=True)

    vals = df_1.state.values.tolist()
    weights = [i - j for i, j in zip(df_1.end.values.tolist(), df_1.start.values.tolist())]

    df = pd.DataFrame({'chromosome': df_1.chromosome.values.tolist(), 'start': df_1.start.values.tolist(),
                       'end': df_1.end.values.tolist(), 'weights': weights, 'state': vals})
    df.to_csv(args.out_dir_plots + '/data/hp1_weights.tsv', sep='\t', index=False)

    for i in range(len(integer_fractional_means)):
        df_2['state'].mask(df_2['state'] == centers[i], integer_fractional_means[i], inplace=True)

    vals = df_2.state.values.tolist()
    weights = [i - j for i, j in zip(df_2.end.values.tolist(), df_2.start.values.tolist())]

    df = pd.DataFrame({'chromosome': df_2.chromosome.values.tolist(), 'start': df_2.start.values.tolist(),
                       'end': df_2.end.values.tolist(), 'weights': weights, 'state': vals})
    df.to_csv(args.out_dir_plots + '/data/hp2_weights.tsv', sep='\t', index=False)
    #########################
    df_1 = pd.read_csv(args.out_dir_plots + '/data/hp1_weights.tsv', sep='\t')
    cent_indices = []
    for index, row in df_1.iterrows():
        for index_cent, row_cent in df_centm.iterrows():
            if row['chromosome'] == row_cent['chr'] and (
                    ((row['start'] - 1) == row_cent['start'] or (row['start']) == row_cent['start']) and row_cent['end'] == row_cent['end']):
                cent_indices.append(index)

    df_1 = df_1.drop(cent_indices)

    vals = df_1.state.values.tolist()
    weights = [i - j for i, j in zip(df_1.end.values.tolist(), df_1.start.values.tolist())]

    hp_1 = weighted_means(vals, weights)
    ########################
    df_2 = pd.read_csv(args.out_dir_plots + '/data/hp2_weights.tsv', sep='\t')
    cent_indices = []
    for index, row in df_2.iterrows():
        for index_cent, row_cent in df_centm.iterrows():
            if row['chromosome'] == row_cent['chr'] and (
                    ((row['start'] - 1) == row_cent['start'] or (row['start']) == row_cent['start']) and row_cent['end'] == row_cent['end']):
                cent_indices.append(index)

    df_2 = df_2.drop(cent_indices)

    vals = df_2.state.values.tolist()
    weights = [i - j for i, j in zip(df_2.end.values.tolist(), df_2.start.values.tolist())]

    hp_2 = weighted_means(vals, weights)

    return hp_1 + hp_2


def average_p_value_genome(args, centers, df_segs_hp1_, df_segs_hp2_, df_hp1, df_hp2):
    df_segs_hp1 = df_segs_hp1_.copy()
    df_segs_hp2 = df_segs_hp2_.copy()

    chroms = get_contigs_list(args.contigs)
    hp_1_values = []
    hp_2_values = []
    for i in range(len(centers)):
        hp_1_values.append([])
        hp_2_values.append([])

    # if args.quick_start:
    #     loh_path = args.quick_start_coverage_path + '/'
    # else:
    #     loh_path = args.out_dir_plots + '/coverage_data/'
    # if args.tumor_phased_vcf and os.path.exists(loh_path + args.genome_name + '_loh_segments.csv'):
    #     df_loh = csv_df_chromosomes_sorter(loh_path + args.genome_name + '_loh_segments.csv', ['chr', 'start', 'end', 'hp'])
    #     indices_loh_hp1 = update_state_with_loh_overlap(df_segs_hp1, df_loh)
    #     indices_loh_hp2 = update_state_with_loh_overlap(df_segs_hp2, df_loh)
    # 
    #     df_segs_hp1 = df_segs_hp1.drop(indices_loh_hp1)
    #     df_segs_hp1 = df_segs_hp1.reset_index(drop=True)
    #
    #     df_segs_hp2 = df_segs_hp2.drop(indices_loh_hp2)
    #     df_segs_hp2 = df_segs_hp2.reset_index(drop=True)

    for index, chrom in enumerate(chroms):
        df_segs_hp_1_updated = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
        df_segs_hp_2_updated = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
        df_hp_1 = df_hp1[df_hp1['chr'] == chrom]
        df_hp_2 = df_hp2[df_hp2['chr'] == chrom]

        df_segs_hp_1_updated_start = df_segs_hp_1_updated.start.values.tolist()
        df_segs_hp_2_updated_start = df_segs_hp_2_updated.start.values.tolist()
        df_segs_hp_1_updated_end = df_segs_hp_1_updated.end.values.tolist()
        df_segs_hp_2_updated_end = df_segs_hp_2_updated.end.values.tolist()
        df_segs_hp_1_updated_state = df_segs_hp_1_updated.state.values.tolist()
        df_segs_hp_2_updated_state = df_segs_hp_2_updated.state.values.tolist()

        df_hp_1_val = df_hp_1.hp1.values.tolist()
        df_hp_2_val = df_hp_2.hp2.values.tolist()

        for j, (start, end) in enumerate(zip(df_segs_hp_1_updated_start, df_segs_hp_1_updated_end)):
            for i in range(len(centers)):
                if df_segs_hp_1_updated_state[j] == centers[i]:
                    hp_1_values[i].extend(df_hp_1_val[start // args.bin_size:end // args.bin_size])

        for j, (start, end) in enumerate(zip(df_segs_hp_2_updated_start, df_segs_hp_2_updated_end)):
            for i in range(len(centers)):
                if df_segs_hp_2_updated_state[j] == centers[i]:
                    hp_2_values[i].extend(df_hp_2_val[start // args.bin_size:end // args.bin_size])

    sample_mean = []
    sample_stdev = []
    for i in range(len(centers)):
        if len(hp_1_values[i] + hp_2_values[i]): #and not args.tumor_phased_vcf:
            sample_mean.append(statistics.median(remove_outliers_iqr(np.array(hp_1_values[i] + hp_2_values[i]))))
            sample_stdev.append(15)#statistics.stdev(remove_outliers_iqr(np.array(hp_1_values[i] + hp_2_values[i]))))
        else:
            sample_mean.append(centers[i])
            sample_stdev.append(5)

    p_value_median = []
    df_segs_hp_1_updated_p_score = []
    df_segs_hp_2_updated_p_score = []

    df_segs_hp_1_updated_weight = []
    df_segs_hp_2_updated_weight = []

    sample_mean = centers

    #diffs = df_segs_hp1.end.values.tolist() - df_segs_hp1.start.values.tolist()

    for index, chrom in enumerate(chroms):

        df_segs_hp_1_updated = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
        df_segs_hp_2_updated = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
        df_hp_1 = df_hp1[df_hp1['chr'] == chrom]
        df_hp_2 = df_hp2[df_hp2['chr'] == chrom]

        df_segs_hp_1_updated_start = df_segs_hp_1_updated.start.values.tolist()
        df_segs_hp_2_updated_start = df_segs_hp_2_updated.start.values.tolist()
        df_segs_hp_1_updated_end = df_segs_hp_1_updated.end.values.tolist()
        df_segs_hp_2_updated_end = df_segs_hp_2_updated.end.values.tolist()

        df_hp_1_val = df_hp_1.hp1.values.tolist()
        df_hp_2_val = df_hp_2.hp2.values.tolist()

        for i, (start,end) in enumerate(zip(df_segs_hp_1_updated_start, df_segs_hp_1_updated_end)):
            if end - start * args.bin_size > 1000000:
                bins = [x for x in df_hp_1_val[start // args.bin_size:end // args.bin_size] if x != 0]
                if bins:
                    seg_mean = statistics.median(bins)
                else:
                    seg_mean = statistics.median(df_hp_1_val[start // args.bin_size:end // args.bin_size])
                sample_mean_init = min(sample_mean, key=lambda x: abs(x - seg_mean))
                index =  sample_mean.index(sample_mean_init)
                z_score =  (seg_mean - sample_mean_init) / 10 #statistics.stdev(remove_outliers_iqr(np.array(df_hp_1_val[start//args.bin_size:end//args.bin_size])))
                p_value = stats.norm.sf(abs(z_score)) * 2
                if p_value > 0:
                    df_segs_hp_1_updated_p_score.append(round(p_value, 7))
                    df_segs_hp_1_updated_weight.append(end-start)

        for i, (start,end) in enumerate(zip(df_segs_hp_2_updated_start, df_segs_hp_2_updated_end)):
            if end-start * args.bin_size > 1000000:
                bins = [x for x in df_hp_2_val[start // args.bin_size:end // args.bin_size] if x != 0]
                if bins:
                    seg_mean = statistics.median(bins)
                else:
                    seg_mean = statistics.median(df_hp_2_val[start // args.bin_size:end // args.bin_size])
                sample_mean_init = min(sample_mean, key=lambda x: abs(x - seg_mean))
                index =  sample_mean.index(sample_mean_init)
                z_score =  (seg_mean - sample_mean_init) / 10 #statistics.stdev(remove_outliers_iqr(np.array(df_hp_2_val[start//args.bin_size:end//args.bin_size])))
                p_value = stats.norm.sf(abs(z_score)) * 2
                if p_value > 0:
                    df_segs_hp_2_updated_p_score.append(round(p_value, 7))
                    df_segs_hp_2_updated_weight.append(end-start)

        p_value_median.append(weighted_means(df_segs_hp_1_updated_p_score + df_segs_hp_2_updated_p_score, weights=df_segs_hp_1_updated_weight + df_segs_hp_2_updated_weight))

    #p_value_median = statistics.mean([weighted_means(df_segs_hp_1_updated_p_score, weights=[i - j for i, j in zip(df_segs_hp1.end.values.tolist(), df_segs_hp1.start.values.tolist())]),
    #                                  weighted_means(df_segs_hp_2_updated_p_score, weights=[i - j for i, j in zip(df_segs_hp2.end.values.tolist(), df_segs_hp2.start.values.tolist())])])

    return statistics.mean(p_value_median)


#https://github.com/jankoslavic/py-tools/tree/master/findpeaks
def find_optimized_normal_peaks(args, data, n, spacing=1, limit=None):
    """Finds peaks in `data` which are of `spacing` width and >=`limit`.

    :param data: values
    :param spacing: minimum spacing to the next peak (should be 1 or more)
    :param limit: peaks should have value greater or equal
    :return:
    """

    ln = data.size
    x = np.zeros(ln+2*spacing)
    x[:spacing] = data[0]-1.e-6
    x[-spacing:] = data[-1]-1.e-6
    x[spacing:spacing+ln] = data
    peak_candidate = np.zeros(ln)
    peak_candidate[:] = True
    for s in range(spacing):
        start = spacing - s - 1
        h_b = x[start : start + ln]  # before
        start = spacing
        h_c = x[start : start + ln]  # central
        start = spacing + s + 1
        h_a = x[start : start + ln]  # after
        peak_candidate = np.logical_and(peak_candidate, np.logical_and(h_c > h_b, h_c > h_a))

    ind = np.argwhere(peak_candidate)
    ind = ind.reshape(ind.size)
    if limit is not None:
        ind = ind[data[ind] > limit]

    # t = np.linspace(0., n, n)
    # plt.plot(t, data)
    # plt.axhline(limit, color='r')
    # plt.plot(t[ind], data[ind], 'ro')
    # plt.title('Peaks: minimum value {limit}, minimum spacing {spacing} points'.format(**{'limit': limit, 'spacing': spacing}))
    # plt.savefig(args.out_dir_plots + '/' + args.genome_name + '_' + "normal_optimized_peak.pdf", format="pdf", bbox_inches="tight")

    if len(ind):
        normals =  [i for i in list(ind)]
    else:
        normals = [0]
        logger.info('No normal peak value detected under confidence [%s] in purity: %s and ploidy: %s, so assuming only first estimated value', limit, args.purity_range, args.ploidy_range)

    return normals


def update_segs_with_normal_optimized(df_hp1, df_hp2, df_segs_hp1, df_segs_hp2, dna_tumor_fraction, args):
    df_segs_hp1_ = df_segs_hp1.copy()
    df_segs_hp2_ = df_segs_hp2.copy()

    df_hp1_ = df_hp1.copy()
    df_hp2_ = df_hp2.copy()

    chroms = get_contigs_list(args.contigs)
    updated_df_segs_hp1 = []
    updated_df_segs_hp2 = []
    updated_df_hp1 = []
    updated_df_hp2 = []

    if args.quick_start:
        loh_path = args.quick_start_coverage_path + '/'
    else:
        loh_path = args.out_dir_plots + '/coverage_data/'
    if args.tumor_phased_vcf and os.path.exists(loh_path + args.genome_name + '_loh_segments.csv'):
        df_loh = csv_df_chromosomes_sorter(loh_path + args.genome_name + '_loh_segments.csv', ['chr', 'start', 'end', 'hp'])

    for index, chrom in enumerate(chroms):
        df_hp1_chrom = df_hp1_[df_hp1_['chr'] == chrom]
        df_hp2_chrom = df_hp2_[df_hp2_['chr'] == chrom]
        df_loh_chrom = df_loh[df_loh['chr'] == chrom]
        haplotype_1_values = df_hp1_chrom.hp1.values.tolist()
        haplotype_2_values = df_hp2_chrom.hp2.values.tolist()

        for j, (starts, ends, hp) in enumerate(zip(df_loh_chrom.start.values.tolist(), df_loh_chrom.end.values.tolist(), df_loh_chrom.hp.values.tolist())):
            for i in range(starts // args.bin_size, ends // args.bin_size):
                if hp == 1:
                   haplotype_1_values[i] = max(haplotype_1_values[i] - dna_tumor_fraction, 0)
                   haplotype_2_values[i] = haplotype_2_values[i] + dna_tumor_fraction
                else:
                   haplotype_1_values[i] = haplotype_1_values[i] + dna_tumor_fraction
                   haplotype_2_values[i] = max(haplotype_2_values[i] - dna_tumor_fraction, 0)

        updated_df_hp1.append(
                pd.DataFrame(list(zip(df_hp1_chrom.chr.values.tolist(), df_hp1_chrom.start.values.tolist(), df_hp1_chrom.end.values.tolist(), haplotype_1_values)),
                             columns=['chr', 'start', 'end', 'hp1']))
        updated_df_hp2.append(
                pd.DataFrame(list(zip(df_hp2_chrom.chr.values.tolist(), df_hp2_chrom.start.values.tolist(), df_hp2_chrom.end.values.tolist(), haplotype_2_values)),
                             columns=['chr', 'start', 'end', 'hp2']))

        df_seg_hp1_chrom = df_segs_hp1_[df_segs_hp1_['chromosome'] == chrom]
        df_seg_hp2_chrom = df_segs_hp2_[df_segs_hp2_['chromosome'] == chrom]

        for idx, s1 in df_seg_hp1_chrom.iterrows():
            sub_list = haplotype_1_values[s1['start'] // args.bin_size:s1['end'] // args.bin_size]
            sub_list = [x for x in sub_list if x > 0.1]
            if len(sub_list):
                median = statistics.median(sub_list)
            else:
                median = 0
            df_seg_hp1_chrom.loc[idx, 'state'] = median

        for idx, s1 in df_seg_hp2_chrom.iterrows():
            sub_list = haplotype_2_values[s1['start'] // args.bin_size:s1['end'] // args.bin_size]
            sub_list = [x for x in sub_list]
            if len(sub_list):
                median = statistics.median(sub_list)
            else:
                median = 0
            df_seg_hp2_chrom.loc[idx, 'state'] = median

        updated_df_segs_hp1.append(df_seg_hp1_chrom)
        updated_df_segs_hp2.append(df_seg_hp2_chrom)

    return pd.concat(updated_df_hp1), pd.concat(updated_df_hp2), pd.concat(updated_df_segs_hp1), pd.concat(updated_df_segs_hp2)


def find_p_values_peaks(p_values):
    def remove_consecutive_duplicates(data):
        if not data:
            return []

        result = [data[0]]
        for i in range(1, len(data)):
            if data[i] != data[i - 1]:
                result.append(data[i])
        return result

    def remove_consecutive_diff_1(data):
        if not data:
            return []

        result = [data[0]]  # Start with the first element
        for i in range(1, len(data)):
            # Only add to result if the difference with the last added value is not 1
            if abs(data[i] - result[-1]) != 1:
                result.append(data[i])
        return result

    def find_peaks(data):
        val = []
        for i in range(1, len(data) - 1):
            if data[i] > data[i - 1] and data[i] > data[i + 1]:
                val.append(data[i])
        return val

    data = remove_consecutive_duplicates(p_values)
    indices = [i for i, val in enumerate(p_values) if val in find_peaks(data)]
    if len(p_values) == 1:
        return [0]
    if p_values[0] > p_values[1]:
        indices = [0] + indices
    indices = remove_consecutive_diff_1(indices)

    return indices


def dna_purity_to_cell_purity(alpha, rho):
    """
    Converts DNA purity to cell purity considering tumor ploidy.

    Parameters:
    alpha (float): DNA purity (proportion between 0 and 1).
    rho (float): Tumor cell ploidy (must be positive).

    Returns:
    float: Cell purity (proportion of tumor cells).

    Raises:
    ValueError: If alpha is not in [0,1] or rho is not positive.
    """
    if not (0 <= alpha <= 1):
        raise ValueError("alpha must be between 0 and 1.")
    if rho <= 0:
        raise ValueError("rho must be positive.")
    denominator = rho + 2 * alpha - alpha * rho
    return (2 * alpha) / denominator


def is_100pct_purity_solution_exist(folder_path):
    """
    Searches a folder for subfolders, parses their names for underscores,
    and returns the name of the subfolder where the value after the underscore
    and before the next underscore is "1.0".

    Args:
        folder_path (str): The path to the folder to search.

    Returns:
        str: The name of the subfolder that matches the criteria, or None if no match is found.
    """
    for subfolder_name in os.listdir(folder_path):
        if os.path.isdir(os.path.join(folder_path, subfolder_name)):
            parts = subfolder_name.split('_')
            if len(parts) > 1 and parts[1] == "1.0":
                return subfolder_name
    return None


def move_100pct_purity_sol(args):
    is_100pct_purity_solution = is_100pct_purity_solution_exist(args.out_dir_plots)
    if is_100pct_purity_solution:
        if os.path.exists(args.out_dir_plots + '/100pct_purity_solution'):
            shutil.rmtree(args.out_dir_plots + '/100pct_purity_solution')
            os.mkdir(args.out_dir_plots + '/100pct_purity_solution')
        else:
            os.mkdir(args.out_dir_plots + '/100pct_purity_solution')
        if os.path.exists(args.out_dir_plots + '/100pct_purity_solution'):
            move_folder(args.out_dir_plots +'/'+is_100pct_purity_solution, args.out_dir_plots + '/100pct_purity_solution')
        else:
            logger.info(f"Error: Source folder '{args.out_dir_plots + '/100pct_purity_solution'}' does not exist.")

