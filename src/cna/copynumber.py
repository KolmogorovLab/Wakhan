import numpy as np
import pandas as pd
import scipy.stats as stats
import math
import statistics

from src.breakpoint.breakpoints import sv_vcf_bps_cn_check
from src.utils.chromosome import get_contigs_list, overlap_check
from src.utils.statistics import remove_outliers_iqr


def detect_first_copy_integers_fractional_cluster_means(args, df_segs_hp1, df_segs_hp2, centers):
    haplotype_1_values_copy = df_segs_hp1.state.values.tolist()
    haplotype_1_start_values_copy = df_segs_hp1.start.values.tolist()
    haplotype_1_end_values_copy = df_segs_hp1.end.values.tolist()

    haplotype_2_values_copy = df_segs_hp2.state.values.tolist()
    haplotype_2_start_values_copy = df_segs_hp2.start.values.tolist()
    haplotype_2_end_values_copy = df_segs_hp2.end.values.tolist()

    centers_count = len(centers)
    diff_lists_1 = [[] for i in range(0, centers_count)]
    diff_lists_2 = [[] for i in range(0, centers_count)]

    #TODO for HP2
    for i, (start, end, value) in enumerate(zip(haplotype_1_start_values_copy, haplotype_1_end_values_copy, haplotype_1_values_copy)):
        k = 0
        for i in range(len(centers)):
            if value == centers[i]:
                k = i
        diff_lists_1[k].append(end-start)

    for i, (start, end, value) in enumerate(zip(haplotype_2_start_values_copy, haplotype_2_end_values_copy, haplotype_2_values_copy)):
        k = 0
        for i in range(len(centers)):
            if value == centers[i]:
                k = i
        diff_lists_2[k].append(end-start)

    diff_lists = diff_lists_1 + diff_lists_2

    #add to calculate total segments length of each centers
    sumsup = []
    for i in range(len(centers)):
        sumsup.append(sum(diff_lists[i]))
    #except center 0, get the highest center, suppose it be a first copy
    max_value = max(sumsup)
    if not sumsup.index(max_value) == 0:
        first_copy = centers[sumsup.index(max_value)]
    else:
        first_copy = centers[sumsup.index(sorted(sumsup, reverse=True)[1])]
        first_copy = centers[2] #TODO hardcoded for test
    integer_centers = []
    if args.without_phasing:
        temp_centers = [0] + [first_copy/2] + [first_copy] + [(c+3)*first_copy/2 for c in range(len(centers))]
        for i in range(len(temp_centers)):
            integer_centers.append(min(centers, key=lambda x:abs(x-temp_centers[i])))
        integer_centers = list(np.unique(integer_centers))
        fractional_centers = sorted(list(set(centers) - set(integer_centers)))
    else: #TODO for HP2
        temp_centers = [0] + [first_copy] + [(c + 2) * first_copy for c in range(len(centers))]
        for i in range(len(temp_centers)):
            integer_centers.append(min(centers, key=lambda x: abs(x - temp_centers[i])))
        integer_centers = list(np.unique(integer_centers))
        fractional_centers = sorted(list(set(centers) - set(integer_centers)))

    return integer_centers, fractional_centers


def integer_fractional_cluster_means(args, df_segs_hp1, df_segs_hp2, centers):
    integer_centers, fractional_centers = detect_first_copy_integers_fractional_cluster_means(args, df_segs_hp1, df_segs_hp2, centers)
    integer_fractional_lables = []
    fractional_centers_ = [[]]

    index_list = [i + 1 for i in [i for i in range(len(centers)) if centers[i] in integer_centers]]
    lists_ = list(map(list, np.split(centers, index_list)))

    for n, cen in enumerate(centers):
        if cen in integer_centers:
            if n == 0:
                integer_fractional_lables.append('0')
            else:
                integer_fractional_lables.append(str(integer_centers.index(cen)))
        else:
            #integer_fractional_lables.append(str(round(cen, 2)))
            for m in range(len(integer_centers)-1):
                if integer_centers[m+1] > cen > integer_centers[m]:
                    integer_fractional_lables.append(str( round( ((cen - integer_centers[m]) / (integer_centers[m+1] - integer_centers[m])) * (m+1 - m) + m, 2)))
                    break
    return integer_fractional_lables


def adjust_first_copy_mean(args, centers, df_segs_hp1, df_segs_hp2, df_hp1, df_hp2):
    chroms = get_contigs_list(args.contigs)
    df_segs_hp1_updated = df_segs_hp1.copy()
    df_segs_hp2_updated = df_segs_hp2.copy()

    hp_1_values = []
    hp_2_values = []
    for i in range(len(centers)):
        hp_1_values.append([])
        hp_2_values.append([])

    for index, chrom in enumerate(chroms):
        df_segs_hp_1_updated = df_segs_hp1_updated[df_segs_hp1_updated['chromosome'] == chrom]
        df_segs_hp_2_updated = df_segs_hp2_updated[df_segs_hp2_updated['chromosome'] == chrom]
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
    for i in range(2):
        sample_mean.append(statistics.mean(remove_outliers_iqr(np.array(hp_1_values[i] + hp_2_values[i]))))
        sample_stdev.append(statistics.stdev(remove_outliers_iqr(np.array(hp_1_values[i] + hp_2_values[i]))))

    if args.tumor_purity and args.tumor_ploidy:
        centers = [int(sample_mean[0]*i) for i in range(1, len(centers))]
    else:
        centers = [centers[0]] + [int(sample_mean[1]*i) for i in range(1, len(centers))]

    return centers


def update_subclonal_means_states(centers, subclonals, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2, args):
    chroms = get_contigs_list(args.contigs)
    updated_df_segs_hp_1 = []
    updated_df_segs_hp_2 = []
    centers[0] = int(centers[0])
    hp_1_values = []
    hp_2_values = []
    for i in range(len(centers)):
        hp_1_values.append([])
        hp_2_values.append([])

    for index, chrom in enumerate(chroms):
        df_segs_hp_1_updated = df_segs_hp1_updated[df_segs_hp1_updated['chromosome'] == chrom]
        df_segs_hp_2_updated = df_segs_hp2_updated[df_segs_hp2_updated['chromosome'] == chrom]
        df_hp_1 = df_hp1[df_hp1['chr'] == chrom]
        df_hp_2 = df_hp2[df_hp2['chr'] == chrom]

        df_segs_hp_1_updated_start = df_segs_hp_1_updated.start.values.tolist()
        df_segs_hp_2_updated_start = df_segs_hp_2_updated.start.values.tolist()
        df_segs_hp_1_updated_end = df_segs_hp_1_updated.end.values.tolist()
        df_segs_hp_2_updated_end = df_segs_hp_2_updated.end.values.tolist()
        df_segs_hp_1_updated_state = df_segs_hp_1_updated.state.values.tolist()
        df_segs_hp_2_updated_state = df_segs_hp_2_updated.state.values.tolist()

        if args.without_phasing:
            df_hp_1_val = df_hp_1.coverage.values.tolist()
            df_hp_2_val = df_hp_2.coverage.values.tolist()
        else:
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
        if len(hp_1_values[i] + hp_2_values[i]) >=2:
            sample_mean.append(statistics.mean(remove_outliers_iqr(np.array(hp_1_values[i] + hp_2_values[i]))))
            sample_stdev.append(statistics.stdev(remove_outliers_iqr(np.array(hp_1_values[i] + hp_2_values[i]))))
        else:
            sample_mean.append(centers[i])
            sample_stdev.append(5)

    #sample_mean = centers

    for index, chrom in enumerate(chroms):
        df_segs_hp_1_updated = df_segs_hp1_updated[df_segs_hp1_updated['chromosome'] == chrom]
        df_segs_hp_2_updated = df_segs_hp2_updated[df_segs_hp2_updated['chromosome'] == chrom]
        df_hp_1 = df_hp1[df_hp1['chr'] == chrom]
        df_hp_2 = df_hp2[df_hp2['chr'] == chrom]

        df_segs_hp_1_updated_start = df_segs_hp_1_updated.start.values.tolist()
        df_segs_hp_2_updated_start = df_segs_hp_2_updated.start.values.tolist()
        df_segs_hp_1_updated_end = df_segs_hp_1_updated.end.values.tolist()
        df_segs_hp_2_updated_end = df_segs_hp_2_updated.end.values.tolist()
        df_segs_hp_1_updated_state = df_segs_hp_1_updated.state.values.tolist()
        df_segs_hp_2_updated_state = df_segs_hp_2_updated.state.values.tolist()
        df_segs_hp_1_updated_depth = []#df_segs_hp_1_updated.depth.values.tolist()
        df_segs_hp_2_updated_depth = []#df_segs_hp_2_updated.depth.values.tolist()

        if args.without_phasing:
            df_hp_1_val = df_hp_1.coverage.values.tolist()
            df_hp_2_val = df_hp_2.coverage.values.tolist()
        else:
            df_hp_1_val = df_hp_1.hp1.values.tolist()
            df_hp_2_val = df_hp_2.hp2.values.tolist()

        df_segs_hp_1_updated_p_score = []
        df_segs_hp_2_updated_p_score = []

        df_segs_hp_1_updated_subclonal_check = []
        df_segs_hp_2_updated_subclonal_check = []

        for i, (start,end) in enumerate(zip(df_segs_hp_1_updated_start, df_segs_hp_1_updated_end)):
            if len(df_hp_1_val[start//args.bin_size:end//args.bin_size]) > 0:
                seg_mean = statistics.median(remove_outliers_iqr(np.array(df_hp_1_val[start//args.bin_size:end//args.bin_size])))
            else:
                seg_mean = 0 #statistics.mean(df_hp_1_val[start // args.bin_size:end // args.bin_size])
            sample_mean_init = min(sample_mean, key=lambda x: abs(x - seg_mean))
            index =  sample_mean.index(sample_mean_init)
            z_score =  (seg_mean - sample_mean_init) / 10 #sample_stdev[index]
            p_value = stats.norm.sf(abs(z_score)) * 2
            df_segs_hp_1_updated_p_score.append(round(p_value, 7))
            if p_value < args.confidence_subclonal_score:
                df_segs_hp_1_updated_state[i] = seg_mean #statistics.median(remove_outliers_iqr(np.array(df_hp_1_val[start//args.bin_size:end//args.bin_size])))
                df_segs_hp_2_updated_subclonal_check.append('Y')
            else:
                df_segs_hp_1_updated_state[i] = df_segs_hp_1_updated_state[i]#min(centers, key=lambda x: abs(x - seg_mean)) #statistics.median(remove_outliers_iqr(np.array(df_hp_1_val[start//args.bin_size:end//args.bin_size])))
                df_segs_hp_2_updated_subclonal_check.append('N')
            df_segs_hp_1_updated_depth.append(seg_mean)

        indices = []
        for i, value in enumerate(df_segs_hp_1_updated_state):
            if value not in centers:
                indices.append(i)
        slices = list(slice_lists(lambda x, y: y - x > 1, sorted(indices)))
        for k, sl in enumerate(slices):
            slices_depth = list(slice_lists(lambda x, y: y - x > 2, sorted([df_segs_hp_1_updated_depth[m] for m in sl])))
            if len(sl) > 1 and not len(slices_depth) > 1:
                up_val = statistics.median(remove_outliers_iqr(np.array(df_segs_hp_1_updated_state[sl[0]:sl[-1]])))
                #up_val = statistics.mean(remove_outliers_iqr(np.array(df_hp_1_val[df_segs_hp_1_updated_start[sl[0]]//args.bin_size:df_segs_hp_1_updated_start[sl[-1]]//args.bin_size])))
                for l in sl:
                    df_segs_hp_1_updated_state[l] = up_val
                    df_segs_hp_1_updated_depth[l] = up_val

        for i, (start,end) in enumerate(zip(df_segs_hp_2_updated_start, df_segs_hp_2_updated_end)):
            if len(df_hp_2_val[start//args.bin_size:end//args.bin_size]) > 0:
                seg_mean = statistics.median(remove_outliers_iqr(np.array(df_hp_2_val[start//args.bin_size:end//args.bin_size])))
            else:
                seg_mean = 0 #statistics.mean(df_hp_2_val[start//args.bin_size:end//args.bin_size])
            sample_mean_init = min(sample_mean, key=lambda x: abs(x - seg_mean))
            index =  sample_mean.index(sample_mean_init)
            z_score =  (seg_mean - sample_mean_init) / 10 #sample_stdev[index]
            p_value = stats.norm.sf(abs(z_score)) * 2
            df_segs_hp_2_updated_p_score.append(round(p_value, 7))
            if p_value < args.confidence_subclonal_score:
                df_segs_hp_2_updated_state[i] = seg_mean #statistics.median(remove_outliers_iqr(np.array(df_hp_2_val[start//args.bin_size:end//args.bin_size])))
                df_segs_hp_2_updated_subclonal_check.append('Y')
            else:
                df_segs_hp_2_updated_state[i] = df_segs_hp_2_updated_state[i] #min(centers, key=lambda x: abs(x - seg_mean)) #statistics.median(remove_outliers_iqr(np.array(df_hp_2_val[start//args.bin_size:end//args.bin_size])))
                df_segs_hp_2_updated_subclonal_check.append('N')
            df_segs_hp_2_updated_depth.append(seg_mean)

        indices = []
        for i, value in enumerate(df_segs_hp_2_updated_state):
            if value not in centers:
                indices.append(i)
        slices = list(slice_lists(lambda x, y: y - x > 1, sorted(indices)))
        for k, sl in enumerate(slices):
            slices_depth = list(slice_lists(lambda x, y: y - x > 2, sorted([df_segs_hp_2_updated_depth[m] for m in sl])))
            if len(sl) > 1 and not len(slices_depth) > 1:
                up_val = statistics.median(remove_outliers_iqr(np.array(df_segs_hp_2_updated_state[sl[0]:sl[-1]])))
                #up_val = statistics.mean(remove_outliers_iqr(np.array(df_hp_2_val[df_segs_hp_2_updated_start[sl[0]] // args.bin_size:df_segs_hp_2_updated_start[sl[-1]] // args.bin_size])))
                for l in sl:
                    df_segs_hp_2_updated_state[l] = up_val
                    df_segs_hp_2_updated_depth[l] = up_val

        updated_df_segs_hp_1.append(pd.DataFrame(list(zip(df_segs_hp_1_updated.chromosome.values.tolist(), df_segs_hp_1_updated.start.values.tolist(),
                                  df_segs_hp_1_updated.end.values.tolist(), [round(i, 2) for i in df_segs_hp_1_updated_depth], [int(i) for i in df_segs_hp_1_updated_state], df_segs_hp_1_updated_p_score)),
                         columns=['chromosome', 'start', 'end', 'depth', 'state', 'p_value']))
        updated_df_segs_hp_2.append(pd.DataFrame(list(zip(df_segs_hp_2_updated.chromosome.values.tolist(), df_segs_hp_2_updated.start.values.tolist(),
                                  df_segs_hp_2_updated.end.values.tolist(), [round(i, 2) for i in df_segs_hp_2_updated_depth], [int(i) for i in df_segs_hp_2_updated_state], df_segs_hp_2_updated_p_score)),
                         columns=['chromosome', 'start', 'end', 'depth', 'state', 'p_value']))

    df_segs_hp1 = pd.concat(updated_df_segs_hp_1)
    df_segs_hp2 = pd.concat(updated_df_segs_hp_2)
    df_segs_hp1.loc[df_segs_hp1['state'] < centers[0], 'state'] = centers[0]
    df_segs_hp2.loc[df_segs_hp2['state'] < centers[0], 'state'] = centers[0]

    return df_segs_hp1, df_segs_hp2


def subclonal_values_adjusted(x, centers):
    if x == -3300.0:
        return x
    if x < 0:
        x = -x
    #centers = centers + [centers[-1] + centers[1]]
    for i in range(len(centers) - 1):
        if centers[i + 1] >= x >= centers[i]:
            return round(i + ((x - centers[i]) / (centers[i + 1] - centers[i])), 2)
        elif x < centers[0]:
            return 0
    return len(centers)-1 #round((len(centers) - 1) + ((x - centers[-1]) / ((centers[-1] + (centers[-1] - centers[-2])) - centers[-1])), 2)


def integers_values_adjusted(x, centers):
    if x == -3300.0:
        return x
    if x < 0:
        x = -x
    #centers = centers + [centers[-1] + centers[1]]
    for i in range(len(centers) - 1):
        if centers[i + 1] >= x >= centers[i]:
            return math.ceil(round(i + ((x - centers[i]) / (centers[i + 1] - centers[i])), 2))
        elif x < centers[0]:
            return 0
    return len(centers)-1 #round((len(centers) - 1) + ((x - centers[-1]) / ((centers[-1] + (centers[-1] - centers[-2])) - centers[-1])), 2)


#def mask_df_states(haplotype_df, centers, integer_fractional_means):
#    for i in range(len(integer_fractional_means)):
#        haplotype_df['state'].mask(haplotype_df['state'] == centers[i], integer_fractional_means[i], inplace=True)
#
#    return haplotype_df


def bins_with_copynumber_states(df_segs_hp1, df_segs_hp2, df_hp1, df_hp2, args, centers, integer_fractional_centers):
    chroms = get_contigs_list(args.contigs)
    df_hp1_updated = []
    df_hp2_updated = []

    for index, chrom in enumerate(chroms):
        df_hp1_chrom = df_hp1[df_hp1['chr'] == chrom]
        df_hp2_chrom = df_hp2[df_hp2['chr'] == chrom]
        df_seg_hp1 = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
        df_seg_hp2 = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]

        seg_coverage_hp1 = df_seg_hp1.state.values.tolist()
        seg_coverage_hp2 = df_seg_hp2.state.values.tolist()

        starts = df_hp1_chrom.start.values.tolist()
        ends = df_hp1_chrom.end.values.tolist()

        hp1_chrom_coverage = df_hp1_chrom.hp1.values.tolist()
        hp2_chrom_coverage = df_hp2_chrom.hp2.values.tolist()

        for i, (start_b, end_b) in enumerate(zip(starts, ends)):
            check, index = overlap_check(start_b, end_b, df_seg_hp1.start.values.tolist(), df_seg_hp1.end.values.tolist())
            if check:
                hp1_chrom_coverage[i] = integer_fractional_centers[centers.index(seg_coverage_hp1[index])]
            else:
                hp1_chrom_coverage[i] = 0
            check, index = overlap_check(start_b, end_b, df_seg_hp2.start.values.tolist(), df_seg_hp2.end.values.tolist())
            if check:
                hp2_chrom_coverage[i] = integer_fractional_centers[centers.index(seg_coverage_hp2[index])]
            else:
                hp2_chrom_coverage[i] = 0
        df_hp1_updated.append(pd.DataFrame(list(zip(df_hp1_chrom.chr.values.tolist(), df_hp1_chrom.start.values.tolist(), df_hp1_chrom.start.values.tolist(), \
                     hp1_chrom_coverage)), columns=['chr', 'start', 'end', 'hp1']))

        df_hp2_updated.append(pd.DataFrame(list(zip(df_hp2_chrom.chr.values.tolist(), df_hp2_chrom.start.values.tolist(), df_hp2_chrom.start.values.tolist(), \
                     hp2_chrom_coverage)), columns=['chr', 'start', 'end', 'hp2']))

    return pd.concat(df_hp1_updated), pd.concat(df_hp2_updated)


#def check_adjust_last_cn_states(cens, df_segs_hp1, df_segs_hp2):
#    df_hp1 = df_segs_hp1.copy()
#    df_hp2 = df_segs_hp2.copy()
#    df_hp1.loc[df_hp1['state'] < cens[0], 'state'] = cens[0]
#    df_hp2.loc[df_hp2['state'] < cens[0], 'state'] = cens[0]

#    return df_hp1, df_hp2


def add_confidence_score_cn_segemnts(centers, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2, args):
    chroms = get_contigs_list(args.contigs)
    updated_df_segs_hp_1 = []
    updated_df_segs_hp_2 = []
    #centers[0] = int(centers[0])
    hp_1_values = []
    hp_2_values = []
    for i in range(len(centers)):
        hp_1_values.append([])
        hp_2_values.append([])

    for index, chrom in enumerate(chroms):
        df_segs_hp_1_updated = df_segs_hp1_updated[df_segs_hp1_updated['chromosome'] == chrom]
        df_segs_hp_2_updated = df_segs_hp2_updated[df_segs_hp2_updated['chromosome'] == chrom]
        df_hp_1 = df_hp1[df_hp1['chr'] == chrom]
        df_hp_2 = df_hp2[df_hp2['chr'] == chrom]

        df_segs_hp_1_updated_start = df_segs_hp_1_updated.start.values.tolist()
        df_segs_hp_2_updated_start = df_segs_hp_2_updated.start.values.tolist()
        df_segs_hp_1_updated_end = df_segs_hp_1_updated.end.values.tolist()
        df_segs_hp_2_updated_end = df_segs_hp_2_updated.end.values.tolist()
        df_segs_hp_1_updated_state = df_segs_hp_1_updated.state.values.tolist()
        df_segs_hp_2_updated_state = df_segs_hp_2_updated.state.values.tolist()

        if args.without_phasing:
            df_hp_1_val = df_hp_1.coverage.values.tolist()
            df_hp_2_val = df_hp_2.coverage.values.tolist()
        else:
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
        if len(hp_1_values[i] + hp_2_values[i]) >=2:
            sample_mean.append(statistics.mean(remove_outliers_iqr(np.array(hp_1_values[i] + hp_2_values[i]))))
            sample_stdev.append(statistics.stdev(remove_outliers_iqr(np.array(hp_1_values[i] + hp_2_values[i]))))
        else:
            sample_mean.append(centers[i])
            sample_stdev.append(5)

    #sample_mean = centers

    for index, chrom in enumerate(chroms):
        df_segs_hp_1_updated = df_segs_hp1_updated[df_segs_hp1_updated['chromosome'] == chrom]
        df_segs_hp_2_updated = df_segs_hp2_updated[df_segs_hp2_updated['chromosome'] == chrom]
        df_hp_1 = df_hp1[df_hp1['chr'] == chrom]
        df_hp_2 = df_hp2[df_hp2['chr'] == chrom]

        df_segs_hp_1_updated_start = df_segs_hp_1_updated.start.values.tolist()
        df_segs_hp_2_updated_start = df_segs_hp_2_updated.start.values.tolist()
        df_segs_hp_1_updated_end = df_segs_hp_1_updated.end.values.tolist()
        df_segs_hp_2_updated_end = df_segs_hp_2_updated.end.values.tolist()
        df_segs_hp_1_updated_state = df_segs_hp_1_updated.state.values.tolist()
        df_segs_hp_2_updated_state = df_segs_hp_2_updated.state.values.tolist()
        df_segs_hp_1_updated_depth = []#df_segs_hp_1_updated.depth.values.tolist()
        df_segs_hp_2_updated_depth = []#df_segs_hp_2_updated.depth.values.tolist()

        if args.without_phasing:
            df_hp_1_val = df_hp_1.coverage.values.tolist()
            df_hp_2_val = df_hp_2.coverage.values.tolist()
        else:
            df_hp_1_val = df_hp_1.hp1.values.tolist()
            df_hp_2_val = df_hp_2.hp2.values.tolist()

        df_segs_hp_1_updated_p_score = []
        df_segs_hp_2_updated_p_score = []

        df_segs_hp_1_updated_subclonal_check = []
        df_segs_hp_2_updated_subclonal_check = []

        for i, (start,end) in enumerate(zip(df_segs_hp_1_updated_start, df_segs_hp_1_updated_end)):
            if len(df_hp_1_val[start//args.bin_size:end//args.bin_size]) > 0:
                seg_mean = statistics.median(remove_outliers_iqr(np.array(df_hp_1_val[start//args.bin_size:end//args.bin_size])))
            else:
                seg_mean = 0 #statistics.mean(df_hp_1_val[start // args.bin_size:end // args.bin_size])
            sample_mean_init = min(sample_mean, key=lambda x: abs(x - seg_mean))
            index =  sample_mean.index(sample_mean_init)
            z_score =  (seg_mean - sample_mean_init) / 10 #sample_stdev[index]
            p_value = stats.norm.sf(abs(z_score)) * 2
            df_segs_hp_1_updated_p_score.append(round(p_value, 7))


        for i, (start,end) in enumerate(zip(df_segs_hp_2_updated_start, df_segs_hp_2_updated_end)):
            if len(df_hp_2_val[start//args.bin_size:end//args.bin_size]) > 0:
                seg_mean = statistics.median(remove_outliers_iqr(np.array(df_hp_2_val[start//args.bin_size:end//args.bin_size])))
            else:
                seg_mean = 0 #statistics.mean(df_hp_2_val[start//args.bin_size:end//args.bin_size])
            sample_mean_init = min(sample_mean, key=lambda x: abs(x - seg_mean))
            index =  sample_mean.index(sample_mean_init)
            z_score =  (seg_mean - sample_mean_init) / 10 #sample_stdev[index]
            p_value = stats.norm.sf(abs(z_score)) * 2
            df_segs_hp_2_updated_p_score.append(round(p_value, 7))

        updated_df_segs_hp_1.append(pd.DataFrame(list(zip(df_segs_hp_1_updated.chromosome.values.tolist(), df_segs_hp_1_updated.start.values.tolist(),
                                  df_segs_hp_1_updated.end.values.tolist(), [round(i, 2) for i in df_segs_hp_1_updated.depth.values.tolist()], df_segs_hp_1_updated.state.values.tolist(), df_segs_hp_1_updated_p_score)),
                         columns=['chromosome', 'start', 'end', 'depth', 'state', 'confidence_value']))
        updated_df_segs_hp_2.append(pd.DataFrame(list(zip(df_segs_hp_2_updated.chromosome.values.tolist(), df_segs_hp_2_updated.start.values.tolist(),
                                  df_segs_hp_2_updated.end.values.tolist(), [round(i, 2) for i in df_segs_hp_2_updated.depth.values.tolist()], df_segs_hp_2_updated.state.values.tolist(), df_segs_hp_2_updated_p_score)),
                         columns=['chromosome', 'start', 'end', 'depth', 'state', 'confidence_value']))

    df_segs_hp1 = pd.concat(updated_df_segs_hp_1)
    df_segs_hp2 = pd.concat(updated_df_segs_hp_2)
    #df_segs_hp1.loc[df_segs_hp1['state'] < centers[0], 'state'] = centers[0]
    #df_segs_hp2.loc[df_segs_hp2['state'] < centers[0], 'state'] = centers[0]

    return df_segs_hp1, df_segs_hp2


def slice_lists(predicate, iterable):
  i, x, size = 0, 0, len(iterable)
  while i < size-1:
    if predicate(iterable[i], iterable[i+1]):
      yield iterable[x:i+1]
      x = i + 1
    i += 1
  yield iterable[x:size]

