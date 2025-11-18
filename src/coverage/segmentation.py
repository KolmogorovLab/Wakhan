import numpy as np
import pandas as pd
import os
import statistics
import ruptures as rpt
from joblib import Parallel, delayed

from src.breakpoint.breakpoints import sv_vcf_bps_cn_check
from src.coverage.binning import update_bins_with_bps, update_bins_with_bps_new
from src.utils.chromosome import csv_df_chromosomes_sorter, df_chromosomes_sorter, get_contigs_list


def update_state_with_loh_overlap(chrom_df, loh_df):
    # Create a mask for rows where state is 0
    mask = chrom_df['state'] == 0

    # Create a subset and reset index to track original indices
    chrom_subset = chrom_df[mask].reset_index()

    # Merge with loh_df on chromosome and chr
    merged = pd.merge(
        chrom_subset,
        loh_df,
        left_on='chromosome',
        right_on='chr',
        suffixes=('_chrom', '_loh')
    )

    # Calculate overlap condition
    merged['overlap'] = (merged['start_chrom'] < merged['end_loh']) & (merged['end_chrom'] > merged['start_loh'])

    # Get the original indices of rows that have at least one overlap
    overlapping_indices = merged[merged['overlap']]['index'].unique()

    return overlapping_indices


def remove_indices(original_list, indices_to_remove):
    indices_to_remove = sorted(set(indices_to_remove), reverse=True)
    for index in indices_to_remove:
        if index < len(original_list):
            original_list.pop(index)
    return original_list


def change_point_detection_means(args, chrom, bps_ids_all, df_chrom, ref_start_values, all_breakpoints, df_centm_chrom, df_loh_chrom):
    df_means_chr = []
    if args.without_phasing:
        means = df_chrom.coverage.values.tolist()
        snps_haplotype1_mean, snps_len, snps_haplotype1_pos = change_point_detection_algo(args.bin_size, means, ref_start_values, args, all_breakpoints, df_centm_chrom, df_loh_chrom)
        snps_pos_start = []
        snps_pos_end = []
        for i in range(len(snps_haplotype1_pos) - 1):
            snps_pos_start.append(snps_haplotype1_pos[i] if snps_haplotype1_pos[i] < 1 else snps_haplotype1_pos[i] + 1)
            snps_pos_end.append(snps_haplotype1_pos[i + 1])

        chr_list = [df_chrom['chr'].iloc[0] for ch in range(len(snps_haplotype1_mean))]
        df_means_chr = pd.DataFrame(list(zip(chr_list, snps_pos_start, snps_pos_end, snps_haplotype1_mean)),
                                             columns=['chromosome', 'start', 'end', 'state'])

        return snps_haplotype1_mean, snps_len, df_means_chr
    else:
        edges_chr = [bps_ids_all[i] + (bps_ids_all[i + 1] if i + 1 < len(bps_ids_all) else []) for i in range(0, len(bps_ids_all), 2)]
        hp_1_breakpoints = []
        hp_2_breakpoints = []
        for j, bp in enumerate(all_breakpoints):
            for i, (a, b, c, hp1, dv, d, e, f, hp2, dv) in enumerate(edges_chr):
                if chrom == a or chrom == d:
                    if (bp == b and hp1 == 0) or (bp == e and hp1 == 0):
                        hp_1_breakpoints.append(bp)
                        hp_2_breakpoints.append(bp)
                        break
                    elif (bp == b and hp1.split('|')[0] == '1'):
                        hp_1_breakpoints.append(bp)
                        break
                    elif (bp == b and hp1.split('|')[0] == '2'):
                        hp_2_breakpoints.append(bp)
                        break
                    elif (bp == e and hp1.split('|')[1] == '1'):
                        hp_1_breakpoints.append(bp)
                        break
                    elif (bp == e and hp1.split('|')[1] == '2'):
                        hp_2_breakpoints.append(bp)
                        break
                    elif (bp == b and hp1.split('|')[0] == '0'):
                        hp_1_breakpoints.append(bp)
                        hp_2_breakpoints.append(bp)
                        break
                    elif (bp == e and hp1.split('|')[1] == '0'):
                        hp_1_breakpoints.append(bp)
                        hp_2_breakpoints.append(bp)
                        break

        snps_haplotype1_mean, snps_haplotype1_len, snps_haplotype1_pos = change_point_detection_algo(args.bin_size, df_chrom.hp1.values.tolist(), ref_start_values, args, sorted(list(set(hp_1_breakpoints))), df_centm_chrom, df_loh_chrom)
        snps_pos_start = []
        snps_pos_end = []
        for i in range(len(snps_haplotype1_pos) - 1):
            snps_pos_start.append(snps_haplotype1_pos[i] if snps_haplotype1_pos[i] < 1 else snps_haplotype1_pos[i] + 1)
            snps_pos_end.append(snps_haplotype1_pos[i + 1])

        chr_list = [df_chrom['chr'].iloc[0] for ch in range(len(snps_haplotype1_mean))]
        df_means_chr.append(pd.DataFrame(list(zip(chr_list, snps_pos_start, snps_pos_end, snps_haplotype1_mean)),
                                             columns=['chromosome', 'start', 'end', 'state']))

        snps_haplotype2_mean, snps_haplotype2_len, snps_haplotype2_pos = change_point_detection_algo(args.bin_size, df_chrom.hp2.values.tolist(), ref_start_values, args, sorted(list(set(hp_2_breakpoints))), df_centm_chrom, df_loh_chrom)
        snps_pos_start = []
        snps_pos_end = []
        for i in range(len(snps_haplotype2_pos) - 1):
            snps_pos_start.append(snps_haplotype2_pos[i] if snps_haplotype2_pos[i] < 1 else snps_haplotype2_pos[i] + 1)
            snps_pos_end.append(snps_haplotype2_pos[i + 1])

        chr_list = [df_chrom['chr'].iloc[0] for ch in range(len(snps_haplotype2_mean))]
        df_means_chr.append(pd.DataFrame(list(zip(chr_list, snps_pos_start, snps_pos_end, snps_haplotype2_mean)),
                                         columns=['chromosome', 'start', 'end', 'state']))

        # cent and loh indices
        fileDir = os.path.dirname(__file__)  # os.path.dirname(os.path.realpath('__file__'))
        cen_coord = os.path.join(fileDir, args.centromere)
        df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
        df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)

        indices_cent_hp1 = update_state_with_loh_overlap(df_means_chr[0], df_centm)
        indices_cent_hp2 = update_state_with_loh_overlap(df_means_chr[1], df_centm)
        if args.quick_start:
            loh_path = args.quick_start_coverage_path + '/'
        else:
            loh_path = args.out_dir_plots + '/coverage_data/'
        if args.tumor_phased_vcf and os.path.exists(loh_path + args.genome_name + '_loh_segments.csv'):
            df_loh = csv_df_chromosomes_sorter(loh_path + args.genome_name + '_loh_segments.csv', ['chr', 'start', 'end', 'hp'])
            indices_loh_hp1 = update_state_with_loh_overlap(df_means_chr[0], df_loh)
            indices_loh_hp2 = update_state_with_loh_overlap(df_means_chr[1], df_loh)

            # remove cent/loh segment
            snps_haplotype1_mean = remove_indices(snps_haplotype1_mean, list(set(indices_cent_hp1[0] + indices_loh_hp1)))
            snps_haplotype1_len = remove_indices(snps_haplotype1_len, list(set(indices_cent_hp1[0] + indices_loh_hp1)))

            snps_haplotype2_mean = remove_indices(snps_haplotype2_mean, list(set(indices_cent_hp2[0] + indices_loh_hp2)))
            snps_haplotype2_len = remove_indices(snps_haplotype2_len, list(set(indices_cent_hp2[0] + indices_loh_hp2)))
        else:
            snps_haplotype1_mean = remove_indices(snps_haplotype1_mean, indices_cent_hp1)
            snps_haplotype1_len = remove_indices(snps_haplotype1_len, indices_cent_hp1)

            snps_haplotype2_mean = remove_indices(snps_haplotype2_mean, indices_cent_hp2)
            snps_haplotype2_len = remove_indices(snps_haplotype2_len, indices_cent_hp2)

        #remove cent segment
        # chroms = get_contigs_list(args.contigs)
        # if not df_chrom['chr'].iloc[0] == chroms[1]:
        #     del snps_haplotype1_mean[cent_index_hp1]
        #     del snps_haplotype2_mean[cent_index_hp2]
        #     del snps_haplotype1_len[cent_index_hp1]
        #     del snps_haplotype2_len[cent_index_hp2]
        return snps_haplotype1_mean + snps_haplotype2_mean, snps_haplotype1_len + snps_haplotype2_len, df_means_chr


def remove_continuous_elements(lst):
    if not lst:
        return lst  # Return empty list if input is empty

    result = [lst[0]]  # Start with the first element
    for i in range(1, len(lst)):
        if lst[i] - lst[i - 1] != 1:  # Check if the difference is not 1
            result.append(lst[i])  # Add to result if difference is not 1
    return result


def change_point_detection_algo(bin_size, hp_data, ref_start_values, args, breakpoints_coordinates, df_centm_chrom, df_loh_chrom):
    breakpoints_coordinates = remove_more_than_one_bp_in_a_single_bin(bin_size, ref_start_values, breakpoints_coordinates)
    if not df_centm_chrom.empty:
        cents = [df_centm_chrom.start.values.tolist()[0], df_centm_chrom.end.values.tolist()[0]]
    else:
        cents = [0,0]
    if args.cpd_internal_segments or args.change_point_detection_for_cna:
        sub_list_data = np.array(hp_data, dtype='int')
        result_cpd_points = parallel_regions_in_cpd(sub_list_data)
        ###################
        #TODO remove, only for sensitive segmentation
        # model = "l2"
        # algo = rpt.Pelt(model=model).fit(sub_list_data)
        # penalty = 5
        # result_cpd_points = algo.predict(pen=penalty)
        ###################
        result_cpd_points = [i*bin_size for i in result_cpd_points]
        breakpoints_coordinates = sorted(list(set([0] + breakpoints_coordinates + result_cpd_points + [ref_start_values[-1]//args.bin_size])))

    cent_indices = []
    for i, bp in enumerate(breakpoints_coordinates):
        if bp >= cents[0] and bp <= cents[1]:
            cent_indices.append(i)
    if cent_indices:
        for index in sorted(list(set(cent_indices)), reverse=True):
            del breakpoints_coordinates[index]

    if args.tumor_phased_vcf:
        if not df_loh_chrom.empty:
            if args.first_copy_breakpoints_filter > 0:
                loh = [0, 0]
            else:
                loh = df_loh_chrom['start'].tolist() + df_loh_chrom['end'].tolist()
        else:
            loh = [0, 0]
        change_points = sorted(list(set(breakpoints_coordinates + cents + loh + [ref_start_values[-1]])))
    else:
        change_points = sorted(list(set(breakpoints_coordinates + cents + [ref_start_values[-1]])))

    change_points = remove_continuous_elements(change_points)
    change_points = [i for i in change_points if i >= 2]

    snps_haplotype_mean = []
    snps_haplotype_len = []
    snps_haplotype_pos = []
    start = 0
    snps_haplotype_pos.append(0)
    for index, point in enumerate(change_points):
        if args.cpd_internal_segments and index == len(change_points):
            break
        sub_list = hp_data[start//bin_size:point//bin_size]
        if sub_list and not point == cents[1]:
            if args.normal_phased_vcf or args.tumor_phased_vcf:
                sub_list = [x for x in sub_list if x != 0]
            if len(sub_list):
                snps_haplotype_mean.append(statistics.median(sub_list))
            else:
                snps_haplotype_mean.append(0)
        else:
            snps_haplotype_mean.append(0)

        snps_haplotype_len.append(len(sub_list))
        snps_haplotype_pos.append(point)
        start = point

    return snps_haplotype_mean, snps_haplotype_len, snps_haplotype_pos


def parallel_regions_in_cpd(signal):
    # Define parameters
    n_jobs = 8  # Number of parallel jobs
    model = "rbf"  # Change point model
    penalty = 10  # Penalty parameter for segmentation
    jump = 25

    # Function to detect change points in a given segment
    def detect_changes(segment, start_idx):
        algo = rpt.Pelt(model=model, jump=jump).fit(segment)
        result = algo.predict(pen=penalty)
        return [idx + start_idx for idx in result[:-1]]  # Exclude the last index

    # Split data into chunks for parallel processing
    chunk_size = len(signal) // n_jobs
    segments = [(signal[i * chunk_size:(i + 1) * chunk_size], i * chunk_size) for i in range(n_jobs)]
    # Run parallel change point detection
    results = Parallel(n_jobs=n_jobs)(delayed(detect_changes)(seg, start) for seg, start in segments)

    # Merge and sort change point indices
    merged_indices = sorted(set([idx for sublist in results for idx in sublist]))

    return merged_indices


def remove_more_than_one_bp_in_a_single_bin(bin_size, ref_start_values, breakpoints_coordinates):
    updated_breakpoints_coordinates = []
    for index,val in enumerate(ref_start_values):
        all_values = []
        for index_bp, bp in enumerate(breakpoints_coordinates):
            if bp >= val and bp <= val+bin_size:
                all_values.append(bp)
        if len(all_values)>1:
            updated_breakpoints_coordinates = updated_breakpoints_coordinates + [all_values[0]]
        elif len(all_values) == 1:
            updated_breakpoints_coordinates = updated_breakpoints_coordinates + all_values
    return updated_breakpoints_coordinates


def merge_adjacent_regions_cn(segarr, args):
    chroms = get_contigs_list(args.contigs)
    dfs = []
    for index, chrom in enumerate(chroms):
        seg = segarr[segarr['chromosome'] == chrom]
        label_groups = seg['state'].ne(seg['state'].shift()).cumsum()
        df = (seg.groupby(label_groups).agg({'chromosome': 'first', 'start': 'min', 'end': 'max', 'depth': 'median', 'state': 'first'}).reset_index(drop=True))
        dfs.append(df)
    out = pd.concat(dfs)
    return out


def merge_adjacent_regions_loh(segarr, args):
    chroms = get_contigs_list(args.contigs)
    dfs = []
    for index, chrom in enumerate(chroms):
        seg = segarr[segarr['chr'] == chrom]
        if len(seg)>1:
            seg['diff'] = seg['start'] - seg['end'].shift().fillna(seg['start'].iloc[0])
            seg['group'] = (seg['diff'] != 1).cumsum()
            df = (seg.groupby(seg['group']).agg({'chr': 'first', 'start': 'min', 'end': 'max'}).reset_index(drop=True))
            dfs.append(df)
    if len(dfs):
        return pd.concat(dfs)
    else:
        return segarr


def merge_adjacent_regions_cn_unphased(segarr, args):
    chroms = get_contigs_list(args.contigs)
    dfs = []
    for index, chrom in enumerate(chroms):
        seg = segarr[segarr['chromosome'] == chrom]
        label_groups = seg['state'].ne(seg['state'].shift()).cumsum()
        df = (seg.groupby(label_groups).agg({'chromosome': 'first', 'start': 'min', 'end': 'max', 'depth': 'first', 'state': 'first'}).reset_index(drop=True))
        dfs.append(df)
    out = pd.concat(dfs)
    return out


def adjust_diversified_segments(centers, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, args):
    chroms = get_contigs_list(args.contigs)
    updated_df_segs_hp1 = []
    updated_df_segs_hp2 = []
    for index, chrom in enumerate(chroms):
        df_chrom_segs_hp1 = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
        haplotype_1_values_copyrnumbers = df_chrom_segs_hp1.state.values.tolist()
        haplotype_1_start_values_copyrnumbers = df_chrom_segs_hp1.start.values.tolist()
        haplotype_1_end_values_copyrnumbers = df_chrom_segs_hp1.end.values.tolist()

        df_chrom_segs_hp2 = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
        haplotype_2_values_copyrnumbers = df_chrom_segs_hp2.state.values.tolist()
        haplotype_2_start_values_copyrnumbers = df_chrom_segs_hp2.start.values.tolist()
        haplotype_2_end_values_copyrnumbers = df_chrom_segs_hp2.end.values.tolist()

        for i, (start,end) in enumerate(zip(haplotype_1_start_values_copyrnumbers, haplotype_1_end_values_copyrnumbers)):
            haplotype_1_values_copyrnumbers[i] = min(centers, key=lambda x:abs(x - haplotype_1_values_copyrnumbers[i]))

        updated_df_segs_hp1.append(pd.DataFrame(list(zip(df_chrom_segs_hp1.chromosome.values.tolist(), df_chrom_segs_hp1.start.values.tolist(),
                df_chrom_segs_hp1.end.values.tolist(),  df_chrom_segs_hp1.state.values.tolist(), haplotype_1_values_copyrnumbers)), columns=['chromosome', 'start', 'end', 'depth', 'state']))

        for i, (start, end) in enumerate(zip(haplotype_2_start_values_copyrnumbers, haplotype_2_end_values_copyrnumbers)):
            haplotype_2_values_copyrnumbers[i] = min(centers, key=lambda x: abs(x - haplotype_2_values_copyrnumbers[i]))

        updated_df_segs_hp2.append(pd.DataFrame(list(zip(df_chrom_segs_hp2.chromosome.values.tolist(), df_chrom_segs_hp2.start.values.tolist(),
                df_chrom_segs_hp2.end.values.tolist(), df_chrom_segs_hp2.state.values.tolist(), haplotype_2_values_copyrnumbers)), columns=['chromosome', 'start', 'end', 'depth', 'state']))

    #return merge_adjacent_regions_cn(pd.concat(updated_df_segs_hp1), args), merge_adjacent_regions_cn(pd.concat(updated_df_segs_hp2), args)
    return pd.concat(updated_df_segs_hp1), pd.concat(updated_df_segs_hp2)


def remove_overlapping_segments(starts, ends):

  if len(starts) != len(ends):
    raise ValueError("Starts and ends lists must have the same length.")

  filtered_starts = []
  filtered_ends = []
  stack = []

  for i in range(len(starts)):
    start, end = starts[i], ends[i]

    # If the stack is empty or the current segment starts after the previous segment ends,
    # push it onto the stack
    if not stack or stack[-1][1] < start:
      stack.append((start, end))

  # Extract filtered starts and ends from the stack
  filtered_starts, filtered_ends = zip(*stack)

  return list(filtered_starts), list(filtered_ends)


def adjust_bps_cn_segments_boundries(args, haplotype_df):
    haplotype_df = haplotype_df.reset_index(drop=True)
    chroms = get_contigs_list(args.contigs)
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)
    cent_starts = df_centm.start.values.tolist()
    cent_ends = df_centm.end.values.tolist()
    cent_chr = df_centm.chr.values.tolist()

    if args.breakpoints:
        _, _, _, bps_ids_all, bps, bps_bnd = sv_vcf_bps_cn_check(args.breakpoints, args)

        for index, row in haplotype_df.iterrows():
            bps_ids = []
            starts = []
            for bp in bps_ids_all:
                if bp[0] == row['chromosome'] and int(bp[1]) >= int(row['start']) and int(bp[1]) <= int(row['end']):
                    bps_ids.append(bp[2])
                    if int(bp[1]) >= int(row['start']) - args.bin_size and int(bp[1]) <= int(row['end']) + args.bin_size:
                        starts.append(int(bp[1]))
            #bps_ids_global.append(bps_ids)
            if len(starts):
                if not row['chromosome'] in cent_chr:
                    continue
                cent_chr_index = cent_chr.index(row['chromosome'])
                if not haplotype_df.loc[index, 'start'] - 1 == cent_starts[cent_chr_index] and not haplotype_df.loc[index, 'end'] == cent_ends[cent_chr_index]:
                    haplotype_df.loc[index, 'end'] = min(starts, key=lambda x: abs(x - int(row['start']))) - 1
        #haplotype_df['svs_breakpoints_ids'] = bps_ids_global

        for index, row in haplotype_df.iterrows():
            if not row['chromosome'] in cent_chr:
                continue
            cent_chr_index = cent_chr.index(row['chromosome'])
            if not haplotype_df.loc[index, 'start'] - 1 == cent_starts[cent_chr_index] and not haplotype_df.loc[index, 'end'] == cent_ends[cent_chr_index]:
                if index < len(haplotype_df) and haplotype_df.loc[index, 'start'] > 0:
                    haplotype_df.loc[index, 'start'] = haplotype_df.loc[index - 1, 'end'] + 1

    return haplotype_df


def split_regions_by_points(starts, ends, values1, values2, split_points):
    split_points = sorted(set(split_points))
    new_starts = []
    new_ends = []
    new_values1 = []
    new_values2 = []

    for start, end, val1, val2 in zip(starts, ends, values1, values2):
        # Get split points within this region
        internal_splits = [pt for pt in split_points if start < pt < end]
        boundaries = [start] + internal_splits + [end]

        for i in range(len(boundaries) - 1):
            new_starts.append(boundaries[i])
            new_ends.append(boundaries[i + 1])
            new_values1.append(val1)
            new_values2.append(val2)

    return new_starts, new_ends, new_values1, new_values2


"""
#TODO: likely bug in this function
def _overlap_check(start, end, starts, ends):
  for i in range(len(starts)):
    if (start < ends[i] and end > starts[i]) or (start <= starts[i] and end >= ends[i]):
      return True
  return False


def merge_regions(starts, ends, values, values_next, loh_starts, loh_ends, threshold=3):
    merged_starts, merged_ends, merged_values, merged_values_2 = [], [], [], []
    i = 0
    while i < len(starts):
        start, end, value = starts[i], ends[i], values[i]
        j = i + 1
        while j < len(starts) and abs(value - values[j]) <= threshold and abs(value - values_next[i]) >= threshold:
            end = max(end, ends[j])
            j += 1
        if not _overlap_check(start, end, loh_starts, loh_ends):
            merged_starts.append(start)
            merged_ends.append(end)
            merged_values.append(sum(values[i:j]) / (j - i))
            merged_values_2.append(sum(values_next[i:j]) / (j - i))
        i = j
    return merged_starts, merged_ends, merged_values, merged_values_2
"""

