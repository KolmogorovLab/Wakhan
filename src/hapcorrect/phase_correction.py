import os
import statistics
import numpy as np
import pandas as pd
from collections import Counter
from itertools import groupby
from src.hapcorrect.utils import csv_df_chromosomes_sorter, df_chromosomes_sorter, write_segments_coverage_snps, get_contigs_list, find_peak_median_without_outliers

import logging
logger = logging.getLogger()

def generate_phasesets_bins(bam, path, bin_size, args):
    return get_phasesets_bins(bam, path, bin_size, args)


def get_phasesets_bins(bam, phasesets, bin_size, args):
    #indices, values, chroms = remove_overlapping_and_small_phasesets(phasesets, bin_size, args)
    dfs = csv_df_chromosomes_sorter(phasesets, ['chr', 'position', 'phasesets'])
    dfs_filtered = remove_overlapping_and_small_phasesets(dfs)

    head, tail = os.path.split(bam)
    bed = []
    for idx, s in dfs_filtered.iterrows():
        bed.append([tail, dfs_filtered.loc[idx, 'chr'], dfs_filtered.loc[idx, 'phaseset_start'], dfs_filtered.loc[idx, 'phaseset_end']])
    return bed


def check_missing_phasesets_original(df1, df2):
    df2 = df2[['chr', 'start', 'end']]
    if len(df1) == 0:
        return df2
    if len(df2) == 0:
        return df1

    def is_overlapping(row, df):
        overlapping = df[
            (df['chr'] == row['chr']) &
            (df['end'] >= row['start']) &
            (df['start'] <= row['end'])
            ]
        return not overlapping.empty

    # Filter df2: keep only those rows that do NOT overlap with any in df1
    non_overlapping_rows = df2[~df2.apply(lambda row: is_overlapping(row, df1), axis=1)]

    # Combine non-overlapping df2 rows into df1
    df1_extended = pd.concat([df1, non_overlapping_rows], ignore_index=True)

    # Sort
    df1_extended = df_chromosomes_sorter(df1_extended, ['chr', 'start', 'end'])

    return df1_extended


def remove_overlapping_and_small_phasesets(df):
    # Step 1: Aggregate to get start and end of each phaseset
    intervals = (
        df.groupby(['chr', 'phasesets'])['position']
        .agg(['min', 'max'])
        .reset_index()
        .rename(columns={'min': 'start', 'max': 'end'})
    )

    # Step 2: Sort by chr, then start position
    intervals.sort_values(['chr', 'start'], inplace=True)

    # Step 3: Process each chromosome separately
    result = []
    for chrom, group in intervals.groupby('chr'):
        # Sort by start position, then by longest duration (end - start) descending
        group = group.copy()
        group['length'] = group['end'] - group['start']
        group.sort_values(by=['start', 'length'], ascending=[True, False], inplace=True)

        merged = []
        for _, row in group.iterrows():
            if not merged:
                merged.append(row)
            else:
                last = merged[-1]
                # Check for overlap
                if row['start'] <= last['end']:
                    # Overlap: keep the longer one
                    if row['length'] > last['length']:
                        merged[-1] = row
                else:
                    merged.append(row)
        for m in merged:
            result.append({'chr': chrom, 'phaseset_start': m['start'], 'phaseset_end': m['end']})

    return pd.DataFrame(result)


def closest(lst):
    s = sorted(set(lst))
    return min([[a, b] for a, b in zip(s, s[1:])], key=lambda x: x[1] - x[0])


def merge_contiguous_indices(indices, haplotype_1_values, haplotype_2_values, ref_start_values, ref_start_values_phasesets, 
                             ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, consider_remaining):
    LARGE_BLOCK = 1000000
    new_starts = []
    new_ends = []
    hp_1_values = []
    hp_2_values = []

    #indices_slices = [sorted(list(set(indices_sub_list))) for indices_sub_list in indices]
    indices_slices = indices

    for i, value in enumerate(indices_slices):
        first_imbalance = _has_imbalance(haplotype_1_values_phasesets[value[0]], haplotype_2_values_phasesets[value[0]])
        first_len = ref_end_values_phasesets[value[0]] - ref_start_values_phasesets[value[0]]

        if len(value) > 1 or (first_imbalance and first_len > LARGE_BLOCK):
        #if len(value) > 1:
            internal_bins = [k for k in ref_start_values if k >= ref_start_values_phasesets[value[0]] and k <= ref_end_values_phasesets[value[-1]]]
            if internal_bins:
                new_starts.append(ref_start_values_phasesets[value[0]])
                new_ends.append(ref_end_values_phasesets[value[-1]])
                hp_1_values.append(median_data(haplotype_1_values[ref_start_values.index(internal_bins[0]):ref_start_values.index(internal_bins[-1])]))
                hp_2_values.append(median_data(haplotype_2_values[ref_start_values.index(internal_bins[0]):ref_start_values.index(internal_bins[-1])]))

    if len(hp_1_values) == 0:
        return ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets

    if consider_remaining:
        indices = []
        for sublist in indices_slices:
            if len(sublist) == 1:
                indices.append(sublist[0])
        haplotype_1_values_phasesets = [i for j, i in enumerate(haplotype_1_values_phasesets) if j in indices] + hp_1_values
        haplotype_2_values_phasesets = [i for j, i in enumerate(haplotype_2_values_phasesets) if j in indices] + hp_2_values
        ref_start_values_phasesets = [i for j, i in enumerate(ref_start_values_phasesets) if j in indices] + new_starts
        ref_end_values_phasesets = [i for j, i in enumerate(ref_end_values_phasesets) if j in indices] + new_ends
    else:
        haplotype_1_values_phasesets = hp_1_values
        haplotype_2_values_phasesets = hp_2_values
        ref_start_values_phasesets = new_starts
        ref_end_values_phasesets = new_ends

    zipped = list(zip(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets))
    zipped_sorted = sorted(zipped, key=lambda x: x[0])
    ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = zip(*zipped_sorted)

    return list(ref_start_values_phasesets), list(ref_end_values_phasesets), list(haplotype_1_values_phasesets), list(haplotype_2_values_phasesets)


def phase_contiguous_cis_trans(breakpoints_additional, indices_merge, haplotype_1_values, haplotype_2_values, ref_start_values, 
                               haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):
    values_ps = []
    for index, value in enumerate(ref_start_values_phasesets):
        values_ps.append([ref_start_values_phasesets[index], ref_end_values_phasesets[index]])

    anti_indices_merge = []
    if len(ref_start_values_phasesets) > 1:
        for i in range(len(ref_start_values_phasesets)-1):
            # Cis = min(|A1 - B1|, |A2 - B2|), Trans = min(|A1-B2|, |A2 - B1|)
            cis_right = min(abs(haplotype_1_values_phasesets[i] - haplotype_1_values_phasesets[i+1]),
                            abs(haplotype_2_values_phasesets[i] - haplotype_2_values_phasesets[i+1]))
            trans_right = min(abs(haplotype_1_values_phasesets[i] - haplotype_2_values_phasesets[i+1]),
                              abs(haplotype_2_values_phasesets[i] - haplotype_1_values_phasesets[i+1]))
            #if trans_right < cis_right: # and ref_start_values_phasesets[i+1] - ref_end_values_phasesets[i] < 1000000:
            #    indices_merge.append(i)
            if haplotype_1_values_phasesets[i] > haplotype_2_values_phasesets[i] and not haplotype_1_values_phasesets[i+1] > haplotype_2_values_phasesets[i+1] or \
               haplotype_1_values_phasesets[i] < haplotype_2_values_phasesets[i] and not haplotype_1_values_phasesets[i + 1] < haplotype_2_values_phasesets[i + 1]:# or \
               #abs(haplotype_1_values_phasesets[i + 1] - haplotype_2_values_phasesets[i + 1]) < 3:
                anti_indices_merge.append(ref_end_values_phasesets[i])

    ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = merge_contiguous_indices(indices_merge, haplotype_1_values, haplotype_2_values, ref_start_values, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, True)

    return haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets



def split_on_singletons_and_gaps(data):
    # Step 1: remove duplicates while preserving order
    seen = set()
    deduped = []
    for x in data:
        if x not in seen:
            deduped.append(x)
            seen.add(x)

    # Step 2: count occurrences in original data
    count = Counter(data)

    # Step 3: split on singleton or gap
    result = []
    group = []

    for i, val in enumerate(deduped):
        if i > 0 and val - deduped[i - 1] > 1:
            result.append(group)
            group = []

        group.append(val)

        if count[val] == 1:
            result.append(group)
            group = []

    if group:
        result.append(group)

    return result


def _has_imbalance(val_1, val_2):
    MIN_HP_DIFF_RATE = 0.3
    MIN_HP_DIFF_ABS = 10
    min_diff = min(MIN_HP_DIFF_RATE * min(val_1, val_2), MIN_HP_DIFF_ABS)
    return abs(val_1 - val_2) > min_diff


def find_indices_to_be_merged(mean_cis_trans_ps, ref_start_values_phasesets, ref_end_values_phasesets, 
                              haplotype_1_values_phasesets, haplotype_2_values_phasesets):

    small_ps = []
    for i, (start,end) in enumerate(zip(ref_start_values_phasesets, ref_end_values_phasesets)):
        if end - start < 200000:
            small_ps.append(i)

    if mean_cis_trans_ps:
       mean_cis_trans_ps = statistics.median(mean_cis_trans_ps)
       #print('mean_cis_trans_ps', mean_cis_trans_ps*1)
    else:
       mean_cis_trans_ps = 5

    ##| A-B | / min(A, B)
    to_merge = []
    if len(ref_start_values_phasesets) > 1:
        cur_merged = [0]
        for i in range(1, len(ref_start_values_phasesets)):
            len_fst = ref_end_values_phasesets[i - 1] - ref_start_values_phasesets[i - 1]
            len_snd = ref_end_values_phasesets[i] - ref_start_values_phasesets[i]
            gap = ref_start_values_phasesets[i] - ref_end_values_phasesets[i - 1]

            if ((haplotype_1_values_phasesets[i - 1] < haplotype_2_values_phasesets[i - 1]) ==
                    (haplotype_1_values_phasesets[i] < haplotype_2_values_phasesets[i])) and \
                _has_imbalance(haplotype_1_values_phasesets[i - 1], haplotype_2_values_phasesets[i - 1]) and \
                _has_imbalance(haplotype_1_values_phasesets[i], haplotype_2_values_phasesets[i]) and \
                gap < (len_fst + len_snd) / 2:

                cur_merged.append(i)
            else:
                to_merge.append(cur_merged)
                cur_merged = [i]

        if cur_merged:
            to_merge.append(cur_merged)
    #print(len(ref_start_values_phasesets), to_merge)

    return to_merge


def updated_means_ps(ref_start_values, haplotype_1_values, haplotype_2_values, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets):
    for i in range(len(ref_start_values_phasesets)-1):
        internal_bins = [k for k in ref_start_values if k >= ref_start_values_phasesets[i + 1] and k <= ref_end_values_phasesets[i + 1]]
        if haplotype_1_values_phasesets[i] > haplotype_2_values_phasesets[i]:
            if internal_bins:
                j = ref_start_values.index(internal_bins[0])
                for l in range(len(internal_bins)):
                    if haplotype_1_values[j] < haplotype_2_values[j]:
                        new_hp2 = haplotype_2_values[j]
                        new_hp1 = haplotype_1_values[j]
                        haplotype_1_values[j] = new_hp2
                        haplotype_2_values[j] = new_hp1
                    j = j + 1
        elif haplotype_1_values_phasesets[i]< haplotype_2_values_phasesets[i]:
            if internal_bins:
                j = ref_start_values.index(internal_bins[0])
                for l in range(len(internal_bins)):
                    if haplotype_1_values[j] > haplotype_2_values[j]:
                        new_hp2 = haplotype_2_values[j]
                        new_hp1 = haplotype_1_values[j]
                        haplotype_1_values[j] = new_hp2
                        haplotype_2_values[j] = new_hp1
                    j = j + 1

        if internal_bins:
            haplotype_1_values_phasesets[i] = median_data([float(x) for x in haplotype_1_values[ref_start_values.index(internal_bins[0]):ref_start_values.index(internal_bins[-1])] if x != 0])
            haplotype_2_values_phasesets[i] = median_data([float(x) for x in haplotype_2_values[ref_start_values.index(internal_bins[0]):ref_start_values.index(internal_bins[-1])] if x != 0])

    return haplotype_1_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values


def phase_blocks_updated_coverage(args, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values, haplotype_2_values):
    values_phasesets_hp1 = []
    values_phasesets_hp2 = []

    for i, val in enumerate(ref_start_values_phasesets):
        values_phasesets_hp1.append(median_data(haplotype_1_values[ref_start_values_phasesets[i]//args.bin_size:ref_end_values_phasesets[i]//args.bin_size]))
        values_phasesets_hp2.append(median_data(haplotype_2_values[ref_start_values_phasesets[i]//args.bin_size:ref_end_values_phasesets[i]//args.bin_size]))

    return ref_start_values_phasesets, ref_end_values_phasesets, values_phasesets_hp1, values_phasesets_hp2

def remove_centromere_phaseblocks(haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                                  ref_start_values_phasesets, ref_end_values_phasesets, centromere_region):
    cen = centromere_region.iloc[0]
    new_median_1 = []
    new_median_2 = []
    new_starts = []
    new_ends = []
    for i, (bs, be) in enumerate(zip(ref_start_values_phasesets, ref_end_values_phasesets)):
        if cen['start'] < bs and be < cen['end']:   #contained
            continue

        new_median_1.append(haplotype_1_values_phasesets[i])
        new_median_2.append(haplotype_2_values_phasesets[i])

        if bs < cen['start'] and cen['start'] < be:
            new_starts.append(bs)
            new_ends.append(cen['start'])
        elif bs < cen['end'] and cen['end'] < be:
            new_starts.append(cen['end'])
            new_ends.append(be)
        else:
            new_starts.append(bs)
            new_ends.append(be)

    return new_median_1, new_median_2, new_starts, new_ends


def phase_flips_cis_trans(chrom, args, breakpoints_additional, haplotype_1_values, haplotype_2_values,
                          ref_start_values, ref_end_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                          ref_start_values_phasesets, ref_end_values_phasesets, internal_ps=False, bins_adjust=False,
                          merge=False, swap_final=False, flank_coverage=False):

    FLANK_LEN = 5000000 // args.bin_size

    def get_block_coverage(idx):
        start_bin, end_bin = ref_start_values_phasesets[idx] // args.bin_size, ref_end_values_phasesets[idx] // args.bin_size
        return haplotype_1_values[start_bin : end_bin + 1], haplotype_2_values[start_bin : end_bin + 1]

    mean_cis_trans_ps = []
    values_ps = []
    for index, value in enumerate(ref_start_values_phasesets):
        values_ps.append([ref_start_values_phasesets[index], ref_end_values_phasesets[index]])

    broken_phasesets = []
    if internal_ps:
    # inside phaseblocks phaseswitch errors
        ref_start_values_phasesets_ = ref_start_values_phasesets
        ref_end_values_phasesets_ = ref_end_values_phasesets
        (broken_phasesets, ref_start_values_phasesets, ref_end_values_phasesets,
         haplotype_1_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values) = \
                 scan_and_update_phaseswitches_inside_phaseblocks(args, chrom, values_ps, haplotype_1_values, haplotype_2_values,
                                                     ref_start_values, ref_end_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                                                     ref_start_values_phasesets, ref_end_values_phasesets)
    if bins_adjust:
        haplotype_1_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values = \
                updated_means_ps(ref_start_values, haplotype_1_values, haplotype_2_values,
                                 ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets)

    indices_merge = []
    if len(ref_start_values_phasesets) > 1:
        for i in range(len(ref_start_values_phasesets)-1):
            prev_cov_1, prev_cov_2 = haplotype_1_values_phasesets[i], haplotype_2_values_phasesets[i]
            next_cov_1, next_cov_2 = haplotype_1_values_phasesets[i + 1], haplotype_2_values_phasesets[i + 1]

            if flank_coverage:
                #print(prev_cov_1, prev_cov_2, next_cov_1, next_cov_2)
                prev_bins_1, prev_bins_2 = get_block_coverage(i)
                next_bins_1, next_bins_2 = get_block_coverage(i + 1)
                prev_cov_1, prev_cov_2 = np.median(prev_bins_1[-FLANK_LEN:]), np.median(prev_bins_2[-FLANK_LEN:])
                next_cov_1, next_cov_2 = np.median(next_bins_1[:FLANK_LEN]), np.median(next_bins_2[:FLANK_LEN])
                #print(prev_cov_1, prev_cov_2, next_cov_1, next_cov_2)

            # Cis = min(|A1 - B1|, |A2 - B2|), Trans = min(|A1-B2|, |A2 - B1|)
            cis_right = min(abs(prev_cov_1 - next_cov_1), abs(prev_cov_2 - next_cov_2))
            trans_right = min(abs(prev_cov_1 - next_cov_2), abs(prev_cov_2 - next_cov_1))
            if flank_coverage:
                cis_right = abs(prev_cov_1 - next_cov_1) * abs(prev_cov_2 - next_cov_2)
                trans_right = abs(prev_cov_1 - next_cov_2) * abs(prev_cov_2 - next_cov_1)
                #print("Cis", cis_right, "trans", trans_right)
                #print("")

            if trans_right < cis_right:
                new_hp2_ps = haplotype_2_values_phasesets[i+1]
                new_hp1_ps = haplotype_1_values_phasesets[i+1]
                haplotype_1_values_phasesets[i+1] = new_hp2_ps
                haplotype_2_values_phasesets[i+1] = new_hp1_ps

                #dict = []
                #dict.append((chrom + '\t' + str(ref_start_values_phasesets[i+1]) + '\t' + str(ref_end_values_phasesets[i+1])))
                #write_segments_coverage_snps(dict, args.genome_name + '_phase_change_segments.csv', args)

                indices_merge.append(i)
                indices_merge.append(i+1)
                #mean_cis_trans_ps.append(abs(new_hp1_ps-new_hp2_ps))

                # Compute RMS deviation from median
                bin_values = np.array(haplotype_1_values[ref_start_values_phasesets[i]//args.bin_size:ref_end_values_phasesets[i]//args.bin_size])
                rmsd_hp1 = np.sqrt(np.mean((bin_values - np.median(bin_values)) ** 2))
                bin_values = np.array(haplotype_2_values[ref_start_values_phasesets[i] // args.bin_size:ref_end_values_phasesets[i] // args.bin_size])
                rmsd_hp2 = np.sqrt(np.mean((bin_values - np.median(bin_values)) ** 2))
                mean_cis_trans_ps.append(statistics.mean([rmsd_hp1, rmsd_hp2]))

                internal_bins = [k for k in ref_start_values if k >= ref_start_values_phasesets[i+1] and k <= ref_end_values_phasesets[i+1]]
                if internal_bins:
                    j = ref_start_values.index(internal_bins[0])
                    for l, bin in enumerate(internal_bins): #for l in range(len(internal_bins)):
                        new_hp2 = haplotype_2_values[j]
                        new_hp1 = haplotype_1_values[j]
                        haplotype_1_values[j] = new_hp2
                        haplotype_2_values[j] = new_hp1

                        dict = []
                        dict.append((chrom + '\t' + str(bin) + '\t' + str(bin + args.bin_size - 1)))
                        write_segments_coverage_snps(dict, args.genome_name + '_phase_change_segments.csv', args)

                        j = j + 1

                    #if swap_final:
                    #    haplotype_1_values, haplotype_2_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets = swap_final_bins_and_ps(chrom, args, ref_start_values, ref_end_values_phasesets, i, j, haplotype_1_values, haplotype_2_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets)
    if internal_ps:
        ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = \
                phase_blocks_updated_coverage(args, ref_start_values_phasesets_, ref_end_values_phasesets_, haplotype_1_values, haplotype_2_values)

    if merge:
        haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = \
                phase_contiguous_cis_trans(breakpoints_additional, indices_merge, haplotype_1_values, haplotype_2_values,
                                           ref_start_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                                           ref_start_values_phasesets, ref_end_values_phasesets)

    return (broken_phasesets, mean_cis_trans_ps, haplotype_1_values, haplotype_2_values,
            haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)


def swap_final_bins_and_ps(chrom, args, ref_start_values, ref_end_values_phasesets, index_ps, index_bin,
                           haplotype_1_values, haplotype_2_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets):
    if index_ps + 2 <= len(haplotype_1_values_phasesets)+1:
        for i in range(index_ps + 2, len(haplotype_1_values_phasesets)):
            new_hp2_ps = haplotype_2_values_phasesets[i]
            new_hp1_ps = haplotype_1_values_phasesets[i]
            haplotype_1_values_phasesets[i] = new_hp2_ps
            haplotype_2_values_phasesets[i] = new_hp1_ps

    internal_bins = [k for k in ref_start_values if k >= ref_end_values_phasesets[index_ps+2] and k <= ref_start_values[-1]]#[k for k in ref_start_values if k >= ref_start_values_phasesets[i + 1] and k <= ref_end_values_phasesets[i + 1]]
    for j, bin in enumerate(internal_bins):
        new_hp2 = haplotype_2_values[j]
        new_hp1 = haplotype_1_values[j]
        haplotype_1_values[j] = new_hp2
        haplotype_2_values[j] = new_hp1

        dict = []
        dict.append((chrom + '\t' + str(bin) + '\t' + str(bin + args.bin_size - 1)))
        write_segments_coverage_snps(dict, args.genome_name + '_phase_change_segments.csv', args)

    return haplotype_1_values, haplotype_2_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets


def scan_and_update_phaseswitches_inside_phaseblocks(args, chrom, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):
    broken_phasesets = []
    first_values = []
    last_values = []
    values_hp1 = []
    values_hp2 = []

    for index, value_ps in enumerate(values_ps):
        # if chrom == 'chr17' and value_ps[0] == 28951555:
        #     print('here')
        internal_bins = [i for i in ref_start_values if i >= value_ps[0] and i <= value_ps[1]]
        if internal_bins:
            sub_list = np.arange(ref_start_values.index(internal_bins[0]), ref_start_values.index(internal_bins[len(internal_bins)-1]), 1).tolist()
            haplotype_1_higher = [list(v) for k, v in groupby(sub_list, lambda sub_list: haplotype_1_values[sub_list] > haplotype_2_values[sub_list]) if k]
            haplotype_2_higher = [list(v) for k, v in groupby(sub_list, lambda sub_list: haplotype_2_values[sub_list] > haplotype_1_values[sub_list]) if k]

            if haplotype_1_higher: last_haplotype_1_index = haplotype_1_higher[-1][-1]
            else: last_haplotype_1_index = 0
            if haplotype_2_higher: last_haplotype_2_index = haplotype_2_higher[-1][-1]
            else: last_haplotype_2_index = 0
            if last_haplotype_1_index > last_haplotype_2_index:
                last_index = last_haplotype_1_index
            elif last_haplotype_1_index < last_haplotype_2_index:
                last_index = last_haplotype_2_index
            else: last_index = 0

            [haplotype_1_higher.remove(ss) for ss in haplotype_1_higher[::-1] if len(ss) < 5]
            [haplotype_2_higher.remove(ss) for ss in haplotype_2_higher[::-1] if len(ss) < 5]

            def is_stable(lst, threshold=5):
                lst = np.array(lst)
                return np.std(lst) < threshold

            means_hp1_updated = []
            means_hp2_updated = []
            if len(haplotype_1_higher) > 0 and len(haplotype_2_higher) > 0:
                # to check if one HP is stable, ignore then and continue
                for mi in range(len(haplotype_1_higher)):
                    means_hp1_updated.append(median_data(haplotype_1_values[haplotype_1_higher[mi][0]:haplotype_1_higher[mi][-1]]))
                    #means_hp1_updated.append(median_data(haplotype_2_values[haplotype_1_higher[mi][0]:haplotype_1_higher[mi][-1]]))

                for mi in range(len(haplotype_2_higher)):
                    means_hp2_updated.append(median_data(haplotype_2_values[haplotype_2_higher[mi][0]:haplotype_2_higher[mi][-1]]))
                    #means_hp2_updated.append(median_data(haplotype_1_values[haplotype_2_higher[mi][0]:haplotype_2_higher[mi][-1]]))

                if (is_stable(means_hp1_updated) and len(means_hp1_updated) >= 1)  or (is_stable(means_hp2_updated) and len(means_hp2_updated) >= 1):
                    if len(haplotype_1_higher) > 1:
                        temp=[]
                        for i in range(len(haplotype_1_higher)-1):
                            if haplotype_1_higher[i+1][0] - haplotype_1_higher[i][-1] < 5:
                                temp.append(i)
                            else: temp.append(None)
                        haplotype_1_higher_list = [list(group) for key, group in groupby(temp, key=lambda x: x == None) if not key]
                        for item in haplotype_1_higher_list:
                            for i in item:
                                haplotype_1_higher[item[0]].extend(haplotype_1_higher[i+1])
                        for item in reversed(haplotype_1_higher_list):
                            for i in reversed(item):
                                haplotype_1_higher.pop(i+1)

                    if len(haplotype_2_higher) > 1:
                        temp=[]
                        for i in range(len(haplotype_2_higher)-1):
                            if haplotype_2_higher[i+1][0] - haplotype_2_higher[i][-1] < 5:
                                temp.append(i)
                            else: temp.append(None)
                        haplotype_2_higher_list = [list(group) for key, group in groupby(temp, key=lambda x: x == None) if not key]
                        for item in haplotype_2_higher_list:
                            for i in item:
                                haplotype_2_higher[item[0]].extend(haplotype_2_higher[i+1])
                        for item in reversed(haplotype_2_higher_list):
                            for i in reversed(item):
                                haplotype_2_higher.pop(i+1)

                    first_value = []
                    last_value = []
                    mean_value_haplotype_1 = []
                    mean_value_haplotype_2 = []

                    if len(haplotype_1_higher):
                        for i in range(len(haplotype_1_higher)):
                            first_value.append(ref_start_values[haplotype_1_higher[i][0]])
                            last_value.append(ref_end_values[haplotype_1_higher[i][-1]])
                            mean_value_haplotype_1.append(median_data(haplotype_1_values[haplotype_1_higher[i][0]:haplotype_1_higher[i][-1]]))
                            mean_value_haplotype_2.append(median_data(haplotype_2_values[haplotype_1_higher[i][0]:haplotype_1_higher[i][-1]]))

                    if len(haplotype_2_higher):
                        for i in range(len(haplotype_2_higher)):
                            first_value.append(ref_start_values[haplotype_2_higher[i][0]])
                            last_value.append(ref_end_values[haplotype_2_higher[i][-1]])
                            mean_value_haplotype_1.append(median_data(haplotype_1_values[haplotype_2_higher[i][0]:haplotype_2_higher[i][-1]]))
                            mean_value_haplotype_2.append(median_data(haplotype_2_values[haplotype_2_higher[i][0]:haplotype_2_higher[i][-1]]))

                    if len(first_value) > 1:
                        sort_function = lambda x: x[0]
                        sort_target = list(zip(first_value, last_value, mean_value_haplotype_1, mean_value_haplotype_2))
                        sort_target.sort(key=sort_function)

                        first_value = [a for a,b,c,d in sort_target]
                        last_value = [b for a,b,c,d in sort_target]
                        mean_value_haplotype_1 = [c for a,b,c,d in sort_target]
                        mean_value_haplotype_2 = [d for a,b,c,d in sort_target]

                    if len(mean_value_haplotype_1) > 1:
                        if inside_phasesets_alternate_hps_threshold_check(mean_value_haplotype_1, mean_value_haplotype_2, ref_start_values, value_ps, haplotype_1_values, haplotype_2_values):
                            broken_phasesets.append(index)
                            first_values.extend([value_ps[0]] + first_value[1:])
                            last_values.extend(last_value[:-1] + [value_ps[1]])
                            values_hp1.extend(mean_value_haplotype_1)
                            values_hp2.extend(mean_value_haplotype_2)
                else:
                    i = ref_start_values.index(internal_bins[0])
                    for j, bin in enumerate(internal_bins):
                        if haplotype_1_values_phasesets[index] > haplotype_2_values_phasesets[index]:
                            if haplotype_1_values[i] < haplotype_2_values[i]:
                                new_hp2 = haplotype_2_values[i]
                                new_hp1 = haplotype_1_values[i]
                                haplotype_1_values[i] = new_hp2
                                haplotype_2_values[i] = new_hp1
                                dict = []
                                dict.append((chrom + '\t' + str(bin) + '\t' + str(bin + args.bin_size - 1)))
                                write_segments_coverage_snps(dict, args.genome_name + '_phase_change_segments.csv', args)
                        elif haplotype_1_values_phasesets[index] < haplotype_2_values_phasesets[index]:
                            if haplotype_1_values[i] > haplotype_2_values[i]:
                                new_hp2 = haplotype_2_values[i]
                                new_hp1 = haplotype_1_values[i]
                                haplotype_1_values[i] = new_hp2
                                haplotype_2_values[i] = new_hp1
                                dict = []
                                dict.append((chrom + '\t' + str(bin) + '\t' + str(bin + args.bin_size - 1)))
                                write_segments_coverage_snps(dict, args.genome_name + '_phase_change_segments.csv', args)
                        i = i + 1
                    haplotype_1_values_phasesets[index] = mean_values(haplotype_1_values, internal_bins[0]//args.bin_size, internal_bins[-1]//args.bin_size)
                    haplotype_2_values_phasesets[index] = mean_values(haplotype_2_values, internal_bins[0]//args.bin_size, internal_bins[-1]//args.bin_size)

    if broken_phasesets:
        ref_start_values_phasesets = [item for i, item in enumerate(ref_start_values_phasesets) if i not in broken_phasesets] + first_values
        ref_end_values_phasesets = [item for i, item in enumerate(ref_end_values_phasesets) if i not in broken_phasesets] + last_values
        haplotype_1_values_phasesets = [item for i, item in enumerate(haplotype_1_values_phasesets) if i not in broken_phasesets] + values_hp1
        haplotype_2_values_phasesets = [item for i, item in enumerate(haplotype_2_values_phasesets) if i not in broken_phasesets] + values_hp2

        zipped = list(zip(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets))
        zipped_sorted = sorted(zipped, key=lambda x: x[0])
        ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = zip(*zipped_sorted)

    return (broken_phasesets, list(ref_start_values_phasesets), list(ref_end_values_phasesets),
            list(haplotype_1_values_phasesets), list(haplotype_2_values_phasesets), haplotype_1_values, haplotype_2_values)


def inside_phasesets_alternate_hps_threshold_check(mean_value_haplotype_1, mean_value_haplotype_2, ref_start_values, value_ps, hp1_values, hp2_values):
    #if abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]) > 3:
    if temp_hps_coverage_allocation_check(ref_start_values, value_ps, hp1_values, hp2_values):
        #for i in range(len(mean_value_haplotype_1)-1):
        #    if abs(mean_value_haplotype_1[i] - mean_value_haplotype_2[i+1]) < 6 or abs(mean_value_haplotype_1[i+1] - mean_value_haplotype_2[i]) < 6:
        return True
    else:
        return False


def temp_hps_coverage_allocation_check(ref_start_values, value_ps, haplotype_1_values, haplotype_2_values):

    internal_bins = [k for k in ref_start_values if k >= value_ps[0] and k <= value_ps[1]]
    temp_hp1_values = haplotype_1_values[ref_start_values.index(internal_bins[0]):ref_start_values.index(internal_bins[-1])]
    temp_hp2_values = haplotype_2_values[ref_start_values.index(internal_bins[0]):ref_start_values.index(internal_bins[-1])]

    for i in range(len(internal_bins)-1):
        if temp_hp2_values[i] > temp_hp1_values[i]:
            new_hp2 = temp_hp2_values[i]
            new_hp1 = temp_hp1_values[i]
            temp_hp1_values[i] = new_hp2
            temp_hp2_values[i] = new_hp1

    if abs(mean_values(temp_hp1_values, 0, len(temp_hp1_values)-1) - mean_values(temp_hp2_values, 0, len(temp_hp1_values)-1)) > 6:
        return True
    else:
        return False

def switch_inter_phaseblocks_bins(chrom, args, ref_start_values, haplotype_1_values, haplotype_2_values, ref_start_values_phasesets, 
                                  ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets):
    inter_phaseblocks_start_phasesets, inter_phaseblocks_end_phasesets, inter_phaseblocks_hp_values = \
            inter_phaseblocks_segments(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets)

    for i, (inter_phaseblocks_start, inter_phaseblocks_end, inter_phaseblocks_hp) in \
            enumerate(zip(inter_phaseblocks_start_phasesets, inter_phaseblocks_end_phasesets, inter_phaseblocks_hp_values)):
        internal_bins = [i for i in ref_start_values if i >= inter_phaseblocks_start and i <= inter_phaseblocks_end]
        if internal_bins:
            i = ref_start_values.index(internal_bins[0])

            for j, bin in enumerate(internal_bins):
                if inter_phaseblocks_hp == 1:
                    if haplotype_1_values[i] < haplotype_2_values[i]:
                        new_hp2 = haplotype_2_values[i]
                        new_hp1 = haplotype_1_values[i]
                        haplotype_1_values[i] = new_hp2
                        haplotype_2_values[i] = new_hp1
                        dict = []
                        dict.append((chrom + '\t' + str(bin) + '\t' + str(bin + args.bin_size - 1)))
                        write_segments_coverage_snps(dict, args.genome_name + '_phase_change_segments.csv', args)
                elif inter_phaseblocks_hp == 2:
                    if haplotype_1_values[i] > haplotype_2_values[i]:
                        new_hp2 = haplotype_2_values[i]
                        new_hp1 = haplotype_1_values[i]
                        haplotype_1_values[i] = new_hp2
                        haplotype_2_values[i] = new_hp1
                        dict = []
                        dict.append((chrom + '\t' + str(bin) + '\t' + str(bin + args.bin_size - 1)))
                        write_segments_coverage_snps(dict, args.genome_name + '_phase_change_segments.csv', args)
                i = i +1

    return haplotype_1_values, haplotype_2_values

def inter_phaseblocks_segments(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets):
    inter_phaseblocks_start_values_phasesets = []
    inter_phaseblocks_end_values_phasesets = []
    inter_phaseblocks_hp_values = []
    for i in range(len(ref_start_values_phasesets)-1):
        if haplotype_1_values_phasesets[i] > haplotype_2_values_phasesets[i] and haplotype_1_values_phasesets[i+1] > haplotype_2_values_phasesets[i+1]:
            inter_phaseblocks_start_values_phasesets.append(ref_end_values_phasesets[i] + 1)
            inter_phaseblocks_end_values_phasesets.append(ref_start_values_phasesets[i+1]-1)
            inter_phaseblocks_hp_values.append(1)
        elif haplotype_1_values_phasesets[i] < haplotype_2_values_phasesets[i] and haplotype_1_values_phasesets[i+1] < haplotype_2_values_phasesets[i+1]:
            inter_phaseblocks_start_values_phasesets.append(ref_end_values_phasesets[i] + 1)
            inter_phaseblocks_end_values_phasesets.append(ref_start_values_phasesets[i+1]-1)
            inter_phaseblocks_hp_values.append(2)

    return inter_phaseblocks_start_values_phasesets, inter_phaseblocks_end_values_phasesets, inter_phaseblocks_hp_values

def phaseblock_flipping_simple_heuristics(chrom, args, is_simple_correction, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):

    values_ps = []
    for index, value in enumerate(ref_start_values_phasesets):
        values_ps.append([ref_start_values_phasesets[index], ref_end_values_phasesets[index]])

    #inside phaseblocks phaseswitch errors
    broken_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values = scan_and_update_phaseswitches_inside_phaseblocks(args, chrom, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                                                     haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
    if len(ref_start_values_phasesets) > 1:
        hp_changed=[]
        for i in range(len(ref_start_values_phasesets) - 1):
            if haplotype_1_values_phasesets[i] > haplotype_2_values_phasesets[i]:
                hp_changed.append([ref_start_values_phasesets[i], ref_end_values_phasesets[i]+1])
                haplotype_1_values_phasesets[i], haplotype_2_values_phasesets[i] = haplotype_2_values_phasesets[i], haplotype_1_values_phasesets[i]

        for i in range(len(ref_start_values)):
            for j in range(len(hp_changed)):
                if (ref_start_values[i] >= hp_changed[j][0] and ref_start_values[i] <= hp_changed[j][1]):
                    new_hp2 = haplotype_2_values[i]
                    new_hp1 = haplotype_1_values[i]
                    haplotype_1_values[i]=new_hp2
                    haplotype_2_values[i]=new_hp1
                elif haplotype_1_values[i] > haplotype_2_values[i]:
                    new_hp2 = haplotype_2_values[i]
                    new_hp1 = haplotype_1_values[i]
                    haplotype_1_values[i]=new_hp2
                    haplotype_2_values[i]=new_hp1
                    break
                else:
                    haplotype_1_values[i]=haplotype_1_values[i]
                    haplotype_2_values[i]=haplotype_2_values[i]
    else:
        for i in range(len(ref_start_values)):
            if haplotype_2_values[i] > haplotype_1_values[i]:
                new_hp2 = haplotype_2_values[i]
                new_hp1 = haplotype_1_values[i]
                haplotype_1_values[i] = new_hp2
                haplotype_2_values[i] = new_hp1

    return haplotype_1_values, haplotype_2_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets
def mean_values(selected_list, start_index, end_index):
    result = []
    for i in range(end_index - start_index):
        try:
            result.append(selected_list[start_index + i])
        except IndexError:
            break
    if result:
        return np.mean(result)
    else:
        return 0.0

def median_data(data):
    data = [float(x) for x in data if x != 0]
    if len(data) > 0:
        return np.median(data)
    else:
        return 0

def reintroduce_broken_phasesets(broken_indices, ref_start_values, haplotype1_mean, haplotype2_mean, phaseset_starts, phaseset_ends, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets):
    still_seperate = []

    for i, ps in enumerate(broken_indices):
        internals = []
        for j, (start,end) in enumerate(zip(ref_start_values_phasesets, ref_end_values_phasesets)):
            if (start >= phaseset_starts[ps] and end <= phaseset_ends[ps]):
                internals.append(j)
        if internals:
            still_seperate.append(internals)

    indices = [i for i in range(len(ref_start_values_phasesets))]
    for i in indices:
        if not i in [item for sublist in still_seperate for item in sublist]:
            still_seperate.append([i])

    ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = merge_contiguous_indices(still_seperate, haplotype1_mean, haplotype2_mean, ref_start_values, ref_start_values_phasesets,
        ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, True)

    return haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets


def update_remaining_phasesets(indices_merge, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets,
    ref_start_values_phasesets_merged, ref_end_values_phasesets_merged, haplotype_1_values_phasesets_merged, haplotype_2_values_phasesets_merged):
    merged = []
    for i, (start, end) in enumerate(zip(ref_start_values_phasesets, ref_end_values_phasesets)):
        for j, (start_merged, end_merged) in enumerate(zip(ref_start_values_phasesets_merged, ref_end_values_phasesets_merged)):
            if start >= start_merged and end <= end_merged:
                merged.append(i)

    merged = sorted(list(set(merged)))

    haplotype_1_values_phasesets = [i for j, i in enumerate(haplotype_1_values_phasesets) if j not in merged] + haplotype_1_values_phasesets_merged
    haplotype_2_values_phasesets = [i for j, i in enumerate(haplotype_2_values_phasesets) if j not in merged] + haplotype_2_values_phasesets_merged
    ref_start_values_phasesets = [i for j, i in enumerate(ref_start_values_phasesets) if j not in merged] + ref_start_values_phasesets_merged
    ref_end_values_phasesets = [i for j, i in enumerate(ref_end_values_phasesets) if j not in merged] + ref_end_values_phasesets_merged

    zipped = list(zip(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets))
    zipped_sorted = sorted(zipped, key=lambda x: x[0])
    ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = zip(*zipped_sorted)

    return list(ref_start_values_phasesets), list(ref_end_values_phasesets), list(haplotype_1_values_phasesets), list(haplotype_2_values_phasesets)

def subtract_intervals(df1, df2):
    result = []


    # Specify the columns to compare
    columns_to_compare = ["chr", "start", "end"]

    # Sort and reset index to ensure order doesn't affect comparison
    df1_sorted = df1[columns_to_compare].sort_values(by=columns_to_compare).reset_index(drop=True)
    df2_sorted = df2[columns_to_compare].sort_values(by=columns_to_compare).reset_index(drop=True)

    # Check equality
    same_values = df1_sorted.equals(df2_sorted)
    if same_values:
        return pd.DataFrame({'chr': [], 'start': [], 'end': []})

    for _, row1 in df1.iterrows():
        chr1, start1, end1 = row1['chr'], row1['start'], row1['end']
        temp_intervals = [(start1, end1)]

        for _, row2 in df2[df2['chr'] == chr1].iterrows():
            start2, end2 = row2['start'], row2['end']
            new_temp = []

            for s, e in temp_intervals:
                if end2 <= s or start2 >= e:
                    # no overlap
                    new_temp.append((s, e))
                else:
                    # partial overlap
                    if s < start2:
                        new_temp.append((s, start2))
                    if end2 < e:
                        new_temp.append((end2, e))
            temp_intervals = new_temp

        for s, e in temp_intervals:
            if s < e:
                result.append({'chr': chr1, 'start': s, 'end': e})

    return pd.DataFrame(result)

def without_phasesets_bins_correction(loh_chrom, args, chrom, ref_start_values, snps_haplotype1_mean, snps_haplotype2_mean, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets):
    ref_start_values_phasesets_new = []
    ref_end_values_phasesets_new = []
    haplotype_1_values_phasesets_new = []
    haplotype_2_values_phasesets_new = []
    starts_loh = loh_chrom.start.values.tolist()
    ends_loh = loh_chrom.end.values.tolist()

    start = False
    end = False

    if len(ref_start_values_phasesets) < 1:
        return snps_haplotype1_mean, snps_haplotype2_mean

    for i in range(len(ref_start_values_phasesets)-1):
        if ref_start_values_phasesets[i+1] - ref_end_values_phasesets[i] > args.bin_size * 2:
            ref_start_values_phasesets_new.append(ref_end_values_phasesets[i]+1)
            ref_end_values_phasesets_new.append(ref_start_values_phasesets[i+1]-1)
            haplotype_1_values_phasesets_new.append(median_data(snps_haplotype1_mean[(ref_end_values_phasesets[i]+1)//args.bin_size:(ref_start_values_phasesets[i+1]-1)//args.bin_size]))
            haplotype_2_values_phasesets_new.append(median_data(snps_haplotype2_mean[(ref_end_values_phasesets[i]+1)//args.bin_size:(ref_start_values_phasesets[i+1]-1)//args.bin_size]))

    if len(ref_start_values_phasesets) > 1:
        if ref_start_values_phasesets[0] > ref_start_values[1] and ref_start_values_phasesets[0] - ref_start_values[0] > args.bin_size * 2:
            start = True
            ref_start_values_phasesets_new.append(ref_start_values[0])
            ref_end_values_phasesets_new.append(ref_start_values_phasesets[0]-1)
            haplotype_1_values_phasesets_new.append(median_data(snps_haplotype1_mean[0:(ref_start_values_phasesets[0]-1)//args.bin_size]))
            haplotype_2_values_phasesets_new.append(median_data(snps_haplotype2_mean[0:(ref_start_values_phasesets[0]-1)//args.bin_size]))

        if ref_end_values_phasesets[-1] < ref_start_values[-1] and ref_start_values[-1] - ref_end_values_phasesets[-1]  > args.bin_size * 2:
            end = True
            ref_start_values_phasesets_new.append(ref_end_values_phasesets[-1]+1)
            ref_end_values_phasesets_new.append(ref_start_values[-1])
            haplotype_1_values_phasesets_new.append(median_data(snps_haplotype1_mean[(ref_end_values_phasesets[-1]+1)//args.bin_size:ref_start_values[-1]//args.bin_size]))
            haplotype_2_values_phasesets_new.append(median_data(snps_haplotype2_mean[(ref_end_values_phasesets[-1]+1)//args.bin_size:ref_start_values[-1]//args.bin_size]))

    haplotype_1_values_phasesets = [i for j, i in enumerate(haplotype_1_values_phasesets)] + haplotype_1_values_phasesets_new
    haplotype_2_values_phasesets = [i for j, i in enumerate(haplotype_2_values_phasesets)] + haplotype_2_values_phasesets_new
    ref_start_values_phasesets = [i for j, i in enumerate(ref_start_values_phasesets)] + ref_start_values_phasesets_new
    ref_end_values_phasesets = [i for j, i in enumerate(ref_end_values_phasesets)] + ref_end_values_phasesets_new
    zipped = list(zip(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets))
    zipped_sorted = sorted(zipped, key=lambda x: x[0])
    ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = zip(*zipped_sorted)

    ref_start_values_phasesets = list(ref_start_values_phasesets)
    ref_end_values_phasesets = list(ref_end_values_phasesets)
    haplotype_1_values_phasesets = list(haplotype_1_values_phasesets)
    haplotype_2_values_phasesets = list(haplotype_2_values_phasesets)

    indices = []
    for i, val in enumerate(ref_start_values_phasesets):
        if start and i == 0:
            continue
        if end and i == len(ref_start_values_phasesets)-1:
            continue
        if val in ref_start_values_phasesets_new:
            indices.append(i)


    for i, val in enumerate(indices):
        if ((haplotype_1_values_phasesets[val-1] > haplotype_2_values_phasesets[val-1] and haplotype_1_values_phasesets[val+1] > haplotype_2_values_phasesets[val+1]) and haplotype_1_values_phasesets[val] < haplotype_2_values_phasesets[val]) or\
                ((haplotype_1_values_phasesets[val-1] < haplotype_2_values_phasesets[val-1] and haplotype_1_values_phasesets[val+1] < haplotype_2_values_phasesets[val+1]) and haplotype_1_values_phasesets[val] > haplotype_2_values_phasesets[val]) or \
                ((abs(haplotype_1_values_phasesets[val-1] - haplotype_2_values_phasesets[val]) < abs(haplotype_1_values_phasesets[val-1] - haplotype_1_values_phasesets[val])) and (abs(haplotype_1_values_phasesets[val+1] - haplotype_2_values_phasesets[val]) < abs(haplotype_1_values_phasesets[val+1] - haplotype_1_values_phasesets[val]))):

            internal_bins = [k for k in ref_start_values if k >= ref_start_values_phasesets[val] and k <= ref_end_values_phasesets[val]]  # [k for k in ref_start_values if k >= ref_start_values_phasesets[i + 1] and k <= ref_end_values_phasesets[i + 1]]
            j = ref_start_values.index(internal_bins[0])
            for l, bin in enumerate(internal_bins):
                new_hp2 = snps_haplotype2_mean[j]
                new_hp1 = snps_haplotype1_mean[j]
                snps_haplotype1_mean[j] = new_hp2
                snps_haplotype2_mean[j] = new_hp1
                j = j + 1

                dict = []
                dict.append((chrom + '\t' + str(bin) + '\t' + str(bin + args.bin_size - 1)))
                write_segments_coverage_snps(dict, args.genome_name + '_phase_change_segments.csv', args)

    if start:
        if (haplotype_1_values_phasesets[0] > haplotype_2_values_phasesets[0] and haplotype_1_values_phasesets[1] < haplotype_2_values_phasesets[1]) or (haplotype_1_values_phasesets[0] < haplotype_2_values_phasesets[0] and haplotype_1_values_phasesets[1] > haplotype_2_values_phasesets[1]):
            val = 0
            internal_bins = [k for k in ref_start_values if k >= ref_start_values_phasesets[val] and k <= ref_end_values_phasesets[val]]  # [k for k in ref_start_values if k >= ref_start_values_phasesets[i + 1] and k <= ref_end_values_phasesets[i + 1]]
            j = ref_start_values.index(internal_bins[0])
            for l, bin in enumerate(internal_bins):
                if not any(start <= ref_start_values[j] <= end for start, end in zip(starts_loh, ends_loh)):
                    new_hp2 = snps_haplotype2_mean[j]
                    new_hp1 = snps_haplotype1_mean[j]
                    snps_haplotype1_mean[j] = new_hp2
                    snps_haplotype2_mean[j] = new_hp1

                    dict = []
                    dict.append((chrom + '\t' + str(bin) + '\t' + str(bin + args.bin_size - 1)))
                    write_segments_coverage_snps(dict, args.genome_name + '_phase_change_segments.csv', args)
                j = j + 1
    if end:
        if (haplotype_1_values_phasesets[-1] > haplotype_2_values_phasesets[-1] and haplotype_1_values_phasesets[-2] < haplotype_2_values_phasesets[-2]) or (haplotype_1_values_phasesets[-1] < haplotype_2_values_phasesets[-1] and haplotype_1_values_phasesets[-2] > haplotype_2_values_phasesets[-2]):
            val = -1
            internal_bins = [k for k in ref_start_values if k >= ref_start_values_phasesets[val] and k <= ref_end_values_phasesets[val]]  # [k for k in ref_start_values if k >= ref_start_values_phasesets[i + 1] and k <= ref_end_values_phasesets[i + 1]]
            j = ref_start_values.index(internal_bins[0])

            for l, bin in enumerate(internal_bins):
                if not any(start <= ref_start_values[j] <= end for start, end in zip(starts_loh, ends_loh)):
                    new_hp2 = snps_haplotype2_mean[j]
                    new_hp1 = snps_haplotype1_mean[j]
                    snps_haplotype1_mean[j] = new_hp2
                    snps_haplotype2_mean[j] = new_hp1


                    dict = []
                    dict.append((chrom + '\t' + str(bin) + '\t' + str(bin + args.bin_size - 1)))
                    write_segments_coverage_snps(dict, args.genome_name + '_phase_change_segments.csv', args)
                j = j + 1
    return snps_haplotype1_mean, snps_haplotype2_mean
