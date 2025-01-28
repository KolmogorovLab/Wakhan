import os
import statistics
import numpy as np
from utils import csv_df_chromosomes_sorter, write_segments_coverage_dict
from vcf_processing import squash_regions
from collections import  defaultdict, Counter
import pysam
import bisect
import pandas as pd
import logging

logger = logging.getLogger()

from utils import get_contigs_list, remove_outliers_iqr
def generate_phasesets_bins(bam, path, bin_size, args):
    return get_phasesets_bins(bam, path, bin_size, args)
def get_phasesets_bins(bam, phasesets, bin_size, args):
    indices, values = remove_overlapping_and_small_phasesets(phasesets, bin_size, args)
    head, tail = os.path.split(bam)
    bed = []
    from utils import get_contigs_list
    chroms = get_contigs_list(args.contigs)
    for ind, chrom in enumerate(chroms) :
        for i in range(0, len(values[ind]), 2): #for i in range(len(values[ind])-1):
            start = values[ind][i]
            end = values[ind][i+1]
            if end - start > args.bin_size*6:
                bed.append([tail, chrom, start, end])
    return bed

def remove_overlapping_and_small_phasesets(phasesets, bin_size, args):
    #dfs = pd.read_csv(phasesets, sep='\t', names=['chr', 'pos', 'ps'])
    dfs = csv_df_chromosomes_sorter(phasesets, ['chr', 'pos', 'ps'])
    values_all = []
    indices_all = []
    from utils import get_contigs_list
    chroms = get_contigs_list(args.contigs)
    #
    for index, chrom in enumerate(chroms):
        df = dfs[dfs['chr'] == chrom]
        unique_ps_by_chr = df.groupby('chr', group_keys=True)['ps'].apply(lambda x: list(np.unique(x)))
        # Find overlapping phase blocks locations
        for unique in unique_ps_by_chr[0]:
            indices_overlapping = []
            indices = np.where(df["ps"] == unique)
            list_new = []
            for i in indices[0]:
                list_new.append(i)
        for i in range(len(list_new) - 1):
            if not list_new[i + 1] - list_new[i] == 1:
                indices_overlapping.append(list_new[i] + 1)
                indices_overlapping.append(list_new[i + 1] - list_new[i] - 1)
        for i in range(0, len(indices_overlapping) - 1, 2):
            df = df.drop(labels=range(indices_overlapping[i], indices_overlapping[i] + indices_overlapping[i + 1]),
                         axis=0)
        final = []
        # Remove overlapping phase blocks
        ps = df.ps.values.tolist()
        pos = df.pos.values.tolist()
        unique_ps_by_chr = df.groupby('chr')['ps'].apply(lambda x: list(np.unique(x)))
        for x in range(len(unique_ps_by_chr[0])):
            if x == len(unique_ps_by_chr[0])-1:
                final.append(pos[ps.index(unique_ps_by_chr[0][x])])
                final.append(pos[len(ps) - 1 - ps[::-1].index(unique_ps_by_chr[0][x])])
            else:
                if (unique_ps_by_chr[0][x + 1] - unique_ps_by_chr[0][x] > bin_size):
                    final.append(pos[ps.index(unique_ps_by_chr[0][x])])
                    final.append(pos[len(ps) - 1 - ps[::-1].index(unique_ps_by_chr[0][x])])

        index = []
        # Remove small (< bin_size) phase blocks
        for i in range(len(final)):
            init = 0
            for m in range(1, final[-1], bin_size):
                init += 1
                if final[i] > m and final[i] < m + bin_size:
                    index.append(init)
        indices_all.append(index)
        values_all.append(final)
    return indices_all, values_all


def closest(lst):
    s = sorted(set(lst))
    return min([[a, b] for a, b in zip(s, s[1:])], key=lambda x: x[1] - x[0])

# def inside_phaseblock_switch_correction_coverage_update(chrom, args, is_simple_correction, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
#                     haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):

def phaseblock_flipping(chrom, args, is_simple_correction, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):

    values_ps = []
    for index, value in enumerate(ref_start_values_phasesets):
        values_ps.append([ref_start_values_phasesets[index], ref_end_values_phasesets[index]])

    #inside phaseblocks phaseswitch errors
    scan_and_update_phaseswitches_inside_phaseblocks(values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                                                     haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
    if is_simple_correction == False:
        diff = []
        values_ps = []
        #ref_start_values_phasesets.sort()
        #ref_end_values_phasesets.sort()

        for index, value in enumerate(ref_start_values_phasesets):
            diff.append(ref_end_values_phasesets[index] - ref_start_values_phasesets[index])
            values_ps.append([ref_start_values_phasesets[index], ref_end_values_phasesets[index]])
        max_value_index = diff.index(max(diff))
        #print(ref_start_values_phasesets[max_value_index], ref_end_values_phasesets[max_value_index])

        #phaseblocks switch errors
        scan_and_update_phaseblocks_switch_errors(chrom, args, max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                        haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

        #case for whole HPs switch, ie., HCC1437 - chr8
        #phase_switch_spanning_haplotypes(chrom, args, max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
        #                haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

        #bins without phaseblocks
        bins_without_phaseblocks(chrom, args, max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                        haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

        #amplified signals
        #
    else:
        if len(ref_start_values_phasesets) > 1:
            hp_changed=[]
            for i in range(len(ref_start_values_phasesets) - 1):
                if haplotype_1_values_phasesets[i] > haplotype_2_values_phasesets[i]:
                    hp_changed.append([ref_start_values_phasesets[i], ref_end_values_phasesets[i]+1])

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
    #haplotype_1_values_updated.extend(haplotype_1_values)
    #haplotype_2_values_updated.extend(haplotype_2_values)


def scan_and_update_phaseswitches_inside_phaseblocks(values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):
    for index, value_ps in enumerate(values_ps):
        if value_ps[0] == 113196847 :#and value_ps[1] == 143633540:# or value_ps[0] == 72108881:#72346221:#132135757:#177448845:#81331913:#72346221:
            print("here")
        #print(value_ps[0])

        internal_bins = [i for i in ref_start_values if i >= value_ps[0] and i <= value_ps[1]]
        if internal_bins:
            from itertools import groupby
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
                    mean_value_haplotype_1.append(statistics.mean(haplotype_1_values[haplotype_1_higher[i][0]:haplotype_1_higher[i][-1]]))
                    mean_value_haplotype_2.append(statistics.mean(haplotype_2_values[haplotype_1_higher[i][0]:haplotype_1_higher[i][-1]]))

            if len(haplotype_2_higher):
                for i in range(len(haplotype_2_higher)):
                    first_value.append(ref_start_values[haplotype_2_higher[i][0]])
                    last_value.append(ref_end_values[haplotype_2_higher[i][-1]])
                    mean_value_haplotype_1.append(statistics.mean(haplotype_1_values[haplotype_2_higher[i][0]:haplotype_2_higher[i][-1]]))
                    mean_value_haplotype_2.append(statistics.mean(haplotype_2_values[haplotype_2_higher[i][0]:haplotype_2_higher[i][-1]]))

            #first_value = sorted(first_value)
            #last_value = sorted(last_value)
            if len(first_value) > 1:
                sort_function = lambda x: x[0]
                sort_target = list(zip(first_value, last_value, mean_value_haplotype_1, mean_value_haplotype_2))
                sort_target.sort(key=sort_function)

                first_value = [a for a,b,c,d in sort_target]
                last_value = [b for a,b,c,d in sort_target]
                mean_value_haplotype_1 = [c for a,b,c,d in sort_target]
                mean_value_haplotype_2 = [d for a,b,c,d in sort_target]

            # if last_value and last_index > ref_end_values.index(last_value[-1]) + 10:
            #     first_value.append(last_value[-1]+1)
            #     last_value.append(ref_end_values[last_index])
            #     mean_value_haplotype_1.append(statistics.mean(haplotype_1_values[ref_start_values.index(first_value[-1]):last_index]))
            #     mean_value_haplotype_2.append(statistics.mean(haplotype_2_values[ref_start_values.index(first_value[-1]):last_index]))
            if (len(mean_value_haplotype_1) == 2  and \
                    (abs(mean_value_haplotype_1[0] - mean_value_haplotype_2[1]) < 6 or abs(mean_value_haplotype_1[1] - mean_value_haplotype_2[0]) < 6))  or \
                    (len(mean_value_haplotype_1) == 3 and \
                     ((abs(mean_value_haplotype_1[0] - mean_value_haplotype_2[1]) < 6 or abs(mean_value_haplotype_1[1] - mean_value_haplotype_2[0]) < 6) or \
                            (abs(mean_value_haplotype_1[1] - mean_value_haplotype_2[2]) < 6 or abs(mean_value_haplotype_1[2] - mean_value_haplotype_2[1]) < 6))):
                start_index = ref_start_values_phasesets.index(value_ps[0])
                ref_start_values_phasesets.pop(start_index)
                ref_end_values_phasesets.pop(start_index)
                haplotype_1_values_phasesets.pop(start_index)
                haplotype_2_values_phasesets.pop(start_index)

                ref_start_values_phasesets[start_index:start_index] = first_value
                ref_end_values_phasesets[start_index:start_index] = last_value
                haplotype_1_values_phasesets[start_index:start_index] = mean_value_haplotype_1
                haplotype_2_values_phasesets[start_index:start_index] = mean_value_haplotype_2

def scan_and_update_phaseblocks_switch_errors(chrom, args, max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):

    PHASESETS_DIFF_THRESHOLD = 12
    #############################################################################################
    # Scan Right-Left (HP1>HP2)
    if haplotype_1_values_phasesets[max_value_index] > haplotype_2_values_phasesets[max_value_index] or \
            haplotype_1_values_phasesets[max_value_index] == haplotype_2_values_phasesets[max_value_index]:
        values_ps_new = values_ps[max_value_index + 1:]
        if not max_value_index == len(ref_start_values_phasesets) - 1:
            for j, value_ps in enumerate(values_ps_new):
                if value_ps[0] == 121986191:  # 101790293 or value_ps[0] == 104622049:
                    print("here")
                if not ref_start_values_phasesets.index(value_ps[0]) == 0 and \
                        not ref_end_values_phasesets.index(value_ps[1]) == len(ref_start_values_phasesets):

                    next_ps_start = max([t for t in ref_start_values if
                                         ref_end_values_phasesets[ref_end_values_phasesets.index(value_ps[1]) - 1] > t])
                    last_ps_start = min([t for t in ref_start_values if ref_start_values_phasesets[
                        ref_start_values_phasesets.index(value_ps[0]) - 1] < t])

                    # if haplotype_1_values[ref_start_values.index(next_ps_start)] > haplotype_2_values[ref_start_values.index(next_ps_start)] and \
                    # haplotype_1_values[ref_start_values.index(last_ps_start)] > haplotype_2_values[ref_start_values.index(last_ps_start)]:

                    # if abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                    #        haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])-1] ) < 5 and \
                    #     abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                    #         haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+1]) < 5:
                    internal_bins = [i for i in ref_start_values if i >= value_ps[0] and i <= value_ps[1]]
                    if internal_bins:
                        i = ref_start_values.index(internal_bins[0])
                        if (((
                                # (haplotype_1_values[ref_start_values.index(last_ps_start)] > haplotype_2_values[ref_start_values.index(last_ps_start)] \
                                #  or haplotype_1_values[ref_start_values.index(last_ps_start)+1] > haplotype_2_values[ref_start_values.index(last_ps_start)+1] \
                                #  or haplotype_1_values[ref_start_values.index(last_ps_start)+2] > haplotype_2_values[ref_start_values.index(last_ps_start)]+2) \
                                (haplotype_1_values[ref_start_values.index(next_ps_start)] > haplotype_2_values[
                                    ref_start_values.index(next_ps_start)] \
                                 or haplotype_1_values[ref_start_values.index(next_ps_start) - 1] > haplotype_2_values[
                                     ref_start_values.index(next_ps_start) - 1] \
                                 or haplotype_1_values[ref_start_values.index(next_ps_start) - 2] > haplotype_2_values[
                                     ref_start_values.index(next_ps_start) - 2] \
                                 or haplotype_1_values[ref_start_values.index(next_ps_start) - 3] > haplotype_2_values[
                                     ref_start_values.index(next_ps_start) - 3] \
                                 or haplotype_1_values[ref_start_values.index(next_ps_start) - 4] > haplotype_2_values[
                                     ref_start_values.index(next_ps_start) - 4])
                                and (haplotype_1_values[i] < haplotype_2_values[i] or haplotype_1_values[i + 1] <
                                     haplotype_2_values[i + 1] or haplotype_1_values[i + 2] < haplotype_2_values[
                                         i + 2] or haplotype_1_values[i + 3] < haplotype_2_values[i + 3] or
                                     haplotype_1_values[i + 4] < haplotype_2_values[i + 4]))) \
                            or (abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                                    haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0]) - 1]) < PHASESETS_DIFF_THRESHOLD)) \
                                and haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] < \
                                haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]:
                            # or haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] > haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]: #3-5 bins check, both  sides of phaseblock
                            dict = []
                            if chrom == 'chr3' and 106100001 in internal_bins:
                                print('here')
                            dict.append((chrom + '\t' + str(internal_bins[0]) + '\t' + str(internal_bins[-1] + args.bin_size-1)))
                            write_segments_coverage_dict(dict, args.genome_name + '_phase_change_segments.csv')

                            new_hp2_ps = haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]
                            new_hp1_ps = haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]
                            haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] = new_hp2_ps
                            haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] = new_hp1_ps

                            for k in range(len(internal_bins)):
                                # if haplotype_1_values[i] < haplotype_2_values[i]:
                                new_hp2 = haplotype_2_values[i]
                                new_hp1 = haplotype_1_values[i]
                                haplotype_1_values[i] = new_hp2
                                haplotype_2_values[i] = new_hp1
                                i = i + 1
        # Scan Left-Right (HP1>HP2)
        if not max_value_index == 0:
            values_ps_new = values_ps[0:max_value_index:]
            for j, value_ps in reversed(list(enumerate(values_ps_new))):
                if value_ps[0] == 121986191:  # 101790293 or value_ps[0] == 104622049:
                    print("here")
                if not ref_end_values_phasesets.index(value_ps[1]) == len(ref_start_values_phasesets):

                    next_ps_start = max([t for t in ref_start_values if
                                         ref_end_values_phasesets[ref_end_values_phasesets.index(value_ps[1]) + 1] > t])
                    last_ps_start = min([t for t in ref_start_values if
                                         ref_start_values_phasesets[
                                             ref_start_values_phasesets.index(value_ps[0]) + 1] < t])

                    # if haplotype_1_values[ref_start_values.index(next_ps_start)] > haplotype_2_values[ref_start_values.index(next_ps_start)] and \
                    # haplotype_1_values[ref_start_values.index(last_ps_start)] > haplotype_2_values[ref_start_values.index(last_ps_start)]:

                    # if abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                    #        haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])-1] ) < 5 and \
                    #     abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                    #         haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+1]) < 5:
                    internal_bins = [i for i in ref_start_values if i >= value_ps[0] and i <= value_ps[1]]
                    if internal_bins:
                        i = ref_start_values.index(internal_bins[0])
                        if (((((haplotype_1_values[ref_start_values.index(last_ps_start)] > haplotype_2_values[
                            ref_start_values.index(last_ps_start)] \
                                or haplotype_1_values[ref_start_values.index(last_ps_start) + 1] > haplotype_2_values[
                                    ref_start_values.index(last_ps_start) + 1] \
                                or haplotype_1_values[ref_start_values.index(last_ps_start) + 2] > haplotype_2_values[
                                    ref_start_values.index(last_ps_start)] + 2 \
                                or haplotype_1_values[ref_start_values.index(last_ps_start) + 3] > haplotype_2_values[
                                    ref_start_values.index(last_ps_start)] + 3 \
                                or haplotype_1_values[ref_start_values.index(last_ps_start) + 4] > haplotype_2_values[
                                    ref_start_values.index(last_ps_start)] + 4) \
                               # and (haplotype_1_values[ref_start_values.index(next_ps_start)] > haplotype_2_values[ref_start_values.index(next_ps_start)]\
                               # or haplotype_1_values[ref_start_values.index(next_ps_start)-1] > haplotype_2_values[ref_start_values.index(next_ps_start)-1] \
                               # or haplotype_1_values[ref_start_values.index(next_ps_start)-2] > haplotype_2_values[ref_start_values.index(next_ps_start)-2])
                               and (haplotype_1_values[i] < haplotype_2_values[i] or haplotype_1_values[i + 1] <
                                    haplotype_2_values[i + 1] or haplotype_1_values[i + 2] < haplotype_2_values[
                                        i + 2] or haplotype_1_values[i + 3] < haplotype_2_values[i + 3] or
                                    haplotype_1_values[i + 4] < haplotype_2_values[i + 4]))) \
                             or (abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                                     haplotype_1_values_phasesets[ref_start_values_phasesets.index(
                                         value_ps[0]) + 1]) < PHASESETS_DIFF_THRESHOLD)) \
                                and haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] <
                                haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]):
                            # or haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] > haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]: #3-5 bins check, both  sides of phaseblock
                            dict = []
                            if chrom == 'chr3' and 106100001 in internal_bins:
                                print('here')
                            dict.append((chrom + '\t' + str(internal_bins[0]) + '\t' + str(internal_bins[-1] + args.bin_size-1)))
                            write_segments_coverage_dict(dict, args.genome_name + '_phase_change_segments.csv')

                            new_hp2_ps = haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]
                            new_hp1_ps = haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]
                            haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] = new_hp2_ps
                            haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] = new_hp1_ps

                            for k in range(len(internal_bins)):
                                # if haplotype_1_values[i] < haplotype_2_values[i]:
                                new_hp2 = haplotype_2_values[i]
                                new_hp1 = haplotype_1_values[i]
                                haplotype_1_values[i] = new_hp2
                                haplotype_2_values[i] = new_hp1
                                i = i + 1
    # Scan Right-Left (HP1<HP2)
    elif haplotype_1_values_phasesets[max_value_index] < haplotype_2_values_phasesets[max_value_index]:
        if not max_value_index == len(ref_start_values_phasesets) - 1:
            values_ps_new = values_ps[max_value_index + 1:]
            for j, value_ps in enumerate(values_ps_new):
                if value_ps[0] == 121986191:  # 101790293 or value_ps[0] == 104622049:
                    print("here")

                if not ref_start_values_phasesets.index(value_ps[0]) == 0 and \
                        not ref_end_values_phasesets.index(value_ps[1]) == len(ref_start_values_phasesets):

                    next_ps_start = max([t for t in ref_start_values if
                                         ref_end_values_phasesets[ref_end_values_phasesets.index(value_ps[1]) - 1] > t])
                    last_ps_start = min([t for t in ref_start_values if ref_start_values_phasesets[
                        ref_start_values_phasesets.index(value_ps[0]) - 1] < t])

                    # if haplotype_1_values[ref_start_values.index(next_ps_start)] > haplotype_2_values[ref_start_values.index(next_ps_start)] and \
                    # haplotype_1_values[ref_start_values.index(last_ps_start)] > haplotype_2_values[ref_start_values.index(last_ps_start)]:

                    # if abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                    #        haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])-1] ) < 5 and \
                    #     abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                    #         haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+1]) < 5:
                    internal_bins = [i for i in ref_start_values if i >= value_ps[0] and i <= value_ps[1]]
                    if internal_bins:
                        i = ref_start_values.index(internal_bins[0])
                        if (((
                                # (haplotype_1_values[ref_start_values.index(last_ps_start)] > haplotype_2_values[ref_start_values.index(last_ps_start)] \
                                #  or haplotype_1_values[ref_start_values.index(last_ps_start)+1] > haplotype_2_values[ref_start_values.index(last_ps_start)+1] \
                                #  or haplotype_1_values[ref_start_values.index(last_ps_start)+2] > haplotype_2_values[ref_start_values.index(last_ps_start)]+2) \
                                (haplotype_1_values[ref_start_values.index(next_ps_start)] < haplotype_2_values[
                                    ref_start_values.index(next_ps_start)] \
                                 or haplotype_1_values[ref_start_values.index(next_ps_start) - 1] < haplotype_2_values[
                                     ref_start_values.index(next_ps_start) - 1] \
                                 or haplotype_1_values[ref_start_values.index(next_ps_start) - 2] < haplotype_2_values[
                                     ref_start_values.index(next_ps_start) - 2] \
                                 or haplotype_1_values[ref_start_values.index(next_ps_start) - 3] < haplotype_2_values[
                                     ref_start_values.index(next_ps_start) - 3] \
                                 or haplotype_1_values[ref_start_values.index(next_ps_start) - 4] < haplotype_2_values[
                                     ref_start_values.index(next_ps_start) - 4])
                                #and (haplotype_1_values[i] > haplotype_2_values[i] or haplotype_1_values[i + 1] >
                                #     haplotype_2_values[i + 1] or haplotype_1_values[i + 2] > haplotype_2_values[
                                #         i + 2] or haplotype_1_values[i + 3] > haplotype_2_values[i + 3] or
                                #     haplotype_1_values[i + 4] > haplotype_2_values[i + 4])
                                and (mean_values(haplotype_1_values, i, i + 4) > mean_values(haplotype_2_values, i, i + 4)) )) \
                            or (abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                                    haplotype_1_values_phasesets[
                                        ref_start_values_phasesets.index(value_ps[0]) - 1]) < PHASESETS_DIFF_THRESHOLD)) \
                                and haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] > \
                                haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]:
                            # or haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] > haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]: #3-5 bins check, both  sides of phaseblock
                            dict = []
                            if chrom == 'chr3' and 106100001 in internal_bins:
                                print('here')
                            dict.append((chrom + '\t' + str(internal_bins[0]) + '\t' + str(internal_bins[-1] + args.bin_size-1)))
                            write_segments_coverage_dict(dict, args.genome_name + '_phase_change_segments.csv')

                            new_hp2_ps = haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]
                            new_hp1_ps = haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]
                            haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] = new_hp2_ps
                            haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] = new_hp1_ps

                            for k in range(len(internal_bins)):
                                # if haplotype_1_values[i] < haplotype_2_values[i]:
                                new_hp2 = haplotype_2_values[i]
                                new_hp1 = haplotype_1_values[i]
                                haplotype_1_values[i] = new_hp2
                                haplotype_2_values[i] = new_hp1
                                i = i + 1
        # Scan Left-Right (HP1<HP2)
        if not max_value_index == 0:
            values_ps_new = values_ps[0:max_value_index:]
            for j, value_ps in reversed(list(enumerate(values_ps_new))):
                if value_ps[0] == 121986191:  # 101790293 or value_ps[0] == 104622049:
                    print("here")
                if not ref_end_values_phasesets.index(value_ps[1]) == len(ref_start_values_phasesets):

                    next_ps_start = max([t for t in ref_start_values if
                                         ref_end_values_phasesets[ref_end_values_phasesets.index(value_ps[1]) + 1] > t])
                    last_ps_start = min([t for t in ref_start_values if
                                         ref_start_values_phasesets[
                                             ref_start_values_phasesets.index(value_ps[0]) + 1] < t])

                    # if haplotype_1_values[ref_start_values.index(next_ps_start)] > haplotype_2_values[ref_start_values.index(next_ps_start)] and \
                    # haplotype_1_values[ref_start_values.index(last_ps_start)] > haplotype_2_values[ref_start_values.index(last_ps_start)]:

                    # if abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                    #        haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])-1] ) < 5 and \
                    #     abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                    #         haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+1]) < 5:
                    internal_bins = [i for i in ref_start_values if i >= value_ps[0] and i <= value_ps[1]]
                    if internal_bins:
                        i = ref_start_values.index(internal_bins[0])
                        if ((((haplotype_1_values[ref_start_values.index(last_ps_start)] < haplotype_2_values[
                            ref_start_values.index(last_ps_start)] \
                               or haplotype_1_values[ref_start_values.index(last_ps_start) + 1] < haplotype_2_values[
                                   ref_start_values.index(last_ps_start) + 1] \
                               or haplotype_1_values[ref_start_values.index(last_ps_start) + 2] < haplotype_2_values[
                                   ref_start_values.index(last_ps_start)] + 2 \
                               or haplotype_1_values[ref_start_values.index(last_ps_start) + 3] < haplotype_2_values[
                                   ref_start_values.index(last_ps_start)] + 3 \
                               or haplotype_1_values[ref_start_values.index(last_ps_start) + 4] < haplotype_2_values[
                                   ref_start_values.index(last_ps_start)] + 4) \
                              # and (haplotype_1_values[ref_start_values.index(next_ps_start)] > haplotype_2_values[ref_start_values.index(next_ps_start)]\
                              # or haplotype_1_values[ref_start_values.index(next_ps_start)-1] > haplotype_2_values[ref_start_values.index(next_ps_start)-1] \
                              # or haplotype_1_values[ref_start_values.index(next_ps_start)-2] > haplotype_2_values[ref_start_values.index(next_ps_start)-2])
                              and (haplotype_1_values[i] > haplotype_2_values[i] or haplotype_1_values[i + 1] >
                                   haplotype_2_values[i + 1] or haplotype_1_values[i + 2] > haplotype_2_values[i + 2] or
                                   haplotype_1_values[i + 3] > haplotype_2_values[i + 3] or haplotype_1_values[i + 4] >
                                   haplotype_2_values[i + 4]))) \
                            or (abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                                    haplotype_1_values_phasesets[
                                        ref_start_values_phasesets.index(value_ps[0]) + 1]) < PHASESETS_DIFF_THRESHOLD)) \
                                and haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] > \
                                haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]:
                            # or haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] > haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]: #3-5 bins check, both  sides of phaseblock
                            dict = []
                            dict.append((chrom + '\t' + str(internal_bins[0]) + '\t' + str(internal_bins[-1] + args.bin_size-1)))
                            write_segments_coverage_dict(dict, args.genome_name + '_phase_change_segments.csv')

                            new_hp2_ps = haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]
                            new_hp1_ps = haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]
                            haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] = new_hp2_ps
                            haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] = new_hp1_ps

                            for k in range(len(internal_bins)):
                                # if haplotype_1_values[i] < haplotype_2_values[i]:
                                new_hp2 = haplotype_2_values[i]
                                new_hp1 = haplotype_1_values[i]
                                haplotype_1_values[i] = new_hp2
                                haplotype_2_values[i] = new_hp1
                                i = i + 1

def bins_without_phaseblocks(chrom, args, max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):

    print("here")

    # if haplotype_1_values_phasesets[max_value_index] > haplotype_2_values_phasesets[max_value_index] or haplotype_1_values_phasesets[max_value_index] == haplotype_2_values_phasesets[max_value_index]:
    #     for j in range(len(values_ps)-1):
    #         internal_bins = [i for i in ref_start_values if i > values_ps[j][1] and i < values_ps[j+1][0]]
    #         if internal_bins:
    #             i = ref_start_values.index(internal_bins[0])
    #             dict = []
    #             if chrom == 'chr3' and 106100001 in internal_bins:
    #                 print('here')
    #             dict.append(
    #                 (chrom + '\t' + str(internal_bins[0]) + '\t' + str(internal_bins[-1] + args.bin_size'] - 1)))
    #             write_segments_coverage_dict(dict, args.genome_name'] + '_phase_change_segments.csv')
    #             for k in range(len(internal_bins)):
    #                 if haplotype_1_values[i] < haplotype_2_values[i]:
    #                     new_hp2 = haplotype_2_values[i]
    #                     new_hp1 = haplotype_1_values[i]
    #                     haplotype_1_values[i] = new_hp2
    #                     haplotype_2_values[i] = new_hp1
    #                 i = i + 1
    # else:
    #     for j in range(len(values_ps) - 1):
    #         internal_bins = [i for i in ref_start_values if i > values_ps[j][1] and i < values_ps[j + 1][0]]
    #         if internal_bins:
    #             i = ref_start_values.index(internal_bins[0])
    #             dict = []
    #             dict.append(
    #                 (chrom + '\t' + str(internal_bins[0]) + '\t' + str(internal_bins[-1] + args.bin_size'] - 1)))
    #             write_segments_coverage_dict(dict, args.genome_name'] + '_phase_change_segments.csv')
    #             for k in range(len(internal_bins)):
    #                 if haplotype_1_values[i] > haplotype_2_values[i]:
    #                     new_hp2 = haplotype_2_values[i]
    #                     new_hp1 = haplotype_1_values[i]
    #                     haplotype_1_values[i] = new_hp2
    #                     haplotype_2_values[i] = new_hp1
    #                 i = i + 1


def phase_switch_spanning_haplotypes(chrom, args, max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):
    MIN_PHASESETS_DIFF_THRESHOLD = 3
    for index, value_ps in enumerate(values_ps):
        #if value_ps[0] == 93800001 and value_ps[1] == 66901398:
        #    print("here")
        if not ref_end_values_phasesets.index(value_ps[1]) == len(ref_start_values_phasesets)-1:
            if (abs(haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                   haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+1])  < MIN_PHASESETS_DIFF_THRESHOLD and \
                abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                    haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+1]) > 25) and \
                abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+1] - \
                        haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0]) + 1]) < MIN_PHASESETS_DIFF_THRESHOLD and \
                abs(haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0]) + 1] - \
                        haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0]) + 1]) < MIN_PHASESETS_DIFF_THRESHOLD and \
                    (not haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] < 1 and \
                    not haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] < 1 and \
                    not haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+1] < 1 and \
                    not haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+1] < 1) :
                internal_bins = [i for i in ref_start_values if i >= ref_start_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+1] and i <= ref_end_values_phasesets[ref_end_values_phasesets.index(value_ps[1])+1]]
                if internal_bins:
                    i = ref_start_values.index(internal_bins[0])
                    new_bins_start = ref_start_values[i:]
                    for k in range(len(new_bins_start)):
                        # if haplotype_1_values[i] < haplotype_2_values[i]:
                        new_hp2 = haplotype_2_values[i]
                        new_hp1 = haplotype_1_values[i]
                        haplotype_1_values[i] = new_hp2
                        haplotype_2_values[i] = new_hp1
                        i = i + 1
                break

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


#Apmlification cases
# or (abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
#         haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0]) + 1] > PHASESETS_DIFF_THRESHOLD) \
#     and abs(haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
#             haplotype_2_values_phasesets[
#                 ref_start_values_phasesets.index(value_ps[0]) + 1] < PHASESETS_DIFF_THRESHOLD // 2)):

def phase_correction_centers(args, segs_hp1_state, segs_hp1_start, segs_hp1_end,  segs_hp2_state, segs_hp2_start, segs_hp2_end, haplotype_1_values, haplotype_2_values):

    diff_hp1 = []
    for index, (start, end, state) in enumerate(zip(segs_hp1_start, segs_hp1_end, segs_hp1_state)):
        diff_hp1.append([end - start, index, state])
    diff_hp1.sort(key=lambda x: x[0])
    max_diff_hp1 = diff_hp1[-1]

    diff_hp2 = []
    for index, (start, end, state) in enumerate(zip(segs_hp2_start, segs_hp2_end, segs_hp2_state)):
        diff_hp2.append([end - start, index, state])
    diff_hp2.sort(key=lambda x: x[0])
    max_diff_hp2 = diff_hp2[-1]
    processed = False
    #if not ( (len(set(segs_hp2_state)) == 1 or any(segs_hp2_state) < 5) or (len(set(segs_hp1_state)) == 1 or any(segs_hp1_state) < 5) ):
    if not max_diff_hp1[2] == max_diff_hp2[2]:
        processed = True
        for index, (start, end, state) in enumerate(zip(segs_hp2_start, segs_hp2_end, segs_hp2_state)):
            if state == max_diff_hp1[2]:
                for i in range(start//args.bin_size,end//args.bin_size):
                    new_hp2 = haplotype_2_values[i]
                    new_hp1 = haplotype_1_values[i]
                    haplotype_1_values[i] = new_hp2
                    haplotype_2_values[i] = new_hp1

        for index, (start, end, state) in enumerate(zip(segs_hp1_start, segs_hp1_end, segs_hp1_state)):
            if state == max_diff_hp2[2]:
                for i in range(start//args.bin_size,end//args.bin_size):
                    new_hp2 = haplotype_2_values[i]
                    new_hp1 = haplotype_1_values[i]
                    haplotype_1_values[i] = new_hp2
                    haplotype_2_values[i] = new_hp1

    return haplotype_1_values, haplotype_2_values, processed

def contiguous_phaseblocks(haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):
    THRESHOLD = 8
    contiguous_hp1_start = []
    contiguous_hp1_end = []
    hp1_values = []
    contiguous_hp1_start.append(ref_start_values_phasesets[0])
    hp1_value = []
    for i in range(len(ref_start_values_phasesets) - 1):
        current_value_hp1 = haplotype_1_values_phasesets[i]
        hp1_value.append(haplotype_1_values_phasesets[i])
        if not haplotype_1_values_phasesets[i+1] - THRESHOLD < current_value_hp1 < haplotype_1_values_phasesets[i+1] + THRESHOLD and \
           abs(haplotype_2_values_phasesets[i] - haplotype_1_values_phasesets[i]) > THRESHOLD:
            contiguous_hp1_start.append(ref_start_values_phasesets[i+1])
            contiguous_hp1_end.append(ref_start_values_phasesets[i+1]-1)
            hp1_values.append(statistics.mean(hp1_value))
            hp1_value.clear()
    contiguous_hp1_end.append(ref_end_values_phasesets[-1])
    if len(contiguous_hp1_start) > 2:
        hp1_values.append(statistics.mean(haplotype_1_values_phasesets[ref_start_values_phasesets.index(contiguous_hp1_start[-2]):ref_start_values_phasesets.index(ref_start_values_phasesets[-1])]))
    else:
        hp1_values.append(0)



    contiguous_hp2_start = []
    contiguous_hp2_end = []
    hp2_values = []
    contiguous_hp2_start.append(ref_start_values_phasesets[0])
    hp2_value = []
    for i in range(len(ref_start_values_phasesets) - 1):
        current_value_hp2 = haplotype_2_values_phasesets[i]
        hp2_value.append(haplotype_2_values_phasesets[i])
        if not haplotype_2_values_phasesets[i+1] - THRESHOLD < current_value_hp2 < haplotype_2_values_phasesets[i+1] + THRESHOLD and \
           abs(haplotype_2_values_phasesets[i] - haplotype_1_values_phasesets[i]) > THRESHOLD:
            contiguous_hp2_start.append(ref_start_values_phasesets[i+1])
            contiguous_hp2_end.append(ref_start_values_phasesets[i+1]-1)
            hp2_values.append(statistics.mean(hp2_value))
            hp2_value.clear()
    contiguous_hp2_end.append(ref_end_values_phasesets[-1])
    if len(contiguous_hp2_start) > 2:
        hp2_values.append(statistics.mean(haplotype_2_values_phasesets[ref_start_values_phasesets.index(contiguous_hp2_start[-2]):ref_start_values_phasesets.index(ref_start_values_phasesets[-1])]))
    else:
        hp2_values.append(0)

    return hp1_values, contiguous_hp1_start, contiguous_hp1_end, hp2_values, contiguous_hp2_start, contiguous_hp2_end

def flip_phaseblocks_contigous(chrom, args, hp1_values, contiguous_hp1_start, contiguous_hp1_end, hp2_values, contiguous_hp2_start, contiguous_hp2_end, ref_start_values, haplotype_1_values_phasesets,
                               haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values, haplotype_2_values):
    THRESHOLD = 12
    diff_hp1 = []
    max_value = []
    for index, value in enumerate(contiguous_hp1_start):
        diff_hp1.append(contiguous_hp1_end[index] - contiguous_hp1_start[index])
    if len(diff_hp1) > 1:
        for i in range(1):
            max_value.append(hp1_values[diff_hp1.index(sorted(diff_hp1, reverse=True)[i])])#diff_hp1.index(max(diff_hp1))
            print(max_value)
        # max_value = sorted(max_value, reverse=True)
        # if max_value[0] - max_value[1] < 8:
        #     max_value.pop(1)
    else:
        max_value.append(hp1_values[diff_hp1.index(diff_hp1[0])])
    print(max_value)
    value_ps_updated = []
    for a in range(len(max_value)):
        for i, val in enumerate(ref_start_values_phasesets):
            value_ps = [ref_start_values_phasesets[i], ref_end_values_phasesets[i]]
            if haplotype_2_values_phasesets[i] - THRESHOLD < max_value[a] < haplotype_2_values_phasesets[i] + THRESHOLD and \
                not haplotype_1_values_phasesets[i] - THRESHOLD < max_value[a] < haplotype_1_values_phasesets[i] + THRESHOLD and \
                    not value_ps[0] in value_ps_updated:

                value_ps_updated.append(value_ps[0])
                internal_bins = [k for k in ref_start_values if k >= value_ps[0] and k <= value_ps[1]]
                if internal_bins:
                    i = ref_start_values.index(internal_bins[0])

                    dict = []
                    dict.append((chrom + '\t' + str(internal_bins[0]) + '\t' + str(internal_bins[-1] + args.bin_size - 1)))
                    write_segments_coverage_dict(dict, args.genome_name + '_phase_change_segments.csv')

                    new_hp2_ps = haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]
                    new_hp1_ps = haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]
                    haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] = new_hp2_ps
                    haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] = new_hp1_ps

                    for j in range(len(internal_bins)):
                        # if haplotype_1_values[i] < haplotype_2_values[i]:
                        new_hp2 = haplotype_2_values[i]
                        new_hp1 = haplotype_1_values[i]
                        haplotype_1_values[i] = new_hp2
                        haplotype_2_values[i] = new_hp1
                        i = i + 1

    diff_hp2 = []
    max_value = []
    for index, value in enumerate(contiguous_hp2_start):
        diff_hp2.append(contiguous_hp2_end[index] - contiguous_hp2_start[index])
    if len(diff_hp2) > 1:
        for i in range(1):
            max_value.append(hp2_values[diff_hp2.index(sorted(diff_hp2, reverse=True)[i])])#diff_hp1.index(max(diff_hp1))
            print(max_value)
        # max_value = sorted(max_value)
        # if max_value[1] - max_value[0] < 8:
        #     max_value.pop(1)
    else:
        max_value.append(hp2_values[diff_hp2.index(diff_hp2[0])])
    print(max_value)
    for b in range(len(max_value)):
        for i, val in enumerate(ref_start_values_phasesets):
            value_ps = [ref_start_values_phasesets[i], ref_end_values_phasesets[i]]
            if haplotype_1_values_phasesets[i] - THRESHOLD < max_value[b] < haplotype_1_values_phasesets[i] + THRESHOLD and \
                    not haplotype_2_values_phasesets[i] - THRESHOLD < max_value[b] < haplotype_2_values_phasesets[i] + THRESHOLD and \
                    not value_ps[0] in value_ps_updated:

                value_ps_updated.append(value_ps[0])
                internal_bins = [k for k in ref_start_values if k >= value_ps[0] and k <= value_ps[1]]
                if internal_bins:
                    i = ref_start_values.index(internal_bins[0])
                    dict = []
                    dict.append((chrom + '\t' + str(internal_bins[0]) + '\t' + str(internal_bins[-1] + args.bin_size - 1)))
                    write_segments_coverage_dict(dict, args.genome_name + '_phase_change_segments.csv')

                    new_hp2_ps = haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]
                    new_hp1_ps = haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]
                    haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] = new_hp2_ps
                    haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] = new_hp1_ps

                    for j in range(len(internal_bins)):
                        # if haplotype_1_values[i] < haplotype_2_values[i]:
                        new_hp2 = haplotype_2_values[i]
                        new_hp1 = haplotype_1_values[i]
                        haplotype_1_values[i] = new_hp2
                        haplotype_2_values[i] = new_hp1
                        i = i + 1

    return haplotype_1_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values
def detect_centromeres(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, hap_1_values, hap_2_values, ref_start_values, bin_size):
    centromere_bins = []
    for i in range(len(ref_start_values)):
        if hap_1_values[i] == 0 and hap_2_values[i] == 0:
            centromere_bins.append(ref_start_values[i])
            centromere_bins.append(ref_start_values[i] + bin_size)

    centromere_region_starts, centromere_region_ends = squash_regions(centromere_bins, bin_size)
    region_starts = []
    region_ends = []
    for i, (start,end) in enumerate(zip(centromere_region_starts, centromere_region_ends)):
        if end - start > bin_size * 6:
            region_starts.append(start)
            region_ends.append(end)
    if len(region_starts):
        ref_start_values_phasesets = ref_start_values_phasesets + region_starts
        ref_end_values_phasesets = ref_end_values_phasesets + region_ends
        haplotype_1_values_phasesets = haplotype_1_values_phasesets + [0 for i in range(len(region_starts))]
        haplotype_2_values_phasesets = haplotype_2_values_phasesets + [0 for i in range(len(region_starts))]

        sort_function = lambda x: x[0]
        sort_target = list(zip(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets))
        sort_target.sort(key=sort_function)

        ref_start_values_phasesets = [a for a, b, c, d in sort_target]
        ref_end_values_phasesets = [b for a, b, c, d in sort_target]
        haplotype_1_values_phasesets = [c for a, b, c, d in sort_target]
        haplotype_2_values_phasesets = [d for a, b, c, d in sort_target]

    return ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets


def rephase_vcf(df, vcf_in, out_vcf):
    chr_list = list(set(df['chr']))
    start_pos = defaultdict(list)
    end_pos = defaultdict(list)
    for seq in chr_list:
        start_pos[seq] = sorted([key for key, val in Counter(df.loc[df['chr'] == seq, 'start']).items() if val%2 == 1])
        end_pos[seq] = sorted([key for key, val in Counter(df.loc[df['chr'] == seq, 'end']).items() if val%2 == 1])
    vcf_in=pysam.VariantFile(vcf_in,"r")
    vcf_out = pysam.VariantFile(out_vcf, 'w', header=vcf_in.header)
    for var in vcf_in:
        sample = var.samples.keys()[0]
        if var.samples[sample].phased:
            strt = bisect.bisect_right(start_pos[var.chrom], var.pos)
            end = bisect.bisect_right(end_pos[var.chrom], var.pos)
            if strt == end + 1:
                (a,b) = var.samples[sample]['GT']
                new_gt = (abs(a-1), abs(b-1))
                var.samples[sample]['GT'] = new_gt
                var.samples[sample].phased = True
        vcf_out.write(var)
    vcf_out.close()


def remove_overlaping_contiguous(chrom, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2):
    overlaping_start_hp2 = []
    overlaping_end_hp2 = []
    overlaping_start_hp1 = []
    overlaping_end_hp1 = []
    for i, (start_hp1,end_hp1) in enumerate(zip(ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1)):
        for j, (start_hp2, end_hp2) in enumerate(zip(ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2)):
            if start_hp2 == 97490471 or start_hp1 == 97490471:
                print('here')
            if (start_hp2 >= start_hp1 and end_hp2 <= end_hp1) and not (start_hp2 == start_hp1 and end_hp2 == end_hp1) :
                overlaping_start_hp2.append(start_hp2)
                overlaping_end_hp2.append(end_hp2)

    for i, (start_hp2,end_hp2) in enumerate(zip(ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2)):
        for j, (start_hp1, end_hp1) in enumerate(zip(ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1)):
            if start_hp2 == 97490471 or start_hp1 == 97490471:
                print('here')
            if (start_hp1 >= start_hp2 and end_hp1 <= end_hp2) and not (start_hp1 == start_hp2 and end_hp1 == end_hp2):
                overlaping_start_hp1.append(start_hp1)
                overlaping_end_hp1.append(end_hp1)

    ref_start_values_phasesets_hp1 = [x for x in ref_start_values_phasesets_hp1 if x not in overlaping_start_hp1]
    ref_start_values_phasesets_hp2 = [x for x in ref_start_values_phasesets_hp2 if x not in overlaping_start_hp2]

    ref_end_values_phasesets_hp1 = [x for x in ref_end_values_phasesets_hp1 if x not in overlaping_end_hp1]
    ref_end_values_phasesets_hp2 = [x for x in ref_end_values_phasesets_hp2 if x not in overlaping_end_hp2]

    print(ref_start_values_phasesets_hp1)
    print(ref_end_values_phasesets_hp1)

    print(ref_start_values_phasesets_hp2)
    print(ref_end_values_phasesets_hp2)

    for i, (start_hp1,end_hp1) in enumerate(zip(ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1)):
        for j, (start_hp2, end_hp2) in enumerate(zip(ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2)):
            if start_hp2 < start_hp1 < end_hp2 and  end_hp1 > end_hp2 :
                print("First: "+str(start_hp1)+':'+str(end_hp1) +' , '+ str(start_hp2)+':'+str(end_hp2))
            if start_hp1 < start_hp2 < end_hp1 and  end_hp2 > end_hp1 :
                print("Second: "+str(start_hp1)+':'+ str(end_hp1) +' , '+ str(start_hp2)+':'+str(end_hp2))

    start_values_phasesets = np.unique(sorted(ref_start_values_phasesets_hp1 + ref_start_values_phasesets_hp2)).tolist()
    print(start_values_phasesets)
    chr_list = [chrom for ch in range(len(start_values_phasesets))]
    df_phasesets_chr = pd.DataFrame(list(zip(chr_list, start_values_phasesets)), columns=['chr', 'start'])

    return df_phasesets_chr

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

def switch_inter_phaseblocks_bins(chrom, args, ref_start_values, haplotype_1_values, haplotype_2_values, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets):
    inter_phaseblocks_start_phasesets, inter_phaseblocks_end_phasesets, inter_phaseblocks_hp_values = inter_phaseblocks_segments(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets)

    for i, (inter_phaseblocks_start, inter_phaseblocks_end, inter_phaseblocks_hp) in enumerate(zip(inter_phaseblocks_start_phasesets, inter_phaseblocks_end_phasesets, inter_phaseblocks_hp_values)):
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
                elif inter_phaseblocks_hp == 2:
                    if haplotype_1_values[i] > haplotype_2_values[i]:
                        new_hp2 = haplotype_2_values[i]
                        new_hp1 = haplotype_1_values[i]
                        haplotype_1_values[i] = new_hp2
                        haplotype_2_values[i] = new_hp1
                i = i +1

    return haplotype_1_values, haplotype_2_values

def fix_inter_cn_phase_switch_errors(args, df_segs_hp_1_, df_segs_hp_2_, df_hp_1_, df_hp_2_):
    chroms = get_contigs_list(args.contigs)
    df_segs_hp_1 = df_segs_hp_1_.copy()
    df_segs_hp_2 = df_segs_hp_2_.copy()
    df_hp_1 = df_hp_1_.copy()
    df_hp_2 = df_hp_2_.copy()

    updated_df_segs_hp1 = []
    updated_df_segs_hp2 = []
    updated_df_hp_1 = []
    updated_df_hp_2 = []

    for index, chrom in enumerate(chroms):
        df_seg_hp_1 = df_segs_hp_1[df_segs_hp_1['chromosome'] == chrom]
        df_seg_hp_2 = df_segs_hp_2[df_segs_hp_2['chromosome'] == chrom]

        hp_1 = df_hp_1[df_hp_1['chr'] == chrom]
        hp_2 = df_hp_2[df_hp_2['chr'] == chrom]

        hp_1_hp1 = hp_1.hp1.values.tolist()
        hp_1_start = hp_1.start.values.tolist()

        hp_2_hp2 = hp_2.hp2.values.tolist()

        haplotype_1_state_copyrnumbers = df_seg_hp_1.state.values.tolist()
        haplotype_1_depth_copyrnumbers = df_seg_hp_1.state.values.tolist()
        haplotype_1_start_values_copyrnumbers = df_seg_hp_1.start.values.tolist()
        haplotype_1_end_values_copyrnumbers = df_seg_hp_1.end.values.tolist()

        haplotype_2_state_copyrnumbers = df_seg_hp_2.state.values.tolist()
        haplotype_2_depth_copyrnumbers = df_seg_hp_2.depth.values.tolist()
        haplotype_2_start_values_copyrnumbers = df_seg_hp_2.start.values.tolist()
        haplotype_2_end_values_copyrnumbers = df_seg_hp_2.end.values.tolist()

        matching_indices_hp1 = []
        matching_indices_hp2 = []
        for i, (start_a, end_a) in enumerate(zip(haplotype_1_start_values_copyrnumbers, haplotype_1_end_values_copyrnumbers)):
            if start_a in haplotype_2_start_values_copyrnumbers and end_a in haplotype_2_end_values_copyrnumbers:
                matching_indices_hp1.append(i)
                matching_indices_hp2.append(haplotype_2_start_values_copyrnumbers.index(start_a))

        for index, (index_hp1, index_hp2) in enumerate(zip(matching_indices_hp1, matching_indices_hp2)):
            if chrom == 'chr10' and ((index_hp1 == 0 and index_hp2 == 0) or (index_hp1 == len(haplotype_1_start_values_copyrnumbers)-1 and index_hp2 == len(haplotype_2_start_values_copyrnumbers)-1)):
                    #CN states change
                    new_state_hp2 = haplotype_2_state_copyrnumbers[index_hp2]
                    new_state_hp1 = haplotype_1_state_copyrnumbers[index_hp1]
                    haplotype_1_state_copyrnumbers[index_hp1] = new_state_hp2
                    haplotype_2_state_copyrnumbers[index_hp2] = new_state_hp1

                    new_depth_hp2 = haplotype_2_depth_copyrnumbers[index_hp2]
                    new_depth_hp1 = haplotype_1_depth_copyrnumbers[index_hp1]
                    haplotype_1_depth_copyrnumbers[index_hp1] = new_depth_hp2
                    haplotype_2_depth_copyrnumbers[index_hp2] = new_depth_hp1

                    #bins change
                    internal_bins = [k for k in hp_1_start if k >= haplotype_2_start_values_copyrnumbers[index_hp2] and k <= haplotype_2_end_values_copyrnumbers[index_hp2]]
                    if internal_bins:
                        l = hp_1_start.index(internal_bins[0])
                        for j in range(len(internal_bins)):
                            new_hp2 = hp_2_hp2[l]
                            new_hp1 = hp_1_hp1[l]
                            hp_1_hp1[l] = new_hp2
                            hp_2_hp2[l] = new_hp1
                            l = l + 1

        updated_df_hp_1.append(pd.DataFrame(list(zip(hp_1.chr.values.tolist(), hp_1.start.values.tolist(),
                                  hp_1.end.values.tolist(), hp_1_hp1)),
                         columns=['chr', 'start', 'end', 'hp1']))
        updated_df_hp_2.append(pd.DataFrame(list(zip(hp_2.chr.values.tolist(), hp_2.start.values.tolist(),
                                  hp_2.end.values.tolist(), hp_2_hp2)),
                         columns=['chr', 'start', 'end', 'hp2']))

        updated_df_segs_hp1.append(pd.DataFrame(list(zip(df_seg_hp_1.chromosome.values.tolist(), df_seg_hp_1.start.values.tolist(),
                                  df_seg_hp_1.end.values.tolist(), haplotype_1_depth_copyrnumbers, haplotype_1_state_copyrnumbers)),
                         columns=['chromosome', 'start', 'end', 'depth', 'state']))
        updated_df_segs_hp2.append(pd.DataFrame(list(zip(df_seg_hp_2.chromosome.values.tolist(), df_seg_hp_2.start.values.tolist(),
                                  df_seg_hp_2.end.values.tolist(), haplotype_2_depth_copyrnumbers, haplotype_2_state_copyrnumbers)),
                         columns=['chromosome', 'start', 'end', 'depth', 'state']))

    return pd.concat(updated_df_segs_hp1), pd.concat(updated_df_segs_hp2), pd.concat(updated_df_hp_1), pd.concat(updated_df_hp_2)

def check_cn_state_phaseset(phaseset, df_segs_hp1_updated_chrom, df_segs_hp2_updated_chrom):
    hp1_state = -1
    hp2_state = -1
    df_hp1_state = df_segs_hp1_updated_chrom.depth.values.tolist()
    df_hp2_state = df_segs_hp2_updated_chrom.depth.values.tolist()

    for i, (start, end) in enumerate(zip(df_segs_hp1_updated_chrom.start.values.tolist(), df_segs_hp1_updated_chrom.end.values.tolist())):
        if start >= phaseset['start'] and end <= phaseset['end']:
            hp1_state = df_hp1_state[i]

    for i, (start, end) in enumerate(zip(df_segs_hp2_updated_chrom.start.values.tolist(), df_segs_hp2_updated_chrom.end.values.tolist())):
        if start >= phaseset['start'] and end <= phaseset['end']:
            hp2_state = df_hp2_state[i]

    return hp1_state, hp2_state

def bins_correction_phaseblocks(args, csv_df_phasesets, df_segs_hp1_updated, df_segs_hp2_updated, df_hp1, df_hp2):

    csv_df_phasesets = csv_df_chromosomes_sorter(args.out_dir_plots+'/data_phasing/' + args.genome_name + '_phasesets.csv', ['chr', 'start', 'end'])

    updated_df_hp_1 = []
    updated_df_hp_2 = []
    chroms = get_contigs_list(args.contigs)
    for i, chrom in enumerate(chroms):
        df_phasesets_chrom = csv_df_phasesets[csv_df_phasesets['chr'] == chrom]
        df_hp1_chrom = df_hp1[df_hp1['chr'] == chrom]
        df_hp2_chrom = df_hp2[df_hp2['chr'] == chrom]

        df_segs_hp1_updated_chrom = df_segs_hp1_updated[df_segs_hp1_updated['chromosome'] == chrom]
        df_segs_hp2_updated_chrom = df_segs_hp2_updated[df_segs_hp2_updated['chromosome'] == chrom]

        hp_1_hp1 = df_hp1_chrom.hp1.values.tolist()
        hp_2_hp2 = df_hp2_chrom.hp2.values.tolist()

        hp_1_start = df_hp1_chrom.start.values.tolist()

        for index, phaseset in df_phasesets_chrom.iterrows():
            hp1_mean_ps, hp2_mean_ps = check_cn_state_phaseset(phaseset, df_segs_hp1_updated_chrom, df_segs_hp2_updated_chrom)
            hp1_mean = statistics.mean(remove_outliers_iqr(np.array(hp_1_hp1[phaseset['start'] // args.bin_size:phaseset['end'] // args.bin_size])))
            hp2_mean = statistics.mean(remove_outliers_iqr(np.array(hp_2_hp2[phaseset['start'] // args.bin_size:phaseset['end'] // args.bin_size])))
            if ((abs(hp1_mean_ps - hp1_mean) > abs(hp1_mean_ps - hp2_mean)) and (not hp1_mean_ps == -1 and not hp2_mean_ps == -1)):
            #if phaseset['chr'] == 'chr18' and phaseset['start'] == 67800001 and phaseset['end'] == 70650001 or \
            #   phaseset['chr'] == 'chr17' and phaseset['start'] == 18532178 and phaseset['end'] == 21159334:

                # bins change
                internal_bins = [k for k in hp_1_start if k >= phaseset['start'] and k <= phaseset['end']]
                if internal_bins:
                    l = hp_1_start.index(internal_bins[0])
                    for j in range(len(internal_bins)):
                        new_hp2 = hp_2_hp2[l]
                        new_hp1 = hp_1_hp1[l]
                        hp_1_hp1[l] = new_hp2
                        hp_2_hp2[l] = new_hp1
                        l = l + 1

        updated_df_hp_1.append(pd.DataFrame(list(zip(df_hp1_chrom.chr.values.tolist(), df_hp1_chrom.start.values.tolist(),
                                  df_hp1_chrom.end.values.tolist(), hp_1_hp1)),
                         columns=['chr', 'start', 'end', 'hp1']))
        updated_df_hp_2.append(pd.DataFrame(list(zip(df_hp2_chrom.chr.values.tolist(), df_hp2_chrom.start.values.tolist(),
                                  df_hp2_chrom.end.values.tolist(), hp_2_hp2)),
                         columns=['chr', 'start', 'end', 'hp2']))

    return pd.concat(updated_df_hp_1), pd.concat(updated_df_hp_2)