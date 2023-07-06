import statistics

import pandas as pd
import numpy as np
import os

def get_phasesets_bins(bam, phasesets, bin_size, arguments):
    indices, values = remove_overlapping_and_small_phasesets(phasesets, bin_size, arguments)
    head, tail = os.path.split(bam)
    bed = []
    from utils import get_contigs_list
    chroms = get_contigs_list(arguments['contigs'])
    for ind, chrom in enumerate(chroms) :
        for i in range(0, len(values[ind]), 2): #for i in range(len(values[ind])-1):
            start = values[ind][i]
            end = values[ind][i+1]
            if end - start > arguments['bin_size']*6:
                bed.append([tail, chrom, start, end])
    return bed

def remove_overlapping_and_small_phasesets(phasesets, bin_size, arguments):
    dfs = pd.read_csv(phasesets, sep='\t', names=['chr', 'pos', 'qual', 'filter', 'ps', 'gt', 'dp', 'vaf'])
    values_all = []
    indices_all = []
    from utils import get_contigs_list
    chroms = get_contigs_list(arguments['contigs'])
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

def phaseblock_flipping(is_simple_correction, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):
    if is_simple_correction == False:
        values_ps = []
        for index, value in enumerate(ref_start_values_phasesets):
            values_ps.append([ref_start_values_phasesets[index], ref_end_values_phasesets[index]])

        #inside phaseblocks phaseswitch errors
        scan_and_update_phaseswitches_inside_phaseblocks(values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                                                         haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
        diff = []
        values_ps = []
        #ref_start_values_phasesets.sort()
        #ref_end_values_phasesets.sort()

        for index, value in enumerate(ref_start_values_phasesets):
            diff.append(ref_end_values_phasesets[index] - ref_start_values_phasesets[index])
            values_ps.append([ref_start_values_phasesets[index], ref_end_values_phasesets[index]])
        max_value_index = diff.index(max(diff))
        print(ref_start_values_phasesets[max_value_index], ref_end_values_phasesets[max_value_index])

        #phaseblocks switch errors
        scan_and_update_phaseblocks_switch_errors(max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                        haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

        #case for whole HPs switch, ie., HCC1437 - chr8
        phase_switch_spanning_haplotypes(max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                        haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

        #bins without phaseblocks
        bins_without_phaseblocks(max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
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
        if value_ps[0] == 42066420 or value_ps[0] == 72346221 or value_ps[0] == 72108881:#72346221:#132135757:#177448845:#81331913:#72346221:
            print("here")
        print(value_ps[0])

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

def scan_and_update_phaseblocks_switch_errors(max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):

    PHASESETS_DIFF_THRESHOLD = 12
    #############################################################################################
    # Scan Right-Left (HP1>HP2)
    if haplotype_1_values_phasesets[max_value_index] > haplotype_2_values_phasesets[max_value_index] or \
            haplotype_1_values_phasesets[max_value_index] == haplotype_2_values_phasesets[max_value_index]:
        values_ps_new = values_ps[max_value_index + 1:]
        if not max_value_index == len(ref_start_values_phasesets) - 1:
            for j, value_ps in enumerate(values_ps_new):
                if value_ps[0] == 93800001:  # 101790293 or value_ps[0] == 104622049:
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
                                    haplotype_1_values_phasesets[
                                        ref_start_values_phasesets.index(value_ps[0]) - 1]) < PHASESETS_DIFF_THRESHOLD)) \
                                and haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] < \
                                haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]:
                            # or haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] > haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]: #3-5 bins check, both  sides of phaseblock
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
                if value_ps[0] == 93800001:  # 101790293 or value_ps[0] == 104622049:
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
                if value_ps[0] == 93800001:  # 101790293:
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
                                and (haplotype_1_values[i] > haplotype_2_values[i] or haplotype_1_values[i + 1] >
                                     haplotype_2_values[i + 1] or haplotype_1_values[i + 2] > haplotype_2_values[
                                         i + 2] or haplotype_1_values[i + 3] > haplotype_2_values[i + 3] or
                                     haplotype_1_values[i + 4] > haplotype_2_values[i + 4]))) \
                            or (abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                                    haplotype_1_values_phasesets[
                                        ref_start_values_phasesets.index(value_ps[0]) - 1]) < PHASESETS_DIFF_THRESHOLD)) \
                                and haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] > \
                                haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]:
                            # or haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] > haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])]: #3-5 bins check, both  sides of phaseblock
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
                if value_ps[0] == 93800001:  # 101790293:
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
                            for k in range(len(internal_bins)):
                                # if haplotype_1_values[i] < haplotype_2_values[i]:
                                new_hp2 = haplotype_2_values[i]
                                new_hp1 = haplotype_1_values[i]
                                haplotype_1_values[i] = new_hp2
                                haplotype_2_values[i] = new_hp1
                                i = i + 1

def bins_without_phaseblocks(max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):
    if haplotype_1_values_phasesets[max_value_index] > haplotype_2_values_phasesets[max_value_index] or haplotype_1_values_phasesets[max_value_index] == haplotype_2_values_phasesets[max_value_index]:
        for j in range(len(values_ps)-1):
            internal_bins = [i for i in ref_start_values if i > values_ps[j][1] and i < values_ps[j+1][0]]
            if internal_bins:
                i = ref_start_values.index(internal_bins[0])
                for k in range(len(internal_bins)):
                    if haplotype_1_values[i] < haplotype_2_values[i]:
                        new_hp2 = haplotype_2_values[i]
                        new_hp1 = haplotype_1_values[i]
                        haplotype_1_values[i] = new_hp2
                        haplotype_2_values[i] = new_hp1
                    i = i + 1
    else:
        for j in range(len(values_ps) - 1):
            internal_bins = [i for i in ref_start_values if i > values_ps[j][1] and i < values_ps[j + 1][0]]
            if internal_bins:
                i = ref_start_values.index(internal_bins[0])
                for k in range(len(internal_bins)):
                    if haplotype_1_values[i] > haplotype_2_values[i]:
                        new_hp2 = haplotype_2_values[i]
                        new_hp1 = haplotype_1_values[i]
                        haplotype_1_values[i] = new_hp2
                        haplotype_2_values[i] = new_hp1
                    i = i + 1


def phase_switch_spanning_haplotypes(max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):
    MIN_PHASESETS_DIFF_THRESHOLD = 3
    for index, value_ps in enumerate(values_ps):
        if value_ps[0] == 93800001 and value_ps[1] == 66901398:
            print("here")
        if not ref_end_values_phasesets.index(value_ps[1]) == len(ref_start_values_phasesets)-1:
            if (abs(haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                   haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+1])  < MIN_PHASESETS_DIFF_THRESHOLD and \
                abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
                    haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+1]) > 25) and \
                abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])+2] - \
                        haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0]) + 1]) < MIN_PHASESETS_DIFF_THRESHOLD and \
                abs(haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0]) + 2] - \
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


#Apmlification cases
# or (abs(haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
#         haplotype_1_values_phasesets[ref_start_values_phasesets.index(value_ps[0]) + 1] > PHASESETS_DIFF_THRESHOLD) \
#     and abs(haplotype_2_values_phasesets[ref_start_values_phasesets.index(value_ps[0])] - \
#             haplotype_2_values_phasesets[
#                 ref_start_values_phasesets.index(value_ps[0]) + 1] < PHASESETS_DIFF_THRESHOLD // 2)):