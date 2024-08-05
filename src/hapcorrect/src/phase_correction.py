import os
import statistics
import numpy as np
import pandas as pd
from hapcorrect.src.utils import csv_df_chromosomes_sorter, write_segments_coverage, merge_regions, get_contigs_list
from hapcorrect.src.process_vcf import squash_regions

def generate_phasesets_bins(bam, path, bin_size, arguments):
    return get_phasesets_bins(bam, path, bin_size, arguments)

def get_phasesets_bins(bam, phasesets, bin_size, arguments):
    indices, values = remove_overlapping_and_small_phasesets(phasesets, bin_size, arguments)
    head, tail = os.path.split(bam)
    bed = []
    chroms = get_contigs_list(arguments['contigs'])
    for ind, chrom in enumerate(chroms) :
        for i in range(0, len(values[ind]), 2): #for i in range(len(values[ind])-1):
            start = values[ind][i]
            end = values[ind][i+1]
            if end - start > arguments['bin_size']*6:
                bed.append([tail, chrom, start, end])
    return bed

def remove_overlapping_and_small_phasesets(phasesets, bin_size, arguments):
    #dfs = pd.read_csv(phasesets, sep='\t', names=['chr', 'pos', 'ps'])
    dfs = csv_df_chromosomes_sorter(phasesets, ['chr', 'pos', 'ps'])
    values_all = []
    indices_all = []
    chroms = get_contigs_list(arguments['contigs'])
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

# def inside_phaseblock_switch_correction_coverage_update(chrom, arguments, is_simple_correction, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
#                     haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):



def phaseblock_flipping(chrom, arguments, is_simple_correction, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):

    values_ps = []
    for index, value in enumerate(ref_start_values_phasesets):
        values_ps.append([ref_start_values_phasesets[index], ref_end_values_phasesets[index]])

    #inside phaseblocks phaseswitch errors
    scan_and_update_phaseswitches_inside_phaseblocks(values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                                                     haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
    if is_simple_correction == False and len(ref_start_values_phasesets):
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
        scan_and_update_phaseblocks_switch_errors(chrom, arguments, max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                        haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

        #case for whole HPs switch, ie., HCC1437 - chr8
        #phase_switch_spanning_haplotypes(chrom, arguments, max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
        #                haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

        #bins without phaseblocks
        bins_without_phaseblocks(chrom, arguments, max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
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
        #if value_ps[0] == 37112432 :#and value_ps[1] == 143633540:# or value_ps[0] == 72108881:#72346221:#132135757:#177448845:#81331913:#72346221:
        #    print("here")
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
            if len(mean_value_haplotype_1) > 1:
                if inside_phasesets_alternate_hps_threshold_check(mean_value_haplotype_1, mean_value_haplotype_2, ref_start_values, value_ps, haplotype_1_values, haplotype_2_values):
                    start_index = ref_start_values_phasesets.index(value_ps[0])
                    ref_start_values_phasesets.pop(start_index)
                    ref_end_values_phasesets.pop(start_index)
                    haplotype_1_values_phasesets.pop(start_index)
                    haplotype_2_values_phasesets.pop(start_index)

                    ref_start_values_phasesets[start_index:start_index] = first_value
                    ref_end_values_phasesets[start_index:start_index] = last_value
                    haplotype_1_values_phasesets[start_index:start_index] = mean_value_haplotype_1
                    haplotype_2_values_phasesets[start_index:start_index] = mean_value_haplotype_2



def scan_and_update_phaseblocks_switch_errors(chrom, arguments, max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):

    PHASESETS_DIFF_THRESHOLD = 12
    #############################################################################################
    # Scan Right-Left (HP1>HP2)
    if haplotype_1_values_phasesets[max_value_index] > haplotype_2_values_phasesets[max_value_index] or \
            haplotype_1_values_phasesets[max_value_index] == haplotype_2_values_phasesets[max_value_index]:
        values_ps_new = values_ps[max_value_index + 1:]
        if not max_value_index == len(ref_start_values_phasesets) - 1:
            for j, value_ps in enumerate(values_ps_new):
                #if value_ps[0] == 121986191:  # 101790293 or value_ps[0] == 104622049:
                #    print("here")
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
                            #if chrom == 'chr3' and 106100001 in internal_bins:
                            #    print('here')
                            dict.append((chrom + '\t' + str(value_ps[0]) + '\t' + str(value_ps[1])))
                            write_segments_coverage(dict, arguments['genome_name'] + '_phase_change_segments.csv')
                            #if chrom == 'chr2' and internal_bins[0] > internal_bins[-1] + arguments['bin_size']-1:
                            #    print('here')
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
                #if value_ps[0] == 121986191:  # 101790293 or value_ps[0] == 104622049:
                #    print("here")
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
                            #if chrom == 'chr3' and 106100001 in internal_bins:
                            #    print('here')
                            dict.append((chrom + '\t' + str(value_ps[0]) + '\t' + str(value_ps[1])))
                            write_segments_coverage(dict, arguments['genome_name'] + '_phase_change_segments.csv')
                            #if chrom == 'chr2' and internal_bins[0] > internal_bins[-1] + arguments['bin_size']-1:
                            #    print('here')
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
                #if value_ps[0] == 121986191:  # 101790293 or value_ps[0] == 104622049:
                #    print("here")

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
                            #if chrom == 'chr3' and 106100001 in internal_bins:
                            #    print('here')
                            dict.append((chrom + '\t' + str(value_ps[0]) + '\t' + str(value_ps[1] )))
                            write_segments_coverage(dict, arguments['genome_name'] + '_phase_change_segments.csv')
                            #if chrom == 'chr2' and internal_bins[0] > internal_bins[-1] + arguments['bin_size']-1:
                            #    print('here')
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
                #if value_ps[0] == 121986191:  # 101790293 or value_ps[0] == 104622049:
                #    print("here")
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
                            dict.append((chrom + '\t' + str(value_ps[0]) + '\t' + str(value_ps[1])))
                            write_segments_coverage(dict, arguments['genome_name'] + '_phase_change_segments.csv')
                            #if chrom == 'chr2' and internal_bins[0] > internal_bins[-1] + arguments['bin_size']-1:
                            #    print('here')
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

def bins_without_phaseblocks(chrom, arguments, max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
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
    #                 (chrom + '\t' + str(internal_bins[0]) + '\t' + str(internal_bins[-1] + arguments['bin_size'] - 1)))
    #             write_segments_coverage(dict, arguments['genome_name'] + '_phase_change_segments.csv')
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
    #                 (chrom + '\t' + str(internal_bins[0]) + '\t' + str(internal_bins[-1] + arguments['bin_size'] - 1)))
    #             write_segments_coverage(dict, arguments['genome_name'] + '_phase_change_segments.csv')
    #             for k in range(len(internal_bins)):
    #                 if haplotype_1_values[i] > haplotype_2_values[i]:
    #                     new_hp2 = haplotype_2_values[i]
    #                     new_hp1 = haplotype_1_values[i]
    #                     haplotype_1_values[i] = new_hp2
    #                     haplotype_2_values[i] = new_hp1
    #                 i = i + 1


def phase_switch_spanning_haplotypes(chrom, arguments, max_value_index, values_ps, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
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

def phase_correction_centers(arguments, segs_hp1_state, segs_hp1_start, segs_hp1_end,  segs_hp2_state, segs_hp2_start, segs_hp2_end, haplotype_1_values, haplotype_2_values):

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
                for i in range(start//arguments['bin_size'],end//arguments['bin_size']):
                    new_hp2 = haplotype_2_values[i]
                    new_hp1 = haplotype_1_values[i]
                    haplotype_1_values[i] = new_hp2
                    haplotype_2_values[i] = new_hp1

        for index, (start, end, state) in enumerate(zip(segs_hp1_start, segs_hp1_end, segs_hp1_state)):
            if state == max_diff_hp2[2]:
                for i in range(start//arguments['bin_size'],end//arguments['bin_size']):
                    new_hp2 = haplotype_2_values[i]
                    new_hp1 = haplotype_1_values[i]
                    haplotype_1_values[i] = new_hp2
                    haplotype_2_values[i] = new_hp1

    return haplotype_1_values, haplotype_2_values, processed

def contiguous_phaseblocks(haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, loh_starts, loh_ends):
    # THRESHOLD = 8
    # contiguous_hp1_start = []
    # contiguous_hp1_end = []
    # hp1_values = []
    # contiguous_hp1_start.append(ref_start_values_phasesets[0])
    # hp1_value = []
    # for i in range(len(ref_start_values_phasesets) - 1):
    #     current_value_hp1 = haplotype_1_values_phasesets[i]
    #     hp1_value.append(haplotype_1_values_phasesets[i])
    #     if not haplotype_1_values_phasesets[i+1] - THRESHOLD < current_value_hp1 < haplotype_1_values_phasesets[i+1] + THRESHOLD:
    #         contiguous_hp1_start.append(ref_start_values_phasesets[i+1])
    #         contiguous_hp1_end.append(ref_start_values_phasesets[i+1]-1)
    #         hp1_values.append(statistics.mean(hp1_value))
    #         hp1_value.clear()
    # contiguous_hp1_end.append(ref_end_values_phasesets[-1])
    # if len(contiguous_hp1_start) > 2:
    #     hp1_values.append(statistics.mean(haplotype_1_values_phasesets[ref_start_values_phasesets.index(contiguous_hp1_start[-2]):ref_start_values_phasesets.index(ref_start_values_phasesets[-1])]))
    # else:
    #     hp1_values.append(0)
    #
    # contiguous_hp2_start = []
    # contiguous_hp2_end = []
    # hp2_values = []
    # contiguous_hp2_start.append(ref_start_values_phasesets[0])
    # hp2_value = []
    # for i in range(len(ref_start_values_phasesets) - 1):
    #     current_value_hp2 = haplotype_2_values_phasesets[i]
    #     hp2_value.append(haplotype_2_values_phasesets[i])
    #     if not haplotype_2_values_phasesets[i+1] - THRESHOLD < current_value_hp2 < haplotype_2_values_phasesets[i+1] + THRESHOLD :
    #         contiguous_hp2_start.append(ref_start_values_phasesets[i+1])
    #         contiguous_hp2_end.append(ref_start_values_phasesets[i+1]-1)
    #         hp2_values.append(statistics.mean(hp2_value))
    #         hp2_value.clear()
    # contiguous_hp2_end.append(ref_end_values_phasesets[-1])
    # if len(contiguous_hp2_start) > 2:
    #     hp2_values.append(statistics.mean(haplotype_2_values_phasesets[ref_start_values_phasesets.index(contiguous_hp2_start[-2]):ref_start_values_phasesets.index(ref_start_values_phasesets[-1])]))
    # else:
    #     hp2_values.append(0)


    ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = merge_regions(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, loh_starts, loh_ends)
    #contiguous_hp2_start, contiguous_hp2_end, hp2_values = merge_regions(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values_phasesets, loh_starts, loh_ends)
    return ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets
    #return hp1_values, contiguous_hp1_start, contiguous_hp1_end, hp2_values, contiguous_hp2_start, contiguous_hp2_end

def flip_phaseblocks_contigous(chrom, arguments, ref_start_values, haplotype_1_values_phasesets,
                               haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values, haplotype_2_values):
    THRESHOLD = 12
    diff_hp1 = []
    max_value = []
    for index, value in enumerate(ref_start_values_phasesets):
        diff_hp1.append(ref_end_values_phasesets[index] - ref_start_values_phasesets[index])
    if len(diff_hp1) > 1:
        for i in range(1):
            max_value.append(haplotype_1_values_phasesets[diff_hp1.index(sorted(diff_hp1, reverse=True)[i])])#diff_hp1.index(max(diff_hp1))
            print(max_value)
        # max_value = sorted(max_value, reverse=True)
        # if max_value[0] - max_value[1] < 8:
        #     max_value.pop(1)
    else:
        if diff_hp1:
            max_value.append(haplotype_1_values_phasesets[diff_hp1.index(diff_hp1[0])])
    print(max_value)
    value_ps_updated = []
    if max_value:
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
                        dict.append((chrom + '\t' + str(value_ps[0]) + '\t' + str(value_ps[1])))
                        write_segments_coverage(dict, arguments['genome_name'] + '_phase_change_segments.csv')
                        #if chrom == 'chr2' and internal_bins[0] > internal_bins[-1] + arguments['bin_size'] - 1:
                        #    print('here')
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
    for index, value in enumerate(ref_start_values_phasesets):
        diff_hp2.append(ref_end_values_phasesets[index] - ref_start_values_phasesets[index])
    if len(diff_hp2) > 1:
        for i in range(1):
            max_value.append(haplotype_2_values_phasesets[diff_hp2.index(sorted(diff_hp2, reverse=True)[i])])#diff_hp1.index(max(diff_hp1))
            print(max_value)
        # max_value = sorted(max_value)
        # if max_value[1] - max_value[0] < 8:
        #     max_value.pop(1)
    else:
        if diff_hp2:
            max_value.append(haplotype_2_values_phasesets[diff_hp2.index(diff_hp2[0])])
    if max_value:
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
                        dict.append((chrom + '\t' + str(value_ps[0]) + '\t' + str(value_ps[1])))
                        write_segments_coverage(dict, arguments['genome_name'] + '_phase_change_segments.csv')
                        #if chrom == 'chr2' and internal_bins[0] > internal_bins[-1] + arguments['bin_size'] - 1:
                        #    print('here')
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

def remove_overlaping_contiguous(chrom, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2):
    overlaping_start_hp2 = []
    overlaping_end_hp2 = []
    overlaping_start_hp1 = []
    overlaping_end_hp1 = []
    for i, (start_hp1,end_hp1) in enumerate(zip(ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1)):
        for j, (start_hp2, end_hp2) in enumerate(zip(ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2)):
            #if start_hp2 == 97490471 or start_hp1 == 97490471:
            #    print('here')
            if (start_hp2 >= start_hp1 and end_hp2 <= end_hp1) and not (start_hp2 == start_hp1 and end_hp2 == end_hp1) :
                overlaping_start_hp2.append(start_hp2)
                overlaping_end_hp2.append(end_hp2)

    for i, (start_hp2,end_hp2) in enumerate(zip(ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2)):
        for j, (start_hp1, end_hp1) in enumerate(zip(ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1)):
            #if start_hp2 == 97490471 or start_hp1 == 97490471:
            #    print('here')
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
    remove_overlaping_small_segments = []
    for i, (start_hp1,end_hp1) in enumerate(zip(ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1)):
        for j, (start_hp2, end_hp2) in enumerate(zip(ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2)):
            if start_hp2 < start_hp1 < end_hp2 and  end_hp1 > end_hp2 :
                if end_hp1 - start_hp1 > end_hp2 - start_hp2:
                    remove_overlaping_small_segments.append(start_hp2)
                else:
                    remove_overlaping_small_segments.append(start_hp1)
                #print("First: "+str(start_hp1)+':'+str(end_hp1) +' , '+ str(start_hp2)+':'+str(end_hp2))
            if start_hp1 < start_hp2 < end_hp1 and  end_hp2 > end_hp1 :
                if end_hp1 - start_hp1 > end_hp2 - start_hp2:
                    remove_overlaping_small_segments.append(start_hp2)
                else:
                    remove_overlaping_small_segments.append(start_hp1)
                #print("Second: "+str(start_hp1)+':'+ str(end_hp1) +' , '+ str(start_hp2)+':'+str(end_hp2))
    final_start_positions = [j for j in sorted(ref_start_values_phasesets_hp1 + ref_start_values_phasesets_hp2) if j not in remove_overlaping_small_segments]

    start_values_phasesets = np.unique(sorted(final_start_positions)).tolist()
    print(start_values_phasesets)
    chr_list = [chrom for ch in range(len(start_values_phasesets))]
    df_phasesets_chr = pd.DataFrame(list(zip(chr_list, start_values_phasesets)), columns=['chr', 'start'])

    return df_phasesets_chr

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

    if abs(mean_values(temp_hp1_values, 0, len(temp_hp1_values)-1) - mean_values(temp_hp2_values, 0, len(temp_hp1_values)-1)) > 9:
        return True
    else:
        return False

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

def switch_inter_phaseblocks_bins(chrom, arguments, ref_start_values, haplotype_1_values, haplotype_2_values, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets):
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
                        dict = []
                        dict.append((chrom + '\t' + str(bin) + '\t' + str(bin + arguments['bin_size'] - 1)))
                        write_segments_coverage(dict, arguments['genome_name'] + '_phase_change_segments.csv')
                elif inter_phaseblocks_hp == 2:
                    if haplotype_1_values[i] > haplotype_2_values[i]:
                        new_hp2 = haplotype_2_values[i]
                        new_hp1 = haplotype_1_values[i]
                        haplotype_1_values[i] = new_hp2
                        haplotype_2_values[i] = new_hp1
                        dict = []
                        dict.append((chrom + '\t' + str(bin) + '\t' + str(bin + arguments['bin_size'] - 1)))
                        write_segments_coverage(dict, arguments['genome_name'] + '_phase_change_segments.csv')
                i = i +1

    return haplotype_1_values, haplotype_2_values
