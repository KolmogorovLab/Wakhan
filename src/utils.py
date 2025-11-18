import csv
import numpy as np
import pandas as pd
import pysam
import scipy.stats as stats
import os
import shutil
import math
import statistics
import ruptures as rpt
import logging
from itertools import takewhile
from joblib import Parallel, delayed
from vcf_parser import VCFParser

logger = logging.getLogger()

from src.coverage.smoothing import smoothing
from src.breakpoint.breakpoints import sv_vcf_bps_cn_check

from src.coverage.binning import update_bins_with_bps, update_bins_with_bps_new
from src.coverage.segmentation import merge_adjacent_regions_loh, remove_overlapping_segments
from src.cna.copynumber import subclonal_values_adjusted, integers_values_adjusted, add_confidence_score_cn_segemnts
from src.utils_tmp.chromosome import csv_df_chromosomes_sorter, df_chromosomes_sorter, get_contigs_list



def write_segments_coverage(coverage_segments, output, args):
    with open(args.out_dir_plots + '/bed_output/' + output, 'a') as fp:
        for items in coverage_segments:
            if not items == None:
                fp.write("%s\n" % items)


def detect_alter_loh_regions(args, event, chrom, ref_ends, haplotype_1_values, haplotype_2_values, unphased_reads_values, starts, ends, switch_hps):
    if ends and ends[-1] > ref_ends[-1]:
        ends[-1] = ref_ends[-1]

    region_starts = []
    region_ends = []
    #print(starts)
    #print(ends)
    for i, (start,end) in enumerate(zip(starts,ends)):
        if end - start > args.hets_loh_seg_size:
            region_starts.append(start)
            region_ends.append(end)

    #print(region_starts)
    #print(region_ends)

    if not args.without_phasing and switch_hps:
        for j, (starts,ends) in enumerate(zip(region_starts, region_ends)):
            #TODO Discuss with Ayse, alternate approach on what HP should be selected for each region
            #if mean_values(haplotype_1_values, starts - 1, starts - 4) > mean_values(haplotype_2_values, starts - 1, starts - 4):
            for i in range(starts//args.bin_size,ends//args.bin_size):
                    haplotype_1_values[i] = haplotype_1_values[i] + haplotype_2_values[i] + unphased_reads_values[i]
                    haplotype_2_values[i] = 0
                    unphased_reads_values[i] = 0
            # else:
            #     for i in range(starts // 50000, ends // 50000):
            #         haplotype_2_values[i] = haplotype_1_values[i] + haplotype_2_values[i] + unphased_reads_values[i]
            #         haplotype_1_values[i] = 0
            #         unphased_reads_values[i] = 0

    return haplotype_1_values, haplotype_2_values, unphased_reads_values, region_starts, region_ends


def loh_regions_phasesets(loh_region_starts, loh_region_ends, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):
    indices = []
    for l, (loh_start, loh_end) in enumerate(zip(loh_region_starts, loh_region_ends)):
        for k, (ps_start, ps_end) in enumerate(zip(ref_start_values_phasesets, ref_end_values_phasesets)):
            if ps_start >= loh_start and ps_end <= loh_end:
                indices.append(k)

    haplotype_1_values_phasesets = [j for i, j in enumerate(haplotype_1_values_phasesets) if i not in indices]
    haplotype_2_values_phasesets = [j for i, j in enumerate(haplotype_2_values_phasesets) if i not in indices]
    ref_start_values_phasesets = [j for i, j in enumerate(ref_start_values_phasesets) if i not in indices]
    ref_end_values_phasesets = [j for i, j in enumerate(ref_end_values_phasesets) if i not in indices]

    return haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets


def loh_regions_events(chrom, region_starts, region_ends, args):
    dict = []
    for i in range(len(region_starts)):
        dict.append((chrom + '\t' + str(region_starts[i]) + '\t' + str(region_ends[i])))
    #write_segments_coverage(dict, args.genome_name + '_loh_segments.bed')
    return dict


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


def collect_loh_centromere_regions(df_segs_hp1_, df_segs_hp2_, centers, integer_fractional_means, args):
    df_segs_hp1 = df_segs_hp1_.copy()
    df_segs_hp2 = df_segs_hp2_.copy()
    for i in range(len(integer_fractional_means)):
        df_segs_hp1['state'].mask(df_segs_hp1['state'] == centers[i], integer_fractional_means[i], inplace=True)
        df_segs_hp2['state'].mask(df_segs_hp2['state'] == centers[i], integer_fractional_means[i], inplace=True)

    chroms = get_contigs_list(args.contigs)
    loh_regions_all = []

    for index, chrom in enumerate(chroms):
        chrs_list = []
        loh_region_starts = []
        loh_region_ends = []

        df_segs_hp_1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
        df_segs_hp_2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]

        hp1_state = df_segs_hp_1_chrom.state.values.tolist()
        hp2_state = df_segs_hp_2_chrom.state.values.tolist()

        hp1_start = df_segs_hp_1_chrom.start.values.tolist()
        hp1_end = df_segs_hp_1_chrom.end.values.tolist()

        hp2_start = df_segs_hp_2_chrom.start.values.tolist()
        hp2_end = df_segs_hp_2_chrom.end.values.tolist()

        # if args.breakpoints and not args.cpd_internal_segments:
        #     for index, (state1, state2, start, end) in enumerate(zip(hp1_state, hp2_state, hp1_start, hp1_end)):
        #         #print(chrom, state1, state2, start, end)
        #         if state1 == 0 or state2 == 0 and end - start > 2000000:
        #             loh_region_starts.append(start)
        #             loh_region_ends.append(end)
        # else:
        for index, (state1, start, end) in enumerate(zip(hp1_state, hp1_start, hp1_end)):
            if state1 == 0 and end - start > 2000000:
                loh_region_starts.append(start)
                loh_region_ends.append(end)

        for index, (state2, start, end) in enumerate(zip(hp2_state, hp2_start, hp2_end)):
            if state2 == 0 and end - start > 2000000:
                if not start in loh_region_starts and not end in loh_region_ends:
                    loh_region_starts.append(start)
                    loh_region_ends.append(end)

        # sort
        if loh_region_starts:
            sort_function = lambda x: x[0]
            sort_target = list(zip(loh_region_starts, loh_region_ends))
            sort_target.sort(key=sort_function)
            loh_region_starts = [a for a, b in sort_target]
            loh_region_ends = [b for a, b in sort_target]

            # remove overlapping segments
            loh_region_starts, loh_region_ends = remove_overlapping_segments(loh_region_starts, loh_region_ends)

        if loh_region_starts:
            chrs_list.extend([chrom for ch in range(len(loh_region_starts))])
            loh_regions_all.append(pd.DataFrame(list(zip(chrs_list, loh_region_starts, loh_region_ends)),
                                                columns=['chr', 'start', 'end']))

    if loh_regions_all:
        return merge_adjacent_regions_loh(pd.concat(loh_regions_all), args)
    else:
        return pd.DataFrame(columns=['chr', 'start', 'end'])


#def extract_centromere_regions(args):
#    fileDir = os.path.dirname(__file__)
#    cen_coord = os.path.join(fileDir, args.centromere)
#    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
#    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)
#
#    return df_centm

# def breakpoints_segments_means(args, csv_df_snps_mean):
#     haplotype_1_segs_dfs = []
#     haplotype_2_segs_dfs = []
#
#     fileDir = os.path.dirname(__file__)  # os.path.dirname(os.path.realpath('__file__'))
#     cen_coord = os.path.join(fileDir, args.centromere)
#     df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
#     df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)
#
#     if args.quick_start:
#         loh_path = args.quick_start_coverage_path + '/'
#     else:
#         loh_path = args.out_dir_plots + '/coverage_data/'
#     if args.tumor_phased_vcf and os.path.exists(loh_path + args.genome_name + '_loh_segments.csv'):
#         df_loh = csv_df_chromosomes_sorter(loh_path + args.genome_name + '_loh_segments.csv', ['chr', 'start', 'end', 'hp'])
#     else:
#         df_loh = pd.DataFrame(columns=['chr', 'start', 'end', 'hp'])
#
#     if args.breakpoints:
#         _, _, _, breakpoints_segemnts, bps, bps_bnd = sv_vcf_bps_cn_check(args.breakpoints, args)
#         df_var_bins, df_var_bins_1 = get_chromosomes_bins(args.target_bam[0], args.bin_size, args)
#
#     chroms = get_contigs_list(args.contigs)
#     for index, chrom in enumerate(chroms):
#         df_centm_chrom = df_centm[df_centm['chr'] == chrom]
#         df_loh_chrom = df_loh[df_loh['chr'] == chrom]
#
#         if args.breakpoints:
#             df_var_bins_chr_1 = df_var_bins_1[df_var_bins_1['chr'] == chrom]
#             all_breakpoints = df_var_bins_chr_1.start.values.tolist()
#         else:
#             all_breakpoints = []
#             breakpoints_segemnts = []
#
#         snps_cpd_means, snps_cpd_lens, df_means_chr = change_point_detection_means(args, chrom, breakpoints_segemnts, csv_df_snps_mean, csv_df_snps_mean.start.tolist(), all_breakpoints, df_centm_chrom, df_loh_chrom)
#         haplotype_1_segs_dfs.append(df_means_chr[0])
#
#         snps_cpd_means, snps_cpd_lens, df_means_chr = change_point_detection_means(args, chrom, breakpoints_segemnts, csv_df_snps_mean, csv_df_snps_mean.start.tolist(), all_breakpoints, df_centm_chrom, df_loh_chrom)
#         haplotype_2_segs_dfs.append(df_means_chr[1])
#
#     return pd.concat(haplotype_1_segs_dfs), pd.concat(haplotype_2_segs_dfs)

# def chunk_range(start, end, num_chunks):
#
#   if num_chunks == 0:
#     raise ValueError("Number of chunks cannot be zero")
#
#   # Calculate the ideal chunk size (might have a remainder)
#   ideal_chunk_size = (end - start) // num_chunks
#
#   # Initialize list to store chunks
#   chunks = []
#   current_start = start
#   last = end
#
#   # Iterate and create chunks
#   for _ in range(num_chunks):
#     # Handle last chunk potentially being larger
#     if _ == num_chunks :
#       end_value = end
#     else:
#       end_value = current_start + ideal_chunk_size
#     if not last == end_value:
#         chunks.append([current_start, end_value-1])
#     else:
#         chunks.append([current_start, end_value])
#     current_start = end_value
#
#   return chunks


# def get_breakpoints(chrom, bp_file_path): #TODO add call in plots
#     break_points = []
#     with open(bp_file_path) as bp_file:
#         next(bp_file)
#         for st in bp_file:
#             # st = '-chr1:2671683|+chr1:2673127,0,0,0,9,10,22'
#             st = st.split(",")
#             #Parse the file
#             chr1 = (st[0].split("|"))[0].split(":")[0][1:]
#             chr2 = (st[0].split("|"))[1].split(":")[0][1:]
#             bp_pos1 = int((st[0].split("|"))[0].split(":")[1])
#             bp_pos2 = int((st[0].split("|"))[1].split(":")[1])
#             val_1=int(st[1])
#             val_2=int(st[4])
#             if ( val_1 == 0 and val_2 >= 3) and (abs(bp_pos1 - bp_pos2) > 1000) and (chr1 == chr2) and (chr2 == chrom):
#                 #1k condition and both connections should be on same chromosome for the moment
#                 mid=bp_pos1 + round((abs(bp_pos1 - bp_pos2) / 2))
#                 break_points.extend([bp_pos1, mid, bp_pos2])
#     return break_points
