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

#from src.coverage.smoothing import smoothing
#from src.breakpoint.breakpoints import sv_vcf_bps_cn_check

#from src.coverage.binning import update_bins_with_bps, update_bins_with_bps_new
#from src.coverage.segmentation import merge_adjacent_regions_loh, remove_overlapping_segments
#from src.cna.copynumber import subclonal_values_adjusted, integers_values_adjusted, add_confidence_score_cn_segemnts
#from src.utils_tmp.chromosome import csv_df_chromosomes_sorter, df_chromosomes_sorter, get_contigs_list

#from src.file_tools.process_bam import get_segments_coverage
#from src.coverage.segmentation import split_regions_by_points
#from src.utils_tmp.chromosome import csv_df_chromosomes_sorter, df_chromosomes_sorter, get_contigs_list


#def mean_values(selected_list, start_index, end_index):
#    result = []
#    for i in range(end_index - start_index):
#        try:
#            result.append(selected_list[start_index + i])
#        except IndexError:
#            break
#    if result:
#        return np.mean(result)
#    else:
#        return 0.0


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


#def loh_regions_events(chrom, region_starts, region_ends, hp):
#    dict = []
#    for i in range(len(region_starts)):
#        dict.append((chrom + '\t' + str(region_starts[i]) + '\t' + str(region_ends[i]) + '\t' + str(hp[i])))
#    #write_segments_coverage(dict, args.genome_name'] + '_loh_segments.bed')
#    return dict
