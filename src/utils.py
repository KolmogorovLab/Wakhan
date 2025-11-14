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

from src.smoothing import smoothing
from src.breakpoints import get_contigs_list, sv_vcf_bps_cn_check

from src.coverage.binning import update_bins_with_bps, update_bins_with_bps_new
from src.coverage.segmentation import merge_adjacent_regions_loh, remove_overlapping_segments
from src.utils_tmp.chromosome import csv_df_chromosomes_sorter, df_chromosomes_sorter


def adjust_extreme_outliers(hp_data):
    for j, val in enumerate(hp_data):
        if val > 500 and j > 0 and j < len(hp_data):
            hp_data[j] = hp_data[j-1]
        elif val > 500 and j == 0:
            hp_data[j] = hp_data[j+1]
    return hp_data


def get_chromosomes_bins(bam_file, bin_size, args):
    bed=[]
    bam_alignment = pysam.AlignmentFile(bam_file)
    headers = bam_alignment.header
    seq_dict = headers['SQ']
    region = [''] * len(seq_dict)
    chrs = [''] * len(seq_dict)
    head, tail = os.path.split(bam_file)
    chroms = get_contigs_list(args.contigs)
    chroms_without_prefix = [str(i).replace( 'chr', '') for i in chroms]
    for i, seq_elem in enumerate(seq_dict):
        region[i] = seq_elem['LN']
        chrs[i] = seq_elem['SN']
        start=0
        end=bin_size
        if (chrs[i] in chroms) or (chrs[i] in chroms_without_prefix):
            for c in range(0,region[i],bin_size):
                if end > region[i]:
                    bed.append([chrs[i], start, region[i]])
                else:
                    bed.append([chrs[i], start, end])
                start=end+1
                end+=bin_size

    _,_,_,_, bps, bps_bnd = sv_vcf_bps_cn_check(args.breakpoints, args)
    bed = update_bins_with_bps(bed, bps, bps_bnd, args, region)
    _, bed_1 = update_bins_with_bps_new(bed, bps, bps_bnd, args, region)

    return bed, bed_1

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


def write_df_csv(df, path):
    fp = open(path, 'a')
    fp.write('#chr: chromosome number\n')
    fp.write('#start: start address gene\n')
    fp.write('#end: end address of gene\n')
    fp.write('#gene: name of gene\n')
    fp.write('#actual_coverage_hp1: actual gene coverage value of HP-1\n')
    fp.write('#actual_coverage_hp2: actual gene coverage value of HP-2\n')
    fp.write('#adjusted_coverage_hp1: adjusted gene coverage value of HP-1\n')
    fp.write('#adjusted_coverage_hp2: adjusted gene coverage value of HP-2\n')
    fp.write('#copynumber_state_hp1: gene copy number state of HP-1\n')
    fp.write('#copynumber_state_hp2: gene copy number state of HP-2\n')
    fp.write('#chr\tstart\tend\tgene\tactual_coverage_hp1\tactual_coverage_hp2\tadjusted_coverage_hp1\tadjusted_coverage_hp2\thp1_state\thp2_state\n')

    df.to_csv(fp, sep='\t', index=False, header=False)


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

def write_segments_coverage_dict(coverage_segments, output, args):
    with open(args.out_dir_plots+'/coverage_data/' + output, 'a') as fp:
        for items in coverage_segments:
            if not items == None:
                fp.write("%s\n" % items)

def write_segments_coverage(coverage_segments, output, args):
    with open(args.out_dir_plots + '/bed_output/' + output, 'a') as fp:
        for items in coverage_segments:
            if not items == None:
                fp.write("%s\n" % items)

def write_header_comments(header, header_comments, output, args):
    with open(args.out_dir_plots + '/bed_output/' + output, 'a') as fp:
        fp.write(header_comments)
        fp.write(header)


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



def mask_df_states(haplotype_df, centers, integer_fractional_means):
    for i in range(len(integer_fractional_means)):
        haplotype_df['state'].mask(haplotype_df['state'] == centers[i], integer_fractional_means[i], inplace=True)

    return haplotype_df
def write_copynumber_segments_csv(df_hp1, df_hp2, haplotype_df_, args, centers, integer_fractional_means, hp, filename, p_value, is_half):

    haplotype_df = haplotype_df_.copy()

    #uniques = sorted(haplotype_df['state'].unique())
    #integer_fractional_means = sorted([i for i in range(0, len(uniques))])
    if not "subclonal" in filename:
        for i in range(len(integer_fractional_means)):
            haplotype_df['state'].mask(haplotype_df['state'] == centers[i], integer_fractional_means[i], inplace=True)
    #haplotype_df = mask_df_states(haplotype_df_copy, centers, integer_fractional_means)

    if 'p_value' in haplotype_df.columns:
        haplotype_df['if_subclonal'] = ['Y' if element < args.confidence_subclonal_score and not state == 0  else 'N' for (element,state) in zip(haplotype_df.p_value.values.tolist(), haplotype_df.state.values.tolist())]
        haplotype_df['state'] =  [subclonal_values_adjusted(element, centers) if sub == 'Y' else integers_values_adjusted(element, centers) for (element, sub) in zip(haplotype_df.state.values.tolist(), haplotype_df.if_subclonal.values.tolist())]
        haplotype_df = haplotype_df.rename(columns={'chromosome': 'chr', 'start': 'start', 'end': 'end', 'depth': 'coverage', 'state': 'copynumber_state', 'p_value': 'confidence', 'if_subclonal': 'if_subclonal'})
    else:
        if hp == 1:
            haplotype_df,_ = add_confidence_score_cn_segemnts(centers, haplotype_df, haplotype_df, df_hp1, df_hp2, args)
        else:
            _,haplotype_df = add_confidence_score_cn_segemnts(centers, haplotype_df, haplotype_df, df_hp1, df_hp2, args)

        haplotype_df = haplotype_df.rename(columns={'chromosome': 'chr', 'start': 'start', 'end': 'end', 'depth':'coverage', 'state':'copynumber_state', 'confidence_value': 'confidence'})

    haplotype_df['coverage'] = haplotype_df['coverage'].apply(lambda x: round(x, 2))
    if args.breakpoints:
        _, _, _, bps_ids_all, bps, bps_bnd = sv_vcf_bps_cn_check(args.breakpoints, args)
        bps_ids_global= []
        for index, row in haplotype_df.iterrows():
            bps_ids = []
            for bp in bps_ids_all:
                if bp[0] == row['chr'] and int(bp[1]) >= int(row['start']) and int(bp[1]) <= int(row['end']):
                    bps_ids.append(bp[2])
            bps_ids_global.append(list(set(bps_ids)))
        haplotype_df['svs_breakpoints_ids'] = bps_ids_global

    if args.without_phasing:
        fp = open(args.out_dir_plots + '/bed_output/' + args.genome_name + filename, 'a')
        fp.write('#chr: chromosome number\n')
        fp.write('#start: start address for CN segment\n')
        fp.write('#end: end address for CN segment\n')
        fp.write('#coverage: median coverage for this segment\n')
        fp.write('#copynumber_state: detected copy number state (integer/fraction)\n')

        if 'confidence' in haplotype_df.columns:
            fp.write('#confidence: confidence score\n')
            if "subclonal" in filename:
                fp.write('#if_subclonal: if entry is subclonal [Y/N]\n')
        if args.breakpoints:
            fp.write('#svs_breakpoints_ids: corresponding structural variations (breakpoints) IDs from VCF file\n')

        if 'confidence' in haplotype_df.columns:
            if args.breakpoints:
                if "subclonal" in filename:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tif_subclonal\tsvs_breakpoints_ids\n')
                else:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tsvs_breakpoints_ids\n')
            else:
                if "subclonal" in filename:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tif_subclonal\n')
                else:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\n')
        else:
            if args.breakpoints:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tsvs_breakpoints_ids\n')
            else:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\n')

        haplotype_df.to_csv(fp, sep='\t', index=False, mode='a', header=False)
    else:
        if is_half:
            if not os.path.isdir(args.out_dir_plots + '/wgd/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output'):
                os.makedirs(args.out_dir_plots + '/wgd/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output')
            fp = open(args.out_dir_plots +'/wgd/'+ str(args.tumor_ploidy) + '_'+ str(args.tumor_purity)  +'_'+ str(p_value) +'/bed_output/' + args.genome_name + filename, 'a')
        else:
            if not os.path.isdir(args.out_dir_plots + '/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output'):
                os.makedirs(args.out_dir_plots + '/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output')
            fp = open(args.out_dir_plots +'/'+ str(args.tumor_ploidy) + '_'+ str(args.tumor_purity)  +'_'+ str(p_value) +'/bed_output/' + args.genome_name + filename, 'a')
        fp.write('#chr: chromosome number\n')
        fp.write('#start: start address for CN segment\n')
        fp.write('#end: end address for CN segment\n')
        fp.write('#coverage: median coverage for this segment\n')
        fp.write('#copynumber_state: detected copy number state (integer/fraction)\n')

        if 'confidence' in haplotype_df.columns:
            fp.write('#confidence: confidence score\n')
            if "subclonal" in filename:
                fp.write('#if_subclonal: if entry is subclonal [Y/N]\n')

        if args.breakpoints:
            fp.write('#svs_breakpoints_ids: corresponding structural variations (breakpoints) IDs from VCF file\n')

        if 'confidence' in haplotype_df.columns:
            if args.breakpoints:
                if "subclonal" in filename:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tif_subclonal\tsvs_breakpoints_ids\n')
                else:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tsvs_breakpoints_ids\n')
            else:
                if "subclonal" in filename:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tif_subclonal\n')
                else:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\n')
        else:
            if args.breakpoints:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tsvs_breakpoints_ids\n')
            else:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\n')

        haplotype_df.to_csv(fp, sep='\t', index=False, mode='a', header=False)

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

def slice_lists(predicate, iterable):
  i, x, size = 0, 0, len(iterable)
  while i < size-1:
    if predicate(iterable[i], iterable[i+1]):
      yield iterable[x:i+1]
      x = i + 1
    i += 1
  yield iterable[x:size]
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


def remove_outliers_iqr(data, threshold=5):
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    lower_bound = q1 - threshold * iqr
    upper_bound = q3 + threshold * iqr
    outliers_removed = data[(data >= lower_bound) & (data <= upper_bound)]
    return outliers_removed

def add_intermediaries(numbers, max_difference):

    if len(numbers) < 2:
        return numbers

    result = [numbers[0]]
    for i in range(1, len(numbers)):
        current_diff = numbers[i] - result[-1]
        if current_diff > max_difference:
            num_intermediaries = current_diff // max_difference
            for j in range(1, num_intermediaries + 1):
                result.append(result[-1] + j + max_difference)
        result.append(numbers[i])

    return result

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


def parse_sv_vcf(path):
    my_parser = VCFParser(infile=path, split_variants=True, check_info = True)
    bp_junctions = [[]]
    for variant in my_parser:
        if ("INV" in variant['ID'] or "DUP" in variant['ID'] or "INS" in variant['ID'] or "DEL" in variant['ID']) and int(variant['info_dict']['SVLEN'][0]) > 50000:
            if int(variant['POS']) < int(variant['POS']) + int(variant['info_dict']['SVLEN'][0]):
                bp_junctions.append([variant['CHROM'], int(variant['POS']), variant['CHROM'], int(variant['POS']) + int(variant['info_dict']['SVLEN'][0])])
            else:
                bp_junctions.append([variant['CHROM'], int(variant['POS']) + int(variant['info_dict']['SVLEN'][0]), variant['CHROM'], int(variant['POS'])])
        elif "BND" in variant['ID'] and not "sBND" in variant['ID']:
            s = variant['ALT']
            for ch in ['[', ']', 'N']:
                if ch in s:
                    s = s.replace(ch, '')
            chr2_id = s.split(':')[0]
            chr2_end = int(s.split(':')[1])
            if variant['CHROM'] == chr2_id:
                if int(variant['POS']) < chr2_end:
                    bp_junctions.append([variant['CHROM'], int(variant['POS']), chr2_id, chr2_end])
                else:
                    bp_junctions.append([variant['CHROM'], chr2_end, chr2_id, int(variant['POS'])])
            else:
                bp_junctions.append([variant['CHROM'], int(variant['POS']), chr2_id, chr2_end])

    df = pd.DataFrame(bp_junctions[1:], columns = ['chr', 'pos', 'chr1', 'pos1'])
    return df.drop_duplicates()


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

def weighted_means(vals, weights):
    if not len(vals) == len(weights) or sum(weights) == 0:
        return 0
    weighted_vals = []
    vals_n_weights = [(vals[i], weights[i]) for i in range(0, len(weights))]
    for tup in vals_n_weights:
        weighted_vals.append(round(tup[0] * tup[1], 3))
    return sum(weighted_vals)/sum(weights)


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

def check_adjust_last_cn_states(cens, df_segs_hp1, df_segs_hp2):
    df_hp1 = df_segs_hp1.copy()
    df_hp2 = df_segs_hp2.copy()
    df_hp1.loc[df_hp1['state'] < cens[0], 'state'] = cens[0]
    df_hp2.loc[df_hp2['state'] < cens[0], 'state'] = cens[0]

    return df_hp1, df_hp2

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


def overlap_check(start, end, starts, ends):
  for i in range(len(starts)):
    #if (start < ends[i] and end > starts[i]) or (start <= starts[i] and end >= ends[i]):
    if start >= starts[i] and start <= ends[i]:
      return True, i
  return False, -1


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
def genes_phase_correction(df_genes, df_segs_hp1, df_segs_hp2, args, centers, integer_fractional_centers):
    chroms = get_contigs_list(args.contigs)
    df_genes_updated = []

    for index, chrom in enumerate(chroms):
        df_gene = df_genes[df_genes['chr'] == chrom]
        df_seg_hp1 = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
        df_seg_hp2 = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]

        seg_coverage_hp1 = df_seg_hp1.depth.values.tolist()
        seg_coverage_hp2 = df_seg_hp2.depth.values.tolist()

        seg_state_hp1 = df_seg_hp1.state.values.tolist()
        seg_state_hp2 = df_seg_hp2.state.values.tolist()

        genes_starts = df_gene.start.values.tolist()
        genes_ends = df_gene.end.values.tolist()
        genes_coverage_hp1 = df_gene.hp1.values.tolist()
        genes_coverage_hp2 = df_gene.hp2.values.tolist()

        gene_hp1_state = [0 for x in range(len(genes_coverage_hp1))]
        gene_hp2_state = [0 for x in range(len(genes_coverage_hp2))]

        for i, (start_b, end_b) in enumerate(zip(genes_starts, genes_ends)):
            hp_1_val = 0
            hp_2_val = 0
            check, index = overlap_check(start_b, end_b, df_seg_hp1.start.values.tolist(), df_seg_hp1.end.values.tolist())
            if check:
                hp_1_val = centers[integers_values_adjusted(seg_state_hp1[index], centers)]
                #gene_hp1_state[i] = integer_fractional_centers[centers.index(min(centers, key=lambda x:abs(x - seg_coverage_hp1[index])))]  #(integer_fractional_centers[centers.index(seg_coverage_hp1[index])])
            check, index = overlap_check(start_b, end_b, df_seg_hp2.start.values.tolist(), df_seg_hp2.end.values.tolist())
            if check:
                hp_2_val = centers[integers_values_adjusted(seg_state_hp2[index], centers)]
                #gene_hp2_state[i] = integer_fractional_centers[centers.index(min(centers, key=lambda x:abs(x - seg_coverage_hp2[index])))]

            # if hp_1_val == 0 and hp_2_val == 0:
            #     gene_hp1_state[i] = integers_values_adjusted(genes_coverage_hp1[i], centers)
            #     gene_hp2_state[i] = integers_values_adjusted(genes_coverage_hp2[i], centers)
            #     continue

            #print(genes_coverage_hp1[i], genes_coverage_hp2[i])
            #print(abs(genes_coverage_hp1[i] -  hp_2_val) + abs(genes_coverage_hp2[i] -  hp_1_val))
            #print(abs(genes_coverage_hp1[i] -  hp_1_val) + abs(genes_coverage_hp2[i] -  hp_2_val))

            # if abs(genes_coverage_hp1[i] -  hp_2_val) < abs(genes_coverage_hp1[i] -  hp_1_val):
            #     new_hp2 = genes_coverage_hp1[i]
            #     new_hp1 = genes_coverage_hp2[i]
            #     genes_coverage_hp2[i] = new_hp2
            #     genes_coverage_hp1[i] = new_hp1

            genes_coverage_hp1[i] = hp_1_val #centers[integers_values_adjusted(genes_coverage_hp1[i], centers)]
            genes_coverage_hp2[i] = hp_2_val #centers[integers_values_adjusted(genes_coverage_hp2[i], centers)]

            gene_hp1_state[i] = integers_values_adjusted(genes_coverage_hp1[i], centers)
            gene_hp2_state[i] = integers_values_adjusted(genes_coverage_hp2[i], centers)

        df_genes_updated.append(pd.DataFrame(list(zip(df_gene.chr.values.tolist(), df_gene.start.values.tolist(), df_gene.end.values.tolist(), \
                     df_gene.gene.values.tolist(), [round(l,2) for l in df_gene.hp1.values.tolist()], [round(l,2) for l in df_gene.hp2.values.tolist()], [round(l,2) for l in genes_coverage_hp1], [round(l,2) for l in genes_coverage_hp2], gene_hp1_state, gene_hp2_state)),
            columns=['chr', 'start', 'end', 'gene', 'actual_coverage_hp1', 'actual_coverage_hp2', 'adjusted_coverage_hp1', 'adjusted_coverage_hp2', 'copynumber_hp1_state', 'copynumber_hp2_state']))

    return pd.concat(df_genes_updated)

def update_genes_phase_corrected_coverage(args, df_segs_hp1, df_segs_hp2, p_value, centers, integer_fractional_centers, is_half):

    df_genes = csv_df_chromosomes_sorter(args.out_dir_plots + '/coverage_data/cancer_genes_coverage.csv', ['chr','start','end','gene', 'hp1', 'hp2'])
    #write_df_csv(df_genes, args.out_dir_plots + '/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output/' + 'cancer_genes_coverage.csv')
    df_genes = genes_phase_correction(df_genes, df_segs_hp1, df_segs_hp2, args, centers, integer_fractional_centers)
    if is_half:
        write_df_csv(df_genes, args.out_dir_plots + '/wgd/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output/' + 'genes_copynumber_states.bed')
    else:
        write_df_csv(df_genes, args.out_dir_plots + '/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output/' + 'genes_copynumber_states.bed')

    return df_genes

def extract_breakpoints_additional(args):
    if args.breakpoints:
        df_var_bins, df_var_bins_1 = get_chromosomes_bins(args.target_bam[0], args.bin_size, args)
        return df_var_bins_1
        # chroms = get_contigs_list(args.contigs)
        # for index, chrom in enumerate(chroms):
        #     df_var_bins_chr_1 = df_var_bins_1[df_var_bins_1['chr'] == chrom]
        #     ref_start_values_1 = df_var_bins_chr_1.start.values.tolist()
    else:
        return None

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

def move_folder(source, destination):
    """Moves a folder from the source path to the destination path.

    Args:
        source: The path of the folder to be moved.
        destination: The path where the folder should be moved.
    """
    try:
        shutil.move(source, destination)
        logger.info(f"Folder '{source}' moved successfully to '{destination}'")
    except FileNotFoundError:
        logger.info(f"Error: Folder '{source}' not found.")
    except FileExistsError:
         logger.info(f"Error: A file or directory with the name '{os.path.basename(source)}' already exists in '{destination}'.")
    except Exception as e:
        logger.info(f"An error occurred: {e}")

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


def centromere_regions_blacklist_bins(args, df_hp1_, df_hp2_, df_segs_hp1_updated_, df_segs_hp2_updated_):
    df_segs_hp1 = df_segs_hp1_updated_.copy()
    df_segs_hp2 = df_segs_hp2_updated_.copy()

    updated_hp1 = []
    updated_hp2 = []

    fileDir = os.path.dirname(__file__)  # os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)

    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        df_centm_chrom = df_centm[df_centm['chr'] == chrom]
        if not df_centm_chrom.empty:
            cents = [df_centm_chrom.start.values.tolist()[0], df_centm_chrom.end.values.tolist()[0]]
        else:
            cents = [0, 0]
        df_segs_hp_1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
        df_segs_hp_2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]

        df_hp_1_chrom = df_hp1_[df_hp1_['chr'] == chrom]
        df_hp_2_chrom = df_hp2_[df_hp2_['chr'] == chrom]
        df_hp_1_chrom = df_hp_1_chrom.reset_index(drop=True)
        df_hp_2_chrom = df_hp_2_chrom.reset_index(drop=True)

        ref_start_values = df_hp_1_chrom.start.values.tolist()

        for idx, s1 in df_segs_hp_1_chrom.iterrows():
            if (df_segs_hp_1_chrom.loc[idx, 'start'] == cents[0] or df_segs_hp_1_chrom.loc[idx, 'start'] == cents[0]+1) and df_segs_hp_1_chrom.loc[idx, 'end'] == cents[1]:# or df_segs_hp_1_chrom.loc[idx+1, 'start'] == cents[0] and df_segs_hp_1_chrom.loc[idx+1, 'end'] == cents[1]:
                internal_bins = [k for k in ref_start_values if k >= s1['start'] and k <= s1['end']]
                if internal_bins:
                    j = ref_start_values.index(internal_bins[0])
                    for l in range(len(internal_bins)):
                        df_hp_1_chrom.loc[j, 'hp1'] = 3300
                        j = j+1

        for idx, s2 in df_segs_hp_2_chrom.iterrows():
            if (df_segs_hp_2_chrom.loc[idx, 'start'] == cents[0] or df_segs_hp_2_chrom.loc[idx, 'start'] == cents[0]+1) and df_segs_hp_2_chrom.loc[idx, 'end'] == cents[1]:# or df_segs_hp_1_chrom.loc[idx+1, 'start'] == cents[0] and df_segs_hp_1_chrom.loc[idx+1, 'end'] == cents[1]:
                internal_bins = [k for k in ref_start_values if k >= s2['start'] and k <= s2['end']]
                if internal_bins:
                    j = ref_start_values.index(internal_bins[0])
                    for l in range(len(internal_bins)):
                        df_hp_2_chrom.loc[j, 'hp2'] = 3300
                        j = j+1

        updated_hp1.append(df_hp_1_chrom)
        updated_hp2.append(df_hp_2_chrom)

    segs_hp1 = pd.concat(updated_hp1)
    segs_hp2 = pd.concat(updated_hp2)

    return segs_hp1, segs_hp2


def centromere_regions_blacklist(args, df_segs_hp1_, df_segs_hp2_):
    df_segs_hp1 = df_segs_hp1_.copy()
    df_segs_hp2 = df_segs_hp2_.copy()

    updated_hp1_segs = []
    updated_hp2_segs = []

    fileDir = os.path.dirname(__file__)  # os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)

    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        df_centm_chrom = df_centm[df_centm['chr'] == chrom]
        if not df_centm_chrom.empty:
            cents = [df_centm_chrom.start.values.tolist()[0], df_centm_chrom.end.values.tolist()[0]]
        else:
            cents = [0, 0]
        df_segs_hp_1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
        df_segs_hp_2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]

        for idx, s1 in df_segs_hp_1_chrom.iterrows():
            if (df_segs_hp_1_chrom.loc[idx, 'start'] == cents[0] or df_segs_hp_1_chrom.loc[idx, 'start'] == cents[0]+1) and df_segs_hp_1_chrom.loc[idx, 'end'] == cents[1]:# or df_segs_hp_1_chrom.loc[idx+1, 'start'] == cents[0] and df_segs_hp_1_chrom.loc[idx+1, 'end'] == cents[1]:
                df_segs_hp_1_chrom.loc[idx, 'depth'] = 3300
                df_segs_hp_1_chrom.loc[idx, 'state'] = 3300

        for idx, s2 in df_segs_hp_2_chrom.iterrows():
            if (df_segs_hp_2_chrom.loc[idx, 'start'] == cents[0] or df_segs_hp_2_chrom.loc[idx, 'start'] == cents[0]+1) and df_segs_hp_2_chrom.loc[idx, 'end'] == cents[1]:# or df_segs_hp_1_chrom.loc[idx+1, 'start'] == cents[0] and df_segs_hp_1_chrom.loc[idx+1, 'end'] == cents[1]:
                df_segs_hp_2_chrom.loc[idx, 'depth'] = 3300
                df_segs_hp_2_chrom.loc[idx, 'state'] = 3300

        updated_hp1_segs.append(df_segs_hp_1_chrom)
        updated_hp2_segs.append(df_segs_hp_2_chrom)

    segs_hp1 = pd.concat(updated_hp1_segs)
    segs_hp2 = pd.concat(updated_hp2_segs)

    return segs_hp1, segs_hp2

def extract_centromere_regions(args):
    fileDir = os.path.dirname(__file__)
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)

    return df_centm

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
