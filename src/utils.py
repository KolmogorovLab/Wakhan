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
from src.cna.copynumber import subclonal_values_adjusted, integers_values_adjusted, add_confidence_score_cn_segemnts
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
