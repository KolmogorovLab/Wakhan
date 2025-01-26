import csv
import numpy as np
import pandas as pd
import pysam
import scipy.stats as stats
import os
import math
import statistics
import logging
import ruptures as rpt

from smoothing import smoothing
#from hmm import call_copynumbers
from extras import get_contigs_list, sv_vcf_bps_cn_check

def get_chromosomes_bins_replica(bam_file, bin_size, args):
    bed=[]
    bam_alignment = pysam.AlignmentFile(bam_file)
    headers = bam_alignment.header
    seq_dict = headers['SQ']
    region = [''] * len(seq_dict)
    chrs = [''] * len(seq_dict)
    head, tail = os.path.split(bam_file)
    chroms = get_contigs_list(args.contigs)
    for i, seq_elem in enumerate(seq_dict):
        region[i] = seq_elem['LN']
        chrs[i] = seq_elem['SN']
        start=0
        end=bin_size
        if chrs[i] in chroms:
            for c in range(0,region[i],bin_size):
                if end > region[i]:
                    bed.append(chrs[i]+'\t'+str(start)+'\t'+str(region[i]))
                else:
                    bed.append(chrs[i]+'\t'+str(start)+'\t'+str(end))
                start=end+1
                end+=bin_size
    return bed

def get_chromosomes_regions(args):
    chroms = get_contigs_list(args.contigs)
    bam_alignment = pysam.AlignmentFile(args.target_bam[0])
    headers = bam_alignment.header
    seq_dict = headers['SQ']
    region = [''] * len(chroms)
    for i, seq_elem in enumerate(seq_dict):
        if seq_elem['SN'] in chroms:
            region[chroms.index(seq_elem['SN'])] = seq_elem['LN']

    return region

def get_chromosomes_bins_bam(bam_file, bin_size, args):
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
                    bed.append([tail, chrs[i], start, region[i]])
                else:
                    bed.append([tail, chrs[i], start, end])
                start=end+1
                end+=bin_size
    return bed

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

def chunk_range(start, end, num_chunks):

  if num_chunks == 0:
    raise ValueError("Number of chunks cannot be zero")

  # Calculate the ideal chunk size (might have a remainder)
  ideal_chunk_size = (end - start) // num_chunks

  # Initialize list to store chunks
  chunks = []
  current_start = start
  last = end

  # Iterate and create chunks
  for _ in range(num_chunks):
    # Handle last chunk potentially being larger
    if _ == num_chunks :
      end_value = end
    else:
      end_value = current_start + ideal_chunk_size
    if not last == end_value:
        chunks.append([current_start, end_value-1])
    else:
        chunks.append([current_start, end_value])
    current_start = end_value

  return chunks


def update_bins_with_bps(bed, bps, bps_bnd, args, region):

    df_bps_segs = []
    chroms = get_contigs_list(args.contigs)
    for i, val in enumerate(chroms):
        bp_junctions = []

        #########################################################################
        bps_bnds_segs = []
        for j, bp in enumerate(bps):
            if bp[0] == val:
                bp_junctions.append([bp[1], bp[2]])

        for j, bps_bnds in enumerate(bps_bnd):
            if bps_bnds[0] == val:
                bps_bnds_segs.append([bps_bnds[1]])
        #########################################################################
        bp_junctions_new = []
        indices_ol = []
        for k in range(len(bp_junctions)-1):
            a = bp_junctions[k][0]
            b = bp_junctions[k][1]
            c = bp_junctions[k+1][0]
            d = bp_junctions[k+1][1]

            if c < b:
                indices_ol.append(k)
                bp_junctions_new.append([a, c-1])

        columns = ['start', 'end']
        df_bps = pd.DataFrame(bp_junctions, columns=columns)

        if bp_junctions_new:
            for index in sorted(list(set(indices_ol)), reverse=True):
                del bp_junctions[index]

            for m, inner in enumerate(bp_junctions_new):
                bp_junctions.append(inner)
            bp_junctions = sorted(bp_junctions)
        #########################################################################
        bnd_bp_list = []
        indices_ol = []
        for l1, bps_bnds_seg in enumerate(bps_bnds_segs):
            for l2, (st,en) in enumerate(bp_junctions):
                if bps_bnds_seg[0] >= st and bps_bnds_seg[0] <= en:
                    if l2 in indices_ol:
                        break
                    indices_ol.append(l2)
                    # if (bps_bnds_seg[0]-1) - st > 100:
                    #     lst1 = chunk_range(st,bps_bnds_seg[0]-1, 100)
                    #     bnd_bp_list.extend(lst1)
                    # else:
                    #     bnd_bp_list.append([st,bps_bnds_seg[0]-1])
                    #
                    # if en - bps_bnds_seg[0] > 100:
                    #     lst2 = chunk_range(bps_bnds_seg[0], en, 100)
                    #     bnd_bp_list.extend(lst2)
                    # else:
                    #     bnd_bp_list.append([bps_bnds_seg[0], en])
                    # break
                    bnd_bp_list.append([st, bps_bnds_seg[0] - 1])
                    bnd_bp_list.append([bps_bnds_seg[0], en])
                    break

        if bnd_bp_list:
            for index in sorted(list(set(indices_ol)), reverse=True):
                del bp_junctions[index]

            for m, inner in enumerate(bnd_bp_list):
                bp_junctions.append(inner)
            bp_junctions = sorted(bp_junctions)

        bp_junctions.append([region[i] - 1, region[i]])
        #########################################################################
        missing_segs = []
        last_val = 0
        for m, inner in enumerate(bp_junctions):
            if inner[0] > last_val:
                missing_segs.append([last_val, inner[0]-1])
            last_val = inner[1] +1

        for m, inner in enumerate(missing_segs):
            bp_junctions.append(inner)
        bp_junctions = sorted(bp_junctions)
        #########################################################################
        bin_size_segs = []
        indices_ol = []
        for l, inner in enumerate(bp_junctions):
            if inner[1] - inner[0] > args.bin_size:
                indices_ol.append(l)
                for i in range(inner[0], inner[1], args.bin_size):
                    # Handle the last bin potentially being smaller
                    next_end = min(i + args.bin_size, inner[1])
                    if not next_end == inner[1]:
                        next_end = next_end-1
                    bin_size_segs.append([i, next_end])

        if bin_size_segs:
            for index in sorted(list(set(indices_ol)), reverse=True):
                del bp_junctions[index]

            for m, inner in enumerate(bin_size_segs):
                bp_junctions.append(inner)
            bp_junctions = sorted(bp_junctions)
        #########################################################################
        indices_ = []
        for l, inner in enumerate(bp_junctions):
            if inner[1] < inner[0]:
                indices_.append(l)
        if indices_:
            for index in sorted(list(set(indices_)), reverse=True):
                del bp_junctions[index]
        #########################################################################
        bp_junct = []
        for l, inner in enumerate(bp_junctions):
            bp_junct.append([val, inner[0], inner[1]])

        columns = ['chr', 'start', 'end']
        df_bps_new = pd.DataFrame(bp_junct, columns=columns)
        df_bps_segs.append(df_bps_new)
        #########################################################################

    df  = pd.concat(df_bps_segs)
    df.drop_duplicates(subset=['chr', 'start', 'end'], keep=False, inplace=True)
    return df

def update_bins_with_bps_new(bed, bps, bps_bnd, args, region):

    df_bps_segs = []
    df_bps_segs_1 = []
    chroms = get_contigs_list(args.contigs)
    for i, val in enumerate(chroms):
        #if val == 'chr20':
        #    print('here')
        bp_junctions = []

        #########################################################################
        bps_bnds_segs = []
        for j, bp in enumerate(bps):
            if bp[0] == val:
                bp_junctions.append([bp[1], bp[2]])

        for j, bps_bnds in enumerate(bps_bnd):
            if bps_bnds[0] == val:
                bps_bnds_segs.append([bps_bnds[1]])
        #########################################################################
        bp_junctions_new = []
        indices_ol = []
        for k in range(len(bp_junctions)-1):
            a = bp_junctions[k][0]
            b = bp_junctions[k][1]
            c = bp_junctions[k+1][0]
            d = bp_junctions[k+1][1]

            if c < b:
                indices_ol.append(k)
                bp_junctions_new.append([a, c-1])

        columns = ['start', 'end']
        df_bps = pd.DataFrame(bp_junctions, columns=columns)

        if bp_junctions_new:
            for index in sorted(list(set(indices_ol)), reverse=True):
                del bp_junctions[index]

            for m, inner in enumerate(bp_junctions_new):
                bp_junctions.append(inner)
            bp_junctions = sorted(bp_junctions)
        #########################################################################
        bps_bnds_segs_ = bps_bnds_segs
        bps_bnds_segs_ = [x for xs in bps_bnds_segs_ for x in xs]
        bps_bnds_segs_ = list(np.unique(bps_bnds_segs_))
        bps_bnds_segs_new = [[]]
        for m, reg in enumerate(bps_bnds_segs_):
            bps_bnds_segs_new.append([reg, reg+1])
        bps_bnds_segs_new = bps_bnds_segs_new[1:]
        #########################################################################
        bnd_bp_list = []
        indices_ol = []
        for l1, bps_bnds_seg in enumerate(bps_bnds_segs):
            for l2, (st,en) in enumerate(bp_junctions):
                if bps_bnds_seg[0] >= st and bps_bnds_seg[0] <= en:
                    if l2 in indices_ol:
                        break
                    indices_ol.append(l2)
                    # if (bps_bnds_seg[0]-1) - st > 100:
                    #     lst1 = chunk_range(st,bps_bnds_seg[0]-1, 100)
                    #     bnd_bp_list.extend(lst1)
                    # else:
                    #     bnd_bp_list.append([st,bps_bnds_seg[0]-1])
                    #
                    # if en - bps_bnds_seg[0] > 100:
                    #     lst2 = chunk_range(bps_bnds_seg[0], en, 100)
                    #     bnd_bp_list.extend(lst2)
                    # else:
                    #     bnd_bp_list.append([bps_bnds_seg[0], en])
                    # break
                    bnd_bp_list.append([st, bps_bnds_seg[0] - 1])
                    bnd_bp_list.append([bps_bnds_seg[0], en])
                    break

        if bnd_bp_list:
            for index in sorted(list(set(indices_ol)), reverse=True):
                del bp_junctions[index]

            for m, inner in enumerate(bnd_bp_list):
                bp_junctions.append(inner)
            bp_junctions = sorted(bp_junctions)

        bp_junctions.append([region[i] - 1, region[i]])
        #########################################################################
        bp_junct_1 = []
        for l, inner in enumerate(bp_junctions):
            bp_junct_1.append([val, inner[0], inner[1]])
        #########################################################################
        for l1, inner in enumerate(bps_bnds_segs_new): #ad-hoc solution to issue in missing ol segemnts
            bp_junct_1.append([val, inner[0], inner[1]])
        bp_junct_1 = sorted(bp_junct_1)
        #########################################################################
        columns = ['chr', 'start', 'end']
        df_bps_new_1 = pd.DataFrame(bp_junct_1, columns=columns)
        df_bps_segs_1.append(df_bps_new_1)
        #########################################################################
        #########################################################################
        missing_segs = []
        last_val = 0
        for m, inner in enumerate(bp_junctions):
            if inner[0] > last_val:
                missing_segs.append([last_val, inner[0]-1])
            last_val = inner[1] +1
        #########################################################################
        # for m, inner in enumerate(missing_segs):
        #     bp_junctions.append(inner)
        # bp_junctions = sorted(bp_junctions)
        #########################################################################
        bin_size_segs = []
        indices_ol = []
        for l, inner in enumerate(missing_segs):
            if inner[1] - inner[0] > args.bin_size:
                indices_ol.append(l)
                for i in range(inner[0], inner[1], args.bin_size):
                    # Handle the last bin potentially being smaller
                    next_end = min(i + args.bin_size, inner[1])
                    if not next_end == inner[1]:
                        next_end = next_end-1
                    bin_size_segs.append([i, next_end])

        if bin_size_segs:
            for index in sorted(list(set(indices_ol)), reverse=True):
                del missing_segs[index]

            for m, inner in enumerate(bin_size_segs):
                bp_junctions.append(inner)
            bp_junctions = sorted(bp_junctions)
        #########################################################################
        indices_ = []
        for l, inner in enumerate(bp_junctions):
            if inner[1] < inner[0]:
                indices_.append(l)
        if indices_:
            for index in sorted(list(set(indices_)), reverse=True):
                del bp_junctions[index]
        #########################################################################
        bp_junct = []
        for l, inner in enumerate(bp_junctions):
            bp_junct.append([val, inner[0], inner[1]])

        columns = ['chr', 'start', 'end']
        df_bps_new = pd.DataFrame(bp_junct, columns=columns)
        df_bps_segs.append(df_bps_new)
        #########################################################################

    df  = pd.concat(df_bps_segs)
    df_1 = pd.concat(df_bps_segs_1)
    df.drop_duplicates(subset=['chr', 'start', 'end'], keep=False, inplace=True)
    df_1.drop_duplicates(subset=['chr', 'start', 'end'], keep=False, inplace=True)
    return df, df_1

def chromosomes_sorter(label):
    from itertools import takewhile
    # Strip "chr" prefix
    chrom = (label[3:] if label.lower().startswith('chr')
             else label)
    if chrom in ('X', 'Y'):
        key = (1000, chrom)
    else:
        # Separate numeric and special chromosomes
        nums = ''.join(takewhile(str.isdigit, chrom))
        chars = chrom[len(nums):]
        nums = int(nums) if nums else 0
        if not chars:
            key = (nums, '')
        elif len(chars) == 1:
            key = (2000 + nums, chars)
        else:
            key = (3000 + nums, chars)
    return key

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

def csv_df_chromosomes_sorter(path, names, sept='\t'):
    dataframe = pd.read_csv(path, sep=sept, names=names)
    dataframe['chr'] = dataframe['chr'].astype(str)
    #if not dataframe['chr'].iloc[0].startswith('chr'):
    #    dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', names[1]], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

def df_chromosomes_sorter(dataframe, names):
    dataframe['chr'] = dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', names[1]], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

def get_breakpoints(chrom, bp_file_path): #TODO add call in plots
    break_points = []
    with open(bp_file_path) as bp_file:
        next(bp_file)
        for st in bp_file:
            # st = '-chr1:2671683|+chr1:2673127,0,0,0,9,10,22'
            st = st.split(",")
            #Parse the file
            chr1 = (st[0].split("|"))[0].split(":")[0][1:]
            chr2 = (st[0].split("|"))[1].split(":")[0][1:]
            bp_pos1 = int((st[0].split("|"))[0].split(":")[1])
            bp_pos2 = int((st[0].split("|"))[1].split(":")[1])
            val_1=int(st[1])
            val_2=int(st[4])
            if ( val_1 == 0 and val_2 >= 3) and (abs(bp_pos1 - bp_pos2) > 1000) and (chr1 == chr2) and (chr2 == chrom):
                #1k condition and both connections should be on same chromosome for the moment
                mid=bp_pos1 + round((abs(bp_pos1 - bp_pos2) / 2))
                break_points.extend([bp_pos1, mid, bp_pos2])
    return break_points

def write_segments_coverage_dict(coverage_segments, output, args):
    with open(args.out_dir_plots+'/data/' + output, 'a') as fp:
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

def seperate_dfs_coverage(args, df, haplotype_1_values_updated, haplotype_2_values_updated, unphased):
    if args.without_phasing:
        return df[['chr', 'start', 'end', 'coverage']].copy()
    else:
        df_hp1 = df[['chr', 'start','end', 'hp1']].copy()
        df_hp2 = df[['chr', 'start','end', 'hp2']].copy()
        df_unphased = df[['chr', 'start','end', 'hp3']].copy()
        df_hp1['hp1'] = haplotype_1_values_updated
        df_hp2['hp2'] = haplotype_2_values_updated
        df_unphased['hp3'] = unphased
        return df_hp1, df_hp2, df_unphased

def flatten(values):
    return [item for sublist in values for item in sublist]

def flatten_smooth(hp1, hp2, unphased):
    hp1 = flatten(hp1)
    hp2 = flatten(hp2)
    unphased = flatten(unphased)
    unphased, hp1, hp2 = smoothing(unphased, hp1, hp2, conv_window_size=15)

    return hp1, hp2, unphased

def get_snps_frquncies_coverage_from_bam(df, chrom):
    df = df[df['chr'] == chrom]
    #df = dict(tuple(df.groupby('hp')))
    haplotype_1_position = df.pos.values.tolist()
    haplotype_1_coverage = df.freq_value_b.values.tolist()
    haplotype_2_position = df.pos.values.tolist()
    haplotype_2_coverage = df.freq_value_a.values.tolist()

    return haplotype_1_position, haplotype_1_coverage, haplotype_2_position, haplotype_2_coverage

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


def is_phasesets_check_simple_heuristics(ref_start_values_phasesets, ref_end_values_phasesets):
    ps_region_starts = []
    ps_region_ends = []
    for i, (start, end) in enumerate(zip(ref_start_values_phasesets, ref_end_values_phasesets)):
        if end - start > 2000000:
            ps_region_starts.append(start)
            ps_region_ends.append(end)

    if len(ps_region_starts) < 5:
        return True
    else:
        return False
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
                cent_chr_index = cent_chr.index(row['chromosome'])
                if not haplotype_df.loc[index, 'start'] - 1 == cent_starts[cent_chr_index] and not haplotype_df.loc[index, 'end'] == cent_ends[cent_chr_index]:
                    haplotype_df.loc[index, 'end'] = min(starts, key=lambda x: abs(x - int(row['start']))) - 1
        #haplotype_df['svs_breakpoints_ids'] = bps_ids_global

        for index, row in haplotype_df.iterrows():
            cent_chr_index = cent_chr.index(row['chromosome'])
            if not haplotype_df.loc[index, 'start'] - 1 == cent_starts[cent_chr_index] and not haplotype_df.loc[index, 'end'] == cent_ends[cent_chr_index]:
                if index < len(haplotype_df) and haplotype_df.loc[index, 'start'] > 0:
                    haplotype_df.loc[index, 'start'] = haplotype_df.loc[index - 1, 'end'] + 1

    return haplotype_df

def mask_df_states(haplotype_df, centers, integer_fractional_means):
    for i in range(len(integer_fractional_means)):
        haplotype_df['state'].mask(haplotype_df['state'] == centers[i], integer_fractional_means[i], inplace=True)

    return haplotype_df
def write_copynumber_segments_csv(haplotype_df_, args, centers, integer_fractional_means, hp, filename, p_value):

    haplotype_df = haplotype_df_.copy()
    #uniques = sorted(haplotype_df['state'].unique())
    #integer_fractional_means = sorted([i for i in range(0, len(uniques))])
    for i in range(len(integer_fractional_means)):
       haplotype_df['state'].mask(haplotype_df['state'] == centers[i], integer_fractional_means[i], inplace=True)
    #haplotype_df = mask_df_states(haplotype_df_copy, centers, integer_fractional_means)

    if 'p_value' in haplotype_df.columns:
        haplotype_df['if_subclonal'] = ['Y' if element < args.confidence_subclonal_score and  not state == 0  else 'N' for (element,state) in zip(haplotype_df.p_value.values.tolist(), haplotype_df.state.values.tolist())]
        haplotype_df['state'] =  [subclonal_values_adjusted(element, centers) if sub == 'Y' else integers_values_adjusted(element, centers) for (element, sub) in zip(haplotype_df.state.values.tolist(), haplotype_df.if_subclonal.values.tolist())]
        haplotype_df = haplotype_df.rename(columns={'chromosome': 'chr', 'start': 'start', 'end': 'end', 'depth': 'coverage', 'state': 'copynumber_state', 'p_value': 'confidence', 'if_subclonal': 'if_subclonal'})
    else:
        haplotype_df = haplotype_df.rename(columns={'chromosome': 'chr', 'start': 'start', 'end': 'end', 'depth':'coverage', 'state':'copynumber_state'})

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
            fp.write('#if_subclonal: if entry is subclonal [Y/N]\n')
        if args.breakpoints:
            fp.write('#svs_breakpoints_ids: corresponding structural variations (breakpoints) IDs from VCF file\n')

        if 'confidence' in haplotype_df.columns:
            if args.breakpoints:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tif_subclonal\tsvs_breakpoints_ids\n')
            else:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tif_subclonal\n')
        else:
            if args.breakpoints:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tsvs_breakpoints_ids\n')
            else:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\n')

        haplotype_df.to_csv(fp, sep='\t', index=False, mode='a', header=False)
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
            fp.write('#if_subclonal: if entry is subclonal [Y/N]\n')

        if args.breakpoints:
            fp.write('#svs_breakpoints_ids: corresponding structural variations (breakpoints) IDs from VCF file\n')

        if 'confidence' in haplotype_df.columns:
            if args.breakpoints:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tif_subclonal\tsvs_breakpoints_ids\n')
            else:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tif_subclonal\n')
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


def change_point_detection_means(args, df_chrom, ref_start_values, ref_start_values_1, df_centm_chrom):
    df_means_chr = []
    if args.without_phasing:
        means = df_chrom.coverage.values.tolist()
        snps_mean, snps_len, snps_pos = change_point_detection_algo(args.bin_size, means, ref_start_values, args, ref_start_values_1, df_centm_chrom)
        snps_pos_start = []
        snps_pos_end = []
        for i in range(len(snps_pos)-1):
            snps_pos_start.append(snps_pos[i] if snps_pos[i] < 1 else snps_pos[i]+1)
            snps_pos_end.append(snps_pos[i+1])
        chr_list = [df_chrom['chr'].iloc[0] for ch in range(len(snps_mean))]
        df_means_chr = pd.DataFrame(list(zip(chr_list, snps_pos_start, snps_pos_end, snps_mean)), columns=['chromosome', 'start', 'end', 'state'])

        return snps_mean, snps_len, df_means_chr
    else:
        df_means_chr = []
        haplotype1_means = df_chrom.hp1.values.tolist()
        haplotype2_means = df_chrom.hp2.values.tolist()
        snps_haplotype1_mean, snps_haplotype1_len, snps_haplotype1_pos = change_point_detection_algo(args.bin_size, df_chrom.hp1.values.tolist(), ref_start_values, args, ref_start_values_1, df_centm_chrom)
        snps_haplotype2_mean, snps_haplotype2_len, snps_haplotype2_pos = change_point_detection_algo(args.bin_size, df_chrom.hp2.values.tolist(), ref_start_values, args, ref_start_values_1, df_centm_chrom)

        cent_index_hp1 = snps_haplotype1_pos.index([df_centm_chrom.start.values.tolist()[0], df_centm_chrom.end.values.tolist()[0]][0]) -1
        cent_index_hp2 = snps_haplotype2_pos.index([df_centm_chrom.start.values.tolist()[0], df_centm_chrom.end.values.tolist()[0]][0]) -1
        #if not cent_index_hp1 == 0:

        if args.breakpoints:
            cent_index_hp1 = cent_index_hp1 + 1
            cent_index_hp2 = cent_index_hp2 + 1
            # snps_haplotype1_mean = []
            # snps_haplotype2_mean = []
            snps_haplotype1_len = []
            snps_haplotype2_len = []
            snps_pos = sorted(list(set(snps_haplotype1_pos + snps_haplotype2_pos)))
            snps_pos_start = list(map(lambda x: x + 1 if x > 0 else x, snps_pos[:-1]))
            snps_pos_end = list(map(lambda x: x if x > 0 else x, snps_pos[1:]))
            for index, (start, end) in enumerate(zip(snps_pos_start, snps_pos_end)):
                if cent_index_hp1 == index:
                    snps_haplotype1_mean[index] = 0
                    snps_haplotype2_mean[index] = 0
                else:
                    snps_haplotype1_mean[index] = statistics.median(remove_outliers_iqr(np.array(haplotype1_means[start//args.bin_size:end//args.bin_size])))
                    snps_haplotype2_mean[index] = statistics.median(remove_outliers_iqr(np.array(haplotype2_means[start//args.bin_size:end//args.bin_size])))
                snps_haplotype1_len.append(len(haplotype1_means[start//args.bin_size:end//args.bin_size]))
                snps_haplotype2_len.append(len(haplotype2_means[start//args.bin_size:end//args.bin_size]))

            chr_list = [df_chrom['chr'].iloc[0] for ch in range(len(snps_haplotype1_mean))]
            df_means_chr.append(pd.DataFrame(list(zip(chr_list, snps_pos_start, snps_pos_end, snps_haplotype1_mean)),
                                        columns=['chromosome', 'start', 'end', 'state']))
            df_means_chr.append(pd.DataFrame(list(zip(chr_list, snps_pos_start, snps_pos_end, snps_haplotype2_mean)),
                                             columns=['chromosome', 'start', 'end', 'state']))
        else:
            snps_pos_start = []
            snps_pos_end = []
            for i in range(len(snps_haplotype1_pos) - 1):
                snps_pos_start.append(snps_haplotype1_pos[i] if snps_haplotype1_pos[i] < 1 else snps_haplotype1_pos[i] + 1)
                snps_pos_end.append(snps_haplotype1_pos[i + 1])
                if cent_index_hp1 == i:
                    snps_haplotype1_mean[i] = 0
            chr_list = [df_chrom['chr'].iloc[0] for ch in range(len(snps_haplotype1_mean))]
            df_means_chr.append(pd.DataFrame(list(zip(chr_list, snps_pos_start, snps_pos_end, snps_haplotype1_mean)),
                                        columns=['chromosome', 'start', 'end', 'state']))
            snps_pos_start = []
            snps_pos_end = []
            for i in range(len(snps_haplotype2_pos) - 1):
                snps_pos_start.append(snps_haplotype2_pos[i] if snps_haplotype2_pos[i] < 1 else snps_haplotype2_pos[i] + 1)
                snps_pos_end.append(snps_haplotype2_pos[i + 1])
                if cent_index_hp2 == i:
                    snps_haplotype2_mean[i] = 0
            chr_list = [df_chrom['chr'].iloc[0] for ch in range(len(snps_haplotype2_mean))]
            df_means_chr.append(pd.DataFrame(list(zip(chr_list, snps_pos_start, snps_pos_end, snps_haplotype2_mean)),
                                             columns=['chromosome', 'start', 'end', 'state']))

        #remove cent segment
        chroms = get_contigs_list(args.contigs)
        if not df_chrom['chr'].iloc[0] == chroms[1]:
            del snps_haplotype1_mean[cent_index_hp1]
            del snps_haplotype2_mean[cent_index_hp2]
            del snps_haplotype1_len[cent_index_hp1]
            del snps_haplotype2_len[cent_index_hp2]
        return snps_haplotype1_mean + snps_haplotype2_mean, snps_haplotype1_len + snps_haplotype2_len, df_means_chr

def change_point_detection_algo(bin_size, hp_data, ref_start_values, args, ref_start_values_1, df_centm_chrom):
    ############################################
    zeros_values = []
    for i in range(len(hp_data)):
        if hp_data[i] < 1:
            zeros_values.append(i)

    from itertools import groupby, count
    groups = groupby(zeros_values, key=lambda item, c=count(): item - next(c))
    tmp = [list(g) for k, g in groups]

    cpd_zeros_points = []
    for i, val in enumerate(tmp):
        if len(val) == 1:
            hp_data[i] = statistics.mean(hp_data[val[0]:val[0] + 2])
        else:
            cpd_zeros_points.append(val[0])
            cpd_zeros_points.append(val[-1])
    ####################################################
    if not args.breakpoints:
        data = np.array(hp_data, dtype='int')  # numpy.clip(data, a_min=0, a_max=1000)
        algo = rpt.Pelt(model="rbf", jump=25).fit(data)
        result = algo.predict(pen=10)
        change_points = [i for i in result if i < len(data)]
        ############################################
        ov_indices = []
        for i, val1 in enumerate(change_points):
            for j, val2 in enumerate(tmp):
                if len(val2) > 2:
                    if val2[-1] >= val1 >= val2[0]:
                        ov_indices.append(i)
        if ov_indices:
            for index in sorted(list(set(ov_indices)), reverse=True):
                del change_points[index]
        ############################################
        change_points = sorted(list(set(change_points + cpd_zeros_points + [len(hp_data) - 1])))
   ############################################
    cent = [df_centm_chrom.start.values.tolist()[0], df_centm_chrom.end.values.tolist()[0]]
   ############################################
    cent_cpds = []
    if args.breakpoints:
        change_points = []
        for i, val in enumerate(ref_start_values_1):
            for j, val1 in enumerate(ref_start_values):
                if val1 > val:
                    change_points.append(j)
                    break
    for k, cen in enumerate(cent):
        for l, cen1 in enumerate(ref_start_values):
            if cen1 > cen:
                cent_cpds.append(l)
                break
    cent_cpds_updated = [i for i in cent_cpds if i > 1]
    change_points = [i for i in change_points if i < cent_cpds[0] or i > cent_cpds[1]]
    change_points = sorted(list(set(change_points + cent_cpds_updated + [len(hp_data) - 1])))
    ############################################
    from smoothing import smooth_triangle
    hp_data_new = []

    if len(hp_data[0:change_points[0]]) > 46:
        hp_data_new.extend(smooth_triangle(hp_data[0:change_points[0]], 15))
    else:
        hp_data_new.extend(hp_data[0:change_points[0]])

    for p in range(1, len(change_points)):
        if change_points[p] - change_points[p - 1] > 46:
            hp_data_new.extend(smooth_triangle(hp_data[change_points[p - 1]:change_points[p]], 15))
        else:
            hp_data_new.extend(hp_data[change_points[p - 1]:change_points[p]])

    if len(hp_data[change_points[-1]:]) > 46:
        hp_data_new.extend(smooth_triangle(hp_data[change_points[-1]:], 15))
    else:
        hp_data_new.extend(hp_data[change_points[-1]:])

    strong_candidates = []
    snps_haplotype_mean = []
    snps_haplotype_pos = []
    snps_haplotype_len = []
    start = 0
    snps_haplotype_pos.append(0)

    #snps_haplotype_pos = []
    #change_points = sorted(set(add_intermediaries([0] + change_points, 40)))
    for index, point in enumerate(change_points):
        sub_list = hp_data_new[start:point]
        if sub_list:
            median = statistics.median(sub_list)
            # TODO make cent_cpds segments as 0
            if (index == 0 and point == cent_cpds[0]):
                snps_haplotype_mean.append(statistics.mean(remove_outliers_iqr(np.array(sub_list))))
            elif point in [cent_cpds[0]] or (index == 0 and cent_cpds[0] == 1):
               snps_haplotype_mean.append(0)
            else:
                snps_haplotype_mean.append(statistics.mean(remove_outliers_iqr(np.array(sub_list))))
        else:
            snps_haplotype_mean.append(0)
        snps_haplotype_len.append(len(sub_list))
        snps_haplotype_pos.append(point * bin_size)
        start = point
        ###################################
        if args.cpd_internal_segments:
            start_sub = change_points[index-1]
            if index > 0 and len(sub_list) > (1000000 // args.bin_size):
                sub_list_data = np.array(sub_list, dtype='int')  # numpy.clip(data, a_min=0, a_max=1000)
                algo_sub_list = rpt.Pelt(model="rbf", jump=15).fit(sub_list_data)
                result_sub_list = algo_sub_list.predict(pen=5)
                change_points_sub_list = [i for i in result_sub_list if i < len(sub_list_data)]
                change_points_sub_list = [x+ change_points[index-1] for x in change_points_sub_list + [len(sub_list)+1]]

                if len(change_points_sub_list) > 2:
                    del snps_haplotype_mean[-1]
                    del snps_haplotype_len[-1]
                    del snps_haplotype_pos[-1]

                    for index_sub, point_sub in enumerate(change_points_sub_list):
                        sub_list_sub = hp_data_new[start_sub:point_sub]
                        if sub_list_sub:
                            snps_haplotype_mean.append(statistics.median(remove_outliers_iqr(np.array(sub_list))))
                        else:
                            snps_haplotype_mean.append(0)
                        snps_haplotype_len.append(len(sub_list_sub))

                        start_sub = point_sub
                        snps_haplotype_pos.append(point_sub * bin_size)
        ###################################

    return snps_haplotype_mean, snps_haplotype_len, snps_haplotype_pos

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
def merge_adjacent_regions_cn(segarr, args):
    chroms = get_contigs_list(args.contigs)
    dfs = []
    for index, chrom in enumerate(chroms):
        seg = segarr[segarr['chromosome'] == chrom]
        label_groups = seg['state'].ne(seg['state'].shift()).cumsum()
        df = (seg.groupby(label_groups).agg({'chromosome': 'first', 'start': 'min', 'end': 'max', 'depth': 'mean', 'state': 'first'}).reset_index(drop=True))
        dfs.append(df)
    out = pd.concat(dfs)
    return out

def merge_adjacent_regions_loh(segarr, args):
    chroms = get_contigs_list(args.contigs)
    dfs = []
    for index, chrom in enumerate(chroms):
        seg = segarr[segarr['chr'] == chrom]
        seg['diff'] = seg['start'] - seg['end'].shift().fillna(seg['start'].iloc[0])
        seg['group'] = (seg['diff'] != 1).cumsum()
        df = (seg.groupby(seg['group']).agg({'chr': 'first', 'start': 'min', 'end': 'max'}).reset_index(drop=True))
        dfs.append(df)
    out = pd.concat(dfs)
    return out

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
def adjust_extreme_outliers(hp_data):
    for j, val in enumerate(hp_data):
        if val > 500 and j > 0 and j < len(hp_data):
            hp_data[j] = hp_data[j-1]
        elif val > 500 and j == 0:
            hp_data[j] = hp_data[j+1]
    return hp_data

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

def get_vafs_from_normal_phased_vcf(df_snps, df_coverages, chroms, args):
    df_final = []

    for index, chrom in enumerate(chroms):
        df = df_snps[df_snps['chr'] == chrom]
        df_coverage = df_coverages[df_coverages['chr'] == chrom]
        starts_pos = df_coverage.start.values.tolist()
        vaf = []
        # df = dict(tuple(df_snps.groupby('hp')))
        haplotype_1_position = df.pos.values.tolist()
        haplotype_1_coverage = df.freq_value_b.values.tolist()
        haplotype_2_position = df.pos.values.tolist()
        haplotype_2_coverage = df.freq_value_a.values.tolist()
        # vaf = a/a+b
        vaf = [round(min(i,j) / (i + j + 0.0000001), 3) for i, j in zip(haplotype_1_coverage, haplotype_2_coverage)]
        #vaf = list(map(lambda x: 1 - x if x > 0.5 else x, vaf))

        snps_het_vaf = []
        for index, pos in enumerate(starts_pos):
            l2 = [i for i in haplotype_1_position if i > pos and i < pos + args.bin_size]
            if l2:
                snps_het_vaf.append(vaf[haplotype_1_position.index(max(l2))])
            else:
                snps_het_vaf.append(0)

        df_final.append(pd.DataFrame(list(zip([df['chr'].iloc[0] for ch in range(len(starts_pos))], starts_pos, snps_het_vaf)), columns=['chr', 'pos', 'vaf']))

    return pd.concat(df_final)

def get_vafs_from_tumor_phased_vcf(df_snps, df_coverages, chroms, args):
    df_final = []

    for index, chrom in enumerate(chroms):
        df = df_snps[df_snps['chr'] == chrom]
        df_coverage = df_coverages[df_coverages['chr'] == chrom]
        starts_pos = df_coverage.start.values.tolist()
        vaf = []
        # df = dict(tuple(df_snps.groupby('hp')))

        snps_df_haplotype1 = df[(df['gt'] == '0|1') | (df['gt'] == '0/1') | (df['gt'] == '1|0') | (df['gt'] == '1/0')]
        snps_df_haplotype1.reindex(snps_df_haplotype1)
        if df.vaf.dtype == object:
            haplotype_1_coverage = [eval(i) for i in df.vaf.str.split(',').str[0].values.tolist()]
        else:
            haplotype_1_coverage = df.vaf.values.tolist()
        haplotype_1_position = snps_df_haplotype1.pos.values.tolist()

        # snps_df_haplotype1 = df[(df['gt'] == '1|0') | (df['gt'] == '1/0')]
        # snps_df_haplotype1.reindex(snps_df_haplotype1)
        # if df.vaf.dtype == object:
        #     haplotype_2_coverage = [eval(i) for i in df.vaf.str.split(',').str[0].values.tolist()]
        # else:
        #     haplotype_2_coverage = df.vaf.values.tolist()
        # haplotype_2_position = snps_df_haplotype1.pos.values.tolist()

        # vaf = a/a+b
        #vaf = [round(min(i,j) / (i + j + 0.0000001), 3) for i, j in zip(haplotype_1_coverage, haplotype_2_coverage)]

        haplotype_1_coverage = list(map(lambda x: 1 - x if x > 0.5 else x, haplotype_1_coverage))

        snps_het_vaf = []
        for index, pos in enumerate(starts_pos):
            l2 = [i for i in haplotype_1_position if i > pos and i < pos + args.bin_size]
            #print(haplotype_1_position.index(min(l2)), ':', haplotype_1_position.index(max(l2)))
            #data = haplotype_1_coverage[haplotype_1_position.index(min(l2)):haplotype_1_position.index(max(l2))]
            if len(l2) == 1:
                snps_het_vaf.append(haplotype_1_coverage[haplotype_1_position.index(l2[0])])
            elif len(l2) > 1:
                snps_het_vaf.append(statistics.mean(haplotype_1_coverage[haplotype_1_position.index(min(l2)):haplotype_1_position.index(max(l2))]))
            else:
                snps_het_vaf.append(0)

        df_final.append(pd.DataFrame(list(zip([df['chr'].iloc[0] for ch in range(len(starts_pos))], starts_pos, snps_het_vaf)), columns=['chr', 'pos', 'vaf']))

    return pd.concat(df_final)

def parse_sv_vcf(path):
    from vcf_parser import VCFParser
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

    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)

    cent_indices = []
    for index, row in df_1.iterrows():
        for index_cent, row_cent in df_centm.iterrows():
            if row['chromosome'] == row_cent['chr'] and (((row['start'] - 1) == row_cent['start'] or (row['start']) == row_cent['start']) and row_cent['end'] == row_cent['end']):
                cent_indices.append(index)
    df_1 = df_1.drop(cent_indices)

    cent_indices = []
    for index, row in df_2.iterrows():
        for index_cent, row_cent in df_centm.iterrows():
            if row['chromosome'] == row_cent['chr'] and (((row['start'] - 1) == row_cent['start'] or (row['start']) == row_cent['start']) and row_cent['end'] == row_cent['end']):
                cent_indices.append(index)
    df_2 = df_2.drop(cent_indices)

    for i in range(len(integer_fractional_means)):
        df_1['state'].mask(df_1['state'] == centers[i], integer_fractional_means[i], inplace=True)

    vals = df_1.state.values.tolist()
    weights = [i-j for i,j in zip(df_1.end.values.tolist(), df_1.start.values.tolist())]

    ploidy_hp1 = weighted_means(vals, weights)

    for i in range(len(integer_fractional_means)):
        df_2['state'].mask(df_2['state'] == centers[i], integer_fractional_means[i], inplace=True)

    vals = df_2.state.values.tolist()
    weights = [i-j for i,j in zip(df_2.end.values.tolist(), df_2.start.values.tolist())]

    ploidy_hp2 = weighted_means(vals, weights)

    return ploidy_hp1 + ploidy_hp2

def weighted_means(vals, weights):
    weighted_vals = []
    vals_n_weights = [(vals[i], weights[i]) for i in range(0, len(weights))]
    for tup in vals_n_weights:
        weighted_vals.append(round(tup[0] * tup[1] / sum(weights), 3))
    return sum(weighted_vals)
def average_p_value_genome(args, centers, df_segs_hp1, df_segs_hp2, df_hp1, df_hp2):
    chroms = get_contigs_list(args.contigs)
    hp_1_values = []
    hp_2_values = []
    for i in range(len(centers)):
        hp_1_values.append([])
        hp_2_values.append([])

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
        if hp_1_values[i] + hp_2_values[i]:
            sample_mean.append(statistics.mean(remove_outliers_iqr(np.array(hp_1_values[i] + hp_2_values[i]))))
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
            if end - start * args.bin_size > 5000000:
                seg_mean = statistics.median(remove_outliers_iqr(np.array(df_hp_1_val[start//args.bin_size:end//args.bin_size])))
                sample_mean_init = min(sample_mean, key=lambda x: abs(x - seg_mean))
                index =  sample_mean.index(sample_mean_init)
                z_score =  (seg_mean - sample_mean_init) / 10 #statistics.stdev(remove_outliers_iqr(np.array(df_hp_1_val[start//args.bin_size:end//args.bin_size])))
                p_value = stats.norm.sf(abs(z_score)) * 2
                if p_value >= 0:
                    df_segs_hp_1_updated_p_score.append(round(p_value, 7))
                    df_segs_hp_1_updated_weight.append(end-start)

        for i, (start,end) in enumerate(zip(df_segs_hp_2_updated_start, df_segs_hp_2_updated_end)):
            if end-start * args.bin_size > 5000000:
                seg_mean = statistics.median(remove_outliers_iqr(np.array(df_hp_2_val[start//args.bin_size:end//args.bin_size])))
                sample_mean_init = min(sample_mean, key=lambda x: abs(x - seg_mean))
                index =  sample_mean.index(sample_mean_init)
                z_score =  (seg_mean - sample_mean_init) / 10 #statistics.stdev(remove_outliers_iqr(np.array(df_hp_2_val[start//args.bin_size:end//args.bin_size])))
                p_value = stats.norm.sf(abs(z_score)) * 2
                if p_value >= 0:
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
    # import matplotlib.pyplot as plt
    # plt.plot(t, data)
    # plt.axhline(limit, color='r')
    # plt.plot(t[ind], data[ind], 'ro')
    # plt.title('Peaks: minimum value {limit}, minimum spacing {spacing} points'.format(**{'limit': limit, 'spacing': spacing}))
    # plt.savefig(args.out_dir_plots + '/' + args.genome_name + '_' + "normal_optimized_peak.pdf", format="pdf", bbox_inches="tight")

    if len(ind):
        normals =  [i for i in list(ind)]
    else:
        normals = [0]
        logging.info('No normal peak value detected under confidence [%.2f], so assuming only first estimated value', args.purity_range, args.ploidy_range, limit)

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

        if args.breakpoints:
            for index, (state1, state2, start, end) in enumerate(zip(hp1_state, hp2_state, hp1_start, hp1_end)):
                #print(chrom, state1, state2, start, end)
                if state1 == 0 or state2 == 0:# and end - start > 2000000:
                    loh_region_starts.append(start)
                    loh_region_ends.append(end)
        else:
            for index, (state1, start, end) in enumerate(zip(hp1_state, hp1_start, hp1_end)):
                if state1 == 0:# and end - start > 2000000:
                    loh_region_starts.append(start)
                    loh_region_ends.append(end)

            for index, (state2, start, end) in enumerate(zip(hp2_state, hp2_start, hp2_end)):
                if state2 == 0:# and end - start > 2000000:
                    if not start in loh_region_starts and not end in loh_region_ends:
                        loh_region_starts.append(start)
                        loh_region_ends.append(end)

            # sort
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

def update_genes_phase_corrected_coverage(args, df_segs_hp1, df_segs_hp2, p_value, centers, integer_fractional_centers):

    df_genes = csv_df_chromosomes_sorter(args.out_dir_plots + '/data_phasing/cancer_genes_coverage.csv', ['chr','start','end','gene', 'hp1', 'hp2'])
    #write_df_csv(df_genes, args.out_dir_plots + '/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output/' + 'cancer_genes_coverage.csv')
    df_genes = genes_phase_correction(df_genes, df_segs_hp1, df_segs_hp2, args, centers, integer_fractional_centers)
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