import statistics

import numpy as np
import pandas as pd
import pysam
import os
from itertools import takewhile
import logging
import csv
import statistics
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from sklearn.neighbors import KernelDensity

logger = logging.getLogger()

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from src.file_tools.process_bam import get_segments_coverage

from src.coverage.segmentation import split_regions_by_points
from src.utils_tmp.chromosome import csv_df_chromosomes_sorter, df_chromosomes_sorter, get_contigs_list


def get_chromosomes_bins_hapcorrect(bam_file, bin_size, args):
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


def write_segments_coverage(coverage_segments, output, args):
    with open(args.out_dir_plots+'/coverage_data/' + output, 'a') as fp:
        for items in coverage_segments:
            if not items == None:
                fp.write("%s\n" % items)


def loh_regions_events(chrom, region_starts, region_ends, hp):
    dict = []
    for i in range(len(region_starts)):
        dict.append((chrom + '\t' + str(region_starts[i]) + '\t' + str(region_ends[i]) + '\t' + str(hp[i])))
    #write_segments_coverage(dict, args.genome_name'] + '_loh_segments.bed')
    return dict


def detect_alter_loh_regions(csv_df_coverage_tumor_chrom, args, event, chrom, ref_ends, haplotype_1_values, haplotype_2_values,
                             unphased_reads_values, starts, ends, switch_hps):
    region_starts = []
    region_ends = []
    hp = []

    if not args.without_phasing:
        histo_unphased_reads_values = csv_df_coverage_tumor_chrom.hp3.values.tolist()
        histo_haplotype_1_values = csv_df_coverage_tumor_chrom.hp1.values.tolist()
        histo_haplotype_2_values = csv_df_coverage_tumor_chrom.hp2.values.tolist()

    if len(ref_ends) == 0:
        return haplotype_1_values, haplotype_2_values, unphased_reads_values, region_starts, region_ends, hp
    if ends and ends[-1] > ref_ends[-1]:
        ends[-1] = ref_ends[-1]

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
            if ends - starts > 20000:
                #TODO Discuss with Ayse, alternate approach on what HP should be selected for each region
                #if mean_values(haplotype_1_values, starts//args.bin_size - 4, starts//args.bin_size - 1) > mean_values(haplotype_2_values, starts//args.bin_size - 4, starts//args.bin_size - 1):
                for i in range(starts//args.bin_size,ends//args.bin_size):
                    haplotype_1_values[i] = histo_unphased_reads_values[i] + histo_haplotype_1_values[i] + histo_haplotype_2_values[i] #haplotype_1_values[i] + haplotype_2_values[i] + unphased_reads_values[i]
                    haplotype_2_values[i] = 0.0001
                    unphased_reads_values[i] = 0
                hp.append(1)
                # else:
                #     for i in range(starts // args.bin_size, ends // args.bin_size):
                #         haplotype_2_values[i] = haplotype_1_values[i] + haplotype_2_values[i] + unphased_reads_values[i]
                #         haplotype_1_values[i] = 0
                #         unphased_reads_values[i] = 0
                #     hp.append(2)

    return haplotype_1_values, haplotype_2_values, unphased_reads_values, region_starts, region_ends, hp


def adjust_loh_cent_phaseblocks(args, loh_region_starts, loh_region_ends, centrom_region, snps_haplotype1_mean, snps_haplotype2_mean,
                                haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):

    #remove centrom phaseblocks
    if not centrom_region.empty:
        cents = [centrom_region.start.values.tolist()[0], centrom_region.end.values.tolist()[0]]
    else:
        cents = [0,0]
    #break phaseblocks at LOH/Cents
    ref_start_values_phasesets,  ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = \
            split_regions_by_points(ref_start_values_phasesets,  ref_end_values_phasesets,
                                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, cents)
    ref_start_values_phasesets,  ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = \
            split_regions_by_points(ref_start_values_phasesets,  ref_end_values_phasesets, haplotype_1_values_phasesets,
                                    haplotype_2_values_phasesets, loh_region_starts+loh_region_ends)

    index_cents = []
    if len(ref_start_values_phasesets) > 1:
        for i in range(len(ref_start_values_phasesets)):
            if ref_start_values_phasesets[i] >= cents[0] and ref_end_values_phasesets[i] <= cents[1]:
                index_cents.append(i)

            for ps, (start,end) in enumerate(zip(loh_region_starts, loh_region_ends)):
                if ref_start_values_phasesets[i] >= start and ref_end_values_phasesets[i] <= end:
                    index_cents.append(i)

    index_cents = list(set(index_cents))

    #extend LOH phaseblocks
    haplotype_1_values_phasesets_loh = []
    haplotype_2_values_phasesets_loh = []
    ref_start_values_phasesets_loh = []
    ref_end_values_phasesets_loh = []

    for ps, (start, end) in enumerate(zip(loh_region_starts, loh_region_ends)):
        haplotype_1_values_phasesets_loh.append(mean_values(snps_haplotype1_mean, start//args.bin_size, end//args.bin_size))
        haplotype_2_values_phasesets_loh.append(mean_values(snps_haplotype2_mean, start//args.bin_size, end//args.bin_size))
        ref_start_values_phasesets_loh.append(start)
        ref_end_values_phasesets_loh.append(end)

    haplotype_1_values_phasesets = [i for j, i in enumerate(haplotype_1_values_phasesets) if j not in index_cents] + haplotype_1_values_phasesets_loh
    haplotype_2_values_phasesets = [i for j, i in enumerate(haplotype_2_values_phasesets) if j not in index_cents] + haplotype_2_values_phasesets_loh
    ref_start_values_phasesets = [i for j, i in enumerate(ref_start_values_phasesets) if j not in index_cents] + ref_start_values_phasesets_loh
    ref_end_values_phasesets = [i for j, i in enumerate(ref_end_values_phasesets) if j not in index_cents] + ref_end_values_phasesets_loh

    zipped = list(zip(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets))
    zipped_sorted = sorted(zipped, key=lambda x: x[0])
    ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = zip(*zipped_sorted)

    return list(haplotype_1_values_phasesets), list(haplotype_2_values_phasesets), list(ref_start_values_phasesets), list(ref_end_values_phasesets)


def mean_values(selected_list, start_index, end_index):
    selected_list = [i for i in selected_list[start_index:end_index] if i > 0]
    if len(selected_list) >= 2:

        return np.median(selected_list)
    elif len(selected_list) == 1:
        return selected_list[0]
    else:
        return 0.0


def loh_regions_phasesets(haplotype_1_values, haplotype_2_values, loh_region_starts, loh_region_ends, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, args):
    indices = []
    for l, (loh_start, loh_end) in enumerate(zip(loh_region_starts, loh_region_ends)):
        for k, (ps_start, ps_end) in enumerate(zip(ref_start_values_phasesets, ref_end_values_phasesets)):
            if (ps_start >= loh_start and ps_end <= loh_end) or (loh_start < ps_end and ps_start < loh_end):
                indices.append(k)

    haplotype_1_values_phasesets = [j for i, j in enumerate(haplotype_1_values_phasesets) if i not in indices]
    haplotype_2_values_phasesets = [j for i, j in enumerate(haplotype_2_values_phasesets) if i not in indices]
    ref_start_values_phasesets = [j for i, j in enumerate(ref_start_values_phasesets) if i not in indices]
    ref_end_values_phasesets = [j for i, j in enumerate(ref_end_values_phasesets) if i not in indices]

    for m, (loh_start, loh_end) in enumerate(zip(loh_region_starts, loh_region_ends)):
        haplotype_1_values_phasesets.append(mean_values(haplotype_1_values, loh_start//args.bin_size, loh_end//args.bin_size))
        haplotype_2_values_phasesets.append(mean_values(haplotype_2_values, loh_start//args.bin_size, loh_end//args.bin_size))
        ref_start_values_phasesets.append(loh_start)
        ref_end_values_phasesets.append(loh_end)

    sort_function = lambda x: x[0]
    sort_target = list(zip(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets))
    sort_target.sort(key=sort_function)

    ref_start_values_phasesets = [a for a, b, c, d in sort_target]
    ref_end_values_phasesets = [b for a, b, c, d in sort_target]
    haplotype_1_values_phasesets = [c for a, b, c, d in sort_target]
    haplotype_2_values_phasesets = [d for a, b, c, d in sort_target]

    return haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets


def update_hp_assignment_loh_segments(args, loh_df, coverage_df):
    updated_loh_segs = []
    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        df_loh_chrom = loh_df[loh_df['chr'] == chrom]
        df_coverage_chrom = coverage_df[coverage_df['chr'] == chrom]
        hp2_cov = df_coverage_chrom.hp2.values.tolist()
        hp2_start = df_coverage_chrom.start.values.tolist()
        hp2_end = df_coverage_chrom.end.values.tolist()
        for idx, seg in df_loh_chrom.iterrows():
            #print(chrom, df_loh_chrom.loc[idx, 'start'], df_loh_chrom.loc[idx, 'end'], hp2_start.index(df_loh_chrom.loc[idx, 'start']), hp2_end.index(df_loh_chrom.loc[idx, 'end']))
            if statistics.mean(hp2_cov[(df_loh_chrom.loc[idx, 'start']//args.bin_size)+1:(df_loh_chrom.loc[idx, 'end']//args.bin_size)-1]) > 0.0001:
                df_loh_chrom.loc[idx, 'hp'] = 2

        updated_loh_segs.append(df_loh_chrom)

    segs_loh = pd.concat(updated_loh_segs)

    return segs_loh
