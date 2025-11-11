import statistics

import numpy as np
import pandas as pd
import pysam
import os
import logging
logger = logging.getLogger()

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from src.extras import get_contigs_list
from src.hapcorrect.src.process_bam import get_segments_coverage

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
                    bed.append([tail, chrs[i], start, region[i]])
                else:
                    bed.append([tail, chrs[i], start, end])
                start=end+1
                end+=bin_size
    return bed

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

def csv_df_chromosomes_sorter(path, names, sept='\t'):
    dataframe = pd.read_csv(path, sep=sept, names=names)
    dataframe['chr'] = dataframe['chr'].astype(str)
    #if not dataframe['chr'].iloc[0].startswith('chr'):
    #    dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', names[1]], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

def df_chromosomes_sorter(dataframe, names, sept='\t'):
    dataframe['chr'] = dataframe['chr'].astype(str)
    #if not dataframe['chr'].iloc[0].startswith('chr'):
    #    dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', names[1]], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

def write_df_csv(df, file_name):
    df.to_csv(file_name, sep='\t', index=False, header=False, mode='w')

def write_segments_coverage(coverage_segments, output, args):
    with open(args.out_dir_plots+'/coverage_data/' + output, 'a') as fp:
        for items in coverage_segments:
            if not items == None:
                fp.write("%s\n" % items)

def write_segments_coverage_snps(coverage_segments, output, args):
    with open(args.out_dir_plots+'/data_phasing/' + output, 'a') as fp:
        for items in coverage_segments:
            if not items == None:
                fp.write("%s\n" % items)

def loh_regions_events(chrom, region_starts, region_ends, hp):
    dict = []
    for i in range(len(region_starts)):
        dict.append((chrom + '\t' + str(region_starts[i]) + '\t' + str(region_ends[i]) + '\t' + str(hp[i])))
    #write_segments_coverage(dict, args.genome_name'] + '_loh_segments.bed')
    return dict

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


def get_snps_frquncies_coverage_from_bam(df, chrom):
    df = df[df['chr'] == chrom]
    #df = dict(tuple(df.groupby('hp')))
    haplotype_1_position = df.pos.values.tolist()
    haplotype_1_coverage = df.freq_value_b.values.tolist()
    haplotype_2_position = df.pos.values.tolist()
    haplotype_2_coverage = df.freq_value_a.values.tolist()

    return haplotype_1_position, haplotype_1_coverage, haplotype_2_position, haplotype_2_coverage


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

def split_regions_by_points(starts, ends, values1, values2, split_points):
    split_points = sorted(set(split_points))
    new_starts = []
    new_ends = []
    new_values1 = []
    new_values2 = []

    for start, end, val1, val2 in zip(starts, ends, values1, values2):
        # Get split points within this region
        internal_splits = [pt for pt in split_points if start < pt < end]
        boundaries = [start] + internal_splits + [end]

        for i in range(len(boundaries) - 1):
            new_starts.append(boundaries[i])
            new_ends.append(boundaries[i + 1])
            new_values1.append(val1)
            new_values2.append(val2)

    return new_starts, new_ends, new_values1, new_values2


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


def extract_centromere_regions(args):
    #fileDir = args.centromere #os.path.dirname(__file__)
    #cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(args.centromere, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)

    return df_centm

def snps_frequencies_chrom_genes(df_snps_frequencies, args):
    df_chroms = []

    df_genes_all = csv_df_chromosomes_sorter(args.cancer_genes, ['chr', 'start', 'end', 'gene'])
    df_empty = pd.DataFrame(columns=['chr', 'start', 'end', 'gene', 'hp1', 'hp2'])
    if args.user_input_genes:
        with open(args.user_input_genes, 'r') as file:
            entries = file.readlines()
        if "\t" in entries[0]:
            entries = [line.rstrip('\n').split('\t')[3] for line in entries]
        else:
            entries = [line.rstrip('\n') for line in entries]
        if args.reference_name:#not entries[0] in ['chm13','grch38']:
            ref = args.reference_name #default
        else:
            ref = 'grch38'
        import csv
        prefix, filename = os.path.split(args.cancer_genes)
        ref_name = prefix + '/' + ref + '_genes.tsv'
        chroms = []
        starts = []
        ends = []
        genes = []
        with open(ref_name, 'r') as file:
            tsv_reader = csv.reader(file, delimiter='\t')
            for row in tsv_reader:
                if row[3] in entries:
                    chroms.append(row[0])
                    starts.append(row[1])
                    ends.append(row[2])
                    genes.append(row[3])

        data = {"chr": chroms, "start": starts, "end": ends, "gene": genes}
        df_user_genes = pd.DataFrame(data)
        if not df_user_genes.empty:
            write_df_csv(df_user_genes, args.out_dir_plots + '/data_phasing/user_genes.csv')
            df_genes_all = csv_df_chromosomes_sorter(args.out_dir_plots + '/data_phasing/user_genes.csv', ['chr', 'start', 'end', 'gene'])

    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        df = df_snps_frequencies[df_snps_frequencies['chr'] == chrom]
        df_genes = df_genes_all[df_genes_all['chr'] == chrom]

        if not df_genes.empty:

            # df = dict(tuple(df_snps.groupby('hp')))
            haplotype_1_position = df.pos.values.tolist()
            haplotype_1_coverage = df.freq_value_b.values.tolist()
            haplotype_2_position = df.pos.values.tolist()
            haplotype_2_coverage = df.freq_value_a.values.tolist()

            snps_haplotype1_mean = []
            for index, (i,j) in enumerate(zip(df_genes.start.values.tolist(), df_genes.end.values.tolist())):
                sub_list = []
                try:
                    sub_list = haplotype_1_coverage[haplotype_1_position.index(
                        min(haplotype_1_position, key=lambda x:abs(x-i))):haplotype_1_position.index(
                        min(haplotype_1_position, key=lambda x:abs(x-j)))]
                except ValueError:
                    logger.info('No Hets pileup found!')
                if sub_list:
                    snps_haplotype1_mean.append(statistics.mean(sub_list))
                else:
                    snps_haplotype1_mean.append(0)

            snps_haplotype2_mean = []
            for index, (i, j) in enumerate(zip(df_genes.start.values.tolist(), df_genes.end.values.tolist())):
                sub_list = []
                try:
                    sub_list = haplotype_2_coverage[haplotype_2_position.index(
                        min(haplotype_2_position, key=lambda x: abs(x - i))):haplotype_2_position.index(
                        min(haplotype_2_position, key=lambda x: abs(x - j)))]
                except ValueError:
                    logger.info('No Hets pileup found!')
                if sub_list:
                    snps_haplotype2_mean.append(statistics.mean(sub_list))
                else:
                    snps_haplotype2_mean.append(0)

            #snps_mean = [round(i + j, 2) for i, j in zip(snps_haplotype1_mean, snps_haplotype2_mean)]
            if len(snps_haplotype1_mean) == 0:
                df_chroms.append(df_empty)
            else:
                df_chroms.append(pd.DataFrame(list(zip(df_genes.chr.values.tolist(), df_genes.start.values.tolist(), df_genes.end.values.tolist(), \
                                         df_genes.gene.values.tolist(), snps_haplotype1_mean, snps_haplotype2_mean)),
                                 columns=['chr', 'start', 'end', 'gene', 'hp1', 'hp2']))
    if len(df_chroms):
        return pd.concat(df_chroms)
    else:
        return df_empty

def genes_segments_list(bam, args):
    head, tail = os.path.split(bam)

    df_genes = csv_df_chromosomes_sorter(args.cancer_genes, ['chr', 'start', 'end', 'gene'])
    starts = df_genes.start.values.tolist()
    ends = df_genes.end.values.tolist()

    bed = []
    for ind, chrom in enumerate(df_genes.chr.values.tolist()):
        bed.append([tail, chrom, starts[ind], ends[ind]])

    return bed
def genes_segments_coverage(genes_coverage, args):

    df_genes = csv_df_chromosomes_sorter(args.cancer_genes, ['chr', 'start', 'end', 'gene'])

    hp1 = []
    hp2 = []
    hp3 = []
    for items in genes_coverage:
        hp1.append(round(float(items.split('\t')[3]), 2))
        hp2.append(round(float(items.split('\t')[4]), 2))
        #hp3.append(round(float(items.split('\t')[5]), 2))

    return pd.DataFrame(list(zip(df_genes.chr.values.tolist(), df_genes.start.values.tolist(), df_genes.end.values.tolist(),
                                 df_genes.gene.values.tolist(), hp1, hp2)),
                         columns=['chr', 'start', 'end', 'gene', 'hp1', 'hp2'])


def mean_values(selected_list, start_index, end_index):
    selected_list = [i for i in selected_list[start_index:end_index] if i > 0]
    if len(selected_list) >= 2:

        return np.median(selected_list)
    elif len(selected_list) == 1:
        return selected_list[0]
    else:
        return 0.0

def infer_missing_phaseblocks(ref_start_values, ref_end_values, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values, bin_size):
    import statistics
    missing_segments_starts = []
    missing_segments_ends = []
    missing_segments_hp1_value = []
    missing_segments_hp2_value = []
    for i in range(len(ref_start_values_phasesets)-1):
        if ref_start_values_phasesets[i +1] - ref_end_values_phasesets[i] > bin_size * 6:
            start = ((ref_end_values_phasesets[i]//bin_size) +1) * bin_size + 1
            end = (ref_start_values_phasesets[i+1]//bin_size) * bin_size
            missing_segments_starts.append(start)#(ref_end_values_phasesets[i] + 1)
            missing_segments_ends.append(end)#(ref_start_values_phasesets[i+1] - 1)
            try:
                missing_segments_hp1_value.append(statistics.mean(list(filter(lambda zer: zer != 0, haplotype_1_values[ref_start_values.index(start):ref_end_values.index(end)]))))
            except:
                missing_segments_hp1_value.append(0)
            try:
                missing_segments_hp2_value.append(statistics.mean(list(filter(lambda zer: zer != 0, haplotype_2_values[ref_start_values.index(start):ref_end_values.index(end)]))))
            except:
                missing_segments_hp2_value.append(0)

    if len(missing_segments_starts):
        ref_start_values_phasesets = ref_start_values_phasesets + missing_segments_starts
        ref_end_values_phasesets = ref_end_values_phasesets + missing_segments_ends
        haplotype_1_values_phasesets = haplotype_1_values_phasesets + missing_segments_hp1_value
        haplotype_2_values_phasesets = haplotype_2_values_phasesets + missing_segments_hp2_value

        sort_function = lambda x: x[0]
        sort_target = list(zip(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets))
        sort_target.sort(key=sort_function)

        ref_start_values_phasesets = [a for a, b, c, d in sort_target]
        ref_end_values_phasesets = [b for a, b, c, d in sort_target]
        haplotype_1_values_phasesets = [c for a, b, c, d in sort_target]
        haplotype_2_values_phasesets = [d for a, b, c, d in sort_target]

    return ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets

def is_phasesets_check_simple_heuristics(ref_start_values_phasesets, ref_end_values_phasesets, args):
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

def overlap_check(start, end, starts, ends):
  for i in range(len(starts)):
    if (start < ends[i] and end > starts[i]) or (start <= starts[i] and end >= ends[i]):
      return True
  return False

def merge_regions_alt(starts, ends, values, values_next, loh_starts, loh_ends, threshold=3):
    merged_starts, merged_ends, merged_values = [], [], []
    for i in range(len(starts)-1):
        while abs(values[i+1] - values[i]) <= threshold and abs(values[i] - values_next[i]) >= threshold:
            if not overlap_check(starts[i], ends[i], loh_starts, loh_ends):
                continue #print('here')
def merge_regions(starts, ends, values, values_next, loh_starts, loh_ends, threshold=3):
    merged_starts, merged_ends, merged_values, merged_values_2 = [], [], [], []
    i = 0
    while i < len(starts):
        start, end, value = starts[i], ends[i], values[i]
        j = i + 1
        while j < len(starts) and abs(value - values[j]) <= threshold and abs(value - values_next[i]) >= threshold:
            end = max(end, ends[j])
            j += 1
        if not overlap_check(start, end, loh_starts, loh_ends):
            merged_starts.append(start)
            merged_ends.append(end)
            merged_values.append(sum(values[i:j]) / (j - i))
            merged_values_2.append(sum(values_next[i:j]) / (j - i))
        i = j
    return merged_starts, merged_ends, merged_values, merged_values_2

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

def get_contigs_list(contigs):
    #chrs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']
    chroms_list_final = []
    chroms = contigs.split(',')
    for chrom in chroms:
        chrom = chrom[len('chr'):] if chrom.startswith('chr') else chrom
        chrom = chrom.split('-')
        if len(chrom) > 1:
            chroms_list_final.extend(list(range(int(chrom[0]), int(chrom[1]) + 1)))
        else:
            chroms_list_final.extend(chrom)

    chroms_list_final = ['chr' + x if chroms[0].startswith('chr') else x for x in map(str, chroms_list_final)]
    return chroms_list_final

def extend_snps_ratios_df(chrom, offset, ref_start_values_updated, snps_het_counts, snps_homo_counts):
    chr_list = [chrom for ch in range(len(ref_start_values_updated))]

    df_snps_ratios_chrom = pd.DataFrame(
        list(zip(chr_list, ref_start_values_updated, snps_het_counts, snps_homo_counts)),
        columns=['chr', 'start', 'hets_ratios', 'homos_ratios'])

    df_snps_ratios_chrom['start_overall'] = df_snps_ratios_chrom['start'].apply(lambda x: x + offset)

    return df_snps_ratios_chrom


def add_breakpoints(args, phasesets_segments, breakpoints):
    bed = []
    head, tail = os.path.split(args.target_bam[0])
    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        original_list_start = []
        original_list_end = []
        df_bps_chrom = breakpoints[breakpoints['chr'] == chrom]
        new_breakpoints = df_bps_chrom.start.values.tolist()

        for k, val in enumerate(phasesets_segments):
            if val[1] == chrom:
                original_list_start.append(val[2])
                original_list_end.append(val[3])

        for i, (start, end) in enumerate(zip(original_list_start, original_list_end)):
            for breakpoint in new_breakpoints:
                if start <= breakpoint <= end and (not breakpoint == start and not breakpoint == end) and (breakpoint - start > args.bin_size * 7 and end - breakpoint > args.bin_size * 7):
                    bed.append([tail, chrom, start, breakpoint])
                    start = breakpoint + 1
            bed.append([tail, chrom, start, end])

    return bed

def find_peak_median_without_outliers(data):
    import numpy as np
    from scipy.stats import gaussian_kde
    # from scipy.signal import find_peaks
    # peaks, _ = find_peaks(data)
    # peak_values = [data[i] for i in peaks]
    # peak_index = peaks[np.argmax(peak_values)]
    # peak_value = data[peak_index]

    if len(data) == 0: #second time almost big bp boundries based PS will be one and will be set to 0, change condition
        return 0
    else:
        return statistics.median(data)

    from scipy.signal import find_peaks
    from sklearn.neighbors import KernelDensity

    #if len(data) < 5:
    #    return statistics.median(data)

    # 2. Remove outliers using IQR method
    # q1 = np.percentile(data, 25)
    # q3 = np.percentile(data, 75)
    # iqr = q3 - q1
    # lower_bound = q1 - 1.5 * iqr
    # upper_bound = q3 + 1.5 * iqr
    # filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
    def remove_outliers_modified_z(data, peak_value, threshold=3.5):
        data = np.array(data)
        median = peak_value #np.median(data)
        mad = np.median(np.abs(data - median))  # Median Absolute Deviation
        if mad == 0:
            return data  # No variation
        modified_z_scores = 0.6745 * (data - median) / mad
        return data[np.abs(modified_z_scores) <= threshold]

    filtered_data = remove_outliers_modified_z(data, peak_value)
    # 3. Compute median of filtered data
    return np.median(filtered_data)

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
