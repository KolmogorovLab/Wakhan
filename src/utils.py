import csv
import numpy as np
import pandas as pd
import pysam
import os
import statistics
from random import randint
import ruptures as rpt

from cnvlib.cmdutil import read_cna
from cnvlib.cnary import CopyNumArray as CNA
from cnvlib import segmentation, coverage, batch, fix, segmetrics, call, scatter
from skgenome import tabio
from smoothing import smoothing
from hmm import call_copynumbers
from extras import get_contigs_list


from pomegranate import HiddenMarkovModel as Model

def get_chromosomes_bins_replica(bam_file, bin_size, arguments):
    bed=[]
    bam_alignment = pysam.AlignmentFile(bam_file)
    headers = bam_alignment.header
    seq_dict = headers['SQ']
    region = [''] * len(seq_dict)
    chrs = [''] * len(seq_dict)
    head, tail = os.path.split(bam_file)
    chroms = get_contigs_list(arguments['contigs'])
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

def get_chromosomes_bins(bam_file, bin_size, arguments):
    bed=[]
    bam_alignment = pysam.AlignmentFile(bam_file)
    headers = bam_alignment.header
    seq_dict = headers['SQ']
    region = [''] * len(seq_dict)
    chrs = [''] * len(seq_dict)
    head, tail = os.path.split(bam_file)
    chroms = get_contigs_list(arguments['contigs'])
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



def update_bins_with_bps(bins):
    size=len(bins)
    new_start=0
    new_end = 0
    indices=[]

    with open('/home/rezkuh/gits/devel/BGA/breakpoints.csv', mode='r') as csvfile:
        segments = csv.reader(csvfile, delimiter='\t')
        #next(segments)
        for s in segments:
            print(s)
            chr_id=s[0]
            bp_start=int(s[1])
            bp_end=int(s[2])
            for index, value in enumerate(bins):
                if (value[0] == chr_id and bp_start > int(value[1]) and bp_start < int(value[2])) or new_start > 0:
                    if new_start==0:
                        new_start=int(value[1])
                    new_end=int(value[2])
                    indices.append(index)
                    if bp_end < new_end:
                        size=len(bins)
                        bins.insert(size, [chr_id, new_start, bp_start - 1])
                        bins.insert(size+1, [chr_id, bp_start, bp_end])
                        bins.insert(size+2, [chr_id, bp_end+1, int(value[2])])
                        break
    for ind in indices:
        del bins[ind]
    return bins

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



# def csv_df_chromosomes_sorter_coverage(path):
#     dataframe = pd.read_csv(path, sep='\t', names=['chr', 'start', 'end', 'coverage'])
#     dataframe['chr'] = dataframe['chr'].astype(str)
#     #if not dataframe['chr'].iloc[0].startswith('chr'):
#     #    dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
#     dataframe.sort_values(by=['chr', 'start'], ascending=[True, True], inplace=True)
#     return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

# def csv_df_chromosomes_sorter_copyratios(path):
#     dataframe = pd.read_csv(path, sep='\t', names=['chr', 'start', 'end', 'gene', 'log2', 'depth', 'probes', 'weight'])
#     dataframe['chr'] = dataframe['chr'].astype(str)
#     #if not dataframe['chr'].iloc[0].startswith('chr'):
#     #    dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
#     dataframe.sort_values(by=['chr', 'start'], ascending=[True, True], inplace=True)
#     return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

# def csv_df_chromosomes_sorter_snps(path, name):
#     dataframe = pd.read_csv(path, sep='\t', names=['chr', 'pos', 'ps'])
#     dataframe['chr'] = dataframe['chr'].astype(str)
#     #if not dataframe['chr'].iloc[0].startswith('chr'):
#     #    dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
#     dataframe.sort_values(by=['chr', 'pos'], ascending=[True, True], inplace=True)
#     return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)
#
# def csv_df_chromosomes_sorter_snps_vcf(path):
#     dataframe = pd.read_csv(path, sep='\t', names=['chr', 'pos', 'qual', 'gt', 'dp', 'vaf'])
#     dataframe['chr'] = dataframe['chr'].astype(str)
#     #if not dataframe['chr'].iloc[0].startswith('chr'):
#     #    dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
#     dataframe.sort_values(by=['chr', 'pos'], ascending=[True, True], inplace=True)
#     return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)
#
# def csv_df_chromosomes_sorter_snps_from_bam(path):
#     dataframe = pd.read_csv(path, sep='\t', names=['chr', 'pos', 'freq_value_a', 'hp_a', 'freq_value_b', 'hp_b'])
#     dataframe['chr'] = dataframe['chr'].astype(str)
#     #if not dataframe['chr'].iloc[0].startswith('chr'):
#     #    dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
#     dataframe.sort_values(by=['chr', 'pos'], ascending=[True, True], inplace=True)
#     return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)
#
# def csv_df_chromosomes_sorter_snps_frequency(path):
#     dataframe = pd.read_csv(path, sep=',', names=['chr', 'start', 'a', 'c', 'g', 't'])
#     dataframe['chr'] = dataframe['chr'].astype(str)
#     #if not dataframe['chr'].iloc[0].startswith('chr'):
#     #    dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
#     dataframe.sort_values(by=['chr', 'start'], ascending=[True, True], inplace=True)
#     return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)
#
# def csv_df_chromosomes_sorter_snps_alts_gts(path):
#     dataframe = pd.read_csv(path, sep='\t', names=['chr', 'start', 'ref', 'alt', 'gt'])
#     dataframe['chr'] = dataframe['chr'].astype(str)
#     #if not dataframe['chr'].iloc[0].startswith('chr'):
#     #    dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
#     dataframe.sort_values(by=['chr', 'start'], ascending=[True, True], inplace=True)
#     return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)


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

def write_segments_coverage_dict(coverage_segments, output):
    with open('data/' + output, 'a') as fp:
        for items in coverage_segments:
            if not items == None:
                fp.write("%s\n" % items)

def write_segments_coverage(coverage_segments, output, arguments):
    with open(arguments['out_dir_plots']+'/bed_output/' + output, 'a') as fp:
        for items in coverage_segments:
            if not items == None:
                fp.write("%s\n" % items)

def write_header_comments(header, header_comments, output, arguments):
    with open(arguments['out_dir_plots']+'/bed_output/' + output, 'a') as fp:
        fp.write(header_comments)
        fp.write(header)

def seperate_dfs_coverage(arguments, df, haplotype_1_values_updated, haplotype_2_values_updated, unphased):
    if arguments['without_phasing']:
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

def apply_copynumbers(csv_df_coverage, depth_values, depth_values1, arguments, snps_cpd_means, snps_cpd_means_collective):

    #depth_values = flatten(haplotype_1_values_updated)
    #depth_values1 = flatten(haplotype_2_values_updated)
    NULL_LOG2_COVERAGE = -20.0

    df_coverage = csv_df_coverage

    depth_values1 = np.clip(depth_values1, a_min=0, a_max=600)
    depth_values = np.clip(depth_values, a_min=0, a_max=600)

    csv_df_coverage.drop(csv_df_coverage[(csv_df_coverage.chr == "chrX") | (csv_df_coverage.chr == "chrY")].index, inplace=True)

    csv_df_coverage = csv_df_coverage.assign(depth=depth_values)
    csv_df_coverage1 = csv_df_coverage.assign(depth=depth_values1)
    #csv_df_coverage = csv_df_coverage[csv_df_coverage['chr'] == 'chr1']
    #csv_df_coverage = csv_df_coverage[csv_df_coverage['chr'] == 'chr7'] #TODO remove it


    csv_df_coverage = csv_df_coverage.assign(log2=0)
    #csv_df_coverage.loc[(csv_df_coverage['depth'] == 0), 'log2'] = NULL_LOG2_COVERAGE
    #ok_idx = csv_df_coverage["depth"] > 0
    #csv_df_coverage.loc[ok_idx, "log2"] = numpy.log2(csv_df_coverage.loc[ok_idx, "depth"])

    csv_df_coverage1 = csv_df_coverage1.assign(log2=0)
    #csv_df_coverage1.loc[(csv_df_coverage1['depth'] == 0), 'log2'] = NULL_LOG2_COVERAGE
    #ok_idx = csv_df_coverage1["depth"] > 0
    #csv_df_coverage1.loc[ok_idx, "log2"] = numpy.log2(csv_df_coverage1.loc[ok_idx, "depth"])

    csv_df_coverage = pd.concat([csv_df_coverage, csv_df_coverage1], ignore_index=True)
    #depth_values.extend(depth_values1)

    from cnvlib import cluster
    #cluster.clusters_tests(depth_values, depth_values1)

    #cluster.hmm_validation()

    # Fill in CNA required columns
    if "gene" in csv_df_coverage:
        csv_df_coverage["gene"] = csv_df_coverage["gene"].fillna("-")
    else:
        csv_df_coverage["gene"] = "-"
    csv_df_coverage.rename(columns={"chr": "chromosome"}, inplace=True)
    if arguments['without_phasing']:
        csv_df_coverage.drop(columns=['coverage'], inplace=True)
    else:
        csv_df_coverage.drop(columns=['hp1', 'hp2', 'hp3'], inplace=True)

    #fasta = '/home/rezkuh/GenData/reference/parts/chr7.fasta'
    #annot = '/home/rezkuh/gits/Wakhan/src/data/refflat.bed'
    #ref_fname, tgt_bed_fname, _ = batch.batch_make_reference([normal_bam], 'data/bins.bed', None, True, fasta, annot, True, 50000, None, None, None, None, "build", 8, False, "wgs", False, )
    #ref_fname = '/home/rezkuh/gits/Wakhan/src/data/colo829/reference.cnn'

    meta = {"sample_id": 'sample'}
    cnarr = CNA(csv_df_coverage, meta)
    log2_key = "log2"
    spread_key = "spread"
    #cnarr, ref_matched = fix.load_adjust_coverages(cnarr, read_cna(ref_fname), True, False, False, False)
    #cnarr.data["log2"] -= ref_matched[log2_key]
    cnarr = fix.apply_weights_replica(cnarr, None, log2_key, spread_key)
    cnarr.center_all(skip_low=False)

    cnarr['log2'] =cnarr['depth']
    #tabio.write(cnarr, "colo829_normal_grch38_md_chr7_haplotagged.cnr")

    #TODO variants = load_het_snps()
    #cnarr = read_cna('data/coverage_cnvkit.cnr')
    #PT8 cnarr.center_all(skip_low=True), skip_low=True, skip_outliers=20

    #segs, states, centers, stdev, cnarr, snps_cpd_means = segmentation.do_segmentation(depth_values, depth_values1, snps_cpd_means, arguments, cnarr, 'hmm', threshold=None, variants=None, skip_low=False, skip_outliers=0,
    #                                    min_weight=0, save_dataframe=False, rscript_path="Rscript", processes=1,
    #                                    smooth_cbs=False)

    segs, states, centers, stdev = call_copynumbers(arguments, cnarr, df_coverage, snps_cpd_means, snps_cpd_means_collective)

    #seg_metrics = segmetrics.do_segmetrics(cnarr, segs, interval_stats=["ci"], alpha=0.5, smoothed=True, skip_low=True,)
    #seg_call = call.do_call(seg_metrics, method="none", filters=["ci"])
    #seg_alltest = segmetrics.do_segmetrics(cnarr, seg_call, location_stats=["p_ttest"], skip_low=True)
    # Finally, assign absolute copy number values to each segment
    #seg_alltest.center_all("median")
    #seg_final = call.do_call(seg_alltest, method="threshold")

    #TODO _log2_ratio_to_absolute
    #tabio.write(segs, "colo829_normal_grch38_md_chr7_haplotagged.cns")

    values = cnarr.as_dataframe(cnarr.data)
    half_values = len(values) // 2

    return values.data.iloc[:half_values, ], segs[0], values.data.iloc[half_values:,], segs[1], states, centers, stdev

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

def detect_alter_loh_regions(arguments, event, chrom, ref_ends, haplotype_1_values, haplotype_2_values, unphased_reads_values, starts, ends, switch_hps):
    if ends and ends[-1] > ref_ends[-1]:
        ends[-1] = ref_ends[-1]

    region_starts = []
    region_ends = []
    #print(starts)
    #print(ends)
    for i, (start,end) in enumerate(zip(starts,ends)):
        if end - start > 1000000:
            region_starts.append(start)
            region_ends.append(end)

    #print(region_starts)
    #print(region_ends)

    if not arguments['without_phasing'] and switch_hps:
        for j, (starts,ends) in enumerate(zip(region_starts, region_ends)):
            #TODO Discuss with Ayse, alternate approach on what HP should be selected for each region
            #if mean_values(haplotype_1_values, starts - 1, starts - 4) > mean_values(haplotype_2_values, starts - 1, starts - 4):
            for i in range(starts//50000,ends//50000):
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
def loh_regions_events(chrom, region_starts, region_ends, arguments):
    dict = []
    for i in range(len(region_starts)):
        dict.append((chrom + '\t' + str(region_starts[i]) + '\t' + str(region_ends[i])))
    #write_segments_coverage(dict, arguments['genome_name'] + '_loh_segments.bed')
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

def detect_first_copy_integers_fractional_cluster_means(arguments, df_segs_hp1, df_segs_hp2, centers):
    haplotype_1_values_copy = df_segs_hp1.state.values.tolist()
    haplotype_1_start_values_copy = df_segs_hp1.start.values.tolist()
    haplotype_1_end_values_copy = df_segs_hp1.end.values.tolist()

    centers_count = len(centers)
    diff_lists = [[] for i in range(0, centers_count)]

    for i, (start, end, value) in enumerate(zip(haplotype_1_start_values_copy, haplotype_1_end_values_copy, haplotype_1_values_copy)):
        k = 0
        for i in range(len(centers)):
            if value == centers[i]:
                k = i
        diff_lists[k].append(end-start)

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
    integer_centers = []
    if arguments['without_phasing']:
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

def write_copynumber_segments_csv(haplotype_df, arguments, centers, integer_fractional_means, hp):
    fp = open(arguments['out_dir_plots']+'/bed_output/' + arguments['genome_name'] + '_copynumbers_segments.bed', 'a')

    for i in range(len(integer_fractional_means)):
        haplotype_df['depth'].mask(haplotype_df['depth'] == i, integer_fractional_means[i], inplace=True)

    haplotype_df = haplotype_df.rename(columns={'chromosome': 'chr', 'start': 'start', 'end': 'end', 'depth':'copynumber_state', 'state':'coverage'})
    if arguments['without_phasing']:
        fp.write('#chr: chromosome number\n')
        fp.write('#start: start address for CN segment\n')
        fp.write('#end: end address for CN segment\n')
        fp.write('#copynumber_state: detected copy number state (integer/fraction)\n')
        fp.write('#coverage: median coverage for this segment\n')

        header = ['chr', 'start', 'end', 'copynumber_state', 'coverage']
        haplotype_df.to_csv(fp, sep='\t', columns=header, index=False, mode='a', header=True)
    else:
        header = ['chr', 'start', 'end', 'copynumber_state', 'coverage', 'haplotype']
        haplotype_df['haplotype'] = hp
        header_enable = False
        if hp == 1:
            fp.write('#chr: chromosome number\n')
            fp.write('#start: start address for CN segment\n')
            fp.write('#end: end address for CN segment\n')
            fp.write('#copynumber_state: detected copy number state (integer/fraction)\n')
            fp.write('#coverage: median coverage for this segment\n')
            fp.write('#haplotype: haplotype number\n')

            header_enable = True
        haplotype_df.to_csv(fp, sep='\t', columns=header, index=False, mode='a', header=header_enable)

def integer_fractional_cluster_means(arguments, df_segs_hp1, df_segs_hp2, centers):
    integer_centers, fractional_centers = detect_first_copy_integers_fractional_cluster_means(arguments, df_segs_hp1, df_segs_hp2, centers)
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


def change_point_detection_means(arguments, df_chrom):
    df_means_chr = []
    if arguments['without_phasing']:
        means = df_chrom.coverage.values.tolist()
        snps_mean, snps_len, snps_pos = change_point_detection_algo(arguments['bin_size'], means)
        snps_pos_start = []
        snps_pos_end = []
        for i in range(len(snps_pos)-1):
            snps_pos_start.append(snps_pos[i] if snps_pos[i] < 1 else snps_pos[i]+1)
            snps_pos_end.append(snps_pos[i+1])
        chr_list = [df_chrom['chr'].iloc[0] for ch in range(len(snps_mean))]
        df_means_chr = pd.DataFrame(list(zip(chr_list, snps_pos_start, snps_pos_end, snps_mean)), columns=['chr', 'start', 'end', 'means'])

        return snps_mean, df_means_chr
    else:
        df_means_chr = []
        haplotype1_means = df_chrom.hp1.values.tolist()
        snps_haplotype1_mean, snps_haplotype1_len, snps_haplotype1_pos = change_point_detection_algo(arguments['bin_size'], haplotype1_means)
        haplotype2_means = df_chrom.hp2.values.tolist()
        snps_haplotype2_mean, snps_haplotype2_len, snps_haplotype2_pos = change_point_detection_algo(arguments['bin_size'], haplotype2_means)
        snps_pos_start = []
        snps_pos_end = []
        for i in range(len(snps_haplotype1_pos) - 1):
            snps_pos_start.append(snps_haplotype1_pos[i] if snps_haplotype1_pos[i] < 1 else snps_haplotype1_pos[i] + 1)
            snps_pos_end.append(snps_haplotype1_pos[i + 1])
        chr_list = [df_chrom['chr'].iloc[0] for ch in range(len(snps_haplotype1_mean))]

        df_means_chr.append(pd.DataFrame(list(zip(chr_list, snps_pos_start, snps_pos_end, snps_haplotype1_mean)),
                                    columns=['chr', 'start', 'end', 'means']))

        snps_pos_start = []
        snps_pos_end = []
        for i in range(len(snps_haplotype2_pos) - 1):
            snps_pos_start.append(snps_haplotype2_pos[i] if snps_haplotype2_pos[i] < 1 else snps_haplotype2_pos[i] + 1)
            snps_pos_end.append(snps_haplotype2_pos[i + 1])
        chr_list = [df_chrom['chr'].iloc[0] for ch in range(len(snps_haplotype2_mean))]

        df_means_chr.append(pd.DataFrame(list(zip(chr_list, snps_pos_start, snps_pos_end, snps_haplotype2_mean)),
                                         columns=['chr', 'start', 'end', 'means']))
        print(snps_haplotype1_mean)
        print(snps_haplotype2_mean)
        return snps_haplotype1_mean + snps_haplotype2_mean, df_means_chr

def change_point_detection_algo(bin_size, haplotype_means):
    data = np.array(haplotype_means, dtype='uint8')  # numpy.clip(haplotype1_means, a_min=1, a_max=300)
    algo = rpt.Pelt(model="rbf", jump=25).fit(data)
    result = algo.predict(pen=10)
    change_points = [i for i in result if i <= len(data)]
    strong_candidates = []
    snps_haplotype_mean = []
    snps_haplotype_pos = []
    snps_haplotype_len = []
    start = 0
    snps_haplotype_pos.append(0)
    for index, point in enumerate(change_points):
        sub_list = haplotype_means[start:point]

        snps_haplotype_mean.append(statistics.median(sub_list))
        snps_haplotype_len.append(len(sub_list))

        # if len(sub_list) < 25:
        #     continue
        # else:
        #     count = 0
        #     sub_list = [int(i) for i in sub_list]
        #     means_hp1 = statistics.median(sub_list)
        #     for i in range(len(sub_list)):
        #         if means_hp1 - 20 <= sub_list[i] <= means_hp1 + 20:
        #             count += 1
        #     if count > len(sub_list) // 3:
        #         snps_haplotype_mean.append(means_hp1)
        #         snps_haplotype_len.append(len(sub_list))
        #
        #     if count > len(sub_list) / 1.1:
        #         strong_candidates.append(int(means_hp1))
        #
        #     for i in range(len(sub_list)):
        #         if sub_list[i] >= means_hp1 + 10 or sub_list[i] <= means_hp1 - 10:
        #             sub_list[i] = randint(int(means_hp1) - 10, int(means_hp1) + 10)
        #     haplotype_means[start:point] = [sub_list[i] for i in range(len(sub_list))]

        start = point + 1
        snps_haplotype_pos.append(point * bin_size)

    return snps_haplotype_mean, snps_haplotype_len, snps_haplotype_pos

def adjust_diversified_segments(centers, snps_cpd_means_df, df_segs_hp1, df_segs_hp2, arguments):
    chroms = get_contigs_list(arguments['contigs'])
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

        if arguments['without_phasing']:
            snps_cpd_means_df_chrom = snps_cpd_means_df[snps_cpd_means_df['chr'] == chrom]
        else:
            snps_cpd_means_df_chrom = snps_cpd_means_df[0][snps_cpd_means_df[0]['chr'] == chrom]

        start_values_cpd = snps_cpd_means_df_chrom.start.values.tolist()
        end_values_cpd = snps_cpd_means_df_chrom.end.values.tolist()
        mean_cpd = snps_cpd_means_df_chrom.means.values.tolist()

        for i, (start,end) in enumerate(zip(haplotype_1_start_values_copyrnumbers, haplotype_1_end_values_copyrnumbers)):
            for j, (start_cpd, end_cpd) in enumerate(zip(start_values_cpd, end_values_cpd)):
                if start >= start_cpd and end <= end_cpd:
                    haplotype_1_values_copyrnumbers[i] = min(centers, key=lambda x:abs(x-mean_cpd[j]))

        updated_df_segs_hp1.append(pd.DataFrame(list(zip(df_chrom_segs_hp1.chromosome.values.tolist(), df_chrom_segs_hp1.start.values.tolist(), df_chrom_segs_hp1.end.values.tolist(),  df_chrom_segs_hp1.depth.values.tolist(), haplotype_1_values_copyrnumbers)),
                                    columns=['chromosome', 'start', 'end', 'depth', 'state']))

        if not arguments['without_phasing']:
            snps_cpd_means_df_chrom = snps_cpd_means_df[1][snps_cpd_means_df[1]['chr'] == chrom]
            start_values_cpd = snps_cpd_means_df_chrom.start.values.tolist()
            end_values_cpd = snps_cpd_means_df_chrom.end.values.tolist()
            mean_cpd = snps_cpd_means_df_chrom.means.values.tolist()

            for i, (start, end) in enumerate(
                    zip(haplotype_2_start_values_copyrnumbers, haplotype_2_end_values_copyrnumbers)):
                for j, (start_cpd, end_cpd) in enumerate(zip(start_values_cpd, end_values_cpd)):
                    if start >= start_cpd and end <= end_cpd:
                        haplotype_2_values_copyrnumbers[i] = min(centers, key=lambda x: abs(x - mean_cpd[j]))

            updated_df_segs_hp2.append(pd.DataFrame(list(zip(df_chrom_segs_hp2.chromosome.values.tolist(), df_chrom_segs_hp2.start.values.tolist(),
                    df_chrom_segs_hp2.end.values.tolist(), df_chrom_segs_hp2.depth.values.tolist(), haplotype_2_values_copyrnumbers)), columns=['chromosome', 'start', 'end', 'depth', 'state']))

    if arguments['without_phasing']:
        return merge_adjacent_regions_cn(pd.concat(updated_df_segs_hp1), arguments), merge_adjacent_regions_cn(pd.concat(updated_df_segs_hp1), arguments)
    else:
        return merge_adjacent_regions_cn(pd.concat(updated_df_segs_hp1), arguments), merge_adjacent_regions_cn(pd.concat(updated_df_segs_hp2), arguments)

def merge_adjacent_regions_cn(segarr, arguments):
    chroms = get_contigs_list(arguments['contigs'])
    dfs = []
    for index, chrom in enumerate(chroms):
        seg = segarr[segarr['chromosome'] == chrom]
        label_groups = seg['state'].ne(seg['state'].shift()).cumsum()
        df = (seg.groupby(label_groups).agg({'chromosome': 'first', 'start': 'min', 'end': 'max', 'depth': 'first', 'state': 'first'}).reset_index(drop=True))
        dfs.append(df)
    out = pd.concat(dfs)
    return out