import csv
import numpy
import pandas as pd
import pysam
import os

from cnvlib.cmdutil import read_cna
from cnvlib.cnary import CopyNumArray as CNA
from cnvlib import segmentation, coverage, batch, fix, segmetrics, call, scatter
from skgenome import tabio
from smoothing import smoothing

from pomegranate import HiddenMarkovModel as Model

from phasing_correction import get_phasesets_bins

def generate_phasesets_bins(bam, path, bin_size, arguments):
    return get_phasesets_bins(bam, path, bin_size, arguments)

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

def get_contigs_list(contigs):
    chroms_list_final = []
    chroms = contigs.split(',')
    for chrom in chroms:
        chrom = chrom[len('chr'):] if chrom.startswith('chr') else chrom
        chrom = chrom.split('-')
        if len(chrom) > 1:
            chroms_list_final.extend(list(range(int(chrom[0]), int(chrom[1]) + 1)))
        else:
            chroms_list_final.extend(chrom)

    chroms_list_final = ['chr' + x if not x.startswith('chr') else x for x in map(str, chroms_list_final)]
    return chroms_list_final

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

def csv_df_chromosomes_sorter(path):
    dataframe = pd.read_csv(path, sep='\t', names=['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
    if not dataframe['chr'].iloc[0].startswith('chr'):
        dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', 'start'], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

def csv_df_chromosomes_sorter_copyratios(path):
    dataframe = pd.read_csv(path, sep='\t', names=['chr', 'start', 'end', 'gene', 'log2', 'depth', 'probes', 'weight'])
    if not dataframe['chr'].iloc[0].startswith('chr'):
        dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', 'start'], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

def csv_df_chromosomes_sorter_snps(path):
    dataframe = pd.read_csv(path, sep='\t', names=['chr', 'pos', 'qual', 'filter', 'ps', 'gt', 'dp', 'vaf'])
    if not dataframe['chr'].iloc[0].startswith('chr'):
        dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', 'pos'], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

def csv_df_chromosomes_sorter_snps_from_bam(path):
    dataframe = pd.read_csv(path, sep='\t', names=['chr', 'pos', 'freq_value_a', 'hp_a', 'freq_value_b', 'hp_b'])
    if not dataframe['chr'].iloc[0].startswith('chr'):
        dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', 'pos'], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

def csv_df_chromosomes_sorter_snps_frequency(path):
    dataframe = pd.read_csv(path, sep=',', names=['chr', 'start', 'a', 'c', 'g', 't'])
    if not dataframe['chr'].iloc[0].startswith('chr'):
        dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', 'start'], ascending=[True, True], inplace=True)
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

def write_segments_coverage(coverage_segments, output):
    with open('data/' + output, 'w') as fp:
        for items in coverage_segments:
            if not items == None:
                fp.write("%s\n" % items)

def seperate_dfs_coverage(df, haplotype_1_values_updated, haplotype_2_values_updated, unphased):

    df_hp1 = df[['chr', 'start','end', 'hp1']].copy()
    df_hp2 = df[['chr', 'start','end', 'hp2']].copy()
    df_unphased = df[['chr', 'start','end', 'hp3']].copy()
    df_hp1['hp1'] = haplotype_1_values_updated
    df_hp2['hp2'] = haplotype_2_values_updated
    df_unphased['hp3'] = unphased
    return df_hp1, df_hp2, df_unphased

def flatten(values):
    return [item for sublist in values for item in sublist]

def apply_copynumbers(csv_df_coverage, depth_values, depth_values1, arguments):

    #depth_values = flatten(haplotype_1_values_updated)
    #depth_values1 = flatten(haplotype_2_values_updated)
    NULL_LOG2_COVERAGE = -20.0

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
    depth_values.extend(depth_values1)

    from cnvlib import cluster
    #cluster.clusters_tests(depth_values, depth_values1)

    #cluster.hmm_validation()

    # Fill in CNA required columns
    if "gene" in csv_df_coverage:
        csv_df_coverage["gene"] = csv_df_coverage["gene"].fillna("-")
    else:
        csv_df_coverage["gene"] = "-"
    csv_df_coverage.rename(columns={"chr": "chromosome"}, inplace=True)
    csv_df_coverage.drop(columns=['hp1', 'hp2', 'hp3'], inplace=True)

    #fasta = '/home/rezkuh/GenData/reference/parts/chr7.fasta'
    #annot = '/home/rezkuh/gits/Wakhan/src/data/refflat.bed'
    #ref_fname, tgt_bed_fname, _ = batch.batch_make_reference([normal_bam], 'data/bins.bed', None, True, fasta, annot, True, 50000, None, None, None, None, "build", 8, False, "wgs", False, )
    ref_fname = '/home/rezkuh/gits/Wakhan/src/data/colo829/reference.cnn'

    meta = {"sample_id": 'sample'}
    cnarr = CNA(csv_df_coverage, meta)
    log2_key = "log2"
    spread_key = "spread"
    #cnarr, ref_matched = fix.load_adjust_coverages(cnarr, read_cna(ref_fname), True, False, False, False)
    #cnarr.data["log2"] -= ref_matched[log2_key]
    cnarr = fix.apply_weights_replica(cnarr, None, log2_key, spread_key)
    cnarr.center_all(skip_low=True)

    cnarr['log2'] =cnarr['depth']
    #tabio.write(cnarr, "colo829_normal_grch38_md_chr7_haplotagged.cnr")

    #TODO variants = load_het_snps()
    #cnarr = read_cna('data/coverage_cnvkit.cnr')
    #PT8 cnarr.center_all(skip_low=True), skip_low=True, skip_outliers=20
    segs = segmentation.do_segmentation(depth_values, arguments, cnarr, 'hmm', threshold=None, variants=None, skip_low=True, skip_outliers=20,
                                        min_weight=0, save_dataframe=False, rscript_path="Rscript", processes=1,
                                        smooth_cbs=False)



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

    return values.data.iloc[:half_values, ], segs[0].as_dataframe(segs[0].data), values.data.iloc[half_values:,], segs[1].as_dataframe(segs[1].data)

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
