import pandas as pd
import pysam
import os
from phasing_correction import get_phasesets_bins
def generate_phasesets_bins(bam, path, bin_size):
    return get_phasesets_bins(bam, path, bin_size)

def get_chromosomes_bins(bam_file, bin_size):
    bed=[]
    bam_alignment = pysam.AlignmentFile(bam_file)
    headers = bam_alignment.header
    seq_dict = headers['SQ']
    region = [''] * len(seq_dict)
    chrs = [''] * len(seq_dict)
    head, tail = os.path.split(bam_file)
    for i, seq_elem in enumerate(seq_dict):
        region[i] = seq_elem['LN']
        chrs[i] = seq_elem['SN']
        start=0
        end=bin_size
        if chrs[i]:
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

def csv_df_chromosomes_sorter(path):
    dataframe = pd.read_csv(path, sep='\t', names=['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
    dataframe.sort_values(by=['chr', 'start'], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

def csv_df_chromosomes_sorter_snps(path):
    dataframe = pd.read_csv(path, sep='\t', names=['chr', 'pos', 'qual', 'filter', 'ps', 'gt', 'dp', 'vaf'])
    dataframe.sort_values(by=['chr', 'pos'], ascending=[True, True], inplace=True)
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
            fp.write("%s\n" % items)