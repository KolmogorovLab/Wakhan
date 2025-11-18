#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from src.__version__ import __version__
from collections import defaultdict
from datetime import datetime
import sys
import os
import pandas as pd
import ast

from src.utils_tmp.chromosome import get_contigs_list, csv_df_chromosomes_sorter
from src.coverage.binning import get_chromosomes_regions

class vcf_format(object):
    __slots__ = (
    'chrom', 'pos', 'ID', 'alt', 'sv_type', 'bps', 'sv_len', 'qual', 'Filter', 'cn_type', 'sample')

    def __init__(self, chrom, pos, ID, alt, sv_type, bps, sv_len, qual, Filter, cn_type, sample):
        self.chrom = chrom
        self.pos = pos
        self.ID = ID
        self.alt = alt
        self.sv_type = sv_type
        self.bps = bps
        self.sv_len = sv_len
        self.qual = qual
        self.Filter = Filter
        self.cn_type = cn_type
        self.sample = sample

    def svlen(self):
        if self.sv_len > 0:
            return f"SVLEN={self.sv_len}"
        else:
            return ""

    def bps_segments(self):
        if len(self.bps) > 0:
            list_bps = ','.join(ast.literal_eval(self.bps))
            return f";BPS={list_bps}"
        else:
            return ""
    def id_data(self):
        return f"wakhan:{self.cn_type}:{self.chrom}:{self.pos}-{self.sv_len+self.pos}" #wakhan:CNLOH:chr1:818023-16526823
    def info(self):
        return f"SVTYPE={self.sv_type};{self.svlen()};END={self.sv_len+self.pos}{self.bps_segments()}"


    def to_vcf(self):
        return f"{self.chrom}\t{self.pos}\t{self.id_data()}\tN\t{self.alt}\t{self.qual}\t{self.Filter}\t{self.info()}\tGT:TCN:CN1:CN2:CNQ1:CNQ2:COV1:COV2\t{self.sample}\n"

#GT:CN1:CN2:CN1Q:CN2Q
class vcf_sample(object):
    __slots__ = ('GT', 'TCN', 'CN1', 'CN2', 'CN1Q', 'CN2Q', 'CN1COV', 'CN2COV')

    def __init__(self, GT, TCN, CN1, CN2, CN1Q, CN2Q, CN1COV, CN2COV):
        self.GT = None
        self.TCN = round(CN1, 2) + round(CN2, 2)
        self.CN1 = round(CN1, 2)
        self.CN2 = round(CN2, 2)
        self.CN1Q = CN1Q
        self.CN2Q = CN2Q
        self.CN1COV = CN1COV
        self.CN2COV = CN2COV

    def call_genotype(self):
        # if self.CN1 > 1 and self.CN2 > 1:
        #     GT = '1/1'
        # elif (self.CN1 > 1 or self.CN2 > 1) and (not self.CN1 == self.CN2):
        #     GT = '0/1'
        # elif (self.CN1 == 0 and self.CN2 > 1) or (self.CN1 > 1 and self.CN2 == 0):
        #     GT = '1/2'
        # else:
        #     GT = '0/0'
        if (round(self.CN1) == 1 and not round(self.CN2) == 1) or (round(self.CN2) == 1 and not round(self.CN1) == 1):
            GT = '0/1'
        elif (round(self.CN1) > 1 and round(self.CN2) > 1) or (round(self.CN1) == 0 and round(self.CN2) == 0):
            GT = '1/1'
        elif round(self.CN1) == 1 and round(self.CN2) == 1:
            GT = '0/0'
        elif (round(self.CN1) == 0 and round(self.CN2) > 1) or (round(self.CN2) == 0 and round(self.CN1) > 1):
            GT = '1/2'
        #print(round(self.CN1), round(self.CN2))
        self.GT = GT

    def sample(self):
        self.call_genotype()
        return f"{self.GT}:{self.TCN}:{self.CN1}:{self.CN2}:{self.CN1Q}:{self.CN2Q}:{self.CN1COV}:{self.CN2COV}"

def db_2_vcf(df):
    vcf_list = []
    for index, seg in df.iterrows():
        if round(seg['state']) == 1 and round(seg['state_2']) == 1:
            continue
        chrom = seg['chr']
        pos= seg['start']
        ID = ''
        sv_type = 'CNV'
        sv_len = seg['end'] - seg['start']
        qual = 1000 #TODO
        Filter = 'PASS'

        if (round(seg['state']) == 0 and round(seg['state_2']) > 1) or (round(seg['state_2']) == 0 and round(seg['state']) > 1):
            cn_type = 'CNLOH'
        elif (round(seg['state']) > 1 and round(seg['state_2']) > 0) or (round(seg['state']) > 0 and round(seg['state_2']) > 1):
            cn_type = 'GAIN'
        elif round(seg['state']) == 1 and round(seg['state_2']) == 1:
            cn_type = 'REF'
        else:
            cn_type = 'LOSS'

        if cn_type == 'LOSS':
            alt_type = '<DEL>'
        elif cn_type == 'GAIN':
            alt_type = '<DUP>'
        elif cn_type == 'REF':
            alt_type = '.'
        else:
            if round(seg['state']) == 0 and round(seg['state_2']) > 1:
                alt_type = '<DEL>,<DUP>'
            else:
                alt_type = '<DUP>,<DEL>'

        alt = alt_type
        bps = seg['bps']
        sample = vcf_sample('', '', seg['state'], seg['state_2'], seg['p_value'], seg['p_value_2'], seg['depth'], seg['depth_2']).sample()
        vcf_list.append(
            vcf_format(chrom, pos, ID, alt, sv_type, bps, sv_len, qual, Filter, cn_type, sample))
    return vcf_list

def write_vcf_header(ref_lengths, outfile, sample_list, type):
    sample = '\t'.join(sample_list)
    outfile.write("##fileformat=VCFv4.2\n")
    outfile.write('##source=Wakhan_v' + __version__ + '\n')
    outfile.write('##CommandLine= ' + " ".join(sys.argv[1:]) + '\n')
    filedate = str(datetime.now()).split(' ')[0]
    outfile.write('##fileDate=' + filedate + '\n')  #
    for chr_id, chr_len in ref_lengths.items():
        outfile.write("##contig=<ID={0},length={1}>\n".format(chr_id, chr_len))  #
    outfile.write('##ALT=<ID=CNV,Description="Copy number variant region">\n')
    outfile.write('##ALT=<ID=DEL,Description="Deletion relative to the reference">\n')
    outfile.write('##ALT=<ID=DUP,Description="Region of elevated copy number relative to the reference">\n')

    outfile.write('##INFO=<ID=REFLEN,Number=1,Type=Integer,Description="Number of REF positions included in this record">\n')
    outfile.write('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
    outfile.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    outfile.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
    outfile.write('##INFO=<ID=HET,Number=0,Type=Flag,Description="Segment is heterogeneous">\n')
    outfile.write('##INFO=<ID=BPS,Number=0,Type=String,Description="Breakpoints covering segment">\n')

    outfile.write('##FILTER=<ID=PASS,Description="All filters passed">\n')

    #GT:CN1:CN2:CN1Q:CN2Q:COV1:COV2
    outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    outfile.write('##FORMAT=<ID=TCN,Number=1,Type=Float,Description="Estimated total copy numbers of segment">\n')

    if type == 'integers':
        outfile.write('##FORMAT=<ID=CN1,Number=1,Type=Float,Description="Estimated haplotype-1 segment copy number">\n')
        outfile.write('##FORMAT=<ID=CN2,Number=1,Type=Float,Description="Estimated haplotype-2 segment copy number">\n')
    else:
        outfile.write('##FORMAT=<ID=CN1,Number=1,Type=Float,Description="Estimated subclonal haplotype-1 segment copy number">\n')
        outfile.write('##FORMAT=<ID=CN2,Number=1,Type=Float,Description="Estimated subclonal haplotype-2 segment copy number">\n')

    outfile.write('##FORMAT=<ID=CNQ1,Number=1,Type=Float,Description="Estimated haplotype-1 segment confidence score">\n')
    outfile.write('##FORMAT=<ID=CNQ2,Number=1,Type=Float,Description="Estimated haplotype-2 segment confidence score">\n')
    outfile.write('##FORMAT=<ID=COV1,Number=1,Type=Float,Description="Estimated haplotype-1 segment coverage value">\n')
    outfile.write('##FORMAT=<ID=COV2,Number=1,Type=Float,Description="Estimated haplotype-2 segment coverage value">\n')

    outfile.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n")

def _sorted_breakpoints(bp_list, ref_lengths):
    db_ls = defaultdict(list)
    dbls = []
    for chr_id in ref_lengths.keys():
        db_ls[chr_id] = []

    for db in bp_list:
        db_ls[db.chrom].append(db)

    for chr_id, dbs in db_ls.items():
        dbls += sorted(dbs, key=lambda x: x.pos)

    return dbls


def write_cna_vcf(vcf_list, outfile, ref_lengths):
    for db in _sorted_breakpoints(vcf_list, ref_lengths):
        outfile.write(db.to_vcf())
    outfile.close()


def write_to_vcf(args, df, vcf_path, type):
    vcf_list = db_2_vcf(df)
    sample_id = 'Sample' #TODO from BAM
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)
    ref_lengths = dict(zip(chroms, regions))
    cna_outfile = open(vcf_path, "w")
    write_vcf_header(ref_lengths, cna_outfile, sample_id, type)
    write_cna_vcf(vcf_list, cna_outfile, ref_lengths)

def split_segments_by_breakpoints(df, breakpoints_list, chr_col='chr', start_col='start', end_col='end', depth_col='depth', state_col='state', p_value_col='p_value', bps_col='bps'):
    """
    Splits segments in a DataFrame based on a list of unique breakpoint positions.
    New rows are created for the resulting sub-segments, inheriting the original state.

    Args:
        df (pd.DataFrame): Input DataFrame with genomic segments.
                           Expected columns: 'chr', 'start', 'end', 'state'.
        breakpoints_list (list): A list of unique integer/float breakpoint positions.
        chr_col (str): Name of the chromosome column.
        start_col (str): Name of the segment start column.
        end_col (str): Name of the segment end column.
        state_col (str): Name of the segment state column.

    Returns:
        pd.DataFrame: A new DataFrame with segments split by breakpoints.
                      The DataFrame will be sorted by chromosome and start position.
    """
    if not all(col in df.columns for col in [chr_col, start_col, end_col, depth_col, state_col, p_value_col, bps_col]):
        raise ValueError(f"Input DataFrame must contain columns: {chr_col}, {start_col}, {end_col}, {depth_col}, {state_col}, {p_value_col}, {bps_col}")

    # Ensure breakpoints are sorted and unique for efficient processing
    unique_breakpoints = sorted(list(set(breakpoints_list)))

    new_segments_data = []

    # Iterate over each original segment in the DataFrame
    for index, row in df.iterrows():
        current_chr = row[chr_col]
        current_start = row[start_col]
        current_end = row[end_col]
        current_depth = row[depth_col]
        current_state = row[state_col]
        current_p_value = row[p_value_col]
        current_bps = row[bps_col]

        # Identify relevant breakpoints within the current segment
        # We include current_start and current_end to define potential split points
        segment_split_points = [current_start]
        for bp in unique_breakpoints:
            if current_start < bp < current_end:
                segment_split_points.append(bp)
        segment_split_points.append(current_end)

        # Sort and remove duplicates from split points
        # This handles cases where breakpoints might coincide with original start/end
        segment_split_points = sorted(list(set(segment_split_points)))

        # Create new sub-segments
        for i in range(len(segment_split_points) - 1):
            new_segment_start = segment_split_points[i]
            new_segment_end = segment_split_points[i+1]

            # Only add if it's a valid, non-zero-length segment
            if new_segment_start < new_segment_end:
                new_segments_data.append({
                    chr_col: current_chr,
                    start_col: new_segment_start,
                    end_col: new_segment_end,
                    depth_col: current_depth,
                    state_col: current_state,
                    p_value_col: current_p_value,
                    bps_col: current_bps
                })

    # Create a new DataFrame from the collected segments
    new_df = pd.DataFrame(new_segments_data, columns=df.columns)

    # Sort the final DataFrame for a clean output
    new_df = new_df.sort_values(by=[chr_col, start_col]).reset_index(drop=True)

    return new_df

def extract_first_from_consecutive_groups(lst):
    if not lst:
        return []

    result = [lst[0]]  # Always include the first element
    for i in range(1, len(lst)):
        if lst[i] != lst[i - 1] + 1:
            result.append(lst[i])
    return result

def read_cn_segments_process_vcf(args, repo, type):
    if type == 'integers':
        fp_hp1 = args.out_dir_plots + '/' + repo + '/bed_output/' + args.genome_name+'_'+ repo + '_copynumbers_segments_HP_1.bed'
        fp_hp2 = args.out_dir_plots + '/' + repo + '/bed_output/' + args.genome_name+'_'+ repo + '_copynumbers_segments_HP_2.bed'
        if not os.path.exists(args.out_dir_plots + '/' + repo + '/vcf_output/'):
            os.mkdir(args.out_dir_plots + '/' + repo + '/vcf_output/')
        output_path = args.out_dir_plots + '/' + repo + '/vcf_output/' + args.genome_name+'_'+ repo + '_wakhan_cna_integers.vcf'
    else:
        fp_hp1 = args.out_dir_plots + '/' + repo + '/bed_output/' + args.genome_name + '_' + repo + '_copynumbers_subclonal_segments_HP_1.bed'
        fp_hp2 = args.out_dir_plots + '/' + repo + '/bed_output/' + args.genome_name + '_' + repo + '_copynumbers_subclonal_segments_HP_2.bed'
        if not os.path.exists(args.out_dir_plots + '/' + repo + '/vcf_output/'):
            os.mkdir(args.out_dir_plots + '/' + repo + '/vcf_output/')
        output_path = args.out_dir_plots + '/' + repo + '/vcf_output/' + args.genome_name + '_' + repo + '_wakhan_cna_subclonals.vcf'

    hp1_segs = pd.read_csv(fp_hp1, sep='\t', header=None, comment='#')
    if args.change_point_detection_for_cna:
        hp1_segs['bps_new'] = None
    if type == 'integers':
        hp1_segs = hp1_segs[[0, 1, 2, 3, 4, 5, 6]]
    else:
        hp1_segs = hp1_segs[[0, 1, 2, 3, 4, 5, 7]]
    hp1_segs.columns = ['chr', 'start', 'end', 'depth', 'state', 'p_value', 'bps']

    hp2_segs = pd.read_csv(fp_hp2, sep='\t', header=None, comment='#')
    if args.change_point_detection_for_cna:
        hp2_segs['bps_new'] = None
    if type == 'integers':
        hp2_segs = hp2_segs[[0, 1, 2, 3, 4, 5, 6]]
    else:
        hp2_segs = hp2_segs[[0, 1, 2, 3, 4, 5, 7]]
    hp2_segs.columns = ['chr', 'start', 'end', 'depth', 'state', 'p_value', 'bps']

    updated_hp1_segs = []
    updated_hp2_segs = []

    fileDir = os.path.dirname(__file__)  # os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)

    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        hp1_segs_chrom = hp1_segs[hp1_segs['chr'] == chrom]
        hp2_segs_chrom = hp2_segs[hp2_segs['chr'] == chrom]
        df_centm_chrom = df_centm[df_centm['chr'] == chrom]
        if not df_centm_chrom.empty:
            cents = [df_centm_chrom.start.values.tolist()[0], df_centm_chrom.end.values.tolist()[0]]
        else:
            cents = [0, 0]

        #hp1_segs_chrom = hp1_segs_chrom.drop(hp1_segs_chrom[hp1_segs_chrom.end == cents[1]].index)
        #hp2_segs_chrom = hp2_segs_chrom.drop(hp2_segs_chrom[hp2_segs_chrom.end == cents[1]].index)

        unique_points = sorted(list(set(hp1_segs_chrom.start.values.tolist()+hp2_segs_chrom.start.values.tolist())))
        unique_points = extract_first_from_consecutive_groups(unique_points)
        hp1_segs_chrom = split_segments_by_breakpoints(hp1_segs_chrom, unique_points)
        hp2_segs_chrom = split_segments_by_breakpoints(hp2_segs_chrom, unique_points)
        hp1_segs_chrom = hp1_segs_chrom.drop_duplicates()
        hp2_segs_chrom = hp2_segs_chrom.drop_duplicates()

        updated_hp1_segs.append(hp1_segs_chrom)
        updated_hp2_segs.append(hp2_segs_chrom)

    segs_hp1 = pd.concat(updated_hp1_segs)
    segs_hp2 = pd.concat(updated_hp2_segs)
    #merge both HPs on states
    if len(segs_hp1) != len(segs_hp2):
        print("DataFrames must have the same length.")
        return None

    segs_hp1['state_2'] = segs_hp2['state'].values
    segs_hp1['p_value_2'] = segs_hp2['p_value'].values
    segs_hp1['depth_2'] = segs_hp2['depth'].values

    mask_to_keep = ~((segs_hp1['state'] == 1) & (segs_hp1['state_2'] == 1))
    filtered_df = segs_hp1[mask_to_keep].reset_index(drop=True)

    write_to_vcf(args, filtered_df, output_path, type)

    return segs_hp1
