#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from src.__version__ import __version__
from collections import defaultdict
from datetime import datetime
import sys
import os
import pandas as pd

from src.utils.chromosome import get_contigs_list
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
            list_bps = self.bps.replace(';', ',')
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

def read_cn_segments_process_vcf(args, repo, type):
    """Reads the already boundary-aligned integer_profile.bed / subclonal_profile.bed
    (written by write_copynumber_segments_csv) and converts it to a CNA VCF. Since
    HP1/HP2 are already merged onto shared segment boundaries in that file, no
    re-splitting is needed here."""
    is_subclonal = type != 'integers'
    bed_filename = 'subclonal_profile.bed' if is_subclonal else 'integer_profile.bed'
    output_filename = 'subclonal_profile.vcf' if is_subclonal else 'integer_profile.vcf'

    columns = ['chr', 'start', 'end', 'hp1_coverage', 'hp1_copynumber_state', 'hp1_confidence',
               'hp2_coverage', 'hp2_copynumber_state', 'hp2_confidence']
    if is_subclonal:
        columns += ['hp1_is_subclonal', 'hp2_is_subclonal']
    if args.breakpoints:
        columns += ['svs_breakpoints_ids']

    segs = pd.read_csv(args.out_dir_plots + '/' + repo + '/' + bed_filename, sep='\t', header=None, comment='#', names=columns, na_filter=False)
    segs = segs.rename(columns={'hp1_coverage': 'depth', 'hp1_copynumber_state': 'state', 'hp1_confidence': 'p_value',
                                'hp2_coverage': 'depth_2', 'hp2_copynumber_state': 'state_2', 'hp2_confidence': 'p_value_2'})
    segs['bps'] = segs['svs_breakpoints_ids'] if 'svs_breakpoints_ids' in segs.columns else ''

    mask_to_keep = ~((segs['state'] == 1) & (segs['state_2'] == 1))
    filtered_df = segs[mask_to_keep].reset_index(drop=True)

    output_path = args.out_dir_plots + '/' + repo + '/' + output_filename
    write_to_vcf(args, filtered_df, output_path, type)

    return segs
