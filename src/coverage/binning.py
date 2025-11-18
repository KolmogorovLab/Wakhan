import numpy as np
import pandas as pd
import pysam
import os

from src.utils.chromosome import get_contigs_list
from src.breakpoint.breakpoints import sv_vcf_bps_cn_check


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
