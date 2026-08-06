import os
import pandas as pd
from intervaltree import IntervalTree, Interval

from src.cna.copynumber import subclonal_values_adjusted, integers_values_adjusted, add_confidence_score_cn_segemnts
from src.breakpoint.breakpoints import sv_vcf_bps_cn_check


def write_df_csv(df, file_name):
    df.to_csv(file_name, sep='\t', index=False, header=False, mode='w')


def write_header_comments(header, header_comments, output, args):
    with open(args.out_dir_plots + '/bed_output/' + output, 'a') as fp:
        fp.write(header_comments)
        fp.write(header)


def _compute_haplotype_segments_df(df_hp1, df_hp2, haplotype_df_, args, centers, integer_fractional_means, hp, is_subclonal, bps_ids_all):
    haplotype_df = haplotype_df_.copy()

    if not is_subclonal:
        for i in range(len(integer_fractional_means)):
            haplotype_df['state'] = haplotype_df['state'].mask(haplotype_df['state'] == centers[i], integer_fractional_means[i])
        # Convert any high-CN states not matched by the mask (extrapolated beyond known centers)
        if len(centers) >= 2:
            haplotype_df['state'] = haplotype_df['state'].apply(
                lambda x: integers_values_adjusted(x, centers) if x > centers[-1] else x
            )

    if 'p_value' in haplotype_df.columns:
        haplotype_df['is_subclonal'] = ['Y' if element < args.confidence_subclonal_score and not state == 0 else 'N'
                                        for (element, state) in zip(haplotype_df.p_value.values.tolist(), haplotype_df.state.values.tolist())]
        haplotype_df['state'] = [subclonal_values_adjusted(element, centers) if sub == 'Y' else integers_values_adjusted(element, centers)
                                 for (element, sub) in zip(haplotype_df.state.values.tolist(), haplotype_df.is_subclonal.values.tolist())]
        haplotype_df = haplotype_df.rename(columns={'chromosome': 'chr', 'depth': 'coverage', 'state': 'copynumber_state', 'p_value': 'confidence'})
    else:
        if hp == 1:
            haplotype_df, _ = add_confidence_score_cn_segemnts(centers, haplotype_df, haplotype_df, df_hp1, df_hp2, args)
        else:
            _, haplotype_df = add_confidence_score_cn_segemnts(centers, haplotype_df, haplotype_df, df_hp1, df_hp2, args)
        haplotype_df = haplotype_df.rename(columns={'chromosome': 'chr', 'depth': 'coverage', 'state': 'copynumber_state', 'confidence_value': 'confidence'})

    haplotype_df['coverage'] = haplotype_df['coverage'].apply(lambda x: round(x, 2))

    if bps_ids_all is not None:
        bps_ids_global = []
        for index, row in haplotype_df.iterrows():
            bps_ids = set()
            for bp in bps_ids_all:
                if bp[0] == row['chr'] and int(bp[1]) >= int(row['start']) and int(bp[1]) <= int(row['end']):
                    bps_ids.add(bp[2])
            bps_ids_global.append(bps_ids)
        haplotype_df['svs_breakpoints_ids'] = bps_ids_global

    return haplotype_df


def merge_haplotype_segments(df_hp1_segs, df_hp2_segs, is_subclonal, has_breakpoints):
    """Boundary-align two independently-segmented haplotype tables into a single
    wide table per genomic fragment. HP1 and HP2 are segmented independently
    (separate change-point detection), so a segment boundary present in only
    one haplotype gets fragmented into the other; both resulting fragments
    keep the parent segment's original coverage/state/confidence. A fragment
    covered by only one haplotype (should not normally happen - centromere/LOH
    boundaries are synchronized between haplotypes upstream) is zero-filled
    for the missing haplotype rather than dropped.
    """
    merged_rows = []
    chroms = sorted(set(df_hp1_segs['chr']).union(df_hp2_segs['chr']))

    for chrom in chroms:
        h1 = df_hp1_segs[df_hp1_segs['chr'] == chrom]
        h2 = df_hp2_segs[df_hp2_segs['chr'] == chrom]

        tree1 = IntervalTree(Interval(row.start, row.end, row) for row in h1.itertuples() if row.end > row.start)
        tree2 = IntervalTree(Interval(row.start, row.end, row) for row in h2.itertuples() if row.end > row.start)

        boundaries = sorted(set(h1['start']).union(h1['end']).union(h2['start']).union(h2['end']))
        for frag_start, frag_end in zip(boundaries[:-1], boundaries[1:]):
            if frag_start >= frag_end:
                continue
            ivs1 = tree1[frag_start]
            ivs2 = tree2[frag_start]
            if not ivs1 and not ivs2:
                continue

            row = {'chr': chrom, 'start': frag_start, 'end': frag_end}
            ids1, ids2 = set(), set()

            if ivs1:
                r1 = next(iter(ivs1)).data
                row['hp1_coverage'] = r1.coverage
                row['hp1_copynumber_state'] = r1.copynumber_state
                row['hp1_confidence'] = r1.confidence
                if is_subclonal:
                    row['hp1_is_subclonal'] = r1.is_subclonal
                if has_breakpoints:
                    ids1 = r1.svs_breakpoints_ids
            else:
                row['hp1_coverage'] = 0
                row['hp1_copynumber_state'] = 0
                row['hp1_confidence'] = 0
                if is_subclonal:
                    row['hp1_is_subclonal'] = 'N'

            if ivs2:
                r2 = next(iter(ivs2)).data
                row['hp2_coverage'] = r2.coverage
                row['hp2_copynumber_state'] = r2.copynumber_state
                row['hp2_confidence'] = r2.confidence
                if is_subclonal:
                    row['hp2_is_subclonal'] = r2.is_subclonal
                if has_breakpoints:
                    ids2 = r2.svs_breakpoints_ids
            else:
                row['hp2_coverage'] = 0
                row['hp2_copynumber_state'] = 0
                row['hp2_confidence'] = 0
                if is_subclonal:
                    row['hp2_is_subclonal'] = 'N'

            if has_breakpoints:
                row['svs_breakpoints_ids'] = ';'.join(sorted(ids1 | ids2))

            merged_rows.append(row)

    columns = ['chr', 'start', 'end', 'hp1_coverage', 'hp1_copynumber_state', 'hp1_confidence',
               'hp2_coverage', 'hp2_copynumber_state', 'hp2_confidence']
    if is_subclonal:
        columns += ['hp1_is_subclonal', 'hp2_is_subclonal']
    if has_breakpoints:
        columns += ['svs_breakpoints_ids']

    return pd.DataFrame(merged_rows, columns=columns)


def _write_segments_header(fp, columns, is_subclonal, has_breakpoints):
    fp.write('#chr: chromosome number\n')
    fp.write('#start: start address for CN segment\n')
    fp.write('#end: end address for CN segment\n')
    fp.write('#hp1_coverage / hp2_coverage: median coverage for this segment, per haplotype\n')
    fp.write('#hp1_copynumber_state / hp2_copynumber_state: detected copy number state (integer/fraction), per haplotype\n')
    fp.write('#hp1_confidence / hp2_confidence: confidence score, per haplotype\n')
    if is_subclonal:
        fp.write('#hp1_is_subclonal / hp2_is_subclonal: if entry is subclonal [Y/N], per haplotype\n')
    if has_breakpoints:
        fp.write('#svs_breakpoints_ids: corresponding structural variations (breakpoints) IDs from VCF file (union across haplotypes)\n')
    fp.write('#' + '\t'.join(columns) + '\n')


def write_copynumber_segments_csv(df_hp1, df_hp2, df_segs_hp1_updated, df_segs_hp2_updated, args, repo,
                                  centers, integer_fractional_means, filename, is_subclonal):
    """Compute per-haplotype CN segments, merge HP1/HP2 onto a shared set of
    boundaries, and write the result to <out_dir_plots>/<repo>/<filename>."""
    bps_ids_all = None
    if args.breakpoints:
        _, _, _, bps_ids_all, bps, bps_bnd = sv_vcf_bps_cn_check(args.breakpoints, args)

    hp1_df = _compute_haplotype_segments_df(df_hp1, df_hp2, df_segs_hp1_updated, args, centers, integer_fractional_means, 1, is_subclonal, bps_ids_all)
    hp2_df = _compute_haplotype_segments_df(df_hp1, df_hp2, df_segs_hp2_updated, args, centers, integer_fractional_means, 2, is_subclonal, bps_ids_all)

    merged_df = merge_haplotype_segments(hp1_df, hp2_df, is_subclonal, args.breakpoints is not None)

    out_dir = args.out_dir_plots + '/' + repo
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    with open(out_dir + '/' + filename, 'w') as fp:
        _write_segments_header(fp, merged_df.columns, is_subclonal, args.breakpoints is not None)
        merged_df.to_csv(fp, sep='\t', index=False, header=False)


def write_unphased_copynumber_segments_csv(df, df_segs_updated, args, centers, integer_fractional_means, filename, is_subclonal):
    """Single-track (--without-phasing) equivalent of write_copynumber_segments_csv:
    there is no HP1/HP2 split to merge, so this just writes the one segment table."""
    bps_ids_all = None
    if args.breakpoints:
        _, _, _, bps_ids_all, bps, bps_bnd = sv_vcf_bps_cn_check(args.breakpoints, args)

    haplotype_df = _compute_haplotype_segments_df(df, df, df_segs_updated, args, centers, integer_fractional_means, 2, is_subclonal, bps_ids_all)

    columns = ['chr', 'start', 'end', 'coverage', 'copynumber_state', 'confidence']
    if is_subclonal:
        columns.append('is_subclonal')
    if bps_ids_all is not None:
        haplotype_df['svs_breakpoints_ids'] = haplotype_df['svs_breakpoints_ids'].apply(lambda ids: ';'.join(sorted(ids)))
        columns.append('svs_breakpoints_ids')

    out_dir = args.out_dir_plots + '/solution'
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    with open(out_dir + '/' + filename, 'w') as fp:
        fp.write('#chr: chromosome number\n')
        fp.write('#start: start address for CN segment\n')
        fp.write('#end: end address for CN segment\n')
        fp.write('#coverage: median coverage for this segment\n')
        fp.write('#copynumber_state: detected copy number state (integer/fraction)\n')
        fp.write('#confidence: confidence score\n')
        if is_subclonal:
            fp.write('#is_subclonal: if entry is subclonal [Y/N]\n')
        if bps_ids_all is not None:
            fp.write('#svs_breakpoints_ids: corresponding structural variations (breakpoints) IDs from VCF file\n')
        fp.write('#' + '\t'.join(columns) + '\n')
        haplotype_df[columns].to_csv(fp, sep='\t', index=False, mode='a', header=False)


def write_segments_coverage_dict(coverage_segments, output, args):
    with open(args.out_dir_plots+'/coverage_data/' + output, 'a') as fp:
        for items in coverage_segments:
            if not items == None:
                fp.write("%s\n" % items)


def write_segments_coverage_snps(coverage_segments, output, args):
    with open(args.out_dir_plots+'/data_phasing/' + output, 'a') as fp:
        for items in coverage_segments:
            if not items == None:
                fp.write("%s\n" % items)
