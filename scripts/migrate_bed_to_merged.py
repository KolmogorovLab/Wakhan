#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Migrate old-format Wakhan per-haplotype CN segment bed files (separate HP1/HP2
files, e.g. "*_copynumbers_segments_HP_1.bed" / "*_HP_2.bed") into the new
merged single-bed format used since the output reorganization
("integer_profile.bed" / "subclonal_profile.bed").

The two haplotypes were segmented independently, so their segment boundaries
don't always line up. This script boundary-aligns them: any start/end present
in only one haplotype fragments the other haplotype's overlapping segment,
and both resulting fragments keep the parent segment's original
coverage/state/confidence. Breakpoint IDs are unioned across haplotypes. A
fragment covered by only one haplotype (should not normally happen) is
zero-filled for the missing haplotype rather than dropped.

Usage:
    python scripts/migrate_bed_to_merged.py \\
        --hp1 sample_2.0_0.87_0.95_copynumbers_segments_HP_1.bed \\
        --hp2 sample_2.0_0.87_0.95_copynumbers_segments_HP_2.bed \\
        --output integer_profile.bed
"""

import argparse
import ast
import sys

import pandas as pd

try:
    from intervaltree import IntervalTree, Interval
except ImportError:
    sys.exit("This script requires the 'intervaltree' package "
              "(pip install intervaltree, or run inside the wakhan conda environment).")


def _read_old_bed(path):
    """Read an old-format per-haplotype bed file, recovering its column names
    from the embedded '#chr\\tstart\\t...' header comment line."""
    header_cols = None
    with open(path) as fp:
        for line in fp:
            if line.startswith('#') and '\t' in line:
                header_cols = line.lstrip('#').rstrip('\n').split('\t')
    if header_cols is None:
        sys.exit(f"Could not find a '#chr\\tstart\\t...' header line in {path}")

    df = pd.read_csv(path, sep='\t', comment='#', header=None, names=header_cols, na_filter=False)

    if 'svs_breakpoints_ids' in df.columns:
        df['svs_breakpoints_ids'] = df['svs_breakpoints_ids'].apply(
            lambda x: set(ast.literal_eval(x)) if isinstance(x, str) and x.strip() else set()
        )

    return df


def merge_haplotype_segments(df_hp1, df_hp2, is_subclonal, has_breakpoints):
    merged_rows = []
    chroms = sorted(set(df_hp1['chr']).union(df_hp2['chr']))

    for chrom in chroms:
        h1 = df_hp1[df_hp1['chr'] == chrom]
        h2 = df_hp2[df_hp2['chr'] == chrom]

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
                    row['hp1_is_subclonal'] = r1.if_subclonal
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
                    row['hp2_is_subclonal'] = r2.if_subclonal
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


def _write_header(fp, columns, is_subclonal, has_breakpoints):
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


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--hp1', required=True, help='Old-format haplotype-1 CN segments bed file')
    parser.add_argument('--hp2', required=True, help='Old-format haplotype-2 CN segments bed file')
    parser.add_argument('--output', required=True,
                         help='Output path for the merged bed file (e.g. integer_profile.bed / subclonal_profile.bed)')
    args = parser.parse_args()

    df_hp1 = _read_old_bed(args.hp1)
    df_hp2 = _read_old_bed(args.hp2)

    is_subclonal = 'if_subclonal' in df_hp1.columns
    has_breakpoints = 'svs_breakpoints_ids' in df_hp1.columns

    if is_subclonal != ('if_subclonal' in df_hp2.columns):
        sys.exit("--hp1/--hp2 disagree on whether this is a subclonal profile (if_subclonal column)")
    if has_breakpoints != ('svs_breakpoints_ids' in df_hp2.columns):
        sys.exit("--hp1/--hp2 disagree on whether breakpoint IDs are present (svs_breakpoints_ids column)")

    merged_df = merge_haplotype_segments(df_hp1, df_hp2, is_subclonal, has_breakpoints)

    with open(args.output, 'w') as fp:
        _write_header(fp, merged_df.columns, is_subclonal, has_breakpoints)
        merged_df.to_csv(fp, sep='\t', index=False, header=False)

    print(f"Wrote {len(merged_df)} merged segments to {args.output}")


if __name__ == '__main__':
    main()
