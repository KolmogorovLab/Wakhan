#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Synchronize the global haplotype orientation of a phased germline VCF with a
Wakhan integer copy-number profile ("integer_profile.bed"), using TUMOR BAM
coverage (not the germline VCF's own depth fields) as the evidence. This does
NOT re-phase anything: phase block boundaries and which variants are phased
vs. unphased are left exactly as they are. It only decides, per phase block,
whether the whole block should be flipped (GT 0|1 <-> 1|0) so that
"haplotype 1" in the VCF consistently means the same physical haplotype as
"haplotype 1" in the Wakhan profile.

Why tumor BAM, not VCF depth: the VCF here is typically called on the normal
sample, so its own AD/DP reflect normal-sample coverage, which carries no CNA
signal. The tumor's allele-specific coverage at these same het SNV positions
is what's needed, so it's read directly from the tumor BAM via
`samtools mpileup`, using the exact conventions (-q 5 -Q 1, no BAQ, no
indels/read-ends) already used elsewhere in Wakhan
(src/file_tools/process_bam.py:process_pileups_parallel).

Algorithm (per phase block, grouped by the VCF's PS tag):

  A phase block is compared against the Wakhan profile one CN segment at a
  time, not as a single aggregate, so that an ordinary CN breakpoint inside
  the block (CN direction legitimately flips between two adjacent segments)
  is not mistaken for a phasing problem.

  0. Extract phased, heterozygous, single-base (true SNP) variants from the
     VCF, grouped by (chrom, PS). GT decides which allele - REF or ALT -
     belongs to HP1 vs HP2; it does not depend on which BAM the depth comes
     from.

  1. For each such SNV, look up TUMOR allele depth via `samtools mpileup`
     restricted to exactly these positions, and split it into
     (hp1_depth, hp2_depth) using the VCF's GT for that variant.

  2. For each Wakhan profile segment overlapping the block (clipped to the
     block span, centromere/blacklisted segments excluded, CN-balanced
     segments excluded since they carry no directional signal):
       - profile_sign = sign(hp1_copynumber_state - hp2_copynumber_state)
       - vcf_sign = majority vote of sign(hp1_depth - hp2_depth) over the
         SNVs in this segment's span (requires >= --min-snvs-per-segment
         usable, non-tied SNVs)
       - agreement = +1 if vcf_sign == profile_sign (this segment says the
         block's current orientation already matches Wakhan's HP1/HP2),
         else -1 (this segment votes to flip)

  3. Aggregate agreement across all segments that produced a vote, weighted
     by the number of SNVs used in each segment's vote:
       - W_keep / W_flip = total SNV-count weight of segments voting
         keep / flip
       - evaluable_fraction = (W_keep + W_flip) / (total usable-depth phased
         SNVs in the block). Below --min-evaluable-fraction, the block is
         skipped - not enough of it produced a usable comparison.
       - minority_fraction = min(W_keep, W_flip) / (W_keep + W_flip). Above
         --mixed-tolerance, the block is "mixed": segments disagree on the
         flip/keep call even after normalizing for each segment's own CN
         direction - a genuine signal of an internal switch error (not just
         an ordinary CN breakpoint). Mixed blocks are, by default, left
         untouched but reported. With --cut-mixed-blocks, the block is split
         at the boundary between adjacent segments whose agreement differs,
         and each resulting sub-block (with a freshly assigned PS) is
         evaluated independently.

  4. Decision: if not mixed and evaluable, flip the block (swap GT for every
     phased variant in it) iff W_flip > W_keep. PS is otherwise left as-is
     (except for new sub-blocks created by --cut-mixed-blocks).

Sanity checks (always on, printed to stderr, non-fatal):
  - depth-extraction coverage: what fraction of candidate SNVs yielded usable
    tumor depth, and why the rest were dropped
  - unexpected-allele fraction: average fraction of tumor pileup depth at
    these SNV sites that matches neither the REF nor any annotated ALT
    allele - high values (>15%) suggest a BAM/reference/VCF mismatch (wrong
    sample, genome build mismatch, coordinate issue)
  - decision confidence: SNV-weighted mean minority fraction across decided
    blocks - real signal on a profile with sustained CN imbalance should be
    lopsided, not close to a coin flip

Requires `samtools` on PATH.

Usage:
    python scripts/sync_vcf_phase_with_cna.py \\
        --vcf normal.phased.vcf.gz \\
        --profile solution_2.0_0.87_0.95/integer_profile.bed \\
        --tumor-bam tumor.bam \\
        --reference genome.fa \\
        --output normal.phased.synced.vcf \\
        --report sync_report.tsv
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict

import pandas as pd


def _log(msg):
    sys.stderr.write(f"[{time.strftime('%H:%M:%S')}] {msg}\n")
    sys.stderr.flush()

try:
    import pysam
except ImportError:
    sys.exit("This script requires 'pysam' (run inside the wakhan conda environment).")

try:
    from intervaltree import IntervalTree, Interval
except ImportError:
    sys.exit("This script requires 'intervaltree' (run inside the wakhan conda environment).")

# Must match CENTROMERE_SENTINEL in src/utils/chromosome.py - marks centromeric /
# blacklisted segments that carry no real CN information.
CENTROMERE_SENTINEL = 3300


def _read_profile_bed(path):
    """Read integer_profile.bed, recovering column names from its embedded
    '#chr\\tstart\\t...' header comment line. Returns {chrom: IntervalTree of
    Interval(start, end, (hp1_copynumber_state, hp2_copynumber_state))}."""
    header_cols = None
    with open(path) as fp:
        for line in fp:
            if line.startswith('#') and '\t' in line:
                header_cols = line.lstrip('#').rstrip('\n').split('\t')
    if header_cols is None:
        sys.exit(f"Could not find a '#chr\\tstart\\t...' header line in {path}")

    df = pd.read_csv(path, sep='\t', comment='#', header=None, names=header_cols, na_filter=False)
    for col in ('hp1_copynumber_state', 'hp2_copynumber_state'):
        if col not in df.columns:
            sys.exit(f"{path} is missing '{col}' - this script expects the merged "
                      f"integer_profile.bed format (per-haplotype CN states already "
                      f"reconstructed), not the subclonal profile.")

    trees = defaultdict(IntervalTree)
    for row in df.itertuples():
        if row.end <= row.start:
            continue
        trees[row.chr].add(Interval(row.start, row.end, (row.hp1_copynumber_state, row.hp2_copynumber_state)))
    return trees


def _is_masked(hp1_state, hp2_state):
    return abs(hp1_state - CENTROMERE_SENTINEL) < 1e-6 or abs(hp2_state - CENTROMERE_SENTINEL) < 1e-6


def _profile_segments_in_range(tree, start, end):
    """Non-masked profile segments overlapping [start, end), clipped, sorted
    by position: list of (seg_start, seg_end, hp1_state, hp2_state)."""
    segs = []
    for iv in tree[start:end]:
        ovl_start, ovl_end = max(iv.begin, start), min(iv.end, end)
        if ovl_end <= ovl_start:
            continue
        hp1_state, hp2_state = iv.data
        if _is_masked(hp1_state, hp2_state):
            continue
        segs.append((ovl_start, ovl_end, hp1_state, hp2_state))
    return sorted(segs, key=lambda s: s[0])


# --- VCF: extract phased het SNVs (no depth here - depth comes from the tumor BAM) ---

def _extract_phased_snvs(vcf_path, sample_arg):
    """Returns (sample, snvs, stats) where snvs is a list of
    (chrom, pos, ps, gt, alleles) for phased, heterozygous, single-base
    (true SNP) variants. alleles = (ref_base, alt_base_1, ...); gt indexes
    into it. stats tallies why candidate variants were excluded."""
    pysam_verbosity = pysam.set_verbosity(0)
    vcf_reader = pysam.VariantFile(vcf_path, 'r')
    sample = sample_arg or list(vcf_reader.header.samples)[0]

    stats = Counter()
    snvs = []
    for var in vcf_reader:
        s = var.samples[sample]
        if not s.phased:
            continue
        ps = s.get('PS')
        if ps is None:
            continue
        gt = s['GT']
        if gt is None or len(gt) != 2 or gt[0] is None or gt[1] is None or gt[0] == gt[1]:
            continue
        stats['phased_het'] += 1

        alleles = (var.ref,) + var.alts
        if gt[0] >= len(alleles) or gt[1] >= len(alleles) or \
           len(alleles[gt[0]]) != 1 or len(alleles[gt[1]]) != 1:
            stats['not_snp'] += 1
            continue

        snvs.append((var.chrom, var.pos, ps, gt, alleles))

    vcf_reader.close()
    pysam.set_verbosity(pysam_verbosity)
    return sample, snvs, stats


# --- Tumor BAM: samtools mpileup at exactly the candidate SNV positions ---

def _write_positions_bed(snvs, path):
    with open(path, 'w') as fp:
        for chrom, pos, ps, gt, alleles in sorted(snvs, key=lambda r: (r[0], r[1])):
            fp.write(f"{chrom}\t{pos - 1}\t{pos}\n")


def _parse_pileup_line(line):
    """Port of src/file_tools/process_vcf.py:parse_pileup_line (copied
    verbatim rather than imported, to keep this script standalone)."""
    chrom, pos, ref, depth, bases, qual = line.strip().split("\t")[:6]
    pos = int(pos)

    i = 0
    clean_bases = ""
    while i < len(bases):
        c = bases[i]
        if c == '^':  # start of a read segment, skip the next character (mapping quality)
            i += 2
        elif c == '$':  # end of a read segment
            i += 1
        elif c in '+-':  # insertion or deletion
            i += 1
            indel_len = ''
            while i < len(bases) and bases[i].isdigit():
                indel_len += bases[i]
                i += 1
            i += int(indel_len)
        else:
            clean_bases += c
            i += 1

    base_map = {'.': ref.upper(), ',': ref.upper(), 'A': 'A', 'a': 'A',
                'C': 'C', 'c': 'C', 'G': 'G', 'g': 'G', 'T': 'T', 't': 'T'}
    converted = [base_map.get(b, 'N') for b in clean_bases if base_map.get(b, 'N') in 'ACGT']
    return chrom, pos, Counter(converted)


def _run_mpileup(chrom, bed_path, reference, bam, min_mapq, min_baseq, out_path):
    cmd = ['samtools', 'mpileup', '-l', bed_path, '-f', reference, '-r', chrom, bam,
           '--no-BAQ', '-q', str(min_mapq), '-Q', str(min_baseq),
           '--no-output-ins', '--no-output-del', '--no-output-ends']
    with open(out_path, 'w') as out_fp:
        subprocess.run(cmd, stdout=out_fp, stderr=subprocess.DEVNULL, check=True)


def _collect_tumor_pileup_counts(snvs, bed_path, reference, bam, min_mapq, min_baseq, threads, tmpdir):
    from concurrent.futures import ThreadPoolExecutor, as_completed

    chroms = sorted(set(chrom for chrom, *_ in snvs))
    out_paths = {chrom: f"{tmpdir}/{chrom}.pileup" for chrom in chroms}

    _log(f"samtools mpileup: {len(chroms)} chromosome(s), {threads} parallel worker(s)")
    n_done = 0
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {pool.submit(_run_mpileup, chrom, bed_path, reference, bam,
                                min_mapq, min_baseq, out_paths[chrom]): chrom for chrom in chroms}
        for future in as_completed(futures):
            chrom = futures[future]
            try:
                future.result()
            except subprocess.CalledProcessError as e:
                sys.exit(f"samtools mpileup failed for {chrom}: {e}")
            n_done += 1
            _log(f"  mpileup done: {chrom} ({n_done}/{len(chroms)})")

    _log("Parsing mpileup output...")
    pileup_counts = {}
    for chrom in chroms:
        with open(out_paths[chrom]) as fp:
            for line in fp:
                _, pos, counts = _parse_pileup_line(line)
                pileup_counts[(chrom, pos)] = counts
    return pileup_counts


def _hp_depths_from_pileup(chrom, pos, gt, alleles, pileup_counts, stats):
    counts = pileup_counts.get((chrom, pos))
    if counts is None:
        stats['no_pileup_data'] += 1
        return None

    depth_by_allele = {i: counts.get(base.upper(), 0) for i, base in enumerate(alleles) if len(base) == 1}
    hp1_depth, hp2_depth = depth_by_allele.get(gt[0]), depth_by_allele.get(gt[1])
    if hp1_depth is None or hp2_depth is None:
        stats['gt_allele_not_in_depths'] += 1
        return None
    if hp1_depth + hp2_depth == 0:
        stats['zero_relevant_depth'] += 1
        return None

    known_bases = {base.upper() for base in alleles if len(base) == 1}
    total_depth = sum(counts.values())
    other_depth = total_depth - sum(counts.get(b, 0) for b in known_bases)
    stats['other_depth_sum'] += other_depth
    stats['total_depth_sum'] += total_depth

    stats['usable'] += 1
    return hp1_depth, hp2_depth


def _snv_signs_in_range(snvs, start, end):
    """signs (list of +1/-1, ties excluded) and total usable-depth SNV count
    (including ties) for snvs whose pos falls in [start, end)."""
    signs = []
    n = 0
    for pos, hp1_depth, hp2_depth in snvs:
        if not (start <= pos < end):
            continue
        n += 1
        if hp1_depth > hp2_depth:
            signs.append(1)
        elif hp2_depth > hp1_depth:
            signs.append(-1)
    return signs, n


def _sign_vote(signs):
    pos = sum(1 for s in signs if s > 0)
    neg = sum(1 for s in signs if s < 0)
    if pos == neg:
        return None
    return 1 if pos > neg else -1


def _evaluate_segments(tree, snvs, start, end, min_snvs_per_segment):
    """Per-segment agreement votes for [start, end), weighted by SNV count
    used in each segment's vote. Returns (votes, n_block_snvs) where votes is
    a list of {'start','end','weight','agreement'} and n_block_snvs is the
    total usable-depth phased SNV count in [start, end)."""
    _, n_block_snvs = _snv_signs_in_range(snvs, start, end)

    votes = []
    for s_start, s_end, hp1_state, hp2_state in _profile_segments_in_range(tree, start, end):
        diff = hp1_state - hp2_state
        if abs(diff) < 1e-9:
            continue  # CN-balanced segment, no directional signal
        profile_sign = 1 if diff > 0 else -1

        signs, n_seg = _snv_signs_in_range(snvs, s_start, s_end)
        if n_seg < min_snvs_per_segment:
            continue
        vcf_sign = _sign_vote(signs)
        if vcf_sign is None:
            continue

        agreement = 1 if vcf_sign == profile_sign else -1
        votes.append({'start': s_start, 'end': s_end, 'weight': n_seg, 'agreement': agreement})

    return votes, n_block_snvs


def _aggregate(votes):
    w_keep = sum(v['weight'] for v in votes if v['agreement'] == 1)
    w_flip = sum(v['weight'] for v in votes if v['agreement'] == -1)
    return w_keep, w_flip


def _cut_points_from_votes(votes):
    points = []
    prev = None
    for v in sorted(votes, key=lambda v: v['start']):
        if prev is not None and v['agreement'] != prev:
            points.append(v['start'])
        prev = v['agreement']
    return points


def _split_range(start, end, cut_points):
    bounds = [start] + sorted(set(p for p in cut_points if start < p < end)) + [end]
    return list(zip(bounds[:-1], bounds[1:]))


def _decide(tree, snvs, start, end, args):
    """Evaluate [start, end) and return a decision dict with 'decision'
    ('flip'/'keep'), 'reason' (None if an actual comparison was made), and
    diagnostic fields for the report."""
    _, n_snvs = _snv_signs_in_range(snvs, start, end)
    row = {'n_snvs': n_snvs}

    if n_snvs < args.min_snvs:
        row.update(decision='keep', reason='too_few_snvs', flip=False)
        return row

    votes, n_block_snvs = _evaluate_segments(tree, snvs, start, end, args.min_snvs_per_segment)
    w_keep, w_flip = _aggregate(votes)
    w_total = w_keep + w_flip

    if w_total == 0 or n_block_snvs == 0:
        row.update(decision='keep', reason='no_profile_data', flip=False)
        return row

    evaluable_fraction = w_total / n_block_snvs
    row['evaluable_fraction'] = round(evaluable_fraction, 3)
    if evaluable_fraction < args.min_evaluable_fraction:
        row.update(decision='keep', reason='insufficient_evaluable_evidence', flip=False)
        return row

    minority_fraction = min(w_keep, w_flip) / w_total
    row.update(w_keep=w_keep, w_flip=w_flip, minority_fraction=round(minority_fraction, 3))

    if minority_fraction > args.mixed_tolerance:
        row.update(decision='keep', reason='mixed_direction', flip=False, votes=votes)
        return row

    flip = w_flip > w_keep
    row.update(decision='flip' if flip else 'keep', reason=None, flip=flip)
    return row


def compute_decisions(blocks, profile_trees, args, report_rows):
    """blocks: {(chrom, ps): [(pos, hp1_depth, hp2_depth), ...]}
    Returns {chrom: IntervalTree of Interval(start, end, {'flip':bool,'new_ps':int|None})}"""
    decisions = defaultdict(IntervalTree)

    for (chrom, ps), snvs in blocks.items():
        block_start = min(p for p, _, _ in snvs)
        block_end = max(p for p, _, _ in snvs) + 1
        tree = profile_trees.get(chrom, IntervalTree())

        sub_ranges = [(block_start, block_end, ps)]
        if args.cut_mixed_blocks:
            top = _decide(tree, snvs, block_start, block_end, args)
            if top['reason'] == 'mixed_direction':
                cuts = _cut_points_from_votes(top['votes'])
                if cuts:
                    pieces = _split_range(block_start, block_end, cuts)
                    sub_ranges = [(s, e, ps if i == 0 else s) for i, (s, e) in enumerate(pieces)]

        for start, end, sub_ps in sub_ranges:
            row = _decide(tree, snvs, start, end, args)
            row.pop('votes', None)
            row.update(chrom=chrom, orig_ps=ps, start=start, end=end, new_ps=sub_ps)
            decisions[chrom].add(Interval(start, end, {'flip': row['flip'], 'new_ps': sub_ps if sub_ps != ps else None}))
            report_rows.append(row)

    return decisions


def _report_extraction_diagnostics(vcf_stats, depth_stats):
    phased_het = vcf_stats.get('phased_het', 0)
    if phased_het == 0:
        sys.stderr.write("WARNING: no phased heterozygous SNVs with a PS tag were found in the VCF.\n")
        return

    not_snp = vcf_stats.get('not_snp', 0)
    n_snp_candidates = phased_het - not_snp
    _log(f"Phased het variants: {phased_het} ({not_snp} not single-base SNPs, skipped)")

    usable = depth_stats.get('usable', 0)
    if n_snp_candidates == 0:
        return
    usable_fraction = usable / n_snp_candidates
    _log(f"SNP candidates: {n_snp_candidates}, usable tumor depth: {usable} ({usable_fraction:.1%})")
    for reason in ('no_pileup_data', 'gt_allele_not_in_depths', 'zero_relevant_depth'):
        count = depth_stats.get(reason, 0)
        if count:
            sys.stderr.write(f"  dropped ({reason}): {count} ({count / n_snp_candidates:.1%})\n")

    if usable_fraction < 0.5:
        sys.stderr.write(
            "WARNING: fewer than 50% of candidate SNPs yielded usable tumor depth. Check that "
            "--tumor-bam and --reference correspond to this VCF (same sample/genome build), and "
            "that the tumor BAM is indexed.\n")

    total_depth_sum = depth_stats.get('total_depth_sum', 0)
    if total_depth_sum > 0:
        other_fraction = depth_stats.get('other_depth_sum', 0) / total_depth_sum
        _log(f"Mean fraction of tumor pileup depth matching neither REF nor ALT: {other_fraction:.1%}")
        if other_fraction > 0.15:
            sys.stderr.write(
                "WARNING: an unusually high fraction of tumor reads at these SNP sites carry a base "
                "that doesn't match REF or ALT. This usually means --tumor-bam/--reference don't "
                "actually correspond to this VCF (wrong sample, mismatched genome build, or a "
                "coordinate/liftover issue) - results are likely unreliable.\n")


def _report_confidence_diagnostics(report_rows):
    """Warn if decisions are, on average, only weakly confident (SNV-weighted
    mean minority fraction close to the coin-flip point). Real, correctly-
    parsed signal on a profile with sustained CN imbalance should produce
    mostly lopsided per-block votes; a genome-wide result hovering near 50/50
    is a sign of noise - either genuinely low-confidence input data, or a
    BAM/reference mismatch that survived the checks above."""
    decided = [r for r in report_rows if r.get('reason') is None]
    if not decided:
        return

    total_w = sum(r['w_keep'] + r['w_flip'] for r in decided)
    if total_w == 0:
        return
    weighted_minority = sum(r['minority_fraction'] * (r['w_keep'] + r['w_flip']) for r in decided) / total_w

    _log(f"Decided blocks: {len(decided)}, SNV-weighted mean minority fraction: {weighted_minority:.3f}")
    if weighted_minority > 0.3:
        sys.stderr.write(
            "WARNING: decisions are only weakly confident on average (high minority fraction across "
            "decided blocks). This can indicate genuinely noisy/low-coverage data, or a BAM/reference "
            "mismatch - inspect --report for per-block detail before trusting the output.\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf', required=True, help='Phased germline VCF (must have PS tags; typically normal-sample)')
    parser.add_argument('--profile', required=True, help='Wakhan integer_profile.bed (merged HP1/HP2 format)')
    parser.add_argument('--tumor-bam', required=True, help='Indexed tumor BAM/CRAM to read allele depth from')
    parser.add_argument('--reference', required=True, help='Reference FASTA (faidx-indexed), required by samtools mpileup')
    parser.add_argument('--output', required=True, help='Output VCF path')
    parser.add_argument('--sample', default=None, help='VCF sample column to use (default: first sample in the VCF)')
    parser.add_argument('--min-mapping-quality', type=int, default=5,
                         help='samtools mpileup -q threshold. Default: 5 (matches Wakhan\'s own SNP pileup step)')
    parser.add_argument('--min-base-quality', type=int, default=1,
                         help='samtools mpileup -Q threshold. Default: 1 (matches Wakhan\'s own SNP pileup step)')
    parser.add_argument('--threads', type=int, default=4, help='Parallel samtools mpileup processes (one per chromosome). Default: 4')
    parser.add_argument('--min-snvs', type=int, default=3,
                         help='Minimum phased het SNVs (with usable tumor depth) required per block before any '
                              'comparison is attempted. Default: 3')
    parser.add_argument('--min-snvs-per-segment', type=int, default=2,
                         help='Minimum phased het SNVs required within a single profile segment for that '
                              "segment's vote to count. Default: 2")
    parser.add_argument('--min-evaluable-fraction', type=float, default=0.5,
                         help='Minimum fraction of a block\'s usable-depth SNVs that must fall in segments '
                              'that produced a vote (CN-imbalanced, non-masked, enough SNVs) before a flip '
                              'decision is made. Default: 0.5')
    parser.add_argument('--mixed-tolerance', type=float, default=0.15,
                         help='Max minority-vote SNV-count fraction still considered a clean majority; '
                              'above this the block is "mixed". Default: 0.15')
    parser.add_argument('--cut-mixed-blocks', action='store_true',
                         help='Split mixed-direction blocks at the boundary between adjacent segments whose '
                              'agreement differs, and evaluate the resulting sub-blocks independently. '
                              'Disabled by default.')
    parser.add_argument('--keep-intermediate', action='store_true',
                         help='Keep the temporary positions BED and mpileup output files instead of deleting them')
    parser.add_argument('--report', default=None, help='Optional path to write a per-block TSV decision report')
    args = parser.parse_args()

    if shutil.which('samtools') is None:
        sys.exit("'samtools' not found on PATH (required for tumor BAM pileup)")

    _log(f"Reading Wakhan profile: {args.profile}")
    profile_trees = _read_profile_bed(args.profile)
    _log(f"  {sum(len(t) for t in profile_trees.values())} CN segments across {len(profile_trees)} chromosome(s)")

    _log(f"Reading phased VCF: {args.vcf}")
    sample, snvs, vcf_stats = _extract_phased_snvs(args.vcf, args.sample)
    if not snvs:
        sys.exit("No phased, heterozygous, single-base SNVs with a PS tag found in the VCF")
    n_blocks_input = len(set((c, ps) for c, _, ps, _, _ in snvs))
    _log(f"  {len(snvs)} candidate SNPs in {n_blocks_input} phase block(s) (sample: {sample})")

    tmpdir_ctx = tempfile.TemporaryDirectory(prefix='sync_vcf_phase_')
    tmpdir = tmpdir_ctx.name
    if args.keep_intermediate:
        tmpdir_ctx = None  # don't auto-clean; leak the directory intentionally
        tmpdir = tempfile.mkdtemp(prefix='sync_vcf_phase_')
        _log(f"Keeping intermediate files in {tmpdir}")

    bed_path = f"{tmpdir}/positions.bed"
    _write_positions_bed(snvs, bed_path)

    _log(f"Reading tumor coverage from BAM: {args.tumor_bam}")
    pileup_counts = _collect_tumor_pileup_counts(
        snvs, bed_path, args.reference, args.tumor_bam,
        args.min_mapping_quality, args.min_base_quality, args.threads, tmpdir)

    if tmpdir_ctx is not None:
        tmpdir_ctx.cleanup()

    _log("Mapping tumor depth to haplotypes...")
    depth_stats = Counter()
    blocks = defaultdict(list)
    for chrom, pos, ps, gt, alleles in snvs:
        depths = _hp_depths_from_pileup(chrom, pos, gt, alleles, pileup_counts, depth_stats)
        if depths is None:
            continue
        blocks[(chrom, ps)].append((pos, depths[0], depths[1]))

    _report_extraction_diagnostics(vcf_stats, depth_stats)

    _log("Evaluating phase blocks against the CNA profile...")
    report_rows = []
    decisions = compute_decisions(blocks, profile_trees, args, report_rows)
    _report_confidence_diagnostics(report_rows)

    _log(f"Writing output VCF: {args.output}")
    vcf_reader = pysam.VariantFile(args.vcf, 'r')
    vcf_out = pysam.VariantFile(args.output, 'w', header=vcf_reader.header)

    n_flipped_variants = 0
    for var in vcf_reader:
        if var.samples[sample].phased and var.samples[sample].get('PS') is not None:
            hits = decisions[var.chrom][var.pos]
            if hits:
                decision = next(iter(hits)).data
                if decision['new_ps'] is not None:
                    var.samples[sample]['PS'] = decision['new_ps']
                if decision['flip']:
                    a, b = var.samples[sample]['GT']
                    var.samples[sample]['GT'] = (b, a)
                    var.samples[sample].phased = True
                    n_flipped_variants += 1
        vcf_out.write(var)

    vcf_reader.close()
    vcf_out.close()

    n_blocks = len(report_rows)
    counts = defaultdict(int)
    for row in report_rows:
        counts[row.get('reason') or row['decision']] += 1
    _log(f"Done: {n_blocks} block(s) processed, {n_flipped_variants} variant(s) flipped")
    for key, count in sorted(counts.items()):
        sys.stderr.write(f"  {key}: {count}\n")

    if args.report:
        cols = ['chrom', 'orig_ps', 'start', 'end', 'new_ps', 'n_snvs', 'evaluable_fraction',
                'w_keep', 'w_flip', 'minority_fraction', 'decision', 'reason']
        df = pd.DataFrame(report_rows)
        df = df[[c for c in cols if c in df.columns] + [c for c in df.columns if c not in cols and c != 'flip']]
        df.to_csv(args.report, sep='\t', index=False)
        _log(f"Wrote per-block report to {args.report}")


if __name__ == '__main__':
    main()
