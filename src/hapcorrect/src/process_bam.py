import pysam
import numpy as np
import subprocess
import pathlib
import os
import pandas
from collections import defaultdict

class ReadSegment(object):
    __slots__ = ("read_start", "read_end", "ref_start", "ref_end", "read_id", "ref_id",
                 "strand", "read_length", 'segment_length', "haplotype", "mapq", "genome_id", 'mismatch_rate',
                 "is_insertion")

    def __init__(self, read_start, read_end, ref_start, ref_end, read_id, ref_id,
                 strand, read_length, segment_length, haplotype, mapq, genome_id, mismatch_rate, is_insertion):
        self.read_start = read_start
        self.read_end = read_end
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_id = read_id
        self.ref_id = ref_id
        self.strand = strand
        self.read_length = read_length
        self.segment_length = segment_length
        self.haplotype = haplotype
        self.mapq = mapq
        self.genome_id = genome_id
        self.mismatch_rate = mismatch_rate
        self.is_insertion = is_insertion

    def __str__(self):
        return "".join(
            ["read_start=", str(self.read_start), " read_end=", str(self.read_end), " ref_start=", str(self.ref_start),
             " ref_end=", str(self.ref_end), " read_id=", str(self.read_id), " ref_id=", str(self.ref_id), " strand=",
             str(self.strand),
             " read_length=", str(self.read_length), " haplotype=", str(self.haplotype),
             " mapq=", str(self.mapq), " genome_id=", str(self.genome_id)])

class GenomicSegment(object):
    __slots__ = "genome_id","haplotype", "ref_id", "dir1", "pos1", "dir2" , "pos2", "coverage", "length_bp"
    def __init__(self, genome_id, haplotype, ref_id, pos1, pos2, coverage, length_bp):
        self.genome_id = genome_id
        self.haplotype = haplotype
        self.ref_id = ref_id
        self.dir1 = '-'
        self.pos1 = pos1
        self.dir2 = '+'
        self.pos2 = pos2
        self.coverage = coverage
        self.length_bp = length_bp
def get_allsegments(segments_by_read_filtered):
    allsegments = []
    for read in segments_by_read_filtered:
        for seg in read:
            if not seg.is_insertion:
                allsegments.append(seg)
    return allsegments
def filter_all_reads(segments_by_read, min_mapq, max_read_error, arguments):
    #MIN_ALIGNED_LENGTH = 5000
    MIN_ALIGNED_RATE = 0.5
    MAX_SEGMENTS = 10
    MIN_SEGMENT_LENGTH = 100

    segments_by_read_filtered = []
    for read_id, segments in segments_by_read.items():
        dedup_segments = []
        segments.sort(key=lambda s: s.read_start)

        for seg in segments:
            if not dedup_segments or dedup_segments[-1].read_start != seg.read_start:
                if seg.is_insertion and seg.mapq > min_mapq:
                    dedup_segments.append(seg)
                elif seg.mapq > min_mapq and seg.segment_length > MIN_SEGMENT_LENGTH and seg.mismatch_rate < max_read_error:
                    dedup_segments.append(seg)

        aligned_len = sum([seg.segment_length for seg in dedup_segments if not seg.is_insertion])
        aligned_ratio = aligned_len / segments[0].read_length
        if aligned_len < arguments['min_aligned_length'] or aligned_ratio < MIN_ALIGNED_RATE or len(segments) > MAX_SEGMENTS:
            continue
        segments_by_read_filtered.append(dedup_segments)

    return segments_by_read_filtered
def get_segment(read, genome_id, sv_size):
    """
    Parses cigar and generate ReadSegment structure with alignment coordinates
    """
    CIGAR_MATCH = [0, 7, 8]
    CIGAR_MM = 8
    CIGAR_DEL = 2
    CIGAR_INS = 1
    CIGAR_CLIP = [4, 5]
    first_clip = True
    read_start = 0
    read_aligned = 0
    read_length = 0
    ref_aligned = 0
    ref_start = read.reference_start
    read_segments = []
    cigar = read.cigartuples
    read_length = np.sum([k[1] for k in cigar if k[0] != 2])
    num_of_mismatch = 0
    strand = '-' if read.is_reverse else '+'
    if read.has_tag('HP'):
        haplotype = read.get_tag('HP')
    else:
        haplotype = 0
    for token in cigar:
        op = token[0]
        op_len = token[1]
        if op in CIGAR_CLIP:
            if first_clip:
                read_start = op_len
        first_clip = False
        if op in CIGAR_MATCH:
            read_aligned += op_len
            ref_aligned += op_len
            if op == CIGAR_MM:
                num_of_mismatch += 1
        if op == CIGAR_DEL:
            if op_len < sv_size:
                ref_aligned += op_len
            elif op_len > sv_size:
                ref_end = ref_start + ref_aligned
                read_end = read_start + read_aligned
                if read.is_reverse:
                    del_start, del_end = read_length - read_end, read_length - read_start
                else:
                    del_start, del_end = read_start, read_end
                mm_rate = num_of_mismatch / read_aligned
                read_segments.append(ReadSegment(del_start, del_end, ref_start, ref_end, read.query_name,
                                                 read.reference_name, strand, read_length, read_aligned, haplotype,
                                                 read.mapping_quality, genome_id, mm_rate, False))
                read_start = read_end + 1
                ref_start = ref_end + op_len + 1
                read_aligned = 0
                ref_aligned = 0
                num_of_mismatch = 0
        if op == CIGAR_INS:
            if op_len < sv_size:
                read_aligned += op_len
            else:
                ins_start = read_start + read_aligned
                read_aligned += op_len
                ins_pos = ref_start + ref_aligned
                ins_end = read_start + read_aligned
                mm_rate = 0
                read_segments.append(ReadSegment(ins_start, ins_end, ins_pos, ins_pos, read.query_name,
                                                 read.reference_name, strand, read_length, op_len, haplotype,
                                                 read.mapping_quality, genome_id, mm_rate, True))
    if ref_aligned != 0:
        ref_end = ref_start + ref_aligned
        read_end = read_start + read_aligned
        mm_rate = num_of_mismatch / read_aligned
        if read.is_reverse:
            read_start, read_end = read_length - read_end, read_length - read_start
        read_segments.append(ReadSegment(read_start, read_end, ref_start, ref_end, read.query_name,
                                         read.reference_name, strand, read_length, read_aligned, haplotype,
                                         read.mapping_quality, genome_id, mm_rate, False))
    return read_segments


def get_all_reads(bam_file, region, genome_id, sv_size):
    """
    Yields set of split reads for each contig separately. Only reads primary alignments
    and infers the split reads from SA alignment tag
    """
    alignments = []
    ref_id, region_start, region_end = region
    aln_file = pysam.AlignmentFile(bam_file, "rb")
    for aln in aln_file.fetch(ref_id, region_start, region_end, multiple_iterators=True):
        if not aln.is_secondary and not aln.is_unmapped:
            new_segment = get_segment(aln, genome_id, sv_size)
            for new_seg in new_segment:
                alignments.append(new_seg)
    return alignments


def get_all_reads_parallel(bam_file, thread_pool, ref_lengths,
                           min_mapq, genome_id, sv_size):
    CHUNK_SIZE = 10000000
    all_reference_ids = [r for r in pysam.AlignmentFile(bam_file, "rb").references]
    fetch_list = []
    for ctg in all_reference_ids:
        ctg_len = ref_lengths[ctg]
        for i in range(0, max(ctg_len // CHUNK_SIZE, 1)):
            reg_start = i * CHUNK_SIZE
            reg_end = (i + 1) * CHUNK_SIZE
            if ctg_len - reg_end < CHUNK_SIZE:
                reg_end = ctg_len
            fetch_list.append((ctg, reg_start, reg_end))

    tasks = [(bam_file, region, genome_id, sv_size) for region in fetch_list]
    parsing_results = None
    # thread_pool = Pool(num_threads)
    parsing_results = thread_pool.starmap(get_all_reads, tasks)
    segments_by_read = defaultdict(list)
    for alignments in parsing_results:
        for aln in alignments:
            segments_by_read[aln.read_id].append(aln)
    return segments_by_read

def haplotype_update_all_bins_parallel(bam_file, thread_pool, bins, bin_size):
    tasks = [(bam_file, region) for region in bins]
    parsing_results = thread_pool.starmap(process_all_reads, tasks)

def process_all_reads(histograms, bam_file, genome_id, region):
    ref_id, region_start, region_end, haplotype_1, haplotype_2 = region

    haplotype_1_coverage = segment_coverage(histograms, genome_id, ref_id, region_start, region_end, 1)
    haplotype_2_coverage = segment_coverage(histograms, genome_id, ref_id, region_start, region_end, 2)

    aln_file = pysam.AlignmentFile(bam_file, "rb")
    out_file = pysam.AlignmentFile(ref_id+'_'+region_start+'.bam', 'wb', aln_file.header)
    for aln in aln_file.fetch(ref_id, region_start, region_end, multiple_iterators=True):
        if not haplotype_1 == haplotype_1_coverage:
            if aln.tag.HP == 1:
                aln.set_tag('HP', 2)
            elif aln.tag.HP == 2:
                aln.set_tag('HP', 1)

            out_file.write(aln)

def segment_coverage(histograms, genome_id, ref_id, ref_start, ref_end, haplotype):
    hist_start = ref_start // COV_WINDOW
    hist_end = ref_end // COV_WINDOW
    cov_list = histograms[(genome_id, haplotype, ref_id)][hist_start : hist_end + 1]
    if not cov_list:
        return 0
    return round(np.mean(cov_list), 2)

def get_segments_coverage(segments, coverage_histograms):
    genomic_segments = []
    for (genome_id, seg_ref, seg_start, seg_end) in segments: # TODO Parallize it through threads pool with tasks
        coverage_hp1 = segment_coverage(coverage_histograms, genome_id, seg_ref, seg_start, seg_end, 1)
        coverage_hp2 = segment_coverage(coverage_histograms, genome_id, seg_ref, seg_start, seg_end, 2)
        coverage_hp0 = segment_coverage(coverage_histograms, genome_id, seg_ref, seg_start, seg_end, 0)
        genomic_segments.append(seg_ref + '\t' + str(seg_start) + '\t' + str(seg_end) + '\t' + str(coverage_hp1)+ '\t' + str(coverage_hp2)+ '\t' + str(coverage_hp0))
    return genomic_segments

COV_WINDOW = 500

def update_coverage_hist(genome_ids, ref_lengths, segments_by_read, min_mapq, max_read_error, arguments):
    NUM_HAPLOTYPES = 3
    segments_by_read_filtered = filter_all_reads(segments_by_read, min_mapq, max_read_error, arguments)
    allsegments = get_allsegments(segments_by_read_filtered)
    coverage_histograms = {}
    for genome_id in genome_ids:
        for chr_id, chr_len in ref_lengths.items():
            for hp in range(0, NUM_HAPLOTYPES):
                coverage_histograms[(genome_id, hp, chr_id)] = [0 for _ in range(chr_len // COV_WINDOW + 1)]
    for read in allsegments:
        hist_start = read.ref_start // COV_WINDOW
        hist_end = read.ref_end // COV_WINDOW
        for i in range(hist_start, hist_end + 1):  ## Check with (hist_start, hist_end + 1)
            coverage_histograms[(read.genome_id, read.haplotype, read.ref_id)][i] += 1

    return coverage_histograms

def get_snps_frequencies(bam, snp_segments, thread_pool):
    tasks = [(bam, region) for region in snp_segments]
    snps_results = None
    snps_results = thread_pool.starmap(compute_snp_frequency, tasks)
    return snps_results

def compute_snp_frequency(bam, region):
    from collections import Counter
    contig, start, ref, alt, gt = region
    bam = pysam.AlignmentFile(bam, 'rb')
    bases = []
    for pileupcolumn in bam.pileup(contig, start - 1, start, truncate=True): #truncate=True #, fastafile=fasta
        #pileupcolumn.set_min_base_quality(0)
        #base = pileupcolumn.get_query_sequences()
        #print('coverage at base %s = %s' % (pileupcolumn.pos, pileupcolumn.n))
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                bases.append(pileupread.alignment.query_sequence[pileupread.query_position])

    if not bases:
        return None

    #print(region)
    acgts = {}
    acgts['A'] = Counter(bases)['A']
    acgts['C'] = Counter(bases)['C']
    acgts['G'] = Counter(bases)['G']
    acgts['T'] = Counter(bases)['T']

    hp = 0
    if gt == '1|0': #ref freqs
        hp = 1
    elif gt == '0|1': #alt freqs
        hp = 2

    ref_value_new = acgts.get(ref)
    alt_value_new = acgts.get(alt)

    return (contig+'\t'+str(start)+'\t'+ref+'\t'+alt+'\t'+str(ref_value_new)+'\t'+str(alt_value_new)+'\t'+str(hp))

def tumor_bam_haplotag(arguments, out_vcf):
    basefile = pathlib.Path(arguments['target_bam'][0]).stem
    output_bam = f"{os.path.join('data', basefile + arguments['genome_name'] + '.rehaplotagged.bam')}"
    whatshap_cmd = ['whatshap', 'haplotag', '--reference', arguments['reference'], out_vcf, arguments['target_bam'][0], '- o', output_bam, '--ignore-read-groups', '--tag-supplementary', '--skip-missing-contigs', '--output-threads', str(arguments['threads'])]

    wh_1 = subprocess.Popen(whatshap_cmd, stdout=subprocess.PIPE)
    wh_1.wait()
    if wh_1.returncode != 0:
        raise ValueError('whatshap haplotag subprocess returned nonzero value: {}'.format(wh_1.returncode))
def process_bam_for_snps_freqs(arguments, thread_pool):
    basefile = pathlib.Path(arguments['target_bam'][0]).stem
    output_bam = f"{os.path.join('data', basefile + '_reduced.bam')}"

    samtools_cmd = ['samtools', 'view', '-@', str(arguments['threads']), '-F', '3844', '-q', '5', '-h', arguments['target_bam'][0]]
    awk_cmd = ['awk', '-v', 'OFS=\t', '{if($0 ~ /^@/){print $0} else {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, "*"}}']
    samtools_bam_cmd = ['samtools', 'view', '-@', str(arguments['threads']), '-Sb', '-o', output_bam]

    st_1 = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE)
    awk = subprocess.Popen(awk_cmd, stdin=st_1.stdout, stdout=subprocess.PIPE)
    st_2 = subprocess.Popen(samtools_bam_cmd, stdin=awk.stdout, stdout=subprocess.PIPE)

    st_1.wait()
    awk.wait()
    st_2.wait()

    if st_1.returncode != 0:
        raise ValueError('samtools view subprocess returned nonzero value: {}'.format(st_1.returncode))
    if awk.returncode != 0:
        raise ValueError('awk subprocess returned nonzero value: {}'.format(awk.returncode))
    if st_2.returncode != 0:
        raise ValueError('samtools view for bam output subprocess returned nonzero value: {}'.format(st_2.returncode))

    basefile = pathlib.Path(arguments['normal_phased_vcf']).stem
    output_csv = basefile + '_het_snps.csv'
    output_csv = f"{os.path.join('data', output_csv)}"

    output_vcf = basefile + '_het_phased_snps.vcf.gz'
    output_vcf = f"{os.path.join('data', output_vcf)}"

    cmd = ['bcftools', 'view', '--threads', str(arguments['threads']),  '--phased', '-g', 'het', '--types', 'snps', arguments['normal_phased_vcf'], '-o', output_vcf]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    cmd = ['bcftools', 'query', '-i', 'GT="het"', '-f',  '%CHROM\t%POS\n', output_vcf,  '-o', output_csv]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    beds = split_file(output_csv, arguments['threads'])
    pileups_outputs = process_pileups(output_bam, arguments['reference'], beds, thread_pool)

    output_pileup = f"{os.path.join('data', arguments['genome_name'] + '_SNPs.csv')}"

    for i in pileups_outputs:
        with open(i, 'r') as content_file:
            content = content_file.read()
        with open(output_pileup, 'a') as target_device:
            target_device.write(content)

    return output_pileup
def split_file(fname, parts):
    out_beds = []
    lines  = pandas.read_csv(fname)
    df_iterator = pandas.read_csv(fname, chunksize=int(len(lines)/parts) + 1, sep='\t', index_col=False, header=None)
    for i, chunk in enumerate(df_iterator):
        fname = f"{os.path.join('data', '%d.bed' % i)}"
        out_beds.append(fname)
        chunk.to_csv(fname, index=False, sep='\t', header=None)
    return out_beds

def process_pileups(bam, ref, input_beds, thread_pool):
    tasks = [(bam, ref, bed) for bed in input_beds]
    pileups_outputs = thread_pool.starmap(process_pileups_parallel, tasks)
    return pileups_outputs

def process_pileups_parallel(bam, ref, bed):
    basefile = pathlib.Path(bed).stem
    output_csv = f"{os.path.join('data', basefile + '_SNPs.csv')}"
    cmd = ['samtools', 'mpileup', '-l', bed, '-f', ref, bam,  '-o', output_csv]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process.wait()

    return output_csv
