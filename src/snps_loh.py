import statistics

import plotly.graph_objects as go
import os
import logging
import pandas as pd

from vcf_processing import get_snps_frquncies, het_homo_snps_gts, vcf_parse_to_csv_for_het_phased_snps_phasesets, snps_mean, cpd_mean, get_snp_segments, get_snps_frquncies_coverage, vcf_parse_to_csv_for_snps, get_snps_counts, get_snps_counts_cn_regions
from utils import csv_df_chromosomes_sorter
from extras import get_contigs_list
from plots import add_scatter_trace_coverage, print_chromosome_html, plots_add_markers_lines, plots_layout_settings,\
whole_genome_combined_df, copy_number_plots_per_chromosome
from utils import apply_copynumbers, detect_alter_loh_regions, loh_regions_events, write_segments_coverage, write_header_comments
def plot_snps_frequencies(arguments, df, df_segs_hp1_w, df_segs_hp2_w, centers, integer_fractional_means):
    if not os.path.isdir(arguments['out_dir_plots']+'/variation_plots'):
        os.mkdir(arguments['out_dir_plots']+'/variation_plots')
    filename = f"{os.path.join(arguments['out_dir_plots'], 'variation_plots/SNPs_LOH_CNs_INDEX.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")

    if arguments['tumor_vcf']:
        output_phasesets_file_path = vcf_parse_to_csv_for_snps(arguments['tumor_vcf'])
        df_snps_in_csv = csv_df_chromosomes_sorter(output_phasesets_file_path, ['chr', 'pos', 'qual', 'gt', 'dp', 'vaf'])

    chroms = get_contigs_list(arguments['contigs'])
    loh_regions_events_all = []
    df_snps_counts_per_cn_region_all = []

    for index, chrom in enumerate(chroms):
        df_chrom = df[df['chr'] == chrom]
        df_segs_hp1_ = df_segs_hp1_w[df_segs_hp1_w['chromosome'] == chrom]
        df_segs_hp2_ = df_segs_hp2_w[df_segs_hp2_w['chromosome'] == chrom]
        ref_start_values = df_chrom.start.values.tolist()
        ref_end_values = df_chrom.end.values.tolist()
        chr = range(len(ref_start_values))

        if arguments['without_phasing']:
            coverage = df_chrom.coverage.values.tolist()
            df_snps_freqs_chr = whole_genome_combined_df(arguments, chrom, chr, ref_start_values, ref_end_values, coverage, coverage, coverage)

        else:
            hp1 = df_chrom.hp1.values.tolist()
            hp2 = df_chrom.hp2.values.tolist()
            hp3 = df_chrom.hp3.values.tolist()
            df_snps_freqs_chr = whole_genome_combined_df(arguments, chrom, chr, ref_start_values, ref_end_values, hp1, hp2, hp3)

        if arguments['tumor_vcf']:
            snps_het_counts, snps_homo_counts, start_values, end_values = get_snps_counts_cn_regions(df_snps_in_csv, chrom, df_segs_hp1_.start.values.tolist(), df_segs_hp1_.end.values.tolist())
            df_snps_counts_per_cn_region_all.append(snps_counts_per_cn_region(snps_het_counts, snps_homo_counts, start_values, end_values, chrom))

            centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends = plot_snps(chrom, index, df_snps_in_csv, html_graphs, arguments, df_chrom)
            #hp1, hp2, hp3 = detect_alter_loh_regions(arguments, 'centromere/no-coverage', chrom, ref_end_values, hp1, hp2, hp3, centromere_region_starts, centromere_region_ends)
            if arguments['without_phasing']:
                hp1, hp2, hp3, loh_region_starts, loh_region_ends = detect_alter_loh_regions(arguments, 'loss-of-heterozygosity', chrom, ref_end_values, coverage, coverage, coverage, loh_region_starts, loh_region_ends, False)
            else:
                hp1, hp2, hp3, loh_region_starts, loh_region_ends = detect_alter_loh_regions(arguments, 'loss-of-heterozygosity', chrom, ref_end_values, hp1, hp2, hp3, loh_region_starts, loh_region_ends, False)
            loh_regions_events_all.extend(loh_regions_events(chrom, loh_region_starts, loh_region_ends, arguments))
        else:
            loh_region_starts = []
            loh_region_ends = []

        copy_number_plots_per_chromosome(centers, integer_fractional_means, ref_start_values, hp1, df_segs_hp1_, hp2, df_segs_hp2_, arguments, chrom, html_graphs, loh_region_starts, loh_region_ends)

    if arguments['tumor_vcf']:
        write_snps_counts_per_cn_region(pd.concat(df_snps_counts_per_cn_region_all), arguments)
        write_header_comments('chr\tstart\tend\n', '#chr: chromosome number\n#start: start address for LOH region\n#end: end address for LOH region\n', arguments['genome_name'] + '_loh_segments.bed', arguments)
        write_segments_coverage(loh_regions_events_all, arguments['genome_name'] + '_loh_segments.bed', arguments)

    html_graphs.write("</body></html>")

def plot_snps(chrom, index, df_snps_in_csv, html_graphs, arguments, df_chrom):
    logging.info('SNPs frequencies plots generation for ' + chrom)
    fig = go.Figure()

    snps_het, snps_homo, snps_het_pos, snps_homo_pos = get_snps_frquncies(df_snps_in_csv, chrom)

    #for debug purpose, all SNPs Freqs
    add_scatter_trace_coverage(fig, snps_het_pos, snps_het, name='Het SNPs Freqs', text=None, yaxis=None,
                               opacity=0.7, color='#E3B448', visibility='legendonly')
    add_scatter_trace_coverage(fig, snps_homo_pos, snps_homo, name='Homo SNPs Freqs', text=None, yaxis=None,
                               opacity=0.7, color='#3A6B35', visibility='legendonly')

    ref_start_values = [i for i in range(0, snps_homo_pos[-1:][0], arguments['bin_size_snps'])]
    ref_start_values_updated, snps_het_counts, snps_homo_counts, centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends = get_snps_frquncies_coverage(
        df_snps_in_csv, chrom, ref_start_values, arguments['bin_size_snps'])

    add_scatter_trace_coverage(fig, ref_start_values_updated, snps_het_counts, name='Het SNPs Ratios', text=None,
                               yaxis=None,
                               opacity=0.7, color='#E3B448')
    add_scatter_trace_coverage(fig, ref_start_values_updated, snps_homo_counts, name='Homo SNPs Ratios', text=None,
                               yaxis=None,
                               opacity=0.7, color='#3A6B35')

    plots_add_markers_lines(fig)
    plots_layout_settings(fig, chrom, arguments, snps_homo_pos[-1:][0], 0)
    fig.update_layout(yaxis_title="<b>SNPs Frequencies</b> (ratios)",)

    plot_snps_counts(chrom, index, ref_start_values, df_snps_in_csv, html_graphs, arguments, df_chrom)

    print_chromosome_html(fig, chrom + '_snps', html_graphs, arguments['out_dir_plots']+'/variation_plots/')
    html_graphs.write("  <object data=\"" + chrom + '_snps' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

    return centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends
def save_snps_counts_per_bin(hets,homos,ref_start, chrom, coverage, arguments, index):
    import pandas as pd
    fp = open(arguments['out_dir_plots']+'/bed_output/' + arguments['genome_name'] + '_snps_counts.bed', 'a')

    if index == 0:
        header = '#chr: chromosome number\n#start: start address for SNPs bin\n#hets_count: number of hetrozygous SNPs\n#homo_counts: number of homozygous SNPs\n#coverage: mean coverage in this bin\nchr\tstart\thets_count\thomos_count\tcoverage\n'
        fp.write(header)
    chr_list = [chrom for ch in range(len(ref_start))]
    snps_counts_df = pd.DataFrame(list(zip(chr_list, ref_start, hets, homos, coverage)), columns=['chr','start','hets_count','homos_count','coverage'])
    snps_counts_df.to_csv(fp, sep='\t', index=False, mode='a', header=False)

def snps_counts_per_cn_region(hets, homos, ref_start, ref_end, chrom):
    import pandas as pd
    chr_list = [chrom for ch in range(len(ref_start))]
    snps_counts_df = pd.DataFrame(list(zip(chr_list, ref_start, ref_end, hets, homos)),
                                columns=['chr', 'start', 'end', 'hets_count', 'homos_count'])
    return snps_counts_df

def write_snps_counts_per_cn_region(df, arguments):
    fp = open(arguments['out_dir_plots']+'/bed_output/' + arguments['genome_name'] + '_snps_counts_cn_regions.bed', 'a')
    fp.write('#chr: chromosome number\n')
    fp.write('#start: start address for CN segment\n')
    fp.write('#end: end address for CN segment\n')
    fp.write('#hets_count: hetrozygous SNPs count in this region\n')
    fp.write('#homos_count: homozygous SNPs count in this region\n')

    header = ['chr','start', 'end', 'hets_count','homos_count']
    df.to_csv(fp, sep='\t', columns=header, index=False, mode='a', header=True)

def coverage_bins(df):
    ref_start_values = df.start.values.tolist()
    ref_end_values = df.end.values.tolist()
    coverage = df.coverage.values.tolist()
    coverage_ = []

    for i in range(0, len(ref_start_values), 20):
        coverage_.append(statistics.mean([j for j in coverage[i:i+20]]))

    return coverage_

def plot_snps_counts(chrom, index, ref_start_values, df_snps_in_csv, html_graphs, arguments, df_chrom):
    fig = go.Figure()

    snps_het_counts, snps_homo_counts, snps_het_pos, snps_homo_pos = get_snps_counts(df_snps_in_csv, chrom, ref_start_values, arguments['bin_size_snps'])

    #BED output SNPs counts per bin
    coverage = coverage_bins(df_chrom)
    save_snps_counts_per_bin(snps_het_counts, snps_homo_counts, ref_start_values, chrom, coverage, arguments, index)

    add_scatter_trace_coverage(fig, snps_het_pos, snps_het_counts, name='Het SNPs counts', text=None, yaxis=None,
                               opacity=1, color='#E3B448')
    add_scatter_trace_coverage(fig, snps_homo_pos, snps_homo_counts, name='Homo SNPs counts', text=None, yaxis=None,
                               opacity=1, color='#3A6B35')

    plots_add_markers_lines(fig)
    plots_layout_settings(fig, chrom, arguments, snps_homo_pos[-1:][0], arguments['cut_threshold_snps_counts'])
    fig.update_layout(yaxis_title="<b>SNPs counts</b> (per bin)", )

    print_chromosome_html(fig, chrom + '_snps_counts', html_graphs, arguments['out_dir_plots']+'/variation_plots/')
    html_graphs.write(
        "  <object data=\"" + chrom + '_snps_counts' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")