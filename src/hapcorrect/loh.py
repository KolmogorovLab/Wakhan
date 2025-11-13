import statistics

import plotly.graph_objects as go
import os
import logging
logger = logging.getLogger()

from src.file_tools.process_vcf import get_snps_frquncies, vcf_parse_to_csv_for_het_phased_snps_phasesets, get_snps_frquncies_coverage, vcf_parse_to_csv_for_snps, get_snps_counts
from src.hapcorrect.utils import csv_df_chromosomes_sorter, loh_regions_phasesets, detect_alter_loh_regions
from src.breakpoints import get_contigs_list
from src.hapcorrect.plots import add_scatter_trace_coverage, print_chromosome_html, plots_add_markers_lines, plots_layout_settings


def detect_loh_centromere_regions(csv_df_coverage_tumor_chrom, chrom, args, centromere_region_starts, centromere_region_ends, 
                                  loh_region_starts, loh_region_ends, ref_start_values, ref_end_values, haplotype_1_values, haplotype_2_values, 
                                  unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, 
                                  ref_start_values_phasesets, ref_end_values_phasesets):
    hp = []
    if args.without_phasing:
        values = haplotype_1_values
        if centromere_region_starts:
            _, _, _, centromere_region_starts, centromere_region_ends, hp = \
                    detect_alter_loh_regions(None, args, 'centromere/no-coverage', chrom, ref_end_values, values, values, values,
                                             centromere_region_starts, centromere_region_ends, True)
        if loh_region_starts:
            _, _, _, loh_region_starts, loh_region_ends, hp = \
                    detect_alter_loh_regions(None, args, 'loss-of-heterozygosity', chrom, ref_end_values, values, values, values,
                                             loh_region_starts, loh_region_ends, True)
    else:
        # if centromere_region_starts:
        #     haplotype_1_values, haplotype_2_values, unphased_reads_values, centromere_region_starts, centromere_region_ends = detect_alter_loh_regions(args, 'centromere/no-coverage', chrom, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values, centromere_region_starts, centromere_region_ends, True)
        #     haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = loh_regions_phasesets(centromere_region_starts, centromere_region_ends, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

        if loh_region_starts:
            haplotype_1_values, haplotype_2_values, unphased_reads_values, loh_region_starts, loh_region_ends, hp = \
                    detect_alter_loh_regions(csv_df_coverage_tumor_chrom, args, 'loss-of-heterozygosity', chrom, ref_end_values,
                                             haplotype_1_values, haplotype_2_values, unphased_reads_values, loh_region_starts, loh_region_ends, True)
            haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = \
                    loh_regions_phasesets(haplotype_1_values, haplotype_2_values, loh_region_starts, loh_region_ends,
                                          haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets,
                                          ref_end_values_phasesets, args)

        return (haplotype_1_values, haplotype_2_values, unphased_reads_values, haplotype_1_values_phasesets,
                haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, loh_region_starts, loh_region_ends, hp)


def plot_snps(args, df_snps_in_csv):
    filename = f"{os.path.join(args.out_dir_plots, 'snps_loh_plots',  'SNPs_DEBUG.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")

    chroms = get_contigs_list(args.contigs)

    for index, chrom in enumerate(chroms):
        plot_snps_freqs_ratios_counts(chrom, index, df_snps_in_csv, html_graphs, args)

    html_graphs.write("</body></html>")

def plot_snps_freqs_ratios_counts(chrom, index, df_snps_in_csv, html_graphs, args):
    logger.info('SNPs frequencies plots generation for ' + chrom)
    fig = go.Figure()

    snps_het, snps_homo, snps_het_pos, snps_homo_pos = get_snps_frquncies(df_snps_in_csv, chrom)

    add_scatter_trace_coverage(fig, snps_het_pos, snps_het, name='Het SNPs Freqs', text=None, yaxis=None,
                               opacity=0.7, color='#E3B448')
    add_scatter_trace_coverage(fig, snps_homo_pos, snps_homo, name='Homo SNPs Freqs', text=None, yaxis=None,
                               opacity=0.7, color='#3A6B35')

    ref_start_values = [i for i in range(0, snps_homo_pos[-1:][0], args.bin_size_snps)]
    ref_start_values_updated, snps_het_counts, snps_homo_counts, centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends = get_snps_frquncies_coverage(
        df_snps_in_csv, chrom, ref_start_values, args.bin_size_snps, args.hets_ratio, args.hets_smooth_window, args)

    add_scatter_trace_coverage(fig, ref_start_values_updated, snps_het_counts, name='Het SNPs Ratios', text=None,
                               yaxis=None,
                               opacity=0.7, color='#E3B448')
    add_scatter_trace_coverage(fig, ref_start_values_updated, snps_homo_counts, name='Homo SNPs Ratios', text=None,
                               yaxis=None,
                               opacity=0.7, color='#3A6B35')

    plots_add_markers_lines(fig)
    plots_layout_settings(fig, chrom, args, snps_homo_pos[-1:][0], 0)
    fig.update_layout(yaxis_title="<b>SNPs Frequencies</b> (ratios)",)

    #plot_snps_counts(chrom, index, ref_start_values, df_snps_in_csv, html_graphs, args)

    print_chromosome_html(fig, chrom + '_snps', html_graphs, args.out_dir_plots+ '/snps_loh_plots')
    html_graphs.write(
        "  <object data=\"" + chrom + '_snps' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

    #return centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends
def plot_snps_counts(chrom, index, ref_start_values, df_snps_in_csv, html_graphs, args):
    fig = go.Figure()

    snps_het_counts, snps_homo_counts, snps_het_pos, snps_homo_pos = get_snps_counts(df_snps_in_csv, chrom, ref_start_values, args.bin_size_snps)

    add_scatter_trace_coverage(fig, snps_het_pos, snps_het_counts, name='Het SNPs counts', text=None, yaxis=None,
                               opacity=1, color='#E3B448')
    add_scatter_trace_coverage(fig, snps_homo_pos, snps_homo_counts, name='Homo SNPs counts', text=None, yaxis=None,
                               opacity=1, color='#3A6B35')

    plots_add_markers_lines(fig)
    plots_layout_settings(fig, chrom, args, snps_homo_pos[-1:][0], 250)
    fig.update_layout(yaxis_title="<b>SNPs counts</b> (per bin)", )

    print_chromosome_html(fig, chrom + '_snps_counts', html_graphs, args.out_dir_plots)
    html_graphs.write(
        "  <object data=\"" + chrom + '_snps_counts' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")
