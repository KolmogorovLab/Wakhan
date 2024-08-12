import statistics

import plotly.graph_objects as go
import os
import logging
import pandas as pd
import numpy as np

from vcf_processing import get_snps_frquncies, het_homo_snps_gts, vcf_parse_to_csv_for_het_phased_snps_phasesets, snps_mean, cpd_mean, get_snp_segments, get_snps_frquncies_coverage, vcf_parse_to_csv_for_snps, get_snps_counts, get_snps_counts_cn_regions, get_snps_frquncies_genome
from utils import csv_df_chromosomes_sorter
from extras import get_contigs_list
from plots import add_scatter_trace_coverage, print_chromosome_html, plots_add_markers_lines, plots_layout_settings,\
whole_genome_combined_df, copy_number_plots_per_chromosome, print_genome_pdf
from utils import detect_alter_loh_regions, loh_regions_events, write_segments_coverage, write_header_comments

def plot_snps_ratios_genome(args):
    if args.tumor_vcf:
        output_phasesets_file_path = vcf_parse_to_csv_for_snps(args.tumor_vcf, args)
        df_snps_in_csv = csv_df_chromosomes_sorter(output_phasesets_file_path, ['chr', 'pos', 'qual', 'gt', 'dp', 'vaf'])

    chroms = get_contigs_list(args.contigs)
    df_snps_ratios = []

    for index, chrom in enumerate(chroms):
        if chrom in chroms:
            logging.info('Collecting SNPs for chromosome ' + chrom)
            snps_het, snps_homo, snps_het_pos, snps_homo_pos = get_snps_frquncies(df_snps_in_csv, chrom)
            if snps_het or snps_homo:
                if snps_homo_pos and snps_het_pos:
                    last_snp_pos = snps_homo_pos[-1:][0] if snps_homo_pos[-1:][0] > snps_het_pos[-1:][0] else \
                    snps_het_pos[-1:][0]
                elif snps_het_pos:
                    last_snp_pos = snps_het_pos[-1:][0]
                elif snps_homo_pos:
                    last_snp_pos = snps_homo_pos[-1:][0]

                ref_start_values = [i for i in range(0, last_snp_pos, args.bin_size_snps)]
                ref_start_values_updated, snps_het_counts, snps_homo_counts, centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends = get_snps_frquncies_coverage(
                    df_snps_in_csv, chrom, ref_start_values, args.bin_size_snps)

                chr_list = [chrom for ch in range(len(ref_start_values_updated))]
                df_snps_ratios_chrom = pd.DataFrame(list(zip(chr_list, ref_start_values_updated, snps_het_counts, snps_homo_counts)),
                    columns=['chr', 'start', 'hets_ratios', 'homos_ratios'])

                df_snps_ratios.append(df_snps_ratios_chrom)

    loh_plots_genome(df_snps_in_csv, pd.concat(df_snps_ratios), args)
def plot_snps_frequencies(args, df, df_segs_hp1_w, df_segs_hp2_w, centers, integer_fractional_means):
    if not os.path.isdir(args.out_dir_plots+'/variation_plots'):
        os.mkdir(args.out_dir_plots+'/variation_plots')
    filename = f"{os.path.join(args.out_dir_plots, 'variation_plots/SNPs_LOH_CNs_INDEX.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")

    if args.tumor_vcf:
        output_phasesets_file_path = vcf_parse_to_csv_for_snps(args.tumor_vcf, args)
        df_snps_in_csv = csv_df_chromosomes_sorter(output_phasesets_file_path, ['chr', 'pos', 'qual', 'gt', 'dp', 'vaf'])

    chroms = get_contigs_list(args.contigs)
    loh_regions_events_all = []
    df_snps_counts_per_cn_region_all = []
    df_snps_ratios = []

    for index, chrom in enumerate(chroms):
        if chrom in chroms:
            df_chrom = df[df['chr'] == chrom]
            df_segs_hp1_ = df_segs_hp1_w[df_segs_hp1_w['chromosome'] == chrom]
            df_segs_hp2_ = df_segs_hp2_w[df_segs_hp2_w['chromosome'] == chrom]
            ref_start_values = df_chrom.start.values.tolist()
            ref_end_values = df_chrom.end.values.tolist()
            chr = range(len(ref_start_values))

            if args.without_phasing:
                coverage = df_chrom.coverage.values.tolist()
                df_snps_freqs_chr = whole_genome_combined_df(args, chrom, chr, ref_start_values, ref_end_values, coverage, coverage, coverage)

            else:
                hp1 = df_chrom.hp1.values.tolist()
                hp2 = df_chrom.hp2.values.tolist()
                hp3 = df_chrom.hp3.values.tolist()
                df_snps_freqs_chr = whole_genome_combined_df(args, chrom, chr, ref_start_values, ref_end_values, hp1, hp2, hp3)

            if args.tumor_vcf:
                snps_het_counts, snps_homo_counts, start_values, end_values = get_snps_counts_cn_regions(df_snps_in_csv, chrom, df_segs_hp1_.start.values.tolist(), df_segs_hp1_.end.values.tolist())
                df_snps_counts_per_cn_region_all.append(snps_counts_per_cn_region(snps_het_counts, snps_homo_counts, start_values, end_values, chrom))

                centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends, df_snps_ratios_chrom = plot_snps(chrom, index, df_snps_in_csv, html_graphs, args, df_chrom)
                df_snps_ratios.append(df_snps_ratios_chrom)
                #hp1, hp2, hp3 = detect_alter_loh_regions(args, 'centromere/no-coverage', chrom, ref_end_values, hp1, hp2, hp3, centromere_region_starts, centromere_region_ends)
                if args.without_phasing:
                    hp1, hp2, hp3, loh_region_starts, loh_region_ends = detect_alter_loh_regions(args, 'loss-of-heterozygosity', chrom, ref_end_values, coverage, coverage, coverage, loh_region_starts, loh_region_ends, False)
                else:
                    hp1, hp2, hp3, loh_region_starts, loh_region_ends = detect_alter_loh_regions(args, 'loss-of-heterozygosity', chrom, ref_end_values, hp1, hp2, hp3, loh_region_starts, loh_region_ends, False)
                loh_regions_events_all.extend(loh_regions_events(chrom, loh_region_starts, loh_region_ends, args))
            else:
                loh_region_starts = []
                loh_region_ends = []

            copy_number_plots_per_chromosome(centers, integer_fractional_means, ref_start_values, hp1, df_segs_hp1_, hp2, df_segs_hp2_, args, chrom, html_graphs, loh_region_starts, loh_region_ends)

    if args.tumor_vcf:
        write_snps_counts_per_cn_region(pd.concat(df_snps_counts_per_cn_region_all), args)
        write_header_comments('chr\tstart\tend\n', '#chr: chromosome number\n#start: start address for LOH region\n#end: end address for LOH region\n', args.genome_name + '_loh_segments.bed', args)
        write_segments_coverage(loh_regions_events_all, args.genome_name + '_loh_segments.bed', args)

        #loh_plots_genome(df_snps_in_csv, pd.concat(df_snps_ratios), args)

    html_graphs.write("</body></html>")


def plot_snps(chrom, index, df_snps_in_csv, html_graphs, args, df_chrom):
    logging.info('SNPs frequencies plots generation for ' + chrom)
    fig = go.Figure()

    snps_het, snps_homo, snps_het_pos, snps_homo_pos = get_snps_frquncies(df_snps_in_csv, chrom)

    if snps_het or snps_homo:
        #for debug purpose, all SNPs Freqs
        add_scatter_trace_coverage(fig, snps_het_pos, snps_het, name='Het SNPs Freqs', text=None, yaxis=None,
                                   opacity=0.7, color='#E3B448', visibility='legendonly')
        add_scatter_trace_coverage(fig, snps_homo_pos, snps_homo, name='Homo SNPs Freqs', text=None, yaxis=None,
                                   opacity=0.7, color='#3A6B35', visibility='legendonly')
        if snps_homo_pos and snps_het_pos:
            last_snp_pos = snps_homo_pos[-1:][0] if snps_homo_pos[-1:][0] > snps_het_pos[-1:][0] else snps_het_pos[-1:][0]
        elif snps_het_pos:
            last_snp_pos = snps_het_pos[-1:][0]
        elif snps_homo_pos:
            last_snp_pos = snps_homo_pos[-1:][0]

        ref_start_values = [i for i in range(0, last_snp_pos, args.bin_size_snps)]
        ref_start_values_updated, snps_het_counts, snps_homo_counts, centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends = get_snps_frquncies_coverage(
            df_snps_in_csv, chrom, ref_start_values, args.bin_size_snps)

        chr_list = [df_chrom['chr'].iloc[0] for ch in range(len(ref_start_values_updated))]
        df_snps_ratios_chrom = pd.DataFrame(list(zip(chr_list, ref_start_values_updated, snps_het_counts, snps_homo_counts)),
                                    columns=['chr', 'start', 'hets_ratios', 'homos_ratios'])

        add_scatter_trace_coverage(fig, ref_start_values_updated, snps_het_counts, name='Het SNPs Ratios', text=None,
                                   yaxis=None,
                                   opacity=0.7, color='#E3B448')
        add_scatter_trace_coverage(fig, ref_start_values_updated, snps_homo_counts, name='Homo SNPs Ratios', text=None,
                                   yaxis=None,
                                   opacity=0.7, color='#3A6B35')

        plots_add_markers_lines(fig)
        plots_layout_settings(fig, chrom, args, snps_homo_pos[-1:][0], 0)
        fig.update_layout(yaxis_title="<b>SNPs Frequencies</b> (ratios)",)

        plot_snps_counts(chrom, index, ref_start_values, df_snps_in_csv, html_graphs, args, df_chrom)

        print_chromosome_html(fig, chrom + '_snps', html_graphs, args.out_dir_plots+'/variation_plots/')
        html_graphs.write("  <object data=\"" + chrom + '_snps' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

        return centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends, df_snps_ratios_chrom
    else:
        return [], [], [], []
def save_snps_counts_per_bin(hets,homos,ref_start, chrom, coverage, args, index):
    import pandas as pd
    fp = open(args.out_dir_plots+'/bed_output/' + args.genome_name + '_snps_counts.bed', 'a')

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

def write_snps_counts_per_cn_region(df, args):
    fp = open(args.out_dir_plots+'/bed_output/' + args.genome_name + '_snps_counts_cn_regions.bed', 'a')
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

def plot_snps_counts(chrom, index, ref_start_values, df_snps_in_csv, html_graphs, args, df_chrom):
    fig = go.Figure()

    snps_het_counts, snps_homo_counts, snps_het_pos, snps_homo_pos = get_snps_counts(df_snps_in_csv, chrom, ref_start_values, args.bin_size_snps)

    if args.without_phasing:
        coverage = coverage_bins(df_chrom)
        save_snps_counts_per_bin(snps_het_counts, snps_homo_counts, ref_start_values, chrom, coverage, args, index)

    #BED output SNPs counts per bin
    add_scatter_trace_coverage(fig, snps_het_pos, snps_het_counts, name='Het SNPs counts', text=None, yaxis=None,
                               opacity=1, color='#E3B448')
    add_scatter_trace_coverage(fig, snps_homo_pos, snps_homo_counts, name='Homo SNPs counts', text=None, yaxis=None,
                               opacity=1, color='#3A6B35')

    plots_add_markers_lines(fig)
    plots_layout_settings(fig, chrom, args, snps_homo_pos[-1:][0], args.cut_threshold_snps_counts)
    fig.update_layout(yaxis_title="<b>SNPs counts</b> (per bin)", )

    print_chromosome_html(fig, chrom + '_snps_counts', html_graphs, args.out_dir_plots+'/variation_plots/')
    html_graphs.write(
        "  <object data=\"" + chrom + '_snps_counts' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

def loh_plots_genome(df_snps_in_csv, df_snps_ratios, args):
    # TODO genome-wide plots

    lengths = []
    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        df_snps_ratios_chrom = df_snps_ratios[df_snps_ratios['chr'] == chrom]
        lengths.append(len(df_snps_ratios_chrom))

    indices = np.arange(0, len(df_snps_ratios)*args.bin_size_snps, args.bin_size_snps, dtype=int)

    fig = go.Figure()

    custom_text_data = df_snps_ratios.start.values.tolist()

    add_scatter_trace_coverage(fig, indices, df_snps_ratios.hets_ratios.values.tolist(), name='Het SNPs Ratios', text=custom_text_data,
                               yaxis=None, opacity=0.7, color='#E3B448', mul_cols=False)
    add_scatter_trace_coverage(fig, indices, df_snps_ratios.homos_ratios.values.tolist(), name='Homo SNPs Ratios', text=custom_text_data,
                               yaxis=None, opacity=0.7, color='#3A6B35', mul_cols=False)

    ###################################
    #snps_het, snps_homo, indices_het_freqs, indices_homo_freqs = get_snps_frquncies_genome(df_snps_in_csv)

    #add_scatter_trace_coverage(fig, indices_het_freqs, snps_het, name='Het SNPs Freqs', text=None, yaxis=None,
    #                           opacity=0.7, color='#E3B448', visibility='legendonly')
    #add_scatter_trace_coverage(fig, indices_homo_freqs, snps_homo, name='Homo SNPs Freqs', text=None, yaxis=None,
    #                           opacity=0.7, color='#3A6B35', visibility='legendonly')

    current = 0
    label_pos = []
    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        current += lengths[index]
        label_pos.append(round(current - (lengths[index] / 2))*args.bin_size) #y0=-10, y1=args.cut_threshold
        fig.add_vline(x=current*args.bin_size,  line_width=1, line_dash="solid", line_color="#D7DBDD")

        if index == 0:
            start_chrom = 0
        else:
            start_chrom += lengths[index-1]
        if index % 2 == 0:
            fig.add_vrect(x0=start_chrom*args.bin_size, x1=current*args.bin_size, fillcolor="#E5E7E9", opacity=0.9, layer="below", line_width=0, )

    fig.add_hline(y=0.7,  line_width=1.5, line_dash="dashdot", line_color="#3A6B35")
    fig.add_hline(y=0.3,  line_width=1.5, line_dash="dashdot", line_color="#E3B448")

    fig.update_layout(
        xaxis=dict(
            tickangle=90,
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        ),
        font=dict(size=18, color="black"))

    plots_layout_settings(fig, 'Genome', args, indices[-1:][0], args.cut_threshold)

    fig.update_yaxes(range=[0,1])
    fig.update_yaxes(title_text="<b>SNPs Frequencies</b> (ratios)")
    fig.update_layout(width=1680, height=600,)
    fig.update_layout(legend=dict(orientation = 'h', xanchor = "center", x = 0.5, y= 1.1))
    if args.pdf_enable:
        print_genome_pdf(fig, args.genome_name+'_loh', args.out_dir_plots)

    fig.write_html(args.out_dir_plots +'/'+ args.genome_name + "_genome_loh.html")