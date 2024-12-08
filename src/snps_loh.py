import statistics

import plotly.graph_objects as go
import os
import logging
import pandas as pd
import numpy as np

from vcf_processing import get_snps_frquncies, het_homo_snps_gts, vcf_parse_to_csv_for_het_phased_snps_phasesets, snps_mean, cpd_mean, get_snp_segments, get_snps_frquncies_coverage, vcf_parse_to_csv_for_snps, get_snps_counts, get_snps_counts_cn_regions, get_snps_frquncies_genome
from utils import csv_df_chromosomes_sorter, get_vafs_from_tumor_phased_vcf
from extras import get_contigs_list
from plots import add_scatter_trace_coverage, print_chromosome_html, plots_add_markers_lines, plots_layout_settings,\
whole_genome_combined_df, copy_number_plots_per_chromosome, print_genome_pdf, add_annotation
from utils import detect_alter_loh_regions, loh_regions_events, write_segments_coverage, write_header_comments, get_vafs_from_normal_phased_vcf, get_chromosomes_regions

def plot_snps_ratios_genome(args, df_snps_in_csv, df_loh_regions):
    if not os.path.isdir(args.out_dir_plots + '/snps_loh_plots'):
        os.makedirs(args.out_dir_plots + '/snps_loh_plots')
    chroms = get_contigs_list(args.contigs)
    df_snps_ratios = []
    offset = 0
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
                    df_snps_in_csv, chrom, ref_start_values, args)

                chr_list = [chrom for ch in range(len(ref_start_values_updated))]

                df_snps_ratios_chrom = pd.DataFrame(list(zip(chr_list, ref_start_values_updated, snps_het_counts, snps_homo_counts)),
                    columns=['chr', 'start', 'hets_ratios', 'homos_ratios'])
                regions = get_chromosomes_regions(args)
                df_snps_ratios_chrom['start_overall'] = df_snps_ratios_chrom['start'].apply(lambda x: x + offset)
                offset += regions[index]

                df_snps_ratios.append(df_snps_ratios_chrom)

    loh_plots_genome(df_snps_in_csv, pd.concat(df_snps_ratios), args, df_loh_regions)

def snps_df_loh(args, thread_pool, df_coverage_data):
    df_coverage = df_coverage_data.copy()
    chroms = get_contigs_list(args.contigs)
    if args.tumor_vcf:
        output_phasesets_file_path = vcf_parse_to_csv_for_snps(args.tumor_vcf, args)
        df_snps_in_csv = csv_df_chromosomes_sorter(output_phasesets_file_path, ['chr', 'pos', 'qual', 'gt', 'dp', 'vaf'])
        df_snps_in_csv = get_vafs_from_tumor_phased_vcf(df_snps_in_csv, df_coverage, chroms, args)
    else:
        get_snp_segments(args, args.target_bam[0], thread_pool)
        if args.phaseblock_flipping_disable and args.quick_start:
            df_snps_frequencies = csv_df_chromosomes_sorter(args.out_dir_plots+'/data/snps_frequencies.csv', ['chr', 'pos', 'freq_value_a', 'hp_a', 'freq_value_b', 'hp_b'])
        else:
            df_snps_frequencies = csv_df_chromosomes_sorter(args.out_dir_plots+'/data_phasing/snps_frequencies.csv', ['chr', 'pos', 'freq_value_a', 'hp_a', 'freq_value_b', 'hp_b'])
        #df_snps_frequencies = df_snps_frequencies.drop(df_snps_frequencies[(df_snps_frequencies.chr == "chrY")].index)
        df_snps_in_csv = get_vafs_from_normal_phased_vcf(df_snps_frequencies, df_coverage, chroms, args)

    return df_snps_in_csv

def plot_snps_frequencies_without_phasing(args, df, df_segs_hp1_w, df_segs_hp2_w, centers, integer_fractional_means, df_snps_in_csv):
    if not os.path.isdir(args.out_dir_plots + '/snps_loh_plots'):
        os.makedirs(args.out_dir_plots + '/snps_loh_plots')

    filename = f"{os.path.join(args.out_dir_plots, 'snps_loh_plots/SNPs_LOH_INDEX.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")

    chroms = get_contigs_list(args.contigs)

    loh_regions_events_all = []
    df_snps_counts_per_cn_region_all = []
    df_snps_ratios = []
    loh_regions_all = []

    for index, chrom in enumerate(chroms):
        if chrom in chroms:
            chrs_list = []
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

            chrs_list.extend([chrom for ch in range(len(loh_region_starts))])
            loh_regions_all.append(pd.DataFrame(list(zip(chrs_list, loh_region_starts, loh_region_ends)), columns=['chr', 'start', 'end']))

            if args.without_phasing:
                copy_number_plots_per_chromosome(centers, integer_fractional_means, ref_start_values, hp1, df_segs_hp1_, hp2, df_segs_hp2_, args, chrom, html_graphs, loh_region_starts, loh_region_ends, None)

    write_snps_counts_per_cn_region(pd.concat(df_snps_counts_per_cn_region_all), args)
    write_header_comments('chr\tstart\tend\n', '#chr: chromosome number\n#start: start address for LOH region\n#end: end address for LOH region\n', args.genome_name + '_loh_segments.bed', args)
    write_segments_coverage(loh_regions_events_all, args.genome_name + '_loh_segments.bed', args)

        #loh_plots_genome(df_snps_in_csv, pd.concat(df_snps_ratios), args)

    html_graphs.write("</body></html>")

    return pd.concat(loh_regions_all)


def variation_plots(args, df, df_segs_hp1_w, df_segs_hp2_w, centers, integer_fractional_means, df_snps_in_csv, df_loh_regions, p_value):

    if not os.path.isdir(args.out_dir_plots + '/' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/variation_plots'):
        os.makedirs(args.out_dir_plots +'/'+ str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/variation_plots')

    filename = f"{os.path.join(args.out_dir_plots, str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value), 'variation_plots/CN_VARIATION_INDEX.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")

    chroms = get_contigs_list(args.contigs)

    loh_regions_events_all = []
    df_snps_counts_per_cn_region_all = []
    df_snps_ratios = []

    for index, chrom in enumerate(chroms):
        if chrom in chroms:
            df_chrom = df[df['chr'] == chrom]
            df_segs_hp1_ = df_segs_hp1_w[df_segs_hp1_w['chromosome'] == chrom]
            df_segs_hp2_ = df_segs_hp2_w[df_segs_hp2_w['chromosome'] == chrom]
            df_loh_region = df_loh_regions[df_loh_regions['chr'] == chrom]
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

            copy_number_plots_per_chromosome(centers, integer_fractional_means, ref_start_values, hp1, df_segs_hp1_, hp2, df_segs_hp2_, args, chrom, html_graphs, df_loh_region.start.values.tolist(), df_loh_region.end.values.tolist(), p_value)

    html_graphs.write("</body></html>")

def plot_snps_frequencies(args, df, df_snps_in_csv):
    if not os.path.isdir(args.out_dir_plots + '/snps_loh_plots'):
        os.makedirs(args.out_dir_plots + '/snps_loh_plots')

    filename = f"{os.path.join(args.out_dir_plots, 'snps_loh_plots/SNPs_LOH_INDEX.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")

    chroms = get_contigs_list(args.contigs)

    loh_regions_events_all = []
    df_snps_counts_per_cn_region_all = []
    df_snps_ratios = []
    loh_regions_all = []

    for index, chrom in enumerate(chroms):
        if chrom in chroms:
            chrs_list = []
            df_chrom = df[df['chr'] == chrom]
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

            centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends, df_snps_ratios_chrom = plot_snps(chrom, index, df_snps_in_csv, html_graphs, args, df_chrom)
            df_snps_ratios.append(df_snps_ratios_chrom)
            #hp1, hp2, hp3 = detect_alter_loh_regions(args, 'centromere/no-coverage', chrom, ref_end_values, hp1, hp2, hp3, centromere_region_starts, centromere_region_ends)
            if args.without_phasing:
                hp1, hp2, hp3, loh_region_starts, loh_region_ends = detect_alter_loh_regions(args, 'loss-of-heterozygosity', chrom, ref_end_values, coverage, coverage, coverage, loh_region_starts, loh_region_ends, False)
            else:
                hp1, hp2, hp3, loh_region_starts, loh_region_ends = detect_alter_loh_regions(args, 'loss-of-heterozygosity', chrom, ref_end_values, hp1, hp2, hp3, loh_region_starts, loh_region_ends, False)
            loh_regions_events_all.extend(loh_regions_events(chrom, loh_region_starts, loh_region_ends, args))

            chrs_list.extend([chrom for ch in range(len(loh_region_starts))])
            loh_regions_all.append(pd.DataFrame(list(zip(chrs_list, loh_region_starts, loh_region_ends)), columns=['chr', 'start', 'end']))

    # write_snps_counts_per_cn_region(pd.concat(df_snps_counts_per_cn_region_all), args)
    write_header_comments('chr\tstart\tend\n', '#chr: chromosome number\n#start: start address for LOH region\n#end: end address for LOH region\n', args.genome_name + '_loh_segments.bed', args)
    write_segments_coverage(loh_regions_events_all, args.genome_name + '_loh_segments.bed', args)

        #loh_plots_genome(df_snps_in_csv, pd.concat(df_snps_ratios), args)

    html_graphs.write("</body></html>")

    return pd.concat(loh_regions_all)


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
            df_snps_in_csv, chrom, ref_start_values, args)

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

        print_chromosome_html(fig, chrom + '_snps', html_graphs, args.out_dir_plots + '/snps_loh_plots/')
        html_graphs.write("  <object data=\"" + chrom + '_snps' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

        return centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends, df_snps_ratios_chrom
    else:
        return [], [], [], []
def save_snps_counts_per_bin(hets,homos,ref_start, chrom, coverage, args, index):
    import pandas as pd
    fp = open(args.out_dir_plots + '/bed_output/' + args.genome_name + '_snps_counts.bed', 'a')

    if index == 0:
        header = '#chr: chromosome number\n#start: start address for SNPs bin\n#hets_count: number of hetrozygous SNPs\n#homo_counts: number of homozygous SNPs\n#coverage: mean coverage in this bin\nchr\tstart\thets_count\thomos_count\tcoverage\n'
        fp.write(header)
    chr_list = [chrom for ch in range(len(ref_start))]
    snps_counts_df = pd.DataFrame(list(zip(chr_list, ref_start, hets, homos, coverage)), columns=['chr','start','hets_count','homos_count','coverage'])
    snps_counts_df.to_csv(fp, sep='\t', index=False, mode='a', header=False)

def snps_counts_per_cn_region(hets, homos, ref_start, ref_end, chrom):
    import pandas as pd
    chr_list = [chrom for ch in range(len(ref_start))]
    snps_counts_df = pd.DataFrame(list(zip(chr_list, ref_start, ref_end, hets, homos)), columns=['chr', 'start', 'end', 'hets_count', 'homos_count'])
    return snps_counts_df

def write_snps_counts_per_cn_region(df, args):
    fp = open(args.out_dir_plots + '/bed_output/' + args.genome_name + '_snps_counts_cn_regions.bed', 'a')
    fp.write('#chr: chromosome number\n')
    fp.write('#start: start address for CN segment\n')
    fp.write('#end: end address for CN segment\n')
    fp.write('#hets_count: hetrozygous SNPs count in this region\n')
    fp.write('#homos_count: homozygous SNPs count in this region\n')
    fp.write('#chr\tstart\tend\thets_count\thomos_count\n')

    df.to_csv(fp, sep='\t', index=False, mode='a', header=False)

def write_loh_regions(df, path, args, p_value_confidence):
    fp = open(args.out_dir_plots + '/' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) +'_'+ str(p_value_confidence) + '/bed_output/' + path, 'a')
    fp.write('#chr: chromosome number\n')
    fp.write('#start: start address for loh segment\n')
    fp.write('#end: end address for loh segment\n')
    fp.write('#chr\tstart\tend\n')

    df.to_csv(fp, sep='\t', index=False, mode='a', header=False)

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

    print_chromosome_html(fig, chrom + '_snps_counts', html_graphs, args.out_dir_plots +'/snps_loh_plots/')
    html_graphs.write(
        "  <object data=\"" + chrom + '_snps_counts' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

def loh_plots_genome(df_snps_in_csv, df_snps_ratios, args, df_loh_regions):
    # TODO genome-wide plots

    lengths = []
    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        df_snps_ratios_chrom = df_snps_ratios[df_snps_ratios['chr'] == chrom]
        lengths.append(len(df_snps_ratios_chrom))

    indices = np.arange(0, len(df_snps_ratios)*args.bin_size_snps, args.bin_size_snps, dtype=int)

    fig = go.Figure()

    custom_text_data = df_snps_ratios.start.values.tolist()

    add_scatter_trace_coverage(fig, df_snps_ratios.start_overall.values.tolist(), df_snps_ratios.hets_ratios.values.tolist(), name='Het SNPs Ratios', text=custom_text_data,
                               yaxis=None, opacity=0.7, color='#E3B448', mul_cols=False)
    add_scatter_trace_coverage(fig, df_snps_ratios.start_overall.values.tolist(), df_snps_ratios.homos_ratios.values.tolist(), name='Homo SNPs Ratios', text=custom_text_data,
                               yaxis=None, opacity=0.7, color='#3A6B35', mul_cols=False)

    ###################################
    #snps_het, snps_homo, indices_het_freqs, indices_homo_freqs = get_snps_frquncies_genome(df_snps_in_csv)

    #add_scatter_trace_coverage(fig, indices_het_freqs, snps_het, name='Het SNPs Freqs', text=None, yaxis=None,
    #                           opacity=0.7, color='#E3B448', visibility='legendonly')
    #add_scatter_trace_coverage(fig, indices_homo_freqs, snps_homo, name='Homo SNPs Freqs', text=None, yaxis=None,
    #                           opacity=0.7, color='#3A6B35', visibility='legendonly')
    regions = get_chromosomes_regions(args)
    current = 0
    label_pos = []
    offset_chroms = 0
    offset_chroms_1 = 0
    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        if not df_loh_regions.empty:
            df_loh_region = df_loh_regions[df_loh_regions['chr'] == chrom]
        current += lengths[index]
        label_pos.append(round(offset_chroms+regions[index]//2)) #y0=-10, y1=args.cut_threshold
        fig.add_vline(x=offset_chroms,  line_width=1, line_dash="solid", line_color="#D7DBDD")
        offset_chroms_1 += regions[index]

        loh_starts = df_loh_region.start.values.tolist()
        loh_ends  = df_loh_region.end.values.tolist()

        # if len(loh_starts):
        #     for i in range(len(df_loh_region.start.values.tolist())):
        #         fig.add_vrect(x0=offset_chroms+loh_starts[i], x1=offset_chroms+loh_ends[i], fillcolor="#2980b9", opacity=0.4, line_width=0, row=2, col=1,)

        if index % 2 == 0:
            fig.add_vrect(x0=offset_chroms, x1=offset_chroms_1, fillcolor="#E5E7E9", opacity=0.9, layer="below", line_width=0, )

        offset_chroms += regions[index]

    fig.add_hline(y=1-args.hets_ratio,  line_width=1.5, line_dash="solid", line_color="#7a7b7b")
    fig.add_hline(y=args.hets_ratio,  line_width=1.5, line_dash="solid", line_color="#7a7b7b")

    fig.update_layout(
        xaxis=dict(
            tickangle=90,
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        ),
        font=dict(size=18, color="black"))

    plots_layout_settings(fig, 'Genome', args, df_snps_ratios.start_overall.values.tolist()[-1], args.cut_threshold)

    ax = 20
    ay = -30
    add_annotation(fig, 960000000 + 250000000 + 120000000, 1, ax, ay, "LOH", '#2980b9')

    fig.add_annotation(text="Het SNPs ratio threshold: " + str(args.hets_ratio), x=960000000 + 250000000 + 120000000, y=0.95, showarrow=False)
    fig.update_yaxes(range=[0,1])
    fig.update_yaxes(range=[0,1])
    fig.update_yaxes(title_text="<b>SNPs Frequencies</b> (ratios)")
    fig.update_layout(width=1680, height=600,)
    fig.update_layout(legend=dict(orientation = 'h', xanchor = "center", x = 0.5, y= 1.1))
    if args.pdf_enable:
        print_genome_pdf(fig, args.genome_name+'_loh', args.out_dir_plots + '/snps_loh_plots')

    fig.write_html(args.out_dir_plots + '/snps_loh_plots' +'/'+ args.genome_name + "_genome_snps_ratio_loh.html")