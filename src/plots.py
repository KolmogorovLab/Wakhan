import numpy
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly
import pandas as pd
import numpy as np
import csv
import os
import math
import statistics
import itertools
import logging
from cnvlib import cluster
from phasing_correction import phaseblock_flipping, phase_correction_centers, contiguous_phaseblocks, detect_centromeres, flip_phaseblocks_contigous, remove_overlaping_contiguous
from smoothing import smoothing
from vcf_processing import get_snps_frquncies_coverage, get_snps_frquncies, het_homo_snps_gts, vcf_parse_to_csv_for_het_phased_snps_phasesets, snps_mean, cpd_mean, get_snp_segments, vcf_parse_to_csv_for_snps
from utils import csv_df_chromosomes_sorter, get_breakpoints, flatten, get_snps_frquncies_coverage_from_bam, apply_copynumbers, detect_alter_loh_regions, is_phasesets_check_simple_heuristics, change_point_detection_means, loh_regions_phasesets
from extras import get_contigs_list

def copy_number_plots_chromosomes(df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, arguments, loh_region_starts, loh_region_ends):
    filename = f"{os.path.join(arguments['out_dir_plots'], 'COPY_NUMBERS.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")

    chroms = get_contigs_list(arguments['contigs'])
    for index, chrom in enumerate(chroms):
        copy_number_plots_per_chromosome(df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, arguments, chrom, html_graphs, loh_region_starts, loh_region_ends)
    html_graphs.write("</body></html>")

def copy_number_plots_per_chromosome(centers, integer_fractional_means, ref_start_values, df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, arguments, chrom, html_graphs, loh_region_starts, loh_region_ends):

    if not chrom == 'chrX' and not chrom == 'chrY':
        logging.info('Plots generation for ' + chrom)
        fig = go.Figure()

        #df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chromosome'] == chrom]
        df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]

        #df_cnr_hp2_chrom = df_cnr_hp2[df_cnr_hp2['chromosome'] == chrom]
        df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]

        haplotype_1_values_cnr = df_cnr_hp1 #df_cnr_hp1_chrom.log2.values.tolist()
        haplotype_2_values_cnr = df_cnr_hp2 #df_cnr_hp2_chrom.log2.values.tolist()
        haplotype_1_start_values_cnr = ref_start_values#df_cnr_hp1_chrom.start.values.tolist()

        haplotype_1_values_copyratios = df_segs_hp1_chrom.state.values.tolist()
        haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
        haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()

        haplotype_2_values_copyratios = df_segs_hp2_chrom.state.values.tolist()
        haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
        haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()

        if arguments['copynumbers_enable']:
            logging.info('copynumbers plots module')
            OFFSET=0
            #haplotype_1_values_cnr = list(np.asarray(haplotype_1_values_cnr) + OFFSET)
            #haplotype_2_values_cnr = list(np.asarray(haplotype_2_values_cnr) - OFFSET)

            if arguments['without_phasing']:
                add_scatter_trace_coverage(fig, haplotype_1_start_values_cnr, haplotype_1_values_cnr, name='Unphased', text=None, yaxis=None, opacity=0.7, color='olive')
                plots_add_markers_lines(fig)
                plot_copynumbers_scatter_lines(arguments, fig,haplotype_1_values_copyratios, haplotype_1_start_values_copyratios, haplotype_1_end_values_copyratios, haplotype_1_values_copyratios, haplotype_1_start_values_copyratios, haplotype_1_end_values_copyratios, 0)
            else:
                haplotype_2_values_cnr = [-x for x in haplotype_2_values_cnr]
                add_scatter_trace_coverage(fig, haplotype_1_start_values_cnr, haplotype_1_values_cnr, name='HP-1', text=None, yaxis=None, opacity=0.7, color='firebrick')
                add_scatter_trace_coverage(fig, haplotype_1_start_values_cnr, haplotype_2_values_cnr, name='HP-2', text=None, yaxis=None, opacity=0.7, color='steelblue')
                plots_add_markers_lines(fig)

                haplotype_2_values_copyratios = [-x for x in haplotype_2_values_copyratios]
                plot_copynumbers_scatter_lines(arguments, fig,haplotype_1_values_copyratios, haplotype_1_start_values_copyratios, haplotype_1_end_values_copyratios, haplotype_2_values_copyratios, haplotype_2_start_values_copyratios, haplotype_2_end_values_copyratios, OFFSET)

            if loh_region_starts:
                for k, (start_loh,end_loh) in enumerate(zip(loh_region_starts, loh_region_ends)):
                    fig.add_vrect(x0=start_loh, x1=end_loh, fillcolor="lightgrey", opacity=0.5, layer="below", line_width=0, )

            plots_layout_settings(fig, chrom, arguments, haplotype_1_end_values_copyratios[-1:][0], arguments['cut_threshold'])
            if arguments['without_phasing']:
                fig.update_yaxes(range = [0, arguments['cut_threshold']])
            else:
                fig.update_yaxes(range=[-(arguments['cut_threshold']), arguments['cut_threshold']])
            #fig.update_yaxes(title='copy ratio (log2)')

            centers_rev = [-x for x in centers[1:]]
            centers_rev.reverse()
            tick_vals = centers_rev + centers

            integer_fractional_means_rev = [x for x in integer_fractional_means[1:]]
            integer_fractional_means_rev.reverse()
            tickt_ext = integer_fractional_means_rev + integer_fractional_means

            fig.update_layout(
                yaxis=dict(
                    tickmode='array',
                    tickvals=[i for i in range(-1000, 1000, 50)],
                    ticktext= [str(abs(i)) for i in range(-1000, 1000, 50)]
                ),
                yaxis2=dict(
                    tickmode='array',
                    tickvals= tick_vals,#[-133, -99, -66, -33, 0.0, 33, 66, 99, 133],#[-99, -33, 0, 33, 99],#[i for i in range(-1000, 1000, 50)],
                    ticktext= tickt_ext#['3-copy', '2-copy', '1-copy', 'neutral','loss', 'neutral', '1-copy', '2-copy', '3-copy'] #[str(abs(i)) for i in range(-1000, 1000, 50)]
                )
            )

        if arguments['pdf_enable']:
            print_chromosome_pdf(fig, chrom, arguments['out_dir_plots'])

        print_chromosome_html(fig, chrom + '_cn', html_graphs, arguments['out_dir_plots']+'/variation_plots/')
        html_graphs.write("  <object data=\"" + chrom + '_cn' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

def plot_copynumbers_scatter_lines(arguments, fig,haplotype_1_values_copyratios, haplotype_1_start_values_copyratios, haplotype_1_end_values_copyratios, haplotype_2_values_copyratios, haplotype_2_start_values_copyratios, haplotype_2_end_values_copyratios, OFFSET):
    haplotype_1_values_copyratios = list(map(float, haplotype_1_values_copyratios))
    haplotype_1_values_copyratios = list(np.asarray(haplotype_1_values_copyratios) + OFFSET)
    haplotype_1_gaps_values = np.full(len(haplotype_1_values_copyratios), 'None')
    haplotype_1_copyratios_values = list(itertools.chain.from_iterable(zip(haplotype_1_values_copyratios, haplotype_1_values_copyratios, haplotype_1_gaps_values)))
    haplotype_1_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_1_start_values_copyratios, haplotype_1_end_values_copyratios, haplotype_1_gaps_values)))

    haplotype_2_values_copyratios = list(map(float, haplotype_2_values_copyratios))
    haplotype_2_values_copyratios = list(np.asarray(haplotype_2_values_copyratios) - OFFSET)
    haplotype_2_gaps_values = np.full(len(haplotype_2_values_copyratios), 'None')
    haplotype_2_copyratios_values = list(itertools.chain.from_iterable(zip(haplotype_2_values_copyratios, haplotype_2_values_copyratios, haplotype_2_gaps_values)))
    haplotype_2_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_2_start_values_copyratios, haplotype_2_end_values_copyratios, haplotype_2_gaps_values)))

    name = "Copynumbers"
    if arguments['without_phasing']:
        add_scatter_trace_copyratios(arguments, fig, ['darkolivegreen'], name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values)
    else:
        add_scatter_trace_copyratios(arguments, fig, ['firebrick','steelblue'], name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values)

def coverage_plots_chromosomes(df, df_phasesets, arguments, thread_pool):
    if not os.path.isdir(arguments['out_dir_plots']+'/coverage_plots'):
        os.mkdir(arguments['out_dir_plots']+'/coverage_plots')
    if not os.path.isdir(arguments['out_dir_plots']+'/bed_output'):
        os.mkdir(arguments['out_dir_plots']+'/bed_output')

    filename = f"{os.path.join(arguments['out_dir_plots'], 'coverage_plots/COVERAGE_INDEX.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")
    haplotype_1_values_updated = []
    haplotype_2_values_updated = []
    hunphased_updated = []
    values_extended = []
    haplotype_1_snps_freqs_updated = []
    haplotype_2_snps_freqs_updated = []

    snps_cpd_points = []
    snps_cpd_points_weights = []

    chr_all = []
    ref_start_values_all = []
    ref_end_values_all = []
    snps_het_all = []
    snps_homo_all = []
    snps_cpd_means_all = []
    df_means_chr_all = []
    df_means_chr_all_ = []
    df_means_chr_all_hp1 = []
    df_means_chr_all_hp2 = []

    if arguments['normal_phased_vcf']:
        get_snp_segments(arguments, arguments['target_bam'][0], thread_pool)
        df_snps = csv_df_chromosomes_sorter('data/snps_frequencies.csv', ['chr', 'pos', 'freq_value_a', 'hp_a', 'freq_value_b', 'hp_b'])
        df_snps = df_snps.drop(df_snps[(df_snps.chr == "chrX") | (df_snps.chr == "chrY")].index)
    if arguments['tumor_vcf']:
        output_phasesets_file_path = vcf_parse_to_csv_for_snps(arguments['tumor_vcf'])
        df_snps_in_csv = csv_df_chromosomes_sorter(output_phasesets_file_path, ['chr', 'pos', 'qual', 'gt', 'dp', 'vaf'])

    chroms = get_contigs_list(arguments['contigs'])

    haplotype_1_segs_dfs = [] #pd.DataFrame()
    haplotype_2_segs_dfs = [] #pd.DataFrame()
    for index, chrom in enumerate(chroms):
        if chrom in chroms:# and (chrom == '1' or chrom == '18'):# and (chrom == 'chr5' or chrom == 'chr16'):
            logging.info('Plots generation for ' + chrom)
            fig = go.Figure()

            df_chrom = df[df['chr'] == chrom]

            if arguments['without_phasing']:
                values = df_chrom.coverage.values.tolist()
                ref_start_values = df_chrom.start.values.tolist()
                ref_end_values = df_chrom.end.values.tolist()
            else:
                df_chrom_phasesets = df_phasesets[df_phasesets['chr'] == chrom]

                #unphased_reads_values = df_chrom.hp3.clip(upper=arguments['cut_threshold']).values.tolist()
                #haplotype_1_values = df_chrom.hp1.clip(upper=arguments['cut_threshold']).values.tolist()
                #haplotype_2_values = df_chrom.hp2.clip(upper=arguments['cut_threshold']).values.tolist()
                unphased_reads_values = df_chrom.hp3.values.tolist()
                haplotype_1_values = df_chrom.hp1.values.tolist()
                haplotype_2_values = df_chrom.hp2.values.tolist()
                ref_start_values = df_chrom.start.values.tolist()
                ref_end_values = df_chrom.end.values.tolist()

                haplotype_1_values_phasesets = df_chrom_phasesets.hp1.clip(upper=arguments['cut_threshold']).values.tolist()
                haplotype_2_values_phasesets = df_chrom_phasesets.hp2.clip(upper=arguments['cut_threshold']).values.tolist()
                ref_start_values_phasesets = df_chrom_phasesets.start.values.tolist()
                ref_end_values_phasesets = df_chrom_phasesets.end.values.tolist()
            ################################################################################
            if arguments['normal_phased_vcf']:
                haplotype_1_values, haplotype_2_values = snps_mean(df_snps, ref_start_values, chrom)
            ################################################################################
            if arguments['tumor_vcf']:
                logging.info('hetrozygous phased snps frequencies coverage module')
                ref_start_values_updated, snps_het_counts, snps_homo_counts, centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends = get_snps_frquncies_coverage(df_snps_in_csv, chrom, ref_start_values, arguments['bin_size'])
                if arguments['without_phasing']:
                    if centromere_region_starts:
                        _,_,_, centromere_region_starts, centromere_region_ends  = detect_alter_loh_regions(arguments, 'centromere/no-coverage', chrom, ref_end_values, values, values, values, centromere_region_starts, centromere_region_ends, True)
                    if loh_region_starts:
                        _,_,_, loh_region_starts, loh_region_ends = detect_alter_loh_regions(arguments, 'loss-of-heterozygosity', chrom, ref_end_values, values, values, values, loh_region_starts, loh_region_ends, True)
                else:
                    if centromere_region_starts:
                        haplotype_1_values, haplotype_2_values, unphased_reads_values, centromere_region_starts, centromere_region_ends  = detect_alter_loh_regions(arguments, 'centromere/no-coverage', chrom, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values, centromere_region_starts, centromere_region_ends, True)
                        haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = loh_regions_phasesets(centromere_region_starts, centromere_region_ends, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

                    if loh_region_starts:
                        haplotype_1_values, haplotype_2_values, unphased_reads_values, loh_region_starts, loh_region_ends = detect_alter_loh_regions(arguments, 'loss-of-heterozygosity', chrom, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values, loh_region_starts, loh_region_ends, True)
                        haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = loh_regions_phasesets(loh_region_starts, loh_region_ends, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
            else:
                loh_region_starts = []
                loh_region_ends = []
            ################################################################################
            if arguments['phaseblock_flipping_enable']:
                logging.info('phaseblock flipping module')
                is_simple_heuristics = True
                if ref_start_values_phasesets:
                    is_simple_heuristics = False
                haplotype_1_values, haplotype_2_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = \
                phaseblock_flipping(chrom, arguments, is_simple_heuristics, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                        haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

                #detect centromeres
                ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = detect_centromeres(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values, ref_start_values, arguments['bin_size'])
                # infer missing phaseblocks
                # ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = infer_missing_phaseblocks(
                #     ref_start_values, ref_end_values, ref_start_values_phasesets, ref_end_values_phasesets,
                #     haplotype_1_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values, arguments['bin_size'])
                # create more contiguous phaseblocks
                haplotype_1_values_phasesets_conti, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, haplotype_2_values_phasesets_conti, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2 = contiguous_phaseblocks(haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

                # flip phaseblocks based on more contiguous phaseblocks
                haplotype_1_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values = flip_phaseblocks_contigous(chrom, arguments, haplotype_1_values_phasesets_conti, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, haplotype_2_values_phasesets_conti, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2, ref_start_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values, haplotype_2_values)
            ################################################################################
            if arguments['smoothing_enable']:
                logging.info('smoothing module')
                haplotype_1_values, haplotype_2_values, unphased_reads_values = smoothing(haplotype_1_values, haplotype_2_values, unphased_reads_values, conv_window_size=5)
            ################################################################################
            #Plots
            if arguments['without_phasing']:
                add_scatter_trace_coverage(fig, ref_start_values, values, name='Unphased', text=None, yaxis=None, opacity=0.7, color='olive')
            else:
                add_scatter_trace_coverage(fig, ref_start_values, haplotype_1_values, name='HP-1', text=None, yaxis=None, opacity=0.7, color='firebrick')
                add_scatter_trace_coverage(fig, ref_start_values, haplotype_2_values, name='HP-2', text=None, yaxis=None, opacity=0.7, color='steelblue')

                if arguments['unphased_reads_coverage_enable']:
                    add_scatter_trace_coverage(fig, ref_start_values, unphased_reads_values, name='Unphased', text=None, yaxis=None, opacity=0.7, color='olive')

            plots_add_markers_lines(fig)
            ################################################################################
            if arguments['phaseblocks_enable']:
                logging.info('phaseblocks plots module')
                gaps_values = np.full(len(haplotype_1_values_phasesets), 'None')
                haplotype_1_phaseblocks_values = list(itertools.chain.from_iterable(zip(haplotype_1_values_phasesets, haplotype_1_values_phasesets, gaps_values)))
                haplotype_2_phaseblocks_values = list(itertools.chain.from_iterable(zip(haplotype_2_values_phasesets, haplotype_2_values_phasesets, gaps_values)))
                phaseblocks_positions = list(itertools.chain.from_iterable(zip(ref_start_values_phasesets, ref_end_values_phasesets, gaps_values)))

                add_scatter_trace_phaseblocks(fig, phaseblocks_positions, haplotype_1_phaseblocks_values, haplotype_2_phaseblocks_values)
            ################################################################################
            if arguments['breakpoints_enable']:
                breakpoints = get_breakpoints(chrom, arguments['breakpoints_file'])
                add_scatter_trace_breakpoints(fig, breakpoints)

            plots_layout_settings(fig, chrom, arguments, ref_end_values[-1:][0], arguments['cut_threshold'])
            ################################################################################
            if arguments['pdf_enable']:
                print_chromosome_pdf(fig, chrom, arguments['out_dir_plots'])
            ################################################################################
            chr = range(len(ref_start_values))
            chr_all.extend([chrom for ch in chr])
            ref_start_values_all.extend(ref_start_values)
            ref_end_values_all.extend(ref_end_values)

            if arguments['without_phasing']:
                df_snps_freqs_chr = whole_genome_combined_df(arguments, chrom, chr, ref_start_values, ref_end_values, values, values, values)
                df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, states, centers, stdev = apply_copynumbers(df_snps_freqs_chr, values, values, arguments, snps_cpd_means, [])
                snps_cpd_points.extend(snps_cpd_means)
                if arguments['enable_debug']:
                    copy_number_plots_per_chromosome(centers, centers, ref_start_values, values, df_segs_hp1, values, df_segs_hp2, arguments, chrom, html_graphs, loh_region_starts, loh_region_ends)
                haplotype_1_segs_dfs.append(df_segs_hp1)
                haplotype_2_segs_dfs.append(df_segs_hp2)
                if arguments['enable_debug']:
                    change_point_detection(values, ref_start_values, ref_end_values, arguments, chrom, html_graphs, 1, color='#6A5ACD')

            else:
                df_snps_freqs_chr = whole_genome_combined_df(arguments, chrom, chr, ref_start_values, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values)

                #change point detection
                snps_cpd_means, df_means_chr = change_point_detection_means(arguments, df_snps_freqs_chr)
                df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, states, centers, stdev = apply_copynumbers(df_snps_freqs_chr, haplotype_1_values, haplotype_2_values, arguments, snps_cpd_means, [])
                snps_cpd_points.extend(snps_cpd_means)

                #TODO For debug - histo_clusters and cpds
                #add_histo_clusters_plot(df_cnr_hp1.log2.values.tolist(), df_cnr_hp2.log2.values.tolist(), states, centers, stdev, arguments, chrom, html_graphs)
                if arguments['enable_debug']:
                    change_point_detection(haplotype_1_values, ref_start_values, ref_end_values, arguments, chrom, html_graphs, 1, color='#6A5ACD')
                    change_point_detection(haplotype_2_values, ref_start_values, ref_end_values, arguments, chrom, html_graphs, 2, color='#2E8B57')
                if arguments['enable_debug']:
                    copy_number_plots_per_chromosome(centers, ref_start_values, haplotype_1_values, df_segs_hp1, haplotype_2_values, df_segs_hp2, arguments, chrom, html_graphs, loh_region_starts, loh_region_ends)

                haplotype_1_segs_dfs.append(df_segs_hp1)
                haplotype_2_segs_dfs.append(df_segs_hp2)

            if arguments['without_phasing']:
                values_extended.extend(values)
            else:
                haplotype_1_values_updated.extend(haplotype_1_values)
                haplotype_2_values_updated.extend(haplotype_2_values)
                hunphased_updated.extend(unphased_reads_values)

            snps_cpd_means_all.extend(snps_cpd_means)
            if arguments['without_phasing']:
                df_means_chr_all.append(df_means_chr)
            else:
                df_means_chr_all_hp1.append(df_means_chr[0])
                df_means_chr_all_hp2.append(df_means_chr[1])

            print_chromosome_html(fig, chrom + '_cov', html_graphs, arguments['out_dir_plots']+'/coverage_plots/')
            html_graphs.write("  <object data=\"" + chrom + '_cov' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

    df_snps_freqs = pd.DataFrame(list(zip(chr_all, ref_start_values_all, ref_end_values_all, haplotype_1_values_updated, haplotype_2_values_updated, hunphased_updated)), columns=['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

    if not arguments['without_phasing']:
        df_means_chr_all_.append(pd.concat(df_means_chr_all_hp1))
        df_means_chr_all_.append(pd.concat(df_means_chr_all_hp2))

    html_graphs.write("</body></html>")
    if arguments['without_phasing']:
        return values_extended, values_extended, values_extended, df_snps_freqs, snps_cpd_points, pd.concat(df_means_chr_all)
    else:
        return haplotype_1_values_updated, haplotype_2_values_updated, hunphased_updated, df_snps_freqs, snps_cpd_points, df_means_chr_all_

def whole_genome_combined_df(arguments, chrom, chr, ref_start_values, ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean, snps_haplotype_mean):
    if arguments['without_phasing']:
        combined_df = pd.DataFrame(list(
            zip([chrom for ch in chr], ref_start_values, ref_end_values, snps_haplotype1_mean)), columns=['chr', 'start', 'end', 'coverage'])
    else:
        combined_df = pd.DataFrame(list(
            zip([chrom for ch in chr], ref_start_values, ref_end_values, snps_haplotype1_mean,
                snps_haplotype2_mean, snps_haplotype2_mean)), columns=['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
    return combined_df
def print_genome_pdf(fig, genome, coverage_plots_path):
    filename = f"{os.path.join(coverage_plots_path, genome + '.pdf')}"
    plotly.io.write_image(fig,  filename, format='pdf')
def print_chromosome_pdf(fig, chrom, coverage_plots_path):
    filename = f"{os.path.join(coverage_plots_path, chrom + '.pdf')}"
    plotly.io.write_image(fig,  filename, format='pdf')

def print_chromosome_html(fig, chrom, html_graphs, coverage_plots_path):
    fname = f"{os.path.join(coverage_plots_path, chrom + '.html')}"
    plotly.offline.plot(fig, filename=fname,auto_open=False)


def plots_genome_coverage(df_hp1, df_hp2, df_unphased, arguments):
    # TODO genome-wide plots

    lengths = []
    chroms = get_contigs_list(arguments['contigs'])
    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_hp1[df_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_hp1)*arguments['bin_size'], arguments['bin_size'], dtype=int)

    fig = go.Figure()
    add_scatter_trace_coverage(fig, indices, df_hp1.hp1.values.tolist(), name='HP-1', text=None,
                               yaxis=None, opacity=0.7, color='firebrick')
    add_scatter_trace_coverage(fig, indices, df_hp2.hp2.values.tolist(), name='HP-2', text=None,
                               yaxis=None, opacity=0.7, color='steelblue')

    if arguments['unphased_reads_coverage_enable']:
        add_scatter_trace_coverage(fig, indices, df_unphased.hp3.values.tolist(), name='Unphased reads', text=None,
                                   yaxis=None, opacity=0.5, color='olive')
    plots_add_markers_lines(fig)

    current = 0
    label_pos = []
    chroms = get_contigs_list(arguments['contigs'])
    for index, chrom in enumerate(chroms):
        current += lengths[index]
        label_pos.append(round(current - (lengths[index] / 2))*arguments['bin_size'])
        fig.add_vline(x=current*arguments['bin_size'], y0=-10, y1=150, line_width=1, line_dash="dashdot", line_color="green")

    fig.update_layout(
        xaxis=dict(
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        ),
        font=dict(size=18, color="black"))

    plots_layout_settings(fig, 'Genome', arguments, indices[-1:][0], arguments['cut_threshold'])
    #fig.update_yaxes(range=[-10, 10])
    #fig.update_yaxes(title='copy ratio (log2)')
    fig.update_layout(width=1280, height=600,)
    fig.update_layout(legend=dict(orientation='h', xanchor="center", x=0.5, y=1.1))

    if arguments['pdf_enable']:
        print_genome_pdf(fig, arguments['genome_name'], arguments['out_dir_plots'])

    fig.write_html(arguments['out_dir_plots'] +'/'+ arguments['genome_name'] + "_genome_coverage.html")

def copy_number_plots_genome(centers, integer_fractional_centers, df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, df_unphased, arguments):
    # TODO genome-wide plots

    lengths = []
    chroms = get_contigs_list(arguments['contigs'])
    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chromosome'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_cnr_hp1)*arguments['bin_size'], arguments['bin_size'], dtype=int)

    fig = go.Figure()
    custom_text_data_hp1 = df_cnr_hp1.start.values.tolist()
    custom_text_data_hp2 = df_cnr_hp2.start.values.tolist()
    if arguments['without_phasing']:
        add_scatter_trace_coverage(fig, indices, df_cnr_hp1.log2.values.tolist(), name='Unphased', text=None,
                                   yaxis=None, opacity=0.7, color='olive', visibility='legendonly')
    else:
        add_scatter_trace_coverage(fig, indices, df_cnr_hp1.log2.values.tolist(), name='HP-1', text=custom_text_data_hp1,
                                   yaxis=None, opacity=0.7, color='firebrick')
        add_scatter_trace_coverage(fig, indices, df_cnr_hp2.log2.values.tolist(), name='HP-2', text=custom_text_data_hp2,
                                   yaxis=None, opacity=0.7, color='steelblue')
        if arguments['unphased_reads_coverage_enable']:
            add_scatter_trace_coverage(fig, indices, df_unphased.hp3.values.tolist(), name='Unphased', text=None,
                                   yaxis=None, opacity=0.7, color='olive', visibility='legendonly')

    plots_add_markers_lines(fig)

    current = 0
    label_pos = []
    chroms = get_contigs_list(arguments['contigs'])
    for index, chrom in enumerate(chroms):
        current += lengths[index]
        label_pos.append(round(current - (lengths[index] / 2))*arguments['bin_size']) #y0=-10, y1=arguments['cut_threshold']
        fig.add_vline(x=current*arguments['bin_size'],  line_width=1.5, line_dash="dashdot", line_color="green")

    fig.update_layout(
        xaxis=dict(
            tickangle=90,
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        ),
        font=dict(size=18, color="black"))

    if arguments['copynumbers_enable']:
        offset = 0
        haplotype_1_start_values = []
        haplotype_1_end_values = []

        haplotype_2_start_values = []
        haplotype_2_end_values = []
        chroms = get_contigs_list(arguments['contigs'])
        for index, chrom in enumerate(chroms):
            if not chrom == 'chrX' and not chrom == 'chrY':
                df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
                df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
                haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
                haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()
                haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
                haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()
                if chrom == arguments['contigs'].split('-')[0]:
                    haplotype_1_start_values.extend(haplotype_1_start_values_copyratios)
                    haplotype_1_end_values.extend(haplotype_1_end_values_copyratios)
                    haplotype_2_start_values.extend(haplotype_2_start_values_copyratios)
                    haplotype_2_end_values.extend(haplotype_2_end_values_copyratios)
                else:
                    offset += lengths[index-1] * arguments['bin_size']
                    haplotype_1_start_values.extend([x + offset for x in haplotype_1_start_values_copyratios])
                    haplotype_1_end_values.extend([x + offset for x in haplotype_1_end_values_copyratios])
                    haplotype_2_start_values.extend([x + offset for x in haplotype_2_start_values_copyratios])
                    haplotype_2_end_values.extend([x + offset for x in haplotype_2_end_values_copyratios])

        haplotype_1_gaps_values = np.full(len(df_segs_hp1.state.values.tolist()), 'None')
        haplotype_1_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp1.state.values.tolist(), df_segs_hp1.state.values.tolist(), haplotype_1_gaps_values)))
        haplotype_1_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_1_start_values, haplotype_1_end_values, haplotype_1_gaps_values)))

        haplotype_2_gaps_values = np.full(len(df_segs_hp2.state.values.tolist()), 'None')
        haplotype_2_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp2.state.values.tolist(), df_segs_hp2.state.values.tolist(), haplotype_2_gaps_values)))
        haplotype_2_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_2_start_values, haplotype_2_end_values, haplotype_2_gaps_values)))

        if arguments['copynumbers_subclonal_enable']:
            search_list = [30.0, 78.0, 130]
            haplotype_1_copyratios_values_normal = [-33.0 if element in search_list else element for element in haplotype_1_copyratios_values]
            search_list = [30.0, 78.0, 130, 'None']
            haplotype_1_copyratios_values_sub = [element if element in search_list else -33.0 for element in haplotype_1_copyratios_values]

            search_list = [30.0, 78.0, 130]
            haplotype_2_copyratios_values_normal = [-33.0 if element in search_list else element for element in
                                                    haplotype_2_copyratios_values]
            search_list = [30.0, 78.0, 130, 'None']
            haplotype_2_copyratios_values_sub = [element if element in search_list else -33.0 for element in
                                                 haplotype_2_copyratios_values]
            OFFSET = arguments['cut_threshold']/150
            haplotype_1_copyratios_values_normal = [x if x == 'None' else x + OFFSET for x in
                                                    haplotype_1_copyratios_values_normal]
            haplotype_2_copyratios_values_normal = [x if x == 'None' else x - OFFSET for x in
                                                    haplotype_2_copyratios_values_normal]
            haplotype_1_copyratios_values_sub = [x if x == 'None' else x + OFFSET for x in
                                                    haplotype_1_copyratios_values_sub]
            haplotype_2_copyratios_values_sub = [x if x == 'None' else x - OFFSET for x in
                                                    haplotype_2_copyratios_values_sub]

            name = "Copynumbers"
            add_scatter_trace_copyratios(arguments, fig, ['firebrick','steelblue'], name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values_normal, haplotype_2_copyratios_values_normal)
            name = "Copynumbers (subclonal)"
            add_scatter_trace_copyratios(arguments, fig, ['#E95F0A','#6EC5E9'], name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values_sub, haplotype_2_copyratios_values_sub)
        else:
            if arguments['without_phasing']:
                OFFSET = 0
                colors = ['darkolivegreen']
            else:
                OFFSET = arguments['cut_threshold']/150
                colors = ['firebrick', 'steelblue']
            haplotype_1_copyratios_values = [x if x == 'None' else x + OFFSET for x in
                                                    haplotype_1_copyratios_values]
            haplotype_2_copyratios_values = [x if x == 'None' else x - OFFSET for x in
                                                    haplotype_2_copyratios_values]
            name = "Copynumbers"

            add_scatter_trace_copyratios(arguments, fig, colors, name, haplotype_1_copyratios_positions,
                                         haplotype_2_copyratios_positions, haplotype_1_copyratios_values,
                                         haplotype_2_copyratios_values)

    plots_layout_settings(fig, 'Genome', arguments, indices[-1:][0], arguments['cut_threshold'])

    fig.update_layout(
        yaxis=dict(
            tickmode='array',
            tickvals=[i for i in range(-1000, 1000, 25)],
            ticktext=[str(abs(i)) for i in range(-1000, 1000, 25)]
        ),
        yaxis2=dict(
            tickmode='array',
            tickvals= centers, #[i for i in range(-1000, 1000, 25)],#sorted(integer_centers + fractional_centers),#[0, 6,12,18], #13,26,39,5
            # [-99, -33, 0, 33, 99],#[i for i in range(-1000, 1000, 50)],
            ticktext= integer_fractional_centers #['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
            # [str(abs(i)) for i in range(-1000, 1000, 50)]
        )
    )
    fig.update_yaxes(range=[-1, arguments['cut_threshold']])
    #fig.update_yaxes(title='coverge (mean depth)')
    #fig.update_yaxes(title_text="<b>secondary</b> yaxis title", secondary_y=True)
    fig.update_layout(width=1480, height=600,)
    fig.update_layout(legend=dict(orientation = 'h', xanchor = "center", x = 0.5, y= 1.1))

    if arguments['pdf_enable']:
        print_genome_pdf(fig, arguments['genome_name']+'_copynumbers', arguments['out_dir_plots'])

    fig.write_html(arguments['out_dir_plots'] +'/'+ arguments['genome_name'] + "_genome_copynumber.html")

    return centers, integer_fractional_centers

def plots_add_markers_lines(fig):
    # style all the traces
    fig.update_traces(
        # hoverinfo="name+x+text+y",
        # line={"width": 0.5},
        marker={"size": 2},
        mode="markers",
        showlegend=True
    )

def add_scatter_trace_coverage(fig, x, y, name, text, yaxis, opacity, color, visibility=True):

    fig.add_trace(go.Scatter(
        #legendgroup="group1",  # this can be any string, not just "group"
        #legendgrouptitle_text="Coverage",
        x=x,
        y=y,
        name=name,
        text=text,
        #yaxis="y5",
        opacity=opacity,
        marker_color=color,
        visible=visibility,
    ))

def add_scatter_trace_phaseblocks(fig, phaseblocks_positions, haplotype_1_phaseblocks_values, haplotype_2_phaseblocks_values):
    fig.add_trace(go.Scatter(
        #legendgroup="group3",  # this can be any string, not just "group"
        #legendgrouptitle_text="Phaseblocks",
        x=phaseblocks_positions,
        y=haplotype_1_phaseblocks_values,
        name="HP-1",
        text=phaseblocks_positions,
        #yaxis="y5",
        line = dict(shape = 'spline', color = 'gray', width= 1, dash = 'solid'),
        mode='lines+markers',
        marker={"size": 5},
        opacity=0.5,
        marker_color=['dimgray', 'darkgray', 'white']*len(phaseblocks_positions),
        showlegend=True,
        marker_symbol='diamond-wide',
        hoverinfo = "x+name+y+text",
        #legendgroup="group2",
        #legendgrouptitle_text="Phaseblocks",
    ))

    fig.add_trace(go.Scatter(
        #legendgroup="group3",
        x=phaseblocks_positions,
        y=haplotype_2_phaseblocks_values,
        name="HP-2",
        text=phaseblocks_positions,
        #yaxis="y5",
        line = dict(shape = 'spline', color = 'green', width= 1, dash = 'solid'),
        mode='lines+markers',
        marker={"size": 5},
        opacity=0.5,
        marker_color=['darkgreen', 'limegreen', 'white']*len(phaseblocks_positions),
        showlegend=True,
        marker_symbol='diamond-wide',
        hoverinfo = "x+name+y+text",
        #legendgroup="group2",
    ))

def add_scatter_trace_copyratios(arguments, fig, colors, name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values):
    if arguments['without_phasing']:
        name_leg = 'Unphased'
    else:
        name_leg = 'HP-1'
    fig.add_trace(go.Scatter(
        #legendgroup="group2",
        #legendgrouptitle_text=name,
        x=haplotype_1_copyratios_positions,
        y= haplotype_1_copyratios_values,
        name=name_leg,
        text=haplotype_1_copyratios_positions,
        yaxis="y2",
        line = dict(shape = 'spline', color = colors[0], width= 5, dash = 'solid'),
        mode='lines',
        #marker={"size": 5},
        opacity=0.9,
        marker_color=[colors[0], colors[0], 'white']*len(haplotype_1_copyratios_positions),
        showlegend=True,
        #marker_symbol='diamond-wide',
        #hoverinfo = "x+name+y+text",
        #legendgroup="group2",
        #legendgrouptitle_text="Phaseblocks",
    ))
    if arguments['without_phasing'] == False:
        fig.add_trace(go.Scatter(
            #legendgroup="group2",
            x=haplotype_2_copyratios_positions,
            y=haplotype_2_copyratios_values,
            name='HP-2',
            text=haplotype_2_copyratios_positions,
            yaxis="y2",
            line = dict(shape = 'spline', color = colors[1], width= 5, dash = 'solid'),
            mode='lines',
            #marker={"size": 5},
            opacity=0.9,
            marker_color=[colors[1], colors[1], 'white']*len(haplotype_2_copyratios_positions),
            showlegend=True,
            #marker_symbol='diamond-wide',
            #hoverinfo = "x+name+y+text",
            #legendgroup="group2",
        ))

def add_scatter_trace_breakpoints(fig, break_points):
    fig.add_trace(go.Scatter(
        x=break_points,
        y=[-1,99,-1]*round((len(break_points)/3)),
        name="BP",
        text=break_points,
        yaxis="y5",
        opacity=0.5,
        mode='lines',
        line=dict(shape='linear', color='dodgerblue', width=1, dash='solid'),
        showlegend=False,
        hoverinfo="x+name+text",
    ))

def plots_layout_settings(fig, chrom, arguments, limit_x, limit_y):
    # Update axes
    fig.update_layout(
        xaxis=dict(
            #autorange=True,
            type="linear",
            showline=True,
            zeroline=True,
            linecolor = "dimgray",
            range=[0,limit_x*1.01]
        ),
        yaxis=dict(
            linecolor="dimgray",
            range=[0, limit_y+1],
            side="left",
            tickfont={"color": "dimgray"},
            tickmode="auto",
            ticks="outside",
            title="<b>Coverage depth</b> (per bin)",
            titlefont={"color": "dimgray"},
            type="linear",
            showline=True,
            zeroline=True,
        ),
        yaxis2=dict(
            linecolor="dimgray",
            range=[1, arguments['cut_threshold'] + 5],
            side="right",
            tickfont={"color": "dimgray"},
            tickmode="auto",
            ticks="outside",
            title="<b>Copies</b> (integers/fractions)",
            titlefont={"color": "dimgray"},
            type="linear",
            # showline=True,
            # zeroline=True,
            anchor="x",
            overlaying="y",
        )
    )
    # Update layout
    fig.update_layout(
        template="plotly_white",
        font_family="Times New Roman"
    )

    fig.update_layout(
        title=chrom,
    )
    #Legend
    fig.update_layout(legend=dict(
        orientation = 'h', xanchor = "center", x = 0.45, y= 1.2, #orientation = 'v', xanchor = "center", x = 1.08, y= .5
    ))
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    fig.update_xaxes(tick0=0.0, rangemode="nonnegative")

    fig.update_layout(legend={'itemsizing': 'constant'})

    fig.update_layout(font_family= "Times New Roman")

    fig.update_layout(
        title={
            'text': chrom + ' - ' +arguments['genome_name'],
            'y':0.96,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family = "Courier New",
        font_color = "dimgray",
        title_font_family = "Times New Roman",
        title_font_color = "red",
        legend_title_font_color = "green",
    )
    #Size
    fig.update_layout(
        width=680,
        height=400,
       )

def plot_bins_ratios(hp1,hp2, arguments):
    import plotly.express as px
    fig = go.Figure()

    numpy.clip(hp1, a_min=1, a_max=300)
    numpy.clip(hp2, a_min=1, a_max=300)

    add_scatter_trace_coverage(fig, hp1/hp2, hp2/hp1, name='HP-1/HP-2 ratio', text=None, yaxis=None,
                               opacity=0.7, color='firebrick')

    fig.write_html(arguments['out_dir_plots'] +'/'+ arguments['genome_name'] + "_genome_ratio.html")

def change_point_detection(data, start, ends, arguments, chrom, html_graphs, hp, color):
    import ruptures as rpt
    fig = go.Figure()
    #starts = [i for i in range(0, len(data), 50000)]
    add_scatter_trace_coverage(fig, start, data, name='HP-'+str(hp), text=None, yaxis=None,
                               opacity=0.7, color=color)

    data = np.array(data, dtype='int') #numpy.clip(data, a_min=0, a_max=1000)
    algo = rpt.Pelt(model="rbf", jump=25).fit(data)
    result = algo.predict(pen=10)
    change_points = [i for i in result if i < len(data)]
    for i, point in enumerate(change_points):
        fig.add_vline(x=point*50000, y0=-10, y1=500, line_width=1, line_dash="dash",
                  line_color=color)

    plots_add_markers_lines(fig)
    plots_layout_settings(fig, chrom, arguments, ends[-1:][0], arguments['cut_threshold'])

    print_chromosome_html(fig, chrom + '_hp_'  + str(hp), html_graphs, arguments['out_dir_plots'])
    html_graphs.write("  <object data=\"" + chrom + '_hp_'  + str(hp)  + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")



def plot_bins_histograms(hp1,hp2, arguments, chrom):
    import plotly.express as px
    from sklearn.cluster import KMeans
    import plotly.figure_factory as ff

    #hp1=numpy.clip(hp1, a_min=1, a_max=300)
    #hp2=numpy.clip(hp2, a_min=1, a_max=300)

    hp1 = [i for i in hp1 if i > 0 and i < 600]
    hp2 = [i for i in hp2 if i > 0 and i < 600]

    hp1 = np.array(hp1, dtype='int')
    hp2 = np.array(hp2, dtype='int')

    ##########################################
    import ruptures as rpt
    algo = rpt.Pelt(model="rbf").fit(hp1)
    result = algo.predict(pen=10)
    change_points = [i for i in result if i < len(hp1)]

    # Plot the time series
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 6))
    plt.scatter([i for i in range(0,len(hp1))], hp1)

    # Plot the abrupt changes
    for i in range(len(change_points)):
        ax.axvline(change_points[i], color='r')

    # Add labels, grid, and legend
    ax.set_xlabel('Time')
    ax.set_ylabel('Signal')
    ax.set_title(chrom)
    ax.legend()
    ax.grid(True)
    ax.set_ylim([0, 160])

    # Show the plot
    plt.show()
    ######################################
    import ruptures as rpt
    algo = rpt.Pelt(model="rbf").fit(hp2)
    result = algo.predict(pen=10)
    change_points = [i for i in result if i < len(hp2)]

    # Plot the time series
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 6))
    plt.scatter([i for i in range(0,len(hp2))], hp2)

    # Plot the abrupt changes
    for i in range(len(change_points)):
        ax.axvline(change_points[i], color='g')

    # Add labels, grid, and legend
    ax.set_xlabel('Time')
    ax.set_ylabel('Signal')
    ax.set_title(chrom)
    ax.legend()
    ax.grid(True)
    ax.set_ylim([0, 160])

    # Show the plot
    plt.show()
    ######################################

    from sklearn.neighbors import KernelDensity
    from numpy import exp
    model = KernelDensity(bandwidth=2, kernel='gaussian')
    values = np.concatenate((hp1, hp2))
    values = values.reshape(-1, 1)
    model.fit(values)
    probabilities_hp1 = model.score_samples(hp1.reshape(-1, 1))
    probabilities_hp1 = exp(probabilities_hp1)
    probabilities_hp2 = model.score_samples(hp2.reshape(-1, 1))
    probabilities_hp2 = exp(probabilities_hp2)

    df = pd.DataFrame(dict(
        Haplotypes=np.concatenate((["HP-1"] * len(hp1), ["HP-2"] * len(hp2))),
        coverage=np.concatenate((hp1, hp2)),
        probabilities = np.concatenate((probabilities_hp1, probabilities_hp2))
    ))

    # df = pd.DataFrame(dict(
    #     Haplotypes=(["HP-1/HP-2"] * len(hp1)),
    #     coverage=(hp1/(hp1+hp2)),
    # ))

    #X=list(zip(hp1,start))

    # df = pd.DataFrame(dict(
    #     Haplotypes=(["HP-1 & HP-2"] * (len(hp1)*2)),
    #     coverage=np.concatenate((hp1, hp2)),
    # ))

    hp1[hp1 != 0]
    hp = list(hp1)
    hp1=hp1.reshape(-1, 1)

    hp2[hp2 != 0]
    hp = list(hp2)
    hp2=hp2.reshape(-1, 1)

    new = np.concatenate((hp1, hp2))#list(zip(hp1, hp2))

    Sum_of_squared_distances = []
    K = range(1, 15)
    for k in K:
        km = KMeans(n_clusters=k)
        km = km.fit(new)
        Sum_of_squared_distances.append(km.inertia_)

    # fig = px.line(x=K, y=Sum_of_squared_distances, markers=True)
    # fig.update_yaxes(title='Sum_of_squared_distances')
    # fig.update_xaxes(title='Number_of_clusters')
    # plotly.io.write_image(fig, "kmeans.pdf", format='pdf')

    from kneed import KneeLocator
    from sklearn.metrics import silhouette_samples, silhouette_score

    kn = KneeLocator(x=K, y=Sum_of_squared_distances, curve='convex', direction='decreasing')
    print(kn.knee)

    ###########################

    ###########################
    km1 = KMeans(n_clusters=kn.knee, n_init="auto", max_iter=100, random_state=0).fit(new)
    labels1 = list(km1.labels_)
    print(km1.cluster_centers_)
    print(silhouette_score(new, labels1))
    ###########################
    km = KMeans(n_clusters=kn.knee+1, n_init="auto", max_iter = 100, random_state=0).fit(new)
    labels = list(km.labels_)
    u_labels = np.unique(labels)
    #centers = list(np.concatenate(km.cluster_centers_))
    print(km.cluster_centers_)
    print(silhouette_score(new, labels))



    stdev = []
    clusters = []
    for i in range(len(u_labels)):
        clusters.append([int(a) for a, b in zip(new, labels) if b == i])
        stdev.append(numpy.std([(a, b) for a, b in zip(new, labels) if b == i]))

    # Group data togetherhist_data = {list: 4} [[array([76.40792133]), array([76.40792133]), array([76.40792133]), array([76.40792133]), array([76.40792133]), array([76.40792133]), array([76.40792133]), array([76.40792133]), array([76.40792133]), array([76.40792133]), array([76.40792133]), array([76.40... View
    hist_data = clusters # [i for i in clusters]
    group_labels = ['Cluster '+str(i) for i in range(len(clusters))]
    fig = ff.create_distplot(hist_data, group_labels, curve_type='normal', # override default 'kde'
    )

    fig = px.histogram(df, x="coverage", y="probabilities", color="Haplotypes", barmode="overlay", marginal="violin")#, log_y=True)

    #fig.update_yaxes(range=[1, 1000])
    #fig.update_xaxes(range=[1, 36])


    fig.update_layout(
        title={
            'text': 'Genome - ' + arguments['genome_name'],
            'y': 0.96,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family="Courier New",
        font_color="dimgray",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="green",
    )

    fig.write_html(arguments['out_dir_plots'] +'/'+ arguments['genome_name'] + "_genome_histogram_"+chrom+".html")

def plot_bins_distributions(hp1,hp2, arguments):
    import plotly.express as px
    from sklearn.cluster import KMeans
    import plotly.figure_factory as ff

    hp1=numpy.clip(hp1, a_min=1, a_max=300)
    hp2=numpy.clip(hp2, a_min=1, a_max=300)

    import plotly.figure_factory as ff
    import numpy as np

    group_labels = ['Group 1', 'Group 2']

    colors = ['slategray', 'magenta']

    # Create distplot with curve_type set to 'normal'
    fig = ff.create_distplot([hp1, hp2], group_labels, bin_size=.5,
                             curve_type='normal',  # override default 'kde'
                             colors=colors)

    # Add title
    fig.update_layout(title_text='Distplot with Normal Distribution')

    fig.write_html(arguments['out_dir_plots'] +'/'+ arguments['genome_name'] + "_genome_histogram_distplot.html")


def add_histo_clusters_plot(depth_values_hp1, depth_values_hp2, labels, means, covar, arguments, chrom, html_graphs):
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import math
    # fig = go.Figure()

    depth_values_hp1 = np.array(depth_values_hp1, dtype='int')
    depth_values_hp2 = np.array(depth_values_hp2, dtype='int')

    depth_values_hp1 = depth_values_hp1.reshape(-1, 1)
    depth_values_hp2 = depth_values_hp2.reshape(-1, 1)
    Y = np.concatenate([depth_values_hp1, depth_values_hp2])

    fig = make_subplots(rows=1, cols=2, shared_yaxes=True, column_widths=[0.8, 0.2], vertical_spacing=0.01,
                        horizontal_spacing=0.01)

    y = np.concatenate(list(Y)).ravel().tolist()
    lst = [i for i in range(0, (len(y) // 2) * 50000, 50000)]
    x = lst + lst

    cdict = {0: '#1f77b4',  # muted blue
             1: '#ff7f0e',  # safety orange
             2: '#2ca02c',  # cooked asparagus green
             3: '#d62728',  # brick red
             4: '#9467bd',  # muted purple
             5: '#8c564b',  # chestnut brown
             6: '#e377c2',  # raspberry yogurt pink
             7: '#7f7f7f',  # middle gray
             8: '#bcbd22',  # curry yellow-green
             9: '#17becf'  # blue-teal
             }
    ldict = {0: 'Cluster_1', 1: 'Cluster_2', 2: 'Cluster_3', 3: 'Cluster_4', 4: 'Cluster_5', \
             5: 'Cluster_6', 6: 'Cluster_7', 7: 'Cluster_8', 8: 'Cluster_9', 9: 'Cluster_10'}

    for g in np.unique(labels):
        ix = [index for index, i in enumerate(x) if labels[index] == g]
        xn = [x[i] for i in ix]
        yn = [y[i] for i in ix]
        fig.add_trace(go.Scatter(x=xn, y=yn, mode='markers', marker=dict(color=cdict[g]), name=ldict[g], opacity=0.7, ),
                      row=1, col=1)

    fig.update_traces(
        # hoverinfo="name+x+text+y",
        # line={"width": 0.5},
        marker={"size": 2},
        mode="markers",
        showlegend=True
    )

    for i in range(len(means)):
        fig.add_hline(y=means[i], line_width=2,
                      line=dict(dash='solid'), line_color=cdict[i], annotation_position="top right", annotation_text="mean_" + str(i + 1))

        fig.add_hline(y=means[i] + covar[i], line_width=1,
                      line=dict(dash='dash'), line_color=cdict[i], annotation_position="top left", annotation_text="stdev_" + str(i + 1))
        fig.add_hline(y=means[i] - covar[i], line_width=1,
                      line=dict(dash='dash'), line_color=cdict[i], annotation_position="top left", annotation_text="stdev_" + str(i + 1))
    fig.update_yaxes(range=[0, 160])

    # fig.add_trace(go.Histogram(y=y, orientation='h', nbinsy=5000,), row=1, col=2)
    fig.add_trace(
        go.Histogram(y=np.concatenate(list(Y)).ravel().tolist(), orientation='h', marker_color='gray', name='Histogram',
                     nbinsy=8000), row=1, col=2)
    fig.update_layout(xaxis2=dict(range=[0, 200]))
    # fig.add_trace(go.Histogram(y=np.concatenate(list(Y[0:len(Y)//2-1])).ravel().tolist(), name='HP-1', orientation='h', nbinsy=8000, marker_color='#6A5ACD'), row=1, col=2)
    # fig.add_trace(go.Histogram(y=np.concatenate(list(Y[len(Y)//2:len(Y)])).ravel().tolist(), name='HP-2', orientation='h', nbinsy=8000, marker_color='#2E8B57'), row=1, col=2)

    # Overlay both histograms
    # fig.update_layout(barmode='overlay')
    # Reduce opacity to see both histograms
    # fig.update_traces(opacity=0.75)
    fig.update_layout(legend={'itemsizing': 'constant'})
    # Update layout
    fig.update_layout(
        template="plotly_white",
        font_family="Times New Roman"
    )
    fig.update_layout(
        title=chrom,
    )
    fig.update_layout(legend=dict(
        orientation = 'h', xanchor = "center", x = 0.5, y= 1.2 #orientation = 'v', xanchor = "center", x = 1.08, y= .5
    ))
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    fig.update_layout(
        width=680,
        height=400,
       )
    fig.update_layout(
        title={
            'text': chrom + ' - ' + arguments['genome_name'],
            'y': 0.96,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family="Courier New",
        font_color="dimgray",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="green",
    )



    print_chromosome_html(fig, chrom + '_clusters_'  + str(len(means)), html_graphs, arguments['out_dir_plots'])
    html_graphs.write("  <object data=\"" + chrom + '_clusters_'  + str(len(means))  + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

    #fig.write_html("coverage_plots/" + chrom + "_cluster_" + str(len(means)) + ".html")
