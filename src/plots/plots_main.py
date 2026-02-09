import numpy
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly
import pandas as pd
import numpy as np
import csv
import os
import shutil
import statistics
import itertools
import logging
import plotly.express as px

logger = logging.getLogger()

#from src.phasing_correction import phaseblock_flipping, phase_correction_centers, contiguous_phaseblocks, detect_centromeres, flip_phaseblocks_contigous, remove_overlaping_contiguous, switch_inter_phaseblocks_bins
from src.coverage.smoothing import smoothing
from src.file_tools.process_vcf_legacy import get_snps_frquncies_coverage, snps_mean, get_snp_segments
from src.file_tools.process_vcf import vcf_parse_to_csv_for_snps
#from src.utils import subclonal_values_adjusted, csv_df_chromosomes_sorter, df_chromosomes_sorter
from src.cna.copynumber import subclonal_values_adjusted
from src.utils.chromosome import csv_df_chromosomes_sorter, df_chromosomes_sorter, get_contigs_list
from src.cna.loh import detect_alter_loh_regions
from src.coverage.binning import get_chromosomes_regions, get_chromosomes_bins
from src.coverage.segmentation import change_point_detection_means
from src.cna.copynumber import bins_with_copynumber_states
from src.breakpoint.breakpoints import sv_vcf_bps_cn_check
from src.breakpoint.breakpoints_arcs import get_all_breakpoints_data
from src.utils.statistics import adjust_extreme_outliers

pd.options.mode.chained_assignment = None

def copy_number_plots_chromosomes(df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, args, loh_region_starts, loh_region_ends, is_half):
    filename = f"{os.path.join(args.out_dir_plots, 'COPY_NUMBERS.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")

    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        copy_number_plots_per_chromosome(df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, args, chrom, html_graphs, loh_region_starts, loh_region_ends, is_half)
    html_graphs.write("</body></html>")

def copy_number_plots_per_chromosome(centers, integer_fractional_means, ref_start_values, df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, args, chrom, html_graphs, loh_region_starts, loh_region_ends, p_value, is_half):

    logger.debug('Plots generation for ' + chrom)
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

    if not args.copynumbers_disable:
        #logger.info('copynumbers plots module')
        OFFSET=0
        #haplotype_1_values_cnr = list(np.asarray(haplotype_1_values_cnr) + OFFSET)
        #haplotype_2_values_cnr = list(np.asarray(haplotype_2_values_cnr) - OFFSET)

        if args.without_phasing:
            add_scatter_trace_coverage(fig, haplotype_1_start_values_cnr, haplotype_1_values_cnr, name='Unphased', text=ref_start_values, yaxis=None, opacity=0.7, color='olive')
            #plots_add_markers_lines(fig)
            plot_copynumbers_scatter_lines(args, fig,haplotype_1_values_copyratios, haplotype_1_start_values_copyratios, haplotype_1_end_values_copyratios, haplotype_1_values_copyratios, haplotype_1_start_values_copyratios, haplotype_1_end_values_copyratios, df_segs_hp1_chrom, df_segs_hp2_chrom, 0)
        else:
            haplotype_2_values_cnr = [-x for x in haplotype_2_values_cnr]
            add_scatter_trace_coverage(fig, haplotype_1_start_values_cnr, haplotype_1_values_cnr, name='HP-1', text=ref_start_values, yaxis=None, opacity=0.7, color='firebrick')
            add_scatter_trace_coverage(fig, haplotype_1_start_values_cnr, haplotype_2_values_cnr, name='HP-2', text=ref_start_values, yaxis=None, opacity=0.7, color='steelblue')
            #plots_add_markers_lines(fig)

            haplotype_2_values_copyratios = [-x for x in haplotype_2_values_copyratios]
            plot_copynumbers_scatter_lines(args, fig,haplotype_1_values_copyratios, haplotype_1_start_values_copyratios, haplotype_1_end_values_copyratios, haplotype_2_values_copyratios, haplotype_2_start_values_copyratios, haplotype_2_end_values_copyratios, df_segs_hp1_chrom, df_segs_hp2_chrom, OFFSET)

        if loh_region_starts:
            for k, (start_loh,end_loh) in enumerate(zip(loh_region_starts, loh_region_ends)):
                fig.add_vrect(x0=start_loh, x1=end_loh, fillcolor="#f1c40f", opacity=0.5, layer="below", line_width=0.5, )

        plots_layout_settings(fig, chrom, args, haplotype_1_end_values_copyratios[-1:][0], args.cut_threshold)
        if args.without_phasing:
            fig.update_yaxes(range = [0, args.cut_threshold])
        else:
            fig.update_yaxes(range=[-(args.cut_threshold), args.cut_threshold])
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

    if args.without_phasing:
        path = args.out_dir_plots + '/snps_loh_plots/'
    else:
        if is_half:
            path = args.out_dir_plots +'/wgd/'+ str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/variation_plots'
        else:
            path = args.out_dir_plots +'/'+ str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/variation_plots'
    if args.pdf_enable:
        print_chromosome_pdf(fig, chrom, path)

    print_chromosome_html(fig, chrom + '_cn', html_graphs, path)
    html_graphs.write("  <object data=\"" + chrom + '_cn' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

def plot_copynumbers_scatter_lines(args, fig,haplotype_1_values_copyratios, haplotype_1_start_values_copyratios, haplotype_1_end_values_copyratios, haplotype_2_values_copyratios, haplotype_2_start_values_copyratios, haplotype_2_end_values_copyratios, df_segs_hp1, df_segs_hp2, OFFSET):
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
    if args.without_phasing:
        add_scatter_trace_copyratios(args, fig, ['darkolivegreen'], name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values, df_segs_hp1, df_segs_hp2, mul_cols=False, row=1)
    else:
        add_scatter_trace_copyratios(args, fig, ['firebrick','steelblue'], name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values, df_segs_hp1, df_segs_hp2, mul_cols=False, row=1)

def coverage_plots_chromosomes(df, df_phasesets, args, thread_pool):
    if not os.path.isdir(args.out_dir_plots+'/coverage_plots'):
        os.mkdir(args.out_dir_plots+'/coverage_plots')
    if args.without_phasing:
        if not os.path.isdir(args.out_dir_plots+'/bed_output'):
            os.mkdir(args.out_dir_plots+'/bed_output')

    filename = f"{os.path.join(args.out_dir_plots, 'coverage_plots/COVERAGE_INDEX.html')}"
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

    if args.phaseblock_flipping_disable and not args.without_phasing:
        if os.path.exists(args.out_dir_plots + '/coverage_data'):
            shutil.rmtree(args.out_dir_plots + '/coverage_data')
            os.mkdir(args.out_dir_plots + '/coverage_data')
        else:
            os.mkdir(args.out_dir_plots + '/coverage_data')
        get_snp_segments(args, args.target_bam[0], thread_pool)
        df_snps = csv_df_chromosomes_sorter(args.out_dir_plots+'/data/snps_frequencies.csv', ['chr', 'pos', 'freq_value_a', 'hp_a', 'freq_value_b', 'hp_b'])
    #     #df_snps = df_snps.drop(df_snps[(df_snps.chr == "chrX") | (df_snps.chr == "chrY")].index)
    #     #df_snps = df_snps.drop(df_snps[(df_snps.chr == "chrY")].index)
    #     output_file_path_snps = vcf_parse_to_csv_for_snps(args.normal_phased_vcf, args)
    # else:
    #     output_file_path_snps = vcf_parse_to_csv_for_snps(args.tumor_phased_vcf, args)
    # snps_het_df = csv_df_chromosomes_sorter(output_file_path_snps, ['chr', 'pos', 'qual', 'gt', 'dp', 'vaf'])

    if args.breakpoints:
        _, _, _, breakpoints_segemnts, bps, bps_bnd = sv_vcf_bps_cn_check(args.breakpoints, args)
        df_var_bins, df_var_bins_1 = get_chromosomes_bins(args.target_bam[0], args.bin_size, args)

    chroms = get_contigs_list(args.contigs)
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)

    if args.quick_start:
        loh_path = args.quick_start_coverage_path + '/'
    else:
        loh_path = args.out_dir_plots + '/coverage_data/'
    if args.tumor_phased_vcf and os.path.exists(loh_path + args.genome_name + '_loh_segments.csv'):
        df_loh = csv_df_chromosomes_sorter(loh_path + args.genome_name + '_loh_segments.csv', ['chr', 'start', 'end', 'hp'])
    else:
        df_loh = pd.DataFrame(columns=['chr', 'start', 'end', 'hp'])

    haplotype_1_segs_dfs = [] #pd.DataFrame()
    haplotype_2_segs_dfs = [] #pd.DataFrame()
    het_snps_df_all = []

    for index, chrom in enumerate(chroms):
        if chrom in chroms:# and (chrom == '1' or chrom == '18'):# and (chrom == 'chr5' or chrom == 'chr16'):
            logger.debug('Plots generation for ' + chrom)
            fig = go.Figure()

            df_chrom = df[df['chr'] == chrom]
            df_centm_chrom = df_centm[df_centm['chr'] == chrom]
            df_loh_chrom = df_loh[df_loh['chr'] == chrom]

            if args.without_phasing:
                values = df_chrom.coverage.values.tolist()
                ref_start_values = df_chrom.start.values.tolist()
                ref_end_values = df_chrom.end.values.tolist()
            else:
                df_chrom_phasesets = df_phasesets[df_phasesets['chr'] == chrom]

                #unphased_reads_values = df_chrom.hp3.clip(upper=args.cut_threshold']).values.tolist()
                #haplotype_1_values = df_chrom.hp1.clip(upper=args.cut_threshold']).values.tolist()
                #haplotype_2_values = df_chrom.hp2.clip(upper=args.cut_threshold']).values.tolist()
                unphased_reads_values = df_chrom.hp3.values.tolist()
                haplotype_1_values = df_chrom.hp1.values.tolist()
                haplotype_2_values = df_chrom.hp2.values.tolist()
                ref_start_values = df_chrom.start.values.tolist()
                ref_end_values = df_chrom.end.values.tolist()

                haplotype_1_values_phasesets = df_chrom_phasesets.hp1.clip(upper=args.cut_threshold).values.tolist()
                haplotype_2_values_phasesets = df_chrom_phasesets.hp2.clip(upper=args.cut_threshold).values.tolist()
                ref_start_values_phasesets = df_chrom_phasesets.start.values.tolist()
                ref_end_values_phasesets = df_chrom_phasesets.end.values.tolist()
                ref_start_values_1 = []
            ################################################################################

            #debug var_bins
            if args.variable_size_bins:
                df_var_bins_chr = df_var_bins[df_var_bins['chr'] == chrom]
                ref_start_values = df_var_bins_chr.start.values.tolist()
                ref_end_values = df_var_bins_chr.end.values.tolist()

            if args.breakpoints:
                df_var_bins_chr_1 = df_var_bins_1[df_var_bins_1['chr'] == chrom]
                ref_start_values_1 = df_var_bins_chr_1.start.values.tolist()
                ref_end_values_1 = df_var_bins_chr_1.end.values.tolist()
            else:
                ref_start_values_1 = []
                breakpoints_segemnts = []
            ################################################################################
            if args.phaseblock_flipping_disable and not args.without_phasing:
                haplotype_1_values, haplotype_2_values = snps_mean(df_snps, ref_start_values, ref_end_values, chrom, args)
            #     # debug var_bins
            #     if args.variable_size_bins:
            #         unphased_reads_values = haplotype_1_values
            ################################################################################
            if not args.without_phasing:
                haplotype_1_values = adjust_extreme_outliers(haplotype_1_values)
                haplotype_2_values = adjust_extreme_outliers(haplotype_2_values)
            ################################################################################
            # if args.tumor_phased_vcf']:
            #     logger.info('hetrozygous phased snps frequencies coverage module')
            #     ref_start_values_updated, snps_het_counts, snps_homo_counts, centromere_region_starts, centromere_region_ends, loh_region_starts, loh_region_ends = get_snps_frquncies_coverage(df_snps_in_csv, chrom, ref_start_values, args.bin_size_snps'])
            #     if snps_het_counts or snps_homo_counts:
            #         if args.without_phasing']:
            #             if centromere_region_starts:
            #                 _,_,_, centromere_region_starts, centromere_region_ends  = detect_alter_loh_regions(args, 'centromere/no-coverage', chrom, ref_end_values, values, values, values, centromere_region_starts, centromere_region_ends, True)
            #             if loh_region_starts:
            #                 _,_,_, loh_region_starts, loh_region_ends = detect_alter_loh_regions(args, 'loss-of-heterozygosity', chrom, ref_end_values, values, values, values, loh_region_starts, loh_region_ends, True)
            #         else:
            #             if centromere_region_starts:
            #                 haplotype_1_values, haplotype_2_values, unphased_reads_values, centromere_region_starts, centromere_region_ends  = detect_alter_loh_regions(args, 'centromere/no-coverage', chrom, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values, centromere_region_starts, centromere_region_ends, True)
            #                 haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = loh_regions_phasesets(centromere_region_starts, centromere_region_ends, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
            #
            #             if loh_region_starts:
            #                 haplotype_1_values, haplotype_2_values, unphased_reads_values, loh_region_starts, loh_region_ends = detect_alter_loh_regions(args, 'loss-of-heterozygosity', chrom, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values, loh_region_starts, loh_region_ends, True)
            #                 haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = loh_regions_phasesets(loh_region_starts, loh_region_ends, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
            # else:
            #     loh_region_starts = []
            #     loh_region_ends = []
            loh_region_starts = []
            loh_region_ends = []
            ################################################################################
            #Raw coverage data plot before phase correction
            # if args.phaseblock_flipping_enable']:
            #     plot_coverage_raw(args, chrom, html_graphs, ref_start_values, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
            ################################################################################
            # if args.phaseblock_flipping_enable']:
            #     logger.info('phaseblock flipping module')
            #     is_simple_heuristics = True
            #     if len(ref_start_values_phasesets) >= 1:
            #         is_simple_heuristics = False
            #     haplotype_1_values, haplotype_2_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = \
            #     phaseblock_flipping(chrom, args, is_simple_heuristics, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
            #             haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
            #
            #     if not is_simple_heuristics:
            #         #detect centromeres
            #         ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = detect_centromeres(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values, ref_start_values, args.bin_size'])
            #         switch_inter_phaseblocks_bins(chrom, args, ref_start_values, haplotype_1_values, haplotype_2_values, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets)
            #
            #         # infer missing phaseblocks
            #         # ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = infer_missing_phaseblocks(ref_start_values, ref_end_values, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values, args.bin_size'])
            #         # create more contiguous phaseblocks
            #         haplotype_1_values_phasesets_conti, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, haplotype_2_values_phasesets_conti, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2 = contiguous_phaseblocks(haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
            #
            #         # flip phaseblocks based on more contiguous phaseblocks
            #         haplotype_1_values_phasesets, haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values = flip_phaseblocks_contigous(chrom, args, haplotype_1_values_phasesets_conti, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, haplotype_2_values_phasesets_conti, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2, ref_start_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values, haplotype_2_values)
            #
            ################################################################################
            if args.smoothing_enable:
                logger.info('smoothing module')
                if args.without_phasing:
                    values, values, unphased_reads_values = smoothing(values, values, values, conv_window_size=5)
                else:
                    haplotype_1_values, haplotype_2_values, unphased_reads_values = smoothing(haplotype_1_values, haplotype_2_values, unphased_reads_values, conv_window_size=5)
            ################################################################################
            #Plots
            if args.without_phasing:
                add_scatter_trace_coverage(fig, ref_start_values, values, name='Unphased', text=ref_start_values, yaxis=None, opacity=0.7, color='olive')
            else:
                add_scatter_trace_coverage(fig, ref_start_values, haplotype_1_values, name='HP-1', text=ref_start_values, yaxis=None, opacity=0.7, color='firebrick')
                add_scatter_trace_coverage(fig, ref_start_values, haplotype_2_values, name='HP-2', text=ref_start_values, yaxis=None, opacity=0.7, color='steelblue')

                if not args.unphased_reads_coverage_disable:
                    add_scatter_trace_coverage(fig, ref_start_values, unphased_reads_values, name='Unphased', text=ref_start_values, yaxis=None, opacity=0.7, color='olive', visibility=False)

            #plots_add_markers_lines(fig)
            ################################################################################
            if args.phaseblocks_enable:
                #logger.info('phaseblocks plots module')
                gaps_values = np.full(len(haplotype_1_values_phasesets), 'None')
                haplotype_1_phaseblocks_values = list(itertools.chain.from_iterable(zip(haplotype_1_values_phasesets, haplotype_1_values_phasesets, gaps_values)))
                haplotype_2_phaseblocks_values = list(itertools.chain.from_iterable(zip(haplotype_2_values_phasesets, haplotype_2_values_phasesets, gaps_values)))
                phaseblocks_positions = list(itertools.chain.from_iterable(zip(ref_start_values_phasesets, ref_end_values_phasesets, gaps_values)))

                add_scatter_trace_phaseblocks(fig, phaseblocks_positions, haplotype_1_phaseblocks_values, haplotype_2_phaseblocks_values)
            ################################################################################
            # if args.breakpoints_enable:
            #     breakpoints = get_breakpoints(chrom, args.breakpoints_file)
            #     add_scatter_trace_breakpoints(fig, breakpoints)
            regions = get_chromosomes_regions(args)
            chroms = get_contigs_list(args.contigs)
            chrom_index = chroms.index(chrom)
            plots_layout_settings(fig, chrom, args, regions[chrom_index], args.cut_threshold)
            ################################################################################
            if args.pdf_enable:
                print_chromosome_pdf(fig, chrom, args.out_dir_plots+'/coverage_plots')
            ################################################################################
            chr = range(len(ref_start_values))
            chr_all.extend([chrom for ch in chr])
            ref_start_values_all.extend(ref_start_values)
            ref_end_values_all.extend(ref_end_values)

            if args.without_phasing:
                df_snps_freqs_chr = whole_genome_combined_df(args, chrom, chr, ref_start_values, ref_end_values, values, values, values)
                # change point detection
                snps_cpd_means, snps_cpd_lens, df_means_chr = change_point_detection_means(args, chrom, breakpoints_segemnts, df_snps_freqs_chr, ref_start_values, ref_start_values_1, df_centm_chrom, df_loh_chrom)
                #df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, states, centers, stdev = apply_copynumbers(df_snps_freqs_chr, values, values, args, snps_cpd_means, [])
                snps_cpd_points.extend(snps_cpd_means)
                snps_cpd_points_weights.extend(snps_cpd_lens)
                #if args.enable_debug']:
                #    copy_number_plots_per_chromosome(centers, centers, ref_start_values, values, df_segs_hp1, values, df_segs_hp2, args, chrom, html_graphs, loh_region_starts, loh_region_ends)
                df_segs_hp1 = df_means_chr
                df_segs_hp2 = df_means_chr
                haplotype_1_segs_dfs.append(df_segs_hp1)
                haplotype_2_segs_dfs.append(df_segs_hp2)
                #if args.enable_debug:
                #    change_point_detection(values, ref_start_values, ref_end_values, args, chrom, html_graphs, 1, color='#6A5ACD')
            else:
                df_snps_freqs_chr = whole_genome_combined_df(args, chrom, chr, ref_start_values, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values)

                #change point detection
                snps_cpd_means, snps_cpd_lens, df_means_chr = change_point_detection_means(args, chrom, breakpoints_segemnts, df_snps_freqs_chr, ref_start_values, ref_start_values_1, df_centm_chrom, df_loh_chrom)
                #df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, states, centers, stdev = apply_copynumbers(df_snps_freqs_chr, haplotype_1_values, haplotype_2_values, args, snps_cpd_means, [])
                snps_cpd_points.extend(snps_cpd_means)
                snps_cpd_points_weights.extend(snps_cpd_lens)

                df_segs_hp1 = df_means_chr[0]
                df_segs_hp2 = df_means_chr[1]

                #TODO For debug - histo_clusters and cpds
                #add_histo_clusters_plot(df_cnr_hp1.log2.values.tolist(), df_cnr_hp2.log2.values.tolist(), states, centers, stdev, args, chrom, html_graphs)
                #if args.enable_debug:
                #    change_point_detection(haplotype_1_values, ref_start_values, ref_end_values, ref_start_values_1, args, chrom, html_graphs, 1, color='#6A5ACD')
                #    change_point_detection(haplotype_2_values, ref_start_values, ref_end_values, ref_start_values_1, args, chrom, html_graphs, 2, color='#2E8B57')
                #if args.enable_debug:
                #    copy_number_plots_per_chromosome(centers, centers, ref_start_values, haplotype_1_values, df_segs_hp1, haplotype_2_values, df_segs_hp2, args, chrom, html_graphs, loh_region_starts, loh_region_ends)

                #from gmm import process_gmm
                #process_gmm(haplotype_1_values, haplotype_2_values)

                haplotype_1_segs_dfs.append(df_segs_hp1)
                haplotype_2_segs_dfs.append(df_segs_hp2)

            if args.without_phasing:
                values_extended.extend(values)
            else:
                haplotype_1_values_updated.extend(haplotype_1_values)
                haplotype_2_values_updated.extend(haplotype_2_values)
                hunphased_updated.extend(unphased_reads_values)

            snps_cpd_means_all.extend(snps_cpd_means)
            if args.without_phasing:
                df_means_chr_all.append(df_means_chr)
            else:
                df_means_chr_all_hp1.append(df_means_chr[0])
                df_means_chr_all_hp2.append(df_means_chr[1])

            print_chromosome_html(fig, chrom + '_cov', html_graphs, args.out_dir_plots+'/coverage_plots/')
            html_graphs.write("  <object data=\"" + chrom + '_cov' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

    if args.without_phasing:
        df_snps_freqs = pd.DataFrame(list(zip(chr_all, ref_start_values_all, ref_end_values_all, values_extended)), columns=['chr', 'start', 'end', 'coverage'])
    else:
        df_snps_freqs = pd.DataFrame(list(zip(chr_all, ref_start_values_all, ref_end_values_all, haplotype_1_values_updated, haplotype_2_values_updated, hunphased_updated)), columns=['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

    if not args.without_phasing:
        df_means_chr_all_.append(pd.concat(df_means_chr_all_hp1))
        df_means_chr_all_.append(pd.concat(df_means_chr_all_hp2))

    html_graphs.write("</body></html>")

    if args.without_phasing:
        return values_extended, values_extended, values_extended, df_snps_freqs, snps_cpd_points, snps_cpd_points_weights, pd.concat(df_means_chr_all)
    else:
        return haplotype_1_values_updated, haplotype_2_values_updated, hunphased_updated, df_snps_freqs, snps_cpd_points, snps_cpd_points_weights, df_means_chr_all_

def breakpoints_segments_means(df, df_phasesets, args, thread_pool):

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

    if args.phaseblock_flipping_disable and not args.without_phasing:
        get_snp_segments(args, args.target_bam[0], thread_pool)
        df_snps = csv_df_chromosomes_sorter(args.out_dir_plots+'/data/snps_frequencies.csv', ['chr', 'pos', 'freq_value_a', 'hp_a', 'freq_value_b', 'hp_b'])

    if args.breakpoints:
        _, _, _, breakpoints_segemnts, bps, bps_bnd = sv_vcf_bps_cn_check(args.breakpoints, args)
        df_var_bins, df_var_bins_1 = get_chromosomes_bins(args.target_bam[0], args.bin_size, args)

    chroms = get_contigs_list(args.contigs)
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)

    if args.quick_start:
        loh_path = args.quick_start_coverage_path + '/'
    else:
        loh_path = args.out_dir_plots + '/coverage_data/'
    if args.tumor_phased_vcf and os.path.exists(loh_path + args.genome_name + '_loh_segments.csv'):
        df_loh = csv_df_chromosomes_sorter(loh_path + args.genome_name + '_loh_segments.csv', ['chr', 'start', 'end', 'hp'])
    else:
        df_loh = pd.DataFrame(columns=['chr', 'start', 'end', 'hp'])

    haplotype_1_segs_dfs = [] #pd.DataFrame()
    haplotype_2_segs_dfs = [] #pd.DataFrame()
    het_snps_df_all = []

    for index, chrom in enumerate(chroms):
        if chrom in chroms:# and (chrom == '1' or chrom == '18'):# and (chrom == 'chr5' or chrom == 'chr16'):
            logger.debug('Mean segmentation for ' + chrom)
            fig = go.Figure()

            df_chrom = df[df['chr'] == chrom]
            df_centm_chrom = df_centm[df_centm['chr'] == chrom]
            df_loh_chrom = df_loh[df_loh['chr'] == chrom]

            if args.without_phasing:
                values = df_chrom.coverage.values.tolist()
                ref_start_values = df_chrom.start.values.tolist()
                ref_end_values = df_chrom.end.values.tolist()
            else:
                df_chrom_phasesets = df_phasesets[df_phasesets['chr'] == chrom]
                unphased_reads_values = df_chrom.hp3.values.tolist()
                haplotype_1_values = df_chrom.hp1.values.tolist()
                haplotype_2_values = df_chrom.hp2.values.tolist()
                ref_start_values = df_chrom.start.values.tolist()
                ref_end_values = df_chrom.end.values.tolist()

                haplotype_1_values_phasesets = df_chrom_phasesets.hp1.clip(upper=args.cut_threshold).values.tolist()
                haplotype_2_values_phasesets = df_chrom_phasesets.hp2.clip(upper=args.cut_threshold).values.tolist()
                ref_start_values_phasesets = df_chrom_phasesets.start.values.tolist()
                ref_end_values_phasesets = df_chrom_phasesets.end.values.tolist()
                ref_start_values_1 = []
            ################################################################################

            #debug var_bins
            if args.variable_size_bins:
                df_var_bins_chr = df_var_bins[df_var_bins['chr'] == chrom]
                ref_start_values = df_var_bins_chr.start.values.tolist()
                ref_end_values = df_var_bins_chr.end.values.tolist()

            if args.breakpoints:
                df_var_bins_chr_1 = df_var_bins_1[df_var_bins_1['chr'] == chrom]
                ref_start_values_1 = df_var_bins_chr_1.start.values.tolist()
                ref_end_values_1 = df_var_bins_chr_1.end.values.tolist()
            else:
                ref_start_values_1 = []
                breakpoints_segemnts = []
            ################################################################################
            if args.phaseblock_flipping_disable and not args.without_phasing:
                haplotype_1_values, haplotype_2_values = snps_mean(df_snps, ref_start_values, ref_end_values, chrom, args)
            ################################################################################
            if not args.without_phasing:
                haplotype_1_values = adjust_extreme_outliers(haplotype_1_values)
                haplotype_2_values = adjust_extreme_outliers(haplotype_2_values)
            ################################################################################
            loh_region_starts = []
            loh_region_ends = []
            ################################################################################
            ################################################################################
            if args.phaseblocks_enable:
                gaps_values = np.full(len(haplotype_1_values_phasesets), 'None')
                haplotype_1_phaseblocks_values = list(itertools.chain.from_iterable(zip(haplotype_1_values_phasesets, haplotype_1_values_phasesets, gaps_values)))
                haplotype_2_phaseblocks_values = list(itertools.chain.from_iterable(zip(haplotype_2_values_phasesets, haplotype_2_values_phasesets, gaps_values)))
                phaseblocks_positions = list(itertools.chain.from_iterable(zip(ref_start_values_phasesets, ref_end_values_phasesets, gaps_values)))

            regions = get_chromosomes_regions(args)
            chroms = get_contigs_list(args.contigs)
            chrom_index = chroms.index(chrom)
            chr = range(len(ref_start_values))
            chr_all.extend([chrom for ch in chr])
            ref_start_values_all.extend(ref_start_values)
            ref_end_values_all.extend(ref_end_values)

            if args.without_phasing:
                df_snps_freqs_chr = whole_genome_combined_df(args, chrom, chr, ref_start_values, ref_end_values, values, values, values)
                snps_cpd_means, snps_cpd_lens, df_means_chr = change_point_detection_means(args, chrom, breakpoints_segemnts, df_snps_freqs_chr, ref_start_values, ref_start_values_1, df_centm_chrom, df_loh_chrom)
                snps_cpd_points.extend(snps_cpd_means)
                snps_cpd_points_weights.extend(snps_cpd_lens)
                df_segs_hp1 = df_means_chr
                df_segs_hp2 = df_means_chr
                haplotype_1_segs_dfs.append(df_segs_hp1)
                haplotype_2_segs_dfs.append(df_segs_hp2)
            else:
                df_snps_freqs_chr = whole_genome_combined_df(args, chrom, chr, ref_start_values, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values)

                #change point detection
                snps_cpd_means, snps_cpd_lens, df_means_chr = change_point_detection_means(args, chrom, breakpoints_segemnts, df_snps_freqs_chr, ref_start_values, ref_start_values_1, df_centm_chrom, df_loh_chrom)
                #df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, states, centers, stdev = apply_copynumbers(df_snps_freqs_chr, haplotype_1_values, haplotype_2_values, args, snps_cpd_means, [])
                snps_cpd_points.extend(snps_cpd_means)
                snps_cpd_points_weights.extend(snps_cpd_lens)

                df_segs_hp1 = df_means_chr[0]
                df_segs_hp2 = df_means_chr[1]

                haplotype_1_segs_dfs.append(df_segs_hp1)
                haplotype_2_segs_dfs.append(df_segs_hp2)

            if args.without_phasing:
                values_extended.extend(values)
            else:
                haplotype_1_values_updated.extend(haplotype_1_values)
                haplotype_2_values_updated.extend(haplotype_2_values)
                hunphased_updated.extend(unphased_reads_values)

            snps_cpd_means_all.extend(snps_cpd_means)
            if args.without_phasing:
                df_means_chr_all.append(df_means_chr)
            else:
                df_means_chr_all_hp1.append(df_means_chr[0])
                df_means_chr_all_hp2.append(df_means_chr[1])

    if args.without_phasing:
        df_snps_freqs = pd.DataFrame(list(zip(chr_all, ref_start_values_all, ref_end_values_all, values_extended)), columns=['chr', 'start', 'end', 'coverage'])
    else:
        df_snps_freqs = pd.DataFrame(list(zip(chr_all, ref_start_values_all, ref_end_values_all, haplotype_1_values_updated, haplotype_2_values_updated, hunphased_updated)), columns=['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

    if not args.without_phasing:
        df_means_chr_all_.append(pd.concat(df_means_chr_all_hp1))
        df_means_chr_all_.append(pd.concat(df_means_chr_all_hp2))

    if args.without_phasing:
        return values_extended, values_extended, values_extended, df_snps_freqs, snps_cpd_points, snps_cpd_points_weights, pd.concat(df_means_chr_all)
    else:
        return haplotype_1_values_updated, haplotype_2_values_updated, hunphased_updated, df_snps_freqs, snps_cpd_points, snps_cpd_points_weights, df_means_chr_all_


def plot_coverage_raw(args, chrom, html_graphs, ref_start_values, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):
    fig = go.Figure()

    add_scatter_trace_coverage(fig, ref_start_values, haplotype_1_values, name='HP-1', text=None, yaxis=None, opacity=0.7, color='firebrick')
    add_scatter_trace_coverage(fig, ref_start_values, haplotype_2_values, name='HP-2', text=None, yaxis=None, opacity=0.7, color='steelblue')

    if not args.unphased_reads_coverage_disable:
        add_scatter_trace_coverage(fig, ref_start_values, unphased_reads_values, name='Unphased', text=None, yaxis=None, opacity=0.7, color='olive')
    #plots_add_markers_lines(fig)

    if args.phaseblocks_enable:
        #logger.info('phaseblocks plots module')
        gaps_values = np.full(len(haplotype_1_values_phasesets), 'None')
        haplotype_1_phaseblocks_values = list(itertools.chain.from_iterable(zip(haplotype_1_values_phasesets, haplotype_1_values_phasesets, gaps_values)))
        haplotype_2_phaseblocks_values = list(itertools.chain.from_iterable(zip(haplotype_2_values_phasesets, haplotype_2_values_phasesets, gaps_values)))
        phaseblocks_positions = list(
            itertools.chain.from_iterable(zip(ref_start_values_phasesets, ref_end_values_phasesets, gaps_values)))

        add_scatter_trace_phaseblocks(fig, phaseblocks_positions, haplotype_1_phaseblocks_values, haplotype_2_phaseblocks_values)
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)
    chrom_index = chroms.index(chrom)
    plots_layout_settings(fig, chrom, args, regions[chrom_index], args.cut_threshold)

    print_chromosome_html(fig, chrom + '_cov_raw', html_graphs, args.out_dir_plots + '/coverage_plots/')
    html_graphs.write(
        "  <object data=\"" + chrom + '_cov_raw' + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")


def whole_genome_combined_df(args, chrom, chr, ref_start_values, ref_end_values, snps_haplotype1_mean, snps_haplotype2_mean, snps_haplotype_mean):
    if args.without_phasing:
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


def plots_genome_coverage(df_hp1, df_hp2, df_unphased, args):
    # TODO genome-wide plots

    lengths = []
    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_hp1[df_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_hp1)*args.bin_size, args.bin_size, dtype=int)

    fig = go.Figure()
    add_scatter_trace_coverage(fig, indices, df_hp1.hp1.values.tolist(), name='HP-1', text=None,
                               yaxis=None, opacity=0.7, color='firebrick')
    add_scatter_trace_coverage(fig, indices, df_hp2.hp2.values.tolist(), name='HP-2', text=None,
                               yaxis=None, opacity=0.7, color='steelblue')

    if not args.unphased_reads_coverage_disable:
        add_scatter_trace_coverage(fig, indices, df_unphased.hp3.values.tolist(), name='Unphased reads', text=None,
                                   yaxis=None, opacity=0.5, color='olive')
    #plots_add_markers_lines(fig)

    current = 0
    label_pos = []
    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        current += lengths[index]
        label_pos.append(round(current - (lengths[index] / 2))*args.bin_size)
        fig.add_vline(x=current*args.bin_size, y0=-10, y1=150, line_width=1, line_dash="dashdot", line_color="green")

    fig.update_layout(
        xaxis=dict(
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        ),
        font=dict(size=18, color="black"))

    plots_layout_settings(fig, 'Genome', args, indices[-1:][0], args.cut_threshold)
    #fig.update_yaxes(range=[-10, 10])
    #fig.update_yaxes(title='copy ratio (log2)')
    fig.update_layout(width=1280, height=600,)
    fig.update_layout(legend=dict(orientation='h', xanchor="center", x=0.5, y=1.1))

    if args.pdf_enable:
        print_genome_pdf(fig, args.genome_name, args.out_dir_plots)

    fig.write_html(args.out_dir_plots +'/'+ args.genome_name + "_genome_coverage.html")

def copy_number_plots_genome_details(centers, integer_fractional_centers, df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, df_unphased, args, x_axis, observed_hist, p_value, is_half):
    # TODO genome-wide plots

    lengths = []
    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    regions = get_chromosomes_regions(args)
    indices = np.arange(0, len(df_cnr_hp1)*args.bin_size, args.bin_size, dtype=int)

    df_cnr_hp1_1 = df_cnr_hp1
    df_1_ = []
    offset = 0
    for i, chrom in enumerate(chroms):
        df_genes_chrom_1 = df_cnr_hp1_1[df_cnr_hp1_1['chr'] == chrom]
        df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
        offset += regions[i]
        df_1_.append(df_genes_chrom_1)
    df_1 = pd.concat(df_1_)

    #fig = go.Figure()
    fig = make_subplots(rows=1, cols=2, shared_yaxes=False, column_widths=[0.80, 0.10], vertical_spacing=0.01,
                        horizontal_spacing=0.08, specs=[[{"secondary_y":True}, {"type":"xy"}]])

    ###########################################################
    # from breakpoints_arcs import get_all_breakpoints_data
    # edges = [[0, 1450001//50000]]
    # interact_strength = [1]
    #
    # labels =  ['0', '1450001', '143700001']#[str(i) for i in indices.tolist()]
    # print(len(labels))
    # L = len(labels)
    #
    # arcs_data = get_all_breakpoints_data(edges, L, interact_strength)
    # fig = go.Figure(data=arcs_data)
    ###########################################################

    custom_text_data_hp1 = df_cnr_hp1.start.values.tolist()
    custom_text_data_hp2 = df_cnr_hp2.start.values.tolist()
    if args.without_phasing:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.coverage.values.tolist(), name='Unphased', text=custom_text_data_hp1,
                                   yaxis=None, opacity=0.7, color='olive', visibility='legendonly', mul_cols=True, row=1)
    else:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.hp1.values.tolist() , name='HP-1', text=custom_text_data_hp1,
                                   yaxis=None, opacity=0.7, color='firebrick', visibility='legendonly', mul_cols=True, row=1)
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp2.hp2.values.tolist(), name='HP-2', text=custom_text_data_hp2,
                                   yaxis=None, opacity=0.7, color='steelblue', visibility='legendonly', mul_cols=True, row=1)
        if not args.unphased_reads_coverage_disable:
            add_scatter_trace_coverage(fig, df_1['start'], df_unphased.hp3.values.tolist(), name='Unphased', text=None,
                                   yaxis=None, opacity=0.7, color='olive', visibility='legendonly', mul_cols=True, row=1)

    #plots_add_markers_lines(fig)

    current = 0
    label_pos = []
    offset_chroms = 0
    offset_chroms_1 = 0
    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        current += lengths[index]
        label_pos.append(round(offset_chroms+regions[index]//2)) #y0=-10, y1=args.cut_threshold
        fig.add_vline(x=offset_chroms,  line_width=1, line_dash="solid", line_color="#D7DBDD")
        offset_chroms_1 += regions[index]

        if index % 2 == 0:
            fig.add_vrect(x0=offset_chroms, x1=offset_chroms_1, fillcolor="#E5E7E9", opacity=0.9, layer="below", line_width=0, )
        offset_chroms += regions[index]

    fig.update_layout(
        xaxis=dict(
            tickangle=90,
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        ),
        font=dict(size=18, color="black"))

    if not args.copynumbers_disable:
        offset_start = 0
        offset_end = 0
        haplotype_1_start_values = []
        haplotype_1_end_values = []

        haplotype_2_start_values = []
        haplotype_2_end_values = []
        chroms = get_contigs_list(args.contigs)
        for index, chrom in enumerate(chroms):
            df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
            df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
            haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
            haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()
            haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
            haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()

            if chrom == args.contigs.split('-')[0]:
                haplotype_1_start_values.extend(haplotype_1_start_values_copyratios)
                haplotype_1_end_values.extend(haplotype_1_end_values_copyratios)
                haplotype_2_start_values.extend(haplotype_2_start_values_copyratios)
                haplotype_2_end_values.extend(haplotype_2_end_values_copyratios)
            else:
                offset_start += regions[index-1]
                #lengths[index-1] * args.bin_size
                haplotype_1_start_values.extend([x + offset_start for x in haplotype_1_start_values_copyratios])
                haplotype_1_end_values.extend([x + offset_start for x in haplotype_1_end_values_copyratios])
                haplotype_2_start_values.extend([x + offset_start for x in haplotype_2_start_values_copyratios])
                haplotype_2_end_values.extend([x + offset_start for x in haplotype_2_end_values_copyratios])
                #offset_start += regions[index-1]

        haplotype_1_gaps_values = np.full(len(df_segs_hp1.state.values.tolist()), 'None')
        haplotype_1_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp1.state.values.tolist(), df_segs_hp1.state.values.tolist(), haplotype_1_gaps_values)))
        haplotype_1_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_1_start_values, haplotype_1_end_values, haplotype_1_gaps_values)))

        haplotype_2_gaps_values = np.full(len(df_segs_hp2.state.values.tolist()), 'None')
        haplotype_2_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp2.state.values.tolist(), df_segs_hp2.state.values.tolist(), haplotype_2_gaps_values)))
        haplotype_2_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_2_start_values, haplotype_2_end_values, haplotype_2_gaps_values)))


        if args.without_phasing:
            OFFSET = 0
            colors = ['darkolivegreen']
        else:
            OFFSET = args.cut_threshold/150
            colors = ['firebrick', 'steelblue']
        haplotype_1_copyratios_values = [x if x == 'None' else x  for x in haplotype_1_copyratios_values]
        haplotype_2_copyratios_values = [x if x == 'None' else x  for x in haplotype_2_copyratios_values]
        name = "Copynumbers"

        add_scatter_trace_copyratios(args, fig, colors, name, haplotype_1_copyratios_positions,
                                     haplotype_2_copyratios_positions, haplotype_1_copyratios_values,
                                     haplotype_2_copyratios_values, df_segs_hp1, df_segs_hp2, mul_cols=True, row=1)

    #############################################################
    # depth_values_hp1 = np.array(df_cnr_hp1.log2.values.tolist(), dtype='int')
    # depth_values_hp2 = np.array(df_cnr_hp2.log2.values.tolist(), dtype='int')
    #
    # depth_values_hp1 = depth_values_hp1.reshape(-1, 1)
    # depth_values_hp2 = depth_values_hp2.reshape(-1, 1)
    # Y = np.concatenate([depth_values_hp1, depth_values_hp2])
    #fig.add_trace(
    #    go.Histogram(y=observed_hist, orientation='h', marker_color='gray', name='Histogram',
    #                  visible="legendonly"), row=1, col=2)
    #fig.update_layout(xaxis2=dict(range=[0, 2000], showticklabels=False))
    fig.add_trace(go.Bar(x=observed_hist, y=x_axis, marker_color='grey', orientation='h',  name='Histogram'), row=1, col=2)
    #fig.add_trace(go.Scatter(x=fit_curve, y=x_axis, marker_color='royalblue', orientation='h', name='opt fit', line=dict(color='royalblue', width=4, dash='dot'),), row=1, col=2)

    #############################################################
    plots_layout_settings(fig, 'Genome', args, len(df_cnr_hp1.start.values.tolist())*args.bin_size, args.cut_threshold)

    centers_rev = [-x for x in centers[1:]]
    centers_rev.reverse()
    tick_vals = centers

    integer_fractional_means_rev = [x for x in integer_fractional_centers[1:]]
    integer_fractional_means_rev.reverse()
    tickt_ext = integer_fractional_centers
    if args.without_phasing:
        fig.update_yaxes(range=[-1, args.cut_threshold])
    else:
        fig.update_yaxes(range=[-1, args.cut_threshold])
    fig.update_layout(
        yaxis=dict(
            tickmode='array',
            tickvals=[i for i in range(0, 1000, 25)],
            ticktext=[str(abs(i)) for i in range(0, 1000, 25)]
        ),
        yaxis2=dict(
            tickmode='array',
            tickvals=tick_vals,
            # [i for i in range(-1000, 1000, 25)],#sorted(integer_centers + fractional_centers),#[0, 6,12,18], #13,26,39,5
            # [-99, -33, 0, 33, 99],#[i for i in range(-1000, 1000, 50)],
            ticktext=tickt_ext  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
            # [str(abs(i)) for i in range(-1000, 1000, 50)]
        )
    )
    #fig.update_yaxes(title='coverge (mean depth)')
    #fig.update_yaxes(title_text="<b>secondary</b> yaxis title", secondary_y=True)
    fig.update_layout(width=1680, height=600,)
    fig.update_layout(legend=dict(orientation = 'h', xanchor = "center", x = 0.5, y= 1.1))
    fig.update_layout(yaxis2={'side': 'right'} )

    if args.without_phasing:
        path_set = args.genome_name + '_genome_copynumbers_details'
    else:
        if is_half:
            path_set = 'wgd/' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) +  '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers_details'
        else:
            path_set = str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) +  '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers_details'

    if args.pdf_enable:
        print_genome_pdf(fig, path_set, args.out_dir_plots)

    fig.write_html(args.out_dir_plots + '/' + path_set + ".html")

    return centers, integer_fractional_centers

def add_annotation(fig, x, y, ax, ay, text, color):
    fig.add_annotation(
        x=x,
        y=y,
        xref="x",
        yref="y",
        text=text,
        showarrow=False,
        font=dict(
            family="Courier New, monospace",
            size=16,
            color="#ffffff"
        ),
        align="center",
        ax=ax,
        ay=ay,
        bordercolor="#c7c7c7",
        borderwidth=2,
        borderpad=4,
        bgcolor=color,
        opacity=0.8
    )

def copy_number_plots_genome_breakpoints_cytos(centers, integer_fractional_centers, df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, df_unphased, args, p_value, loh_regions, het_snps_bafs_means, is_half):
    # TODO genome-wide plots

    if args.breakpoints:
        coords_single, coords_single_chr, coords, coords_chr, _,_ = sv_vcf_bps_cn_check(args.breakpoints, args)

    lengths = []
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)

    #chroms = ['chr1', 'chr5','chr8']

    df_cnr_hp1_ = []
    df_segs_hp1_ = []
    df_cnr_hp2_ = []
    df_segs_hp2_ = []
    df_genes_1_ = []
    df_genes_2_ = []
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    #cancer_genes = os.path.join(fileDir, args.cancer_genes)
    df_genes = csv_df_chromosomes_sorter(args.cancer_genes, ['chr','start','end','gene'])
    genestart_1 = []
    genestart_2 = []
    last_len = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        df_segs_hp1_.append(df_segs_hp1[df_segs_hp1['chromosome'] == chrom])
        df_cnr_hp2_.append(df_cnr_hp2[df_cnr_hp2['chr'] == chrom])
        df_segs_hp2_.append(df_segs_hp2[df_segs_hp2['chromosome'] == chrom])

    df_cnr_hp1_1 = df_cnr_hp1
    df_1_ = []
    offset = 0
    for i, chrom in enumerate(chroms):
        df_genes_chrom_1 = df_cnr_hp1_1[df_cnr_hp1_1['chr'] == chrom]
        df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
        offset += regions[i]
        df_1_.append(df_genes_chrom_1)
    df_1 = pd.concat(df_1_)

    if args.breakpoints:
        offset = 0
        for i, chrom in enumerate(chroms):
            new_len = len(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
            if i > 0:
                for j, pair in enumerate(coords):
                    if chrom == pair[0]:
                        coords[j][1] = pair[1] + offset

            if i > 0:
                for k, pair in enumerate(coords_single):
                    if chrom == pair[0]:
                        coords_single[k][1] = pair[1] + offset

            last_len += new_len
            offset += regions[i]

        n = len(coords)
        coords = [x[1:] for x in coords]
        coords = [coords[i] + (coords[i + 1] if i + 1 < n else []) for i in range(0, n, 2)]
        coords = list(map(sorted, coords))
        last_len = 0
        coords_single = [j for i in [x[1:2] for x in coords_single] for j in i]

        arcs_data = get_all_breakpoints_data(coords, coords_chr, 75, args.breakpoints, args)
        arcs_data1 = get_all_breakpoints_data(coords, coords_chr, -75, args.breakpoints, args)

    df_cnr_hp1 = pd.concat(df_cnr_hp1_)
    df_segs_hp1 = pd.concat(df_segs_hp1_)

    df_cnr_hp2 = pd.concat(df_cnr_hp2_)
    df_segs_hp2 = pd.concat(df_segs_hp2_)

    df_genes_1 = df_genes[df_genes['start'] < df_genes['end']]
    df_genes_2 = df_genes[df_genes['start'] > df_genes['end']]

    last_len = 0
    offset = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        new_len = len(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        if not chrom.startswith('chr'):
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == 'chr' + chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == 'chr' + chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        else:
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        last_len += new_len
        offset += regions[i]

    df_genes_1 = pd.concat(df_genes_1_)
    df_genes_2 = pd.concat(df_genes_2_)

    genename_1 = df_genes_1['gene'].values.tolist()
    genename_2 = df_genes_2['gene'].values.tolist()

    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_cnr_hp1)*args.bin_size, args.bin_size, dtype=int)

    ###########################################################
    if args.without_phasing:
        row_heights = [220, 160, 150, 20, 40]
    else:
        row_heights = [220, 320, 20, 150, 40]

    #fig = make_subplots(rows=2, cols=1, shared_yaxes=False, shared_xaxes=True, vertical_spacing=0.02, horizontal_spacing=0.02)
    fig = make_subplots(rows=5, cols=1, shared_yaxes=False, shared_xaxes='columns',  vertical_spacing=0.01, row_heights=row_heights,
                        horizontal_spacing=0.02, specs=[[{"type":"xy"}], [{"secondary_y":True}], [{"type":"xy"}], [{"type":"xy"}], [{"type":"xy"}]]) #https://community.plotly.com/t/can-subplot-support-multiple-y-axes/38891/20
    # #############################################################
    if args.breakpoints:
        for i in range(len(arcs_data)):
            fig.add_trace(go.Scatter(x=arcs_data[i][0],
                                     y=arcs_data[i][1],
                                     name=arcs_data[i][2],
                                     mode=arcs_data[i][3],
                                     line=arcs_data[i][4],
                                     yaxis="y",
                                     opacity=0.5,
                                     text = arcs_data[i][5],
                                     customdata = [arcs_data[i][6]],
                                     hovertemplate= '<br><b>Breakpoint info: </b>' + arcs_data[i][6],
                                     showlegend=arcs_data[i][7]), row=1, col=1,
                         )
    # #############################################################

    # #############################################################
    custom_text_data_hp1 = df_cnr_hp1.start.values.tolist()
    custom_text_data_hp2 = df_cnr_hp2.start.values.tolist()
    if args.without_phasing:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.coverage.values.tolist(), name='Unphased coverage',
                                   text=custom_text_data_hp1,
                                   yaxis="y2", opacity=0.7, color='olive', visibility=True, mul_cols=True)
    else:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.hp1.values.tolist(), name='HP-1 coverage',
                                   text=custom_text_data_hp1,
                                   yaxis=None, opacity=0.7, color='firebrick', visibility=True, mul_cols=True)
        add_scatter_trace_coverage(fig, df_1['start'], [-x for x in df_cnr_hp2.hp2.values.tolist()], name='HP-2 coverage',
                                   text=custom_text_data_hp2,
                                   yaxis=None, opacity=0.7, color='steelblue', visibility=True, mul_cols=True)

    add_scatter_trace_coverage(fig, df_1['start'], het_snps_bafs_means.vaf.values.tolist(), name='BAF',
                                   text=het_snps_bafs_means.pos.values.tolist(),
                                   yaxis="y5", opacity=0.7, color='olive', visibility='legendonly', mul_cols=True, row=4, baf=True)
        # if not args.unphased_reads_coverage_disable:
        #     add_scatter_trace_coverage(fig, indices, df_unphased.hp3.values.tolist(), name='Unphased', text=None,
        #                                yaxis=None, opacity=0.7, color='olive', visibility='legendonly', mul_cols=True)

    fig.add_trace(go.Scatter(x=df_genes_1['start'], y=[0.5]*len(df_genes_1),#[1, 2, 3, 4, 5] * (len(df_genes_1)//5),
                             mode='markers',
                             text = genename_1,
                             customdata = genestart_1,
                             #hovertext=df_genes['gene'],
                             hovertemplate=
                             '<br><b>Gene</b>: %{text}'+
                             '<br><b>Pos</b>: %{customdata}<br>',
                             marker=dict(
                                 symbol="y-left",
                                 color="#3A6B35",
                                 size=6,
                                 line=dict(width=1, color="#7F7F7F"),
                             ),
                             yaxis="y6",
                             name= 'GeneInfo',
                             showlegend=False,), row=5, col=1,)

    fig.add_hline(y=2,  line_width=1, line_dash="solid", line_color="black", row=3, col=1,)
    fig.add_hline(y=-(args.cut_threshold + 5),  line_width=1, line_dash="solid", line_color="black", row=2, col=1,)

    #cytobands
    fileDir = os.path.dirname(__file__)
    if '13' in args.reference_name:
        cytos_file = os.path.join(fileDir, "annotations/chm13.cytobands.bed")
    else:
        cytos_file = os.path.join(fileDir, "annotations/grch38.cytobands.bed")
    df = pd.read_csv(cytos_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'gieStain'])
    #df = df_chromosomes_sorter(df, names=['chr', 'start', 'end', 'name', 'gieStain'])
    plot_cytobands_genome_wide(fig, df)

    current = 0
    label_pos = []
    offset_chroms = 0
    offset_chroms_1 = 0
    chroms = get_contigs_list(args.contigs)
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])

    for index, chrom in enumerate(chroms):
        loh_starts = []
        if not loh_regions.empty:
            df_loh_region = loh_regions[loh_regions['chr'] == chrom]
            loh_starts = df_loh_region.start.values.tolist()
            loh_ends = df_loh_region.end.values.tolist()
        if not df_centm.empty:
            df_cent_region = df_centm[df_centm['chr'] == chrom]
            cent_starts = df_cent_region.start.values.tolist()
            cent_ends = df_cent_region.end.values.tolist()

        current += lengths[index]
        label_pos.append(round(offset_chroms+regions[index]//2))

        fig.add_vline(x=offset_chroms,  line_width=1, line_dash="solid", line_color="#D7DBDD")

        offset_chroms_1 += regions[index]

        if index % 2 == 0:
            fig.add_vrect(x0=offset_chroms, x1=offset_chroms_1, fillcolor="#E5E7E9", opacity=0.9, layer="below", line_width=0, )

        if len(loh_starts) and args.loh_enable:
            for i in range(len(df_loh_region.start.values.tolist())):
                fig.add_vrect(x0=offset_chroms+loh_starts[i], x1=offset_chroms+loh_ends[i], fillcolor="#2980b9", opacity=0.3, layer="below", line_width=0, row=2, col=1,)
        if len(cent_starts) and args.loh_enable:
            fig.add_vrect(x0=offset_chroms+cent_starts[0], x1=offset_chroms+cent_ends[0], fillcolor="#7e1f14", opacity=0.3, layer="above", line_width=0, row=2, col=1,)

        offset_chroms += regions[index]

    if not args.copynumbers_disable:
        offset_start = 0
        haplotype_1_start_values = []
        haplotype_1_end_values = []

        haplotype_2_start_values = []
        haplotype_2_end_values = []
        for index, chrom in enumerate(chroms):
            df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
            df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
            haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
            haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()
            haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
            haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()
            if chrom == args.contigs.split('-')[0]:
                haplotype_1_start_values.extend(haplotype_1_start_values_copyratios)
                haplotype_1_end_values.extend(haplotype_1_end_values_copyratios)
                haplotype_2_start_values.extend(haplotype_2_start_values_copyratios)
                haplotype_2_end_values.extend(haplotype_2_end_values_copyratios)
            else:
                offset_start += regions[index-1]#lengths[index-1] * args.bin_size
                haplotype_1_start_values.extend([x + offset_start for x in haplotype_1_start_values_copyratios])
                haplotype_1_end_values.extend([x + offset_start for x in haplotype_1_end_values_copyratios])
                haplotype_2_start_values.extend([x + offset_start for x in haplotype_2_start_values_copyratios])
                haplotype_2_end_values.extend([x + offset_start for x in haplotype_2_end_values_copyratios])

        haplotype_1_gaps_values = np.full(len(df_segs_hp1.state.values.tolist()), 'None')
        haplotype_1_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp1.state.values.tolist(), df_segs_hp1.state.values.tolist(), haplotype_1_gaps_values)))
        haplotype_1_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_1_start_values, haplotype_1_end_values, haplotype_1_gaps_values)))

        haplotype_2_gaps_values = np.full(len(df_segs_hp2.state.values.tolist()), 'None')
        haplotype_2_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp2.state.values.tolist(), df_segs_hp2.state.values.tolist(), haplotype_2_gaps_values)))
        haplotype_2_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_2_start_values, haplotype_2_end_values, haplotype_2_gaps_values)))


        if args.without_phasing:
            OFFSET = 0
            colors = ['darkolivegreen']
        else:
            OFFSET = args.cut_threshold/150
            colors = ['firebrick', 'steelblue']
        haplotype_1_copyratios_values = [x if x == 'None' else x  for x in haplotype_1_copyratios_values]
        haplotype_2_copyratios_values = [x if x == 'None' else -(x) for x in haplotype_2_copyratios_values]
        name = "Copynumbers"

        add_scatter_trace_copyratios(args, fig, colors, name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values, df_segs_hp1, df_segs_hp2, mul_cols=True)

    centers_rev = [-x for x in centers]
    centers_rev.reverse()
    integer_fractional_means_rev = [x for x in integer_fractional_centers]
    integer_fractional_means_rev.reverse()
    if args.without_phasing:
        tick_vals = centers
        tickt_ext = integer_fractional_centers
        tickvals = [i for i in range(0, 1000, 25)]
        ticktext = [str(abs(i)) for i in range(0, 1000, 25)]
        yaxis2_3_range = [0, args.cut_threshold + 5]
        plot_height = 850 + 150 + 20
        legend_y = 1.090
    else:
        tick_vals = centers_rev + centers
        tickt_ext = integer_fractional_means_rev + integer_fractional_centers
        tickvals = [i for i in range(-1000, 1000, 25)]
        ticktext = [str(abs(i)) for i in range(-1000, 1000, 25)]
        yaxis2_3_range = [-(args.cut_threshold + 5), args.cut_threshold + 5]
        plot_height = 850 + 150 + 130 + 40 + 20
        legend_y = 1.05
    # #############################################################
    # #############################################################
    #fig.update_yaxes(range=[-1, args.cut_threshold])
    fig.update_layout(yaxis=dict(title="<b>Breakpoints</b>", range=[0, 75], showticklabels = False, showgrid=False, zeroline=False),
                      yaxis2=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis3=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis5=dict(title="<b>B-allele frequency</b>", range=[0, 0.6], showticklabels=True, showgrid=False, zeroline=False),
                      yaxis6=dict(title="<b>Genes</b>", range=[0, 1], showticklabels = False, showgrid=False, zeroline=True, zerolinewidth=2, zerolinecolor='black'),
                      yaxis4=dict(title="<b>Cytos</b>", range=[0, 1], showticklabels=False, showgrid=False, zeroline=True, zerolinewidth=2, zerolinecolor='black'),

                      xaxis=dict(showspikes=True, tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], showticklabels = False, showgrid=False, zeroline=False),
                      xaxis2=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], zeroline=True, zerolinewidth=1, zerolinecolor='black', showgrid=False,),
                      xaxis3=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size],showgrid=False,),
                      xaxis5=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist()) * args.bin_size], showgrid=False,),
                      xaxis6=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], showgrid=False, zeroline=True, zerolinewidth=1, zerolinecolor='black'),
                      xaxis4=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], showgrid=False, zeroline=True, zerolinewidth=1, zerolinecolor='black'))

    fig.update_layout(
        yaxis2=dict(
            title="<b>Coverage depth</b> (per bin)",
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext
        ),
        yaxis3=dict(
            title="<b>Copies</b> (integers)",
            zeroline=True, zerolinewidth=1, zerolinecolor='black',
            tickmode='array',
            tickvals=tick_vals,
            ticktext=tickt_ext  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
        ),
    )
    # fig.update_xaxes(
    #     tickmode='array',
    #     tickvals=label_pos_chrms,
    #     ticktext=label_chrms  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
    # )

    # fig.update_layout(
    #     hovermode="x unified",
    #     legend_traceorder="normal")

    fig.update_xaxes(
        #xaxis2=dict(
            tickangle=90,
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        #),
        #font=dict(size=18, color="black")
    )
    ax = 20
    ay = -30
    add_annotation(fig, 960000000, 75, ax, ay, "DEL", '#CF0759')
    add_annotation(fig, 960000000+250000000, 75, ax, ay, "INV", '#2830DE')
    add_annotation(fig, 960000000 + 250000000 + 250000000, 75, ax, ay, "INS", '#e0cf03')
    add_annotation(fig, 960000000 + 250000000 + 250000000 + 250000000, 75, ax, ay, "BND", '#737373')
    add_annotation(fig, 960000000 + 250000000 + 250000000 + 250000000 + 250000000, 75, ax, ay, "DUP", '#178117')

    if args.loh_enable:
        add_annotation(fig, 960000000 + 250000000 + 250000000, 0, ax, ay, "LOH", '#2980b9')
        #fig.add_annotation(text="Het SNPs ratio threshold: " + str(args.hets_ratio) , x = 960000000 + 250000000 + 250000000, y=args.cut_threshold, showarrow=False, row=2, col=1)

    # Update layout
    fig.update_layout(
        template="plotly_white",
        font_family="Times New Roman"
    )

    fig.update_layout(
        title=chrom,
    )
    # Legend
    fig.update_layout(legend=dict(
        orientation='h', xanchor="center", x=0.475, y=legend_y,  # orientation = 'v', xanchor = "center", x = 1.08, y= .5
    ))
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    #fig.update_xaxes(tick0=0.0, rangemode="nonnegative")
    fig.update_xaxes(spikedash='dashdot', spikemode='across', spikesnap='cursor', showspikes=True)


    fig.update_layout(legend={'itemsizing': 'constant'})

    fig.update_layout(font_family="Times New Roman")
    if args.without_phasing:
        genome_tile = args.genome_name
    else:
        genome_tile = args.genome_name + "<br>" + "<span style='color:blue'>Ploidy: </span>" + str(args.tumor_ploidy) + "     " + "<span style='color:blue'>Cellular tumor fraction: </span>" + str(args.tumor_purity) + "     " + "<span style='color:blue'>Confidence: </span>" + str(p_value)
    fig.update_layout(
        title={
            'text': genome_tile,
            'y': 0.98,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family="Courier New",
        font_color="dimgray",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="green",
    )
    fig.update_layout(
        width=1380,
        height=plot_height,
       )
    if args.without_phasing:
        path_set = args.genome_name + '_genome_copynumbers_breakpoints'
    else:
        if is_half:
            path_set = 'wgd/'+str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers_breakpoints'
        else:
            path_set = str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers_breakpoints'

    if args.pdf_enable:
        print_genome_pdf(fig, path_set, args.out_dir_plots)

    fig.write_html(args.out_dir_plots + '/' + path_set + ".html")

def copy_number_plots_genome_breakpoints(centers, integer_fractional_centers, df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, df_unphased, args, p_value, loh_regions, het_snps_bafs_means, is_half):
    # TODO genome-wide plots

    if args.breakpoints:
        coords_single, coords_single_chr, coords, coords_chr, _,_ = sv_vcf_bps_cn_check(args.breakpoints, args)

    lengths = []
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)

    #chroms = ['chr1', 'chr5','chr8']

    df_cnr_hp1_ = []
    df_segs_hp1_ = []
    df_cnr_hp2_ = []
    df_segs_hp2_ = []
    df_genes_1_ = []
    df_genes_2_ = []
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    #cancer_genes = os.path.join(fileDir, args.cancer_genes)
    df_genes = csv_df_chromosomes_sorter(args.cancer_genes, ['chr','start','end','gene'])
    genestart_1 = []
    genestart_2 = []
    last_len = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        df_segs_hp1_.append(df_segs_hp1[df_segs_hp1['chromosome'] == chrom])
        df_cnr_hp2_.append(df_cnr_hp2[df_cnr_hp2['chr'] == chrom])
        df_segs_hp2_.append(df_segs_hp2[df_segs_hp2['chromosome'] == chrom])

    df_cnr_hp1_1 = df_cnr_hp1
    df_1_ = []
    offset = 0
    for i, chrom in enumerate(chroms):
        df_genes_chrom_1 = df_cnr_hp1_1[df_cnr_hp1_1['chr'] == chrom]
        df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
        offset += regions[i]
        df_1_.append(df_genes_chrom_1)
    df_1 = pd.concat(df_1_)

    if args.breakpoints:
        offset = 0
        for i, chrom in enumerate(chroms):
            new_len = len(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
            if i > 0:
                for j, pair in enumerate(coords):
                    if chrom == pair[0]:
                        coords[j][1] = pair[1] + offset

            if i > 0:
                for k, pair in enumerate(coords_single):
                    if chrom == pair[0]:
                        coords_single[k][1] = pair[1] + offset

            last_len += new_len
            offset += regions[i]

        n = len(coords)
        coords = [x[1:] for x in coords]
        coords = [coords[i] + (coords[i + 1] if i + 1 < n else []) for i in range(0, n, 2)]
        coords = list(map(sorted, coords))
        last_len = 0
        coords_single = [j for i in [x[1:2] for x in coords_single] for j in i]

        arcs_data = get_all_breakpoints_data(coords, coords_chr, 75, args.breakpoints, args)
        arcs_data1 = get_all_breakpoints_data(coords, coords_chr, -75, args.breakpoints, args)

    df_cnr_hp1 = pd.concat(df_cnr_hp1_)
    df_segs_hp1 = pd.concat(df_segs_hp1_)

    df_cnr_hp2 = pd.concat(df_cnr_hp2_)
    df_segs_hp2 = pd.concat(df_segs_hp2_)

    df_genes_1 = df_genes[df_genes['start'] < df_genes['end']]
    df_genes_2 = df_genes[df_genes['start'] > df_genes['end']]

    last_len = 0
    offset = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        new_len = len(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        if not chrom.startswith('chr'):
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == 'chr' + chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == 'chr' + chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        else:
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        last_len += new_len
        offset += regions[i]

    df_genes_1 = pd.concat(df_genes_1_)
    df_genes_2 = pd.concat(df_genes_2_)

    genename_1 = df_genes_1['gene'].values.tolist()
    genename_2 = df_genes_2['gene'].values.tolist()

    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_cnr_hp1)*args.bin_size, args.bin_size, dtype=int)

    ###########################################################
    if args.without_phasing:
        row_heights = [220, 160, 150, 40]
    else:
        row_heights = [220, 320, 150, 40]

    #fig = make_subplots(rows=2, cols=1, shared_yaxes=False, shared_xaxes=True, vertical_spacing=0.02, horizontal_spacing=0.02)
    fig = make_subplots(rows=4, cols=1, shared_yaxes=False, shared_xaxes='columns',  vertical_spacing=0.01, row_heights=row_heights,
                        horizontal_spacing=0.02, specs=[[{"type":"xy"}], [{"secondary_y":True}], [{"type":"xy"}], [{"type":"xy"}]]) #https://community.plotly.com/t/can-subplot-support-multiple-y-axes/38891/20
    # #############################################################
    if args.breakpoints:
        for i in range(len(arcs_data)):
            fig.add_trace(go.Scatter(x=arcs_data[i][0],
                                     y=arcs_data[i][1],
                                     name=arcs_data[i][2],
                                     mode=arcs_data[i][3],
                                     line=arcs_data[i][4],
                                     yaxis="y",
                                     opacity=0.5,
                                     text = arcs_data[i][5],
                                     customdata = [arcs_data[i][6]],
                                     hovertemplate= '<br><b>Breakpoint info: </b>' + arcs_data[i][6],
                                     showlegend=arcs_data[i][7]), row=1, col=1,
                         )
    # #############################################################

    # #############################################################
    custom_text_data_hp1 = df_cnr_hp1.start.values.tolist()
    custom_text_data_hp2 = df_cnr_hp2.start.values.tolist()
    if args.without_phasing:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.coverage.values.tolist(), name='Unphased coverage',
                                   text=custom_text_data_hp1,
                                   yaxis="y2", opacity=0.7, color='olive', visibility=True, mul_cols=True)
    else:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.hp1.values.tolist(), name='HP-1 coverage',
                                   text=custom_text_data_hp1,
                                   yaxis=None, opacity=0.7, color='firebrick', visibility=True, mul_cols=True)
        add_scatter_trace_coverage(fig, df_1['start'], [-x for x in df_cnr_hp2.hp2.values.tolist()], name='HP-2 coverage',
                                   text=custom_text_data_hp2,
                                   yaxis=None, opacity=0.7, color='steelblue', visibility=True, mul_cols=True)

    add_scatter_trace_coverage(fig, df_1['start'], het_snps_bafs_means.vaf.values.tolist(), name='BAF',
                                   text=het_snps_bafs_means.pos.values.tolist(),
                                   yaxis="y4", opacity=0.7, color='olive', visibility='legendonly', mul_cols=True, row=3, baf=True)
        # if not args.unphased_reads_coverage_disable:
        #     add_scatter_trace_coverage(fig, indices, df_unphased.hp3.values.tolist(), name='Unphased', text=None,
        #                                yaxis=None, opacity=0.7, color='olive', visibility='legendonly', mul_cols=True)

    fig.add_trace(go.Scatter(x=df_genes_1['start'], y=[0.5]*len(df_genes_1),#[1, 2, 3, 4, 5] * (len(df_genes_1)//5),
                             mode='markers',
                             text = genename_1,
                             customdata = genestart_1,
                             #hovertext=df_genes['gene'],
                             hovertemplate=
                             '<br><b>Gene</b>: %{text}'+
                             '<br><b>Pos</b>: %{customdata}<br>',
                             marker=dict(
                                 symbol="y-left",
                                 color="#3A6B35",
                                 size=6,
                                 line=dict(width=1, color="#7F7F7F"),
                             ),
                             yaxis="y5",
                             name= 'GeneInfo',
                             showlegend=False,), row=4, col=1,)

    fig.add_hline(y=2,  line_width=1, line_dash="solid", line_color="black", row=3, col=1,)
    fig.add_hline(y=-(args.cut_threshold + 5),  line_width=1, line_dash="solid", line_color="black", row=2, col=1,)

    current = 0
    label_pos = []
    offset_chroms = 0
    offset_chroms_1 = 0
    chroms = get_contigs_list(args.contigs)
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])

    for index, chrom in enumerate(chroms):
        loh_starts = []
        if not loh_regions.empty:
            df_loh_region = loh_regions[loh_regions['chr'] == chrom]
            loh_starts = df_loh_region.start.values.tolist()
            loh_ends = df_loh_region.end.values.tolist()
        if not df_centm.empty:
            df_cent_region = df_centm[df_centm['chr'] == chrom]
            cent_starts = df_cent_region.start.values.tolist()
            cent_ends = df_cent_region.end.values.tolist()

        current += lengths[index]
        label_pos.append(round(offset_chroms+regions[index]//2))

        fig.add_vline(x=offset_chroms,  line_width=1, line_dash="solid", line_color="#D7DBDD")

        offset_chroms_1 += regions[index]

        if index % 2 == 0:
            fig.add_vrect(x0=offset_chroms, x1=offset_chroms_1, fillcolor="#E5E7E9", opacity=0.9, layer="below", line_width=0, )

        if len(loh_starts) and args.loh_enable:
            for i in range(len(df_loh_region.start.values.tolist())):
                fig.add_vrect(x0=offset_chroms+loh_starts[i], x1=offset_chroms+loh_ends[i], fillcolor="#2980b9", opacity=0.3, layer="below", line_width=0, row=2, col=1,)
        if len(cent_starts) and args.loh_enable:
            fig.add_vrect(x0=offset_chroms+cent_starts[0], x1=offset_chroms+cent_ends[0], fillcolor="#7e1f14", opacity=0.3, layer="above", line_width=0, row=2, col=1,)

        offset_chroms += regions[index]

    if not args.copynumbers_disable:
        offset_start = 0
        haplotype_1_start_values = []
        haplotype_1_end_values = []

        haplotype_2_start_values = []
        haplotype_2_end_values = []
        for index, chrom in enumerate(chroms):
            df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
            df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
            haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
            haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()
            haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
            haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()
            if chrom == args.contigs.split('-')[0]:
                haplotype_1_start_values.extend(haplotype_1_start_values_copyratios)
                haplotype_1_end_values.extend(haplotype_1_end_values_copyratios)
                haplotype_2_start_values.extend(haplotype_2_start_values_copyratios)
                haplotype_2_end_values.extend(haplotype_2_end_values_copyratios)
            else:
                offset_start += regions[index-1]#lengths[index-1] * args.bin_size
                haplotype_1_start_values.extend([x + offset_start for x in haplotype_1_start_values_copyratios])
                haplotype_1_end_values.extend([x + offset_start for x in haplotype_1_end_values_copyratios])
                haplotype_2_start_values.extend([x + offset_start for x in haplotype_2_start_values_copyratios])
                haplotype_2_end_values.extend([x + offset_start for x in haplotype_2_end_values_copyratios])

        haplotype_1_gaps_values = np.full(len(df_segs_hp1.state.values.tolist()), 'None')
        haplotype_1_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp1.state.values.tolist(), df_segs_hp1.state.values.tolist(), haplotype_1_gaps_values)))
        haplotype_1_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_1_start_values, haplotype_1_end_values, haplotype_1_gaps_values)))

        haplotype_2_gaps_values = np.full(len(df_segs_hp2.state.values.tolist()), 'None')
        haplotype_2_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp2.state.values.tolist(), df_segs_hp2.state.values.tolist(), haplotype_2_gaps_values)))
        haplotype_2_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_2_start_values, haplotype_2_end_values, haplotype_2_gaps_values)))


        if args.without_phasing:
            OFFSET = 0
            colors = ['darkolivegreen']
        else:
            OFFSET = args.cut_threshold/150
            colors = ['firebrick', 'steelblue']

        search_list = centers
        haplotype_1_copyratios_values = [-3300.0 if element not in search_list else element for element in haplotype_1_copyratios_values]
        haplotype_1_copyratios_values = [x if x == 'None' else x for x in haplotype_1_copyratios_values]

        haplotype_2_copyratios_values = [-3300.0 if element not in search_list else element for element in haplotype_2_copyratios_values]
        haplotype_2_copyratios_values = [x if x == 'None' else -x for x in haplotype_2_copyratios_values]

        #haplotype_1_copyratios_values = [x if x == 'None' else x  for x in haplotype_1_copyratios_values]
        #haplotype_2_copyratios_values = [x if x == 'None' else -(x) for x in haplotype_2_copyratios_values]
        name = "Copynumbers"

        add_scatter_trace_copyratios(args, fig, colors, name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values, df_segs_hp1, df_segs_hp2, mul_cols=True, centers=centers)

    centers_rev = [-x for x in centers]
    centers_rev.reverse()
    integer_fractional_means_rev = [x for x in integer_fractional_centers]
    integer_fractional_means_rev.reverse()
    if args.without_phasing:
        tick_vals = centers
        tickt_ext = integer_fractional_centers
        tickvals = [i for i in range(0, 1000, 25)]
        ticktext = [str(abs(i)) for i in range(0, 1000, 25)]
        yaxis2_3_range = [0, args.cut_threshold + 5]
        plot_height = 850 + 150 + 20 + 20
        legend_y = 1.085
    else:
        tick_vals = centers_rev + centers
        tickt_ext = integer_fractional_means_rev + integer_fractional_centers
        tickvals = [i for i in range(-1000, 1000, 25)]
        ticktext = [str(abs(i)) for i in range(-1000, 1000, 25)]
        yaxis2_3_range = [-(args.cut_threshold + 5), args.cut_threshold + 5]
        plot_height = 850 + 150 + 130+ 10
        legend_y = 1.06
    # #############################################################
    # #############################################################
    #fig.update_yaxes(range=[-1, args.cut_threshold])
    fig.update_layout(yaxis=dict(title="<b>Breakpoints</b>", range=[0, 75], showticklabels = False, showgrid=False, zeroline=False),
                      yaxis2=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis3=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis4=dict(title="<b>B-allele frequency</b>", range=[0, 0.6], showticklabels=True, showgrid=False, zeroline=False),
                      yaxis5=dict(title="<b>Genes</b>", range=[0, 1], showticklabels = False, showgrid=False, zeroline=True, zerolinewidth=2, zerolinecolor='black'),

                      xaxis=dict(showspikes=True, tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], showticklabels = False, showgrid=False, zeroline=False),
                      xaxis2=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], zeroline=True, zerolinewidth=1, zerolinecolor='black', showgrid=False,),
                      xaxis3=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size],showgrid=False,),
                      xaxis4=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist()) * args.bin_size], showgrid=False,),
                      xaxis5=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], showgrid=False, zeroline=True, zerolinewidth=1, zerolinecolor='black'))

    fig.update_layout(
        yaxis2=dict(
            title="<b>Coverage depth</b> (per bin)",
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext
        ),
        yaxis3=dict(
            title="<b>Copies</b> (integers/fractions)",
            zeroline=True, zerolinewidth=1, zerolinecolor='black',
            tickmode='array',
            tickvals=tick_vals,
            ticktext=tickt_ext  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
        ),
    )
    # fig.update_xaxes(
    #     tickmode='array',
    #     tickvals=label_pos_chrms,
    #     ticktext=label_chrms  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
    # )

    # fig.update_layout(
    #     hovermode="x unified",
    #     legend_traceorder="normal")

    fig.update_xaxes(
        #xaxis2=dict(
            tickangle=90,
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        #),
        #font=dict(size=18, color="black")
    )

    fig.add_annotation(text="DEL",xref="paper",yref="paper",x=0.25,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#CF0759',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)
    fig.add_annotation(text="INV",xref="paper",yref="paper",x=0.30,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#2830DE',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)
    fig.add_annotation(text="INS",xref="paper",yref="paper",x=0.36,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#e0cf03',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)
    fig.add_annotation(text="BND",xref="paper",yref="paper",x=0.42,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#737373',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)
    fig.add_annotation(text="DUP",xref="paper",yref="paper",x=0.47,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#178117',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)

    fig.add_annotation(text="LOH Regions",xref="paper",yref="paper",x=0.60,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#2980b9',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)
    fig.add_annotation(text="Centromeres",xref="paper",yref="paper",x=0.75,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#7e1f14',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)

    # Update layout
    fig.update_layout(
        template="plotly_white",
        font_family="Times New Roman"
    )

    fig.update_layout(
        title=chrom,
    )
    # Legend
    fig.update_layout(legend=dict(
        orientation='h', xanchor="center", x=0.475, y=legend_y,  # orientation = 'v', xanchor = "center", x = 1.08, y= .5
    ))
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    #fig.update_xaxes(tick0=0.0, rangemode="nonnegative")


    fig.update_layout(legend={'itemsizing': 'constant'})

    fig.update_layout(font_family="Times New Roman")
    if args.without_phasing:
        genome_tile = args.genome_name
    else:
        genome_tile = args.genome_name + "<br>" + "<span style='color:blue'>Ploidy: </span>" + str(args.tumor_ploidy) + "     " + "<span style='color:blue'>Cellular tumor fraction: </span>" + str(args.tumor_purity) + "     " + "<span style='color:blue'>Confidence: </span>" + str(p_value)
    fig.update_layout(
        title={
            'text': genome_tile,
            'y': 0.98,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family="Courier New",
        font_color="dimgray",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="green",
    )
    fig.update_layout(
        width=1380,
        height=plot_height,
       )
    if args.without_phasing:
        path_set = args.genome_name + '_genome_copynumbers_breakpoints'
    else:
        if is_half:
            path_set = 'wgd/'+str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers_breakpoints'
        else:
            path_set = str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers_breakpoints'

    if args.pdf_enable:
        print_genome_pdf(fig, path_set, args.out_dir_plots)

    fig.write_html(args.out_dir_plots + '/' + path_set + ".html")

def copy_number_plots_genome(centers, integer_fractional_centers, df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, df_unphased, args, p_value, loh_regions, het_snps_bafs_means, is_half):
    # TODO genome-wide plots
    lengths = []
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)

    df_cnr_hp1_ = []
    df_segs_hp1_ = []
    df_cnr_hp2_ = []
    df_segs_hp2_ = []
    df_genes_1_ = []
    df_genes_2_ = []
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    #cancer_genes = os.path.join(fileDir, args.cancer_genes)
    df_genes = csv_df_chromosomes_sorter(args.cancer_genes, ['chr','start','end','gene'])
    genestart_1 = []
    genestart_2 = []
    last_len = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        df_segs_hp1_.append(df_segs_hp1[df_segs_hp1['chromosome'] == chrom])
        df_cnr_hp2_.append(df_cnr_hp2[df_cnr_hp2['chr'] == chrom])
        df_segs_hp2_.append(df_segs_hp2[df_segs_hp2['chromosome'] == chrom])

    df_cnr_hp1_1 = df_cnr_hp1
    df_1_ = []
    offset = 0
    for i, chrom in enumerate(chroms):
        df_genes_chrom_1 = df_cnr_hp1_1[df_cnr_hp1_1['chr'] == chrom]
        df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
        offset += regions[i]
        df_1_.append(df_genes_chrom_1)
    df_1 = pd.concat(df_1_)

    df_cnr_hp1 = pd.concat(df_cnr_hp1_)
    df_segs_hp1 = pd.concat(df_segs_hp1_)

    df_cnr_hp2 = pd.concat(df_cnr_hp2_)
    df_segs_hp2 = pd.concat(df_segs_hp2_)

    df_genes_1 = df_genes[df_genes['start'] < df_genes['end']]
    df_genes_2 = df_genes[df_genes['start'] > df_genes['end']]

    last_len = 0
    offset = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        new_len = len(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        if not chrom.startswith('chr'):
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == 'chr' + chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == 'chr' + chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        else:
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        last_len += new_len
        offset += regions[i]

    df_genes_1 = pd.concat(df_genes_1_)
    df_genes_2 = pd.concat(df_genes_2_)

    genename_1 = df_genes_1['gene'].values.tolist()
    genename_2 = df_genes_2['gene'].values.tolist()

    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_cnr_hp1)*args.bin_size, args.bin_size, dtype=int)

    ###########################################################
    if args.without_phasing:
        row_heights = [160, 150, 40]
    else:
        row_heights = [320, 150, 40]

    #fig = make_subplots(rows=2, cols=1, shared_yaxes=False, shared_xaxes=True, vertical_spacing=0.02, horizontal_spacing=0.02)
    fig = make_subplots(rows=3, cols=1, shared_yaxes=False, shared_xaxes='columns',  vertical_spacing=0.01, row_heights=row_heights,
                        horizontal_spacing=0.02, specs=[[{"secondary_y":True}], [{"type":"xy"}], [{"type":"xy"}]]) #https://community.plotly.com/t/can-subplot-support-multiple-y-axes/38891/20
    # #############################################################

    # #############################################################
    custom_text_data_hp1 = df_cnr_hp1.start.values.tolist()
    custom_text_data_hp2 = df_cnr_hp2.start.values.tolist()
    if args.without_phasing:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.coverage.values.tolist(), name='Unphased coverage',
                                   text=custom_text_data_hp1,
                                   yaxis="y", opacity=0.7, color='olive', visibility=True, mul_cols=True, row=1)
    else:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.hp1.values.tolist(), name='HP-1 coverage',
                                   text=custom_text_data_hp1,
                                   yaxis="y", opacity=0.7, color='firebrick', visibility=True, mul_cols=True, row=1)
        add_scatter_trace_coverage(fig, df_1['start'], [-x for x in df_cnr_hp2.hp2.values.tolist()], name='HP-2 coverage',
                                   text=custom_text_data_hp2,
                                   yaxis="y", opacity=0.7, color='steelblue', visibility=True, mul_cols=True, row=1)

    add_scatter_trace_coverage(fig, df_1['start'], het_snps_bafs_means.vaf.values.tolist(), name='BAF',
                                   text=het_snps_bafs_means.pos.values.tolist(),
                                   yaxis="y3", opacity=0.7, color='olive', visibility='legendonly', mul_cols=True, row=2, baf=True)
        # if not args.unphased_reads_coverage_disable:
        #     add_scatter_trace_coverage(fig, indices, df_unphased.hp3.values.tolist(), name='Unphased', text=None,
        #                                yaxis=None, opacity=0.7, color='olive', visibility='legendonly', mul_cols=True)

    fig.add_trace(go.Scatter(x=df_genes_1['start'], y=[0.5]*len(df_genes_1),#[1, 2, 3, 4, 5] * (len(df_genes_1)//5),
                             mode='markers',
                             text = genename_1,
                             customdata = genestart_1,
                             #hovertext=df_genes['gene'],
                             hovertemplate=
                             '<br><b>Gene</b>: %{text}'+
                             '<br><b>Pos</b>: %{customdata}<br>',
                             marker=dict(
                                 symbol="y-left",
                                 color="#3A6B35",
                                 size=6,
                                 line=dict(width=1, color="#7F7F7F"),
                             ),
                             yaxis="y4",
                             name= 'GeneInfo',
                             showlegend=False,), row=3, col=1,)

    fig.add_hline(y=2,  line_width=1, line_dash="solid", line_color="black", row=2, col=1,)
    fig.add_hline(y=-(args.cut_threshold + 5),  line_width=1, line_dash="solid", line_color="black", row=1, col=1,)

    current = 0
    label_pos = []
    offset_chroms = 0
    offset_chroms_1 = 0
    chroms = get_contigs_list(args.contigs)
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])

    for index, chrom in enumerate(chroms):
        loh_starts = []
        if not loh_regions.empty:
            df_loh_region = loh_regions[loh_regions['chr'] == chrom]
            loh_starts = df_loh_region.start.values.tolist()
            loh_ends = df_loh_region.end.values.tolist()
        if not df_centm.empty:
            df_cent_region = df_centm[df_centm['chr'] == chrom]
            cent_starts = df_cent_region.start.values.tolist()
            cent_ends = df_cent_region.end.values.tolist()

        current += lengths[index]
        label_pos.append(round(offset_chroms+regions[index]//2))

        fig.add_vline(x=offset_chroms,  line_width=1, line_dash="solid", line_color="#D7DBDD")

        offset_chroms_1 += regions[index]

        if index % 2 == 0:
            fig.add_vrect(x0=offset_chroms, x1=offset_chroms_1, fillcolor="#E5E7E9", opacity=0.9, layer="below", line_width=0, )

        if len(loh_starts) and args.loh_enable:
            for i in range(len(df_loh_region.start.values.tolist())):
                fig.add_vrect(x0=offset_chroms+loh_starts[i], x1=offset_chroms+loh_ends[i], fillcolor="#2980b9", opacity=0.3, layer="below", line_width=0, row=1, col=1,)
        if len(cent_starts) and args.loh_enable:
            fig.add_vrect(x0=offset_chroms+cent_starts[0], x1=offset_chroms+cent_ends[0], fillcolor="#7e1f14", opacity=0.3, layer="above", line_width=0, row=1, col=1,)

        offset_chroms += regions[index]

    if not args.copynumbers_disable:
        offset_start = 0
        haplotype_1_start_values = []
        haplotype_1_end_values = []

        haplotype_2_start_values = []
        haplotype_2_end_values = []
        for index, chrom in enumerate(chroms):
            df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
            df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
            haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
            haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()
            haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
            haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()
            if chrom == args.contigs.split('-')[0]:
                haplotype_1_start_values.extend(haplotype_1_start_values_copyratios)
                haplotype_1_end_values.extend(haplotype_1_end_values_copyratios)
                haplotype_2_start_values.extend(haplotype_2_start_values_copyratios)
                haplotype_2_end_values.extend(haplotype_2_end_values_copyratios)
            else:
                offset_start += regions[index-1]#lengths[index-1] * args.bin_size
                haplotype_1_start_values.extend([x + offset_start for x in haplotype_1_start_values_copyratios])
                haplotype_1_end_values.extend([x + offset_start for x in haplotype_1_end_values_copyratios])
                haplotype_2_start_values.extend([x + offset_start for x in haplotype_2_start_values_copyratios])
                haplotype_2_end_values.extend([x + offset_start for x in haplotype_2_end_values_copyratios])

        haplotype_1_gaps_values = np.full(len(df_segs_hp1.state.values.tolist()), 'None')
        haplotype_1_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp1.state.values.tolist(), df_segs_hp1.state.values.tolist(), haplotype_1_gaps_values)))
        haplotype_1_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_1_start_values, haplotype_1_end_values, haplotype_1_gaps_values)))

        haplotype_2_gaps_values = np.full(len(df_segs_hp2.state.values.tolist()), 'None')
        haplotype_2_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp2.state.values.tolist(), df_segs_hp2.state.values.tolist(), haplotype_2_gaps_values)))
        haplotype_2_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_2_start_values, haplotype_2_end_values, haplotype_2_gaps_values)))


        if args.without_phasing:
            OFFSET = 0
            colors = ['darkolivegreen']
        else:
            OFFSET = args.cut_threshold/150
            colors = ['firebrick', 'steelblue']
        haplotype_1_copyratios_values = [x if x == 'None' else x  for x in haplotype_1_copyratios_values]
        haplotype_2_copyratios_values = [x if x == 'None' else -(x) for x in haplotype_2_copyratios_values]
        name = "Copynumbers"

        add_scatter_trace_copyratios(args, fig, colors, name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values, df_segs_hp1, df_segs_hp2, mul_cols=True, row=1)

    centers_rev = [-x for x in centers]
    centers_rev.reverse()
    integer_fractional_means_rev = [x for x in integer_fractional_centers]
    integer_fractional_means_rev.reverse()
    if args.without_phasing:
        tick_vals = centers
        tickt_ext = integer_fractional_centers
        tickvals = [i for i in range(0, 1000, 25)]
        ticktext = [str(abs(i)) for i in range(0, 1000, 25)]
        yaxis2_3_range = [0, args.cut_threshold + 5]
        plot_height = 630 + 150 + 20
        legend_y = 1.085
    else:
        tick_vals = centers_rev + centers
        tickt_ext = integer_fractional_means_rev + integer_fractional_centers
        tickvals = [i for i in range(-1000, 1000, 25)]
        ticktext = [str(abs(i)) for i in range(-1000, 1000, 25)]
        yaxis2_3_range = [-(args.cut_threshold + 5), args.cut_threshold + 5]
        plot_height = 630 + 150 + 130
        legend_y = 1.06
    # #############################################################
    # #############################################################
    #fig.update_yaxes(range=[-1, args.cut_threshold])
    fig.update_layout(
                      yaxis=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis2=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis3=dict(title="<b>B-allele frequency</b>", range=[0, 0.6], showticklabels=True, showgrid=False, zeroline=False),
                      yaxis4=dict(title="<b>Genes</b>", range=[0, 1], showticklabels = False, showgrid=False, zeroline=True, zerolinewidth=2, zerolinecolor='black'),

                      xaxis=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], zeroline=True, zerolinewidth=1, zerolinecolor='black', showgrid=False,),
                      xaxis2=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size],showgrid=False,),
                      xaxis3=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist()) * args.bin_size], showgrid=False,),
                      xaxis4=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], showgrid=False, zeroline=True, zerolinewidth=1, zerolinecolor='black'))

    fig.update_layout(
        yaxis=dict(
            title="<b>Coverage depth</b> (per bin)",
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext
        ),
        yaxis2=dict(
            title="<b>Copies</b> (integers)",
            zeroline=True, zerolinewidth=1, zerolinecolor='black',
            tickmode='array',
            tickvals=tick_vals,
            ticktext=tickt_ext  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
        ),
    )
    # fig.update_xaxes(
    #     tickmode='array',
    #     tickvals=label_pos_chrms,
    #     ticktext=label_chrms  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
    # )

    # fig.update_layout(
    #     hovermode="x unified",
    #     legend_traceorder="normal")

    fig.update_xaxes(
        #xaxis2=dict(
            tickangle=90,
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        #),
        #font=dict(size=18, color="black")
    )
    ax = 20
    ay = -30

    if args.loh_enable:
        add_annotation(fig, 960000000 + 250000000 + 250000000, 0, ax, ay, "LOH regions", '#2980b9')
        add_annotation(fig, 960000000 + 250000000 + 250000000 + 250000000 + 250000000, 0, ax, ay, "Centromeres", '#7e1f14')
        #fig.add_annotation(text="Het SNPs ratio threshold: " + str(args.hets_ratio) , x = 960000000 + 250000000 + 250000000, y=args.cut_threshold, showarrow=False, row=2, col=1)

    # Update layout
    fig.update_layout(
        template="plotly_white",
        font_family="Times New Roman"
    )

    fig.update_layout(
        title=chrom,
    )
    # Legend
    fig.update_layout(legend=dict(
        orientation='h', xanchor="center", x=0.475, y=legend_y,  # orientation = 'v', xanchor = "center", x = 1.08, y= .5
    ))
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    #fig.update_xaxes(tick0=0.0, rangemode="nonnegative")


    fig.update_layout(legend={'itemsizing': 'constant'})

    fig.update_layout(font_family="Times New Roman")
    if args.without_phasing:
        genome_tile = args.genome_name
    else:
        genome_tile = args.genome_name + "<br>" + "<span style='color:blue'>Ploidy: </span>" + str(args.tumor_ploidy) + "     " + "<span style='color:blue'>Cellular tumor fraction: </span>" + str(args.tumor_purity) + "     " + "<span style='color:blue'>Confidence: </span>" + str(p_value)
    fig.update_layout(
        title={
            'text': genome_tile,
            'y': 0.98,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family="Courier New",
        font_color="dimgray",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="green",
    )
    fig.update_layout(
        width=1380,
        height=plot_height,
       )
    if args.without_phasing:
        path_set = args.genome_name + '_genome_copynumbers'
    else:
        if is_half:
            path_set = 'wgd/'+str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers'
        else:
            path_set = str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers'

    if args.pdf_enable:
        print_genome_pdf(fig, path_set, args.out_dir_plots)

    fig.write_html(args.out_dir_plots + '/' + path_set + ".html")

def copy_number_plots_genome_breakpoints_subclonal_cytos(centers, integer_fractional_centers, df_cnr_hp1, df_segs_hp1_, df_cnr_hp2, df_segs_hp2_, df_unphased, args, p_value, loh_regions, het_snps_bafs_means, is_half):
    # TODO genome-wide plots
    df_segs_hp1 = df_segs_hp1_.copy()
    df_segs_hp2 = df_segs_hp2_.copy()

    coords_single, coords_single_chr, coords, coords_chr, _,_ = sv_vcf_bps_cn_check(args.breakpoints, args)

    lengths = []
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)

    #chroms = ['chr1', 'chr5','chr8']

    df_cnr_hp1_ = []
    df_segs_hp1_ = []
    df_cnr_hp2_ = []
    df_segs_hp2_ = []
    df_genes_1_ = []
    df_genes_2_ = []
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    #cancer_genes = os.path.join(fileDir, args.cancer_genes)
    df_genes = csv_df_chromosomes_sorter(args.cancer_genes, ['chr','start','end','gene'])
    genestart_1 = []
    genestart_2 = []
    last_len = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        df_segs_hp1_.append(df_segs_hp1[df_segs_hp1['chromosome'] == chrom])
        df_cnr_hp2_.append(df_cnr_hp2[df_cnr_hp2['chr'] == chrom])
        df_segs_hp2_.append(df_segs_hp2[df_segs_hp2['chromosome'] == chrom])

    #coords_chr =  [['chr1',223550001, 'HP-1'],['chr6',60250001, 'HP-1'], ['chr6',27900001, 'HP-1'],['chr15',22200001, 'HP-1'], ['chr1', 144850001, 'HP-2'], ['chr9', 68300001, 'HP-2']]
    #coords = [['chr1',223550001],['chr6',60250001], ['chr6',27900001],['chr15',22200001], ['chr1', 144850001], ['chr9', 68300001]]
    #coords_single = [['chr1',22355007, 22355008], ['chr2', 123444555, 123444556], ['chr19', 1355007, 1355008]]
    #coords_single_chr = [['chr1', 22355007,22355008, 'DEL'], ['chr2', 123444555,123444556, 'DUP'], ['chr19', 1355007,1355008, 'DEL']]

    offset = 0
    for i, chrom in enumerate(chroms):
        new_len = len(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        # [[9:64526327,1:9769900],
        # [6:91468745,1:9769559]]
        if i > 0:
            for j, pair in enumerate(coords):
                if chrom == pair[0]:
                    coords[j][1] = pair[1] + offset#(last_len * (args.bin_size))

        if i > 0:
            for k, pair in enumerate(coords_single):
                if chrom == pair[0]:
                    coords_single[k][1] = pair[1] + offset  # (last_len * (args.bin_size))

        last_len += new_len
        offset += regions[i]

    df_cnr_hp1_1 = df_cnr_hp1
    df_1_ = []
    offset = 0
    for i, chrom in enumerate(chroms):
        df_genes_chrom_1 = df_cnr_hp1_1[df_cnr_hp1_1['chr'] == chrom]
        df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
        offset += regions[i]
        df_1_.append(df_genes_chrom_1)
    df_1 = pd.concat(df_1_)

    n = len(coords)
    coords = [x[1:] for x in coords]
    coords = [coords[i] + (coords[i + 1] if i + 1 < n else []) for i in range(0, n, 2)]
    coords = list(map(sorted, coords))

    last_len = 0
    coords_single = [j for i in [x[1:2] for x in coords_single] for j in i]

    df_cnr_hp1 = pd.concat(df_cnr_hp1_)
    df_segs_hp1 = pd.concat(df_segs_hp1_)

    df_cnr_hp2 = pd.concat(df_cnr_hp2_)
    df_segs_hp2 = pd.concat(df_segs_hp2_)

    df_genes_1 = df_genes[df_genes['start'] < df_genes['end']]
    df_genes_2 = df_genes[df_genes['start'] > df_genes['end']]

    last_len = 0
    offset = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        new_len = len(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        if not chrom.startswith('chr'):
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == 'chr' + chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == 'chr' + chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        else:
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        last_len += new_len
        offset += regions[i]

    df_genes_1 = pd.concat(df_genes_1_)
    df_genes_2 = pd.concat(df_genes_2_)

    genename_1 = df_genes_1['gene'].values.tolist()
    genename_2 = df_genes_2['gene'].values.tolist()

    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_cnr_hp1)*args.bin_size, args.bin_size, dtype=int)

    ###########################################################

    arcs_data = get_all_breakpoints_data(coords, coords_chr, 75, args.breakpoints, args)
    arcs_data1 = get_all_breakpoints_data(coords, coords_chr, -75, args.breakpoints, args)
    # fig = go.Figure(data=arcs_data)

    #fig = go.Figure()
    #fig = go.Figure().set_subplots(rows=2, cols=1)
    if args.without_phasing:
        row_heights = [220, 160, 150, 40]
    else:
        row_heights = [220, 320, 20, 150, 40]

    #fig = make_subplots(rows=2, cols=1, shared_yaxes=False, shared_xaxes=True, vertical_spacing=0.02, horizontal_spacing=0.02)
    fig = make_subplots(rows=5, cols=1, shared_yaxes=False, shared_xaxes='columns',  vertical_spacing=0.01, row_heights=row_heights,
                        horizontal_spacing=0.02, specs=[[{"type":"xy"}], [{"secondary_y":True}], [{"type":"xy"}], [{"type":"xy"}], [{"type":"xy"}]]) #https://community.plotly.com/t/can-subplot-support-multiple-y-axes/38891/20
    # #############################################################
    for i in range(len(arcs_data)):
        fig.add_trace(go.Scatter(x=arcs_data[i][0],
                                 y=arcs_data[i][1],
                                 name=arcs_data[i][2],
                                 mode=arcs_data[i][3],
                                 line=arcs_data[i][4],
                                 yaxis="y",
                                 opacity=0.5,
                                 text = arcs_data[i][5],
                                 hovertemplate= '<br><b>Breakpoint info: </b>' + arcs_data[i][6],
                                 showlegend=arcs_data[i][7]), row=1, col=1,
                     )
    # #############################################################

    # #############################################################
    custom_text_data_hp1 = df_cnr_hp1.start.values.tolist()
    custom_text_data_hp2 = df_cnr_hp2.start.values.tolist()
    if args.without_phasing:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.coverage.values.tolist(), name='Unphased coverage',
                                   text=custom_text_data_hp1,
                                   yaxis="y2", opacity=0.7, color='olive', visibility=True, mul_cols=True)
    else:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.hp1.values.tolist(), name='HP-1 coverage',
                                   text=custom_text_data_hp1,
                                   yaxis=None, opacity=0.7, color='firebrick', visibility=True, mul_cols=True)
        add_scatter_trace_coverage(fig, df_1['start'], [-x for x in df_cnr_hp2.hp2.values.tolist()], name='HP-2 coverage',
                                   text=custom_text_data_hp2,
                                   yaxis=None, opacity=0.7, color='steelblue', visibility=True, mul_cols=True)
        # if not args.unphased_reads_coverage_disable:
        #     add_scatter_trace_coverage(fig, indices, df_unphased.hp3.values.tolist(), name='Unphased', text=None,
        #                                yaxis=None, opacity=0.7, color='olive', visibility='legendonly', mul_cols=True)

    add_scatter_trace_coverage(fig, df_1['start'], het_snps_bafs_means.vaf.values.tolist(), name='BAF',
                                   text=het_snps_bafs_means.pos.values.tolist(),
                                   yaxis="y5", opacity=0.7, color='olive', visibility='legendonly', mul_cols=True, row=4, baf=True)

    fig.add_trace(go.Scatter(x=df_genes_1['start'], y=[0.5] * len(df_genes_1), #[1, 2, 3, 4, 5] * (len(df_genes_1)//5),
                             mode='markers',
                             text = genename_1,
                             customdata = genestart_1,
                             #hovertext=df_genes['gene'],
                             hovertemplate=
                             '<br><b>Gene</b>: %{text}'+
                             '<br><b>Pos</b>: %{customdata}<br>',
                             marker=dict(
                                 symbol="y-left",
                                 color="#3A6B35",
                                 size=6,
                                 line=dict(width=1, color="#7F7F7F"),
                             ),
                             yaxis="y6",
                             name= 'GeneInfo',
                             showlegend=False,), row=5, col=1,)

    # fig.add_trace(go.Scatter(x=df_genes_2['start'], y=[1, 2, 3, 4, 5] * (len(df_genes_2)//5),
    #                          mode='markers',
    #                          text = genename_2,
    #                          customdata = genestart_2,
    #                          #hovertext=df_genes['gene'],
    #                          hovertemplate=
    #                          '<br><b>Gene</b>: %{text}'+
    #                          '<br><b>Pos</b>: %{customdata}<br>',
    #                          marker=dict(
    #                              symbol="y-right",
    #                              color="#E3B448",
    #                              size=6,
    #                              line=dict(width=1, color="#7F7F7F"),
    #                          ),
    #                          yaxis="y4",
    #                          name= 'GeneInfo',
    #                          showlegend=False,), row=3, col=1,)

    #cytobands
    fileDir = os.path.dirname(__file__)
    if '13' in args.reference_name:
        cytos_file = os.path.join(fileDir, "annotations/chm13.cytobands.bed")
    else:
        cytos_file = os.path.join(fileDir, "annotations/grch38.cytobands.bed")
    df = pd.read_csv(cytos_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'gieStain'])
    #df = df_chromosomes_sorter(df, names=['chr', 'start', 'end', 'name', 'gieStain'])
    plot_cytobands_genome_wide(fig, df)

    current = 0
    label_pos = []
    label_pos_chrms = []
    label_chrms = []



    #coords_lines = [j for i in coords for j in i] + coords_single
    #for i, co in enumerate(coords_single):
        # if coords_single_chr[i][3] == 'DEL':
        #     fig.add_vline(x=co,  line_width=1, line_dash="solid", opacity=0.3, line_color="#CF0759", row=1, col=1,)
        # elif coords_single_chr[i][3] == 'DUP' or coords_single_chr[i][3] == 'INS':
        #     fig.add_vline(x=co,  line_width=1, line_dash="solid", opacity=0.3, line_color="#178117", row=1, col=1,)
        # else:

        #fig.add_vline(x=co,  line_width=1, line_dash="solid", opacity=0.3, line_color="#2B40A0", row=1, col=1,)

    fig.add_hline(y=2,  line_width=1, line_dash="solid", line_color="black", row=3, col=1,)
    fig.add_hline(y=-(args.cut_threshold + 5),  line_width=1, line_dash="solid", line_color="black", row=2, col=1,)

    # for index, chrom in enumerate(chroms):
    #     current += lengths[index]
    #     label_pos.append(round(current - (lengths[index] / 2))*args.bin_size) #y0=-10, y1=args.cut_threshold
    #
    #     label_pos_chrms.append(round(current - (lengths[index]))*args.bin_size)
    #     label_pos_chrms.append(round(current - (lengths[index]*.75)) * args.bin_size)
    #     label_pos_chrms.append(round(current - (lengths[index] / 2))*args.bin_size)
    #     label_pos_chrms.append(round(current - (lengths[index]*.25)) * args.bin_size)
    #
    #     label_chrms.append('0')
    #     label_chrms.append(str((round((lengths[index] * .25) * args.bin_size)) // 1000000) + 'M')
    #     label_chrms.append(str(((lengths[index]//2)*args.bin_size)// 1000000)+'M')
    #     label_chrms.append(str((round((lengths[index] * .75) * args.bin_size)) // 1000000) + 'M')
    #
    #     fig.add_vline(x=current*args.bin_size,  line_width=2, line_dash="solid", line_color="black", row=2, col=1,)

    current = 0
    label_pos = []
    offset_chroms = 0
    offset_chroms_1 = 0
    chroms = get_contigs_list(args.contigs)
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])

    for index, chrom in enumerate(chroms):
        loh_starts = []
        if not loh_regions.empty:
            df_loh_region = loh_regions[loh_regions['chr'] == chrom]
            loh_starts = df_loh_region.start.values.tolist()
            loh_ends = df_loh_region.end.values.tolist()
        if not df_centm.empty:
            df_cent_region = df_centm[df_centm['chr'] == chrom]
            cent_starts = df_cent_region.start.values.tolist()
            cent_ends = df_cent_region.end.values.tolist()

        current += lengths[index]
        label_pos.append(round(offset_chroms+regions[index]//2))

        fig.add_vline(x=offset_chroms,  line_width=1, line_dash="solid", line_color="#D7DBDD")

        #label_chrms.append('0')
        #label_chrms.append(str(((lengths[index]//2) + regions[index])// 1000000)+'M')

        offset_chroms_1 += regions[index]

        if index % 2 == 0:
            fig.add_vrect(x0=offset_chroms, x1=offset_chroms_1, fillcolor="#E5E7E9", opacity=0.9, layer="below", line_width=0, )

        if len(loh_starts) and args.loh_enable:
            for i in range(len(df_loh_region.start.values.tolist())):
                fig.add_vrect(x0=offset_chroms+loh_starts[i], x1=offset_chroms+loh_ends[i], fillcolor="#2980b9", opacity=0.3, layer="below", line_width=0, row=2, col=1,)
        if len(cent_starts) and args.loh_enable:
            fig.add_vrect(x0=offset_chroms+cent_starts[0], x1=offset_chroms+cent_ends[0], fillcolor="#7e1f14", opacity=0.3, layer="above", line_width=0, row=2, col=1,)

        offset_chroms += regions[index]

    if not args.copynumbers_disable:
        offset_start = 0
        haplotype_1_start_values = []
        haplotype_1_end_values = []

        haplotype_2_start_values = []
        haplotype_2_end_values = []
        for index, chrom in enumerate(chroms):
            df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
            df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
            haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
            haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()
            haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
            haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()
            if chrom == args.contigs.split('-')[0]:
                haplotype_1_start_values.extend(haplotype_1_start_values_copyratios)
                haplotype_1_end_values.extend(haplotype_1_end_values_copyratios)
                haplotype_2_start_values.extend(haplotype_2_start_values_copyratios)
                haplotype_2_end_values.extend(haplotype_2_end_values_copyratios)
            else:
                offset_start += regions[index-1]#lengths[index-1] * args.bin_size
                haplotype_1_start_values.extend([x + offset_start for x in haplotype_1_start_values_copyratios])
                haplotype_1_end_values.extend([x + offset_start for x in haplotype_1_end_values_copyratios])
                haplotype_2_start_values.extend([x + offset_start for x in haplotype_2_start_values_copyratios])
                haplotype_2_end_values.extend([x + offset_start for x in haplotype_2_end_values_copyratios])

        haplotype_1_gaps_values = np.full(len(df_segs_hp1.state.values.tolist()), 'None')
        haplotype_1_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp1.state.values.tolist(), df_segs_hp1.state.values.tolist(), haplotype_1_gaps_values)))
        haplotype_1_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_1_start_values, haplotype_1_end_values, haplotype_1_gaps_values)))

        haplotype_2_gaps_values = np.full(len(df_segs_hp2.state.values.tolist()), 'None')
        haplotype_2_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp2.state.values.tolist(), df_segs_hp2.state.values.tolist(), haplotype_2_gaps_values)))
        haplotype_2_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_2_start_values, haplotype_2_end_values, haplotype_2_gaps_values)))

        if args.copynumbers_subclonal_enable:

            search_list = centers #[centers[i] for i in [integer_fractional_centers.index(x) for x in integer_fractional_centers if not float(x).is_integer()]] #[27.025, 42.025]
            search_list_none = search_list + ['None']

            haplotype_1_copyratios_values_normal = [-3300.0 if element not in search_list else element for element in haplotype_1_copyratios_values]
            haplotype_1_copyratios_values_sub = [element if element not in search_list_none else -3300.0 for element in haplotype_1_copyratios_values]

            haplotype_2_copyratios_values_normal = [-3300.0 if element not in search_list else element for element in haplotype_2_copyratios_values]
            haplotype_2_copyratios_values_sub = [element if element not in search_list_none else -3300.0 for element in haplotype_2_copyratios_values]
            OFFSET = args.cut_threshold/150

            haplotype_1_copyratios_values_normal = [x if x == 'None' else x  for x in haplotype_1_copyratios_values_normal]
            haplotype_2_copyratios_values_normal = [x if x == 'None' else -x  for x in haplotype_2_copyratios_values_normal]
            haplotype_1_copyratios_values_sub = [x if (x == 'None' or x == -3300.0) else  x for x in haplotype_1_copyratios_values_sub]
            haplotype_2_copyratios_values_sub = [x if (x == 'None' or x == -3300.0) else  -x  for x in haplotype_2_copyratios_values_sub]

            if args.without_phasing:
                colors_cn = ['darkolivegreen']
                colors_sub = ['#57a456']
            else:
                colors_cn = ['firebrick', 'steelblue']
                colors_sub = ['#E95F0A', '#6EC5E9']
            name = "Copynumbers"
            add_scatter_trace_copyratios(args, fig, colors_cn, name, haplotype_1_copyratios_positions,
                                         haplotype_2_copyratios_positions, haplotype_1_copyratios_values_normal,
                                         haplotype_2_copyratios_values_normal, df_segs_hp1, df_segs_hp2, mul_cols=True, centers=centers)
            name = "Copynumbers (subclonal)"
            add_scatter_trace_copyratios(args, fig, colors_sub, name, haplotype_1_copyratios_positions,
                                         haplotype_2_copyratios_positions, haplotype_1_copyratios_values_sub,
                                         haplotype_2_copyratios_values_sub, df_segs_hp1, df_segs_hp2, mul_cols=True, centers=centers)
        else:
            if args.without_phasing:
                OFFSET = 0
                colors = ['darkolivegreen']
            else:
                OFFSET = args.cut_threshold/150
                colors = ['firebrick', 'steelblue']
            haplotype_1_copyratios_values = [x if x == 'None' else x  for x in haplotype_1_copyratios_values]
            haplotype_2_copyratios_values = [x if x == 'None' else -(x) for x in haplotype_2_copyratios_values]
            name = "Copynumbers"

            add_scatter_trace_copyratios(args, fig, colors, name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values, df_segs_hp1.start.values.tolist(), df_segs_hp1.end.values.tolist(), df_segs_hp2.start.values.tolist(), df_segs_hp2.end.values.tolist(), mul_cols=True)

    # #############################################################
    # #############################################################
    centers_rev = [-x for x in centers]
    centers_rev.reverse()
    integer_fractional_means_rev = [x for x in integer_fractional_centers]
    integer_fractional_means_rev.reverse()
    if args.without_phasing:
        tick_vals = centers
        tickt_ext = integer_fractional_centers
        tickvals = [i for i in range(0, 1000, 25)]
        ticktext = [str(abs(i)) for i in range(0, 1000, 25)]
        yaxis2_3_range = [0, args.cut_threshold + 5]
        plot_height = 850 + 150 + 20
        legend_y = 1.085
    else:
        tick_vals = centers_rev + centers
        tickt_ext = integer_fractional_means_rev + integer_fractional_centers
        tickvals = [i for i in range(-1000, 1000, 25)]
        ticktext = [str(abs(i)) for i in range(-1000, 1000, 25)]
        yaxis2_3_range = [-(args.cut_threshold + 5), args.cut_threshold + 5]
        plot_height = 850 + 150 + 130 + 40 + 20
        legend_y = 1.05
    # #############################################################
    # #############################################################
    #fig.update_yaxes(range=[-1, args.cut_threshold])
    fig.update_layout(yaxis=dict(title="<b>Breakpoints</b>", range=[0, 75], showticklabels = False, showgrid=False, zeroline=False),
                      yaxis2=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis3=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis5=dict(title="<b>B-allele frequency</b>", range=[0, 0.6], showticklabels=True, showgrid=False, zeroline=False),
                      yaxis6=dict(title="<b>Genes</b>", range=[0, 1], showticklabels = False, showgrid=False, zeroline=True, zerolinewidth=2, zerolinecolor='black'),
                      yaxis4=dict(title="<b>Cytos</b>", range=[0, 1], showticklabels=False, showgrid=False, zeroline=True, zerolinewidth=2, zerolinecolor='black'),

                      xaxis=dict(showspikes=True, tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], showticklabels = False, showgrid=False, zeroline=False),
                      xaxis2=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], zeroline=True, zerolinewidth=1, zerolinecolor='black', showgrid=False,),
                      xaxis3=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size],showgrid=False,),
                      xaxis5=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist()) * args.bin_size], showgrid=False,),
                      xaxis6=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], showgrid=False, zeroline=True, zerolinewidth=1, zerolinecolor='black'),
                      xaxis4=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist()) * args.bin_size], showgrid=False, zeroline=True, zerolinewidth=1, zerolinecolor='black'))

    fig.update_layout(
        yaxis2=dict(
            title="<b>Coverage depth</b> (per bin)",
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext
        ),
        yaxis3=dict(
            title="<b>Copies</b> (integers/fractions)",
            zeroline=True, zerolinewidth=1, zerolinecolor='black',
            tickmode='array',
            tickvals=tick_vals,
            ticktext=tickt_ext  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
        ),
    )

    fig.update_xaxes(
        #xaxis2=dict(
            tickangle=90,
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        #),
        #font=dict(size=18, color="black")
    )
    ax = 20
    ay = -30
    add_annotation(fig, 960000000, 75, ax, ay, "DEL", '#CF0759')
    add_annotation(fig, 960000000 + 250000000, 75, ax, ay, "INV", '#2830DE')
    add_annotation(fig, 960000000 + 250000000 + 250000000, 75, ax, ay, "INS", '#e0cf03')
    add_annotation(fig, 960000000 + 250000000 + 250000000 + 250000000, 75, ax, ay, "BND", '#737373')
    add_annotation(fig, 960000000 + 250000000 + 250000000 + 250000000 + 250000000, 75, ax, ay, "DUP", '#178117')

    if args.loh_enable:
        add_annotation(fig, 960000000 + 250000000 + 250000000, 0, ax, ay, "LOH", '#2980b9')
        #fig.add_annotation(text="Het SNPs ratio threshold: " + str(args.hets_ratio), x=960000000 + 250000000 + 250000000, y=args.cut_threshold, showarrow=False, row=2, col=1)

    # Update layout
    fig.update_layout(
        template="plotly_white",
        font_family="Times New Roman"
    )

    fig.update_layout(
        title=chrom,
    )
    # Legend
    fig.update_layout(legend=dict(
        orientation='h', xanchor="center", x=0.45, y=legend_y,  # orientation = 'v', xanchor = "center", x = 1.08, y= .5
    ))
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    #fig.update_xaxes(tick0=0.0, rangemode="nonnegative")
    fig.update_xaxes(spikedash='longdashdot', spikemode='across', spikesnap='cursor', showspikes=True)

    fig.update_layout(legend={'itemsizing': 'constant'})

    fig.update_layout(font_family="Times New Roman")

    if args.without_phasing:
        genome_tile = args.genome_name
    else:
        genome_tile = args.genome_name + "<br>" + "<span style='color:blue'>Ploidy: </span>" + str(
            args.tumor_ploidy) + "     " + "<span style='color:blue'>Cellular tumor fraction: </span>" + str(
            args.tumor_purity) + "     " + "<span style='color:blue'>Confidence: </span>" + str(p_value)
    fig.update_layout(
        title={
            'text': genome_tile,
            'y': 0.98,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family="Courier New",
        font_color="dimgray",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="green",
    )
    fig.update_layout(
        width=1380,
        height=plot_height,
    )
    if args.without_phasing:
        path_set = args.genome_name + '_genome_copynumbers_breakpoints_subclonal'
    else:
        if is_half:
            path_set = 'wgd/'+str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers_breakpoints_subclonal'
        else:
            path_set = str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers_breakpoints_subclonal'

    if args.pdf_enable:
        print_genome_pdf(fig, path_set, args.out_dir_plots)

    fig.write_html(args.out_dir_plots + '/' + path_set + ".html")

def copy_number_plots_genome_breakpoints_subclonal(centers, integer_fractional_centers, df_cnr_hp1, df_segs_hp1_, df_cnr_hp2, df_segs_hp2_, df_unphased, args, p_value, loh_regions, het_snps_bafs_means, is_half):
    # TODO genome-wide plots
    df_segs_hp1 = df_segs_hp1_.copy()
    df_segs_hp2 = df_segs_hp2_.copy()

    coords_single, coords_single_chr, coords, coords_chr, _,_ = sv_vcf_bps_cn_check(args.breakpoints, args)

    lengths = []
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)

    #chroms = ['chr1', 'chr5','chr8']

    df_cnr_hp1_ = []
    df_segs_hp1_ = []
    df_cnr_hp2_ = []
    df_segs_hp2_ = []
    df_genes_1_ = []
    df_genes_2_ = []
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    #cancer_genes = os.path.join(fileDir, args.cancer_genes)
    df_genes = csv_df_chromosomes_sorter(args.cancer_genes, ['chr','start','end','gene'])
    genestart_1 = []
    genestart_2 = []
    last_len = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        df_segs_hp1_.append(df_segs_hp1[df_segs_hp1['chromosome'] == chrom])
        df_cnr_hp2_.append(df_cnr_hp2[df_cnr_hp2['chr'] == chrom])
        df_segs_hp2_.append(df_segs_hp2[df_segs_hp2['chromosome'] == chrom])

    #coords_chr =  [['chr1',223550001, 'HP-1'],['chr6',60250001, 'HP-1'], ['chr6',27900001, 'HP-1'],['chr15',22200001, 'HP-1'], ['chr1', 144850001, 'HP-2'], ['chr9', 68300001, 'HP-2']]
    #coords = [['chr1',223550001],['chr6',60250001], ['chr6',27900001],['chr15',22200001], ['chr1', 144850001], ['chr9', 68300001]]
    #coords_single = [['chr1',22355007, 22355008], ['chr2', 123444555, 123444556], ['chr19', 1355007, 1355008]]
    #coords_single_chr = [['chr1', 22355007,22355008, 'DEL'], ['chr2', 123444555,123444556, 'DUP'], ['chr19', 1355007,1355008, 'DEL']]

    offset = 0
    for i, chrom in enumerate(chroms):
        new_len = len(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        # [[9:64526327,1:9769900],
        # [6:91468745,1:9769559]]
        if i > 0:
            for j, pair in enumerate(coords):
                if chrom == pair[0]:
                    coords[j][1] = pair[1] + offset#(last_len * (args.bin_size))

        if i > 0:
            for k, pair in enumerate(coords_single):
                if chrom == pair[0]:
                    coords_single[k][1] = pair[1] + offset  # (last_len * (args.bin_size))

        last_len += new_len
        offset += regions[i]

    df_cnr_hp1_1 = df_cnr_hp1
    df_1_ = []
    offset = 0
    for i, chrom in enumerate(chroms):
        df_genes_chrom_1 = df_cnr_hp1_1[df_cnr_hp1_1['chr'] == chrom]
        df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
        offset += regions[i]
        df_1_.append(df_genes_chrom_1)
    df_1 = pd.concat(df_1_)

    n = len(coords)
    coords = [x[1:] for x in coords]
    coords = [coords[i] + (coords[i + 1] if i + 1 < n else []) for i in range(0, n, 2)]
    coords = list(map(sorted, coords))

    last_len = 0
    coords_single = [j for i in [x[1:2] for x in coords_single] for j in i]

    df_cnr_hp1 = pd.concat(df_cnr_hp1_)
    df_segs_hp1 = pd.concat(df_segs_hp1_)

    df_cnr_hp2 = pd.concat(df_cnr_hp2_)
    df_segs_hp2 = pd.concat(df_segs_hp2_)

    df_genes_1 = df_genes[df_genes['start'] < df_genes['end']]
    df_genes_2 = df_genes[df_genes['start'] > df_genes['end']]

    last_len = 0
    offset = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        new_len = len(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        if not chrom.startswith('chr'):
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == 'chr' + chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == 'chr' + chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        else:
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        last_len += new_len
        offset += regions[i]

    df_genes_1 = pd.concat(df_genes_1_)
    df_genes_2 = pd.concat(df_genes_2_)

    genename_1 = df_genes_1['gene'].values.tolist()
    genename_2 = df_genes_2['gene'].values.tolist()

    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_cnr_hp1)*args.bin_size, args.bin_size, dtype=int)

    ###########################################################

    arcs_data = get_all_breakpoints_data(coords, coords_chr, 75, args.breakpoints, args)
    arcs_data1 = get_all_breakpoints_data(coords, coords_chr, -75, args.breakpoints, args)
    # fig = go.Figure(data=arcs_data)

    #fig = go.Figure()
    #fig = go.Figure().set_subplots(rows=2, cols=1)
    if args.without_phasing:
        row_heights = [220, 160, 150, 40]
    else:
        row_heights = [220, 320, 150, 40]

    #fig = make_subplots(rows=2, cols=1, shared_yaxes=False, shared_xaxes=True, vertical_spacing=0.02, horizontal_spacing=0.02)
    fig = make_subplots(rows=4, cols=1, shared_yaxes=False, shared_xaxes='columns',  vertical_spacing=0.01, row_heights=row_heights,
                        horizontal_spacing=0.02, specs=[[{"type":"xy"}], [{"secondary_y":True}], [{"type":"xy"}], [{"type":"xy"}]]) #https://community.plotly.com/t/can-subplot-support-multiple-y-axes/38891/20
    # #############################################################
    for i in range(len(arcs_data)):
        fig.add_trace(go.Scatter(x=arcs_data[i][0],
                                 y=arcs_data[i][1],
                                 name=arcs_data[i][2],
                                 mode=arcs_data[i][3],
                                 line=arcs_data[i][4],
                                 yaxis="y",
                                 opacity=0.5,
                                 text = arcs_data[i][5],
                                 hovertemplate= '<br><b>Breakpoint info: </b>' + arcs_data[i][6],
                                 showlegend=arcs_data[i][7]), row=1, col=1,
                     )
    # #############################################################

    # #############################################################
    custom_text_data_hp1 = df_cnr_hp1.start.values.tolist()
    custom_text_data_hp2 = df_cnr_hp2.start.values.tolist()
    if args.without_phasing:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.coverage.values.tolist(), name='Unphased coverage',
                                   text=custom_text_data_hp1,
                                   yaxis="y2", opacity=0.7, color='olive', visibility=True, mul_cols=True)
    else:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.hp1.values.tolist(), name='HP-1 coverage',
                                   text=custom_text_data_hp1,
                                   yaxis=None, opacity=0.7, color='firebrick', visibility=True, mul_cols=True)
        add_scatter_trace_coverage(fig, df_1['start'], [-x for x in df_cnr_hp2.hp2.values.tolist()], name='HP-2 coverage',
                                   text=custom_text_data_hp2,
                                   yaxis=None, opacity=0.7, color='steelblue', visibility=True, mul_cols=True)
        # if not args.unphased_reads_coverage_disable:
        #     add_scatter_trace_coverage(fig, indices, df_unphased.hp3.values.tolist(), name='Unphased', text=None,
        #                                yaxis=None, opacity=0.7, color='olive', visibility='legendonly', mul_cols=True)

    add_scatter_trace_coverage(fig, df_1['start'], het_snps_bafs_means.vaf.values.tolist(), name='BAF',
                                   text=het_snps_bafs_means.pos.values.tolist(),
                                   yaxis="y4", opacity=0.7, color='olive', visibility='legendonly', mul_cols=True, row=3, baf=True)

    fig.add_trace(go.Scatter(x=df_genes_1['start'], y=[0.5] * len(df_genes_1), #[1, 2, 3, 4, 5] * (len(df_genes_1)//5),
                             mode='markers',
                             text = genename_1,
                             customdata = genestart_1,
                             #hovertext=df_genes['gene'],
                             hovertemplate=
                             '<br><b>Gene</b>: %{text}'+
                             '<br><b>Pos</b>: %{customdata}<br>',
                             marker=dict(
                                 symbol="y-left",
                                 color="#3A6B35",
                                 size=6,
                                 line=dict(width=1, color="#7F7F7F"),
                             ),
                             yaxis="y5",
                             name= 'GeneInfo',
                             showlegend=False,), row=4, col=1,)

    # fig.add_trace(go.Scatter(x=df_genes_2['start'], y=[1, 2, 3, 4, 5] * (len(df_genes_2)//5),
    #                          mode='markers',
    #                          text = genename_2,
    #                          customdata = genestart_2,
    #                          #hovertext=df_genes['gene'],
    #                          hovertemplate=
    #                          '<br><b>Gene</b>: %{text}'+
    #                          '<br><b>Pos</b>: %{customdata}<br>',
    #                          marker=dict(
    #                              symbol="y-right",
    #                              color="#E3B448",
    #                              size=6,
    #                              line=dict(width=1, color="#7F7F7F"),
    #                          ),
    #                          yaxis="y4",
    #                          name= 'GeneInfo',
    #                          showlegend=False,), row=3, col=1,)
    current = 0
    label_pos = []
    label_pos_chrms = []
    label_chrms = []

    #coords_lines = [j for i in coords for j in i] + coords_single
    #for i, co in enumerate(coords_single):
        # if coords_single_chr[i][3] == 'DEL':
        #     fig.add_vline(x=co,  line_width=1, line_dash="solid", opacity=0.3, line_color="#CF0759", row=1, col=1,)
        # elif coords_single_chr[i][3] == 'DUP' or coords_single_chr[i][3] == 'INS':
        #     fig.add_vline(x=co,  line_width=1, line_dash="solid", opacity=0.3, line_color="#178117", row=1, col=1,)
        # else:

        #fig.add_vline(x=co,  line_width=1, line_dash="solid", opacity=0.3, line_color="#2B40A0", row=1, col=1,)

    fig.add_hline(y=2,  line_width=1, line_dash="solid", line_color="black", row=3, col=1,)
    fig.add_hline(y=-(args.cut_threshold + 5),  line_width=1, line_dash="solid", line_color="black", row=2, col=1,)

    # for index, chrom in enumerate(chroms):
    #     current += lengths[index]
    #     label_pos.append(round(current - (lengths[index] / 2))*args.bin_size) #y0=-10, y1=args.cut_threshold
    #
    #     label_pos_chrms.append(round(current - (lengths[index]))*args.bin_size)
    #     label_pos_chrms.append(round(current - (lengths[index]*.75)) * args.bin_size)
    #     label_pos_chrms.append(round(current - (lengths[index] / 2))*args.bin_size)
    #     label_pos_chrms.append(round(current - (lengths[index]*.25)) * args.bin_size)
    #
    #     label_chrms.append('0')
    #     label_chrms.append(str((round((lengths[index] * .25) * args.bin_size)) // 1000000) + 'M')
    #     label_chrms.append(str(((lengths[index]//2)*args.bin_size)// 1000000)+'M')
    #     label_chrms.append(str((round((lengths[index] * .75) * args.bin_size)) // 1000000) + 'M')
    #
    #     fig.add_vline(x=current*args.bin_size,  line_width=2, line_dash="solid", line_color="black", row=2, col=1,)

    current = 0
    label_pos = []
    offset_chroms = 0
    offset_chroms_1 = 0
    chroms = get_contigs_list(args.contigs)
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])

    for index, chrom in enumerate(chroms):
        loh_starts = []
        if not loh_regions.empty:
            df_loh_region = loh_regions[loh_regions['chr'] == chrom]
            loh_starts = df_loh_region.start.values.tolist()
            loh_ends = df_loh_region.end.values.tolist()
        if not df_centm.empty:
            df_cent_region = df_centm[df_centm['chr'] == chrom]
            cent_starts = df_cent_region.start.values.tolist()
            cent_ends = df_cent_region.end.values.tolist()

        current += lengths[index]
        label_pos.append(round(offset_chroms+regions[index]//2))

        fig.add_vline(x=offset_chroms,  line_width=1, line_dash="solid", line_color="#D7DBDD")

        #label_chrms.append('0')
        #label_chrms.append(str(((lengths[index]//2) + regions[index])// 1000000)+'M')

        offset_chroms_1 += regions[index]

        if index % 2 == 0:
            fig.add_vrect(x0=offset_chroms, x1=offset_chroms_1, fillcolor="#E5E7E9", opacity=0.9, layer="below", line_width=0, )

        if len(loh_starts) and args.loh_enable:
            for i in range(len(df_loh_region.start.values.tolist())):
                fig.add_vrect(x0=offset_chroms+loh_starts[i], x1=offset_chroms+loh_ends[i], fillcolor="#2980b9", opacity=0.3, layer="below", line_width=0, row=2, col=1,)
        if len(cent_starts) and args.loh_enable:
            fig.add_vrect(x0=offset_chroms+cent_starts[0], x1=offset_chroms+cent_ends[0], fillcolor="#7e1f14", opacity=0.3, layer="above", line_width=0, row=2, col=1,)

        offset_chroms += regions[index]

    if not args.copynumbers_disable:
        offset_start = 0
        haplotype_1_start_values = []
        haplotype_1_end_values = []

        haplotype_2_start_values = []
        haplotype_2_end_values = []
        for index, chrom in enumerate(chroms):
            df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
            df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
            haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
            haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()
            haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
            haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()
            if chrom == args.contigs.split('-')[0]:
                haplotype_1_start_values.extend(haplotype_1_start_values_copyratios)
                haplotype_1_end_values.extend(haplotype_1_end_values_copyratios)
                haplotype_2_start_values.extend(haplotype_2_start_values_copyratios)
                haplotype_2_end_values.extend(haplotype_2_end_values_copyratios)
            else:
                offset_start += regions[index-1]#lengths[index-1] * args.bin_size
                haplotype_1_start_values.extend([x + offset_start for x in haplotype_1_start_values_copyratios])
                haplotype_1_end_values.extend([x + offset_start for x in haplotype_1_end_values_copyratios])
                haplotype_2_start_values.extend([x + offset_start for x in haplotype_2_start_values_copyratios])
                haplotype_2_end_values.extend([x + offset_start for x in haplotype_2_end_values_copyratios])

        haplotype_1_gaps_values = np.full(len(df_segs_hp1.state.values.tolist()), 'None')
        haplotype_1_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp1.state.values.tolist(), df_segs_hp1.state.values.tolist(), haplotype_1_gaps_values)))
        haplotype_1_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_1_start_values, haplotype_1_end_values, haplotype_1_gaps_values)))

        haplotype_2_gaps_values = np.full(len(df_segs_hp2.state.values.tolist()), 'None')
        haplotype_2_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp2.state.values.tolist(), df_segs_hp2.state.values.tolist(), haplotype_2_gaps_values)))
        haplotype_2_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_2_start_values, haplotype_2_end_values, haplotype_2_gaps_values)))

        if args.copynumbers_subclonal_enable:

            search_list = centers #[centers[i] for i in [integer_fractional_centers.index(x) for x in integer_fractional_centers if not float(x).is_integer()]] #[27.025, 42.025]
            search_list_none = search_list + ['None']

            haplotype_1_copyratios_values_normal = [-3300.0 if element not in search_list else element for element in haplotype_1_copyratios_values]
            haplotype_1_copyratios_values_sub = [element if element not in search_list_none else -3300.0 for element in haplotype_1_copyratios_values]

            haplotype_2_copyratios_values_normal = [-3300.0 if element not in search_list else element for element in haplotype_2_copyratios_values]
            haplotype_2_copyratios_values_sub = [element if element not in search_list_none else -3300.0 for element in haplotype_2_copyratios_values]
            OFFSET = args.cut_threshold/150

            haplotype_1_copyratios_values_normal = [x if x == 'None' else x  for x in haplotype_1_copyratios_values_normal]
            haplotype_2_copyratios_values_normal = [x if x == 'None' else -x  for x in haplotype_2_copyratios_values_normal]
            haplotype_1_copyratios_values_sub = [x if (x == 'None' or x == -3300.0) else  x for x in haplotype_1_copyratios_values_sub]
            haplotype_2_copyratios_values_sub = [x if (x == 'None' or x == -3300.0) else  -x  for x in haplotype_2_copyratios_values_sub]

            if args.without_phasing:
                colors_cn = ['darkolivegreen']
                colors_sub = ['#57a456']
            else:
                colors_cn = ['firebrick', 'steelblue']
                colors_sub = ['#E95F0A', '#6EC5E9']
            name = "Copynumbers"
            add_scatter_trace_copyratios(args, fig, colors_cn, name, haplotype_1_copyratios_positions,
                                         haplotype_2_copyratios_positions, haplotype_1_copyratios_values_normal,
                                         haplotype_2_copyratios_values_normal, df_segs_hp1, df_segs_hp2, mul_cols=True, centers=centers)
            name = "Copynumbers (subclonal)"
            add_scatter_trace_copyratios(args, fig, colors_sub, name, haplotype_1_copyratios_positions,
                                         haplotype_2_copyratios_positions, haplotype_1_copyratios_values_sub,
                                         haplotype_2_copyratios_values_sub, df_segs_hp1, df_segs_hp2, mul_cols=True, centers=centers)
        else:
            if args.without_phasing:
                OFFSET = 0
                colors = ['darkolivegreen']
            else:
                OFFSET = args.cut_threshold/150
                colors = ['firebrick', 'steelblue']
            haplotype_1_copyratios_values = [x if x == 'None' else x  for x in haplotype_1_copyratios_values]
            haplotype_2_copyratios_values = [x if x == 'None' else -(x) for x in haplotype_2_copyratios_values]
            name = "Copynumbers"

            add_scatter_trace_copyratios(args, fig, colors, name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values, df_segs_hp1.start.values.tolist(), df_segs_hp1.end.values.tolist(), df_segs_hp2.start.values.tolist(), df_segs_hp2.end.values.tolist(), mul_cols=True)

    # #############################################################
    # #############################################################
    centers_rev = [-x for x in centers]
    centers_rev.reverse()
    integer_fractional_means_rev = [x for x in integer_fractional_centers]
    integer_fractional_means_rev.reverse()
    if args.without_phasing:
        tick_vals = centers
        tickt_ext = integer_fractional_centers
        tickvals = [i for i in range(0, 1000, 25)]
        ticktext = [str(abs(i)) for i in range(0, 1000, 25)]
        yaxis2_3_range = [0, args.cut_threshold + 5]
        plot_height = 850 + 150 + 20 + 20
        legend_y = 1.085
    else:
        tick_vals = centers_rev + centers
        tickt_ext = integer_fractional_means_rev + integer_fractional_centers
        tickvals = [i for i in range(-1000, 1000, 25)]
        ticktext = [str(abs(i)) for i in range(-1000, 1000, 25)]
        yaxis2_3_range = [-(args.cut_threshold + 5), args.cut_threshold + 5]
        plot_height = 850 + 150 + 130 + 20
        legend_y = 1.06
    # #############################################################
    # #############################################################
    #fig.update_yaxes(range=[-1, args.cut_threshold])
    fig.update_layout(yaxis=dict(title="<b>Breakpoints</b>", range=[0, 75], showticklabels = False, showgrid=False, zeroline=False),
                      yaxis2=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis3=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis4=dict(title="<b>B-allele frequency</b>", range=[0, 0.6], showticklabels=True, showgrid=False, zeroline=False),
                      yaxis5=dict(title="<b>Genes</b>", range=[0, 1], showticklabels = False, showgrid=False, zeroline=True, zerolinewidth=2, zerolinecolor='black'),

                      xaxis=dict(showspikes=True, tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], showticklabels = False, showgrid=False, zeroline=False),
                      xaxis2=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], zeroline=True, zerolinewidth=1, zerolinecolor='black', showgrid=False,),
                      xaxis3=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size],showgrid=False,),
                      xaxis4=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist()) * args.bin_size], showgrid=False,),
                      xaxis5=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], showgrid=False, zeroline=True, zerolinewidth=1, zerolinecolor='black'))

    fig.update_layout(
        yaxis2=dict(
            title="<b>Coverage depth</b> (per bin)",
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext
        ),
        yaxis3=dict(
            title="<b>Copies</b> (integers/fractions)",
            zeroline=True, zerolinewidth=1, zerolinecolor='black',
            tickmode='array',
            tickvals=tick_vals,
            ticktext=tickt_ext  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
        ),
    )

    fig.update_xaxes(
        #xaxis2=dict(
            tickangle=90,
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        #),
        #font=dict(size=18, color="black")
    )
    fig.add_annotation(text="DEL",xref="paper",yref="paper",x=0.25,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#CF0759',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)
    fig.add_annotation(text="INV",xref="paper",yref="paper",x=0.30,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#2830DE',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)
    fig.add_annotation(text="INS",xref="paper",yref="paper",x=0.36,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#e0cf03',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)
    fig.add_annotation(text="BND",xref="paper",yref="paper",x=0.42,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#737373',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)
    fig.add_annotation(text="DUP",xref="paper",yref="paper",x=0.47,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#178117',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)

    fig.add_annotation(text="LOH Regions",xref="paper",yref="paper",x=0.60,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#2980b9',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)
    fig.add_annotation(text="Centromeres",xref="paper",yref="paper",x=0.75,y=1.03,showarrow=False,font=dict(size=17, color="white"),bgcolor='#7e1f14',bordercolor="#c7c7c7",borderwidth=2,borderpad=4, opacity=0.7,)

    # Update layout
    fig.update_layout(
        template="plotly_white",
        font_family="Times New Roman"
    )

    fig.update_layout(
        title=chrom,
    )
    # Legend
    fig.update_layout(legend=dict(
        orientation='h', xanchor="center", x=0.45, y=legend_y,  # orientation = 'v', xanchor = "center", x = 1.08, y= .5
    ))
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    #fig.update_xaxes(tick0=0.0, rangemode="nonnegative")

    fig.update_layout(legend={'itemsizing': 'constant'})

    fig.update_layout(font_family="Times New Roman")

    if args.without_phasing:
        genome_tile = args.genome_name
    else:
        genome_tile = args.genome_name + "<br>" + "<span style='color:blue'>Ploidy: </span>" + str(
            args.tumor_ploidy) + "     " + "<span style='color:blue'>Cellular tumor fraction: </span>" + str(
            args.tumor_purity) + "     " + "<span style='color:blue'>Confidence: </span>" + str(p_value)
    fig.update_layout(
        title={
            'text': genome_tile,
            'y': 0.98,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family="Courier New",
        font_color="dimgray",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="green",
    )
    fig.update_layout(
        width=1380,
        height=plot_height,
    )
    if args.without_phasing:
        path_set = args.genome_name + '_genome_copynumbers_breakpoints_subclonal'
    else:
        if is_half:
            path_set = 'wgd/'+str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers_breakpoints_subclonal'
        else:
            path_set = str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers_breakpoints_subclonal'

    if args.pdf_enable:
        print_genome_pdf(fig, path_set, args.out_dir_plots)

    fig.write_html(args.out_dir_plots + '/' + path_set + ".html")

def copy_number_plots_genome_subclonal(centers, integer_fractional_centers, df_cnr_hp1, df_segs_hp1_, df_cnr_hp2, df_segs_hp2_, df_unphased, args, p_value, loh_regions, het_snps_bafs_means, is_half):
    # TODO genome-wide plots
    df_segs_hp1 = df_segs_hp1_.copy()
    df_segs_hp2 = df_segs_hp2_.copy()

    lengths = []
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)

    #chroms = ['chr1', 'chr5','chr8']

    df_cnr_hp1_ = []
    df_segs_hp1_ = []
    df_cnr_hp2_ = []
    df_segs_hp2_ = []
    df_genes_1_ = []
    df_genes_2_ = []
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    #cancer_genes = os.path.join(fileDir, args.cancer_genes)
    df_genes = csv_df_chromosomes_sorter(args.cancer_genes, ['chr','start','end','gene'])
    genestart_1 = []
    genestart_2 = []
    last_len = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        df_segs_hp1_.append(df_segs_hp1[df_segs_hp1['chromosome'] == chrom])
        df_cnr_hp2_.append(df_cnr_hp2[df_cnr_hp2['chr'] == chrom])
        df_segs_hp2_.append(df_segs_hp2[df_segs_hp2['chromosome'] == chrom])

    df_cnr_hp1_1 = df_cnr_hp1
    df_1_ = []
    offset = 0
    for i, chrom in enumerate(chroms):
        df_genes_chrom_1 = df_cnr_hp1_1[df_cnr_hp1_1['chr'] == chrom]
        df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
        offset += regions[i]
        df_1_.append(df_genes_chrom_1)
    df_1 = pd.concat(df_1_)

    df_cnr_hp1 = pd.concat(df_cnr_hp1_)
    df_segs_hp1 = pd.concat(df_segs_hp1_)

    df_cnr_hp2 = pd.concat(df_cnr_hp2_)
    df_segs_hp2 = pd.concat(df_segs_hp2_)

    df_genes_1 = df_genes[df_genes['start'] < df_genes['end']]
    df_genes_2 = df_genes[df_genes['start'] > df_genes['end']]

    last_len = 0
    offset = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        new_len = len(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        if not chrom.startswith('chr'):
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == 'chr' + chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == 'chr' + chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        else:
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        last_len += new_len
        offset += regions[i]

    df_genes_1 = pd.concat(df_genes_1_)
    df_genes_2 = pd.concat(df_genes_2_)

    genename_1 = df_genes_1['gene'].values.tolist()
    genename_2 = df_genes_2['gene'].values.tolist()

    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_cnr_hp1)*args.bin_size, args.bin_size, dtype=int)

    ###########################################################
    if args.without_phasing:
        row_heights = [160, 150, 40]
    else:
        row_heights = [320, 150, 40]

    # fig = make_subplots(rows=2, cols=1, shared_yaxes=False, shared_xaxes=True, vertical_spacing=0.02, horizontal_spacing=0.02)
    fig = make_subplots(rows=3, cols=1, shared_yaxes=False, shared_xaxes='columns', vertical_spacing=0.01, row_heights=row_heights,
                        horizontal_spacing=0.02, specs=[[{"secondary_y": True}], [{"type": "xy"}], [{"type": "xy"}]])  # https://community.plotly.com/t/can-subplot-support-multiple-y-axes/38891/20
    # #############################################################

    # #############################################################

    # #############################################################
    custom_text_data_hp1 = df_cnr_hp1.start.values.tolist()
    custom_text_data_hp2 = df_cnr_hp2.start.values.tolist()
    if args.without_phasing:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.coverage.values.tolist(), name='Unphased coverage',
                                   text=custom_text_data_hp1,
                                   yaxis="y", opacity=0.7, color='olive', visibility=True, mul_cols=True, row=1)
    else:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.hp1.values.tolist(), name='HP-1 coverage',
                                   text=custom_text_data_hp1,
                                   yaxis="y", opacity=0.7, color='firebrick', visibility=True, mul_cols=True, row=1)
        add_scatter_trace_coverage(fig, df_1['start'], [-x for x in df_cnr_hp2.hp2.values.tolist()], name='HP-2 coverage',
                                   text=custom_text_data_hp2,
                                   yaxis="y", opacity=0.7, color='steelblue', visibility=True, mul_cols=True, row=1)
        # if not args.unphased_reads_coverage_disable:
        #     add_scatter_trace_coverage(fig, indices, df_unphased.hp3.values.tolist(), name='Unphased', text=None,
        #                                yaxis=None, opacity=0.7, color='olive', visibility='legendonly', mul_cols=True)
    add_scatter_trace_coverage(fig, df_1['start'], het_snps_bafs_means.vaf.values.tolist(), name='BAF',
                                   text=het_snps_bafs_means.pos.values.tolist(),
                                   yaxis="y3", opacity=0.7, color='olive', visibility='legendonly', mul_cols=True, row=2, baf=True)

    fig.add_trace(go.Scatter(x=df_genes_1['start'], y=[0.5] * len(df_genes_1), #[1, 2, 3, 4, 5] * (len(df_genes_1)//5),
                             mode='markers',
                             text = genename_1,
                             customdata = genestart_1,
                             #hovertext=df_genes['gene'],
                             hovertemplate=
                             '<br><b>Gene</b>: %{text}'+
                             '<br><b>Pos</b>: %{customdata}<br>',
                             marker=dict(
                                 symbol="y-left",
                                 color="#3A6B35",
                                 size=6,
                                 line=dict(width=1, color="#7F7F7F"),
                             ),
                             yaxis="y4",
                             name= 'GeneInfo',
                             showlegend=False,), row=3, col=1,)

    fig.add_hline(y=2,  line_width=1, line_dash="solid", line_color="black", row=2, col=1,)
    fig.add_hline(y=-(args.cut_threshold + 5),  line_width=1, line_dash="solid", line_color="black", row=1, col=1,)

    current = 0
    label_pos = []
    offset_chroms = 0
    offset_chroms_1 = 0
    chroms = get_contigs_list(args.contigs)
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])

    for index, chrom in enumerate(chroms):
        loh_starts = []
        if not loh_regions.empty:
            df_loh_region = loh_regions[loh_regions['chr'] == chrom]
            loh_starts = df_loh_region.start.values.tolist()
            loh_ends = df_loh_region.end.values.tolist()
        if not df_centm.empty:
            df_cent_region = df_centm[df_centm['chr'] == chrom]
            cent_starts = df_cent_region.start.values.tolist()
            cent_ends = df_cent_region.end.values.tolist()

        current += lengths[index]
        label_pos.append(round(offset_chroms+regions[index]//2))

        fig.add_vline(x=offset_chroms,  line_width=1, line_dash="solid", line_color="#D7DBDD")

        offset_chroms_1 += regions[index]

        if index % 2 == 0:
            fig.add_vrect(x0=offset_chroms, x1=offset_chroms_1, fillcolor="#E5E7E9", opacity=0.9, layer="below", line_width=0, )

        if len(loh_starts) and args.loh_enable:
            for i in range(len(df_loh_region.start.values.tolist())):
                fig.add_vrect(x0=offset_chroms+loh_starts[i], x1=offset_chroms+loh_ends[i], fillcolor="#2980b9", opacity=0.3, layer="below", line_width=0, row=1, col=1,)
        if len(cent_starts) and args.loh_enable:
            fig.add_vrect(x0=offset_chroms+cent_starts[0], x1=offset_chroms+cent_ends[0], fillcolor="#7e1f14", opacity=0.3, layer="above", line_width=0, row=1, col=1,)

        offset_chroms += regions[index]

    if not args.copynumbers_disable:
        offset_start = 0
        haplotype_1_start_values = []
        haplotype_1_end_values = []

        haplotype_2_start_values = []
        haplotype_2_end_values = []
        for index, chrom in enumerate(chroms):
            df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
            df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
            haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
            haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()
            haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
            haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()
            if chrom == args.contigs.split('-')[0]:
                haplotype_1_start_values.extend(haplotype_1_start_values_copyratios)
                haplotype_1_end_values.extend(haplotype_1_end_values_copyratios)
                haplotype_2_start_values.extend(haplotype_2_start_values_copyratios)
                haplotype_2_end_values.extend(haplotype_2_end_values_copyratios)
            else:
                offset_start += regions[index-1]#lengths[index-1] * args.bin_size
                haplotype_1_start_values.extend([x + offset_start for x in haplotype_1_start_values_copyratios])
                haplotype_1_end_values.extend([x + offset_start for x in haplotype_1_end_values_copyratios])
                haplotype_2_start_values.extend([x + offset_start for x in haplotype_2_start_values_copyratios])
                haplotype_2_end_values.extend([x + offset_start for x in haplotype_2_end_values_copyratios])

        haplotype_1_gaps_values = np.full(len(df_segs_hp1.state.values.tolist()), 'None')
        haplotype_1_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp1.state.values.tolist(), df_segs_hp1.state.values.tolist(), haplotype_1_gaps_values)))
        haplotype_1_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_1_start_values, haplotype_1_end_values, haplotype_1_gaps_values)))

        haplotype_2_gaps_values = np.full(len(df_segs_hp2.state.values.tolist()), 'None')
        haplotype_2_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp2.state.values.tolist(), df_segs_hp2.state.values.tolist(), haplotype_2_gaps_values)))
        haplotype_2_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_2_start_values, haplotype_2_end_values, haplotype_2_gaps_values)))

        if args.copynumbers_subclonal_enable:

            search_list = centers #[centers[i] for i in [integer_fractional_centers.index(x) for x in integer_fractional_centers if not float(x).is_integer()]] #[27.025, 42.025]
            search_list_none = search_list + ['None']

            haplotype_1_copyratios_values_normal = [-3300.0 if element not in search_list else element for element in haplotype_1_copyratios_values]
            haplotype_1_copyratios_values_sub = [element if element not in search_list_none else -3300.0 for element in haplotype_1_copyratios_values]

            haplotype_2_copyratios_values_normal = [-3300.0 if element not in search_list else element for element in haplotype_2_copyratios_values]
            haplotype_2_copyratios_values_sub = [element if element not in search_list_none else -3300.0 for element in haplotype_2_copyratios_values]
            OFFSET = args.cut_threshold/150

            haplotype_1_copyratios_values_normal = [x if x == 'None' else x  for x in haplotype_1_copyratios_values_normal]
            haplotype_2_copyratios_values_normal = [x if x == 'None' else -x  for x in haplotype_2_copyratios_values_normal]
            haplotype_1_copyratios_values_sub = [x if x == 'None' else x  for x in haplotype_1_copyratios_values_sub]
            haplotype_2_copyratios_values_sub = [x if x == 'None' else -x  for x in haplotype_2_copyratios_values_sub]

            if args.without_phasing:
                colors_cn = ['darkolivegreen']
                colors_sub = ['#57a456']
            else:
                colors_cn = ['firebrick', 'steelblue']
                colors_sub = ['#E95F0A', '#6EC5E9']
            name = "Copynumbers"
            add_scatter_trace_copyratios(args, fig, colors_cn, name, haplotype_1_copyratios_positions,
                                         haplotype_2_copyratios_positions, haplotype_1_copyratios_values_normal,
                                         haplotype_2_copyratios_values_normal, df_segs_hp1, df_segs_hp2, mul_cols=True, row=1)
            name = "Copynumbers (subclonal)"
            add_scatter_trace_copyratios(args, fig, colors_sub, name, haplotype_1_copyratios_positions,
                                         haplotype_2_copyratios_positions, haplotype_1_copyratios_values_sub,
                                         haplotype_2_copyratios_values_sub, df_segs_hp1, df_segs_hp2, mul_cols=True, row=1)
        else:
            if args.without_phasing:
                OFFSET = 0
                colors = ['darkolivegreen']
            else:
                OFFSET = args.cut_threshold/150
                colors = ['firebrick', 'steelblue']
            haplotype_1_copyratios_values = [x if x == 'None' else x  for x in haplotype_1_copyratios_values]
            haplotype_2_copyratios_values = [x if x == 'None' else -(x) for x in haplotype_2_copyratios_values]
            name = "Copynumbers"

            add_scatter_trace_copyratios(args, fig, colors, name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values, df_segs_hp1.start.values.tolist(), df_segs_hp1.end.values.tolist(), df_segs_hp2.start.values.tolist(), df_segs_hp2.end.values.tolist(), mul_cols=True, row=1)

    # #############################################################
    # #############################################################
    centers_rev = [-x for x in centers]
    centers_rev.reverse()
    integer_fractional_means_rev = [x for x in integer_fractional_centers]
    integer_fractional_means_rev.reverse()
    if args.without_phasing:
        tick_vals = centers
        tickt_ext = integer_fractional_centers
        tickvals = [i for i in range(0, 1000, 25)]
        ticktext = [str(abs(i)) for i in range(0, 1000, 25)]
        yaxis2_3_range = [0, args.cut_threshold + 5]
        plot_height = 630 + 150 + 20
        legend_y = 1.085
    else:
        tick_vals = centers_rev + centers
        tickt_ext = integer_fractional_means_rev + integer_fractional_centers
        tickvals = [i for i in range(-1000, 1000, 25)]
        ticktext = [str(abs(i)) for i in range(-1000, 1000, 25)]
        yaxis2_3_range = [-(args.cut_threshold + 5), args.cut_threshold + 5]
        plot_height = 630 + 150 + 130
        legend_y = 1.06
    # #############################################################
    # #############################################################
    # fig.update_yaxes(range=[-1, args.cut_threshold])
    fig.update_layout(
                      yaxis=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis2=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis3=dict(title="<b>B-allele frequency</b>", range=[0, 0.6], showticklabels=True, showgrid=False, zeroline=False),
                      yaxis4=dict(title="<b>Genes</b>", range=[0, 1], showticklabels = False, showgrid=False, zeroline=True, zerolinewidth=2, zerolinecolor='black'),

                      xaxis=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], zeroline=True, zerolinewidth=1, zerolinecolor='black', showgrid=False,),
                      xaxis2=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size],showgrid=False,),
                      xaxis3=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist()) * args.bin_size], showgrid=False,),
                      xaxis4=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], showgrid=False, zeroline=True, zerolinewidth=1, zerolinecolor='black'))

    fig.update_layout(
        yaxis=dict(
            title="<b>Coverage depth</b> (per bin)",
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext
        ),
        yaxis2=dict(
            title="<b>Copies</b> (integers/fractions)",
            zeroline=True, zerolinewidth=1, zerolinecolor='black',
            tickmode='array',
            tickvals=tick_vals,
            ticktext=tickt_ext  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
        ),
    )

    fig.update_xaxes(
        #xaxis2=dict(
            tickangle=90,
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        #),
        #font=dict(size=18, color="black")
    )
    ax = 20
    ay = -30
    if args.loh_enable:
        add_annotation(fig, 960000000 + 250000000 + 250000000, 0, ax, ay, "LOH", '#2980b9')
        add_annotation(fig, 960000000 + 250000000 + 250000000 + 250000000 + 250000000, 0, ax, ay, "Centromeres", '#7e1f14')
        #fig.add_annotation(text="Het SNPs ratio threshold: " + str(args.hets_ratio), x=960000000 + 250000000 + 250000000, y=args.cut_threshold, showarrow=False, row=2, col=1)

    # Update layout
    fig.update_layout(
        template="plotly_white",
        font_family="Times New Roman"
    )

    fig.update_layout(
        title=chrom,
    )
    # Legend
    fig.update_layout(legend=dict(
        orientation='h', xanchor="center", x=0.45, y=legend_y,  # orientation = 'v', xanchor = "center", x = 1.08, y= .5
    ))
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    #fig.update_xaxes(tick0=0.0, rangemode="nonnegative")

    fig.update_layout(legend={'itemsizing': 'constant'})

    fig.update_layout(font_family="Times New Roman")

    if args.without_phasing:
        genome_tile = args.genome_name
    else:
        genome_tile = args.genome_name + "<br>" + "<span style='color:blue'>Ploidy: </span>" + str(
            args.tumor_ploidy) + "     " + "<span style='color:blue'>Cellular tumor fraction: </span>" + str(
            args.tumor_purity) + "     " + "<span style='color:blue'>Confidence: </span>" + str(p_value)
    fig.update_layout(
        title={
            'text': genome_tile,
            'y': 0.98,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family="Courier New",
        font_color="dimgray",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="green",
    )
    fig.update_layout(
        width=1380,
        height=plot_height,
    )
    if args.without_phasing:
        path_set = args.genome_name + '_genome_copynumbers_subclonal'
    else:
        if is_half:
            path_set = 'wgd/'+str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers_subclonal'
        else:
            path_set = str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genome_copynumbers_subclonal'

    if args.pdf_enable:
        print_genome_pdf(fig, path_set, args.out_dir_plots)

    fig.write_html(args.out_dir_plots + '/' + path_set + ".html")

def genes_copy_number_plots_genome(df_genes, centers, integer_fractional_centers, df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, df_unphased, args, p_value, loh_regions, het_snps_bafs_means, is_half):
    # TODO genome-wide plots
    lengths = []
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)

    df_cnr_hp1_ = []
    df_segs_hp1_ = []
    df_cnr_hp2_ = []
    df_segs_hp2_ = []
    df_genes_1_ = []
    df_genes_2_ = []

    genestart_1 = []
    genestart_2 = []
    last_len = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        df_segs_hp1_.append(df_segs_hp1[df_segs_hp1['chromosome'] == chrom])
        df_cnr_hp2_.append(df_cnr_hp2[df_cnr_hp2['chr'] == chrom])
        df_segs_hp2_.append(df_segs_hp2[df_segs_hp2['chromosome'] == chrom])

    df_cnr_hp1_1 = df_cnr_hp1
    df_1_ = []
    offset = 0
    for i, chrom in enumerate(chroms):
        df_genes_chrom_1 = df_cnr_hp1_1[df_cnr_hp1_1['chr'] == chrom]
        df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
        offset += regions[i]
        df_1_.append(df_genes_chrom_1)
    df_1 = pd.concat(df_1_)

    df_cnr_hp1 = pd.concat(df_cnr_hp1_)
    df_segs_hp1 = pd.concat(df_segs_hp1_)

    df_cnr_hp2 = pd.concat(df_cnr_hp2_)
    df_segs_hp2 = pd.concat(df_segs_hp2_)

    df_genes_1 = df_genes[df_genes['start'] < df_genes['end']]
    df_genes_2 = df_genes[df_genes['start'] > df_genes['end']]

    last_len = 0
    offset = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        new_len = len(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        if not chrom.startswith('chr'):
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == 'chr' + chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == 'chr' + chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        else:
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        last_len += new_len
        offset += regions[i]

    df_genes_1 = pd.concat(df_genes_1_)
    df_genes_2 = pd.concat(df_genes_2_)

    genename_1 = df_genes_1['gene'].values.tolist()
    genename_2 = df_genes_2['gene'].values.tolist()

    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_cnr_hp1)*args.bin_size, args.bin_size, dtype=int)

    ###########################################################
    row_heights = [420, 40]

    #fig = make_subplots(rows=2, cols=1, shared_yaxes=False, shared_xaxes=True, vertical_spacing=0.02, horizontal_spacing=0.02)
    fig = make_subplots(rows=2, cols=1, shared_yaxes=False, shared_xaxes='columns',  vertical_spacing=0.01, row_heights=row_heights,
                        horizontal_spacing=0.02, specs=[[{"secondary_y":True}], [{"type":"xy"}]]) #https://community.plotly.com/t/can-subplot-support-multiple-y-axes/38891/20
    # #############################################################

    # #############################################################
    custom_text_data_hp1 = df_cnr_hp1.start.values.tolist()
    custom_text_data_hp2 = df_cnr_hp2.start.values.tolist()
    if args.without_phasing:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.coverage.values.tolist(), name='Unphased coverage',
                                   text=custom_text_data_hp1,
                                   yaxis="y", opacity=0.7, color='olive', visibility=True, mul_cols=True, row=1)
    else:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.hp1.values.tolist(), name='HP-1 coverage',
                                   text=custom_text_data_hp1,
                                   yaxis="y", opacity=0.7, color='firebrick', visibility=True, mul_cols=True, row=1)
        add_scatter_trace_coverage(fig, df_1['start'],[-x for x in df_cnr_hp2.hp2.values.tolist()], name='HP-2 coverage',
                                   text=custom_text_data_hp2,
                                   yaxis="y", opacity=0.7, color='steelblue', visibility=True, mul_cols=True, row=1)

    # add_scatter_trace_coverage(fig, df_1['start'], het_snps_bafs_means.vaf.values.tolist(), name='BAF',
    #                                text=het_snps_bafs_means.pos.values.tolist(),
    #                                yaxis="y3", opacity=0.7, color='olive', visibility='legendonly', mul_cols=True, row=2)
        # if not args.unphased_reads_coverage_disable:
        #     add_scatter_trace_coverage(fig, indices, df_unphased.hp3.values.tolist(), name='Unphased', text=None,
        #                                yaxis=None, opacity=0.7, color='olive', visibility='legendonly', mul_cols=True)

    fig.add_trace(go.Scatter(x=df_genes_1['start'], y=[0.5]*len(df_genes_1),#[1, 2, 3, 4, 5] * (len(df_genes_1)//5),
                             mode='markers',
                             text = genename_1,
                             customdata = genestart_1,
                             #hovertext=df_genes['gene'],
                             hovertemplate=
                             '<br><b>Gene</b>: %{text}'+
                             '<br><b>Pos</b>: %{customdata}<br>',
                             marker=dict(
                                 symbol="y-left",
                                 color="#3A6B35",
                                 size=6,
                                 line=dict(width=1, color="#7F7F7F"),
                             ),
                             yaxis="y3",
                             name= 'GeneInfo',
                             showlegend=False,), row=2, col=1,)

    fig.add_hline(y=2,  line_width=1, line_dash="solid", line_color="black", row=2, col=1,)
    fig.add_hline(y=-(args.cut_threshold + 5),  line_width=1, line_dash="solid", line_color="black", row=1, col=1,)

    current = 0
    label_pos = []
    offset_chroms = 0
    offset_chroms_1 = 0
    chroms = get_contigs_list(args.contigs)
    fileDir = os.path.dirname(__file__) #os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])

    fig.add_trace(go.Bar(x=df_genes_1['start'], y=[x for x in df_genes.hp1.values.tolist()], name='HP-1', marker={'color': '#E3B448'}, width=2000000), row = 1, col = 1)
    #fig.add_bar(x=df_genes_1['start'], y=[-(x/2) for x in df_genes.hp2.values.tolist()], name='HP-2', marker={'color': '#3A6B35'}, width=2000000)
    fig.add_trace(go.Bar(x=df_genes_1['start'], y=[-(x) for x in df_genes.hp2.values.tolist()], name='HP-2', marker={'color': '#3A6B35'}, width=2000000), row = 1, col = 1)

    #fig.update_layout(barmode='stack')

    for index, chrom in enumerate(chroms):
        df_gene = df_genes[df_genes['chr'] == chrom]
        genes_starts = df_gene.start.values.tolist()
        genes_ends = df_gene.end.values.tolist()
        genes_name = df_gene.gene.values.tolist()
        genes_coverage_hp1 = df_gene.hp1.values.tolist()
        genes_coverage_hp2 = df_gene.hp2.values.tolist()

        loh_starts = []
        if not loh_regions.empty:
            df_loh_region = loh_regions[loh_regions['chr'] == chrom]
            loh_starts = df_loh_region.start.values.tolist()
            loh_ends = df_loh_region.end.values.tolist()
        if not df_centm.empty:
            df_cent_region = df_centm[df_centm['chr'] == chrom]
            cent_starts = df_cent_region.start.values.tolist()
            cent_ends = df_cent_region.end.values.tolist()

        current += lengths[index]
        label_pos.append(round(offset_chroms+regions[index]//2))

        fig.add_vline(x=offset_chroms,  line_width=1, line_dash="solid", line_color="#D7DBDD")

        offset_chroms_1 += regions[index]

        if index % 2 == 0:
            fig.add_vrect(x0=offset_chroms, x1=offset_chroms_1, fillcolor="#E5E7E9", opacity=0.9, layer="below", line_width=0, )

        # if len(loh_starts) and args.loh_enable:
        #     for i in range(len(df_loh_region.start.values.tolist())):
        #         fig.add_vrect(x0=offset_chroms+loh_starts[i], x1=offset_chroms+loh_ends[i], fillcolor="#2980b9", opacity=0.3, layer="below", line_width=0, row=1, col=1,)
        #
        # if len(cent_starts) and args.loh_enable:
        #     fig.add_vrect(x0=offset_chroms+cent_starts[0], x1=offset_chroms+cent_ends[0], fillcolor="#7e1f14", opacity=0.3, layer="above", line_width=0, row=1, col=1,)

        # if len(genes_starts):
        #     for i in range(len(genes_starts)):
        #         if genes_coverage[i] > 0:
        #             #fig.add_vrect(x0=offset_chroms+genes_starts[i], x1=offset_chroms+genes_ends[i], y0=0, y1=genes_coverage[i], fillcolor="#108c0b",  layer="below", line_width=1, row=1, col=1,)
        #             fig.add_shape(
        #                 plotly.graph_objects.layout.Shape(
        #                     type="rect",
        #                     xref="x",
        #                     yref="y",
        #                     x0=offset_chroms+genes_starts[i],
        #                     y0 = 0,
        #                     x1 = offset_chroms+genes_ends[i],
        #                     y1 = genes_coverage[i],
        #                     name = 'name',
        #                     fillcolor="#0c8702",
        #                     line = dict(color="#c7c7c7", width=0.01, dash='solid'),
        #                     ),
        #                     row = 1,
        #                     col = 1,
        #                     )
        #             fig.add_annotation(x=offset_chroms+genes_starts[i], y=genes_coverage[i],
        #                                text=genes_name[i],
        #                                showarrow=False,
        #                                yshift=10)
        #         else:
        #             fig.add_shape(
        #                 plotly.graph_objects.layout.Shape(
        #                     type="rect",
        #                     xref="x",
        #                     yref="y",
        #                     x0=offset_chroms+genes_starts[i],
        #                     y0 = 0,
        #                     x1 = offset_chroms+genes_ends[i],
        #                     y1 = centers[1],
        #                     name = 'name',
        #                     fillcolor="#d30000",
        #                     line = dict(color="#c7c7c7", width=0.01, dash='solid'),
        #                     ),
        #                     row = 1,
        #                     col = 1,
        #                     )
        #             fig.add_annotation(x=offset_chroms+genes_starts[i], y=centers[1],
        #                                text=genes_name[i],
        #                                showarrow=False,
        #                                yshift=10)

        offset_chroms += regions[index]

    if not args.copynumbers_disable:
        offset_start = 0
        haplotype_1_start_values = []
        haplotype_1_end_values = []

        haplotype_2_start_values = []
        haplotype_2_end_values = []
        for index, chrom in enumerate(chroms):
            df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
            df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
            haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
            haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()
            haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
            haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()
            if chrom == args.contigs.split('-')[0]:
                haplotype_1_start_values.extend(haplotype_1_start_values_copyratios)
                haplotype_1_end_values.extend(haplotype_1_end_values_copyratios)
                haplotype_2_start_values.extend(haplotype_2_start_values_copyratios)
                haplotype_2_end_values.extend(haplotype_2_end_values_copyratios)
            else:
                offset_start += regions[index-1]#lengths[index-1] * args.bin_size
                haplotype_1_start_values.extend([x + offset_start for x in haplotype_1_start_values_copyratios])
                haplotype_1_end_values.extend([x + offset_start for x in haplotype_1_end_values_copyratios])
                haplotype_2_start_values.extend([x + offset_start for x in haplotype_2_start_values_copyratios])
                haplotype_2_end_values.extend([x + offset_start for x in haplotype_2_end_values_copyratios])

        haplotype_1_gaps_values = np.full(len(df_segs_hp1.state.values.tolist()), 'None')
        haplotype_1_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp1.state.values.tolist(), df_segs_hp1.state.values.tolist(), haplotype_1_gaps_values)))
        haplotype_1_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_1_start_values, haplotype_1_end_values, haplotype_1_gaps_values)))

        haplotype_2_gaps_values = np.full(len(df_segs_hp2.state.values.tolist()), 'None')
        haplotype_2_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp2.state.values.tolist(), df_segs_hp2.state.values.tolist(), haplotype_2_gaps_values)))
        haplotype_2_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_2_start_values, haplotype_2_end_values, haplotype_2_gaps_values)))


        if args.without_phasing:
            OFFSET = 0
            colors = ['darkolivegreen']
        else:
            OFFSET = args.cut_threshold/150
            colors = ['firebrick', 'steelblue']
        haplotype_1_copyratios_values = [x if x == 'None' else x  for x in haplotype_1_copyratios_values]
        haplotype_2_copyratios_values = [x if x == 'None' else -x for x in haplotype_2_copyratios_values]
        name = "Copynumbers"

        add_scatter_trace_copyratios(args, fig, colors, name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values, df_segs_hp1, df_segs_hp2, mul_cols=True, row=1, visibility='legendonly')

    centers_rev = [-x for x in centers]
    centers_rev.reverse()
    integer_fractional_means_rev = [x for x in integer_fractional_centers]
    integer_fractional_means_rev.reverse()

    tick_vals = centers + centers_rev
    tickt_ext = integer_fractional_centers + integer_fractional_means_rev
    tickvals = [i for i in range(-1000, 1000, 25)]
    ticktext = [str(abs(i)) for i in range(-1000, 1000, 25)]
    yaxis2_3_range = [-(args.cut_threshold + 5), args.cut_threshold + 5]
    plot_height = 420 + 40 + 15
    legend_y = 1.12

    # #############################################################
    # #############################################################
    #fig.update_yaxes(range=[-1, args.cut_threshold])
    fig.update_layout(
                      yaxis=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis2=dict(range=yaxis2_3_range, showgrid=False,),
                      yaxis3=dict(title="<b>Genes</b>", range=[0, 1], showticklabels = False, showgrid=False, zeroline=True, zerolinewidth=2, zerolinecolor='black'),

                      xaxis=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], zeroline=True, zerolinewidth=1, zerolinecolor='black', showgrid=False,),
                      xaxis2=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size],showgrid=False,),
                      xaxis3=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist())*args.bin_size], showgrid=False, zeroline=True, zerolinewidth=1, zerolinecolor='black'))

    fig.update_layout(
        yaxis=dict(
            title="<b>Coverage depth</b> (per bin)",
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext
        ),
        yaxis2=dict(
            title="<b>Copies</b> (integers/fractions)",
            zeroline=True, zerolinewidth=1, zerolinecolor='black',
            tickmode='array',
            tickvals=tick_vals,
            ticktext=tickt_ext  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
        ),
    )
    # fig.update_xaxes(
    #     tickmode='array',
    #     tickvals=label_pos_chrms,
    #     ticktext=label_chrms  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
    # )

    # fig.update_layout(
    #     hovermode="x unified",
    #     legend_traceorder="normal")

    fig.update_xaxes(
        #xaxis2=dict(
            tickangle=90,
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        #),
        #font=dict(size=18, color="black")
    )
    ax = 20
    ay = -30

    # if args.loh_enable:
    #     add_annotation(fig, 960000000 + 250000000 + 250000000, args.cut_threshold, ax, ay, "LOSS", '#d30000')
    #     add_annotation(fig, 960000000 + 250000000 + 250000000 + 250000000, args.cut_threshold, ax, ay, "GAIN", '#0c8702')

        #fig.add_annotation(text="Het SNPs ratio threshold: " + str(args.hets_ratio) , x = 960000000 + 250000000 + 250000000, y=args.cut_threshold, showarrow=False, row=2, col=1)

    # Update layout
    fig.update_layout(
        template="plotly_white",
        font_family="Times New Roman"
    )

    fig.update_layout(
        title=chrom,
    )
    # Legend
    fig.update_layout(legend=dict(
        orientation='h', xanchor="center", x=0.475, y=legend_y,  # orientation = 'v', xanchor = "center", x = 1.08, y= .5
    ))
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    #fig.update_xaxes(tick0=0.0, rangemode="nonnegative")


    fig.update_layout(legend={'itemsizing': 'constant'})

    fig.update_layout(font_family="Times New Roman")
    if args.without_phasing:
        genome_tile = args.genome_name
    else:
        genome_tile = args.genome_name + "<br>" + "<span style='color:blue'>Ploidy: </span>" + str(args.tumor_ploidy) + "     " + "<span style='color:blue'>Cellular tumor fraction: </span>" + str(args.tumor_purity) + "     " + "<span style='color:blue'>Confidence: </span>" + str(p_value)
    fig.update_layout(
        title={
            'text': genome_tile,
            'y': 0.92,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family="Courier New",
        font_color="dimgray",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="green",
    )
    fig.update_layout(
        width=1380,
        height=plot_height,
       )
    if args.without_phasing:
        path_set = args.genome_name + '_genes_genome_copynumbers'
    else:
        if is_half:
            path_set = 'wgd/'+str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genes_genome_copynumbers'
        else:
            path_set = str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genes_genome_copynumbers'

    if args.pdf_enable:
        print_genome_pdf(fig, path_set, args.out_dir_plots)

    fig.write_html(args.out_dir_plots + '/' + path_set + ".html")

def genes_plots_genome(df_genes, centers, integer_fractional_centers, df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, df_unphased, args, p_value, loh_regions, het_snps_bafs_means, is_half):
    # TODO genome-wide plots
    lengths = []
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)

    df_cnr_hp1_ = []
    df_segs_hp1_ = []
    df_cnr_hp2_ = []
    df_segs_hp2_ = []
    df_genes_1_ = []
    df_genes_2_ = []

    genestart_1 = []
    genestart_2 = []
    gene_cn_hp1 = []
    gene_cn_hp2 = []
    last_len = 0
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        df_segs_hp1_.append(df_segs_hp1[df_segs_hp1['chromosome'] == chrom])
        df_cnr_hp2_.append(df_cnr_hp2[df_cnr_hp2['chr'] == chrom])
        df_segs_hp2_.append(df_segs_hp2[df_segs_hp2['chromosome'] == chrom])

    df_cnr_hp1_1 = df_cnr_hp1
    df_1_ = []
    offset = 0
    for i, chrom in enumerate(chroms):
        df_genes_chrom_1 = df_cnr_hp1_1[df_cnr_hp1_1['chr'] == chrom]
        df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
        offset += regions[i]
        df_1_.append(df_genes_chrom_1)
    df_1 = pd.concat(df_1_)

    df_cnr_hp1 = pd.concat(df_cnr_hp1_)
    df_segs_hp1 = pd.concat(df_segs_hp1_)

    df_cnr_hp2 = pd.concat(df_cnr_hp2_)
    df_segs_hp2 = pd.concat(df_segs_hp2_)

    df_genes_1 = df_genes[df_genes['start'] < df_genes['end']]
    df_genes_2 = df_genes[df_genes['start'] > df_genes['end']]

    last_len = 0
    offset = 0
    #chrs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']
    for i, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chr'] == chrom]
        if df_cnr_hp1_chrom.empty:
            continue
        df_cnr_hp1_.append(df_cnr_hp1_chrom)
        new_len = len(df_cnr_hp1_chrom)
        if not chrom.startswith('chr'):
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == 'chr' + chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())
            gene_cn_hp1.extend(df_genes_chrom_1['copynumber_hp1_state'].values.tolist())
            gene_cn_hp2.extend(df_genes_chrom_1['copynumber_hp2_state'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == 'chr' + chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        else:
            df_genes_chrom_1 = df_genes_1[df_genes_1['chr'] == chrom]
            genestart_1.extend(df_genes_chrom_1['start'].values.tolist())
            gene_cn_hp1.extend(df_genes_chrom_1['copynumber_hp1_state'].values.tolist())
            gene_cn_hp2.extend(df_genes_chrom_1['copynumber_hp2_state'].values.tolist())

            df_genes_chrom_2 = df_genes_2[df_genes_2['chr'] == chrom]
            genestart_2.extend(df_genes_chrom_2['start'].values.tolist())
            if i > 0:
                df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
            df_genes_1_.append(df_genes_chrom_1)

            if i > 0:
                df_genes_chrom_2['start'] = df_genes_chrom_2['start'].apply(lambda x: x + offset)
            df_genes_2_.append(df_genes_chrom_2)
        last_len += new_len
        offset += regions[i]

    df_genes_1 = pd.concat(df_genes_1_)
    df_genes_2 = pd.concat(df_genes_2_)

    genename_1 = df_genes_1['gene'].values.tolist()
    genename_2 = df_genes_2['gene'].values.tolist()

    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_cnr_hp1) * args.bin_size, args.bin_size, dtype=int)

    ###########################################################
    row_heights = [420, 40]

    # fig = make_subplots(rows=2, cols=1, shared_yaxes=False, shared_xaxes=True, vertical_spacing=0.02, horizontal_spacing=0.02)
    fig = make_subplots(rows=2, cols=1, shared_yaxes=False, shared_xaxes='columns', vertical_spacing=0.01,
                        row_heights=row_heights,
                        horizontal_spacing=0.02, specs=[[{"secondary_y": True}], [
            {"type": "xy"}]])  # https://community.plotly.com/t/can-subplot-support-multiple-y-axes/38891/20
    # #############################################################

    # #############################################################
    custom_text_data_hp1 = df_cnr_hp1.start.values.tolist()
    custom_text_data_hp2 = df_cnr_hp2.start.values.tolist()
    if args.without_phasing:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.coverage.values.tolist(), name='Unphased coverage',
                                   text=custom_text_data_hp1,
                                   yaxis="y", opacity=0.7, color='olive', visibility=True, mul_cols=True, row=1)
    else:
        add_scatter_trace_coverage(fig, df_1['start'], df_cnr_hp1.hp1.values.tolist(), name='HP-1 coverage',
                                   text=custom_text_data_hp1,
                                   yaxis="y", opacity=0.7, color='firebrick', visibility=True, mul_cols=True,
                                   row=1)
        add_scatter_trace_coverage(fig, df_1['start'], [x for x in df_cnr_hp2.hp2.values.tolist()], name='HP-2 coverage',
                                   text=custom_text_data_hp2,
                                   yaxis="y", opacity=0.7, color='steelblue', visibility=True, mul_cols=True,
                                   row=1)

    # add_scatter_trace_coverage(fig, df_1['start'], het_snps_bafs_means.vaf.values.tolist(), name='BAF',
    #                                text=het_snps_bafs_means.pos.values.tolist(),
    #                                yaxis="y3", opacity=0.7, color='olive', visibility='legendonly', mul_cols=True, row=2)
    # if not args.unphased_reads_coverage_disable:
    #     add_scatter_trace_coverage(fig, indices, df_unphased.hp3.values.tolist(), name='Unphased', text=None,
    #                                yaxis=None, opacity=0.7, color='olive', visibility='legendonly', mul_cols=True)

    fig.add_trace(
        go.Scatter(x=df_genes_1['start'], y=[0.5] * len(df_genes_1),  # [1, 2, 3, 4, 5] * (len(df_genes_1)//5),
                   mode='markers',
                   text=genename_1,
                   customdata=genestart_1,
                   # hovertext=df_genes['gene'],
                   hovertemplate=
                   '<br><b>Gene</b>: %{text}' +
                   '<br><b>Pos</b>: %{customdata}<br>',
                   marker=dict(
                       symbol="y-left",
                       color="#3A6B35",
                       size=6,
                       line=dict(width=1, color="#7F7F7F"),
                   ),
                   yaxis="y3",
                   name='GeneInfo',
                   showlegend=False, ), row=2, col=1, )

    fig.add_hline(y=2, line_width=1, line_dash="solid", line_color="black", row=2, col=1, )
    fig.add_hline(y=-(args.cut_threshold + 5), line_width=1, line_dash="solid", line_color="black", row=1, col=1, )

    current = 0
    label_pos = []
    offset_chroms = 0
    offset_chroms_1 = 0
    chroms = get_contigs_list(args.contigs)
    fileDir = os.path.dirname(__file__)  # os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])

    fig.add_trace(go.Bar(x=df_genes_1['start'], y=[x for x in df_genes.adjusted_coverage_hp1.values.tolist()], text=genename_1,
                   customdata=gene_cn_hp1,
                   hovertemplate=
                   '<br><b>Gene</b>: %{text}' +
                   '<br><b>CN</b>: %{customdata}<br>', name='HP-1', marker={'color': '#E3B448'}, width=2000000), row=1, col=1)
    fig.add_bar(x=df_genes_1['start'], y=[(x) for x in df_genes.adjusted_coverage_hp2.values.tolist()], text=genename_1,
                   customdata=gene_cn_hp2,
                   hovertemplate=
                   '<br><b>Gene</b>: %{text}' +
                   '<br><b>CN</b>: %{customdata}<br>',name='HP-2', marker={'color': '#3A6B35'}, width=2000000)
    #fig.add_trace(go.Bar(x=df_genes_1['start'], y=[-(x) for x in df_genes.hp2.values.tolist()], name='HP-2',
    #                     marker={'color': '#3A6B35'}, width=2000000), row=1, col=1)

    fig.update_layout(barmode='stack')

    for index, chrom in enumerate(chroms):
        df_gene = df_genes[df_genes['chr'] == chrom]
        genes_starts = df_gene.start.values.tolist()
        genes_ends = df_gene.end.values.tolist()
        genes_name = df_gene.gene.values.tolist()
        genes_coverage_hp1 = df_gene.adjusted_coverage_hp1.values.tolist()
        genes_coverage_hp2 = df_gene.adjusted_coverage_hp2.values.tolist()

        loh_starts = []
        if not loh_regions.empty:
            df_loh_region = loh_regions[loh_regions['chr'] == chrom]
            loh_starts = df_loh_region.start.values.tolist()
            loh_ends = df_loh_region.end.values.tolist()
        if not df_centm.empty:
            df_cent_region = df_centm[df_centm['chr'] == chrom]
            cent_starts = df_cent_region.start.values.tolist()
            cent_ends = df_cent_region.end.values.tolist()

        current += lengths[index]
        label_pos.append(round(offset_chroms + regions[index] // 2))

        fig.add_vline(x=offset_chroms, line_width=1, line_dash="solid", line_color="#D7DBDD")

        offset_chroms_1 += regions[index]

        if index % 2 == 0:
            fig.add_vrect(x0=offset_chroms, x1=offset_chroms_1, fillcolor="#E5E7E9", opacity=0.9, layer="below",
                          line_width=0, )

        # if len(loh_starts) and args.loh_enable:
        #     for i in range(len(df_loh_region.start.values.tolist())):
        #         fig.add_vrect(x0=offset_chroms+loh_starts[i], x1=offset_chroms+loh_ends[i], fillcolor="#2980b9", opacity=0.3, layer="below", line_width=0, row=1, col=1,)
        #
        # if len(cent_starts) and args.loh_enable:
        #     fig.add_vrect(x0=offset_chroms+cent_starts[0], x1=offset_chroms+cent_ends[0], fillcolor="#7e1f14", opacity=0.3, layer="above", line_width=0, row=1, col=1,)

        # if len(genes_starts):
        #     for i in range(len(genes_starts)):
        #         if genes_coverage[i] > 0:
        #             #fig.add_vrect(x0=offset_chroms+genes_starts[i], x1=offset_chroms+genes_ends[i], y0=0, y1=genes_coverage[i], fillcolor="#108c0b",  layer="below", line_width=1, row=1, col=1,)
        #             fig.add_shape(
        #                 plotly.graph_objects.layout.Shape(
        #                     type="rect",
        #                     xref="x",
        #                     yref="y",
        #                     x0=offset_chroms+genes_starts[i],
        #                     y0 = 0,
        #                     x1 = offset_chroms+genes_ends[i],
        #                     y1 = genes_coverage[i],
        #                     name = 'name',
        #                     fillcolor="#0c8702",
        #                     line = dict(color="#c7c7c7", width=0.01, dash='solid'),
        #                     ),
        #                     row = 1,
        #                     col = 1,
        #                     )
        #             fig.add_annotation(x=offset_chroms+genes_starts[i], y=genes_coverage[i],
        #                                text=genes_name[i],
        #                                showarrow=False,
        #                                yshift=10)
        #         else:
        #             fig.add_shape(
        #                 plotly.graph_objects.layout.Shape(
        #                     type="rect",
        #                     xref="x",
        #                     yref="y",
        #                     x0=offset_chroms+genes_starts[i],
        #                     y0 = 0,
        #                     x1 = offset_chroms+genes_ends[i],
        #                     y1 = centers[1],
        #                     name = 'name',
        #                     fillcolor="#d30000",
        #                     line = dict(color="#c7c7c7", width=0.01, dash='solid'),
        #                     ),
        #                     row = 1,
        #                     col = 1,
        #                     )
        #             fig.add_annotation(x=offset_chroms+genes_starts[i], y=centers[1],
        #                                text=genes_name[i],
        #                                showarrow=False,
        #                                yshift=10)

        offset_chroms += regions[index]

    if not args.copynumbers_disable:
        offset_start = 0
        haplotype_1_start_values = []
        haplotype_1_end_values = []

        haplotype_2_start_values = []
        haplotype_2_end_values = []
        for index, chrom in enumerate(chroms):
            df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
            df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
            haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
            haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()
            haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
            haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()
            if chrom == args.contigs.split('-')[0]:
                haplotype_1_start_values.extend(haplotype_1_start_values_copyratios)
                haplotype_1_end_values.extend(haplotype_1_end_values_copyratios)
                haplotype_2_start_values.extend(haplotype_2_start_values_copyratios)
                haplotype_2_end_values.extend(haplotype_2_end_values_copyratios)
            else:
                offset_start += regions[index - 1]  # lengths[index-1] * args.bin_size
                haplotype_1_start_values.extend([x + offset_start for x in haplotype_1_start_values_copyratios])
                haplotype_1_end_values.extend([x + offset_start for x in haplotype_1_end_values_copyratios])
                haplotype_2_start_values.extend([x + offset_start for x in haplotype_2_start_values_copyratios])
                haplotype_2_end_values.extend([x + offset_start for x in haplotype_2_end_values_copyratios])

        haplotype_1_gaps_values = np.full(len(df_segs_hp1.state.values.tolist()), 'None')
        haplotype_1_copyratios_values = list(itertools.chain.from_iterable(
            zip(df_segs_hp1.state.values.tolist(), df_segs_hp1.state.values.tolist(), haplotype_1_gaps_values)))
        haplotype_1_copyratios_positions = list(itertools.chain.from_iterable(
            zip(haplotype_1_start_values, haplotype_1_end_values, haplotype_1_gaps_values)))

        haplotype_2_gaps_values = np.full(len(df_segs_hp2.state.values.tolist()), 'None')
        haplotype_2_copyratios_values = list(itertools.chain.from_iterable(
            zip(df_segs_hp2.state.values.tolist(), df_segs_hp2.state.values.tolist(), haplotype_2_gaps_values)))
        haplotype_2_copyratios_positions = list(itertools.chain.from_iterable(
            zip(haplotype_2_start_values, haplotype_2_end_values, haplotype_2_gaps_values)))

        if args.without_phasing:
            OFFSET = 0
            colors = ['darkolivegreen']
        else:
            OFFSET = args.cut_threshold / 150
            colors = ['firebrick', 'steelblue']
        haplotype_1_copyratios_values = [x if x == 'None' else x for x in haplotype_1_copyratios_values]
        haplotype_2_copyratios_values = [x if x == 'None' else x for x in haplotype_2_copyratios_values]
        name = "Copynumbers"

        add_scatter_trace_copyratios(args, fig, colors, name, haplotype_1_copyratios_positions,
                                     haplotype_2_copyratios_positions, haplotype_1_copyratios_values,
                                     haplotype_2_copyratios_values, df_segs_hp1, df_segs_hp2, mul_cols=True, row=1,
                                     visibility='legendonly')

    centers_rev = [-x for x in centers]
    centers_rev.reverse()
    integer_fractional_means_rev = [x for x in integer_fractional_centers]
    integer_fractional_means_rev.reverse()

    tick_vals = centers
    tickt_ext = integer_fractional_centers
    tickvals = [i for i in range(0, 1000, 25)]
    ticktext = [str(abs(i)) for i in range(0, 1000, 25)]
    yaxis2_3_range = [0, args.cut_threshold + 5]
    plot_height = 420 + 40 + 15
    legend_y = 1.12

    # #############################################################
    # #############################################################
    # fig.update_yaxes(range=[-1, args.cut_threshold])
    fig.update_layout(
        yaxis=dict(range=yaxis2_3_range, showgrid=False, ),
        yaxis2=dict(range=yaxis2_3_range, showgrid=False, ),
        yaxis3=dict(title="<b>Genes</b>", range=[0, 1], showticklabels=False, showgrid=False, zeroline=True,
                    zerolinewidth=2, zerolinecolor='black'),

        xaxis=dict(tick0=0.0, rangemode="nonnegative", range=[0, len(df_cnr_hp1.start.values.tolist()) * args.bin_size],
                   zeroline=True, zerolinewidth=1, zerolinecolor='black', showgrid=False, ),
        xaxis2=dict(tick0=0.0, rangemode="nonnegative",
                    range=[0, len(df_cnr_hp1.start.values.tolist()) * args.bin_size], showgrid=False, ),
        xaxis3=dict(tick0=0.0, rangemode="nonnegative",
                    range=[0, len(df_cnr_hp1.start.values.tolist()) * args.bin_size], showgrid=False, zeroline=True,
                    zerolinewidth=1, zerolinecolor='black'))

    fig.update_layout(
        yaxis=dict(
            title="<b>Coverage depth</b> (per bin)",
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext
        ),
        yaxis2=dict(
            title="<b>Copies</b> (integers/fractions)",
            zeroline=True, zerolinewidth=1, zerolinecolor='black',
            tickmode='array',
            tickvals=tick_vals,
            ticktext=tickt_ext  # ['loss'] + [str(i) + '_copy' for i in range(1,len(centers))]
        ),
    )

    fig.update_xaxes(
        # xaxis2=dict(
        tickangle=90,
        tickmode='array',  # change 1
        tickvals=label_pos,  # change 2
        ticktext=chroms,  # change 3
        # ),
        # font=dict(size=18, color="black")
    )
    ax = 20
    ay = -30

    # if args.loh_enable:
    #     add_annotation(fig, 960000000 + 250000000 + 250000000, args.cut_threshold, ax, ay, "LOSS", '#d30000')
    #     add_annotation(fig, 960000000 + 250000000 + 250000000 + 250000000, args.cut_threshold, ax, ay, "GAIN", '#0c8702')

    # fig.add_annotation(text="Het SNPs ratio threshold: " + str(args.hets_ratio) , x = 960000000 + 250000000 + 250000000, y=args.cut_threshold, showarrow=False, row=2, col=1)

    # Update layout
    fig.update_layout(
        template="plotly_white",
        font_family="Times New Roman"
    )

    fig.update_layout(
        title=chrom,
    )
    # Legend
    fig.update_layout(legend=dict(
        orientation='h', xanchor="center", x=0.475, y=legend_y,
        # orientation = 'v', xanchor = "center", x = 1.08, y= .5
    ))
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    # fig.update_xaxes(tick0=0.0, rangemode="nonnegative")

    fig.update_layout(legend={'itemsizing': 'constant'})

    fig.update_layout(font_family="Times New Roman")
    if args.without_phasing:
        genome_tile = args.genome_name
    else:
        genome_tile = args.genome_name + "<br>" + "<span style='color:blue'>Ploidy: </span>" + str(
            args.tumor_ploidy) + "     " + "<span style='color:blue'>Cellular tumor fraction: </span>" + str(
            args.tumor_purity) + "     " + "<span style='color:blue'>Confidence: </span>" + str(p_value)
    fig.update_layout(
        title={
            'text': genome_tile,
            'y': 0.92,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family="Courier New",
        font_color="dimgray",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="green",
    )
    fig.update_layout(
        width=1380,
        height=plot_height,
    )
    if args.without_phasing:
        path_set = args.genome_name + '_genes_genome'
    else:
        if is_half:
            path_set = 'wgd/'+str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genes_genome'
        else:
            path_set = str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_genes_genome'

    if args.pdf_enable:
        print_genome_pdf(fig, path_set, args.out_dir_plots)

    fig.write_html(args.out_dir_plots + '/' + path_set + ".html")

def heatmap_copy_number_plots_genome(df_genes, centers, integer_fractional_centers, df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, df_unphased, args, p_value, loh_regions, het_snps_bafs_means, is_half):
    # TODO genome-wide plots
    lengths = []
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)

    df_cnr_hp1_ = []
    df_segs_hp1_ = []
    df_cnr_hp2_ = []
    df_segs_hp2_ = []
    df_genes_1_ = []
    df_genes_2_ = []

    genestart_1 = []
    genestart_2 = []
    last_len = 0
    df_cnr_hp1, df_cnr_hp2 = bins_with_copynumber_states(df_segs_hp1, df_segs_hp2, df_cnr_hp1, df_cnr_hp2, args, centers, integer_fractional_centers)

    for i, chrom in enumerate(chroms):
        df_cnr_hp1_.append(df_cnr_hp1[df_cnr_hp1['chr'] == chrom])
        df_segs_hp1_.append(df_segs_hp1[df_segs_hp1['chromosome'] == chrom])
        df_cnr_hp2_.append(df_cnr_hp2[df_cnr_hp2['chr'] == chrom])
        df_segs_hp2_.append(df_segs_hp2[df_segs_hp2['chromosome'] == chrom])

    df_cnr_hp1_1 = df_cnr_hp1
    df_1_ = []
    offset = 0
    for i, chrom in enumerate(chroms):
        df_genes_chrom_1 = df_cnr_hp1_1[df_cnr_hp1_1['chr'] == chrom]
        df_genes_chrom_1['start'] = df_genes_chrom_1['start'].apply(lambda x: x + offset)
        offset += regions[i]
        df_1_.append(df_genes_chrom_1)
    df_1 = pd.concat(df_1_)

    df_cnr_hp1 = pd.concat(df_cnr_hp1_)
    df_segs_hp1 = pd.concat(df_segs_hp1_)

    df_cnr_hp2 = pd.concat(df_cnr_hp2_)
    df_segs_hp2 = pd.concat(df_segs_hp2_)

    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    fig = make_subplots(rows=2, cols=1, vertical_spacing=0.1, shared_xaxes=True, horizontal_spacing=0.02, row_heights=[0.30, 0.30])

    cmap = [
        (r, c)
        for r, c in zip(
            np.repeat(np.linspace(0, 1, len(px.colors.sequential.RdBu) + 1), 2)[1:],
            np.repeat(px.colors.sequential.RdBu, 2),
        )
    ]
    HP1 = go.Heatmap(
        z=[df_cnr_hp1.hp1.values.tolist()],
        x=[i for i in range(0, len(df_cnr_hp1.hp1.values.tolist()), 1)],#df_1['start'],
        y=['HP-1'],
        colorscale=cmap,
        # colorscale=[
        #     [0, '#1c46e3'],
        #     [1, '#e3281c']],
        # colorbar_tickvals=np.arange(0, 8, 1),
        # colorbar_tickmode='array',
        #colorscale='Bluered',
        colorbar=dict(y=0.8, len=.5),
    )
    HP2 = go.Heatmap(
        z=[df_cnr_hp2.hp2.values.tolist()],
        x=[i for i in range(0, len(df_cnr_hp2.hp2.values.tolist()), 1)],#df_1['start'],
        y=['HP-2'],
        colorscale=cmap,
        #colorscale='Bluered',
        colorbar=dict(y=.2, len=.5),
        # colorscale=[[0, "rgb(166,206,227)"],
        #             [0.25, "rgb(31,120,180)"],
        #             [0.50, "rgb(178,223,138)"],
        #             [0.75, "rgb(51,160,44)"],
        #             [1, "rgb(227,26,28)"]],
        # colorscale=[
        #     [0, '#1c46e3'],
        #     [1, '#e3281c']],
        # colorbar_tickvals=np.arange(0, 8, 1), #np.arange(0, integer_fractional_centers[-1], 1),
        # colorbar_tickmode='array',
    )



    fig.add_trace(HP1, row=1, col=1)
    fig.add_trace(HP2, row=2, col=1)

    current = 0
    label_pos = []
    for index, chrom in enumerate(chroms):
        if index % 2 == 0:
            # fig.add_shape(type="rect",
            #               xref="x", yref="y",
            #               x0=current, y0=0,
            #               x1=current+lengths[index], y1=2,
            #
            #               line=dict(
            #                   color="#E5E7E9",
            #                   width=0.5,
            #               ),
            #               fillcolor="#E5E7E9")

            fig.add_vrect(x0=current, x1=current+lengths[index], fillcolor="#E5E7E9", opacity=0.8, layer="below", line_width=5, )

        label_pos.append(round(current + lengths[index] // 2))
        current += lengths[index]

    fig.update_xaxes(
        # xaxis2=dict(
        tickangle=90,
        tickmode='array',  # change 1
        tickvals=label_pos,  # change 2
        ticktext=chroms,  # change 3
        # ),
        # font=dict(size=18, color="black")
    )

    plot_height = 420 + 40 + 15
    legend_y = 1.12
    # Update layout
    fig.update_layout(
        template="plotly_white",
        font_family="Times New Roman"
    )

    fig.update_layout(
        title=chrom,
    )
    # Legend
    fig.update_layout(legend=dict(
        orientation='h', xanchor="center", x=0.475, y=legend_y,
        # orientation = 'v', xanchor = "center", x = 1.08, y= .5
    ))
    fig.update_layout(margin=dict(l=15, r=15, b=15, pad=10))
    # fig.update_xaxes(tick0=0.0, rangemode="nonnegative")

    fig.update_layout(legend={'itemsizing': 'constant'})

    fig.update_layout(font_family="Times New Roman")
    if args.without_phasing:
        genome_tile = args.genome_name
    else:
        genome_tile = args.genome_name + "<br>" + "<span style='color:blue'>Ploidy: </span>" + str(
            args.tumor_ploidy) + "     " + "<span style='color:blue'>Cellular tumor fraction: </span>" + str(
            args.tumor_purity) + "     " + "<span style='color:blue'>Confidence: </span>" + str(p_value)
    fig.update_layout(
        title={
            'text': genome_tile,
            'y': 0.92,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family="Courier New",
        font_color="dimgray",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="green",
    )
    fig.update_layout(
        width=1380,
        height=plot_height,
    )

    if args.without_phasing:
        path_set = args.genome_name + '_heatmap_genome'
    else:
        if is_half:
            path_set = 'wgd/'+str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_heatmap_genome'
        else:
            path_set = str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '/' + args.genome_name + '_' + str(args.tumor_ploidy) + '_' + str(args.tumor_purity) + '_' + str(p_value) + '_heatmap_genome'

    if args.pdf_enable:
        print_genome_pdf(fig, path_set, args.out_dir_plots)

    fig.write_html(args.out_dir_plots + '/' + path_set + ".html")

def plots_add_markers_lines(fig):
    # style all the traces
    fig.update_traces(
        # hoverinfo="name+x+text+y",
        # line={"width": 0.5},
        marker={"size": 2},
        mode="markers",
        showlegend=True
    )

def add_scatter_trace_coverage(fig, x, y, name, text, yaxis, opacity, color, visibility=True, mul_cols=False, row=2, baf=None):
    if any(value < 0 for value in y):
        xy = [round(-y2, 2) for y2 in y]
    else:
        xy = [round(y2, 2) for y2 in y]
    if baf:
        ht = '<br><b>Position</b>: %{text}' + '<br><b>BAF</b>: %{customdata}<br>'
    else:
        ht = '<br><b>Position</b>: %{text}' + '<br><b>Coverage</b>: %{customdata}<br>'
    if mul_cols:
        fig.add_trace(go.Scatter(
            customdata=xy,
            # legendgroup="group1",  # this can be any string, not just "group"
            # legendgrouptitle_text="Coverage",
            hovertemplate= ht,
            x=x,
            y=y,
            name=name,
            text=text,
            yaxis=yaxis,
            opacity=opacity,
            marker_color=color,
            visible=visibility,
            marker={"size": 2},
            mode="markers",
        ), row = row, col = 1)
    else:
        fig.add_trace(go.Scatter(
            # legendgroup="group1",  # this can be any string, not just "group"
            # legendgrouptitle_text="Coverage",
            customdata=xy,
            hovertemplate=ht,
            x=x,
            y=y,
            name=name,
            text=text,
            # yaxis="y5",
            opacity=opacity,
            marker_color=color,
            visible=visibility,
            marker={"size": 2},
            mode="markers",
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
        visible='legendonly',
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
        visible='legendonly',
        opacity=0.5,
        marker_color=['darkgreen', 'limegreen', 'white']*len(phaseblocks_positions),
        showlegend=True,
        marker_symbol='diamond-wide',
        hoverinfo = "x+name+y+text",
        #legendgroup="group2",
    ))

def add_scatter_trace_copyratios(args, fig, colors, name, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values, df_segs_hp1, df_segs_hp2, mul_cols, row=2, centers=[], visibility=True):
    haplotype_1_copyratios_positions_start = df_segs_hp1.start.values.tolist()
    haplotype_1_copyratios_positions_end = df_segs_hp1.end.values.tolist()
    haplotype_2_copyratios_positions_start = df_segs_hp2.start.values.tolist()
    haplotype_2_copyratios_positions_end = df_segs_hp2.end.values.tolist()

    hp1_chr_info = [j for i in  [[x,y, None] for (x,y) in zip(haplotype_1_copyratios_positions_start, haplotype_1_copyratios_positions_end)] for j in i]
    hp2_chr_info = [j for i in  [[x,y, None] for (x,y) in zip(haplotype_2_copyratios_positions_start, haplotype_2_copyratios_positions_end)] for j in i]

    if 'confidence_value' in df_segs_hp1.columns:
        haplotype_1_copyratios_positions_p_values = [j for i in  [[x,y, None] for (x,y) in zip(df_segs_hp1.confidence_value.values.tolist(), df_segs_hp1.confidence_value.values.tolist())] for j in i]
        haplotype_2_copyratios_positions_p_values = [j for i in  [[x,y, None] for (x,y) in zip(df_segs_hp2.confidence_value.values.tolist(), df_segs_hp2.confidence_value.values.tolist())] for j in i]
        haplotype_1_copyratios_positions_cn_values = [subclonal_values_adjusted(x, centers) for x in haplotype_1_copyratios_values] #[j for i in  [[x,y, None] for (x,y) in zip(df_segs_hp1.state.values.tolist(), df_segs_hp1.state.values.tolist())] for j in i]
        haplotype_2_copyratios_positions_cn_values = [subclonal_values_adjusted(x, centers) for x in haplotype_2_copyratios_values] #[j for i in  [[x,y, None] for (x,y) in zip(df_segs_hp2.state.values.tolist(), df_segs_hp2.state.values.tolist())] for j in i]
        df_cn_sub_info_1 = pd.DataFrame(list(zip(haplotype_1_copyratios_positions_p_values, haplotype_1_copyratios_positions_cn_values)), columns=['p_1', 'cn_1'])
        df_cn_sub_info_2 = pd.DataFrame(list(zip(haplotype_2_copyratios_positions_p_values, haplotype_2_copyratios_positions_cn_values)), columns=['p_2', 'cn_2'])
        ht = '<br><b>Position</b>: %{text}' + '<br><b>CN</b>: %{customdata[1]}' + '<br><b>confidence</b>: %{customdata[0]}<br>'
    elif 'p_value' in df_segs_hp1.columns:
        haplotype_1_copyratios_positions_p_values = [j for i in  [[x,y, None] for (x,y) in zip(df_segs_hp1.p_value.values.tolist(), df_segs_hp1.p_value.values.tolist())] for j in i]
        haplotype_2_copyratios_positions_p_values = [j for i in  [[x,y, None] for (x,y) in zip(df_segs_hp2.p_value.values.tolist(), df_segs_hp2.p_value.values.tolist())] for j in i]
        haplotype_1_copyratios_positions_cn_values = [subclonal_values_adjusted(x, centers) for x in haplotype_1_copyratios_values] #[j for i in  [[subclonal_values_adjusted(x, centers), subclonal_values_adjusted(y, centers), None] for (x,y) in zip(haplotype_1_copyratios_values, haplotype_1_copyratios_values)] for j in i]
        haplotype_2_copyratios_positions_cn_values = [subclonal_values_adjusted(x, centers) for x in haplotype_2_copyratios_values] #[j for i in  [[subclonal_values_adjusted(x, centers), subclonal_values_adjusted(y, centers), None] for (x,y) in zip(haplotype_2_copyratios_values, haplotype_2_copyratios_values)] for j in i]
        df_cn_sub_info_1 = pd.DataFrame(list(zip(haplotype_1_copyratios_positions_p_values, haplotype_1_copyratios_positions_cn_values)), columns=['p_1', 'cn_1'])
        df_cn_sub_info_2 = pd.DataFrame(list(zip(haplotype_2_copyratios_positions_p_values, haplotype_2_copyratios_positions_cn_values)), columns=['p_2', 'cn_2'])
        ht = '<br><b>Position</b>: %{text}' + '<br><b>CN</b>: %{customdata[1]}' + '<br><b>confidence</b>: %{customdata[0]}<br>'
    else:
        haplotype_1_copyratios_positions_p_values = []
        haplotype_2_copyratios_positions_p_values = []
        df_cn_sub_info_1 = []
        df_cn_sub_info_2 = []
        ht = '<br><b>Position</b>: %{text}' + '<br><b>CN</b>: %{y}<br>'

    if mul_cols:
        if args.without_phasing:
            name_leg = 'Unphased'
        else:
            name_leg = 'HP-1 CNA'
        fig.add_trace(go.Scatter(
            #legendgroup="group2",
            #legendgrouptitle_text=name,
            hovertemplate = ht,
            customdata = df_cn_sub_info_1,
            x=haplotype_1_copyratios_positions,
            y= haplotype_1_copyratios_values,
            name=name_leg,
            text=hp1_chr_info,
            yaxis="y3",
            line = dict(shape = 'spline', color = colors[0], width= 5, dash = 'solid'),
            mode='lines',
            #marker={"size": 5},
            opacity=0.9,
            visible=visibility,
            marker_color=[colors[0], colors[0], 'white']*len(haplotype_1_copyratios_positions),
            showlegend=True,
            #marker_symbol='diamond-wide',
            #hoverinfo = "x+name+y+text",
            #legendgroup="group2",
            #legendgrouptitle_text="Phaseblocks",
        ), row = row, col = 1, secondary_y=True)
        if args.without_phasing == False:
            fig.add_trace(go.Scatter(
                #legendgroup="group2",
                customdata=df_cn_sub_info_2,
                hovertemplate= ht,
                x=haplotype_2_copyratios_positions,
                y=haplotype_2_copyratios_values,
                name='HP-2 CNA',
                text=hp2_chr_info,
                yaxis="y3",
                line = dict(shape = 'spline', color = colors[1], width= 5, dash = 'solid'),
                mode='lines',
                visible=visibility,
                #marker={"size": 5},
                opacity=0.9,
                marker_color=[colors[1], colors[1], 'white']*len(haplotype_2_copyratios_positions),
                showlegend=True,
                #marker_symbol='diamond-wide',
                #hoverinfo = "x+name+y+text",
                #legendgroup="group2",
            ), row = row, col = 1, secondary_y=True)
    else:
        if args.without_phasing:
            name_leg = 'Unphased'
        else:
            name_leg = 'HP-1 CNA'
        fig.add_trace(go.Scatter(
            #legendgroup="group2",
            #legendgrouptitle_text=name,
            hovertemplate=
            '<br><b>Position</b>: %{text}' +
            '<br><b>CN</b>: %{y}<br>',
            x=haplotype_1_copyratios_positions,
            y= haplotype_1_copyratios_values,
            name=name_leg,
            text=hp1_chr_info,
            yaxis="y2",
            visible=visibility,
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
        if args.without_phasing == False:
            fig.add_trace(go.Scatter(
                #legendgroup="group2",
                hovertemplate=
                '<br><b>Position</b>: %{text}' +
                '<br><b>CN</b>: %{y}<br>',
                x=haplotype_2_copyratios_positions,
                y=haplotype_2_copyratios_values,
                name='HP-2 CNA',
                text=hp2_chr_info,
                yaxis="y2",
                visible=visibility,
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

def plots_layout_settings(fig, chrom, args, limit_x, limit_y):
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
            #title="<b>Coverage depth</b> (per bin)",
            #titlefont={"color": "dimgray"},
            title=dict(text="<b>Coverage depth</b> (per bin)", font={"color": "dimgray"}),
            type="linear",
            showline=True,
            zeroline=True,
        ),
        yaxis2=dict(
            linecolor="dimgray",
            range=[1, args.cut_threshold + 5],
            side="right",
            tickfont={"color": "dimgray"},
            tickmode="auto",
            ticks="outside",
            #title="<b>Copies</b> (integers/fractions)",
            #titlefont={"color": "dimgray"},
            title=dict(text="<b>Copies</b> (integers/fractions)", font={"color": "dimgray"}),
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
            'text': chrom + ' - ' +args.genome_name,
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

# def flat_with_diff(A, atol=3):
#     runs = np.split(A, np.where(np.abs(np.diff(A)) >= atol)[0] + 1)
#     return [list(x) for x in runs if len(x) > 1]

# def change_point_detection(hp_data, start, ends, ref_start_values_1, args, chrom, html_graphs, hp, color):
#     import ruptures as rpt
#     fig = go.Figure()
# 
#     #starts = [i for i in range(0, len(data), 50000)]
#     add_scatter_trace_coverage(fig, start, hp_data, name='HP-'+str(hp), text=None, yaxis=None, opacity=0.7, color=color)
# 
#     ############################################
#     # temp_start = []
#     # extra_cpd_points = []
#     # start_pos = False
#     # for k in range(5, len(start)-5):
#     #     if start_pos == False and statistics.mean(hp_data[k-5:k]) > abs(statistics.mean(hp_data[k-4:k-1]) - hp_data[k])  > statistics.mean(hp_data[k-5:k]):
#     #         temp_start.append(k)
#     #         start_pos = True
#     #     elif start_pos and statistics.mean(hp_data[k-5:k]) > abs(statistics.mean(hp_data[k+1:k+5]) - hp_data[k])  > statistics.mean(hp_data[k:k+5]) :
#     #         temp_start.append(k)
#     #         extra_cpd_points.append(temp_start[0])
#     #         extra_cpd_points.append(temp_start[-1])
#     #         temp_start = []
#     #         start_pos = False
# 
#     ##############################################
#     cpd_peaks_points = []
#     # temp_start = []
#     # for k in range(1, len(start) - 1):
#     #     if abs(hp_data[k-1] - hp_data[k]) > 30:
#     #         temp_start.append(k)
#     #         for l in range(k+1,len(start) - 1):
#     #             if abs(hp_data[k] - hp_data[l]) > 30:
#     #                 continue
#     #             else:
#     #                 if l > k:
#     #                     temp_start.append(l)
#     #                 break
#     # temp_start = sorted(list(set(temp_start)))
#     # temp_start_groups = flat_with_diff(temp_start)
#     #
#     # for i, val in enumerate(temp_start_groups):
#     #     if len(val) > 5:
#     #         cpd_peaks_points.append(val[0])
#     #         cpd_peaks_points.append(val[-1])
# 
#     ############################################
#     zeros_values = []
#     for i in range(len(hp_data)):
#         if hp_data[i] < 1:
#             zeros_values.append(i)
# 
#     from itertools import groupby, count
#     groups = groupby(zeros_values, key=lambda item, c=count(): item - next(c))
#     tmp = [list(g) for k, g in groups]
# 
#     cpd_zeros_points = []
#     for i, val in enumerate(tmp):
#         if len(val) == 1:
#             hp_data[i] = statistics.mean(hp_data[val[0]:val[0]+2])
#         else:
#             cpd_zeros_points.append(val[0])
#             cpd_zeros_points.append(val[-1])
#     ####################################################
#     data = np.array(hp_data, dtype='int') #numpy.clip(data, a_min=0, a_max=1000)
#     algo = rpt.Pelt(model="rbf", jump=25).fit(data)
#     result = algo.predict(pen=10)
#     change_points = [i for i in result if i < len(data)]
#     ############################################
#     ov_indices = []
#     for i, val1 in enumerate(change_points):
#         for j, val2 in enumerate(tmp):
#             if len(val2) > 2:
#                 if val2[-1]>= val1 >= val2[0]:
#                     ov_indices.append(i)
#     if ov_indices:
#         for index in sorted(list(set(ov_indices)), reverse=True):
#             del change_points[index]
#     ############################################
#     change_points = sorted(list(set(change_points + cpd_zeros_points + cpd_peaks_points)))
#     ############################################
#     ############################################
#     from src.smoothing import smooth_triangle
#     hp_data_new = []
# 
#     if len(hp_data[0:change_points[0]]) > 46:
#         hp_data_new.extend(smooth_triangle(hp_data[0:change_points[0]], 15))
#     else:
#         hp_data_new.extend(hp_data[0:change_points[0]])
# 
#     for p in range(1, len(change_points)):
#         if change_points[p] - change_points[p - 1] > 46:
#             hp_data_new.extend(smooth_triangle(hp_data[change_points[p - 1]:change_points[p]], 15))
#         else:
#             hp_data_new.extend(hp_data[change_points[p - 1]:change_points[p]])
# 
#     if len(hp_data[change_points[-1]:]) > 46:
#         hp_data_new.extend(smooth_triangle(hp_data[change_points[-1]:], 15))
#     else:
#         hp_data_new.extend(hp_data[change_points[-1]:])
#     ############################################
#     add_scatter_trace_coverage(fig, start, hp_data_new, name='HP-' + str(hp), text=None, yaxis=None, opacity=0.7, color='grey')
#     ############################################
# 
#     #model = CUSUM(k=1., h=2., burnin=50, mu=0., sigma=1.)
#     #model = EWMA(r=0.15, L=2.4, burnin=50, mu=0., sigma=1.)
#     #model = TwoSample(statistic="Lepage", threshold=3.1)
# 
#     # model.process(data)
#     # change_points = model.changepoints
# 
#     if args.variable_size_bins:
#         start_pos = 0
#         snps_haplotype_mean = []
#         snps_haplotype_len = []
#         snps_haplotype_pos = []
#         snps_haplotype_pos.append(0)
#         for index, point in enumerate(change_points):
#             i = start[point-1]
# 
#             sub_list = hp_data[start_pos:i]
#             if sub_list:
#                 snps_haplotype_mean.append(statistics.median(sub_list))
#                 snps_haplotype_len.append(len(sub_list))
#                 snps_haplotype_pos.append(i)
#             start_pos = point + 1
# 
#         for i, point in enumerate(snps_haplotype_pos):
#             fig.add_vline(x=point, line_width=1, line_dash="dash", line_color=color)
#     else:
#         for i, point in enumerate(change_points):
#             fig.add_vline(x=point*args.bin_size, line_width=1, line_dash="dash", line_color=color)
# 
#     #plots_add_markers_lines(fig)
#     plots_layout_settings(fig, chrom, args, ends[-1:][0], args.cut_threshold)
# 
#     print_chromosome_html(fig, chrom + '_hp_'  + str(hp), html_graphs, args.out_dir_plots)
#     html_graphs.write("  <object data=\"" + chrom + '_hp_'  + str(hp)  + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")
# 



# def add_histo_clusters_plot(depth_values_hp1, depth_values_hp2, labels, means, covar, args, chrom, html_graphs):
#     from plotly.subplots import make_subplots
#     import plotly.graph_objects as go
#     import math
#     # fig = go.Figure()
# 
#     depth_values_hp1 = np.array(depth_values_hp1, dtype='int')
#     depth_values_hp2 = np.array(depth_values_hp2, dtype='int')
# 
#     depth_values_hp1 = depth_values_hp1.reshape(-1, 1)
#     depth_values_hp2 = depth_values_hp2.reshape(-1, 1)
#     Y = np.concatenate([depth_values_hp1, depth_values_hp2])
# 
#     fig = make_subplots(rows=1, cols=2, shared_yaxes=True, column_widths=[0.8, 0.2], vertical_spacing=0.01,
#                         horizontal_spacing=0.01)
# 
#     y = np.concatenate(list(Y)).ravel().tolist()
#     lst = [i for i in range(0, (len(y) // 2) * args.bin_size, args.bin_size)]
#     x = lst + lst
# 
#     cdict = {0: '#1f77b4',  # muted blue
#              1: '#ff7f0e',  # safety orange
#              2: '#2ca02c',  # cooked asparagus green
#              3: '#d62728',  # brick red
#              4: '#9467bd',  # muted purple
#              5: '#8c564b',  # chestnut brown
#              6: '#e377c2',  # raspberry yogurt pink
#              7: '#7f7f7f',  # middle gray
#              8: '#bcbd22',  # curry yellow-green
#              9: '#17becf'  # blue-teal
#              }
#     ldict = {0: 'Cluster_1', 1: 'Cluster_2', 2: 'Cluster_3', 3: 'Cluster_4', 4: 'Cluster_5', \
#              5: 'Cluster_6', 6: 'Cluster_7', 7: 'Cluster_8', 8: 'Cluster_9', 9: 'Cluster_10'}
# 
#     for g in np.unique(labels):
#         ix = [index for index, i in enumerate(x) if labels[index] == g]
#         xn = [x[i] for i in ix]
#         yn = [y[i] for i in ix]
#         fig.add_trace(go.Scatter(x=xn, y=yn, mode='markers', marker=dict(color=cdict[g]), name=ldict[g], opacity=0.7, ),
#                       row=1, col=1)
# 
#     fig.update_traces(
#         # hoverinfo="name+x+text+y",
#         # line={"width": 0.5},
#         marker={"size": 2},
#         mode="markers",
#         showlegend=True
#     )
# 
#     for i in range(len(means)):
#         fig.add_hline(y=means[i], line_width=2,
#                       line=dict(dash='solid'), line_color=cdict[i], annotation_position="top right", annotation_text="mean_" + str(i + 1))
# 
#         fig.add_hline(y=means[i] + covar[i], line_width=1,
#                       line=dict(dash='dash'), line_color=cdict[i], annotation_position="top left", annotation_text="stdev_" + str(i + 1))
#         fig.add_hline(y=means[i] - covar[i], line_width=1,
#                       line=dict(dash='dash'), line_color=cdict[i], annotation_position="top left", annotation_text="stdev_" + str(i + 1))
#     fig.update_yaxes(range=[0, 160])
# 
#     # fig.add_trace(go.Histogram(y=y, orientation='h', nbinsy=5000,), row=1, col=2)
#     fig.add_trace(
#         go.Histogram(y=np.concatenate(list(Y)).ravel().tolist(), orientation='h', marker_color='gray', name='Histogram',
#                      nbinsy=8000), row=1, col=2)
#     fig.update_layout(xaxis2=dict(range=[0, 200]))
#     # fig.add_trace(go.Histogram(y=np.concatenate(list(Y[0:len(Y)//2-1])).ravel().tolist(), name='HP-1', orientation='h', nbinsy=8000, marker_color='#6A5ACD'), row=1, col=2)
#     # fig.add_trace(go.Histogram(y=np.concatenate(list(Y[len(Y)//2:len(Y)])).ravel().tolist(), name='HP-2', orientation='h', nbinsy=8000, marker_color='#2E8B57'), row=1, col=2)
# 
#     # Overlay both histograms
#     # fig.update_layout(barmode='overlay')
#     # Reduce opacity to see both histograms
#     # fig.update_traces(opacity=0.75)
#     fig.update_layout(legend={'itemsizing': 'constant'})
#     # Update layout
#     fig.update_layout(
#         template="plotly_white",
#         font_family="Times New Roman"
#     )
#     fig.update_layout(
#         title=chrom,
#     )
#     fig.update_layout(legend=dict(
#         orientation = 'h', xanchor = "center", x = 0.5, y= 1.2 #orientation = 'v', xanchor = "center", x = 1.08, y= .5
#     ))
#     fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
#     fig.update_layout(
#         width=680,
#         height=400,
#        )
#     fig.update_layout(
#         title={
#             'text': chrom + ' - ' + args.genome_name,
#             'y': 0.96,
#             'x': 0.5,
#             'xanchor': 'center',
#             'yanchor': 'top'},
# 
#         font_family="Courier New",
#         font_color="dimgray",
#         title_font_family="Times New Roman",
#         title_font_color="red",
#         legend_title_font_color="green",
#     )
# 
#     print_chromosome_html(fig, chrom + '_clusters_'  + str(len(means)), html_graphs, args.out_dir_plots)
#     html_graphs.write("  <object data=\"" + chrom + '_clusters_'  + str(len(means))  + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")
# 
#     #fig.write_html("coverage_plots/" + chrom + "_cluster_" + str(len(means)) + ".html")

def plot_ploidy_purity_p_values(args, ploidy_values, purity_values, p_values):
    def plot_heatmap_with_nan_fill(p_values, ploidy, purity, title="Heatmap Plot"):
        """
        Plots a heatmap with NaN values in p_values filled with the lowest p_value.

        Args:
            p_values: List or NumPy array of p-values.
            ploidy: List or NumPy array of ploidy values.
            purity: List or NumPy array of purity values.
            title: Title of the plot.
        """

        # Create a DataFrame
        df = pd.DataFrame({'p_value': p_values, 'ploidy': ploidy, 'purity': purity})

        # Find the lowest non-NaN p-value
        min_p_value = df['p_value'].min()

        # Fill NaN values with the lowest p-value
        df['p_value_filled'] = df['p_value'].fillna(min_p_value)

        # Create the heatmap
        '''
        fig = go.Figure(data=go.Heatmap(
            z=df['p_value_filled'],
            x=df['ploidy'],
            y=df['purity'],
            colorscale='Viridis',  # Choose a colorscale
            zmin=df['p_value_filled'].min(), # set zmin and zmax for consistent color scaling.
            zmax=df['p_value_filled'].max()
        ))
        '''
        pivot_table = df.pivot_table(values='p_value', index='purity', columns='ploidy')
        # Find the lowest and highest p-values
        lowest_p_value = df['p_value'].min()
        highest_p_value = df['p_value'].max()

        # Fill NaN values with the lowest p-value
        pivot_table = pivot_table.fillna(lowest_p_value)

        # Text annotations, initialize with empty strings
        # annotations = np.empty_like(data, dtype=str)

        # Specify the row and column index for annotation
        # row_index = 2
        # col_index = 3

        # Annotation text
        annotations = [
            dict(x=10, y=10, bgcolor='white', text="Top Left", showarrow=False),  # Top-left cell
            dict(x=20, y=20, bgcolor='white', text="Middle Right", showarrow=True),  # Middle-right cell
            dict(x=39, y=30, bgcolor='white', text="Bottom Center", showarrow=False)  # Bottom-center cell
        ]
        # Add the annotation text
        # annotations[row_index, col_index] = "Important Value"

        fig = go.Figure(data=go.Heatmap(
            z=pivot_table,
            x=pivot_table.columns,
            y=pivot_table.index,
            xgap=0.41, ygap=0.41,
            # text=annotations,
            # texttemplate="%{text}", # Display the text
            # textfont={"size":10},
            colorbar=dict(title='p-value'),
            colorscale='BrBG'  # 'Cividis_r'#'Viridis'  # Choose a suitable colorscale
        ))

        fig.update_layout(
            width=800,
            height=600,
            # title=title,
            xaxis_title="Ploidy",
            yaxis_title="Purity",
            xaxis=dict(type='category'),  # Treat ploidy as categorical
            yaxis=dict(type='category'),  # Treat purity as categorical
        )
        fig.update_layout(
            {
                "plot_bgcolor": "rgba(0, 0, 0, 0)",  # make the background transparent
                # 'paper_bgcolor': 'rgba(0, 0, 0, 0)'
                "margin": {"l": 0, "r": 0, "t": 0, "b": 0}
            }
        )
        fig.update_layout(
            template="plotly_white",
            font_family="Times New Roman"
        )
        #TODO Add annotations
        # fig.update_layout(
        #     annotations=annotations
        # )

        # fig.add_shape(
        # type="rect",
        # xref=pivot_table.columns,
        # yref=pivot_table.index,
        # x0=4,  # Start x coordinate (column index 1)
        # y0=10,  # Start y coordinate (row index 1)
        # x1=2,  # End x coordinate (column index 2)
        # y1=20,  # End y coordinate (row index 2)
        # line=dict(color="red", width=3)  # Customize the rectangle's appearance
        # )

        #TODO
        # fig.add_shape(type="rect",
        #               x0=-0.5, y0=1.5, x1=3.5, y1=5.5,
        #               line=dict(color="blue", width=4),
        #               )

        # fig.update_layout(
        # bargap=1,        # Gap between bars of adjacent categories (as a fraction of the bar width)
        # bargroupgap=1,   # Gap between bars within the same category (as a fraction of the bar width)
        # barmode='group'    # Display bars in groups
        # )
        path_set = args.genome_name + '_heatmap_ploidy_purity.html'
        if args.pdf_enable:
            print_genome_pdf(fig, path_set, args.out_dir_plots)

        fig.write_html(args.out_dir_plots + '/' + path_set)

    """
    if (purity_values[0] < float(args.purity_range.split('-')[0]) or purity_values[0] > float(args.purity_range.split('-')[1])) or (
            ploidy_values[0] < float(args.ploidy_range.split('-')[0]) or ploidy_values[0] > float(args.ploidy_range.split('-')[1])):
        del purity_values[0]
        del ploidy_values[0]
        del p_values[0]
    """

    purity = [round(num, 4) for num in purity_values]
    ploidy = [round(num, 4) for num in ploidy_values]

    plot_heatmap_with_nan_fill(p_values, ploidy, purity)

def plot_cytobands_genome_wide(fig, cytobands_df, title="Cytobands Whole Genome"):
    """
    Plots cytobands with labels and hover information.

    Args:
        cytobands_df: DataFrame with columns 'chrom', 'start', 'end', 'name', 'gieStain'.
        title: Title of the plot.

    Returns:
        A Plotly Figure object.
    """

    cytobands_df = cytobands_df.copy()
    cytobands_df['chrom'] = pd.Categorical(cytobands_df['chrom'], categories=[f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"], ordered=True)
    cytobands_df = cytobands_df.sort_values(by=['chrom', 'start'])

    x_positions = []
    shapes = []
    hover_texts = []
    annotations = []  # Added annotations list
    current_x = 0

    for _, row in cytobands_df.iterrows():
        chrom, start, end, name, gie_stain = row['chrom'], row['start'], row['end'], row['name'], row['gieStain']
        width = end - start

        color = {
            'gneg': 'white',
            'gpos25': 'lightgray',
            'gpos50': 'gray',
            'gpos75': 'darkgray',
            'gpos100': 'black',
            'acen': 'red',
            'stalk': 'darkred',
            'gvar': 'lightgreen',
        }.get(gie_stain, 'white')

        shapes.append(
            dict(
                xref="x4",
                yref="y4",
                type="rect",
                x0=current_x,
                y0=0,
                x1=current_x + width,
                y1=1,
                fillcolor=color,
                line=dict(color='lightgray', width=0.5),
            )
        )
        hover_texts.append(f"<br>{name}<br>{chrom}:{start/1000000:.2f}M-{end/1000000:.2f}M")
        x_positions.append(current_x + width / 2)

        # Add annotation for each cytoband
        # annotations.append(
        #     dict(
        #         x=current_x + width / 2,
        #         y=0.5,
        #         xref="x",
        #         yref="y",
        #         text=f"<br>{name}<br>{chrom}:{start/1000000:.2f}M-{end/1000000:.2f}M", #f"{chrom}<br>{name}<br>{start:,}:{end:,}",
        #         showarrow=False,
        #         font=dict(size=4),
        #         textangle=-90,
        #         xanchor="center",
        #         yanchor="middle",
        #     )
        # )

        current_x += width

    # fig.update_layout(
    #     shapes=shapes,
    #     #hovermode="closest",
    #     #annotations=annotations,  # Add annotations to the layout
    # )
    for shape in shapes:
        fig.add_shape(shape, row=3, col=1)

    fig.add_trace(go.Scatter(x=x_positions, y=[0.5]*len(x_positions), yaxis="y4", showlegend=False, marker=dict(symbol="y-left", color="#5a5c5a", size=6, line=dict(width=1, color="#7F7F7F"),), hovertemplate = '%{text}', name= 'CytoBandsInfo', text=hover_texts, mode = 'markers'), row=3, col=1)

    return fig


