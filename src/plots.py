import numpy
import plotly.graph_objects as go
import plotly
import pandas as pd
import numpy as np
import csv
import os
import math
import statistics
import itertools
import logging
from phasing_correction import phaseblock_flipping
from smoothing import smoothing
from vcf_processing import get_snps_frquncies_coverage, vcf_parse_to_csv_for_het_phased_snps_phasesets
from utils import csv_df_chromosomes_sorter_snps, get_breakpoints, flatten

chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']#, 'chrX', 'chrY']
def copy_number_log2_ratios_plots_chromosomes(df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, arguments):
    filename = f"{os.path.join(arguments['out_dir_plots'], 'COPY_NUMBER_LOG2_RATIOS.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")

    for index, chrom in enumerate(chroms):
        if not chrom == 'chrX' and not chrom == 'chrY':
            logging.info('Plots generation for ' + chrom)
            fig = go.Figure()

            df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chromosome'] == chrom]
            df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]

            df_cnr_hp2_chrom = df_cnr_hp2[df_cnr_hp2['chromosome'] == chrom]
            df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]

            haplotype_1_values_cnr = df_cnr_hp1_chrom.log2.values.tolist()
            haplotype_2_values_cnr = df_cnr_hp2_chrom.log2.values.tolist()
            haplotype_1_start_values_cnr = df_cnr_hp1_chrom.start.values.tolist()

            haplotype_1_values_copyratios = df_segs_hp1_chrom.data.state.values.tolist()
            haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
            haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()

            haplotype_2_values_copyratios = df_segs_hp2_chrom.data.state.values.tolist()
            haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
            haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()

            if arguments['copyratios_enable']:
                logging.info('copyratios plots module')
                OFFSET=.25
                haplotype_1_values_cnr = list(np.asarray(haplotype_1_values_cnr) + OFFSET)
                haplotype_2_values_cnr = list(np.asarray(haplotype_2_values_cnr) - OFFSET)
                add_scatter_trace_coverage(fig, haplotype_1_start_values_cnr, haplotype_1_values_cnr, name='HP-1', text=None, yaxis=None, opacity=0.7, color='firebrick')
                add_scatter_trace_coverage(fig, haplotype_1_start_values_cnr, haplotype_2_values_cnr, name='HP-2', text=None, yaxis=None, opacity=0.7, color='steelblue')
                plots_add_markers_lines(fig)

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

                add_scatter_trace_copyratios(fig, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values)

                plots_layout_settings(fig, chrom, arguments, haplotype_1_end_values_copyratios[-1:][0])
                #fig.update_yaxes(range = [-10,10])
                #fig.update_yaxes(title='copy ratio (log2)')

            if arguments['pdf_enable']:
                print_chromosome_pdf(fig, chrom, arguments['out_dir_plots'])

            print_chromosome_html(fig, chrom+'_cnr', html_graphs, arguments['out_dir_plots'])

    html_graphs.write("</body></html>")

def coverage_plots_chromosomes(df, df_phasesets, arguments):
    filename = f"{os.path.join(arguments['out_dir_plots'], 'COVERAGE.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")
    haplotype_1_values_updated = [[]]
    haplotype_2_values_updated = [[]]
    hunphased_updated = [[]]

    for index, chrom in enumerate(chroms):
        if not chrom == 'chrX' and not chrom == 'chrY':
            logging.info('Plots generation for ' + chrom)
            fig = go.Figure()

            df_chrom = df[df['chr'] == chrom]
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

            if arguments['phaseblock_flipping_enable']:
                logging.info('phaseblock flipping module')
                is_simple_correction = False
                if len(ref_start_values_phasesets) < 5:
                    is_simple_correction = True
                haplotype_1_values, haplotype_2_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = \
                phaseblock_flipping(is_simple_correction, haplotype_1_values, haplotype_2_values, ref_start_values, ref_end_values, \
                        haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

            if arguments['smoothing_enable']:
                logging.info('smoothing module')
                unphased_reads_values, haplotype_1_values, haplotype_2_values = smoothing(unphased_reads_values, haplotype_1_values, haplotype_2_values, conv_window_size=15)

            add_scatter_trace_coverage(fig, ref_start_values, haplotype_1_values, name='HP-1', text=None, yaxis=None, opacity=0.7, color='firebrick')
            add_scatter_trace_coverage(fig, ref_start_values, haplotype_2_values, name='HP-2', text=None, yaxis=None, opacity=0.7, color='steelblue')

            if arguments['unphased_reads_coverage_enable']:
                add_scatter_trace_coverage(fig, ref_start_values, unphased_reads_values, name='Unphased reads', text=None, yaxis=None, opacity=0.7, color='olive')

            if arguments['het_phased_snps_freq_enable']:
                logging.info('hetrozygous phased snps frequencies coverage module')
                output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(arguments['phased_vcf_snps_freqs'])
                csv_df_snps = csv_df_chromosomes_sorter_snps(output_phasesets_file_path)
                haplotype_1_snps_freqs, haplotype_2_snps_freqs = get_snps_frquncies_coverage(csv_df_snps, chrom, ref_start_values, arguments['bin_size'])
                add_scatter_trace_coverage(fig, ref_start_values, haplotype_1_snps_freqs, name='HP-1 SNPs Freqs', text=None, yaxis=None,
                                           opacity=0.7, color=None)
                add_scatter_trace_coverage(fig, ref_start_values, haplotype_2_snps_freqs, name='HP-2 SNPs Freqs', text=None, yaxis=None,
                                           opacity=0.7, color=None)
            plots_add_markers_lines(fig)

            if arguments['phaseblocks_enable']:
                logging.info('phaseblocks plots module')
                gaps_values = np.full(len(haplotype_1_values_phasesets), 'None')
                haplotype_1_phaseblocks_values = list(itertools.chain.from_iterable(zip(haplotype_1_values_phasesets, haplotype_1_values_phasesets, gaps_values)))
                haplotype_2_phaseblocks_values = list(itertools.chain.from_iterable(zip(haplotype_2_values_phasesets, haplotype_2_values_phasesets, gaps_values)))
                phaseblocks_positions = list(itertools.chain.from_iterable(zip(ref_start_values_phasesets, ref_end_values_phasesets, gaps_values)))

                add_scatter_trace_phaseblocks(fig, phaseblocks_positions, haplotype_1_phaseblocks_values, haplotype_2_phaseblocks_values)

            if arguments['breakpoints_enable']:
                breakpoints = get_breakpoints(chrom, arguments['breakpoints_file'])
                add_scatter_trace_breakpoints(fig, breakpoints)

            plots_layout_settings(fig, chrom, arguments, ref_end_values[-1:][0])

            if arguments['pdf_enable']:
                print_chromosome_pdf(fig, chrom, arguments['out_dir_plots'])

            print_chromosome_html(fig, chrom, html_graphs, arguments['out_dir_plots'])

            haplotype_1_values_updated.append(haplotype_1_values)
            haplotype_2_values_updated.append(haplotype_2_values)
            hunphased_updated.append(unphased_reads_values)

    #plot_bins_histograms(flatten(haplotype_1_values_updated), flatten(haplotype_2_values_updated), ref_start_values, arguments)
    #plot_bins_ratios(flatten(haplotype_1_values_updated), flatten(haplotype_2_values_updated), arguments)
    html_graphs.write("</body></html>")

    return haplotype_1_values_updated, haplotype_2_values_updated, hunphased_updated

def print_genome_pdf(fig, genome, coverage_plots_path):
    filename = f"{os.path.join(coverage_plots_path, genome + '.pdf')}"
    plotly.io.write_image(fig,  filename, format='pdf')
def print_chromosome_pdf(fig, chrom, coverage_plots_path):
    filename = f"{os.path.join(coverage_plots_path, chrom + '.pdf')}"
    plotly.io.write_image(fig,  filename, format='pdf')

def print_chromosome_html(fig, chrom, html_graphs, coverage_plots_path):
    fname = f"{os.path.join(coverage_plots_path, chrom + '.html')}"
    plotly.offline.plot(fig, filename=fname,auto_open=False)
    html_graphs.write("  <object data=\""+chrom+'.html'+"\" width=\"700\" height=\"420\"></object>"+"\n")

def plots_genome_coverage(df_hp1, df_hp2, df_unphased, arguments):
    # TODO genome-wide plots

    lengths = []
    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_hp1[df_hp1['chr'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_hp1)*50000, 50000, dtype=int)

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

    for index, chrom in enumerate(chroms):
        current += lengths[index]
        label_pos.append(round(current - (lengths[index] / 2))*50000)
        fig.add_vline(x=current*50000, y0=-10, y1=150, line_width=1, line_dash="dashdot", line_color="green")

    fig.update_layout(
        xaxis=dict(
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        ),
        font=dict(size=18, color="black"))

    plots_layout_settings(fig, 'Genome', arguments, indices[-1:][0])
    #fig.update_yaxes(range=[-10, 10])
    #fig.update_yaxes(title='copy ratio (log2)')
    fig.update_layout(width=1080, height=400,)

    if arguments['pdf_enable']:
        print_genome_pdf(fig, arguments['genome_name'], arguments['out_dir_plots'])

    fig.write_html(arguments['out_dir_plots'] +'/'+ arguments['genome_name'] + "_genome_coverage.html")

def plots_genome(df_cnr_hp1, df_segs_hp1, df_cnr_hp2, df_segs_hp2, arguments):
    # TODO genome-wide plots

    lengths = []
    for index, chrom in enumerate(chroms):
        df_cnr_hp1_chrom = df_cnr_hp1[df_cnr_hp1['chromosome'] == chrom]
        lengths.append(len(df_cnr_hp1_chrom))

    indices = np.arange(0, len(df_cnr_hp1)*50000, 50000, dtype=int)

    fig = go.Figure()
    add_scatter_trace_coverage(fig, indices, df_cnr_hp1.log2.values.tolist(), name='HP-1', text=None,
                               yaxis=None, opacity=0.7, color='firebrick')
    add_scatter_trace_coverage(fig, indices, df_cnr_hp2.log2.values.tolist(), name='HP-2', text=None,
                               yaxis=None, opacity=0.7, color='steelblue')
    plots_add_markers_lines(fig)

    current = 0
    label_pos = []

    for index, chrom in enumerate(chroms):
        current += lengths[index]
        label_pos.append(round(current - (lengths[index] / 2))*50000)
        fig.add_vline(x=current*50000, y0=-10, y1=150, line_width=1, line_dash="dashdot", line_color="green")

    fig.update_layout(
        xaxis=dict(
            tickmode='array',  # change 1
            tickvals=label_pos,  # change 2
            ticktext=chroms,  # change 3
        ),
        font=dict(size=18, color="black"))

    if arguments['copyratios_enable']:
        offset = 0
        haplotype_1_start_values = []
        haplotype_1_end_values = []

        haplotype_2_start_values = []
        haplotype_2_end_values = []
        for index, chrom in enumerate(chroms):
            if not chrom == 'chrX' and not chrom == 'chrY':
                df_segs_hp1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
                df_segs_hp2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]
                haplotype_1_start_values_copyratios = df_segs_hp1_chrom.start.values.tolist()
                haplotype_1_end_values_copyratios = df_segs_hp1_chrom.end.values.tolist()
                haplotype_2_start_values_copyratios = df_segs_hp2_chrom.start.values.tolist()
                haplotype_2_end_values_copyratios = df_segs_hp2_chrom.end.values.tolist()
                if chrom == 'chr1':
                    haplotype_1_start_values.extend(haplotype_1_start_values_copyratios)
                    haplotype_1_end_values.extend(haplotype_1_end_values_copyratios)
                    haplotype_2_start_values.extend(haplotype_2_start_values_copyratios)
                    haplotype_2_end_values.extend(haplotype_2_end_values_copyratios)
                else:
                    offset += lengths[index-1] * 50000
                    haplotype_1_start_values.extend([x + offset for x in haplotype_1_start_values_copyratios])
                    haplotype_1_end_values.extend([x + offset for x in haplotype_1_end_values_copyratios])
                    haplotype_2_start_values.extend([x + offset for x in haplotype_2_start_values_copyratios])
                    haplotype_2_end_values.extend([x + offset for x in haplotype_2_end_values_copyratios])

        haplotype_1_gaps_values = np.full(len(df_segs_hp1.data.state.values.tolist()), 'None')
        haplotype_1_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp1.data.state.values.tolist(), df_segs_hp1.data.state.values.tolist(), haplotype_1_gaps_values)))
        haplotype_1_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_1_start_values, haplotype_1_end_values, haplotype_1_gaps_values)))

        haplotype_2_gaps_values = np.full(len(df_segs_hp2.data.state.values.tolist()), 'None')
        haplotype_2_copyratios_values = list(itertools.chain.from_iterable(zip(df_segs_hp2.data.state.values.tolist(), df_segs_hp2.data.state.values.tolist(), haplotype_2_gaps_values)))
        haplotype_2_copyratios_positions = list(itertools.chain.from_iterable(zip(haplotype_2_start_values, haplotype_2_end_values, haplotype_2_gaps_values)))

        add_scatter_trace_copyratios(fig, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values)

    plots_layout_settings(fig, 'Genome', arguments, indices[-1:][0])
    fig.update_yaxes(range=[-1, arguments['cut_threshold']])
    #fig.update_yaxes(title='copy ratio (log2)')
    fig.update_layout(width=1080, height=400,)

    if arguments['pdf_enable']:
        print_genome_pdf(fig, arguments['genome_name']+'_copynumbers', arguments['out_dir_plots'])

    fig.write_html(arguments['out_dir_plots'] +'/'+ arguments['genome_name'] + "_genome_copynumber.html")

def plots_add_markers_lines(fig):
    # style all the traces
    fig.update_traces(
        # hoverinfo="name+x+text+y",
        # line={"width": 0.5},
        marker={"size": 2},
        mode="markers",
        showlegend=True
    )

def add_scatter_trace_coverage(fig, x, y, name, text, yaxis, opacity, color):

    fig.add_trace(go.Scatter(
        x=x,
        y=y,
        name=name,
        #text=text,
        yaxis="y5",
        opacity=opacity,
        marker_color=color,
    ))

def add_scatter_trace_phaseblocks(fig, phaseblocks_positions, haplotype_1_phaseblocks_values, haplotype_2_phaseblocks_values):
    fig.add_trace(go.Scatter(
        x=phaseblocks_positions,
        y=haplotype_1_phaseblocks_values,
        name="HP-1 Phaseblocks",
        text=phaseblocks_positions,
        yaxis="y5",
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
        x=phaseblocks_positions,
        y=haplotype_2_phaseblocks_values,
        name="HP-2 Phaseblocks",
        text=phaseblocks_positions,
        yaxis="y5",
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

def add_scatter_trace_copyratios(fig, haplotype_1_copyratios_positions, haplotype_2_copyratios_positions, haplotype_1_copyratios_values, haplotype_2_copyratios_values):
    fig.add_trace(go.Scatter(
        x=haplotype_1_copyratios_positions,
        y= haplotype_1_copyratios_values,
        name="HP-1 Copy-numbers",
        text=haplotype_1_copyratios_positions,
        yaxis="y5",
        line = dict(shape = 'spline', color = 'gray', width= 5, dash = 'solid'),
        mode='lines',
        #marker={"size": 5},
        opacity=0.5,
        marker_color=['dimgray', 'darkgray', 'white']*len(haplotype_1_copyratios_positions),
        showlegend=True,
        #marker_symbol='diamond-wide',
        #hoverinfo = "x+name+y+text",
        #legendgroup="group2",
        #legendgrouptitle_text="Phaseblocks",
    ))

    fig.add_trace(go.Scatter(
        x=haplotype_2_copyratios_positions,
        y=haplotype_2_copyratios_values,
        name="HP-2 Copy-numbers",
        text=haplotype_2_copyratios_positions,
        yaxis="y5",
        line = dict(shape = 'spline', color = 'green', width= 5, dash = 'solid'),
        mode='lines',
        #marker={"size": 5},
        opacity=0.5,
        marker_color=['darkgreen', 'limegreen', 'white']*len(haplotype_2_copyratios_positions),
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

def plots_layout_settings(fig, chrom, arguments, limit):
    # Update axes
    fig.update_layout(
        xaxis=dict(
            #autorange=True,
            type="linear",
            showline=True,
            zeroline=True,
            linecolor = "dimgray",
            range=[0,limit*1.01]
        ),
        yaxis5=dict(
            linecolor="dimgray",
            range=[1, arguments['cut_threshold']+5],
            side="left",
            tickfont={"color": "dimgray"},
            tickmode="auto",
            ticks="outside",
            title="coverage (meandepth)",
            titlefont={"color": "dimgray"},
            type="linear",
            showline=True,
            zeroline=True,
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
        orientation = 'h', xanchor = "center", x = 0.5, y= 1.2
    ))
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    fig.update_xaxes(tick0=0.0, rangemode="nonnegative")

    fig.update_layout(legend={'itemsizing': 'constant'})

    fig.update_layout(font_family= "Times New Roman")

    fig.update_layout(
        title={
            'text': arguments['genome_name'] + ' - ' +chrom,
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

def plot_bins_histograms(hp1,hp2, start, arguments):
    import plotly.express as px
    from sklearn.cluster import KMeans
    import plotly.figure_factory as ff

    hp1=numpy.clip(hp1, a_min=1, a_max=300)
    hp2=numpy.clip(hp2, a_min=1, a_max=300)

    df = pd.DataFrame(dict(
        Haplotypes=np.concatenate((["HP-1"] * len(hp1), ["HP-2"] * len(hp2))),
        coverage=np.concatenate((hp1, hp2))
    ))
    # df = pd.DataFrame(dict(
    #     Haplotypes=(["HP-1/HP-2"] * len(hp1)),
    #     coverage=(hp1/(hp1+hp2)),
    # ))


    #X=list(zip(hp1,start))



    df = pd.DataFrame(dict(
        Haplotypes=(["HP-1 & HP-2"] * (len(hp1)*2)),
        coverage=np.concatenate((hp1, hp2)),
    ))

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
    kn = KneeLocator(x=K, y=Sum_of_squared_distances, curve='convex', direction='decreasing')
    print(kn.knee)

    km = KMeans(n_clusters=kn.knee+1, n_init=25, max_iter = 600, random_state=0)
    km = km.fit(new)
    labels = list(km.labels_)
    u_labels = np.unique(labels)
    centers = list(np.concatenate(km.cluster_centers_))

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

    #fig = px.histogram(df, x="coverage", color="Haplotypes", barmode="overlay", marginal="violin")#, log_y=True)


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

    fig.write_html(arguments['out_dir_plots'] +'/'+ arguments['genome_name'] + "_genome_histogram.html")
