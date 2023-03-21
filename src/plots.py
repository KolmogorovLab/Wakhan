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
from vcf_processing import get_snps_frquncies_coverage
from utils import csv_df_chromosomes_sorter_snps

chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']#, 'chrX', 'chrY']
def coverage_plots_chromosomes(df, df_phasesets, arguments, output_phasesets_file_path):
    filename = f"{os.path.join(arguments['out_dir_plots'], 'DASHBOARD.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")

    for index, chrom in enumerate(chroms):
        logging.info('Plots generation for ' + chrom)
        fig = go.Figure()

        df_chrom = df[df['chr'] == chrom]
        df_chrom_phasesets = df_phasesets[df_phasesets['chr'] == chrom]

        unphased_reads_values = df_chrom.hp3.values.tolist()
        haplotype_1_values = df_chrom.hp1.values.tolist()
        haplotype_2_values = df_chrom.hp2.values.tolist()
        ref_start_values = df_chrom.start.values.tolist()
        #ref_end_values = df_chrom.end.values.tolist()

        haplotype_1_values_phasesets = df_chrom_phasesets.hp1.values.tolist()
        haplotype_2_values_phasesets = df_chrom_phasesets.hp2.values.tolist()
        ref_start_values_phasesets = df_chrom_phasesets.start.values.tolist()
        ref_end_values_phasesets = df_chrom_phasesets.end.values.tolist()

        if arguments['phaseblock_flipping_enable']:
            logging.info('phaseblock flipping module')
            haplotype_1_values, haplotype_2_values = \
            phaseblock_flipping(haplotype_1_values, haplotype_2_values, ref_start_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

        if arguments['smoothing_enable']:
            logging.info('smoothing module')
            unphased_reads_values, haplotype_1_values, haplotype_2_values = smoothing(unphased_reads_values, haplotype_1_values, haplotype_2_values, conv_window_size=15)

        add_scatter_trace_coverage(fig, ref_start_values, unphased_reads_values, name='Unphased reads', text=None, yaxis=None, opacity=None, color='lightgreen')
        add_scatter_trace_coverage(fig, ref_start_values, haplotype_1_values, name='HP-1', text=None, yaxis=None, opacity=None, color='lightpink')
        add_scatter_trace_coverage(fig, ref_start_values, haplotype_2_values, name='HP-2', text=None, yaxis=None, opacity=None, color='cornflowerblue')

        if arguments['het_phased_snps_freq_enable']:
            logging.info('hetrozygous phased snps frequencies coverage module')
            csv_df_snps = csv_df_chromosomes_sorter_snps(output_phasesets_file_path)
            haplotype_1_snps_freqs, haplotype_2_snps_freqs = get_snps_frquncies_coverage(csv_df_snps, chrom, ref_start_values, arguments['bin_size'])
            add_scatter_trace_coverage(fig, ref_start_values, haplotype_1_snps_freqs, name='HP-1 SNPs Freqs', text=None, yaxis=None,
                                       opacity=None, color=None)
            add_scatter_trace_coverage(fig, ref_start_values, haplotype_2_snps_freqs, name='HP-2 SNPs Freqs', text=None, yaxis=None,
                                       opacity=None, color=None)

        plots_add_markers_lines(fig)

        if arguments['phaseblocks_enable']:
            logging.info('phaseblocks plots module')
            gaps_values = np.full(len(haplotype_1_values_phasesets), 'None')
            haplotype_1_phaseblocks_values = list(itertools.chain.from_iterable(zip(haplotype_1_values_phasesets, haplotype_1_values_phasesets, gaps_values)))
            haplotype_2_phaseblocks_values = list(itertools.chain.from_iterable(zip(haplotype_2_values_phasesets, haplotype_2_values_phasesets, gaps_values)))
            phaseblocks_positions = list(itertools.chain.from_iterable(zip(ref_start_values_phasesets, ref_end_values_phasesets, gaps_values)))

            add_scatter_trace_phaseblocks(fig, phaseblocks_positions, haplotype_1_phaseblocks_values, haplotype_2_phaseblocks_values)

        plots_layout_settings(fig, chrom, arguments)

        if arguments['pdf_enable']:
            print_chromosome_pdf(fig, chrom, arguments['out_dir_plots'])

        print_chromosome_html(fig, chrom, html_graphs, arguments['out_dir_plots'])

    html_graphs.write("</body></html>")
def print_chromosome_pdf(fig, chrom, coverage_plots_path):
    filename = f"{os.path.join(coverage_plots_path, chrom + '.pdf')}"
    plotly.io.write_image(fig,  filename, format='pdf')

def print_chromosome_html(fig, chrom, html_graphs, coverage_plots_path):
    fname = f"{os.path.join(coverage_plots_path, chrom + '.html')}"
    plotly.offline.plot(fig, filename=fname,auto_open=False)
    html_graphs.write("  <object data=\""+chrom+'.html'+"\" width=\"700\" height=\"420\"></object>"+"\n")
def coverage_plots_genome(df):
    # TODO genome-wide plots
    print('')

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
        opacity=0.7,
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
def plots_layout_settings(fig, chrom, arguments):
    # Update axes
    fig.update_layout(
        xaxis=dict(
            autorange=True,
            type="linear",
            showline=True,
            zeroline=True,
            linecolor = "dimgray",
        ),
        yaxis5=dict(
            linecolor="dimgray",
            range=[0, arguments['cut_threshold']],
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



