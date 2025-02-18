import plotly.graph_objects as go
import plotly
import numpy as np
import os
import itertools

import logging
logger = logging.getLogger()

from src.hapcorrect.src.utils import get_contigs_list, get_chromosomes_regions

def plot_coverage_data(html_graphs, args, chrom, ref_start_values, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, sufix):
    fig = go.Figure()

    add_scatter_trace_coverage(fig, ref_start_values, haplotype_1_values, name='HP-1', text=ref_start_values, yaxis=None,
                               opacity=0.7, color='firebrick')
    add_scatter_trace_coverage(fig, ref_start_values, haplotype_2_values, name='HP-2', text=ref_start_values, yaxis=None,
                               opacity=0.7, color='steelblue')

    if not args.unphased_reads_coverage_disable:
        add_scatter_trace_coverage(fig, ref_start_values, unphased_reads_values, name='Unphased', text=ref_start_values,
                                   yaxis=None, opacity=0.7, color='olive')
    plots_add_markers_lines(fig)

    if args.phaseblocks_enable:
        gaps_values = np.full(len(haplotype_1_values_phasesets), 'None')
        haplotype_1_phaseblocks_values = list(
            itertools.chain.from_iterable(zip(haplotype_1_values_phasesets, haplotype_1_values_phasesets, gaps_values)))
        haplotype_2_phaseblocks_values = list(
            itertools.chain.from_iterable(zip(haplotype_2_values_phasesets, haplotype_2_values_phasesets, gaps_values)))
        phaseblocks_positions = list(
            itertools.chain.from_iterable(zip(ref_start_values_phasesets, ref_end_values_phasesets, gaps_values)))

        add_scatter_trace_phaseblocks(fig, phaseblocks_positions, haplotype_1_phaseblocks_values,
                                      haplotype_2_phaseblocks_values)
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)
    chrom_index = chroms.index(chrom)
    plots_layout_settings(fig, chrom, args, regions[chrom_index], args.cut_threshold)

    if args.pdf_enable:
        print_chromosome_pdf(fig, chrom, args.out_dir_plots+'/phasing_output')

    print_chromosome_html(fig, chrom + '_' + sufix, html_graphs, args.out_dir_plots+'/phasing_output')
    html_graphs.write(
        "  <object data=\"" + chrom + '_' + sufix + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")


def plot_coverage_data_after_correction(html_graphs, args, chrom, ref_start_values, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2, sufix, loh_region_starts, loh_region_ends):
    fig = go.Figure()
    add_scatter_trace_coverage(fig, ref_start_values, haplotype_1_values, name='HP-1', text=None, yaxis=None,
                               opacity=0.7, color='firebrick')
    add_scatter_trace_coverage(fig, ref_start_values, haplotype_2_values, name='HP-2', text=None, yaxis=None,
                               opacity=0.7, color='steelblue')

    if not args.unphased_reads_coverage_disable:
        add_scatter_trace_coverage(fig, ref_start_values, unphased_reads_values, name='Unphased', text=None,
                                   yaxis=None, opacity=0.7, color='olive')
    plots_add_markers_lines(fig)

    if args.phaseblocks_enable:
        gaps_values_hp1 = np.full(len(haplotype_1_values_phasesets), 'None')
        gaps_values_hp2 = np.full(len(haplotype_2_values_phasesets), 'None')
        haplotype_1_phaseblocks_values = list(
            itertools.chain.from_iterable(zip(haplotype_1_values_phasesets, haplotype_1_values_phasesets, gaps_values_hp1)))
        haplotype_2_phaseblocks_values = list(
            itertools.chain.from_iterable(zip(haplotype_2_values_phasesets, haplotype_2_values_phasesets, gaps_values_hp2)))
        phaseblocks_positions_hp1 = list(
            itertools.chain.from_iterable(zip(ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, gaps_values_hp1)))
        phaseblocks_positions_hp2 = list(
            itertools.chain.from_iterable(zip(ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2, gaps_values_hp2)))

        add_scatter_trace_phaseblocks_seperate(fig, phaseblocks_positions_hp1, phaseblocks_positions_hp2, haplotype_1_phaseblocks_values,
                                      haplotype_2_phaseblocks_values)

    if loh_region_starts:
        for k, (start_loh, end_loh) in enumerate(zip(loh_region_starts, loh_region_ends)):
            fig.add_vrect(x0=start_loh, x1=end_loh, fillcolor="lightgrey", opacity=0.5, layer="below", line_width=0, )
    regions = get_chromosomes_regions(args)
    chroms = get_contigs_list(args.contigs)
    chrom_index = chroms.index(chrom)
    plots_layout_settings(fig, chrom, args, regions[chrom_index], args.cut_threshold)

    if args.pdf_enable:
        print_chromosome_pdf(fig, chrom, args.out_dir_plots+'/phasing_output')

    print_chromosome_html(fig, chrom + '_' + sufix, html_graphs, args.out_dir_plots+'/phasing_output')
    html_graphs.write(
        "  <object data=\"" + chrom + '_' + sufix + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

def change_point_detection(data, start, ends, args, chrom, html_graphs, hp, color):
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
        fig.add_vline(x=point*args.bin_size, y0=-10, y1=500, line_width=1, line_dash="dash",
                  line_color=color)

    plots_add_markers_lines(fig)
    plots_layout_settings(fig, chrom, args, ends[-1:][0], args.cut_threshold)

    print_chromosome_html(fig, chrom + '_hp_'  + str(hp), html_graphs, args.out_dir_plots+'/phasing_output')
    html_graphs.write("  <object data=\"" + chrom + '_hp_'  + str(hp)  + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")


def print_chromosome_pdf(fig, chrom, coverage_plots_path):
    filename = f"{os.path.join(coverage_plots_path, chrom + '.pdf')}"
    plotly.io.write_image(fig,  filename, format='pdf')

def print_chromosome_html(fig, chrom, html_graphs, coverage_plots_path):
    fname = f"{os.path.join(coverage_plots_path, chrom + '.html')}"
    plotly.offline.plot(fig, filename=fname,auto_open=False)

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
            title="<b>Coverage</b> (mean depth)",
            titlefont={"color": "dimgray"},
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
            title="<b>Copynumbers</b> (integers)",
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

def add_scatter_trace_phaseblocks_seperate(fig, phaseblocks_positions_hp1, phaseblocks_positions_hp2, haplotype_1_phaseblocks_values, haplotype_2_phaseblocks_values):
    fig.add_trace(go.Scatter(
        #legendgroup="group3",  # this can be any string, not just "group"
        #legendgrouptitle_text="Phaseblocks",
        x=phaseblocks_positions_hp1,
        y=haplotype_1_phaseblocks_values,
        name="HP-1",
        text=phaseblocks_positions_hp1,
        #yaxis="y5",
        line = dict(shape = 'spline', color = 'gray', width= 1, dash = 'solid'),
        mode='lines+markers',
        marker={"size": 5},
        opacity=0.5,
        marker_color=['dimgray', 'darkgray', 'white']*len(phaseblocks_positions_hp1),
        showlegend=True,
        marker_symbol='diamond-wide',
        hoverinfo = "x+name+y+text",
        #legendgroup="group2",
        #legendgrouptitle_text="Phaseblocks",
    ))

    fig.add_trace(go.Scatter(
        #legendgroup="group3",
        x=phaseblocks_positions_hp2,
        y=haplotype_2_phaseblocks_values,
        name="HP-2",
        text=phaseblocks_positions_hp2,
        #yaxis="y5",
        line = dict(shape = 'spline', color = 'green', width= 1, dash = 'solid'),
        mode='lines+markers',
        marker={"size": 5},
        opacity=0.5,
        marker_color=['darkgreen', 'limegreen', 'white']*len(phaseblocks_positions_hp2),
        showlegend=True,
        marker_symbol='diamond-wide',
        hoverinfo = "x+name+y+text",
        #legendgroup="group2",
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

def slice_list_sums(input):
    res, last = [[]], None
    for x in sorted(input):
        if last is None or abs(last - x) <= 1.5:
            res[-1].append(x)
        else:
            res.append([x])
        last = x
    first = [res[0][0]]
    print(res)

    return first + [sum(sub_list) / len(sub_list) for sub_list in res[1:]]

def loh_plots_genome(df_snps_ratios, args, df_loh_regions):
    # TODO genome-wide plots
    if not os.path.isdir(args.out_dir_plots + '/snps_loh_plots'):
        os.makedirs(args.out_dir_plots + '/snps_loh_plots')

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

    #ax = 20
    #ay = -30
    #add_annotation(fig, 960000000 + 250000000 + 120000000, 1, ax, ay, "LOH", '#2980b9')

    fig.add_annotation(text="Het SNPs ratio threshold: " + str(args.hets_ratio), x=960000000 + 250000000 + 120000000, y=0.95, showarrow=False)
    fig.update_yaxes(range=[0,1])
    fig.update_yaxes(range=[0,1])
    fig.update_yaxes(title_text="<b>SNPs Frequencies</b> (ratios)")
    fig.update_layout(width=1680, height=600,)
    fig.update_layout(legend=dict(orientation = 'h', xanchor = "center", x = 0.5, y= 1.1))
    #if args.pdf_enable:
    #    print_genome_pdf(fig, args.genome_name+'_loh', args.out_dir_plots + '/snps_loh_plots')

    fig.write_html(args.out_dir_plots + '/snps_loh_plots' +'/'+ args.genome_name + "_genome_snps_ratio_loh.html")