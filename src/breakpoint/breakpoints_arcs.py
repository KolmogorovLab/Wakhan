import numpy as np
import pandas as pd
import logging
from vcf_parser import VCFParser

logger = logging.getLogger()

from src.utils.chromosome import get_contigs_list, df_chromosomes_sorter

def get_all_breakpoints_data(edges, edges_chr, height, path, args):
    _logging_level = logger.level
    logger.setLevel(logging.CRITICAL)
    my_parser = VCFParser(infile=path, split_variants=True, check_info=True)
    logger.setLevel(_logging_level)

    chroms = get_contigs_list(args.contigs)
    bp_junctions_inv = [[]]
    bp_junctions_ins = [[]]
    bp_junctions_dup = [[]]
    bp_junctions_del = [[]]
    bp_junctions_sbnd = [[]]
    for variant in my_parser:
        if not variant['CHROM'] in chroms:
            continue
        if 'SVLEN' not in variant['info_dict']: #rare bug in Severus output
            continue
        if variant['info_dict']['SVTYPE'][0] == 'INV':
            bp_junctions_inv.append([variant['CHROM'], int(variant['POS'])])
        elif variant['info_dict']['SVTYPE'][0] == 'sBND':
            bp_junctions_sbnd.append([variant['CHROM'], int(variant['POS'])])
        elif (variant['info_dict']['SVTYPE'][0] == 'DUP') and int(variant['info_dict']['SVLEN'][0]) > args.breakpoints_min_length:
            bp_junctions_dup.append([variant['CHROM'], int(variant['POS'])])
        elif (variant['info_dict']['SVTYPE'][0] == 'INS') and int(variant['info_dict']['SVLEN'][0]) > args.breakpoints_min_length:
            bp_junctions_ins.append([variant['CHROM'], int(variant['POS'])])
        elif variant['info_dict']['SVTYPE'][0] == 'DEL' and int(variant['info_dict']['SVLEN'][0]) > args.breakpoints_min_length:
            bp_junctions_del.append([variant['CHROM'], int(variant['POS'])])

    #keys = sorted(set(interact_strength))
    #widths = [0.5 + k * 0.25 for k in range(5)] + [2 + k * 0.25 for k in range(4)] + [3, 3.25, 3.75, 4.25, 5, 5.25, 7]
    #d = dict(zip(keys, widths))
    #nwidths = [d[val] for val in interact_strength]

    data = []
    tooltips = []  # list of strings to be displayed when hovering the mouse over the middle of the circle arcs
    xx = []
    yy = []
    #X = list(range(L))  # node x-coordinates

    #9:64526327 - 1:9769900
    #6:91468745 - 1:9769559

    n = len(edges_chr)
    hover_info = []
    colors = []
    widths = []
    names_conv = []
    edges_chr = [edges_chr[i] + (edges_chr[i + 1] if i + 1 < n else []) for i in range(0, n, 2)]
    #BNDs and INVs
    for i, (a,b,c,hp1, dv, d,e,f,hp2, dv) in enumerate(edges_chr):
        if ('_2' in c and c.replace('_2', '_1') in names_conv) or ('_1' in c and c.replace('_1', '_2') in names_conv):
            hp1 = '|'.join(reversed(hp1.split('|')))
            hp2 = '|'.join(reversed(hp2.split('|')))
        names_conv.append(c)

        #sort coordinates
        data_df = {'chr': [a, d], 'start': [b, e], 'hps': [hp1, hp2]}
        df = pd.DataFrame(data_df)
        df_sorted = df_chromosomes_sorter(df, ['chr','start'])
        a = df_sorted.iloc[0,0]
        b = df_sorted.iloc[0,1]
        d = df_sorted.iloc[1,0]
        e = df_sorted.iloc[1,1]
        hp = df_sorted.iloc[0,2]

        if a == d:
            nr = 35
        else:
            nr = 75

        if [a,b] in bp_junctions_inv:
            colors.append('#2830DE')
            hover_info.append('HP='+str(hp)+ '- ' + 'INV - ' + a + ':' + str(b) + '-' + d + ':' + str(e))
        elif [a,b] in bp_junctions_dup:
            colors.append('#178117')
            hover_info.append('HP='+str(hp)+ '- ' + 'DUP - ' + a + ':' + str(b) + '-' + d + ':' + str(e))
        elif [a,b] in bp_junctions_sbnd:
            colors.append('#7f8c8d')
            hover_info.append('HP='+str(hp)+ '- ' + 'sBND - ' + a + ':' + str(b) + '-' + d + ':' + str(e))
        elif [a,b] in bp_junctions_ins:
            colors.append('#e0cf03')
            hover_info.append('HP='+str(hp)+ '- ' + 'INS - ' + a + ':' + str(b) + '-' + d + ':' + str(e))
        elif [a,b] in bp_junctions_del:
            colors.append('#CF0759')
            hover_info.append('HP='+str(hp)+ '- ' + 'DEL - ' + a + ':' + str(b) + '-' + d + ':' + str(e))
        else:
            colors.append('#737373')
            hover_info.append('HP='+str(hp)+ '- ' + 'BND - ' + a + ':' + str(b) + '-' + d + ':' + str(e))

    for i, (j, k) in enumerate(edges):
        b0 = [j, 0.0]
        b2 = [k, 0.0]
        b1 = get_b1(b0, b2)

        a = dim_plus_1([b0, b1, b2], [1, height/(k-j), 1])
        pts = Rational_Bezier_curve(a, nr)
        xx.append(pts[nr // 2][0])  # abscissa of the middle point on the computed arc
        yy.append(pts[nr // 2][1])  # ordinate of the same point
        x, y = zip(*pts)

        data.append([x, y, '',  'lines', dict(width=1.2, color=colors[i], shape='spline'), hover_info[i], hover_info[i], False])

    return data
#Foolowing Bezier curve code is adopted from https://notebook.community/empet/Plotly-plots/Arc-diagram-Force-Awakens
def get_b1(b0, b2):
    # b0, b1 list of x, y coordinates
    if len(b0) != len(b2) != 2:
        raise ValueError('b0, b1 must be lists of two elements')
    b1 = 0.5 * (np.asarray(b0)+np.asarray(b2))+\
         0.5 * np.array([0,1.0]) * np.sqrt(3) * np.linalg.norm(np.array(b2)-np.array(b0))
    return b1.tolist()

def dim_plus_1(b, w):#lift the points b0, b1, b2 to 3D points a0, a1, a2 (see Gallier book)
    #b is a list of 3 lists of 2D points, i.e. a list of three 2-lists
    #w is a list of numbers (weights) of len equal to the len of b
    if not isinstance(b, list) or  not isinstance(b[0], list):
        raise ValueError('b must be a list of three 2-lists')
    if len(b) != len(w)   != 3:
        raise ValueError('the number of weights must be  equal to the nr of points')
    else:
        a = np.array([point + [w[i]] for (i, point) in enumerate(b)])
        a[1, :2] *= w[1]
        return a
def Bezier_curve(bz, nr): #the control point coordinates are passed in a list bz=[bz0, bz1, bz2]
    # bz is a list of three 2-lists
    # nr is the number of points to be computed on each arc
    t = np.linspace(0, 1, nr)
    #for each parameter t[i] evaluate a point on the Bezier curve with the de Casteljau algorithm
    N = len(bz)
    points = [] # the list of points to be computed on the Bezier curve
    for i in range(nr):
        aa = np.copy(bz)
        for r in range(1, N):
            aa[:N-r,:] = (1-t[i]) * aa[:N-r,:] + t[i] * aa[1:N-r+1,:]  # convex combination of points
        points.append(aa[0,:])
    return np.array(points)

def Rational_Bezier_curve(a, nr):
    discrete_curve = Bezier_curve(a, nr )
    return [p[:2]/p[2] for p in discrete_curve]


