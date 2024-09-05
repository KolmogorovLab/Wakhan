import numpy as np

def get_all_breakpoints_data(edges, edges_chr, height, path):
    from vcf_parser import VCFParser
    my_parser = VCFParser(infile=path, split_variants=True, check_info=True)
    bp_junctions_inv = [[]]
    bp_junctions_dup_ins = [[]]
    bp_junctions_del = [[]]
    for variant in my_parser:
        if variant['CHROM'] == 'chrY':
            continue
        if "INV" in variant['ID']:
            bp_junctions_inv.append([variant['CHROM'], int(variant['POS'])])
        elif ("DUP" in variant['ID'] or "INS" in variant['ID']) and int(variant['info_dict']['SVLEN'][0]) > 10000:
            bp_junctions_dup_ins.append([variant['CHROM'], int(variant['POS'])])
        elif "DEL" in variant['ID'] and int(variant['info_dict']['SVLEN'][0]) > 10000:
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
    edges_chr = [edges_chr[i] + (edges_chr[i + 1] if i + 1 < n else []) for i in range(0, n, 2)]
    #BNDs and INVs
    for i, (a,b,c,d,e,f) in enumerate(edges_chr):
        hover_info.append(a+':'+str(b)+'-'+d+':'+str(e))

        if a == d:
            nr = 35
        else:
            nr = 75
        # if c == 1 and f == 1:
        #     colors.append('#DC3A3A')
        # elif c == 2 and f == 2:
        #     colors.append('#1D28C3')
        # else:
        #     colors.append('#737373')
        if [a,b] in bp_junctions_inv:
            colors.append('#2830DE')
        elif [a,b] in bp_junctions_dup_ins:
            colors.append('#178117')
        elif [a,b] in bp_junctions_del:
            colors.append('#CF0759')
        else:
            colors.append('#737373')

    #edges = [[0, 1000000], [0, 100000000], [0, 190000000], [358000001, 234000000], [68754445, 345000000], [6577777, 462000000]]
    #edges = [[358000001, 234000000]]

    for i, (j, k) in enumerate(edges):
        #print(edges_chr[i])
        # if j < k:
        #     tooltips.append(f'interactions({labels[j]}, {labels[k]})={interact_strength[i]}')
        # else:
        #     tooltips.append(f'interactions({labels[k]}, {labels[j]})={interact_strength[i]}')
        b0 = [j, 0.0]
        b2 = [k, 0.0]
        b1 = get_b1(b0, b2)

        a = dim_plus_1([b0, b1, b2], [1, height/(k-j), 1])
        pts = Rational_Bezier_curve(a, nr)
        xx.append(pts[nr // 2][0])  # abscissa of the middle point on the computed arc
        yy.append(pts[nr // 2][1])  # ordinate of the same point
        x, y = zip(*pts)

        # data.append(dict(type='scatter',
        #                  x=x,
        #                  y=y,
        #                  name='',
        #                  mode='lines',
        #                  line=dict(width=1, color='#6b8aca', shape='spline'),
        #                  hoverinfo='none',
        #                  showlegend=False
        #                  )
        #             )

        data.append([x, y, '',  'lines', dict(width=1.2, color=colors[i], shape='spline'), hover_info[i], '<br><b>Breakpoint info</b>: %{text}<br>', False])
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


