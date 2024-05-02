import numpy as np

def get_all_breakpoints_data(edges, edges_chr, height):
    #keys = sorted(set(interact_strength))
    #widths = [0.5 + k * 0.25 for k in range(5)] + [2 + k * 0.25 for k in range(4)] + [3, 3.25, 3.75, 4.25, 5, 5.25, 7]
    #d = dict(zip(keys, widths))
    #nwidths = [d[val] for val in interact_strength]

    data = []
    tooltips = []  # list of strings to be displayed when hovering the mouse over the middle of the circle arcs
    xx = []
    yy = []
    #X = list(range(L))  # node x-coordinates
    nr = 75


    #9:64526327 - 1:9769900
    #6:91468745 - 1:9769559

    n = len(edges_chr)
    hover_info = []
    colors = []
    widths = []
    edges_chr = [edges_chr[i] + (edges_chr[i + 1] if i + 1 < n else []) for i in range(0, n, 2)]
    for i, (a,b,c,d,e,f) in enumerate(edges_chr):
        hover_info.append(a+':'+str(b)+'-'+d+':'+str(e))
        if c == 1 and f == 1:
            colors.append('#DC3A3A')
        elif c == 2 and f == 2:
            colors.append('#1D28C3')
        else:
            colors.append('#737373')



    #edges = [[0, 1000000], [0, 100000000], [0, 190000000], [358000001, 234000000], [68754445, 345000000], [6577777, 462000000]]
    #edges = [[358000001, 234000000]]

    for i, (j, k) in enumerate(edges):
        print(edges_chr[i])
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

        data.append([x, y, '',  'lines', dict(width=1, color=colors[i], shape='spline'), hover_info[i], '<br><b>Breakpoint info</b>: %{text}<br>', False])
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

class bps_sample(object):
    __slots__ = ('bp_1_id', 'bp_1_pos', 'bp_1_hp', 'bp_2_id', 'bp_2_pos', 'bp_2_hp', 'status')
    def __init__(self, bp_1_id, bp_1_pos, bp_1_hp, bp_2_id, bp_2_pos, bp_2_hp, status):
        self.bp_1_id = bp_1_id
        self.bp_1_pos = bp_1_pos
        self.bp_2_id = bp_2_id
        self.bp_2_pos = bp_2_pos
        self.bp_1_hp = bp_1_hp
        self.bp_2_hp = bp_2_hp
        self.status = False

    def call_bp_cn(self):
        #call CN on BP and check if overlaps

        #if no_overlaps:
        #    self.status = False
        print('call CN on BP and check if overlaps')
        self.status = True

    def bp_cn_overlaps(self):
        self.call_bp_cn()
        if self.status == None:
            return False
        else:
            return f"{self.bp_1_id}:{self.bp_1_pos}:{self.bp_1_hp}-{self.bp_2_id}:{self.bp_2_pos}:{self.bp_2_hp}"

def sv_vcf_bps_cn_check(path, df_segs_hp1, df_segs_hp2):
    #########################################
    from vcf_parser import VCFParser
    my_parser = VCFParser(infile=path, split_variants=True, check_info=True)
    from collections import defaultdict
    sample_list = defaultdict(list)
    sample_single_list = defaultdict(list)
    bp_junctions = [[]]
    bp_junctions_chr = [[]]
    for variant in my_parser:
        if "INV" in variant['ID'] or "INS" in variant['ID']:
            continue
        elif "BND" in variant['ID']:
            s = variant['ALT']
            for ch in ['[', ']', 'N']:
                if ch in s:
                    s = s.replace(ch, '')
            chr2_id = s.split(':')[0]
            chr2_end = int(s.split(':')[1])
            hp = 0
            if 'HP' in variant['info_dict']:
                hp = int(variant['info_dict']['HP'][0])




            # for k, v in variant['info_dict'].items():
            #     if k == 'CHR2':
            #         chr2_id = v[0]
            #     elif k == 'END':
            #         chr2_end = v[0]
            if chr2_id == '' or chr2_id == 'chrX' or chr2_id == 'chrY' or variant['CHROM'] == 'chrX' or variant['CHROM'] == 'chrY' or \
                    (variant['CHROM'] == chr2_id and chr2_end - int(variant['POS']) < 50000):
                continue
            else:
                sample_list[variant['ID']] = bps_sample(variant['CHROM'], int(variant['POS']), hp, chr2_id, chr2_end, hp, False)
        else:
            hp = 0
            if 'HP' in variant['info_dict']:
                hp = int(variant['info_dict']['HP'][0])

            if abs((int(variant['info_dict']['SVLEN'][0]) + int(variant['POS'])) - int(variant['POS'])) < 50000 or \
                    variant['CHROM'] == 'chrX' or variant['CHROM'] == 'chrY':
                continue
            sample_single_list[variant['ID']] = bps_sample(variant['CHROM'], int(variant['POS']), hp, variant['CHROM'],  int(variant['info_dict']['SVLEN'][0]) + int(variant['POS']) + 1, hp, False )

    for j, dict in enumerate(sample_list.items()):
        if dict[1].bp_1_id == dict[1].bp_2_id and dict[1].bp_1_pos == dict[1].bp_2_pos:
            continue
        bp_junctions_chr.append([dict[1].bp_1_id, dict[1].bp_1_pos, dict[1].bp_1_hp])
        bp_junctions_chr.append([dict[1].bp_2_id, dict[1].bp_2_pos, dict[1].bp_2_hp])

    for j, dict in enumerate(sample_list.items()):
        if dict[1].bp_1_id == dict[1].bp_2_id and dict[1].bp_1_pos == dict[1].bp_2_pos:
            continue
        bp_junctions.append([dict[1].bp_1_id, dict[1].bp_1_pos])
        bp_junctions.append([dict[1].bp_2_id, dict[1].bp_2_pos])

    return bp_junctions[1:], bp_junctions_chr[1:]
    #########################################