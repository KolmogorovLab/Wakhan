import os
import pandas as pd
from itertools import takewhile


def centromere_regions_blacklist(args, df_segs_hp1_, df_segs_hp2_):
    df_segs_hp1 = df_segs_hp1_.copy()
    df_segs_hp2 = df_segs_hp2_.copy()

    updated_hp1_segs = []
    updated_hp2_segs = []

    fileDir = os.path.dirname(__file__)  # os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)

    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        df_centm_chrom = df_centm[df_centm['chr'] == chrom]
        if not df_centm_chrom.empty:
            cents = [df_centm_chrom.start.values.tolist()[0], df_centm_chrom.end.values.tolist()[0]]
        else:
            cents = [0, 0]
        df_segs_hp_1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
        df_segs_hp_2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]

        for idx, s1 in df_segs_hp_1_chrom.iterrows():
            if (df_segs_hp_1_chrom.loc[idx, 'start'] == cents[0] or df_segs_hp_1_chrom.loc[idx, 'start'] == cents[0]+1) and df_segs_hp_1_chrom.loc[idx, 'end'] == cents[1]:# or df_segs_hp_1_chrom.loc[idx+1, 'start'] == cents[0] and df_segs_hp_1_chrom.loc[idx+1, 'end'] == cents[1]:
                df_segs_hp_1_chrom.loc[idx, 'depth'] = 3300
                df_segs_hp_1_chrom.loc[idx, 'state'] = 3300

        for idx, s2 in df_segs_hp_2_chrom.iterrows():
            if (df_segs_hp_2_chrom.loc[idx, 'start'] == cents[0] or df_segs_hp_2_chrom.loc[idx, 'start'] == cents[0]+1) and df_segs_hp_2_chrom.loc[idx, 'end'] == cents[1]:# or df_segs_hp_1_chrom.loc[idx+1, 'start'] == cents[0] and df_segs_hp_1_chrom.loc[idx+1, 'end'] == cents[1]:
                df_segs_hp_2_chrom.loc[idx, 'depth'] = 3300
                df_segs_hp_2_chrom.loc[idx, 'state'] = 3300

        updated_hp1_segs.append(df_segs_hp_1_chrom)
        updated_hp2_segs.append(df_segs_hp_2_chrom)

    segs_hp1 = pd.concat(updated_hp1_segs)
    segs_hp2 = pd.concat(updated_hp2_segs)

    return segs_hp1, segs_hp2


def centromere_regions_blacklist_bins(args, df_hp1_, df_hp2_, df_segs_hp1_updated_, df_segs_hp2_updated_):
    df_segs_hp1 = df_segs_hp1_updated_.copy()
    df_segs_hp2 = df_segs_hp2_updated_.copy()

    updated_hp1 = []
    updated_hp2 = []

    fileDir = os.path.dirname(__file__)  # os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)

    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        df_centm_chrom = df_centm[df_centm['chr'] == chrom]
        if not df_centm_chrom.empty:
            cents = [df_centm_chrom.start.values.tolist()[0], df_centm_chrom.end.values.tolist()[0]]
        else:
            cents = [0, 0]
        df_segs_hp_1_chrom = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
        df_segs_hp_2_chrom = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]

        df_hp_1_chrom = df_hp1_[df_hp1_['chr'] == chrom]
        df_hp_2_chrom = df_hp2_[df_hp2_['chr'] == chrom]
        df_hp_1_chrom = df_hp_1_chrom.reset_index(drop=True)
        df_hp_2_chrom = df_hp_2_chrom.reset_index(drop=True)

        ref_start_values = df_hp_1_chrom.start.values.tolist()

        for idx, s1 in df_segs_hp_1_chrom.iterrows():
            if (df_segs_hp_1_chrom.loc[idx, 'start'] == cents[0] or df_segs_hp_1_chrom.loc[idx, 'start'] == cents[0]+1) and df_segs_hp_1_chrom.loc[idx, 'end'] == cents[1]:# or df_segs_hp_1_chrom.loc[idx+1, 'start'] == cents[0] and df_segs_hp_1_chrom.loc[idx+1, 'end'] == cents[1]:
                internal_bins = [k for k in ref_start_values if k >= s1['start'] and k <= s1['end']]
                if internal_bins:
                    j = ref_start_values.index(internal_bins[0])
                    for l in range(len(internal_bins)):
                        df_hp_1_chrom.loc[j, 'hp1'] = 3300
                        j = j+1

        for idx, s2 in df_segs_hp_2_chrom.iterrows():
            if (df_segs_hp_2_chrom.loc[idx, 'start'] == cents[0] or df_segs_hp_2_chrom.loc[idx, 'start'] == cents[0]+1) and df_segs_hp_2_chrom.loc[idx, 'end'] == cents[1]:# or df_segs_hp_1_chrom.loc[idx+1, 'start'] == cents[0] and df_segs_hp_1_chrom.loc[idx+1, 'end'] == cents[1]:
                internal_bins = [k for k in ref_start_values if k >= s2['start'] and k <= s2['end']]
                if internal_bins:
                    j = ref_start_values.index(internal_bins[0])
                    for l in range(len(internal_bins)):
                        df_hp_2_chrom.loc[j, 'hp2'] = 3300
                        j = j+1

        updated_hp1.append(df_hp_1_chrom)
        updated_hp2.append(df_hp_2_chrom)

    segs_hp1 = pd.concat(updated_hp1)
    segs_hp2 = pd.concat(updated_hp2)

    return segs_hp1, segs_hp2


def extract_centromere_regions(args):
    fileDir = os.path.dirname(__file__)  # os.path.dirname(os.path.realpath('__file__'))
    cen_coord = os.path.join(fileDir, args.centromere)
    df_centm = csv_df_chromosomes_sorter(cen_coord, ['chr', 'start', 'end'])
    df_centm['start'].mask(df_centm['start'] == 1, 0, inplace=True)
    return df_centm


def get_contigs_list(contigs):
    #chrs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']
    chroms_list_final = []
    chroms = contigs.split(',')
    for chrom in chroms:
        chrom = chrom[len('chr'):] if chrom.startswith('chr') else chrom
        chrom = chrom.split('-')
        if len(chrom) > 1:
            chroms_list_final.extend(list(range(int(chrom[0]), int(chrom[1]) + 1)))
        else:
            chroms_list_final.extend(chrom)

    chroms_list_final = ['chr' + x if chroms[0].startswith('chr') else x for x in map(str, chroms_list_final)]
    return chroms_list_final


def chromosomes_sorter(label):
    # Strip "chr" prefix
    chrom = (label[3:] if label.lower().startswith('chr')
             else label)
    if chrom in ('X', 'Y'):
        key = (1000, chrom)
    else:
        # Separate numeric and special chromosomes
        nums = ''.join(takewhile(str.isdigit, chrom))
        chars = chrom[len(nums):]
        nums = int(nums) if nums else 0
        if not chars:
            key = (nums, '')
        elif len(chars) == 1:
            key = (2000 + nums, chars)
        else:
            key = (3000 + nums, chars)
    return key


def csv_df_chromosomes_sorter(path, names, sept='\t'):
    dataframe = pd.read_csv(path, sep=sept, names=names, header=None)
    dataframe['chr'] = dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', names[1]], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)


def df_chromosomes_sorter(dataframe, names, sept='\t'):
    dataframe['chr'] = dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', names[1]], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)


#TODO: likely bug in this function
def overlap_check(start, end, starts, ends):
  for i in range(len(starts)):
    #if (start < ends[i] and end > starts[i]) or (start <= starts[i] and end >= ends[i]):
    if start >= starts[i] and start <= ends[i]:
      return True, i
  return False, -1
