import pandas as pd
from itertools import takewhile


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
