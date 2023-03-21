import pandas as pd
import numpy as np
import os

chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
          'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']  # , 'chrX', 'chrY']

def get_phasesets_bins(bam, phasesets, bin_size):
    indices, values = remove_overlapping_and_small_phasesets(phasesets, bin_size)
    head, tail = os.path.split(bam)
    bed = []
    for ind, chrom in enumerate(chroms) :
        for i in range(0, len(values[ind])-1, 2): #for i in range(len(values[ind])-1):
            start = values[ind][i]
            end = values[ind][i+1]
            bed.append([tail, chrom, start, end])
    return bed

def remove_overlapping_and_small_phasesets(phasesets, bin_size):
    dfs = pd.read_csv(phasesets, sep='\t', names=['chr', 'pos', 'qual', 'filter', 'ps', 'gt', 'dp', 'vaf'])
    values_all = []
    indices_all = []
    for index, chrom in enumerate(chroms):
        df = dfs[dfs['chr'] == chrom]
        unique_ps_by_chr = df.groupby('chr', group_keys=True)['ps'].apply(lambda x: list(np.unique(x)))
        # Find overlapping phase blocks locations
        for unique in unique_ps_by_chr[0]:
            indices_overlapping = []
            indices = np.where(df["ps"] == unique)
            list_new = []
            for i in indices[0]:
                list_new.append(i)
        for i in range(len(list_new) - 1):
            if not list_new[i + 1] - list_new[i] == 1:
                indices_overlapping.append(list_new[i] + 1)
                indices_overlapping.append(list_new[i + 1] - list_new[i] - 1)
        for i in range(0, len(indices_overlapping) - 1, 2):
            df = df.drop(labels=range(indices_overlapping[i], indices_overlapping[i] + indices_overlapping[i + 1]),
                         axis=0)
        final = []
        # Remove overlapping phase blocks
        ps = df.ps.values.tolist()
        pos = df.pos.values.tolist()
        unique_ps_by_chr = df.groupby('chr')['ps'].apply(lambda x: list(np.unique(x)))
        for x in range(len(unique_ps_by_chr[0]) - 1):
            if (unique_ps_by_chr[0][x + 1] - unique_ps_by_chr[0][x] > bin_size):
                final.append(pos[ps.index(unique_ps_by_chr[0][x])])
                final.append(pos[len(ps) - 1 - ps[::-1].index(unique_ps_by_chr[0][x])])
        index = []
        # Remove small (< bin_size) phase blocks
        for i in range(len(final)):
            init = 0
            for m in range(1, final[-1], bin_size):
                init += 1
                if final[i] > m and final[i] < m + bin_size:
                    index.append(init)

        indices_all.append(index)
        values_all.append(final)
    return indices_all, values_all

def phaseblock_flipping(haplotype_1_values, haplotype_2_values, ref_start_values, \
                    haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets):
    haplotype_block_changed = []
    for i in range(len(ref_start_values_phasesets) - 1):
        if haplotype_1_values_phasesets[i] > haplotype_2_values_phasesets[i]:
            haplotype_block_changed.append([ref_start_values_phasesets[i], ref_end_values_phasesets[i] + 1])

    for i in range(len(ref_start_values)):
        for j in range(len(haplotype_block_changed)): #TODO update it
            if (ref_start_values[i] >= haplotype_block_changed[j][0] and ref_start_values[i] <= haplotype_block_changed[j][1]):
                new_hp2 = haplotype_2_values[i]
                new_hp1 = haplotype_1_values[i]
                haplotype_1_values[i] = new_hp2
                haplotype_2_values[i] = new_hp1
            elif haplotype_1_values[i] > haplotype_2_values[i]:
                new_hp2 = haplotype_2_values[i]
                new_hp1 = haplotype_1_values[i]
                haplotype_1_values[i] = new_hp2
                haplotype_2_values[i] = new_hp1
                break
            else:
                haplotype_1_values[i] = haplotype_1_values[i]
                haplotype_2_values[i] = haplotype_2_values[i]

    return haplotype_1_values, haplotype_2_values
    #haplotype_1_values_updated.extend(haplotype_1_values)
    #haplotype_2_values_updated.extend(haplotype_2_values)