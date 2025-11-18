import csv
import pandas as pd
import os
import statistics
import logging

logger = logging.getLogger()

from src.cna.copynumber import subclonal_values_adjusted, integers_values_adjusted, add_confidence_score_cn_segemnts
from src.utils.chromosome import csv_df_chromosomes_sorter, df_chromosomes_sorter, overlap_check
from src.utils.chromosome import get_contigs_list


def write_df_csv_header(df, path):
    fp = open(path, 'a')
    fp.write('#chr: chromosome number\n')
    fp.write('#start: start address gene\n')
    fp.write('#end: end address of gene\n')
    fp.write('#gene: name of gene\n')
    fp.write('#actual_coverage_hp1: actual gene coverage value of HP-1\n')
    fp.write('#actual_coverage_hp2: actual gene coverage value of HP-2\n')
    fp.write('#adjusted_coverage_hp1: adjusted gene coverage value of HP-1\n')
    fp.write('#adjusted_coverage_hp2: adjusted gene coverage value of HP-2\n')
    fp.write('#copynumber_state_hp1: gene copy number state of HP-1\n')
    fp.write('#copynumber_state_hp2: gene copy number state of HP-2\n')
    fp.write('#chr\tstart\tend\tgene\tactual_coverage_hp1\tactual_coverage_hp2\tadjusted_coverage_hp1\tadjusted_coverage_hp2\thp1_state\thp2_state\n')

    df.to_csv(fp, sep='\t', index=False, header=False)


def genes_phase_correction(df_genes, df_segs_hp1, df_segs_hp2, args, centers, integer_fractional_centers):
    chroms = get_contigs_list(args.contigs)
    df_genes_updated = []

    for index, chrom in enumerate(chroms):
        df_gene = df_genes[df_genes['chr'] == chrom]
        df_seg_hp1 = df_segs_hp1[df_segs_hp1['chromosome'] == chrom]
        df_seg_hp2 = df_segs_hp2[df_segs_hp2['chromosome'] == chrom]

        seg_coverage_hp1 = df_seg_hp1.depth.values.tolist()
        seg_coverage_hp2 = df_seg_hp2.depth.values.tolist()

        seg_state_hp1 = df_seg_hp1.state.values.tolist()
        seg_state_hp2 = df_seg_hp2.state.values.tolist()

        genes_starts = df_gene.start.values.tolist()
        genes_ends = df_gene.end.values.tolist()
        genes_coverage_hp1 = df_gene.hp1.values.tolist()
        genes_coverage_hp2 = df_gene.hp2.values.tolist()

        gene_hp1_state = [0 for x in range(len(genes_coverage_hp1))]
        gene_hp2_state = [0 for x in range(len(genes_coverage_hp2))]

        for i, (start_b, end_b) in enumerate(zip(genes_starts, genes_ends)):
            hp_1_val = 0
            hp_2_val = 0
            check, index = overlap_check(start_b, end_b, df_seg_hp1.start.values.tolist(), df_seg_hp1.end.values.tolist())
            if check:
                hp_1_val = centers[integers_values_adjusted(seg_state_hp1[index], centers)]
                #gene_hp1_state[i] = integer_fractional_centers[centers.index(min(centers, key=lambda x:abs(x - seg_coverage_hp1[index])))]  #(integer_fractional_centers[centers.index(seg_coverage_hp1[index])])
            check, index = overlap_check(start_b, end_b, df_seg_hp2.start.values.tolist(), df_seg_hp2.end.values.tolist())
            if check:
                hp_2_val = centers[integers_values_adjusted(seg_state_hp2[index], centers)]
                #gene_hp2_state[i] = integer_fractional_centers[centers.index(min(centers, key=lambda x:abs(x - seg_coverage_hp2[index])))]

            # if hp_1_val == 0 and hp_2_val == 0:
            #     gene_hp1_state[i] = integers_values_adjusted(genes_coverage_hp1[i], centers)
            #     gene_hp2_state[i] = integers_values_adjusted(genes_coverage_hp2[i], centers)
            #     continue

            #print(genes_coverage_hp1[i], genes_coverage_hp2[i])
            #print(abs(genes_coverage_hp1[i] -  hp_2_val) + abs(genes_coverage_hp2[i] -  hp_1_val))
            #print(abs(genes_coverage_hp1[i] -  hp_1_val) + abs(genes_coverage_hp2[i] -  hp_2_val))

            # if abs(genes_coverage_hp1[i] -  hp_2_val) < abs(genes_coverage_hp1[i] -  hp_1_val):
            #     new_hp2 = genes_coverage_hp1[i]
            #     new_hp1 = genes_coverage_hp2[i]
            #     genes_coverage_hp2[i] = new_hp2
            #     genes_coverage_hp1[i] = new_hp1

            genes_coverage_hp1[i] = hp_1_val #centers[integers_values_adjusted(genes_coverage_hp1[i], centers)]
            genes_coverage_hp2[i] = hp_2_val #centers[integers_values_adjusted(genes_coverage_hp2[i], centers)]

            gene_hp1_state[i] = integers_values_adjusted(genes_coverage_hp1[i], centers)
            gene_hp2_state[i] = integers_values_adjusted(genes_coverage_hp2[i], centers)

        df_genes_updated.append(pd.DataFrame(list(zip(df_gene.chr.values.tolist(), df_gene.start.values.tolist(), df_gene.end.values.tolist(), \
                     df_gene.gene.values.tolist(), [round(l,2) for l in df_gene.hp1.values.tolist()], [round(l,2) for l in df_gene.hp2.values.tolist()], [round(l,2) for l in genes_coverage_hp1], [round(l,2) for l in genes_coverage_hp2], gene_hp1_state, gene_hp2_state)),
            columns=['chr', 'start', 'end', 'gene', 'actual_coverage_hp1', 'actual_coverage_hp2', 'adjusted_coverage_hp1', 'adjusted_coverage_hp2', 'copynumber_hp1_state', 'copynumber_hp2_state']))

    return pd.concat(df_genes_updated)


def update_genes_phase_corrected_coverage(args, df_segs_hp1, df_segs_hp2, p_value, centers, integer_fractional_centers, is_half):

    df_genes = csv_df_chromosomes_sorter(args.out_dir_plots + '/coverage_data/cancer_genes_coverage.csv', ['chr','start','end','gene', 'hp1', 'hp2'])
    #write_df_csv_header(df_genes, args.out_dir_plots + '/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output/' + 'cancer_genes_coverage.csv')
    df_genes = genes_phase_correction(df_genes, df_segs_hp1, df_segs_hp2, args, centers, integer_fractional_centers)
    if is_half:
        write_df_csv_header(df_genes, args.out_dir_plots + '/wgd/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output/' + 'genes_copynumber_states.bed')
    else:
        write_df_csv_header(df_genes, args.out_dir_plots + '/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output/' + 'genes_copynumber_states.bed')

    return df_genes


def snps_frequencies_chrom_genes(df_snps_frequencies, args):
    df_chroms = []

    df_genes_all = csv_df_chromosomes_sorter(args.cancer_genes, ['chr', 'start', 'end', 'gene'])
    df_empty = pd.DataFrame(columns=['chr', 'start', 'end', 'gene', 'hp1', 'hp2'])
    if args.user_input_genes:
        with open(args.user_input_genes, 'r') as file:
            entries = file.readlines()
        if "\t" in entries[0]:
            entries = [line.rstrip('\n').split('\t')[3] for line in entries]
        else:
            entries = [line.rstrip('\n') for line in entries]
        if args.reference_name:#not entries[0] in ['chm13','grch38']:
            ref = args.reference_name #default
        else:
            ref = 'grch38'
        prefix, filename = os.path.split(args.cancer_genes)
        ref_name = prefix + '/' + ref + '_genes.tsv'
        chroms = []
        starts = []
        ends = []
        genes = []
        with open(ref_name, 'r') as file:
            tsv_reader = csv.reader(file, delimiter='\t')
            for row in tsv_reader:
                if row[3] in entries:
                    chroms.append(row[0])
                    starts.append(row[1])
                    ends.append(row[2])
                    genes.append(row[3])

        data = {"chr": chroms, "start": starts, "end": ends, "gene": genes}
        df_user_genes = pd.DataFrame(data)
        if not df_user_genes.empty:
            write_df_csv(df_user_genes, args.out_dir_plots + '/data_phasing/user_genes.csv')
            df_genes_all = csv_df_chromosomes_sorter(args.out_dir_plots + '/data_phasing/user_genes.csv', ['chr', 'start', 'end', 'gene'])

    chroms = get_contigs_list(args.contigs)
    for index, chrom in enumerate(chroms):
        df = df_snps_frequencies[df_snps_frequencies['chr'] == chrom]
        df_genes = df_genes_all[df_genes_all['chr'] == chrom]

        if not df_genes.empty:

            # df = dict(tuple(df_snps.groupby('hp')))
            haplotype_1_position = df.pos.values.tolist()
            haplotype_1_coverage = df.freq_value_b.values.tolist()
            haplotype_2_position = df.pos.values.tolist()
            haplotype_2_coverage = df.freq_value_a.values.tolist()

            snps_haplotype1_mean = []
            for index, (i,j) in enumerate(zip(df_genes.start.values.tolist(), df_genes.end.values.tolist())):
                sub_list = []
                try:
                    sub_list = haplotype_1_coverage[haplotype_1_position.index(
                        min(haplotype_1_position, key=lambda x:abs(x-i))):haplotype_1_position.index(
                        min(haplotype_1_position, key=lambda x:abs(x-j)))]
                except ValueError:
                    logger.info('No Hets pileup found!')
                if sub_list:
                    snps_haplotype1_mean.append(statistics.mean(sub_list))
                else:
                    snps_haplotype1_mean.append(0)

            snps_haplotype2_mean = []
            for index, (i, j) in enumerate(zip(df_genes.start.values.tolist(), df_genes.end.values.tolist())):
                sub_list = []
                try:
                    sub_list = haplotype_2_coverage[haplotype_2_position.index(
                        min(haplotype_2_position, key=lambda x: abs(x - i))):haplotype_2_position.index(
                        min(haplotype_2_position, key=lambda x: abs(x - j)))]
                except ValueError:
                    logger.info('No Hets pileup found!')
                if sub_list:
                    snps_haplotype2_mean.append(statistics.mean(sub_list))
                else:
                    snps_haplotype2_mean.append(0)

            #snps_mean = [round(i + j, 2) for i, j in zip(snps_haplotype1_mean, snps_haplotype2_mean)]
            if len(snps_haplotype1_mean) == 0:
                df_chroms.append(df_empty)
            else:
                df_chroms.append(pd.DataFrame(list(zip(df_genes.chr.values.tolist(), df_genes.start.values.tolist(), df_genes.end.values.tolist(), \
                                         df_genes.gene.values.tolist(), snps_haplotype1_mean, snps_haplotype2_mean)),
                                 columns=['chr', 'start', 'end', 'gene', 'hp1', 'hp2']))
    if len(df_chroms):
        return pd.concat(df_chroms)
    else:
        return df_empty


def genes_segments_list(bam, args):
    head, tail = os.path.split(bam)

    df_genes = csv_df_chromosomes_sorter(args.cancer_genes, ['chr', 'start', 'end', 'gene'])
    starts = df_genes.start.values.tolist()
    ends = df_genes.end.values.tolist()

    bed = []
    for ind, chrom in enumerate(df_genes.chr.values.tolist()):
        bed.append([tail, chrom, starts[ind], ends[ind]])

    return bed


def genes_segments_coverage(genes_coverage, args):

    df_genes = csv_df_chromosomes_sorter(args.cancer_genes, ['chr', 'start', 'end', 'gene'])

    hp1 = []
    hp2 = []
    hp3 = []
    for items in genes_coverage:
        hp1.append(round(float(items.split('\t')[3]), 2))
        hp2.append(round(float(items.split('\t')[4]), 2))
        #hp3.append(round(float(items.split('\t')[5]), 2))

    return pd.DataFrame(list(zip(df_genes.chr.values.tolist(), df_genes.start.values.tolist(), df_genes.end.values.tolist(),
                                 df_genes.gene.values.tolist(), hp1, hp2)),
                         columns=['chr', 'start', 'end', 'gene', 'hp1', 'hp2'])

