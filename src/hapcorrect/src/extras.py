import logging
logger = logging.getLogger()

def get_contigs_list(contigs):
    chroms_list_final = []
    #chrs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']
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