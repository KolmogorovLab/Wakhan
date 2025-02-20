import logging

logger = logging.getLogger()

def get_contigs_list(contigs):
    chrs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']
    chroms_list_final = []
    chroms = contigs.split(',')
    for chrom in chroms:
        chrom = chrom[len('chr'):] if chrom.startswith('chr') else chrom
        chrom = chrom.split('-')
        if len(chrom) > 1:
            chroms_list_final.extend(list(range(int(chrom[0]), int(chrom[1]) + 1)))
        else:
            chroms_list_final.extend(chrom)

    chroms_list_final = ['chr' + x if chroms[0].startswith('chr') and x in chrs else x for x in map(str, chroms_list_final)]
    return chroms_list_final

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

#def sv_vcf_bps_cn_check(path, df_segs_hp1, df_segs_hp2):
def sv_vcf_bps_cn_check(path, args):
    #########################################
    from vcf_parser import VCFParser
    my_parser = VCFParser(infile=path, split_variants=True, check_info=True)
    from collections import defaultdict
    sample_list = defaultdict(list)
    sample_single_list = defaultdict(list)
    chroms = get_contigs_list(args.contigs)

    bp_junctions = [[]]
    bp_junctions_bnd = [[]]

    bp_junctions_chr = [[]]
    bp_junctions_values = [[]]
    bp_junctions_single = [[]]
    bp_junctions_single_chr = [[]]

    for variant in my_parser:
        #if ("INV" in variant['ID'] and variant['info_dict']['DETAILED_TYPE'] == ['reciprocal_inv']) or  "INS" in variant['ID']:
        if variant['info_dict']['SVTYPE'][0] == 'INV' or variant['info_dict']['SVTYPE'][0] == 'DUP' or ((variant['info_dict']['SVTYPE'][0] == 'INS' or variant['info_dict']['SVTYPE'][0] == 'DEL') and int(variant['info_dict']['SVLEN'][0]) > args.breakpoints_min_length):
            if not variant['CHROM'] in chroms:
                continue
            hp = 0
            if 'HP' in variant['info_dict']:
                hp = int(variant['info_dict']['HP'][0])

            # sample_list[variant['ID']] = bps_sample(variant['CHROM'], int(variant['POS']), hp, variant['CHROM'], int(variant['POS'])+1, hp, False)
            # sample_list[variant['ID']] = bps_sample(chr2_id, chr2_end, hp, chr2_id, chr2_end+1, hp, False)
            bp_junctions_bnd.append([variant['CHROM'], int(variant['POS'])])
            bp_junctions_bnd.append([variant['CHROM'], int(variant['POS'])+int(variant['info_dict']['SVLEN'][0])])
            sample_list[variant['ID']] = bps_sample(variant['CHROM'], int(variant['POS']), variant['ID'], variant['CHROM'], int(variant['POS'])+int(variant['info_dict']['SVLEN'][0]), variant['ID'], False)

        elif variant['info_dict']['SVTYPE'][0] == 'sBND' or 'sBND' in variant['ID']:
            hp = 0
            if 'HP' in variant['info_dict']:
                hp = int(variant['info_dict']['HP'][0])
            bp_junctions_bnd.append([variant['CHROM'], int(variant['POS'])])
            bp_junctions_bnd.append([variant['CHROM'], int(variant['POS']) + 1])
            sample_list[variant['ID']] = bps_sample(variant['CHROM'], int(variant['POS']), variant['ID'], variant['CHROM'], int(variant['POS']) + 1, variant['ID'], False)

        elif variant['info_dict']['SVTYPE'][0] == 'BND':
            s = variant['ALT']
            for ch in ['[', ']', 'N']:
                if ch in s:
                    s = s.replace(ch, '')
            chr2_id = s.split(':')[0]
            chr2_end = int(s.split(':')[1])

            if not variant['CHROM'] in chroms or not chr2_id in chroms:
                continue
            hp = 0
            if 'HP' in variant['info_dict']:
                hp = int(variant['info_dict']['HP'][0])

            bp_junctions_bnd.append([variant['CHROM'], int(variant['POS'])])
            bp_junctions_bnd.append([chr2_id, chr2_end])
            sample_list[variant['ID']] = bps_sample(variant['CHROM'], int(variant['POS']), variant['ID'], chr2_id, chr2_end, variant['ID'], False)

    bp_junctions = sorted(bp_junctions[1:], key=lambda x: (x[0], x[1]))
    bp_junctions_bnd = sorted(bp_junctions_bnd[1:], key=lambda x: (x[0], x[1]))

    for j, dict in enumerate(sample_list.items()):
        if dict[1].bp_1_id == dict[1].bp_2_id and dict[1].bp_1_pos == dict[1].bp_2_pos:
            continue
        bp_junctions_chr.append([dict[1].bp_1_id, dict[1].bp_1_pos, dict[1].bp_1_hp])
        bp_junctions_chr.append([dict[1].bp_2_id, dict[1].bp_2_pos, dict[1].bp_2_hp])

    for j, dict in enumerate(sample_list.items()):
        if dict[1].bp_1_id == dict[1].bp_2_id and dict[1].bp_1_pos == dict[1].bp_2_pos:
            continue
        bp_junctions_values.append([dict[1].bp_1_id, dict[1].bp_1_pos])
        bp_junctions_values.append([dict[1].bp_2_id, dict[1].bp_2_pos])

    for j, dict in enumerate(sample_single_list.items()):
        if dict[1].bp_1_id == dict[1].bp_2_id and dict[1].bp_1_pos == dict[1].bp_2_pos:
            continue
        bp_junctions_single.append([dict[1].bp_1_id, dict[1].bp_1_pos, dict[1].bp_2_pos])
        bp_junctions_single_chr.append([dict[1].bp_1_id, dict[1].bp_1_pos, dict[1].bp_2_pos, dict[1].bp_1_hp])

    return bp_junctions_single[1:], bp_junctions_single_chr[1:], bp_junctions_values[1:], bp_junctions_chr[1:], bp_junctions, bp_junctions_bnd