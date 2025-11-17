import csv
import os

from src.cna.copynumber import subclonal_values_adjusted, integers_values_adjusted, add_confidence_score_cn_segemnts
from src.breakpoint.breakpoints import sv_vcf_bps_cn_check


def write_df_csv(df, file_name):
    df.to_csv(file_name, sep='\t', index=False, header=False, mode='w')


def write_header_comments(header, header_comments, output, args):
    with open(args.out_dir_plots + '/bed_output/' + output, 'a') as fp:
        fp.write(header_comments)
        fp.write(header)

def write_copynumber_segments_csv(df_hp1, df_hp2, haplotype_df_, args, centers, integer_fractional_means, hp, filename, p_value, is_half):

    haplotype_df = haplotype_df_.copy()

    #uniques = sorted(haplotype_df['state'].unique())
    #integer_fractional_means = sorted([i for i in range(0, len(uniques))])
    if not "subclonal" in filename:
        for i in range(len(integer_fractional_means)):
            haplotype_df['state'].mask(haplotype_df['state'] == centers[i], integer_fractional_means[i], inplace=True)
    #haplotype_df = mask_df_states(haplotype_df_copy, centers, integer_fractional_means)

    if 'p_value' in haplotype_df.columns:
        haplotype_df['if_subclonal'] = ['Y' if element < args.confidence_subclonal_score and not state == 0  else 'N' for (element,state) in zip(haplotype_df.p_value.values.tolist(), haplotype_df.state.values.tolist())]
        haplotype_df['state'] =  [subclonal_values_adjusted(element, centers) if sub == 'Y' else integers_values_adjusted(element, centers) for (element, sub) in zip(haplotype_df.state.values.tolist(), haplotype_df.if_subclonal.values.tolist())]
        haplotype_df = haplotype_df.rename(columns={'chromosome': 'chr', 'start': 'start', 'end': 'end', 'depth': 'coverage', 'state': 'copynumber_state', 'p_value': 'confidence', 'if_subclonal': 'if_subclonal'})
    else:
        if hp == 1:
            haplotype_df,_ = add_confidence_score_cn_segemnts(centers, haplotype_df, haplotype_df, df_hp1, df_hp2, args)
        else:
            _,haplotype_df = add_confidence_score_cn_segemnts(centers, haplotype_df, haplotype_df, df_hp1, df_hp2, args)

        haplotype_df = haplotype_df.rename(columns={'chromosome': 'chr', 'start': 'start', 'end': 'end', 'depth':'coverage', 'state':'copynumber_state', 'confidence_value': 'confidence'})

    haplotype_df['coverage'] = haplotype_df['coverage'].apply(lambda x: round(x, 2))
    if args.breakpoints:
        _, _, _, bps_ids_all, bps, bps_bnd = sv_vcf_bps_cn_check(args.breakpoints, args)
        bps_ids_global= []
        for index, row in haplotype_df.iterrows():
            bps_ids = []
            for bp in bps_ids_all:
                if bp[0] == row['chr'] and int(bp[1]) >= int(row['start']) and int(bp[1]) <= int(row['end']):
                    bps_ids.append(bp[2])
            bps_ids_global.append(list(set(bps_ids)))
        haplotype_df['svs_breakpoints_ids'] = bps_ids_global

    if args.without_phasing:
        fp = open(args.out_dir_plots + '/bed_output/' + args.genome_name + filename, 'a')
        fp.write('#chr: chromosome number\n')
        fp.write('#start: start address for CN segment\n')
        fp.write('#end: end address for CN segment\n')
        fp.write('#coverage: median coverage for this segment\n')
        fp.write('#copynumber_state: detected copy number state (integer/fraction)\n')

        if 'confidence' in haplotype_df.columns:
            fp.write('#confidence: confidence score\n')
            if "subclonal" in filename:
                fp.write('#if_subclonal: if entry is subclonal [Y/N]\n')
        if args.breakpoints:
            fp.write('#svs_breakpoints_ids: corresponding structural variations (breakpoints) IDs from VCF file\n')

        if 'confidence' in haplotype_df.columns:
            if args.breakpoints:
                if "subclonal" in filename:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tif_subclonal\tsvs_breakpoints_ids\n')
                else:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tsvs_breakpoints_ids\n')
            else:
                if "subclonal" in filename:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tif_subclonal\n')
                else:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\n')
        else:
            if args.breakpoints:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tsvs_breakpoints_ids\n')
            else:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\n')

        haplotype_df.to_csv(fp, sep='\t', index=False, mode='a', header=False)
    else:
        if is_half:
            if not os.path.isdir(args.out_dir_plots + '/wgd/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output'):
                os.makedirs(args.out_dir_plots + '/wgd/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output')
            fp = open(args.out_dir_plots +'/wgd/'+ str(args.tumor_ploidy) + '_'+ str(args.tumor_purity)  +'_'+ str(p_value) +'/bed_output/' + args.genome_name + filename, 'a')
        else:
            if not os.path.isdir(args.out_dir_plots + '/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output'):
                os.makedirs(args.out_dir_plots + '/' + str(args.tumor_ploidy) + '_'+ str(args.tumor_purity) +'_'+ str(p_value) +'/bed_output')
            fp = open(args.out_dir_plots +'/'+ str(args.tumor_ploidy) + '_'+ str(args.tumor_purity)  +'_'+ str(p_value) +'/bed_output/' + args.genome_name + filename, 'a')
        fp.write('#chr: chromosome number\n')
        fp.write('#start: start address for CN segment\n')
        fp.write('#end: end address for CN segment\n')
        fp.write('#coverage: median coverage for this segment\n')
        fp.write('#copynumber_state: detected copy number state (integer/fraction)\n')

        if 'confidence' in haplotype_df.columns:
            fp.write('#confidence: confidence score\n')
            if "subclonal" in filename:
                fp.write('#if_subclonal: if entry is subclonal [Y/N]\n')

        if args.breakpoints:
            fp.write('#svs_breakpoints_ids: corresponding structural variations (breakpoints) IDs from VCF file\n')

        if 'confidence' in haplotype_df.columns:
            if args.breakpoints:
                if "subclonal" in filename:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tif_subclonal\tsvs_breakpoints_ids\n')
                else:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tsvs_breakpoints_ids\n')
            else:
                if "subclonal" in filename:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\tif_subclonal\n')
                else:
                    fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tconfidence\n')
        else:
            if args.breakpoints:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\tsvs_breakpoints_ids\n')
            else:
                fp.write('#chr\tstart\tend\tcoverage\tcopynumber_state\n')

        haplotype_df.to_csv(fp, sep='\t', index=False, mode='a', header=False)

