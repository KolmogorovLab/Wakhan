import statistics


def infer_missing_phaseblocks(ref_start_values, ref_end_values, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets,
                              haplotype_2_values_phasesets, haplotype_1_values, haplotype_2_values, bin_size):
    missing_segments_starts = []
    missing_segments_ends = []
    missing_segments_hp1_value = []
    missing_segments_hp2_value = []
    for i in range(len(ref_start_values_phasesets)-1):
        if ref_start_values_phasesets[i +1] - ref_end_values_phasesets[i] > bin_size * 6:
            start = ((ref_end_values_phasesets[i]//bin_size) +1) * bin_size + 1
            end = (ref_start_values_phasesets[i+1]//bin_size) * bin_size
            missing_segments_starts.append(start)#(ref_end_values_phasesets[i] + 1)
            missing_segments_ends.append(end)#(ref_start_values_phasesets[i+1] - 1)
            try:
                missing_segments_hp1_value.append(statistics.mean(list(filter(lambda zer: zer != 0, haplotype_1_values[ref_start_values.index(start):ref_end_values.index(end)]))))
            except:
                missing_segments_hp1_value.append(0)
            try:
                missing_segments_hp2_value.append(statistics.mean(list(filter(lambda zer: zer != 0, haplotype_2_values[ref_start_values.index(start):ref_end_values.index(end)]))))
            except:
                missing_segments_hp2_value.append(0)

    if len(missing_segments_starts):
        ref_start_values_phasesets = ref_start_values_phasesets + missing_segments_starts
        ref_end_values_phasesets = ref_end_values_phasesets + missing_segments_ends
        haplotype_1_values_phasesets = haplotype_1_values_phasesets + missing_segments_hp1_value
        haplotype_2_values_phasesets = haplotype_2_values_phasesets + missing_segments_hp2_value

        sort_function = lambda x: x[0]
        sort_target = list(zip(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets))
        sort_target.sort(key=sort_function)

        ref_start_values_phasesets = [a for a, b, c, d in sort_target]
        ref_end_values_phasesets = [b for a, b, c, d in sort_target]
        haplotype_1_values_phasesets = [c for a, b, c, d in sort_target]
        haplotype_2_values_phasesets = [d for a, b, c, d in sort_target]

    return ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets


def is_phasesets_check_simple_heuristics(ref_start_values_phasesets, ref_end_values_phasesets, args):
    ps_region_starts = []
    ps_region_ends = []
    for i, (start, end) in enumerate(zip(ref_start_values_phasesets, ref_end_values_phasesets)):
        if end - start > 2000000:
            ps_region_starts.append(start)
            ps_region_ends.append(end)

    if len(ps_region_starts) < 5:
        return True
    else:
        return False

