This sub-directory includes some cell lines sample phasing with switching errors and SomaticCNA updated phasing correction.

Default directory for each cell line contains coverage plots for each chromosome output direct from haplotagged BAMs while arbitrary directory contains coverage plots after fixing phase-switch errors through this tool.

Bin size for these cell lines coverage plots is set 50K, which can be adjusted as per requirements through `--bin-size` parameter.

DASHBOARD.html provides a whole genome view for each chromosome. 

Unphased reads coverage, phaseblocks for both haplotypes, smoothing could be enabled through `--unphased-reads-coverage-enable`, `--phaseblocks-enable`, `--smoothing-enable` parameters.

Phase flipping to correct phaseblocks/switch error for haplotypes is enabled throguh `--phaseblock-flipping-enable`. 