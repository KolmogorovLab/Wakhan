# SomaticCNA

###### Note: This repository is under extensive updates. Currently, only haplotagged BAM file coverage plots module is working.

Requirements:
* Python3
* Plotly
* Samtools
* Pysam
* Bcftools
* Tabix

## Installation
```
git clone https://github.com/KolmogorovLab/SomaticCNA.git
cd SomaticCNA/src
```

## Usage
```
python3 main.py --target-bam </home/rezkuh/GenData/COLO829/colo829_tumor_grch38_md_chr7:78318498-78486891_haplotagged.bam> --out-dir-plots <coverage_plots> --genome-name <Cellline> --phased-vcf <HCC1437BL.phased.vcf.gz>
```

## Required parameters

* `--target-bam` path to one or multiple target bam files (must be indexed)
  
* `--out-dir-plots` path to output coverage plots

* `--genome-name` genome cellline/sample name to be displayed on plots

* `--phased-vcf` phased VCF file for the corresponding haplotagged BAM
  
## Optional parameters
  
* `--phaseblock-flipping-enable` enabling phaseblock flipping in coverage plots
  
* `--smoothing-enable` enabling smoothing in coverage plots

* `--phaseblocks-enable` enabling phaseblocks display in coverage plots

* `--het-phased-snps-freq-enable` enabling hetrozygous phased snps frequencies in coverage plots

* `--breakpoints-enable` enabling breakpoints in coverage plots
  
