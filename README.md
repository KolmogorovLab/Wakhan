# Wakhan

###### Note: This repository is under extensive updates.

A tool to analyze haplotype-specific chromosome-scale somatic copy number aberrations and aneuploidy using long reads (Oxford Nanopore, PacBio). 
Wakhan takes long-read alignment and phased heterozygous variants as input, and first uses extends the phased blocks, taking
advantage of the CNA differences between the haplotypes. Wakhan then generates inetractive haplotype-specific coverage plots.    

#### Breakpoints/SVs based segmentation and Copy numbers estimation:
<img width="1373" alt="plots_example" src="examples/images/1437.png">

[//]: # (#### LOH detection, phasing correction and CopyNumbers profiling &#40;COLO357&#41;:)

[//]: # (<p float="left">)

[//]: # (  <img src="examples/images/loh1.png" width="270" />)

[//]: # (  <img src="examples/images/loh2.png" width="270" /> )

[//]: # (  <img src="examples/images/loh3.png" width="270" />)

[//]: # (</p>)

[//]: # (#### Copy number segmentation:)

[//]: # (<img width="1373" alt="plots_example" src="examples/images/HG008_CN.png">)

[//]: # ()
[//]: # (#### LOH detection)

[//]: # (<img width="1373" alt="plots_example" src="examples/images/HG008_LOH.png">)

[//]: # (## Installation &#40;individual packages through conda and pip&#41;)

[//]: # (```)

[//]: # (git clone https://github.com/KolmogorovLab/Wakhan.git)

[//]: # (cd Wakhan/)

[//]: # (conda create --name Wakhan python=3.8)

[//]: # (conda activate Wakhan)

[//]: # (pip install pysam pyfaidx numpy pandas plotly scikit-learn==1.2.0 scipy==1.9.2 ruptures vcf_parser)

[//]: # (conda install -c bioconda samtools bcftools)

[//]: # (cd src/)

[//]: # (```)

## Installation (enabling through conda environment)
```
git clone https://github.com/KolmogorovLab/Wakhan.git
cd Wakhan/
conda env create -f environment.yml -n Wakhan
conda activate Wakhan
cd src/
```

## Usage

### Tumor-Normal Mode (requires normal phased VCF)
```
python main.py --threads <4> --reference <ref.fa>  --target-bam <data.tumor.bam>  --normal-phased-vcf <data.normal_phased.vcf.gz>  --genome-name <cellline/dataset name> --out-dir-plots <genome_abc_output> --breakpoints <severus-sv-VCF>
```

### Tumor-only (requires tumor phased/haplotagged BAM and phased VCF)
```
python main.py --threads <4> --reference <ref.fa>  --target-bam <data.tumor_haplotagged.bam>  --tumor-vcf <data.tumor_phased.vcf.gz> --genome-name <cellline/dataset name> --out-dir-plots <genome_abc_output> --breakpoints <severus-sv-VCF>
```

##### Breakpoints/Structural variations or change point detection algo for copy number model

Wakhan accepts [Severus](https://github.com/KolmogorovLab/Severus) or any other structural variant caller VCF as breakpoints with param `--breakpoints` inputs to detect copy number changes and this option is highly recommended. 
However, if `--breakpoints` option is not used, `--change-point-detection-for-cna` should be used instead to use change point detection algorithm [ruptures](https://centre-borelli.github.io/ruptures-docs/) alternatively.

##### Tumor-Normal mixture and purity/ploidy estimation

User can input both `--ploidy-range` [default: 1.5-5.5 -> [min-max]] and `--purity-range` [default: 0.5-1.0 -> [min-max] to inform copy number model about normal contamination in tumor to estimate copy number states correctly.

##### Genes/copy number annotations

By default, Wakhan uses [COSMIC](https://cancer.sanger.ac.uk/cosmic) cancer genes to display corresponding copy number states in `<genome_name>_<ploidy>_<purity>_<confidence>_genes_genome.html` file. 
User can also input path through param `--user-input-genes` to custom input genes/subset of genes bed file these genes will be used in plots instead of default COSMIC cancer genes.
grch38 reference genes will be use as default, user can input alternate (i.e, chm13) `--reference-name` to change to T2T-CHM13 instead. 

##### Quick-run if coverage/pileup data is already available

Wakhan produces reads coverage `coverage.csv` (bin-size based reads coverage) and phasesets reads coverage `coverage_ps.csv` data, phase-corrected coverage `phase_corrected_coverage.csv` (as well as tumor BAM pileup `pileup_SNPs.csv` in case Tumor/normal mode) and stores in directory `coverage_data` inside `--out-dir-plots` location.
If this data has already been generated in a previous Wakhan run, user can rerun the Wakhan with additionally passing `--quick-start` and `--quick-start-coverage-path <path to coverage_data directory -> e.g., /home/rezkuh/data/1437/coverage_data>` in addition to required params in above example runs.
This will save runtime significantly by not invoking coverage and pileup methods. 

[//]: # (## Note)

[//]: # (If for some reason you have already generated coverage and pileup data from Wakhan &#40;as pileup takes some time&#41; and want to rerun the tool, you can avoid generating coverage/pileup data again by copying this data and using it again:)

[//]: # (1. Copy `coverage.csv`, `coverage_ps.csv` and `<Your genome name>_SNPs.csv` files from `data/` output dir &#40;it should be in src&#41; of your current run to some separate directory i.e.,  `/home/abc/dry_run_data`.)

[//]: # (2. Then run again the tool with adding this command additional to what you use already: `--dryrun True` `--dryrun-path` <This is the path where you copied CSVs files in step-1, ie. like, `/home/abc/dry_run_data/`>)

## Examples
Few cell lines arbitrary phase-switch correction and copy number estimation output with coverage profile is included in the [examples](https://github.com/KolmogorovLab/Wakhan/tree/main/examples) directory. 

## Required parameters
* `--reference` Reference file path

* `--target-bam` path to target bam files (must be indexed)
  
* `--out-dir-plots` path to output coverage plots

* `--genome-name` genome cellline/sample name to be displayed on plots

* `--normal-phased-vcf` normal phased VCF file to generate het SNPs frequncies pileup for tumor BAM (if tumor-only mode, use phased `--tumor-vcf` instead)

* `--tumor-vcf` phased VCF is required in tumor-only mode


## Optional parameters
* `--breakpoints` For segmentation to use in CN estimation, structural variations/breakpoints VCF file is required

* `--breakpoints-min-length` To adjust breakpoints min length to be included in copy number analysis [default: 10000]

* `--cpd-internal-segments` For change point detection algo on internal segments after breakpoint/cpd algo for more precise segmentation.

* `--copynumbers-subclonal-enable` Enabling subclonal/fractional copy number states in plots

* `--loh-enable` Enabling LOH regions display in CN plots

* `--phaseblock-flipping-disable` disabling phaseblock flipping if traget tumor BAM doesn't need phase-correction (default: enabled)

* `--phaseblocks-enable` enabling phaseblocks display in coverage plots

* `--contigs` List of contigs (chromosomes, default: chr1-22) to be included in the plots [e.g., chr1-22,X,Y]

* `--cut-threshold` Maximum cut threshold for coverage (readdepth) plots [default: 100]

* `--centromere` Path to centromere annotations BED file [default: annotations/grch38.cen_coord.curated.bed]

* `--cancer-genes` Path to Cancer Genes in TSV format to display in CNA plots [default: annotations/CancerGenes.tsv]

* `--pdf-enable` Enabling PDF output for plots

Wakhan can also be used in case phasing is not good in input tumor or analysis is being performed without considering phasing:

* `--without-phasing` enable it if CNA analysis is being performed without phasing in conjunction with `--phaseblock-flipping-disable` with all other required parameters as mentioned in example command

Here is a sample copy number/breakpoints output plot without phasing.
<img width="1373" alt="plots_example" src="examples/images/C15.png">

## Output produced
Based on best confidence scores, tumor purity and ploidy values are calculated and copy number analysis is performed. 
Each subfolder in output directory represents best <`ploidy`>_<`purity`>_<`confidence`> values.

* `<genome-name>_genome_copynumber.html` Genome-wide copy number plots with coverage information on same axis
* `<genome-name>_copynumber_breakpoints.html` Genome-wide copy number plots with coverage information on opposite axis, additionally breakpoints and genes annotations 
* `<genome-name>_copynumber_breakpoints_subclonal.html` Genome-wide subclonal/fractional copy number plots with coverage information on opposite axis, additionally breakpoints and genes annotations (`--copynumbers-subclonal-enable`)
* `bed_output` It contains copy numbers segments in bed format
* `variation_plots` Copy number chromosomes-scale plots with segmentation, coverage and LOH

Following are coverage and SNPs/LOH plots and bed directories in output folder, independent of CNA analysis

* `snps_loh_plots` SNPs and SNPs ratios plots with LOH representation in chromosomes-scale and genome-wide
* `<genome-name>_genome_loh.html` Genome-wide LOH plot
* `bed_output` It contains LOH segments in bed format
* `coverage_plots` Haplotype specific coverage plots for chromosomes with option for unphased coverage
* `phasing_output` Phase-switch error correction plots and phase corrected VCF file (*rephased.vcf.gz)


## Prerequisite
This tool requires haplotagged tumor BAM and phased VCF in case tumor-only mode and normal phased VCF in case tumor-normal mode. This can be done through any phasing tools like Margin, Whatshap and Longphase. 
Following commands could be helpful for phasing VCFs and haplotagging BAMs.

#### For normal/tumor pair:
```
# ClairS phase and haplotag both normal and tumor samples
singularity run clairs_latest.sif /opt/bin/run_clairs --threads 56 --phase_tumor True --use_whatshap_for_final_output_haplotagging --use_whatshap_for_final_output_phasing --tumor_bam_fn normal.bam --normal_bam_fn tumor.bam --ref ref.fasta --output_dir clairS --platform ont_r10
```
or
```
# Phase normal sample
pepper_margin_deepvariant call_variant -b normal.bam -f ref.fasta -o pepper/output -t 56 --ont_r9_guppy5_sup -p pepper --phased_output

# Haplotag tumor sample with normal phased VCF (phased.vcf.gz) output from previous step
whatshap haplotag --ignore-read-groups phased.vcf.gz tumor.bam  --reference ref.fasta -o tumor_whatshap_haplotagged.bam
```
#### For tumor only:
```
# Phase and haplotag tumor sample
singularity run clair3_latest.sif /opt/bin/run_clair3.sh --use_whatshap_for_final_output_haplotagging --use_whatshap_for_final_output_phasing --bam_fn=tumor.bam --ref_fn=ref.fasta --threads=56 --platform=ont --model_path=r941_prom_sup_g5014 --output=clair3 --enable_phasing
```
or
```
# Phase and haplotag tumor sample
pepper_margin_deepvariant call_variant -b tumor.bam -f ref.fasta -o pepper/output -t 56 --ont_r9_guppy5_sup -p pepper --phased_output
```