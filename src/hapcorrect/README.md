# hapcorrect

###### Note: This repository is under extensive updates.

hapcorrect takes long-read alignment and phased heterozygous variants as input, and corrects the haplotypes phase-switching errors in BAMs around phased blocks as well as inside phase-blocks.


#### Phasing errors:
<img width="1373" alt="plots_example" src="src/images/1.png">

#### Phasing errors correction:
<img width="1373" alt="plots_example" src="src/images/2.png">

## Installation (individual packages through conda and pip)
```
git clone https://github.com/KolmogorovLab/hapcorrect.git
cd hapcorrect/
conda create --name hapcorrect python=3.8
conda activate hapcorrect
pip install pysam pyfaidx numpy pandas Bio scipy ReportLab==3.6.12 matplotlib py plotly hmmlearn tqdm>=4.24.0 numba>=0.43.0 nose pomegranate==0.14.8 scikit-genome==0.0.1 scikit-learn==1.2.0 scipy==1.9.2  matplotlib==3.6.2 kneed ruptures
conda install -c bioconda samtools bcftools
cd src/
```

## Usage 
hapcorrect works with both tumor/normal pair and tumor-only data. 
To generate input for hapcorrect please follow the [following](https://github.com/KolmogorovLab/Wakhan/tree/main?tab=readme-ov-file#prerequisite) instructions.


### Tumor-Normal pair (requires tumor BAM and normal phased VCF)
```
python main.py --threads <4> --reference <ref.fa>  --target-bam <data.tumor.bam>  --normal-phased-vcf <data.normal_phased.vcf.gz>  --genome-name <cellline/dataset name> --cut-threshold <150> --out-dir-plots <genome_abc_output>
```
## Optional parameters
* `--tumor-vcf` If tumor VCF is also provided, it helps hapcorrect to detect LOH and enables better phase correction.
* `--rephase-normal-vcf` Set true to rephase normal VCF with hapcorrect phasing correction (default: disabled)
* `--rehaplotag-tumor-bam` Set true to rehaplotag the tumor BAM with new rephased VCF from above step (default: disabled)


### Tumor-only (requires tumor phased/haplotagged BAM and tumor phased VCF)
```
python main.py --threads <4> --reference <ref.fa>  --target-bam <data.tumor_haplotagged.bam>  --tumor-vcf <data.tumor_phased.vcf.gz>  --genome-name <cellline/dataset name> --cut-threshold <150> --out-dir-plots <genome_abc_output>
```
## Optional parameters
* `--rephase-tumor-vcf` Set true to rephase tumor VCF with hapcorrect phasing correction (default: disabled)

## Note:
In some cases, phaseblocks are not good enough for `hapcorrect` to correct phase-switch errors (too small phaseblocks), in that scenario, user can set `enable-simple-heuristics` parameter to apply simple heuristics which will assign higher coverage bins values to HP-1 and lower bins coverage values to HP-2.  