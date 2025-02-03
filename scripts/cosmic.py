#!/usr/bin/env python3

import pandas as pd

#Complete COSMIC academic/research purpose cancer census genes set (Cosmic_CancerGeneCensus_v101_GRCh38.tsv) could be downloaded from [COSMIC](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v101/cancergenecensus).

df = pd.read_csv("Cosmic_CancerGeneCensus_Tsv_v101_GRCh38/Cosmic_CancerGeneCensus_v101_GRCh38.tsv", sep="\t", header=0)
df = df[['CHROMOSOME','GENOME_START','GENOME_STOP','GENE_SYMBOL']]
df['CHROMOSOME'] = 'chr' + df['CHROMOSOME']
df.to_csv('cosmic_genes.tsv', sep='\t', index=False, header=False)