# DPM: Directional integration and pathway enrichment analysis for multi-omics data

## Supplementary Figure 3. Directional integration of genes and pathways jointly upregulated or downregulated in HOXA10-AS knockdown and overexpression experiments.

## Overview

We analysed the genes and pathways prioritised by directional integration of transcriptomics (RNA-seq) data from HOXA10-AS lncRNA knockdown (KD) and overexpression (OE) experiments in GBM cells from our earlier [study](https://pubmed.ncbi.nlm.nih.gov/34686327/)

## Input Data

``` markdown
./input_data/
    |- 060821_hoxa10as_merged_oe_kd.rsav  # differential expression results
    |- CCDS_current.txt # dataset of protein-coding genes, downloaded January 16 2023
    |- hsapiens.GO_BP_REACTOME.name.gmt # gmt file of GO biological processes and molecular pathways of Reactome, downloaded March 27 2023
```

## Reproducibility

``` bash
Rscript 001_enrichment_map.R
```
