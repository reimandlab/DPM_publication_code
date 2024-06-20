# DPM: Directional integration and pathway enrichment analysis for multi-omics data

## Supplementary Figure 2. Examples of directionally penalised pathways regulated by HOXA10-AS knockdown and overexpression

## Overview

We analysed the genes and pathways prioritised by directional integration of transcriptomics (RNA-seq) data from HOXA10-AS lncRNA knockdown (KD) and overexpression (OE) experiments in GBM cells from our earlier [study](https://pubmed.ncbi.nlm.nih.gov/34686327/)

## Input Data

``` markdown
./input_data/
    |- enriched_pathways.tsv # enriched pathways from the KD and OE experiments
    |- hoxa10_fc.csv # fold change values from the KD and OE experiments
    |- hoxa10_fdr.csv # fdr values from the KD and OE experiments
    |- hsapiens.GO_BP_REACTOME.name.gmt # gmt file of GO and Reactome terms
```

## Reproducibility

First run the scripts in the `Figure_3` directory, particularly `003_enrichment_map_panel_e.R` to obtain the enriched pathways. Then run the analysis sequentially: 

``` bash
python 001a-gene_pathway_dotplots.py
Rscript 001b-gene_pathway_dotplots.R

```





