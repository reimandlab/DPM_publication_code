# DPM: Directional integration and pathway enrichment analysis for multi-omics data 
## Figure 3. Directional integration of transcriptomics data from functional experiments of HOXA10-AS lncRNA in GBM cells. 

## Overview

We analysed the genes and pathways prioritised by directional integration of transcriptomics (RNA-seq) data from HOXA10-AS lncRNA knockdown (KD) and overexpression (OE) experiments in GBM cells from our earlier [study](https://pubmed.ncbi.nlm.nih.gov/34686327/)

## Input Data

``` markdown
./input_data/
    |- 060821_hoxa10as_merged_oe_kd.rsav  # differential expression results
    |- CCDS_current.txt # dataset of protein-coding genes, downloaded January 16 2023
    |- 'Census_allWed Jan 18 16_45_05 2023.tsv' # cancer genes of the COSMIC Cancer Gene Census database, downloaded January 18 2023
    |- hsapiens.GO_BP_REACTOME.name.gmt # gmt file of GO biological processes and molecular pathways of Reactome, downloaded March 27 2023
```

## Reproducibility

Run the analysis sequentially:

``` bash
Rscript 001_gene_scatterplot_panel_b.R
Rscript 002a-dotplots_panel_c.py
Rscript 002b_dotplots_panel_c.R
Rscript 003_enrichment_map_panel_e.R
Rscript 004_venn_diagram_panel_d.R
Rscript 005_dot_plot_panel_f.R
```
