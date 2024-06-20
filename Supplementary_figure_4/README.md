# DPM: Directional integration and pathway enrichment analysis for multi-omics data

## Supplementary Figure 4. Integrating transcriptomic and proteomic datasets with cancer patient survival in 10 cancer types to find prognostic biomarkers and enriched pathways.

## Overview

We integrated quantitative proteomics (isobaric label quantitation analysis with orbitrap) and transcriptomic (RNA-seq) data of cancer samples with patient survival information obtained 
from the CPTAC-3 and TCGA PanCanAtlas projects. This dataset included 1,140 cancer samples of ten cancer types: pancreatic, ovarian, colorectal, breast, kidney, head & neck, and endometrial cancers, two subtypes of lung cancer, and GBM.

We used the combined dataset assembled by Zhang et al. (2022) that included transcriptomics data for 15,424 genes and proteomics data for ~10,000 genes that varied between cancer types. This dataset was preprocessed and represented as standard deviations from cohort median values. See [cancer-proteomics-compendium-n2002](https://github.com/chadcreighton/cancer-proteomics-compendium-n2002) to download the compendium data files.

## Input Data

First run the `Figure_4` script `000_survival_analysis.R`, which will produce a `survival_analysis_results.csv` file for each cancer type. This file is used as input for this figure.

## Reproducibility

First run the `Figure_4` script `000_survival_analysis.R`and then run the below scripts sequentially:

``` bash
Rscript 001_gene_scatterplots_panel_a.R
Rscript 002_venn_diagrams_panel_b.R
```

