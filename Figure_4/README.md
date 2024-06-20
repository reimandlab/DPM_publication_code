# Figure 4 Analyses

## Overview

We integrated quantitative proteomics (isobaric label quantitation analysis with orbitrap) and transcriptomic (RNA-seq) data of cancer samples with patient survival information obtained 
from the CPTAC-3 and TCGA PanCanAtlas projects. This dataset included 1,140 cancer samples of ten cancer types: pancreatic, ovarian, colorectal, breast, kidney, head & neck, and endometrial cancers, two subtypes of lung cancer, and GBM. The main analysis focused on ovarian cancer (OV).

We used the combined dataset assembled by Zhang et al. (2022) that included transcriptomics data for 15,424 genes and proteomics data for ~10,000 genes that varied between cancer types. This dataset was preprocessed and represented as standard deviations from cohort median values. See [cancer-proteomics-compendium-n2002](https://github.com/chadcreighton/cancer-proteomics-compendium-n2002) to download the compendium data files.

## Input Data

``` markdown
./input_data/
    |- cptac_clinical # folder containing 7 cancer types with CPTAC patient ids
        |- clinical_and_survival_information.txt # CPTAC patient clinical and survival data
        |- *_normalized_mRNA.csv # normalized transcriptomics data
        |- *_normalized_total_protein.csv # normalized proteomics data
    |- tcga_clinical # folder containing 3 cancer types with TCGA patient ids
        |- clinical_and_survival_information.txt # TCGA patient clinical and survival data
        |- *_normalized_mRNA.csv # normalized transcriptomics data
        |- *_normalized_total_protein.csv # normalized proteomics data
    |- alias_corrected_names.csv # file used to convert gene alias names to full official names
    |- *_rename_numeric_genes.tsv # file used to convert numeric gene names to their full official names
```

## Reproducibility

Run the analysis sequentially:

``` bash
Rscript 000_survival_analysis.R
Rscript 001_gene_scatterplot_panel_b.R
Rscript 002_gene_log2hr_plot_panel_c.R
Rscript 003_correlation_scatterplots_panel_e.R
Rscript 004_enrichment_map_panel_g.R
Rscript 005_venn_diagram_panel_f.R
Rscript 006_dot_plot_panel_h.R
```

The Panel D script to create the Kaplan-Meier plots can be found in `000_survival_analysis.R` lines 148 to 160.