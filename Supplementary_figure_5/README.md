# DPM: Directional integration and pathway enrichment analysis for multi-omics data

## Supplementary Figure 5. Directional analysis of fibroblast growth factor receptor (FGFR) pathway in IDH mutant gliomas.

## Overview

We integrated quantitative proteomics (isobaric label quantitation analysis with orbitrap), transcriptomic (RNA-seq), and DNA methylation data of glioblastoma samples from the CPTAC-3 and TCGA PanCanAtlas projects.

We used the proteomics dataset assembled by Zhang et al. (2022) that included proteomics data for ~10,000 genes. This dataset was preprocessed and represented as standard deviations from cohort median values. See [cancer-proteomics-compendium-n2002](https://github.com/chadcreighton/cancer-proteomics-compendium-n2002) to download the compendium data files.

The transcriptomics and DNA methylation data were obtained using [TCGA biolinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html). The DNA methylation data was represented as beta values. Methlyation probes were associated with gene promoters using the [Improved DNA Methylation Array Probe Annotation](https://gdc.cancer.gov/content/improved-dna-methylation-array-probe-annotation) from the GDC. 

## Input Data

``` markdown
./input_data/
    |- intermediate_files # directory containing files that are generated during the analysis
    |- figure5_fcs.tsv # file of fcs from the different data sources
    |- figure5_pvals.tsv # file of pvalues from the different data sources
    |- gbm_tpm.tsv # tpm expression data for the GBM samples
    |- GBM-CPTAC_normalized_total_protein.csv # protein expression data 
    |- hsapiens.GO_BP_REACTOME.name.gmt # gmt file of GO and Reactome terms
    |- README_EPICv2.hg38.manifest.gencode.v41.txt # file describing how to download the EPIC probe associations with gene promoters
    |- README_tcga_gbm_methylation.txt # file describing how to download the TCGA GBM methylation data


```

## Reproducibility

Run the analysis sequentially:

``` bash
python 001a-gene_pathway_dotplots.py
python 002a-tcga_expression_boxplots.py
python 002c-tcga_methylation_boxplots.py
python 002e-cptac_protein_boxplots.py

```


