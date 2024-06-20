# DPM: Directional integration and pathway enrichment analysis for multi-omics data

## Figure 5. Integrating transcriptomic, proteomic, and DNA methylation profiles of IDH-mutant gliomas.

## Overview

We integrated quantitative proteomics (isobaric label quantitation analysis with orbitrap), transcriptomic (RNA-seq), and DNA methylation data of glioblastoma samples from the CPTAC-3 and TCGA PanCanAtlas projects.

We used the proteomics dataset assembled by Zhang et al. (2022) that included proteomics data for ~10,000 genes. This dataset was preprocessed and represented as standard deviations from cohort median values. See [cancer-proteomics-compendium-n2002](https://github.com/chadcreighton/cancer-proteomics-compendium-n2002) to download the compendium data files.

The transcriptomics and DNA methylation data were obtained using [TCGA biolinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html). The DNA methylation data was represented as beta values. Methlyation probes were associated with gene promoters using the [Improved DNA Methylation Array Probe Annotation](https://gdc.cancer.gov/content/improved-dna-methylation-array-probe-annotation) from the GDC. 

## Input Data

``` markdown
./input_data/
    |- intermediate_files # directory containing files that are generated during the analysis
        |- figure5_fcs.tsv # file of FCs from each data source
        |- figure5_merged_pvalues.tsv # file of p-values merged using DPM or Brown's method
        |- figure5_pvals.tsv # file of p-values from each data source
        |- protein_degs.tsv # differential protein expression results
        |- tcga_gbm_degs.tsv # differential gene expression results
        |- tcga_gbm_methylation_degs.tsv # differential gene promoter methylation results
    |- cgc_v98_29062023.tsv # file of CGC genes
    |- GBM-CPTAC_normalized_total_protein.csv # file of normalized protein expression data
    |- hgnc_gene_annotations.tsv # file describing protein-coding genes
    |- hsapiens.GO_BP_REACTOME.name.gmt # gmt file of GO and Reactome terms
    |- README_EPICv2.hg38.manifest.gencode.v41.txt # file describing how to download the EPIC probe associations with gene promoters
    |- README_tcga_gbm_methylation.txt # file describing how to download the TCGA GBM methylation data
    |- tcga_gbm_counts.tsv # file of TCGA GBM RNA-seq data (forward strand raw STAR counts)


```

## Reproducibility

Run the analysis sequentially, with the exception of running "006-panel_f-enrichment_map.R" before the "005" scripts:

``` bash
python 000a-gbm_rna_dgea.py
python 000b-gbm_diff_methlyation.py
Rscript 001a-pvalue_merging.R
python 002a-panel_b-gene_venn_diagrams.py
python 003a-panel_c-genes_lineplot.py
python 004a-panel_d-heatmap.py
Rscript 006-panel_f-enrichment_map.R
python 005a-panel_e-common_pathways_venn_diagram.py
python 007a-panel_g-gene_pathway_dotplots.py
Rscript 008-panel_h-themes_plot.R

```


