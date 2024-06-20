# DPM: Directional integration and pathway enrichment analysis for multi-omics data

## Supplementary Figure 6. Validating the integrative analysis of IDH-mutant gliomas in an independent dataset of cancer samples.

## Overview

We integrated quantitative proteomics (isobaric label quantitation analysis with orbitrap), transcriptomic (RNA-seq), and DNA methylation data of glioblastoma samples from the publication Oh et al. (2022) and the Glioma longitudinal analysis (GLASS) cohort. The proteomics dataset was downloaded from https://proteomecentral.proteomexchange.org/ accession number PXD015545.

Input datasets from the GLASS project require additional approval by GLASS and can be downloaded from the Synapse database at https://www.synapse.org:
- TPM gene expression data: syn31121291
- Methylationd data: syn23594913
- Clinical data: syn31121219

Methlyation probes were associated with gene promoters using the [Improved DNA Methylation Array Probe Annotation](https://gdc.cancer.gov/content/improved-dna-methylation-array-probe-annotation) from the GDC.


## Input Data

``` markdown
./input_data/
    |- intermediate_files # directory containing files that are generated during the analysis
        |- glass_dgea_methylation.tsv # differentially methylated gene promoter results
        |- glass_dgea.tsv # file of differentially expressed gene results
        |- glass_expression.tsv # file of raw expression values for each gene and sample
        |- glass_fcs.tsv # fcs of each gene from each data source
        |- glass_pvals.tsv # pvalues of each gene from each data source
        |- merged_pvalues.tsv # pvalues merged using dpm and brown's method
        |- validation_proteomics_dgea.tsv # file of differentially expressed proteins from the validation dataset
    |- hgnc_gene_annotations.tsv # file describing protein-coding genes
    |- hsapiens.GO_BP_REACTOME.name.gmt # gmt file of GO and Reactome terms
    |- README_EPICv2.hg38.manifest.gencode.v41.txt # file describing how to download the EPIC probe associations with gene promoters
    |- README_glass_clinical.txt # file describing how to download the GLASS clinical data
    |- README_glass_expression.txt # file describing how to download the GLASS gene expression data
    |- README_glass_methylation.txt # file describing how to download the GLASS methylation data
    |- validation_proteomics_clinical.tsv # file of clinical data for the validation proteomics dataset
    |- validation_proteomics_expression.tsv # file of proteomics expression data from the validation dataset


```

## Reproducibility

Run the analysis sequentially:

``` bash
python 000a-glass_dgea.py
python 000b-glass_methylation.py
python 000c-proteomics_dpea.py
Rscript 001-pvalue_merging.R
python 001a-glass_venn_diagrams.py
python 002a-genes_lineplot.py
Rscript 003-enrichment_map.R

```


