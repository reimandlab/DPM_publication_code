# DPM: Directional integration and pathway enrichment analysis for multi-omics data

## Source code and input files

This repository contains the source code and most input files required to generate results and figures of the manuscript:

Directional integration and pathway enrichment analysis for multi-omics data.

Mykhaylo Slobodyanyuk \^, Alexander T. Bahcheli \^, Zoe P. Klein, Masroor Bayati, Lisa J. Strug, J�ri Reimand.

Nature Communications (2024) (in press; June 2024).

\^ - co-first authors; \@ - correspondence: [juri.reimand\@utoronto.ca](mailto:juri.reimand@utoronto.ca)

The materials are organized into sub-directories for each Figure that includes source code and input datasets. The R package ActivePathways v2.0.3 was used. A Zenodo archive of the R package is available at DOI [10.5281/zenodo.12118089](https://doi.org/10.5281/zenodo.12118089)

A few input datasets in Figure 5 require access from the GLASS consortium and can be retrieved from the Synapse database at <https://www.synapse.org> (accession numbers: mRNA: syn31121291; clinical: syn31121219; methylation: syn23594913). A few very large datasets from TCGA are not provided as part of this repository since these are publicly available and their duplication is unnecessary. Download instructions for these files are shown within individual files.

### Source code and package versions

The source code is in R (version 4.3.1) and python (version 3.10.13). 
The following R package versions were used:


``` bash
ggpubr v0.6.0
survival v3.5-5
eulerr v7.0.2
ActivePathways v2.0.3
dplyr v1.1.2         
data.table v1.14.8
ggplot2 v3.4.2
optparse v1.7.3
forcats v1.0.0
gridExtra v2.3
ggrepel v0.9.5
gplots v3.1.3.1
RColorBrewer v1.1.3
gtools v3.9.5
```


The following python package versions were used:


``` bash
pandas v2.2.0
numpy v1.26.4
scipy v1.12.0
statsmodels v0.14.1
```

