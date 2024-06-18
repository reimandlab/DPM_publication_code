#### Figure 4E
### Scatterplots of mRNA and protein expression of ACTN4 and PIK3R4 with Spearman
### correlation coefficients and P-values

library("ggplot2")
library("data.table")
library('dplyr')
library('ggpubr')

# load rna data
rna <- read.table("input_data/tcga_clinical/OV-CPTAC-TCGA_normalized_mRNA.csv", sep =",", header=TRUE)

# load protein data
protein <- read.table("input_data/tcga_clinical/OV-CPTAC-TCGA_normalized_total_protein.csv", sep=",", header=TRUE)

# lets create a scatterplot with correlation information for ACTN4 and PIK3R4 RNA vs Protein (for each patient)
scatterplot_function <- function(gene_name){
      # filter RNA and Protein data for the gene of interest
      rna_gene <- rna[rna$Gene == gene_name,]
      protein_gene <- protein[protein$Gene == gene_name,]
      
      # identify common patients between the two datasets
      rna_gene <- rna_gene[names(rna_gene) %in% names(protein_gene)]
      protein_gene <- protein_gene[names(protein_gene) %in% names(rna_gene)]
      
      df_scatterplot <- data.frame(rna_data = as.numeric(rna_gene),protein_data = as.numeric(protein_gene))
      df_scatterplot <- df_scatterplot[complete.cases(df_scatterplot),]
      
      pdf(paste(paste("figure4_panel_e_scatterplot",gene_name,sep = "_"),"pdf",sep = "."),width = 8, height = 7)
      print(ggscatter(df_scatterplot, x = "rna_data", y = "protein_data", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.coef.size = 5, cor.method = "spearman",
                      xlab = "RNA expression (Z-score)",
                      ylab = "Protein expression (Z-score)", title = gene_name) + 
                  theme_bw(base_size=11) +
                  theme(
                        plot.title = element_text(size=20,hjust = 0.5),
                        axis.title.x = element_text(size=18,margin = unit(c(2, 0, 0, 0), "mm")),
                        axis.title.y = element_text(size=18,margin = unit(c(0,4,0,0), "mm")),
                        axis.text = element_text(size = 18),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.text = element_text(size=24))
      )
      dev.off()
}
scatterplot_function("ACTN4")
scatterplot_function("PIK3R4")
