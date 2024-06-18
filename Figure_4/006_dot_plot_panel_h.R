#### Figure 4H
### Dot plots of significant genes involved in directionally penalized pathways.
### These pathways are only enriched with Brown's method, but lost when using DPM.

library(ggplot2)
library(dplyr)

# load the enriched pathways from panel G
load("panelG_Rdata.Rdata")

# load known cancer genes from COSMIC Cancer Gene Census
cancer_genes <- read.table("../Figure_3/input_data/Census_allWed Jan 18 16_45_05 2023.tsv",sep="\t",header=TRUE)

# Order genes by merged P-value
rownames(df) <- df$gene
brown_scores <- data.frame(row.names = rownames(df),rna = df$rna_pval,protein = df$protein_pval)
brown_scores <- as.matrix(brown_scores)
brown_scores[is.na(brown_scores)] <- 1
brown_merged <- merge_p_values(brown_scores, "Brown")
ordered_brown_merged <- brown_merged[order(brown_merged)]

# Order pathways by FDR
res_brown <- res_brown[order(res_brown$adjusted_p_val),]
res_dpm <- res_dpm[order(res_dpm$adjusted_p_val),]

# find the pathways enriched by Brown's method and lost with DPM
res_brown <- res_brown[!res_brown$term_id %in% res_dpm$term_id,]

# create a dotplot for each pathway
pdf('dotplot_panelH.pdf',height = 1.8, width = 14)
for (i in 1:length(res_brown$term_id)){
      
      #order by the genes merged P-value
      gene_names <- unlist(res_brown[i,]$overlap)
      gene_names <- names(ordered_brown_merged)[names(ordered_brown_merged) %in% gene_names]
      
      gene_names_col <- rep(gene_names,2)
      dataset_contribution <- rep(c("rna","protein"),each=length(gene_names))

      pval_rna <- df[gene_names,]$rna_pval
      pval_protein <- df[gene_names,]$protein_pval
      pvals <- c(pval_rna,pval_protein)

      log2hr_rna <- df[gene_names,]$rna_log2hr
      log2hr_protein <- df[gene_names,]$protein_log2hr
      log2hr <- c(log2hr_rna,log2hr_protein)

      df_bubble <- data.frame(genes = gene_names_col,
                              dataset = dataset_contribution,
                              pval = as.numeric(pvals),
                              log2hr = as.numeric(log2hr))
      df_bubble$pval <- -log10(df_bubble$pval)

      # limit the color scale so the log2hr direction is clearer
      if (any(df_bubble[!is.na(df_bubble$log2hr),]$log2hr >= 2)){
            df_bubble[!is.na(df_bubble$log2hr) & df_bubble$log2hr >= 2,]$log2hr <- 2
      } 
      
      if (any(df_bubble[!is.na(df_bubble$log2hr),]$log2hr <= -2)){
            df_bubble[!is.na(df_bubble$log2hr) & df_bubble$log2hr <= -2, ]$log2hr <- -2      
      }
      
      df_bubble$dataset <- factor(df_bubble$dataset, levels=unique(df_bubble$dataset))
      df_bubble$genes <- factor(df_bubble$genes, levels=unique(df_bubble$genes))
      colnames(df_bubble) <- c("genes","dataset","-log10pval","log2HR")

      df_bubble$cancer <- "black"
      df_bubble$cancer[df_bubble$genes %in% cancer_genes$Gene.Symbol] <- "red"
      
      tit <- sym("-log10pval")
      print(ggplot() + theme_bw(base_size=7.5) +
                  geom_point(data = df_bubble, aes(x=genes, y=dataset, size= !!tit, fill=log2HR), color='black', pch=21)+
                  scale_size(limits = c(0.5,6.1),range = c(4,12))  +
                  scale_fill_gradient2(low="blue", high="red", mid = "white",limits = c(-2,2)) + 
                  xlab("") + ylab("") +
                  ggtitle(paste(res_brown$term_name[i],res_brown$term_id[i],
                         paste("Evidence",res_brown$evidence[i],sep=":"),
                         paste("Overlap Genes", length(gene_names),sep=":"),
                         paste("Size",res_brown$term_size[i],sep = ":"),
                         paste("FDR",formatC(res_brown$adjusted_p_val[i],format="e",digits = 2),
                               sep=":"),sep=" | ")) + 
                  theme(plot.title = element_text(hjust=1,size=8),
                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 10,colour=df_bubble$cancer),
                        legend.text = element_text(size = 12),
                        axis.text.y = element_text(size = 12),
                        axis.title = element_blank()))
}
dev.off()
