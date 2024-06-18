#### Figure 3F
### Dot plots of significant genes involved in directionally penalized pathways.
### These pathways are only enriched with Brown's method, but lost when using DPM.

library(ggplot2)
library(dplyr)

# load the enriched pathways from panel E
load("input_data/panelE_Rdata.Rdata")

# load known cancer genes from COSMIC Cancer Gene Census
cancer_genes <- read.table("input_data/Census_allWed Jan 18 16_45_05 2023.tsv",sep="\t",header=TRUE)

# Order genes by merged P-value
rownames(df_oe_kd) <- df_oe_kd$gene_name
brown_scores <- data.frame(row.names = rownames(df_oe_kd),overexpression = df_oe_kd$OE_P.Value,knockdown = df_oe_kd$KD_P.Value)
brown_scores <- as.matrix(brown_scores)
brown_scores[is.na(brown_scores)] <- 1
brown_merged <- merge_p_values(brown_scores, "Brown")
ordered_brown_merged <- brown_merged[order(brown_merged)]

# Order pathways by FWER
res_brown <- res_brown[order(res_brown$adjusted_p_val),]
res_dpm <- res_dpm[order(res_dpm$adjusted_p_val),]

# find the pathways enriched by Brown's method and lost with DPM
res_brown <- res_brown[!res_brown$term_id %in% res_dpm$term_id,]

# create a dotplot for each pathway
pdf('dotplot_panelF.pdf',height = 1.8, width = 14)
for (i in 1:length(res_brown$term_id)){
      # order by the genes merged P-value
      gene_names <- unlist(res_brown[i,]$overlap)
      gene_names <- names(ordered_brown_merged)[names(ordered_brown_merged) %in% gene_names]
      
      gene_names_col <- rep(gene_names,2)
      dataset_contribution <- rep(c("Overexpression","Knockdown"),each=length(gene_names))

      pval_oe <- df_oe_kd[gene_names,]$OE_P.Value
      pval_kd <- df_oe_kd[gene_names,]$KD_P.Value
      pvals <- c(pval_oe,pval_kd)
      
      fc_oe <- df_oe_kd[gene_names,]$OE_logFC
      fc_kd <- df_oe_kd[gene_names,]$KD_logFC
      fcs <- c(fc_oe,fc_kd)
      
      df_bubble <- data.frame(genes = gene_names_col,dataset = dataset_contribution,pval = pvals, log2fc = fcs)
      df_bubble$pval <- -log10(df_bubble$pval)
      print(df_bubble)
      
      # limit the color scale so the log2fc direction is clearer
      if (any(df_bubble[!is.na(df_bubble$log2fc),]$log2fc >= 2)){
            df_bubble[!is.na(df_bubble$log2fc) & df_bubble$log2fc >= 2,]$log2fc <- 2
      } 
      
      if (any(df_bubble[!is.na(df_bubble$log2fc),]$log2fc <= -2)){
            df_bubble[!is.na(df_bubble$log2fc) & df_bubble$log2fc <= -2, ]$log2fc <- -2      
      }
      
      df_bubble$dataset <- factor(df_bubble$dataset, levels=unique(df_bubble$dataset))
      df_bubble$genes <- factor(df_bubble$genes, levels=unique(df_bubble$genes))
      colnames(df_bubble) <- c("genes","dataset","-log10pval","log2fc")

      df_bubble$cancer <- "black"
      df_bubble$cancer[df_bubble$genes %in% cancer_genes$Gene.Symbol] <- "red"
     
      tit <- sym("-log10pval")
      print(ggplot() + theme_bw(base_size=7.5) +
            geom_point(data = df_bubble, aes(x=genes, y=dataset, size= !!tit, fill=log2fc), color='black', pch=21)+
            scale_size(limits = c(1,9),range = c(4,12))  +
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

