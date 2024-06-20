#### Supplementary Figure 4A
### Scatter plot of merged P-values of OS associations from directional and non-directional analyses

library("ggplot2")
library("data.table")
library("dplyr")
library("ggrastr")
library("ActivePathways")

# integrate results for the 9 remaining cancer types
cancer_types <- c("GBM","HNSC","LUAD","LUSC","PAAD","Renal","UCEC","BREAST","CRC") 
df_scatterplot <- data.frame(brown_results = integer(0), dpm_results = integer(0),cancer = character(0))
for (type in cancer_types){
      # load the survival analysis results
      df <- read.table(paste0("../Figure_4/",type,"/survival_analysis_results.csv"), sep=",", header=TRUE)
      
      scores <- data.frame(row.names = df$gene, rna = df$rna_pval, protein = df$protein_pval)
      scores <- as.matrix(scores)
      scores[is.na(scores)] <- 1
      scores_dir <- data.frame(row.names = df$gene, rna = df$rna_log2hr, protein = df$protein_log2hr)
      scores_dir <- as.matrix(scores_dir)
      scores_dir[is.na(scores_dir)] <- 0
      constraints_vector <- c(1,1)
      
      brown_merged <- merge_p_values(scores=scores, method="Brown")
      dpm_merged <- merge_p_values(
            scores=scores,
            method="DPM",
            scores_direction=scores_dir,
            constraints_vector=constraints_vector
      )
      
      # add the DPM results and Brown results for each cancer type to a dataframe for plotting
      df_integrated <- data.frame(brown_results = -log10(brown_merged), dpm_results = -log10(dpm_merged))
      df_integrated$cancer <- type
      df_scatterplot <- rbind(df_scatterplot, df_integrated)
}


df_scatterplot = df_scatterplot %>%
  mutate(cancer = factor(cancer))

# remove one gene outlier from UCEC for plotting
df_scatterplot[df_scatterplot$cancer == "UCEC" & df_scatterplot $brown > 50,]
df_scatterplot <- df_scatterplot[rownames(df_scatterplot) != "ZNF348",]

# plot
pdf(file = "supplementary_fig4_panel_a_scatterplot.pdf",   
    width = 8, 
    height = 7)
g_plot <- ggplot(df_scatterplot) + 
      geom_point(size = 2.4,shape = 19,
                 aes(brown_results, dpm_results,
                     color = ifelse(brown_results <= 1.301,"gray",
                                    ifelse(dpm_results > 1.301,"#1F449C","#F05039")))) +
      labs(title="", x ="Brown P (-log10)", y = "DPM P (-log10)") + 
      geom_hline(yintercept=1.301, linetype='dashed', col = 'black', size = 0.5)+
      geom_vline(xintercept = 1.301, linetype = "dashed", col = "black", size = 0.5) +
      facet_wrap( ~ cancer, ncol = 3, scales = "free") +
      theme(plot.title = element_text(size=23,hjust = 0.5),
            axis.title.x = element_text(size=18,margin = unit(c(2, 0, 0, 0), "mm")),
            axis.title.y = element_text(size=18,margin = unit(c(0,4,0,0), "mm")),
            axis.text = element_text(size = 16),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))+
      geom_abline(size = 0.5, slope=1,intercept = 0) +
      scale_color_identity()
print(g_plot)
dev.off()
