#### Figure 4B
### Scatter plot of merged P-values of OS associations from directional and non-directional analyses in OV

library("ggplot2")
library("data.table")
library("dplyr")
library("ActivePathways")

# load the survival analysis results
df <- read.table("OV/survival_analysis_results.csv", sep=",", header=TRUE)

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

# prepare a scatterplot where DPM results are on the y-axis and Brown on the x-axis
df_scatterplot <- data.frame(brown_results = -log10(brown_merged), dpm_results = -log10(dpm_merged))

pdf(file = paste("figure4_panel_b_scatterplot",".pdf",sep = ""),   
    width = 8, 
    height = 7)
print(ggplot(df_scatterplot) +
            geom_point(size = 2.4,shape = 19,aes(brown_results, dpm_results,color = ifelse(brown_results <= 1.301,"gray",ifelse(dpm_results > 1.301,"#1F449C","#F05039")))) +
            labs(title="", x ="Brown P (-log10)", y = "DPM P (-log10)") + 
            geom_hline(yintercept=1.301, linetype='dashed', col = 'black', size = 0.5)+
            geom_vline(xintercept = 1.301, linetype = "dashed", col = "black", size = 0.5) +
            theme(
                  plot.title = element_text(size=23,hjust = 0.5),
                  axis.title.x = element_text(size=18,margin = unit(c(2, 0, 0, 0), "mm")),
                  axis.title.y = element_text(size=18,margin = unit(c(0,4,0,0), "mm")),
                  axis.text = element_text(size = 16),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))+
            geom_abline(size = 0.5, slope=1,intercept = 0) +
            scale_color_identity())
dev.off()
