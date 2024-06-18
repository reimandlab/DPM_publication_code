#### Figure 3B
### Compare the directionally prioritised to the directionally penalised genes using DPM vs Brown's method. 

library("ggplot2")
library("data.table")
library('dplyr')
library('ActivePathways')

# load the differential expression results from HOXA10-AS KD vs WT and
# HOXA10-AS overexpression vs WT
load("input_data/060821_hoxa10as_merged_oe_kd.rsav")
df_oe_kd <- df[,c(2,3,6,7,8,11,12)]

### lets analyze the inverse constraints_vector of c(-1,1)
scores <- data.frame(row.names = df_oe_kd$gene_name, overexpression = df_oe_kd$OE_P.Value, knockdown = df_oe_kd$KD_P.Value)
scores <- as.matrix(scores)
scores[is.na(scores)] <- 1

scores_dir <- data.frame(row.names = df_oe_kd$gene_name, overexpression = df_oe_kd$OE_logFC, knockdown = df_oe_kd$KD_logFC)
scores_dir <- as.matrix(scores_dir)
scores_dir[is.na(scores_dir)] <- 0
constraints_vector <- c(-1,1)

brown_merged <- merge_p_values(scores=scores, method="Brown")
dpm_merged <- merge_p_values(
      scores=scores,
      method="DPM",
      scores_direction=scores_dir,
      constraints_vector=constraints_vector
      )

# prepare a scatterplot where DPM results are on the y-axis and Brown on the x-axis
df_scatterplot <- data.frame(brown_results = -log10(brown_merged), dpm_results = -log10(dpm_merged))

pdf(file = paste("figure3_panel_b_scatterplot",".pdf",sep = ""),   
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