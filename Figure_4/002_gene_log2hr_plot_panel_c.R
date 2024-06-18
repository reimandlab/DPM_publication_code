#### Figure 4C
### log-transformed HR values of top 100 genes prioritised or penalised by DPM

library("ggplot2")
library("data.table")
library('dplyr')
library('ggpubr')

# load the survival analysis results
df <- read.table("OV/survival_analysis_results.csv", sep=",", header=TRUE)
rownames(df) <- df$gene

# obtain the merged brown and DPM p-values
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

grouped_plotter <- function(x,plot_name){
      
      # identify significant DPM genes with a merged P-value < 0.05
      significant_dpm_merged <- dpm_merged[dpm_merged <= 0.05]
      
      # find genes with conflicting direction
      s_direction <- scores_dir/abs(scores_dir)
      disagreement <- rownames(s_direction[s_direction[,1] != s_direction[,2],])
      disagreement <- disagreement[!is.na(disagreement)]
      
      # identify the top 10 and top x penalized genes by DPM
      penalized <- brown_merged[brown_merged <= 0.05 & !names(brown_merged) %in% names(significant_dpm_merged)]
      penalized <- head(penalized[order(penalized)],x)
      penalized_10 <- head(penalized,10)

      # identify the top 10 and top x prioritized genes by DPM
      agreement <- dpm_merged[!names(dpm_merged) %in% disagreement & dpm_merged <= 0.05]
      top_prioritized <- head(agreement[order(agreement)],x)
      top_prioritized_10 <- head(top_prioritized,10)

      # prepare a dataframe for plotting
      gene_names_col <- c(names(top_prioritized),names(top_prioritized),names(penalized),names(penalized))
      dataset_contribution <- rep(rep(c("rna","protein"),each=x),times=2)
      label <- c(rep(c("Prioritized"),each=x*2),rep(c("Penalized"),each=x*2))
      
      log2hr_rna_penalized <- df[names(penalized),]$rna_log2hr
      log2hr_protein_penalized <- df[names(penalized),]$protein_log2hr
      log2hr_rna_prioritize <- df[names(top_prioritized),]$rna_log2hr
      log2hr_protein_prioritize <- df[names(top_prioritized),]$protein_log2hr
      log2hr <- c(log2hr_rna_prioritize,log2hr_protein_prioritize,log2hr_rna_penalized,log2hr_protein_penalized)

      df_bubble <- data.frame(
            genes = gene_names_col,
            dataset = dataset_contribution,
            log2hr = log2hr,
            label = label
            )
      
      df_bubble$colors <- "not10"
      df_bubble[df_bubble$genes %in% c(names(penalized_10),names(top_prioritized_10)),]$colors <- "top10"
      
      df_bubble <- df_bubble %>%
            filter(!is.na(log2hr)) %>%
            mutate(true_label = case_when(
                  log2hr < 0 & label == "Prioritized" ~ "Prioritized - Low Risk",
                  log2hr > 0 & label == "Prioritized" ~ "Prioritized - High Risk",
                  label == "Penalized" ~ "Penalized",
                  TRUE ~ "Neither"),
                  dataset = factor(dataset, levels=unique(dataset)),
                  genes = factor(genes, levels=unique(genes)),
                  label = factor(label, levels=unique(label)),
                  true_label = factor(true_label, levels = unique(true_label)),
                  colors = factor(colors, levels = unique(colors))
            )

      # plot
      plt = ggplot(df_bubble, aes(dataset, log2hr)) +
            geom_line(data = subset(df_bubble,colors=="not10"),aes(group = genes,color=colors), alpha = 0.15, linewidth = 0.2) + 
            geom_line(data = subset(df_bubble,colors=="top10"),aes(group = genes,color=colors), alpha = 0.9, linewidth = 0.8) +
            geom_hline(yintercept=0, linetype='solid', col = 'black', size = 1.0) + 
            geom_point(data = subset(df_bubble,colors=="not10"),aes(color=colors),size = 2.5, alpha = 1) + 
            geom_point(data = subset(df_bubble,colors=="top10"),aes(color=colors),size = 2.5, alpha = 1) +
            facet_grid(cols = vars(true_label)) +
            scale_color_manual(values = c("not10" = "#4d4d4d","top10" = "#ef3b2c")) +
            scale_y_continuous("log2HR") +
            scale_x_discrete(NULL) +
            theme_bw() +
            theme(    
                  plot.title = element_text(size = 20,margin = unit(c(0,0,3,0),"mm")),
                  plot.caption = element_text(size = 15),
                  plot.subtitle = element_text(size = 16),
                  axis.title.x = element_text(size = 18,margin = unit(c(2,0,0,0),"mm")),
                  axis.title.y = element_text(size = 18, margin = unit(c(0,4,0,0),"mm")),
                  axis.text.x = element_text(size = 15,angle = 90, hjust = 1, vjust=0.5, color = "black"),
                  axis.text.y = element_text(size = 15, color = "black"),
                  legend.title = element_text(size = 16),
                  legend.text = element_text(size = 15)
            )
      ggsave(plot=plt,filename=plot_name, useDingbats=FALSE,height = 4.5,width = 10)

}

grouped_plotter(x=100,"figure4_panel_c_plot.pdf")
