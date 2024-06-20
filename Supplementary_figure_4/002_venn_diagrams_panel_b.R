#### Supplementary Figure 4B
### Venn diagram of enriched pathways of OS associations with mRNA and protein levels from 
### directional and non-directional analyses (ActivePathways, false discovery rate (FDR) < 0.05)

library('eulerr')
library("data.table")
library("dplyr")
library("ActivePathways")

# perform pathway enrichment analysis for the 9 remaining cancer types
cancer_types <- c("GBM","HNSC","LUAD","LUSC","PAAD","Renal","UCEC","BREAST","CRC") 
for (type in cancer_types){
      # load the survival analysis results
      df <- read.table(paste0("../Figure_4/",type,"/survival_analysis_results.csv"), sep=",", header=TRUE)
      
      # set the statistical background set to genes where transcriptomic and proteomic measurements were available
      scores <- data.frame(row.names = df$gene, rna = df$rna_pval, protein = df$protein_pval)
      ap_background <- rownames(scores[!is.na(scores$protein),])
      scores <- as.matrix(scores)
      scores[is.na(scores)] <- 1
      scores_dir <- data.frame(row.names = df$gene, rna = df$rna_log2hr, protein = df$protein_log2hr)
      scores_dir <- as.matrix(scores_dir)
      scores_dir[is.na(scores_dir)] <- 0
      constraints_vector <- c(1,1)

      # gmt file location
      fname_GMT = "../Figure_3/input_data/hsapiens.GO_BP_REACTOME.name.gmt"

      # run ActivePathways with Brown's merging method
      res_brown <- ActivePathways(scores=scores,
                                  background=ap_background,
                                  merge_method="Brown",
                                  gmt=fname_GMT,
                                  cytoscape_file_tag=paste0("Brown_",type),
                                  custom_colors=c("#fdb863","#b2abd2"),
                                  color_integrated_only="#bababa",
                                  geneset_filter=c(10,750),
                                  correction_method="fdr",
                                  cutoff=0.10,
                                  significant=0.05
      )
      
      # run ActivePathways with DPM
      res_dpm <- ActivePathways(scores=scores,
                                scores_direction=scores_dir,
                                constraints_vector=constraints_vector,
                                background=ap_background,
                                merge_method="DPM",
                                gmt=fname_GMT,
                                cytoscape_file_tag=paste0("DPM_",type),
                                custom_colors=c("#fdb863","#b2abd2"),
                                color_integrated_only="#bababa",
                                geneset_filter=c(10,750),
                                correction_method="fdr",
                                cutoff=0.10,
                                significant=0.05
      )

      # find the unique DPM, unique Brown, and shared pathways between the two methods
      unique_dpm <- length(res_dpm$term_id[!res_dpm$term_id %in% res_brown$term_id])
      unique_brown <- length(res_brown$term_id[!res_brown$term_id %in% res_dpm$term_id])
      shared_pathways <- length(res_brown$term_id[res_brown$term_id %in% res_dpm$term_id])
      
      print(type)
      print(paste0("DPM unique pathways: ", unique_dpm))
      print(paste0("Brown unique pathways: ", unique_brown))
      print(paste0("Shared pathways: " ,shared_pathways))
     
      # plot the venn diagram
      pdf(file = paste0("venn_pathways_",type,".pdf"),   
          width = 8, 
          height = 8)
      fit1 <- euler(c("Non-directional pathways" = unique_brown , "Directional pathways" = unique_dpm, 
                      "Non-directional pathways&Directional pathways" = shared_pathways))
      
      print(plot(fit1,
           fills = list(fill=c("#B3C28A","#2ca25f","#bababa",alpha = 0.2)))
      )
      dev.off()
}   