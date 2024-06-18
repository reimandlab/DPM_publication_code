#### Figure 4G
### Enrichment map of pathways and processes with OS associations (FDR < 0.05)

library("data.table")
library("dplyr")
library("ActivePathways")

# load the survival analysis results
df <- read.table("OV/survival_analysis_results.csv", sep=",", header=TRUE)

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
                            cytoscape_file_tag="Brown_",
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
                          cytoscape_file_tag="DPM_",
                          custom_colors=c("#fdb863","#b2abd2"),
                          color_integrated_only="#bababa",
                          geneset_filter=c(10,750),
                          correction_method="fdr",
                          cutoff=0.10,
                          significant=0.05
)

# save results
save(df, ap_background, res_brown, res_dpm, file = "panelG_Rdata.Rdata")
export_as_CSV(res_brown,"res_brown.csv")
export_as_CSV(res_dpm,"res_dpm.csv")

### aggregate the res_brown pathways with the res_dpm pathways, while 
### keeping track of lost/gained/maintained pathways between the two methods. 

# 1) Write the aggregated subgroups txt file
col_colors <- c("#fdb863","#b2abd2","#bababa")
tests <- c('rna','protein','combined')
all_pathways <- data.frame(term_id = c(res_brown$term_id,res_dpm$term_id))
evidence <- append(res_brown$evidence, res_dpm$evidence)
all_pathways$evidence <- evidence
head(all_pathways)

# retain the evidence contribution from one method (contributions for shared pathways will be the same).
sub_df = all_pathways %>% group_by(term_id) %>% summarise(evidence = list(evidence))
sub_df$evidence <- lapply(sub_df$evidence,'[[',1)
head(sub_df)

col_significance <- sub_df
evidence_columns = do.call(rbind, lapply(col_significance$evidence,
                                         function(x) 0+(tests %in% x)))
colnames(evidence_columns) = tests
col_significance = cbind(col_significance[,"term_id"], evidence_columns)
head(col_significance)

# check for lost/gained/maintained pathways between methods
lostp <- res_brown$term_id[!res_brown$term_id %in% res_dpm$term_id]
gainedp <- res_dpm$term_id[!res_dpm$term_id %in% res_brown$term_id]
sharedp <- res_brown$term_id[res_brown$term_id %in% res_dpm$term_id]
col_significance$directional_impact <- 0
col_significance[col_significance$term_id %in% lostp,]$directional_impact <- 1
col_significance[col_significance$term_id %in% gainedp,]$directional_impact <- 2


col_significance <- as.data.table(col_significance)
instruct.str <- paste('piechart:',
                      ' attributelist="', 
                      paste(tests, collapse=','),
                      '" colorlist="', 
                      paste(col_colors, collapse=','), 
                      '" showlabels=FALSE', sep='')
col_significance[, "instruct" := instruct.str]
utils::write.table(col_significance, 
                   file=paste0("COMBINED__", "subgroups.txt"), 
                   row.names=FALSE, 
                   sep="\t", 
                   quote=FALSE)


# 2) Write the aggregated pathways txt file
df_txtpathways <- data.frame(term_id = c(res_brown$term_id,res_dpm$term_id),
                             term_name = c(res_brown$term_name,res_dpm$term_name),
                             adjusted_p_val = c(res_brown$adjusted_p_val,res_dpm$adjusted_p_val))

sub_df2 = df_txtpathways %>% group_by(term_id,term_name) %>% summarise(adjusted_p_val = min(adjusted_p_val))
sub_df2
pathways_txt <- as.data.table(sub_df2)

utils::write.table(pathways_txt, 
                   file=paste0("COMBINED__", "pathways.txt"), 
                   row.names=FALSE, 
                   sep="\t", 
                   quote=FALSE)


# 3) Write the aggregated gmt file
gmt_main <- read.GMT("../Figure_3/input_data/hsapiens.GO_BP_REACTOME.name.gmt")
gmt_main <- gmt_main[pathways_txt$term_id]
write.GMT(gmt_main,paste0("COMBINED__","pathways.gmt"))
