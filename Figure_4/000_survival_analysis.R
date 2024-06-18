#### Figure 4
### Integrating ovarian cancer transcriptomes and proteomes with patient survival
### information for pathway and biomarker analyses

library('survival')
library("ggplot2")
library("data.table")
library("dplyr")

load_and_preprocess <- function(clinical_directory, rna_directory, protein_directory){
      # load clinical data, rename patient ids and remove duplicate entries
      clinical_data <- read.table(clinical_directory, sep="\t", header=TRUE)
      clinical_data$case_submitter_id <- unlist(lapply(clinical_data$case_submitter_id, function(x) {gsub("-", ".", x)}))
      clinical_data <- clinical_data[order(clinical_data[,'case_submitter_id']),]
      clinical_data <- clinical_data[!duplicated(clinical_data$case_submitter_id),]
      
      # further preprocess the clinical data for survival analysis
      clinical_data <- clinical_data[!clinical_data$vital_status == "Not Reported",]
      clinical_data[clinical_data$vital_status == "Alive",]$vital_status <- 0
      clinical_data[clinical_data$vital_status == "Dead",]$vital_status <- 1
      clinical_data[is.na(clinical_data$days_to_death),]$days_to_death <- clinical_data[
            is.na(clinical_data$days_to_death), "days_to_last_follow_up"
      ]
      clinical_data$days_to_last_follow_up <- NULL
      clinical_data$days_to_last_known_disease_status <- NULL
      
      # the normalized rna and protein expression data was acquired from:
      # https://github.com/chadcreighton/cancer-proteomics-compendium-n2002
      # these datasets contain a "Gene" column with several numeric names, dates, and gene aliases
      # that need to be corrected
      gene_aliases <- read.table("input_data/alias_corrected_names.csv", sep=",", header=TRUE)
      rna_numeric_genes <- read.table("input_data/rna_rename_numeric_genes.tsv", sep=",", header=TRUE)
      protein_numeric_genes <- read.table("input_data/protein_rename_numeric_genes.tsv", sep=",", header=TRUE)
            
      # load rna data, fix aliases and numeric gene names
      rna_data <- read.table(rna_directory, sep =",", header=TRUE)
      rna_gene_aliases <- gene_aliases[gene_aliases$initial_alias %in% rna_data$Gene,]
      rna_data$Gene[match(rna_gene_aliases$initial_alias, rna_data$Gene)] <- rna_gene_aliases$name
      rna_data$Gene[match(rna_numeric_genes$Gene, rna_data$Gene)] <- rna_numeric_genes$correct_name
      
      # remove replicate gene entries by taking the median across replicates
      rna <- as.data.frame(rna_data %>% group_by(Gene) %>% summarise(across(everything(),list(median))))
      rownames(rna) <- rna$Gene
      rna$Gene <- NULL
      colnames(rna) <- unlist(lapply(colnames(rna),function(x) substring(x,0,nchar(x)-2)))

      # load protein data, fix aliases and numeric gene names
      protein_data <- read.table(protein_directory,sep=",",header=TRUE)
      protein_gene_aliases <- gene_aliases[gene_aliases$initial_alias %in% protein_data$Gene,]
      protein_data$Gene[match(protein_gene_aliases$initial_alias, protein_data$Gene)] <- protein_gene_aliases$name 
      protein_data$Gene[match(protein_numeric_genes$Gene, protein_data$Gene)] <- protein_numeric_genes$correct_name
      
      # remove replicate genes entries by taking the median across replicates
      protein <- as.data.frame(protein_data %>% group_by(Gene) %>% summarise(across(everything(),list(median))))
      rownames(protein) <- protein$Gene
      protein$Gene <- NULL
      colnames(protein) <- unlist(lapply(colnames(protein),function(x) substring(x,0,nchar(x)-2)))

      return(list(clinical_data, rna, protein))
}


survival_analysis <- function(gene, dataset, survival_info, omic_name, cancer){
      ## Median dichotomized expression survival analysis. Returns 
      ## the p-value and log2hr of a gene.
      
      # filter survival and expression data for shared patient ids
      df_si <- survival_info[survival_info$case_submitter_id %in% colnames(dataset),]
      dataset <- dataset[,colnames(dataset) %in% df_si$case_submitter_id]
      
      # acquire the expression information for the gene of interest
      expression <- dataset[rownames(dataset) == gene,]
      
      # remove NA entries
      n_expression <- names(expression)[!is.na(expression)]
      expression <- expression[!is.na(expression)]
      names(expression) <- n_expression
      
      # obtain the median expression and split patients into two groups
      m_exp <- median(as.numeric(expression))
      expression[expression < m_exp] <- FALSE
      expression[expression != FALSE] <- TRUE

      # prepare survival data for coxph model
      df_expression <- data.frame(patient = names(expression), group = as.numeric(expression))
      df_si <- df_si[df_si$case_submitter_id %in% df_expression$patient,]
      
      df_si$group <- df_expression$group[match(df_si$case_submitter_id,df_expression$patient)]
      df_si$group <- as.numeric(df_si$group)
      df_si$vital_status <- as.numeric(df_si$vital_status)
      df_si$days_to_death <- as.numeric(df_si$days_to_death)
      df_si$days_to_birth <- as.integer(df_si$days_to_birth/-365)
      
      if (any(c("Stage IB","Stage IA","Stage I") %in% df_si$ajcc_pathologic_stage)){
            df_si[df_si$ajcc_pathologic_stage %in% c("Stage IB","Stage IA","Stage I"),]$ajcc_pathologic_stage <- 1
      }
      if (any(c("Stage IIB","Stage IIA","Stage II") %in% df_si$ajcc_pathologic_stage)){
            df_si[df_si$ajcc_pathologic_stage %in% c("Stage IIB","Stage IIA","Stage II"),]$ajcc_pathologic_stage <- 2
      }
      if (any(c("Stage IIIB","Stage IIIA","Stage III") %in% df_si$ajcc_pathologic_stage)){
            df_si[df_si$ajcc_pathologic_stage %in% c("Stage IIIB","Stage IIIA","Stage III"),]$ajcc_pathologic_stage <- 3
      }
      if (any(c("Stage IVB","Stage IVA","Stage IV") %in% df_si$ajcc_pathologic_stage)){
            df_si[df_si$ajcc_pathologic_stage %in% c("Stage IVB","Stage IVA","Stage IV"),]$ajcc_pathologic_stage <- 4
      }
      if (any("Not Reported" %in% df_si$ajcc_pathologic_stage)){
            df_si[df_si$ajcc_pathologic_stage == "Not Reported",]$ajcc_pathologic_stage <- 0
      }
      if (any(NA %in% df_si$ajcc_pathologic_stage)){
            df_si[is.na(df_si$ajcc_pathologic_stage),]$ajcc_pathologic_stage <- 0
      }
      df_si$ajcc_pathologic_stage <- as.numeric(df_si$ajcc_pathologic_stage)
      
      df_si <- df_si[,c("days_to_death","vital_status","group","ajcc_pathologic_stage","days_to_birth",
                                                "gender")]
      
      # for each gene ensure there are at least 10 patients in each median dichotomized group
      g1_len <- length(df_si[df_si$group == 0,]$group)
      g2_len <- length(df_si[df_si$group == 1,]$group)
      if (g1_len >= 10 & g2_len >= 10){
            
            if (length(unique(df_si$gender)) == 2){
                  so <- 'Surv(days_to_death, vital_status) ~ ajcc_pathologic_stage + days_to_birth + gender'
                  fo <- as.formula(so)
                  cox.mod_ho <- coxph(fo, data = df_si)
                  
                  sa <- 'Surv(days_to_death, vital_status) ~ ajcc_pathologic_stage + days_to_birth + gender + group'
                  fa <- as.formula(sa)
                  cox.mod_ha <- coxph(fa, data = df_si) 
                  
                  output_anova <- anova(cox.mod_ho,cox.mod_ha,test='chisq')
                  pval <- output_anova[,4][2]
                  hr <- as.numeric(coef(summary(cox.mod_ha))[,2][4])
                  return(c(gene, pval, log2(hr)))
            } else {
                  so <- 'Surv(days_to_death, vital_status) ~ ajcc_pathologic_stage + days_to_birth'
                  fo <- as.formula(so)
                  cox.mod_ho <- coxph(fo, data = df_si)
                  
                  sa <- 'Surv(days_to_death, vital_status) ~ ajcc_pathologic_stage + days_to_birth + group'
                  fa <- as.formula(sa)
                  cox.mod_ha <- coxph(fa, data = df_si)
                  
                  output_anova <- anova(cox.mod_ho,cox.mod_ha,test='chisq')
                  pval <- output_anova[,4][2]
                  hr <- as.numeric(coef(summary(cox.mod_ha))[,2][3])
                  
                  # Figure 4D prepare KM plot for ACTN4 and PIK3R4 in OV cancer
                  if (gene %in% c("ACTN4","PIK3R4") & cancer == "OV"){
                        kmsurvival <- survfit(Surv(df_si$days_to_death,df_si$vital_status) ~ df_si$group)
                        pdf(paste(paste("KM",gene,omic_name,sep = "_"),"pdf",sep = "."),width = 7, height = 7)
                        par(mar=c(5,6,4,1)+.1)
                        plot(kmsurvival, col = 1:2, lwd=4,mark.time = TRUE, xlab = "Time (days)",
                             ylab = "Survival Probability",cex.axis=1.8,cex.lab = 2,
                             main = paste(gene,omic_name,sep = " "), cex.main = 2)
                        legend("bottomleft",legend=c("Low","High"),col=1:2,lty = 1,lwd = 4,horiz=FALSE,
                               bty="n",cex = 1.8)
                        legend("topright",legend=c(paste("p =",formatC(pval,format="e",digits = 3),sep = " "),paste("log2(HR) =",round(log2(hr),5),sep = " ")),
                               bty = "n", cex = 1.2)
                        dev.off()
                  }
                  return(c(gene, pval, log2(hr)))
            }
            
      } else{
            return(c(gene,NA,NA))
      }
}


main <- function(cancer_type, cohort, subcohort_name){
      
      # create a directory for each cancer type
      dir.create(cancer_type)
    
      # load and preprocess the clinical, rna, and protein datasets
      filename_dataset <- paste("input_data/", cohort, "/", cancer_type, "-" ,subcohort_name, "_", sep = "")
      loaded_data <- load_and_preprocess(
            clinical_directory = paste("input_data", cohort, "clinical_and_survival_information.txt", sep="/"),
            rna_directory = paste(filename_dataset,"normalized_mRNA.csv", sep =""),
            protein_directory = paste(filename_dataset, "normalized_total_protein.csv", sep = "")
            )
      survival_info <- loaded_data[[1]]
      rna <- loaded_data[[2]]
      protein <- loaded_data[[3]]

      # perform median dichotomized survival coxph analysis for each
      # gene and acquire p-value and log2hr value
      protein_pval <- c()
      protein_hazard <- c()
      for (x in rownames(protein)){
            survival_result <- survival_analysis(
                  gene=x,
                  dataset=protein,
                  survival_info=survival_info,
                  omic_name="Protein",
                  cancer=cancer_type
                  )
            protein_pval <- c(protein_pval, survival_result[2])
            protein_hazard <- c(protein_hazard, survival_result[3])
      }
      protein_summary <- data.frame(gene = rownames(protein), protein_pval = protein_pval, protein_log2hr = protein_hazard,
                                    protein_fdr = p.adjust(as.numeric(protein_pval), method = "fdr"))
      protein_summary <- protein_summary[order(protein_summary$protein_pval),]

      # perform median dichotomized survival coxph analysis for each
      # gene and acquire p-value and log2hr value
      rna_pval <- c()
      rna_hazard <- c()
      for (x in rownames(rna)){
            survival_result <- survival_analysis(
                  gene=x,
                  dataset=rna,
                  survival_info=survival_info,
                  omic_name="RNA",
                  cancer=cancer_type)
            rna_pval <- c(rna_pval, survival_result[2])
            rna_hazard <- c(rna_hazard, survival_result[3])
      }
      rna_summary <- data.frame(gene = rownames(rna), rna_pval = rna_pval, rna_log2hr = rna_hazard,
                                rna_fdr = p.adjust(as.numeric(rna_pval), method = "fdr"))
      rna_summary <- rna_summary[order(rna_summary$rna_pval),]

      # join the two datasets
      df_combined <- full_join(rna_summary, protein_summary, by = "gene")
      print(head(df_combined))
      write.csv(df_combined,paste0(cancer_type,"/survival_analysis_results.csv"), quote=F, row.names= F)
}

# perform coxph survival analysis using RNA and Protein median dichotomized expression along with
# clinical covariates across 10 cancer types
cancer_types <- c("GBM","HNSC","LUAD","LUSC","PAAD","Renal","UCEC","BREAST","CRC","OV")     
for (type in cancer_types){
      if (type %in% c("BREAST","CRC","OV")){
            main(
                  cancer_type=type,
                  cohort="tcga_clinical",
                  subcohort_name="CPTAC-TCGA"
                  )
      } else {
            main(
                  cancer_type=type,
                  cohort="cptac_clinical",
                  subcohort_name="CPTAC"
                  )
      }}
