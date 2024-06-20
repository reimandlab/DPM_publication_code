# merge P-values from different sources using directional and non-directional methods

library(optparse)
library(ActivePathways)
library(dplyr)

# options list for parser options
option_list <- list(
    make_option(c("-a","--gene_expression_file"), type="character", default=NULL,
            help="",
            dest="gene_expression_file"),
    make_option(c("-b","--protein_expression_file"), type="character", default=NULL,
            help="",
            dest="protein_expression_file"),
    make_option(c("-c","--methylation_file"), type="character", default=NULL,
            help="",
            dest="methylation_file"),
    make_option(c("-d","--combined_pval_file"), type="character", default=NULL,
            help="",
            dest="combined_pval_file"),
    make_option(c("-e","--combined_fc_file"), type="character", default=NULL,
            help="",
            dest="combined_fc_file"),
    make_option(c("-f","--merged_pvalues_file"), type="character", default=NULL,
            help="",
            dest="merged_pvalues_file")
            )

parser <- OptionParser(usage = "%prog -i data.csv etc.", option_list=option_list)

# get the parse information and assign to groups in the opt variable
opt = parse_args(parser)

plot_theme = function(...) {
    theme_bw() +
    theme(    
        plot.title = element_text(size = 22),
        plot.caption = element_text(size = 12),
        plot.subtitle = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 12,
                angle = 90, hjust = 1, vjust=0.5, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        ...
    )
}


opt = list()
opt$gene_expression_file = '~/input_data/intermediate_files/tcga_gbm_degs.tsv'
opt$protein_expression_file = '~/input_data/intermediate_files/protein_degs.tsv'
opt$methylation_file = '~/input_data/intermediate_files/tcga_gbm_methylation_degs.tsv'
opt$combined_pval_file = '~/output_data/intermediate_files/figure5_pvals.tsv'
opt$combined_fc_file = '~/output_data/intermediate_files/figure5_fcs.tsv'
opt$merged_pvalues_file = '~/output_data/intermediate_files/figure5_merged_pvalues.tsv'




# load the dfs
gene_expression = read.csv(opt$gene_expression_file, sep='\t', row.names=1)
protein_expression = read.csv(opt$protein_expression_file, sep='\t', row.names=1)
methylation = read.csv(opt$methylation_file, sep='\t', row.names=1)


# Identify common row names across all data frames
common_rows <- Reduce(intersect, list(rownames(gene_expression), rownames(protein_expression), rownames(methylation)))

# subset to common_rows
gene_expression = gene_expression[common_rows,]
protein_expression = protein_expression[common_rows,]
methylation = methylation[common_rows,]

# separate into p_value and fold_change dfs
pval_df = gene_expression[,'pvalue',drop=FALSE]
colnames(pval_df) = 'rna'
pval_df$protein = protein_expression$p_value
pval_df$methylation = methylation$pvalue


fc_df = gene_expression[,'FC',drop=FALSE]
colnames(fc_df) = 'rna'
fc_df$protein = protein_expression$fold_change
fc_df$methylation = methylation$FC

fc_df[is.na(fc_df)] = 1

# log-transform the fc
fc_df = log2(fc_df)


fc_df[is.na(fc_df)] = 0
pval_df[is.na(pval_df)] = 1

# merge pvals using brown's method
browns_df = merge_p_values(as.matrix(pval_df), method="Brown")

# repeat using DPM
dpm_df = merge_p_values(as.matrix(pval_df), method='DPM', scores_direction = as.matrix(fc_df), constraints_vector = c(1,1,-1))

# combine and save to file
res_df = data.frame(Brown = browns_df, DPM = dpm_df)
write.table(res_df, file=opt$merged_pvalues_file, sep='\t', quote=FALSE)


# save of dfs to files
write.table(pval_df, file=opt$combined_pval_file, sep='\t', quote=FALSE)
write.table(fc_df, file=opt$combined_fc_file, sep='\t', quote=FALSE)





