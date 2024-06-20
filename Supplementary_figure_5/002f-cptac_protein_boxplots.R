# visualzie the protein expression of genes of interest in the CPTAC cohort

library(optparse)
library(ggplot2)
library(forcats)
library(gridExtra)


# options list for parser options
option_list <- list(
    make_option(c("-a","--figure_stats_file"), type="character", default=NULL,
            help="",
            dest="figure_stats_file"),
    make_option(c("-b","--figure_data_file"), type="character", default=NULL,
            help="",
            dest="figure_data_file"),
    make_option(c("-c","--figure_file"), type="character", default=NULL,
            help="",
            dest="figure_file")
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





dotplot = function(gene, raw_input_df, stats_df){


# subset to sample
raw_input_df = raw_input_df[raw_input_df$gene == gene,]
stats_df = stats_df[stats_df$gene == gene,]

wildtype = sum(raw_input_df$mutation_status == 'wildtype')
mutant = sum(raw_input_df$mutation_status == 'mutant')


p = ggplot(raw_input_df, aes(x = mutation_status, y = value, fill = mutation_status)) + plot_theme() +
geom_boxplot() +
geom_jitter(alpha = 0.35, pch=21) + 

labs(x = paste0("N wt=", wildtype, ", N mut =", mutant), y = "Expression (Z-score)") +
ggtitle(paste0(gene, ", p=", signif(stats_df$protein, 3)))

# theme(legend.title = element_text(size = 10))


print(p)
    
return()

}




pdf(opt$figure_file, width = 5, height = 5)

# load dfs
input_df = read.csv(opt$figure_data_file, sep='\t')
stats_df = read.csv(opt$figure_stats_file, sep='\t')


# refactor 
input_df$mutation_status = factor(input_df$mutation_status, levels = c('wildtype', 'mutant'))

# create plots
lapply(unique(input_df$gene), dotplot, input_df, stats_df)


dev.off()


print(opt$figure_file)



