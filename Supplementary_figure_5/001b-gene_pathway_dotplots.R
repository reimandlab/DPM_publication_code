# create dotplots of the contribution of each gene to the pathway

library(optparse)
library(forcats)
library(ggplot2)

# options list for parser options
option_list <- list(
    make_option(c("-a","--figure_data_file"), type="character", default=NULL,
            help="",
            dest="figure_data_file"),
    make_option(c("-b","--figure_stats_file"), type="character", default=NULL,
            help="",
            dest="figure_stats_file"),
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





dot_plot = function(pathway, input_df){

# subset to pathway
sub_df = input_df[input_df$gmt_key == pathway,]

# order factors
sub_df$gene = fct_reorder(sub_df$gene, sub_df$log_p, max, .desc=TRUE)

fill_values = c('firebrick', 'dodgerblue')
names(fill_values) = c('up', 'down')


p = ggplot(data=sub_df, aes(x = gene, y = sample, fill = log2fc, size = log_p)) + plot_theme() +
geom_point(pch=21) +

# scale_fill_manual(values = fill_values) +
scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, limits=c(-2,2)) +

ggtitle(paste0(unique(sub_df$gmt_key), '-', signif(unique(sub_df$adjusted_p_val), 3))) + ylab('') + xlab('') +

guides(fill = guide_colourbar(title="Average FC \n(log2)", ticks.colour = "black", frame.colour = 'black')) +
theme(strip.placement = "outside", strip.background = element_blank(),
    legend.position="right", plot.title = element_text(size = 12))


print(p)

return()

}



pdf(opt$figure_file, width=7, height = 3)

p_limit = 8
fc_limit = 2

# loads dfs
input_df = read.csv(opt$figure_data_file, sep='\t')

# subset to pathways of interest
input_df = input_df[grepl('fibroblast', input_df$gmt_key),]

# pvalue limits
input_df$log_p = -log10(input_df$pvalue)
input_df$log_p[input_df$log_p > p_limit] = p_limit

# fc limits
input_df$log2fc[input_df$log2fc > fc_limit] = fc_limit
input_df$log2fc[input_df$log2fc < -fc_limit] = -fc_limit


# factor pathways according to pval
input_df$gmt_key = fct_reorder(input_df$gmt_key, input_df$adjusted_p_val, .desc=TRUE)


# lapply(levels(input_df$gmt_key), dot_plot, input_df)
lapply(unique(input_df$gmt_key), dot_plot, input_df)

dev.off()


print(opt$figure_file)









