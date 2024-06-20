# create venn diagrams of overlapping genes

library(optparse)
library(ggplot2)
library(forcats)
library(gridExtra)
library(eulerr)

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






# plot venn diagram
venn_plot = function(condition, input_df){

# subset to condition
sub_df = input_df[input_df$upregulated == condition,]
if (condition == 'none'){
    sub_df = input_df
}

# create list for each condition
data_list = list(Protein = sub_df$gene[sub_df$source == 'protein'], RNA = sub_df$gene[sub_df$source == 'rna'], Methylation = sub_df$gene[sub_df$source == 'methylation'])
venn = euler(data_list)
    
# plot
p = plot(venn, counts = TRUE, quantities = TRUE, fills = list(fill=c("#b2abd2","#fdb863","#66c2a5","#8dd3c7","#ffffb3","#80b1d3","#d9d9d9",alpha=0.2)))

print(p)
}


# load df
input_df = read.csv(opt$figure_data_file, sep='\t')

pdf(opt$figure_file, width = 7, height = 7)

conditions = c("True", "False", "none")
lapply(conditions, venn_plot, input_df)

dev.off()

print(opt$figure_file)



