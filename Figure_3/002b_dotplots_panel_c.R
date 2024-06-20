# visualize a dot-plot comparing KD and OE genes P-values and FC

library(optparse)
library(ggplot2)
library(forcats)


# options list for parser options
option_list <- list(
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



goi = c('CPED1', 'FAT1', 'CRYM', 'CACNA1H', 'NEGR1', 'B4GALT6')



dotplot = function(input_df){



# create plot
p = ggplot(input_df, aes(x = condition, y = gene, fill = logfc, size = pvalue)) + 

geom_point() +

scale_fill_gradient2(low = "blue", mid = "white", high = "firebrick", midpoint = 0)


print(p)

}




pdf(opt$figure_file, width = 8, height = 7)


# load df and convert
input_df = read.csv(opt$figure_data_file, sep='\t', row.names=1)

fc_limit = 2

# convert
input_df$pvalue = -log10(input_df$pvalue)

# subset to genes of interest
input_df = input_df[input_df$gene %in% goi,]


# create plots
dotplot(input_df)


dev.off()


print(opt$figure_file)



