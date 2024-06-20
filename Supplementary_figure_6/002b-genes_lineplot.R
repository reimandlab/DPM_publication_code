# create a lineplot of genes' pvalues using the two methods

library(optparse)
library(ggplot2)
library(forcats)
library(gridExtra)


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


opt = list(
    figure_data_file = '~/input_data/intermediate_files/genes_lineplot.tsv',
    figure_file = '~/input_data/intermediate_files/genes_lineplot.pdf'
)



lineplot = function(input_df, pval_cutoff){

# create plot
p = ggplot(input_df, aes(x = Brown, y = DPM)) + 
geom_point(size = 2.4,shape = 19, aes(color = ifelse(Brown <= pval_cutoff, "gray", ifelse(DPM > pval_cutoff, "#1F449C", "#F05039")))) +

labs(title="", x ="Brown's P (-log10)", y = "DPM P (-log10)") +

geom_hline(yintercept=pval_cutoff, linetype='dashed', col = 'black', linewidth = 0.5)+
geom_vline(xintercept = pval_cutoff, linetype = "dashed", col = "black", linewidth = 0.5) + 

geom_abline(linewidth = 0.5, slope=1,intercept = 0) +
scale_color_identity() +

ylim(0,15) + xlim(0,15) +

coord_cartesian() +

theme(plot.title = element_text(size=23,hjust = 0.5),
    axis.title.x = element_text(size=18,margin = unit(c(2, 0, 0, 0), "mm")),
    axis.title.y = element_text(size=18,margin = unit(c(0,4,0,0), "mm")),
    axis.text = element_text(size = 16),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")) 

print(p)

}






pdf(opt$figure_file, width = 8, height = 7)


# load df and convert
input_df = read.csv(opt$figure_data_file, sep='\t', row.names=1)
input_df = -log10(input_df)

# log10 cutoff for P-value is 0.05
pval_cutoff = 1.301


# create plots
lineplot(input_df, pval_cutoff)


dev.off()


print(opt$figure_file)



