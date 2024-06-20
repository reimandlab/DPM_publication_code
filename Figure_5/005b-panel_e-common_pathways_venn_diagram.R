# create a venn diagram of the common pathways between the two merging methods

library(optparse)
library(ggplot2)
library(data.table)
library(dplyr)
library(eulerr)

# options list for parser options
option_list <- list(
    make_option(c("-a","--figure_data_file"), type="character", default=NULL,
            help="",
            dest="figure_data_file"),
    make_option(c("-b","--figure_file"), type="character", default=NULL,
            help="",
            dest="figure_file")
            )

parser <- OptionParser(usage = "%prog -i data.csv etc.", option_list=option_list)

# get the parse information and assign to groups in the opt variable
opt = parse_args(parser)


# load df
input_df = read.csv(opt$figure_data_file, sep='\t')


pdf(opt$figure_file, width = 7, height = 7)


# create list for each condition
data_list = list(Brown = input_df$term_id[input_df$source == 'brown'], DPM = input_df$term_id[input_df$source == 'dpm'])
venn = euler(data_list)

print(data_list)
    
# plot
p = plot(venn, counts = TRUE, quantities = TRUE, fills = list(fill=c("#B3C28A","#2ca25f","#bababa",alpha=0.2)))

print(p)



dev.off()




