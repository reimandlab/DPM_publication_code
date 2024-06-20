# visualize the directional agreements and disagreement between data sources in a heatmap

library(optparse)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(gtools)


# options list for parser options
option_list <- list(
    make_option(c("-a","--figure_data_file"), type="character", default=NULL,
            help="",
            dest="figure_data_file"),
    make_option(c("-b","--figure_file"), type="character", default=NULL,
            help="",
            dest="figure_file")
            )

# get the parse information and assign to groups in the opt variable
parser <- OptionParser(usage = "%prog -i data.csv etc.", option_list=option_list)
opt = parse_args(parser)


sub_df = read.csv(opt$figure_data_file, sep='\t', row.names=1)
rownames(sub_df) = rownames(sub_df)

pdf(opt$figure_file, width=30, height=30)



# colors
myCol <- colorRampPalette(c('#301abc', 'white', '#d50000'))(100)


# create plot matrices
fc_mat = sub_df[,c('rna', 'protein', 'methylation')]
rownames(fc_mat) = rownames(sub_df)

# add direction 
site_types_cols = sub_df$direction


# define distance metric
hr = hclust(dist(fc_mat, method = 'euclidean'), method="complete") 

heatmap.2(as.matrix(fc_mat), trace = "none", 
Colv = FALSE,
Rowv = as.dendrogram(hr),
col = myCol, 
key=TRUE,
main = "Complete euclidean distance", 
dendrogram = "row",
RowSideColors = site_types_cols,
labRow = sub_df$significant,
colRow = site_types_cols,
cexCol=0.75,
cexRow=1,
notecol = "black")


site_types_cols = sub_df$cancer_gene


# define distance metric
hr = hclust(dist(fc_mat, method = 'euclidean'), method="complete") 

heatmap.2(as.matrix(fc_mat), trace = "none", 
Colv = FALSE,
Rowv = as.dendrogram(hr),
col = myCol, 
key=TRUE,
main = "Complete euclidean distance", 
dendrogram = "row",
RowSideColors = site_types_cols,
labRow = sub_df$significant,
colRow = site_types_cols,
cexCol=0.75,
cexRow=1,
notecol = "black")



dev.off()


print(opt$figure)









