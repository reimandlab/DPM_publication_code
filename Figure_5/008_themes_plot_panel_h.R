#### Figure 5H 
### Functional themes from the discovery dataset (TCGA, CPTAC) and
### validation dataset (GLASS, Oh et al.) are compared

library("ggplot2")
library('dplyr')

# validation themes from the annotated enrichment map
validation_themes <- c("Platelet aggregation", "Immune system response",
                       "Defense response","Cell-cell adhesion",
                       "Extracellular matrix", "MAPK and ERK signal transduction",
                       "Lymphocyte activation", "Cytokine production",
                       "Cell motility", "Proteolysis", "Response to hypoxia",
                       "Organismal development", "Tissue homeostasis",
                       "Actin filament", "Hemostasis",
                       "Wound healing", "Angiogenesis","Morphogenesis",
                       "Collagen metabolism", "Cell surface signaling",
                       "Carbohydrate catabolism","Carbohydrate biosynthesis",
                       "Basement membrane organization","Phagocytosis",
                       "Gland morphogenesis","Response to ion",
                       "Apoptosis", "Non-integrin membrane-ECM interactions",
                       "Integrin-mediated signaling pathway", "Response to nutrient",
                       "Cell junction organization", "Cellular component organization",
                       "Regulation of cell proliferation",
                       "Anatomical structure", "Phosphorous metabolism", 
                       "Response to endogenous stimulus",
                       "Immune cell apoptosis", "Cell killing")

# discovery themes from the annotated enrichment map in figure 5F
discovery_themes <- c("Taxis", 
                 "Organismal development",
                 "Morphogenesis", 
                 "Response to ion","Immune cell apoptosis",
                 "Response to hypoxia", "Extracellular matrix",
                 "Angiogenesis","Response to growth factor",
                 "Detoxification","Actin filament","Apoptosis",
                 "Cell motility","Cell-cell adhesion","Epithelial cell differentiation",
                 "Cytokine production", "ROS metabolism",
                 "Trans-synaptic signaling","Gliogenesis","Osteoblast differentiation",
                 "Fatty acid transport", "Response to mechanical stimulus",
                 "Regulation of cell proliferation",
                 "Cell fate commitment", "Phosphorous metabolism",
                 "FGFR signaling")

# identify all unique themes, themes unique to each dataset, and common themes
unique_themes <- unique(c(validation_themes,discovery_themes))
unique_validation <- validation_themes[!validation_themes %in% discovery_themes]
unique_discovery <- discovery_themes[!discovery_themes %in% validation_themes]
common <- validation_themes[validation_themes %in% discovery_themes]

# add themes to a dataframe and filter
df_themes <- data.frame(themes = rep(unique_themes, each = 2),
                        contribution =  c(rep(c("DISC", "VAL"), 50))
                        )
index1 <- rownames(df_themes[
      df_themes$contribution == "DISC" & df_themes$themes %in% unique_validation,
      ])
index2 <- rownames(df_themes[
      df_themes$contribution == "VAL" & df_themes$themes %in% unique_discovery,
      ])
df_themes <- df_themes[!rownames(df_themes) %in% c(index1,index2),]
row.names(df_themes) <- NULL

# add quantity of nodes in each theme by visual inspection
quantity_nodes <- c(1, 27, 10, 2, 16, 4, 10, 18, 21, 2, 13, 5, 10,
                    9, 3, 2, 13, 6, 2, 2, 5, 5, 7, 8, 13, 3, 9, 1,
                    3, 1, 2, 1, 1, 1, 13, 1, 3, 1, 1, 1, 1, 1, 1,
                    1, 4, 2, 2, 14, 1, 4, 2, 1, 12, 5, 2, 4, 1, 1,
                    1, 1, 1, 1, 1, 2
)
df_themes$nodes <- quantity_nodes

# prepare dataframe for plotting
df_themes <- df_themes %>%
      mutate(
            themes = factor(themes, levels=c(df_common,df_absent$themes, df_absent2$themes)),
            contribution = factor(contribution, levels=c("DISC","VAL")),
      )

# plot
pdf(file = paste("panelH_themes_plot",".pdf",sep = ""),   
    width = 12, 
    height = 3.5)
print(ggplot(df_themes, aes(x=themes, y=contribution, size=nodes)) + theme_bw() +
            geom_point(shape=15, color = '#0570b0') + 
            labs(title="", x ="Pathway Theme", y = "Dataset") + 
            scale_size_continuous(range = c(1, 7),breaks = c(1, 10, 20), labels = c("1", "10", "20")) + 
            theme(
                  plot.title = element_text(size=18,hjust = 0.5),
                  axis.title.x = element_text(size=12,margin = unit(c(0, 0, 0, 0), "mm")),
                  axis.title.y = element_text(size=12,margin = unit(c(0,4,0,0), "mm")),
                  axis.text.x = element_text(angle=45,vjust=1,hjust=1, size = 10),
                  axis.text.y = element_text(angle=0, size=10),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))
)

dev.off()
