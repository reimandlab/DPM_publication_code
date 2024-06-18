#### Figure 3D
### Venn diagram of enriched pathways found with directional and non-directional analyses

library('eulerr')

# load the enriched pathways from panel E
load("input_data/panelE_Rdata.Rdata")

# find the unique DPM, unique Brown, and shared pathways between the two methods
unique_dpm <- length(res_dpm$term_id[!res_dpm$term_id %in% res_brown$term_id])
unique_brown <- length(res_brown$term_id[!res_brown$term_id %in% res_dpm$term_id])
shared_pathways <- length(res_brown$term_id[res_brown$term_id %in% res_dpm$term_id])

# plot the venn diagram
fit1 <- euler(c("Non-directional pathways" = unique_brown , "Directional pathways" = unique_dpm, 
                "Non-directional pathways&Directional pathways" = shared_pathways))

plot(fit1,
     fills = list(fill=c("#B3C28A","#2ca25f","#bababa",alpha = 0.2)))
