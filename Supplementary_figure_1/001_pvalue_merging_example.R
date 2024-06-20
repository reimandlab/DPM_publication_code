#### Supplementary Figure 1
### Example of directional P-value merging (DPM) in three input omics datasets. 

library(ActivePathways)

# prepare pvalue matrix
pvals <- matrix(
      data=rep(c(0.05,0.1,0.20),4),
      nrow=4,
      ncol=3,
      byrow=TRUE
      )
colnames(pvals) <- c("Dataset 1", "Dataset 2", "Dataset 3")

# prepare directions matrix
directions <- matrix(
      data=c(1,1,0,-1,-1,0,1,-1,0,-1,1,0),
      nrow=4,
      ncol=3,
      byrow=TRUE
      )
colnames(directions) <- c("Dataset 1", "Dataset 2", "Dataset 3")

# set contraints vector such that Dataset 3 contains no directional information,
# while Dataset 1 and Dataset 2 have direct directional relationship
constraints_vector <- c(1,1,0)

# merge data using DPM
merge_p_values(
      scores=pvals,
      method="DPM",
      scores_direction=directions,
      constraints_vector=constraints_vector
)
