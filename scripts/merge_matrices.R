# Takes as input two lineages to be merged in one single matrix
# The columns that don't intersect are filled with zeroes
# This can be chained to combine even larger matrices

create_mutations_set_lineages <- function(df1, df2) {
  # Select column names a.k.a mutations from each set
  mutations1 <- colnames(df1)
  mutations2 <- colnames(df2)
  # Compute the union of these sets
  mutation_set <- union(mutations1, mutations2)
  # Compute the set difference between the whole set of mutations and each lineage's mutation set
  diff1 <- setdiff(mutation_set, mutations1)
  diff2 <- setdiff(mutation_set, mutations2)
  # For each mutation in diff1 add a column of zeroes, in this way we make the matrices homogeneous
  for(mutation in diff1) {
    # Binding of the column of zeroes
    df1 <- cbind(df1, integer(nrow(df1)))
    # We rename the variable
    names(df1)[length(names(df1))] <- mutation
  }
  for(mutation in diff2) {
    df2 <- cbind(df2, integer(nrow(df2)))
    names(df2)[length(names(df2))] <- mutation
  }
  # We sort the columns alphabetically to make the matrices fully homogeneous
  df1 <- df1[,sort(names(df1))]
  df2 <- df2[,sort(names(df2))]
  # All of this allows us to use rbind
  return(rbind(df1,df2))
}

# Creation of binary matrix with 2 lineages 

setwd("~/Documents/GitHub/AppliedStatisticsProject2022/data/lineages")
binmatrix_b117 <- read.csv("B.1.1.7/binmatrix_b117.csv")
binmatrix_ay4 <- read.csv("AY.4/binmatrix_ay4.csv")
binmatrix <- create_mutations_set_lineages(binmatrix_ay4, binmatrix_b117)
remove(binmatrix_ay4, binmatrix_b117)
gc()

write.csv(data.frame(binmatrix),file = "fullmatrix.csv", row.names = FALSE)
remove(binmatrix)
gc()