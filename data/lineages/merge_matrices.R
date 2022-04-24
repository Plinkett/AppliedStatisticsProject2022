# setwd("~/Documents/GitHub/AppliedStatisticsProject2022/data/lineages")
# b117 <- read.csv("B.1.1.7/binmatrix_b117.csv")
# ay4 <- read.csv("AY.4/binmatrix_ay4.csv")
# ba1 <- read.csv("BA.1/binmatrix_ba1.csv")
# ba11 <- read.csv("BA.1.1/binmatrix_ba11.csv")
# ay44 <- read.csv("AY.44/binmatrix_ay44.csv")
# ay3 <- read.csv("AY.3/binmatrix_ay3.csv")
# ay25 <- read.csv("AY.25/binmatrix_ay25.csv")
# ay43 <- read.csv("AY.43/binmatrix_ay43.csv")
# ay122 <- read.csv("AY.122/binmatrix_ay122.csv")
# Takes as input 2 binary matrices, the ones we built from the original dataset
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

# Example of union

finalmatrix <- create_mutations_set_lineages(binmatrix1, binmatrix2)
write.csv(data.frame(finalmatrix),file = "finalmatrix.csv", row.names = FALSE)
#matrice <- create_mutations_set_lineages(b117,ay4)
#matrice <- create_mutations_set_lineages(matrice,ba1)
#remove(b117, ay4, ba1)
#gc()
#matrice <- create_mutations_set_lineages(ay122, matrice)
#matrice <- create_mutations_set_lineages(matrice, ay25)
#matrice <- create_mutations_set_lineages(matrice, ay3)
#remove(ay122, ay25, ay3)
#gc()
#matrice <- create_mutations_set_lineages(matrice, ay43)
#matrice <- create_mutations_set_lineages(matrice, ay44)
#remove(ay43, ay44)
#gc()
#matrice <- create_mutations_set_lineages(matrice, ba11)
#ay251 <- read.csv("AY.25.1/binmatrix_ay251.csv")
#b16172
#b16172 <- read.csv("B.1.617.2/binmatrix_b16172.csv")
#matrice <- create_mutations_set_lineages(matrice, ay251)
#remove(ay251)
#gc()
#matrice <- create_mutations_set_lineages(matrice, b16172)
#remove(b16172)
#gc()
#View(matrice)
#write.csv(data.frame(matrice),file = "fullmatrix.csv", row.names = FALSE)
