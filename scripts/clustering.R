library(tidyverse)
library(broom)
library(palmerpenguins)

# Perform clustering on our data
# Idea: shuffle the data and select 1000 samples from 6 lineages
#       after selecting those samples perform k-means clustering (k = 6)
#       and, if possible, hierarchical clustering. On the whole data?
#       Let's try the principal components

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

################################################################

# For B.1.1.7
b117 <- read_csv("Documents/GitHub/AppliedStatisticsProject2022/data/lineages/B.1.1.7/binmatrix_b117.csv")
b117shuffle <- b117[sample(nrow(b117)),]
b117shuffle <- b117shuffle[sample(nrow(b117)),]

# Select 10k

b117_topk <- b117shuffle[1:10000,-c(2,3)]
remove(b117, b117shuffle)
gc()

# For AY.4
ay4 <- read_csv("Documents/GitHub/AppliedStatisticsProject2022/data/lineages/AY.4/binmatrix_ay4.csv")
ay4shuffle <- ay4[sample(nrow(ay4)),]
ay4shuffle <- ay4[sample(nrow(ay4)),]

ay4_topk <- ay4shuffle[1:10000,-c(2,3)]
remove(ay4, ay4shuffle)
gc()

# For BA.1
ba1 <- read_csv("Documents/GitHub/AppliedStatisticsProject2022/data/lineages/BA.1/binmatrix_ba1.csv")
ba1shuffle <- ba1[sample(nrow(ba1)),]
ba1shuffle <- ba1[sample(nrow(ba1)),]

ba1_topk <- ba1shuffle[1:10000,-c(2,3)]
remove(ba1, ba1shuffle)
gc()


# For BA.1.1
ba11 <- read_csv("Documents/GitHub/AppliedStatisticsProject2022/data/lineages/BA.1.1/binmatrix_ba11.csv")
ba11shuffle <- ba11[sample(nrow(ba11)),]
ba11shuffle <- ba11[sample(nrow(ba11)),]

ba11_topk <- ba11shuffle[1:10000,-c(2,3)]
remove(ba11, ba11shuffle)
gc()

# For AY.103
#ay103 <- read_csv("Documents/GitHub/AppliedStatisticsProject2022/data/lineages/AY.103/binmatrix_ay103.csv")
#ay103shuffle <- ay103[sample(nrow(ay103)),]
#ay103shuffle <- ay103[sample(nrow(ay103)),]

#ay103_topk <- ay103shuffle[1:10000,-c(2,3)]
#remove(ay103, ay103shuffle)
#gc()

# For AY.44
#ay44 <- read_csv("Documents/GitHub/AppliedStatisticsProject2022/data/lineages/AY.44/binmatrix_ay44.csv")
#ay44shuffle <- ay44[sample(nrow(ay44)),]
#ay44shuffle <- ay44[sample(nrow(ay44)),]

#ay44_topk <- ay44shuffle[1:10000,-c(2,3)]
#remove(ay44, ay44shuffle)
#gc()


# Merge the datasets into one
binmatrix <- create_mutations_set_lineages(ay4_topk, b117_topk)
#binmatrix <- create_mutations_set_lineages(binmatrix, ay103_topk)
binmatrix <- create_mutations_set_lineages(binmatrix, ba11_topk)
binmatrix <- create_mutations_set_lineages(binmatrix, ba1_topk)
#binmatrix <- create_mutations_set_lineages(binmatrix, ay44_topk)

forpc <- binmatrix[,-2]
pc <- princomp(forpc, score =T)
summary(pc)

# Plot the scree plot 

variance <- pc$sdev^2 / sum(pc$sdev^2)
variance
library(ggplot2)

var_explained_df <- data.frame(PC= paste0("PC",1:9),
                               var_explained=(pc$sdev[1:9])^2/sum((pc$sdev[1:9])^2))
# With ggplot2
X11()
var_explained_df %>%
  ggplot(aes(x=PC,y=var_explained, group=1))+
  geom_point(size=4)+
  geom_line()+
  labs(title="Scree plot: PCA on scaled data")

# Plain old R (professor Arnone's code)
X11()
plot(cumsum(pc$sd^2)/sum(pc$sd^2), type='b', axes=F, xlab='number of components', 
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(binmatrix),labels=1:ncol(binmatrix),las=2)
graphics.off()

# Plot along first 2 PCs for first 4 lineages
x11()
par(mfrow = c(1,3))
plot(pc$scores[,1], pc$scores[,2], xlab="PC1", ylab="PC2", main="Scores with no jitter", col=factor(binmatrix$lineage))
abline(h=0, v=0, lty=2, col='black')
plot(pc$scores[,1] + rnorm(40000, sd=0.005), pc$scores[,2] + rnorm(40000, sd=0.01), xlab="PC1", ylab="PC2", main ="Scores with gaussian jitter (sd = 0.01)", col=factor(binmatrix$lineage))
abline(h=0, v=0, lty=2, col='black')
plot(pc$scores[,1] + rnorm(40000, sd=0.01), pc$scores[,2] + rnorm(40000, sd=0.1), xlab="PC1", ylab="PC2", main ="Scores with gaussian jitter (sd = 0.1)", col=factor(binmatrix$lineage))
abline(h=0, v=0, lty=2, col='black')

# The plots with jitter give us an idea of the sizes of the clusters before-hand

X11()
plot(pc$scores[,1] + rnorm(40000, sd=0), pc$scores[,2] + rnorm(40000, sd=0), xlab="PC1", ylab="PC2", main ="Scores with gaussian jitter (sd = 0.1)", col=factor(binmatrix$lineage))
abline(h=0, v=0, lty=2, col='black')
legend('bottomright',unique(binmatrix$lineage),col=1:6,pch=15)
X11()
barplot(pc$loadings[which(abs(pc$loadings[,1]) > 0.1),1], main="PC 1 Loadings Plot", las=2)

length(which(abs(pc$loadings[,1]) > 0.005))


# What if we look at one lineage? 

# For B.1.1.7
b117 <- read_csv("Documents/GitHub/AppliedStatisticsProject2022/data/lineages/AY.4/binmatrix_ay4.csv")
b117shuffle <- b117[sample(nrow(b117)),]
b117shuffle <- b117shuffle[sample(nrow(b117)),]

# Select 60k sample from B.1.1.7 (English variant)

b117_topk <- b117shuffle[1:60000,-c(2,3)]
remove(b117, b117shuffle)
gc()

forpc <- b117_topk[,-1]
pc <- princomp(forpc, score =T)
summary(pc)

# scree plot
X11()
plot(cumsum(pc$sd^2)/sum(pc$sd^2), type='b', axes=F, xlab='number of components', 
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(b117_topk),labels=1:ncol(b117_topk),las=2)
graphics.off()

X11()
plot(pc$scores[,1] + rnorm(60000, sd=0), pc$scores[,2] + rnorm(60000, sd=0), xlab="PC1", ylab="PC2", main ="Scores with gaussian jitter (sd = 0.1)")
abline(h=0, v=0, lty=2, col='black')


# Hierarchical clustering with small number of samples
b117 <- read.csv("Documents/GitHub/AppliedStatisticsProject2022/data/lineages/B.1.1.7/binmatrix_b117.csv")
b117shuffle <- b117[sample(nrow(b117)),]
b117shuffle <- b117shuffle[sample(nrow(b117)),]
b117_topk <- b117shuffle[1:10000,-c(2,3)]
remove(b117, b117shuffle)
gc()

forpc <- b117_topk[,-1]
pc <- princomp(forpc, score =T)
summary(pc)

b117_topk <- b117_topk[,-1]
distance_mat <- dist(b117_topk, method = 'euclidean')
Hierar_cl <- hclust(distance_mat, method = "ward.D")
X11()
plot(Hierar_cl)
abline(h = 110, col = "green")
fit <- cutree(Hierar_cl, k = 3 )

# What do I see for AY.4?

ay4 <- read.csv("Documents/GitHub/AppliedStatisticsProject2022/data/lineages/AY.4/binmatrix_ay4.csv")
ay4shuffle <- ay4[sample(nrow(ay4)),]
ay4shuffle <- ay4[sample(nrow(ay4)),]

ay4_topk <- ay4shuffle[1:10000,-c(2,3)]
remove(ay4, ay4shuffle)
gc()

ay4_topk <- ay4_topk[,-1]
distance_mat <- dist(ay4_topk, method = 'euclidean')
Hierar_cl <- hclust(distance_mat, method = "ward.D")
X11()
plot(Hierar_cl)
abline(h = 110, col = "green")
fit <- cutree(Hierar_cl, k = 3 )

# One-way MANOVA
fit <- manova(as.matrix(iris4) ~ species.name)
summary.manova(fit,test="Wilks")