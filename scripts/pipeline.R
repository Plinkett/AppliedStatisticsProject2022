# library
library(plyr)
library(arules)
library(arulesViz)


#########################
# PCA on Lineage

# Cambiare il nome del lineage sul quale stai lavorando 
lineage <- read.csv("~/Documents/GitHub/AppliedStatisticsProject2022/data/lineages/B.1.1.7/binmatrix_b117.csv")
lineageshuffle <- lineage[sample(nrow(lineage)),]
lineageshuffle <- lineageshuffle[sample(nrow(lineage)),]
lineage_topk <- lineageshuffle[1:10000,-c(2,3)]
remove(lineage, lineageshuffle)
gc()

# Perform PCA 
forpc <- lineage_topk[,-1]
pc <- princomp(forpc, score =T)
summary(pc)

# Scree plot of PCs
X11()
plot(cumsum(pc$sd^2)/sum(pc$sd^2), type='b', axes=F, xlab='number of components', 
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=c(0.8, 0.5), lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(lineage_topk),labels=1:ncol(lineage_topk),las=2)

# Plot of the scores (pairs)
scores <- pc$scores[,1:5]
# If you want the non-noisy data skip this for loop
for(i in 1:ncol(scores))
  scores[,i] <- scores[,i] + rnorm(nrow(scores), sd=0.01)
X11()
pdf(file="~/test.pdf")
par(mar=c(10,5,10,5))
pairs(scores, main ="Pairs plot PCA with jitter (sd = 0.05) - lineage")
dev.off()
# hierarchical clustering
lineage_clust <- lineage_topk[,-1]
# Compute distance matrix
distance_mat <- dist(lineage_clust, method = 'euclidean')
# Compute hierchical clustering with ward linkage
Hierar_cl <- hclust(distance_mat, method = "ward.D")
X11()
# The results are promising and it seems that there are 5 well defined groups. 
# Is this result consistent?
plot(Hierar_cl, xlab="sample", ylab="height", main = "Lineage Dendrogram")
rect.hclust(Hierar_cl, k=8)

fit <- cutree(Hierar_cl, k = 8)
# Replot the pairs plot and color using clustering criterion
scores <- pc$scores[,1:5]
# If you want the non-noisy data skip this for loop
for(i in 1:ncol(scores))
  scores[,i] <- scores[,i] + rnorm(nrow(scores), sd=0.05)
X11()
pairs(scores, col = fit, main ="Pairs plot PCA with jitter (sd = 0.05) - lineage")

# Why are the clusters different?
# simple frequency computation function
mutation_frequencies_bis <- function(binmatrix) {
  freq <-  as.data.frame(matrix(nrow=ncol(binmatrix),ncol=2))
  colnames(freq) <- c("mutation", "frequency")
  for(i in 1:ncol(binmatrix)) {
    freq[i,1] <- colnames(binmatrix)[i]
    freq[i,2] <- sum(binmatrix[,i] == 1)/nrow(binmatrix)
  }
  return(freq)
}

freq <- mutation_frequencies_bis(forpc)
freq <- freq[order(freq$frequency, decreasing = T),] 

X11()
par(mar=c(10,5,10,5))
barplot(height = freq$frequency, names = freq$mutation, las = 2, ylab="mutation percentage", main="Barplot of mutation frequencies in lineage")

sum(fit == 1)  
sum(fit == 2)
sum(fit == 3)
sum(fit == 4)
sum(fit == 5)
sum(fit == 6)
sum(fit == 7)
sum(fit == 8)

# Barplot for cluster 1
cluster1 <- forpc[fit == 1,]
cluster1
freq1 <- mutation_frequencies_bis(cluster1)
freq1 <- freq1[match(freq$mutation, freq1$mutation),]
#freq1 <- freq1[order(freq1$frequency, decreasing = T),] 
X11()
par(mar=c(10,5,10,5))
barplot(height = freq1$frequency, names = freq1$mutation, las = 2, ylab="mutation percentage", main="Barplot cluster 1 in lineage")

# Barplot for cluster 3
cluster3 <- forpc[fit == 3,]
cluster3
freq3 <- mutation_frequencies_bis(cluster3)
freq3 <- freq3[match(freq$mutation, freq3$mutation),]
#freq1 <- freq1[order(freq1$frequency, decreasing = T),] 
X11()
par(mar=c(10,5,10,5))
barplot(height = freq3$frequency, names = freq3$mutation, las = 2, ylab="mutation percentage", main="Barplot cluster 3 in lineage")

# Association rule mining on the interesting mutations seen on clusters
# dataframe for creating transactions
# Extract randomly 10k samples and extract rules
lineage <- read.csv("~/Documents/GitHub/AppliedStatisticsProject2022/data/lineages/B.1.1.7/binmatrix_b117.csv")
lineageshuffle <- lineage[sample(nrow(lineage)),]
lineageshuffle <- lineageshuffle[sample(nrow(lineage)),]
lineage_topk <- lineageshuffle[1:10000,-c(2,3)]
remove(lineage, lineageshuffle)
gc()

lineage_topk <- lineage_topk[,-1]
df_for_rules <- matrix(nrow=nrow(lineage_topk),ncol=ncol(lineage_topk))
names <- colnames(lineage_topk)

for(i in 1:nrow(lineage_topk)) {
  for(j in 1:length(names)) {
    if(lineage_topk[i,j] == 1)
      df_for_rules[i,j] = names[j] 
  }
}

df_for_rules <- as.data.frame(df_for_rules)
write.table(df_for_rules, file="~/transactions_10k_lineage.csv", row.names = FALSE, col.names = FALSE, sep=",")

transactions <- read.transactions("~/transactions_10k_lineage.csv", sep=",", rm.duplicates=TRUE)
rules <- apriori(transactions, parameter = list(supp = 0.001, conf = 0.7, maxlen = 3))
rules <- sort(subset(rules, subset = (lift >= 1.1 | lift <= 0.9)), by="lift")
View(inspect( subset( rules, subset = rhs %pin% "ORF1a.Q3966R" )))
