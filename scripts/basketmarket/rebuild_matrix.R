
# Creation of matrix for easy transaction access and rule mining
# In this case we take the two biggest lineages (AY.4 and B.1.1.7) and generate
# a .csv for the 
library(plyr)
library(arules)
library(arulesViz)


# CURRENT PLAN
# Filter common mutations < 95%
# And then generate rules, further filtering for rules with for
# lift < 0.9 or lift > 1.1
# Maximum length? 5
# Perform chi square test on pairs



filterXtoYpercent <- function(frequencies, totalsamples, x, y) {
  frequencies[,3] <- frequencies[,2] / totalsamples
  frequencies <- subset(frequencies, frequencies[,3] > x & frequencies[,3] < y)
  frequencies <- frequencies[,-3]
  frequencies <- as.data.frame(frequencies)
  frequencies[,1] <- gsub("\\:", "\\.", frequencies[,1])
  return(frequencies)
}

calculate_frequencies_binmatrix <- function(binmatrix) {
    frequencies <- matrix(nrow = (ncol(binmatrix) - 3), ncol = 2)
    for(i in 4:ncol(binmatrix)) {
      frequencies[i-3,1] <- colnames(binmatrix)[i]
      frequencies[i-3,2] <- sum(binmatrix[,i] == 1)
    }  
    frequencies <- as.data.frame(frequencies)
    frequencies[,2] <- as.numeric(frequencies[,2])
    return(frequencies)
}
setwd("~/Documents/GitHub/AppliedStatisticsProject2022/data/lineages")
binmatrix <- read.csv("B.1.1.7/binmatrix_b117.csv", sep=",")
frequencies <- calculate_frequencies_binmatrix(binmatrix)
freq_for_rules <- filterXtoYpercent(frequencies, nrow(binmatrix), 0.005, 0.95)

attach(freq_for_rules)
#freq_for_rules <- freq_for_rules[order(-t.freq),]
detach(freq_for_rules)

mutationsToKeep <- as.data.frame(freq_for_rules[,1])
binmatrix <- binmatrix[,-1:-3]
totalMutations <- as.data.frame(names(binmatrix))
mutationsToDrop <- setdiff(t(totalMutations), t(mutationsToKeep))
mutationsToDrop <- as.character(mutationsToDrop)
binmatrix <- binmatrix[,!(names(binmatrix) %in% mutationsToDrop)]
df_for_rules <- matrix(nrow=nrow(binmatrix),ncol=ncol(binmatrix))
names <- colnames(binmatrix)



for(i in 1:nrow(binmatrix)) {
  for(j in 1:length(names)) {
    if(binmatrix[i,j] == 1)
      df_for_rules[i,j] = names[j] 
  }
}

df_for_rules <- as.data.frame(df_for_rules)
write.table(df_for_rules, file="~/final_rules.csv", row.names = FALSE, col.names = FALSE, sep=",")
remove(list = ls())
gc()
# Just checking if it works
df_for_rules <- read.csv("~/final_rules.csv")









###################################################
# Rebulding balanced dataset
###################################################

# b117, ay4, ay103, ba1, ay3
# 100k per lineage

b117 <- read_csv("B.1.1.7/binmatrix_b117.csv")
b117 <- b117[sample(nrow(b117)),]
b117 <- b117[1:10000,]

ay4 <- read_csv("AY.4/binmatrix_ay4.csv")
ay4 <- ay4[sample(nrow(ay4)),]
ay4 <- ay4[1:10000,]

ba1 <- read_csv("BA.1/binmatrix_ba1.csv")
ba1 <- ba1[sample(nrow(ba1)),]
ba1 <- ba1[1:10000,]

ay103 <- read_csv("AY.103/binmatrix_ay103.csv")
ay103 <- ay103[sample(nrow(ay103)),]
ay103 <- ay103[1:100000,]

binmatrix <- create_mutations_set_lineages(b117, ay4)
binmatrix <- create_mutations_set_lineages(binmatrix, ba1)
binmatrix <- create_mutations_set_lineages(binmatrix, ay103)
remove(b117, ay4, ba1, ay103)
gc()
binmatrix <- binmatrix[,-1:-3]

