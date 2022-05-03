# We already have the binary matrix and such
# For association rule mining we will filter the mutations that appear more than
# 90% of the time and those that appear less than 5% of the time.
library(plyr)
library(arules)
library(arulesViz)

filterXtoYpercent <- function(frequencies, totalsamples, x, y) {
  frequencies[,3] <- frequencies[,2] / totalsamples
  frequencies <- subset(frequencies, frequencies[,3] > x & frequencies[,3] < y)
  frequencies <- frequencies[,-3]
  frequencies <- as.data.frame(frequencies)
  frequencies[,1] <- gsub("\\:", "\\.", frequencies[,1])
  return(frequencies)
}


dfLineage <- subset(dfcovid, dfcovid$pango_lineage == "AY.3")
freq_for_rules <- filterXtoYpercent(frequencies, nrow(binmatrix), 0.01, 0.95)
attach(freq_for_rules)
#freq_for_rules <- freq_for_rules[order(-t.freq),]
detach(freq_for_rules)

mutationsToKeep <- as.data.frame(freq_for_rules[,1])
tempbin <- binmatrix[,-1:-3]
totalMutations <- as.data.frame(names(tempbin))
mutationsToDrop <- setdiff(t(totalMutations), t(mutationsToKeep))
mutationsToDrop <- as.character(mutationsToDrop)
tempbin <- tempbin[,!(names(tempbin) %in% mutationsToDrop)]
df_for_rules <- matrix(nrow=nrow(binmatrix),ncol=nrow(freq_for_rules))
names <- colnames(tempbin)

for(i in 1:nrow(binmatrix)) {
  for(j in 1:length(names)) {
    if(tempbin[i,j] == 1)
      df_for_rules[i,j] = names[j] 
  }
}

df_for_rules <- as.data.frame(df_for_rules)

write.table(df_for_rules,file="testing.csv",row.names = FALSE, col.names=FALSE, sep=",")
transactions <- read.transactions("testing.csv", sep=",", rm.duplicates=TRUE)
itemFrequencyPlot(transactions)
rules <- apriori(transactions, parameter=list(supp = 0.01, conf = 0.5))
inspect(rules)
# ON D.2!!!!
X11()
forPlot <- rev(tail(sort(itemFrequency(transactions)),40))
barplot(forPlot, las=2, cex.names=0.8)

# ORF1b:K709N

rules <- apriori(data = transactions, parameter =list(sup=0.01, conf=0.7))
