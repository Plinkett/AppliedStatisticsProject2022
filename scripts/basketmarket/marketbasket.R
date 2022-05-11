# Association rule mining in R
# https://towardsdatascience.com/association-rule-mining-in-r-ddf2d044ae50 mostly stolen from here
# https://stackoverflow.com/questions/41303266/r-arules-how-to-remove-certain-itemsets-from-lhs-rhs
# https://www.geeksforgeeks.org/association-rule-mining-in-r-programming/
# https://www.kirenz.com/post/2020-05-14-r-association-rule-mining/
# https://www.datacamp.com/community/tutorials/market-basket-analysis-r

# For 3D plots (hopefully)
# https://rpubs.com/swagatamdas03/507581
# https://irapoenya.wordpress.com/2020/10/04/rstudio-tutorial-an-intro-to-3d-plots/
# https://rdrr.io/rforge/arulesViz/man/plot.html
# https://github.com/mhahsler/arulesViz
# https://towardsdatascience.com/association-rules-2-aa9a77241654
# 
library(arules)
library(arulesViz)

# FOR A SINGLE LINEAGE
# !!!!!! Eseguire quando la directory di partenza è /data/lineages!
transactions <- read.transactions("AY.103/transactions_ay103.csv", sep=",", rm.duplicates=TRUE)
transactions <- read.transactions("B.1.1.7/transactions_b117.csv", sep=",", rm.duplicates=TRUE)
# Test with the big 4 lineages into one transaction file
transactions <- read.transactions("transaction_ay4_b117_ba1_ay103.csv", sep=",", rm.duplicates=TRUE)

#--------------------------------------------------------------------------------------------
# Build and save rules for AY.4
transactions <- read.transactions("AY.4/transactions_ay4.csv", sep=",", rm.duplicates=TRUE)
summary(transactions)
items = rev(tail(sort(itemFrequency(transactions)), 15))
X11()
barplot(items, las=2, cex.names=0.8, col="gold")
dev.off()
rules0001ay4 <- apriori(transactions, parameter = list(supp = 0.001, conf = 0.5))
rules001ay4 <- apriori(transactions, parameter = list(supp = 0.01, conf = 0.5))
# Exclude "ORF7.aV82A", "S.T95I", "S.G142D" from rules
rules <- apriori(transactions, parameter = list(supp = 0.001, conf = 0.5), 
                 appearance = list(none = c("ORF7a.V82A", "S.T95I", "S.G142D"), default = "both"))
# Very few rules...
# So: consider mutations > 95% and < 5%, BUT consider a support of 0.1%
write(rules0001ay4, file = "AY.4/rules_ay4_0001.csv", sep = ",", quote = TRUE, row.names = FALSE)
write(rules001ay4, file = "AY.4/rules_ay4_001.csv", sep = ",", quote = TRUE, row.names = FALSE)
write.PMML(rules0001ay4, file = "AY.4/rules_ay4_0001.xml")
write.PMML(rules001ay4, file = "AY.4/rules_ay4_001.xml")

# Same thing for B.1.1.7
transactions <- read.transactions("B.1.1.7/transactions_b117.csv", sep=",", rm.duplicates=TRUE)
rules0001b117 <- apriori(transactions, parameter = list(supp = 0.001, conf = 0.5))
rules001b117 <- apriori(transactions, parameter = list(supp = 0.01, conf = 0.5))
write(rules0001b117, file = "B.1.1.7/rules_b117_0001.csv", sep = ",", quote = TRUE, row.names = FALSE)
write(rules001b117, file = "B.1.1.7/rules_b117_001.csv", sep = ",", quote = TRUE, row.names = FALSE)
write.PMML(rules0001b117, file = "B.1.1.7/rules_b117_0001.xml")
write.PMML(rules001b117, file = "B.1.1.7/rules_b117_001.xml")
# BA.1
transactions <- read.transactions("BA.1/transactions_ba1.csv", sep=",", rm.duplicates=TRUE)
rules0001ba1 <- apriori(transactions, parameter = list(supp = 0.001, conf = 0.5, maxlen = 15))
rules001ba1 <- apriori(transactions, parameter = list(supp = 0.01, conf = 0.5, maxlen = 15))
write(rules0001ba1, file = "BA.1/rules_ba1_0001.csv", sep = ",", quote = TRUE, row.names = FALSE)
write(rules001ba1, file = "BA.1/rules_ba1_001.csv", sep = ",", quote = TRUE, row.names = FALSE)
write.PMML(rules0001ba1, file = "BA.1/rules_ba1_0001.xml")
write.PMML(rules001ba1, file = "BA.1/rules_ba1_001.xml")

# AY.103
transactions <- read.transactions("AY.103/transactions_ay103.csv", sep=",", rm.duplicates=TRUE)
rules0001ay103 <- apriori(transactions, parameter = list(supp = 0.001, conf = 0.5))
rules001ay103 <- apriori(transactions, parameter = list(supp = 0.01, conf = 0.5))
write(rules0001ay103, file = "AY.103/rules_ay103_0001.csv", sep = ",", quote = TRUE, row.names = FALSE)
write(rules001ay103, file = "AY.103/rules_ay103_001.csv", sep = ",", quote = TRUE, row.names = FALSE)
write.PMML(rules0001ay103, file = "AY.103/rules_ay103_0001.xml")
write.PMML(rules001ay103, file = "AY.103/rules_ay103_0001.xml")

rules0001tot <- c(rules0001ba1, rules0001ay103, rules0001b117, rules0001ay4)
rules001tot <- c(rules001ba1, rules001ay103, rules001b117, rules001ay4)
rules <- rules001tot
# SAVE ALL RULES

write(rules0001tot, file = "ay4_b117_ba1_ay103_0001_rules.csv", sep = ",", quote = TRUE, row.names = FALSE)
write.PMML(rules0001tot, file = "ay4_b117_ba1_ay103_0001_rules.xml")
write(rules001tot, file = "ay4_b117_ba1_ay103_001_rules.csv", sep = ",", quote = TRUE, row.names = FALSE)
write.PMML(rules001tot, file = "ay4_b117_ba1_ay103_001_rules.xml")



test <- read.PMML("AY.103/rules_ay103_001.csv")

# EXPERIMENTING WITH PLOTS
#-------------------------------------------------------------------------------------------
# Density is ~ 0.04, this means that the matrix is sparse
summary(transactions)
items = rev(tail(sort(itemFrequency(transactions)), 15))
X11()
barplot(items, las=2, cex.names=0.8, col="gold")
dev.off()
# For AY.103 -> take out S.G142D
# For B.1.1.7 -> take out N.G204R
rules <- apriori(transactions, parameter = list(supp = 0.01, conf = 0.5, maxlen = 15))
inspect(rules)
# Avoid mutations by using "none" attribute, I'm excluding the most frequent mutation(s).
# For B.1.1.7
rules <- apriori(transactions, parameter = list(supp = 0.01, conf = 0.5, maxlen = 15), 
                 appearance = list(none = c("N.G204R"), default = "both"))
# For AY.103
rules <- apriori(transactions, parameter = list(supp = 0.001, conf = 0.5, maxlen = 15), 
                 appearance = list(none = c("S.G142D"), default = "both"))

# I don't quite know how to read this...
X11()
plot(rules, method = "graph", measure = "confidence", shading = "lift")
dev.off()

# Scatter plot, hue is chosen by lift (most significant I would say)
X11()
plot(rules)
dev.off()

# Scatter plot, hue is chosen by confidence (could be misleading)
X11()
plot(rules, measure = "confidence")
dev.off()

# Two-key plot
X11()
plot(rules, method = "two-key plot")
dev.off()

# Interactive two-key plot
X11()
plot(rules, engine = "plotly")
dev.off()

# Parallel coordinate plot
X11()
plot(rules, method="paracoord")
dev.off()

# Graph network
subrules <- head(rules, n = 36, by = list("lift","confidence"))
inspect(subrules)
X11()
plot(rules, method = "graph",  engine = "htmlwidget")
graphics.off()





#----------------------------------------------------------------------------
# FOR HANDLING THE 4 BIG LINEAGES ALL AT ONCE
#----------------------------------------------------------------------------
transactions <- read.transactions("transaction_ay4_b117_ba1_ay103.csv", sep=",", rm.duplicates=TRUE)
# With a support of 0.001 (over 1.5 mln samples) we consider 
rules <- apriori(transactions, parameter = list(supp = 0.001, conf = 0.5, maxlen = 15)) 
# Let's see the most common mutations
items = rev(tail(sort(itemFrequency(transactions)), 15))
X11()
barplot(items, las=2, cex.names=0.8, col="gold")
dev.off()

# Now we ignore mutations that appear more than 20% of the time
rules <- apriori(transactions, parameter = list(supp = 0.007, conf = 0.5, maxlen = 15), 
                 appearance = list(none = c("S.G142D","ORF7a.V82A","S.T95I"), default = "both"))

write(rules, file = "rules_filtered.csv", sep = ",", quote = TRUE, row.names = FALSE)

# Scatter plot, hue is chosen by lift (most significant I would say)
X11()
plot(rules)
dev.off()

# Scatter plot, hue is chosen by confidence (could be misleading)
X11()
plot(rules, measure = "confidence")
dev.off()

# Two-key plot
X11()
plot(rules, method = "two-key plot")
dev.off()

# Interactive two-key plot
X11()
plot(rules, engine = "plotly")
dev.off()

# Parallel coordinate plot
X11()
plot(rules, method="paracoord")
dev.off()

# Graph network
subrules <- head(rules, n = 36, by = list("lift","confidence"))
inspect(subrules)
X11()
plot(rules, method = "graph",  engine = "htmlwidget")
graphics.off()

#-------------------------------------------------------------------------------
# FOR CHECKING VALIDITY OF RULES
# We read the rules 
rulescsv <- read.csv("ay4_b117_ba1_ay103_001_rules.csv")
head(rulescsv)

clean_rules <- function(rule) {
  rule <- gsub('[{}]', '', rule)
  rule <- gsub(' => ', ',', rule)
  return(rule)
}
result <- t(as.data.frame(lapply(rulescsv[,1], clean_rules)))
rulescsv[,"filtered"] <- result[,1]

#Cleaning
remove(result,rules)
# We could start, for example, with AY.3
# Then we will put this inside a for loop
setwd("~/Documents/GitHub/AppliedStatisticsProject2022/data/lineages")
test_lineage <- read.csv("B.1.2/binmatrix_b12.csv")
# First single out applicable rules
# 
app_rules <- as.data.frame(matrix(nrow = 0, ncol = ncol(rulescsv)))
for(i in 1:nrow(rulescsv)) {
  tokens <- strsplit(rulescsv$filtered[i], ",") 
  applicable <- TRUE
  if(!(tokens[[1]][1] == "")) {
    for(j in 1:length(tokens[[1]]))
      if(!(tokens[[1]][j] %in% colnames(test_lineage))) 
        applicable <- FALSE
    if(applicable == TRUE)
      app_rules <- rbind(app_rules, rulescsv[i,])
  }
}

app_rules[,"count_in_lineage"] = NA
app_rules[,"success_rate"] = NA

for(i in 1:nrow(app_rules)) {
  tokens <- strsplit(app_rules$filtered[i], ",") 
  lhs_length <- length(tokens[[1]]) - 1
  rhs_mutation <- tokens[[1]][lhs_length]
  # Since we can't have a variable number of conditions we will just use
  # subsets...
  # Initialization of subset 
  lhsIndices <- vector(mode="numeric", length=lhs_length)
  for(j in 1:lhs_length) {
    lhsIndices[j] <- which(colnames(test_lineage) == tokens[[1]][j])
  }
  lhs_subset <- subset(test_lineage, test_lineage[,lhsIndices[1]] == 1)
  rhs <- which(colnames(test_lineage) == tokens[[1]][lhs_length + 1])
  if(lhs_length > 1) {
    for(j in 2:lhs_length) {
      lhs_subset <- subset(lhs_subset, lhs_subset[,lhsIndices[j]] == 1)
    }
  }
  app_rules$count_in_lineage[i] <- nrow(lhs_subset)
  app_rules$success_rate[i] <- sum(lhs_subset[,rhs] == 1) / nrow(lhs_subset)
  # print(nrow(lhs_subset))
  #print(sum(lhs_subset[,rhs] == 1))
  #print(sum(lhs_subset[,rhs] == 0))
}
app_rules <- app_rules[,-7]
write.csv(data.frame(app_rules), file = "AY.103/apprules_ay103_001.csv",row.names=FALSE)


test <- rulescsv[10000,7]
test
tokens <- strsplit(test, ",") 
tokens
length(tokens[[1]])
rhs_mutation
# Rhs mutation
# 
which(colnames(test_lineage)=="date")
sum(test_lineage$M.I82T == 1)







#-------------------------------------------------------------------------------


# RANDOM NOTES
# Apriori algorithm, extracts association rules exploiting the idea that
# the relative support of X has to be more than or equal than the relative
# support of a set Y that includes X. On the other hand the relative support
# of a set Y that includes X is less than or equal to that of X.

# This realization helps us prune branches from our prefix tree and not even
# look at associations that will surely be under our minimum threshold anyway

# What's its output? The association rules found above a minimum threshold 
# defined by us.

# Keep in mind that arules and its algorithms receive as input data in
# the "transactions" format and NOT the dataframe format.

data("groceries")
class(groceries)

# The inspect() method helps us explore the format of the dataset
# Each row is a transaction
inspect(head(groceries,2))

# Support -> In how many transactions does the item appear? 
# Confidence ->
grocery_rules <- apriori(groceries, parameter = list(support = 0.01, confidence = 0.5))

# Let's try to identify read the output of apriori
# minval -> minimum value of the support to even be considere (given by us at the beginning)
# 
inspect(head(sort(grocery_rules, by = "confidence"), 3))

# We can play with the "support" and "confidence" parameters to obtain 
# meaningful rules

# By specifying rhs  (right-hand-side) whole milk we want the rules
# that produce the wholemilk product

# N.B. the one below produces not very meaning rules because of the very low support to begin with
# Why? Suppose I have an item that is almost never bought, its support will be very low. Consider
# that whole milk is the most frequent item so it is very likely that in the few instances  in which
# this infrequent item appears whole milk also appears. This could lead you to generate a rule like this:
# {infrequent item} => {whol milk}, with a VERY high confidence because of the reasons above
# I guess we could remove the most frequent items from the right-hand side

wholemilk_rules <- apriori(data=groceries, parameter=list (supp=0.001,conf = 0.08), appearance = list (rhs="whole milk"))


wholemilk_rules <- apriori(data=groceries, parameter=list (supp=0.01,conf = 0.5), appearance = list (rhs="whole milk"))
# We get the top 10 rules (by confidence) ending in wholemilk
inspect(head(sort(wholemilk_rules, by = "confidence"),10))
# Higher confidence => stronger rules

# How can we plot this?
itemFrequencyPlot(groceries, topN = 10)

# Now let's look at transforming dataframes into transactions
data("AdultUCI")
class(AdultUCI)
View(AdultUCI)
str(AdultUCI)

# In order to transition to transactions all the columns must in the "factor" form
AdultUCI <- lapply(AdultUCI, function(x){as.factor(x)})
str(AdultUCI)

# This still gives us issues...
data <- as(AdultUCI, "transactions")

groceries <- read.transactions("Market_Basket_Optimisation.csv",sep=",",rm.duplicates=TRUE)
str(groceries)
rules <- apriori(data=groceries, parameter=list(support=0.01,confidence=0.2))
itemFrequencyPlot(groceries, topN=10)
inspect(head(sort(rules,by="confidence")))
X11()
plot(rules, method="graph",measure="confidence",shading="lift")
