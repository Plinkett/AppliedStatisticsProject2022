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

# Strategy for generating rules
# Don't filter and use high support, rules will be overwhelmingly about the most common
# mutations. There will also be many rules (limit length to 3 or 4). Filter the rules by lift
# 1.1 might be a good choice. Keep those, I went from 100k to 2k rules.
# Then reconsider the transactions and filter out the most common mutations (above 80%). 
# Now generate rules from those transactions. Put the rules together 

#--------------------------------------------------------------------------------------------
# Build and save rules for AY.4
transactions <- read.transactions("AY.4/transaction_ay4_005_95_filter.csv", sep=",", rm.duplicates=TRUE)
summary(transactions)
items = rev(tail(sort(itemFrequency(transactions)), 91))

X11()
barplot(items, las=2, cex.names=0.8, col="gold")
dev.off()
# FINAL PARAMETERS ON OVERALL MATRIX!!!!!
# Minimum support of 0.01
# Minimum confidence of 0.8
# Maximum length of 4
rules <- apriori(transactions, parameter = list(supp = 0.005, conf = 0.8, maxlen = 4))

# We want rules with lift >= 1.1
rulesay4 <- sort(subset(rules, subset = (lift >= 1.1 | lift <= 0.9)), by="lift") 
rulesay4

write(rulesay4, file = "AY.4/rules_ay4_filtered.csv", sep = ",", quote = TRUE, row.names = FALSE)
write.PMML(rulesay4, file = "AY.4/rules_ay4.xml")

# Same thing for B.1.1.7
transactions <- read.transactions("B.1.1.7/transaction_b117_no_filter.csv", sep=",", rm.duplicates=TRUE)
items = rev(tail(sort(itemFrequency(transactions)), 128))
X11()
barplot(items, las=2, cex.names=0.8, col="gold")
dev.off()
rules <- apriori(transactions, parameter = list(supp = 0.005, conf = 0.8, maxlen = 4))
rulesb117 <- sort(subset(rules, subset = (lift >= 1.1 | lift <= 0.9)), by="lift") 
rulesb117

write(rulesb117, file = "B.1.1.7/rules_b117_filtered.csv", sep = ",", quote = TRUE, row.names = FALSE)
write.PMML(rulesb117, file = "B.1.1.7/rules_b117.xml")

# BA.1
transactions <- read.transactions("BA.1/transactions_ba1.csv", sep=",", rm.duplicates=TRUE)
items = rev(tail(sort(itemFrequency(transactions)), 128))
X11()
barplot(items, las=2, cex.names=0.8, col="gold")
dev.off()


rules <- apriori(transactions, parameter = list(supp = 0.005, conf = 0.8, maxlen = 4))
write(rules, file = "BA.1/rules_ba1.csv", sep = ",", quote = TRUE, row.names = FALSE)
write.PMML(rules, file = "BA.1/rules_ba1.xml")

# AY.103
transactions <- read.transactions("AY.103/transaction_ay103_no_filter.csv", sep=",", rm.duplicates=TRUE)
items = rev(tail(sort(itemFrequency(transactions)), 128))
X11()
barplot(items, las=2, cex.names=0.8, col="gold")
dev.off()

rules <- apriori(transactions, parameter = list(supp = 0.005, conf = 0.8, maxlen = 4))
rulesay103 <- sort(subset(rules, subset = (lift >= 1.1 | lift <= 0.9)), by="lift") 
rulesay103

write(rulesay103, file = "AY.103/rules_ay103_filtered.csv", sep = ",", quote = TRUE, row.names = FALSE)
write.PMML(rulesay103, file = "AY.103/rules_ay103.xml")

rulesay103 <- read.csv("AY.103/rules_ay103_filtered.csv", sep=",")
rulesb117 <- read.csv("B.1.1.7/rules_b117_filtered.csv", sep=",")
rulesay4 <- read.csv("AY.4/rules_ay4_filtered.csv", sep=",")
rulesTot <- rbind(rulesay4,rulesb117,rulesay103)

write.csv(rulesTot, file = "rules_ay4_b117_ay103.csv", row.names = FALSE)
rulesTot <- read.csv("rules_ay4_b117_ay103.csv")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PLOT FOR RULES TOT

support <- jitter(rulesTot$support,2)
confidence <- jitter(rulesTot$confidence, 2)
lift <- jitter(rulesTot$lift, 2)

X11()
plot(support, log(lift), pch=21)

#-------------------------------------------------------------------------------------------
# Density is ~ 0.04, this means that the matrix is sparse
summary(transactions)
items = rev(tail(sort(itemFrequency(transactions)), 100))
X11()
barplot(items, las=2, cex.names=0.8, col="gold")
dev.off()
# For AY.103 -> take out S.G142D
# For B.1.1.7 -> take out N.G204R
rules <- apriori(transactions, parameter = list(supp = 0.9, conf = 0.9, maxlen = 3))
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


#-------------------------------------------------------------------------------
# FOR CHECKING VALIDITY OF RULES
# We read the rules and clean them

# INIZIA QUAAAAAAAA

setwd("~/Documents/GitHub/AppliedStatisticsProject2022/data/lineages")
rulescsv <- read.csv("rules_ay4_b117_ay103.csv")
head(rulescsv)

clean_rules <- function(rule) {
  rule <- gsub('[{}]', '', rule)
  rule <- gsub(' => ', ',', rule)
  return(rule)
}

result <- as.data.frame(lapply(rulescsv[,1], clean_rules))
result <- as.data.frame(t(result))
rulescsv[,"filtered"] <- result[,1]

toRemove <- numeric()
for(i in 1:nrow(rulescsv)) {
  rule_str <- rulescsv[i,1]
  if(substring(rule_str, 1, 2) == "{}")
     toRemove <- append(toRemove, i) 
}

if(length(toRemove) != 0)
  rulescsv <- rulescsv[-toRemove,]
#Cleaning
remove(result,rules)
# We could start, for example, with AY.44
# Then we will put this inside a for loop

# CAMBIATE IL LINEAGE A QUELLO CHE VOLETE
test_lineage <- read.csv("B.1.2/binmatrix_b12.csv")
# First single out applicable rules
# TO MODIFY!!!!!
app_rules <- as.data.frame(matrix(nrow = 0, ncol = ncol(rulescsv)))
app_rules_none <- as.data.frame(matrix(nrow = 0, ncol = ncol(rulescsv)))

for(i in 1:nrow(rulescsv)) {
  tokens <- strsplit(rulescsv$filtered[i], ",") 
  applicable <- TRUE
  if(!(tokens[[1]][1] == "")) {
    for(j in 1:(length(tokens[[1]])-1))
      if(!(tokens[[1]][j] %in% colnames(test_lineage))) 
        applicable <- FALSE
    if(applicable == TRUE)
      if(tokens[[1]][length(tokens[[1]])] %in% colnames(test_lineage)) 
        app_rules <- rbind(app_rules, rulescsv[i,])
      else
        app_rules_none <- rbind(app_rules_none, rulescsv[i,])
  }
}

nrules <- nrow(app_rules)
app_rules <- rbind(app_rules, app_rules_none)
app_rules[,"count_in_lineage"] <- 0
app_rules[,"success_rate"] <- 0

for(i in 1:nrules) {
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
}
app_rules <- app_rules[,-7]

# SALVATE QUAAAAAAAAAAAAAA
write.csv(data.frame(app_rules), file = "B.1.2/apprules_b12.csv",row.names=FALSE)

# Swapping and putting ALL other mutations in the rhs.
# See success rates.


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