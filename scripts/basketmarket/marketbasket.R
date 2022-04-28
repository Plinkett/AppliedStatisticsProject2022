# Association rule mining in R
# https://towardsdatascience.com/association-rule-mining-in-r-ddf2d044ae50 mostly stolen from here
# https://stackoverflow.com/questions/41303266/r-arules-how-to-remove-certain-itemsets-from-lhs-rhs
# https://www.geeksforgeeks.org/association-rule-mining-in-r-programming/
# https://www.kirenz.com/post/2020-05-14-r-association-rule-mining/
# https://www.datacamp.com/community/tutorials/market-basket-analysis-r

library(arules)
library(arulesViz)

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

data("Groceries")
class(Groceries)

# The inspect() method helps us explore the format of the dataset
# Each row is a transaction
inspect(head(Groceries,2))

# Support -> In how many transactions does the item appear? 
# Confidence ->
grocery_rules <- apriori(Groceries, parameter = list(support = 0.01, confidence = 0.5))

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

wholemilk_rules <- apriori(data=Groceries, parameter=list (supp=0.001,conf = 0.08), appearance = list (rhs="whole milk"))


wholemilk_rules <- apriori(data=Groceries, parameter=list (supp=0.01,conf = 0.5), appearance = list (rhs="whole milk"))
# We get the top 10 rules (by confidence) ending in wholemilk
inspect(head(sort(wholemilk_rules, by = "confidence"),10))
# Higher confidence => stronger rules

# How can we plot this?
itemFrequencyPlot(Groceries, topN = 10)




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










