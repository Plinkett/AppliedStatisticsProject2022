library(arules)
library(arulesViz)
#https://stackoverflow.com/questions/7714677/scatterplot-with-too-many-points
#https://stackoverflow.com/questions/38962015/10-000-point-3d-scatter-plots-in-python-with-quick-rendering
setwd("~/Documents/GitHub/AppliedStatisticsProject2022/data/lineages")
transactions <- read.transactions("AY.4/transaction_ay4_005_95_filter.csv", sep=",", rm.duplicates=TRUE)
summary(transactions)

items = rev(tail(sort(itemFrequency(transactions)), 91))
X11()
barplot(items, las=2, cex.names=0.8, col="gold")
dev.off()

rules <- apriori(transactions, parameter = list(supp = 0.0005, conf = 0.8, maxlen = 3))
# We want rules with lift >= 1.1
rules <- sort(subset(rules, subset = (lift >= 1.1 | lift <= 0.9)), by="lift") 
rules

write(rules, file = "~/Documents/GitHub/AppliedStatisticsProject2022/scripts/basketmarket/rules_testing/rules.csv",
      sep = ",", quote = TRUE, row.names = FALSE)
rules <- read.csv("~/Documents/GitHub/AppliedStatisticsProject2022/scripts/basketmarket/rules_testing/rules.csv")

plot_ly(x=rules$suport, y=rules$confidence, z=rules$lift, type="scatter3d", mode="ù
        markers")

# Two-key plot
X11()
plot(rules, method = "two-key plot")
dev.off()

X11()
plot(rules)
dev.off()


set.seed(417)
library(plotly)
temp <- rnorm(100, mean=30, sd=5)
pressure <- rnorm(100)
dtime <- 1:100
plot_ly(x=rules$support, y=rules$confidence, z=rules$lift, type="scatter3d", mode="markers")
