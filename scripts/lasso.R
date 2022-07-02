library(glmnet)

# This approach DOESN'T WORK because of non gaussianity of data and discreteness

q3 <- read_csv("Documents/GitHub/AppliedStatisticsProject2022/data/lineages/Q.3/binmatrix_q3.csv")
View(binmatrix_q3)


# response variable
y <- q3$N.D3L
x <- data.matrix(q3[,c("N.R203K","N.G204R","N.S235F", "ORF1a.T1001I", "ORF1a.A1708D",
                       "ORF1a.I2230T", "ORF1b.P314L", "ORF1b.D2499Y", "ORF3a.W131C",
                       "ORF8.R52I", "ORF8.Q27.", "ORF8.R52I", "ORF8.K68.", "ORF8.Y73C",
                       "S.N501Y", "S.A570D", "S.D614G", "S.P681H", "S.T716I", "S.S982A",
                       "S.A1020S", "S.D1118H", "M.V70L", "N.Q389H", "ORF1a.M1586I",
                       "ORF1b.K1383R")])
# Perform cross validation on alpha (10 fold)
# x <- scale(x, center = T, scale = T)
# y <- scale(y, center = T, scale = T)
cv_model <- cv.glmnet(x, y, alpha = 1)
best_lambda <- cv_model$lambda.min
best_lambda
# Look at the plot of lambda
plot(cv_model) 

# We fit the model with that lambda
best_lambda
best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
# Coefficients
coefs <- coef(best_model)
coefs

# Extract significant features
summary(best_model)
betas <- as.data.frame(as.matrix(best_model$beta))
betas


rownames(betas)[which(betas$s0 != 0)] 

fit <- lm(N.D3L ~ N.R203K + N.G204R + N.S235F + ORF1a.T1001I + ORF3a.W131C
          + ORF8.Q27. + S.D614G + S.A1020S + M.V70L + N.Q389H + ORF1a.M1586I
          + ORF1b.K1383R, data = q3)

fit <- lm(N.D3L ~ N.R203K + N.G204R + N.S235F + ORF1a.T1001I + ORF3a.W131C
          + ORF8.Q27. + S.D614G + S.A1020S + M.V70L + ORF1a.M1586I
          + ORF1b.K1383R, data = q3)

fit <- lm(N.D3L ~ N.R203K + N.G204R + N.S235F + ORF1a.T1001I + ORF3a.W131C
          + ORF8.Q27. + S.D614G + S.A1020S + ORF1a.M1586I
          + ORF1b.K1383R, data = q3)

fit <- lm(N.D3L ~ N.R203K + N.G204R + N.S235F + ORF1a.T1001I + ORF3a.W131C
          + S.D614G + S.A1020S + ORF1a.M1586I
          + ORF1b.K1383R, data = q3)

fit <- lm(N.D3L ~ N.R203K + N.G204R + N.S235F
          + S.D614G + S.A1020S + ORF1a.M1586I
          + ORF1b.K1383R, data = q3)

summary(fit)


summary(best_model$beta)
y_predicted <- predict(best_model, s = best_lambda, newx = x)
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)
rsq <- 1 - sse/sst
rsq
