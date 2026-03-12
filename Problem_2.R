#TO DO:
# DO some EDA
#Do we need to do quadratic effects?
# ======

# Load necessary libraries ====
library(boot)
library(ggplot2)
library(rpart)
library(randomForest)
library(caret)
library(rpart.plot)
library(rattle)
library(pROC)
library(ipred)
library(devtools)
library(readxl)
library(lars)
library(glmnet)
library(plotly)

set.seed(1717)  # For reproducibility

### Restructure Data

covid_dataset <- read_excel("Immunologic profiles of patients with COVID-19.xlsx")
covid_dataset <- rename(covid_dataset, Severity = Severirty)


## Remove missing values
missing_vals_index <- which(!complete.cases((covid_dataset)))
missing_vals_index #Integer(0) so there are no missing values

dim(covid_dataset) #76 x 32
str(covid_dataset)

# Prepare data
covid_metrics <- covid_dataset[,-1]  # Remove patient number
covid_predictors <- model.matrix(Severity ~ ., data = covid_metrics)[, -1] #Remove the intercept (X)
covid_response <- as.factor(covid_dataset$Severity) #(Y)

# Check dimensions
cat("Dataset dimensions:", dim(covid_dataset), "\n")
cat("Number of predictors:", ncol(covid_predictors), "\n")
cat("Number of observations:", length(covid_response), "\n")

#Splitting the data into training (75%) and testing (25%)
train_index <- createDataPartition(covid_response, p = 0.75, list = FALSE)
predictors_train <- covid_predictors[train_index, ]
predictors_test <- covid_predictors[-train_index, ]
response_train <- covid_response[train_index]
response_test <- covid_response[-train_index]

# Standardize predictors (glmnet does this by default with standardize=TRUE)

# Grid search ====
cat("\n=== GRID SEARCH FOR OPTIMAL ALPHA ===\n")

# Define alpha grid
alpha_grid <- seq(0, 1, by = 0.01)

# Create container for results
results_cv10 <- data.frame()
results_cv20 <- data.frame()

#Need to use type.measure = "deviance" because not enough results for auc
# Perform grid search for CV = 10 fold ====
for (alpha_val in alpha_grid) {
  set.seed(1717)
  # Cross-validation for current alpha
  cv_model <- cv.glmnet(predictors_train, response_train, 
                        alpha = alpha_val, 
                        nfolds = 10,
                        type.measure = "deviance",
                        keep = TRUE,
                        family = "binomial", #Severity is either mild or severe
                        standardize = TRUE)
  
  # Store results
  results_cv10 <- rbind(results_cv10, data.frame(
    alpha = alpha_val,
    lambda_min = cv_model$lambda.min,
    lambda_1se = cv_model$lambda.1se,
    cvm_min = min(cv_model$cvm),
    cvm_1se = cv_model$cvm[which(cv_model$lambda == cv_model$lambda.1se)],
    nzero_min = cv_model$nzero[which(cv_model$lambda == cv_model$lambda.min)],
    nzero_1se = cv_model$nzero[which(cv_model$lambda == cv_model$lambda.1se)]
  ))
}

# Perform grid search for CV fold 20 ====
for (alpha_val in alpha_grid) {
  set.seed(1717)
  # Cross-validation for current alpha
  cv_model <- cv.glmnet(predictors_train, response_train, 
                        alpha = alpha_val, 
                        nfolds = 20,
                        type.measure = "deviance",
                        keep = TRUE,
                        family = "binomial", #Severity is either mild or severe
                        standardize = TRUE,
                        grouped = FALSE) #Adding grouped = FALSE removes warnings for CV20 because it creates an error matrix at the observation level of all folds to summarize results
  
  # Store results
  results_cv20 <- rbind(results_cv20, data.frame(
    alpha = alpha_val,
    lambda_min = cv_model$lambda.min,
    lambda_1se = cv_model$lambda.1se,
    cvm_min = min(cv_model$cvm),
    cvm_1se = cv_model$cvm[which(cv_model$lambda == cv_model$lambda.1se)],
    nzero_min = cv_model$nzero[which(cv_model$lambda == cv_model$lambda.min)],
    nzero_1se = cv_model$nzero[which(cv_model$lambda == cv_model$lambda.1se)]
  ))
}

# Print results
print(results_cv10)
print(results_cv20)

# Find best alpha based on minimum CV error
best_alpha_min10 <- results_cv10$alpha[which.min(results_cv10$cvm_min)]
best_alpha_1se10 <- results_cv10$alpha[which.min(results_cv10$cvm_1se)]

best_alpha_min20 <- results_cv20$alpha[which.min(results_cv20$cvm_min)]
best_alpha_1se20 <- results_cv20$alpha[which.min(results_cv20$cvm_1se)]

#Best alpha for CV 10
cat("\nBest alpha CV10 (lambda.min):", best_alpha_min10, "\n")
cat("Best alpha CV10 (lambda.1se):", best_alpha_1se10, "\n")

#Best alpha for CV 20
cat("\nBest alpha CV20 (lambda.min):", best_alpha_min20, "\n")
cat("Best alpha CV20 (lambda.1se):", best_alpha_1se20, "\n")

# Visualize alpha grid search - CV10
ggplot(results_cv10, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), size = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), size = 1.2) +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 3) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 3) +
  labs(title = "CV Error vs Alpha (10 fold)",
       x = "Alpha (0=Ridge, 1=Lasso)",
       y = "CV Error (binomial deviance)",
       color = "Lambda Type") +
  theme_minimal()

# Visualize alpha grid search - CV20
ggplot(results_cv20, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), size = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), size = 1.2) +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 3) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 3) +
  labs(title = "CV Error vs Alpha (20 fold)",
       x = "Alpha (0=Ridge, 1=Lasso)",
       y = "CV Error (binomial deviance)",
       color = "Lambda Type") +
  theme_minimal()

par(mfrow = c(1,2))
plot(results_cv10$alpha, results_cv10$cvm_min,
     main = "CV10: Error vs Alpha",
     xlab = "Alpha", ylab = "CV Error",
     type = "l")
abline(v = best_alpha_min10, col = "red", lty = 2)

plot(results_cv20$alpha, results_cv20$cvm_min,
     main = "CV20: Error vs Alpha",
     xlab = "Alpha", ylab = "CV Error",
     type = "l")
abline(v = best_alpha_min20, col = "red", lty = 2)
par(mfrow = c(1,1))

# Comparing CV10 and CV20 schemes  
# need one representative model from each scheme 
# decided on lambda.min instead of lambda.1se, as it keeps more predictors and prioritizes accuracy over simplicity 

# Extract lambda values
best_row_cv10 <- which.min(results_cv10$cvm_min)
best_lambda_cv10 <- results_cv10$lambda_min[best_row_cv10]

best_row_cv20 <- which.min(results_cv20$cvm_min)
best_lambda_cv20 <- results_cv20$lambda_min[best_row_cv20]

cat("Best lambda (CV10):", best_lambda_cv10, "\n")
cat("Best alpha (CV10):", best_alpha_min10, "\n")
cat("CV10 error:", min(results_cv10$cvm_min), "\n")
cat("CV10 - Non-zero coefficients:", 
    results_cv10$nzero_min[which.min(results_cv10$cvm_min)], "\n")

cat("Best lambda (CV20):", best_lambda_cv20, "\n")
cat("Best alpha (CV20):", best_alpha_min20, "\n")
cat("CV20 error:", min(results_cv20$cvm_min), "\n")
cat("CV20 - Non-zero coefficients:", 
    results_cv20$nzero_min[which.min(results_cv20$cvm_min)], "\n")

# Fit final models
enet_cv10_model_min <- glmnet(predictors_train, response_train,
                              family = "binomial",
                              alpha = best_alpha_min10,
                              lambda = best_lambda_cv10)
enet_cv20_model_min <- glmnet(predictors_train, response_train,
                              family = "binomial",
                              alpha = best_alpha_min20,
                              lambda = best_lambda_cv20)

# CV10 predictions
prds.train.cv10 <- predict(enet_cv10_model_min, newx = predictors_train,
                           type = "response")[,1]
prds.test.cv10 <- predict(enet_cv10_model_min, newx = predictors_test,
                          type = "response")[,1]

# CV20 predictions
prds.train.cv20 <- predict(enet_cv20_model_min, newx = predictors_train,
                           type = "response")[,1]
prds.test.cv20 <- predict(enet_cv20_model_min, newx = predictors_test,
                          type = "response")[,1]


# ROC CURVES 
auc.train.cv10 <- roc(response_train, prds.train.cv10)
auc.test.cv10 <- roc(response_test, prds.test.cv10)
auc.train.cv20 <- roc(response_train, prds.train.cv20)
auc.test.cv20 <- roc(response_test, prds.test.cv20)

par(mfrow = c(2,2))
plot(auc.train.cv10, main = "CV10 - Train")
plot(auc.test.cv10, main = "CV10 - Test")
plot(auc.train.cv20, main = "CV20 - Train")
plot(auc.test.cv20, main = "CV20 - Test")
par(mfrow = c(1,1))

cat("=== AUC Values ===\n")
cat("CV10 - Train AUC:", auc(auc.train.cv10), "\n")
cat("CV10 - Test AUC:", auc(auc.test.cv10), "\n")
cat("CV20 - Train AUC:", auc(auc.train.cv20), "\n")
cat("CV20 - Test AUC:", auc(auc.test.cv20), "\n")

# Calculate optimal classification threshold using min-max approach ====
# Finds threshold that maximises the minimum of sensitivity and specificity
snsp.train.cv10 <- cbind(auc.train.cv10$sensitivities, auc.train.cv10$specificities)
indx.cv10 <- which.max(apply(snsp.train.cv10, 1, min))
cutoff.train.cv10 <- auc.train.cv10$thresholds[indx.cv10]

snsp.test.cv10 <- cbind(auc.test.cv10$sensitivities, auc.test.cv10$specificities)
indx.cv10.test <- which.max(apply(snsp.test.cv10, 1, min))
cutoff.test.cv10 <- auc.test.cv10$thresholds[indx.cv10.test]

snsp.train.cv20 <- cbind(auc.train.cv20$sensitivities, auc.train.cv20$specificities)
indx.cv20 <- which.max(apply(snsp.train.cv20, 1, min))
cutoff.train.cv20 <- auc.train.cv20$thresholds[indx.cv20]

snsp.test.cv20 <- cbind(auc.test.cv20$sensitivities, auc.test.cv20$specificities)
indx.cv20.test <- which.max(apply(snsp.test.cv20, 1, min))
cutoff.test.cv20 <- auc.test.cv20$thresholds[indx.cv20.test]


cat("=== Optimal Thresholds ===\n")
cat("CV10 - Train cutoff:", cutoff.train.cv10, "\n")
cat("CV10 - Test cutoff:", cutoff.test.cv10, "\n")
cat("CV20 - Train cutoff:", cutoff.train.cv20, "\n")
cat("CV20 - Test cutoff:", cutoff.test.cv20, "\n")

# Visualize ROC curves with optimal threshold marked (blue dotted lines) ====
par(mfrow = c(2,2))

plot(auc.train.cv10, main = "CV10 - Train")
abline(h = snsp.train.cv10[indx.cv10, 1], 
       v = snsp.train.cv10[indx.cv10, 2], 
       col = 'blue', lty = 2)

plot(auc.test.cv10, main = "CV10 - Test")
abline(h = snsp.test.cv10[indx.cv10.test, 1], 
       v = snsp.test.cv10[indx.cv10.test, 2], 
       col = 'blue', lty = 2)

plot(auc.train.cv20, main = "CV20 - Train")
abline(h = snsp.train.cv20[indx.cv20, 1], 
       v = snsp.train.cv20[indx.cv20, 2], 
       col = 'blue', lty = 2)

plot(auc.test.cv20, main = "CV20 - Test")
abline(h = snsp.test.cv20[indx.cv20.test, 1], 
       v = snsp.test.cv20[indx.cv20.test, 2], 
       col = 'blue', lty = 2)

par(mfrow = c(1,1))

# Evaluate model performance at optimal threshold ====
# Using training cutoff applied to both train and test sets
conf.mat.train.cv10 <- table(y=response_train, yhat=as.numeric(prds.train.cv10 > cutoff.train.cv10))
conf.mat.test.cv10 <- table(y=response_test, yhat=as.numeric(prds.test.cv10 > cutoff.train.cv10))
conf.mat.train.cv20 <- table(y=response_train, yhat=as.numeric(prds.train.cv20 > cutoff.train.cv20))
conf.mat.test.cv20 <- table(y=response_test, yhat=as.numeric(prds.test.cv20 > cutoff.train.cv20))

# Sensitivity, specificity, accuracy function
sn.sp <- function(mat){
  sn <- mat[2,2]/sum(mat[2,])
  sp <- mat[1,1]/sum(mat[1,])
  acc <- (mat[1,1] + mat[2,2])/sum(mat)
  return(unlist(list(sensitivity=sn, specificity=sp, accuracy=acc)))
}

# Print results
cat("=== CV10 Train ===\n"); sn.sp(conf.mat.train.cv10)
cat("=== CV10 Test ===\n"); sn.sp(conf.mat.test.cv10)
cat("=== CV20 Train ===\n"); sn.sp(conf.mat.train.cv20)
cat("=== CV20 Test ===\n"); sn.sp(conf.mat.test.cv20)

# Variable selection - identify predictors selected by each model ====
coef.cv10 <- coef(enet_cv10_model_min)[-1, 1]
coef.cv10.nonzero <- coef.cv10[coef.cv10 != 0]

# Rank by absolute value
coef.cv10.ranked <- sort(abs(coef.cv10.nonzero), decreasing = TRUE)

cat("=== CV10 Selected Predictors ===\n")
print(coef.cv10.ranked)
cat("Number of predictors selected:", length(coef.cv10.ranked), "\n")

# CV20 variable selection ====
coef.cv20 <- coef(enet_cv20_model_min)[-1, 1]
coef.cv20.nonzero <- coef.cv20[coef.cv20 != 0]

# Rank by absolute value
coef.cv20.ranked <- sort(abs(coef.cv20.nonzero), decreasing = TRUE)

cat("=== CV20 Selected Predictors ===\n")
print(coef.cv20.ranked)
cat("Number of predictors selected:", length(coef.cv20.ranked), "\n")

# Compare predictor selection between CV10 and CV20 models ====
shared <- intersect(names(coef.cv10.nonzero), names(coef.cv20.nonzero))
cat("=== Cytokines selected by both models ===\n")
print(shared)

# Examine direction and magnitude of shared predictors ====
# Positive coefficient = associated with severe, negative = associated with mild
cat("=== CV10 shared predictor coefficients ===\n")
print(coef.cv10.nonzero[shared])
cat("=== CV20 shared predictor coefficients ===\n")
print(coef.cv20.nonzero[shared])

# Visualize relationship between age and COVID-19 severity ====
boxplot(as.numeric(as.character(covid_dataset$AGE)) ~ covid_dataset$Severity,
        main = "Age vs COVID-19 Severity",
        xlab = "Severity",
        ylab = "Age",
        col = c("lightblue", "lightyellow"))




