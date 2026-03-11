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
cat("\nBest alpha (lambda.min):", best_alpha_min10, "\n")
cat("Best alpha (lambda.1se):", best_alpha_1se10, "\n")



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

# how many predictor left? ====
# Extract lambda values
enet_lambda_min <- cv_enet$lambda.min
enet_lambda_1se <- cv_enet$lambda.1se

cat("Elastic Net - Lambda min:", enet_lambda_min, "\n")
cat("Elastic Net - Lambda 1se:", enet_lambda_1se, "\n")

# Fit final models
enet_model_min <- glmnet(X_train_scaled, y_train, 
                         alpha = 0.5, 
                         lambda = enet_lambda_min)

enet_model_1se <- glmnet(X_train_scaled, y_train, 
                         alpha = 0.5, 
                         lambda = enet_lambda_1se)

