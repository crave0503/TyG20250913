# --------------------------
# 0. Environment Preparation and Data Loading
# --------------------------
library(tidyverse)
library(caret)
library(pROC)
library(SuperLearner)
library(doParallel)
library(caretEnsemble)
#install.packages("caretEnsemble") # Uncomment to install if the package is not available
# --------------------------
# Corrected Machine Learning Modeling Code
# --------------------------
library(tidyverse)
library(caret)
library(pROC)
library(doParallel)
library(survey)
set.seed(42)  # Fix global random seed for reproducibility
cl <- makeCluster(detectCores()-1)  # Enable parallel computing for acceleration
registerDoParallel(cl)
# Load data
#cleaned_data<-read.csv("cleaned_CAD_data.csv",) # Uncomment to load data
# Alternative way to read CSV file
#cleaned_data1 <- read.csv("cleaned_CAD_data.csv", stringsAsFactors = FALSE)
# Optional: Verify data loading success by checking basic information
head(cleaned_data)  # View the first 6 rows of data
dim(cleaned_data)   # Check data dimensions (number of rows and columns)
str(cleaned_data)   # Check data structure (variable types and values)
# --------------------------
# 1. Data Preprocessing and Variable Selection (Univariate Only)
# --------------------------
# 1.1 Univariate analysis function (weighted logistic regression)
univariate_analysis <- function(data, outcome, covariates, design) {
  significant_vars <- c()
  for (var in covariates) {
    formula <- as.formula(paste(outcome, "~", var))
    model <- tryCatch(
      svyglm(formula, design = design, family = quasibinomial()),
      error = function(e) NULL  # Return NULL if model training fails
    )
    if (!is.null(model)) {
      p_value <- summary(model)$coefficients[2, 4]
      if (p_value < 0.05) significant_vars <- c(significant_vars, var)  # Select variables with p<0.05
    }
  }
  return(significant_vars)
}

# 1.2 Execute univariate variable selection
outcome <- "CVD_CAD"  # Target variable (Coronary Artery Disease)
excluded_vars <- c("SEQN", "TyG_quartile")  # Variables to exclude (e.g., ID, derived quartiles)
covariates <- setdiff(colnames(cleaned_data), c(outcome, excluded_vars))  # All potential predictors

# Split data into training (80%) and testing (20%) sets
trainIndex <- createDataPartition(cleaned_data[[outcome]], p = 0.8, list = FALSE)
trainData <- cleaned_data[trainIndex, ]
testData <- cleaned_data[-trainIndex, ]

# Create weighted survey design object (for complex sampling data, e.g., NHANES)
train_design <- svydesign(
  id = ~SDMVPSU,        # Primary sampling unit
  strata = ~SDMVSTRA,   # Stratification variable
  weights = ~WTMEC2YR,  # Sampling weight
  data = trainData,
  nest = TRUE           # Nested sampling design
)

# Get significant variables from univariate analysis
selected_vars <- univariate_analysis(trainData, outcome, covariates, train_design)
cat("Univariate significant variables:", selected_vars, "\n")
# 0. Data Preprocessing Check --------------------------------------------------------
# Confirm target variable format (critical fix for modeling)
# Ensure target variable is a factor with correct levels
cleaned_data$CVD_CAD <- factor(cleaned_data$CVD_CAD, 
                               levels = c("No", "Yes"), 
                               labels = c("No", "Yes"))

# Example expected output (class distribution):
#  No  Yes 
# 8000  132 
# Handle class imbalance if present (e.g., using SMOTE)
#  Apply SMOTE for oversampling minority class

# Load required packages
library(tidyverse)
library(caret)
library(pROC)
library(doParallel)
library(DMwR)

# 1. Data Splitting (Enforce 8:2 Ratio)------------------------------------------------
set.seed(1010)  # Fix seed for consistent data split
trainIndex <- createDataPartition(
  cleaned_data$CVD_CAD, 
  p = 0.8, 
  list = FALSE,
  times = 1  # Ensure only one split is generated
)
trainData <- cleaned_data[trainIndex, ]
testData <- cleaned_data[-trainIndex, ]
colnames(trainData)  # Check variable names in training data


# Extract features (X) and labels (Y)
train_X <- trainData %>% select(-CVD_CAD)  # Predictors for training set
train_Y <- trainData$CVD_CAD              # Labels for training set
test_X <- testData %>% select(-CVD_CAD)    # Predictors for testing set
test_Y <- testData$CVD_CAD                # Labels for testing set

# Define columns to retain (key predictors)
columns_to_keep <- c("RIDAGEYR", "RIAGENDR", "LBDLDL", "TyG",  "BPXSY1", "BPXDI1", "SMQ020")
# Select specified columns for training predictors
train_X <- train_X %>% select(all_of(columns_to_keep))
# Select specified columns for testing predictors
test_X <- test_X %>% select(all_of(columns_to_keep))

# Check class distribution in training set (before balancing)
print(table(train_Y))

# Apply SMOTE to balance training set class distribution
set.seed(42)
balanced_data <- SMOTE(CVD_CAD ~ ., data = cbind(train_X, CVD_CAD = train_Y))
train_X <- balanced_data[, -ncol(balanced_data)]  # Update training predictors after SMOTE
train_Y <- balanced_data$CVD_CAD                  # Update training labels after SMOTE
# Check class distribution in training set (after balancing)
print(table(train_Y))
# Convert all factor variables to numeric (required for modeling)
train_X <- train_X %>% mutate_if(is.factor, as.numeric)
test_X <- test_X %>% mutate_if(is.factor, as.numeric)

# Add KNN-specific data preprocessing (after SMOTE)
# KNN requires standardized features (mean=0, SD=1) and removal of near-zero variance variables
knn_preProc <- preProcess(train_X, method = c("center", "scale", "nzv"))
train_X_knn <- predict(knn_preProc, train_X)  # Preprocessed predictors for KNN training
test_X_knn <- predict(knn_preProc, test_X)    # Preprocessed predictors for KNN testing
# Ensure factor variables are correctly converted (KNN requires numeric input)
train_X_knn <- train_X_knn %>% mutate(across(where(is.factor), as.numeric))
test_X_knn <- test_X_knn %>% mutate(across(where(is.factor), as.numeric))


# 2. Cross-Validation Configuration (Critical Fix)-----------------------------------------------
ctrl <- trainControl(
  method = "cv",               # 10-fold cross-validation
  number = 10,
  classProbs = TRUE,           # Enable probability prediction (required for ROC/AUC calculation)
  summaryFunction = twoClassSummary,  # Use two-class evaluation metrics (sensitivity, specificity, ROC)
  savePredictions = "final",   # Save final cross-validation predictions
  allowParallel = TRUE,        # Enable parallel computing for cross-validation
  verboseIter = TRUE           # Print training progress (iteration details)
)

# 3. Model Parameter Correction ---------------------------------------------------------
# Define list of base models with optimized parameters
models <- list(
  "RF" = list(
    method = "rf",             # Random Forest
    tuneGrid = data.frame(mtry = c(2, sqrt(ncol(train_X)))), # Tune mtry (2 vs sqrt(number of predictors))
    preProcess = c("center", "scale")  # Standardize features before modeling
  ),
  "XGB" = list(
    method = "xgbTree",        # XGBoost
    tuneGrid = expand.grid(
      nrounds = 100,           # Number of boosting rounds
      max_depth = c(4, 6),     # Tune tree depth (4 vs 6)
      eta = c(0.1, 0.3),       # Tune learning rate (0.1 vs 0.3)
      gamma = 0,               # Minimum loss reduction required for splitting
      colsample_bytree = 0.8,  # Fraction of features used per tree
      min_child_weight = 1,    # Minimum sum of instance weights (hessian) per leaf
      subsample = 0.8          # Fraction of samples used per tree
    ),
    preProcess = c("center", "scale")  # Standardize features before modeling
  ),
  "SVM" = list(
    method = "svmRadial",      # SVM with radial kernel
    tuneLength = 5,            # Tune 5 hyperparameter combinations
    preProcess = c("center", "scale")  # SVM requires feature standardization
  ),
  "GLM" = list(
    method = "glm",            # Generalized Linear Model (logistic regression)
    family = binomial(),       # Binomial family for binary classification
    preProcess = c("center", "scale")  # Standardize features before modeling
  ),
  "KNN" = list(
    method = "knn",            # K-Nearest Neighbors
    tuneGrid = expand.grid(
      k = seq(5, 25, by = 5)   # Tune k value (5, 10, 15, 20, 25)
    ),
    preProcess = NULL,  # Already preprocessed separately (see knn_preProc above)
    metric = "ROC"      # Optimize model based on ROC/AUC
  )
)
# 4. Parallel Model Training (Add Error Handling)-------------------------------------------
# Initialize empty list to store trained models
#model_list <- list()
#for (model_name in names(models)) {
#  cat("\n--- Training:", model_name, "---\n")
#  tryCatch({
#    model_list[[model_name]] <- train(
#      x = train_X,
#      y = train_Y,
#      method = models[[model_name]]$method,
#      trControl = ctrl,
#      tuneGrid = models[[model_name]]$tuneGrid,
#      preProcess = models[[model_name]]$preProcess,
#      metric = "ROC"
#    )
#    print(model_list[[model_name]])  # Print model summary
#  }, error = function(e) {
#    message("Model ", model_name, " training failed: ", e$message)
#  })
#}
# Revised parallel training (use KNN-specific preprocessed data for KNN)
model_list <- list()
for (model_name in names(models)) {
  cat("\n--- Training:", model_name, "---\n")
  
  # Use KNN-specific preprocessed data for KNN model; use regular data for other models
  if(model_name == "KNN") {
    train_data <- train_X_knn
    test_data <- test_X_knn
  } else {
    train_data <- train_X
    test_data <- test_X
  }
  
  tryCatch({
    model_list[[model_name]] <- train(
      x = train_data,
      y = train_Y,
      method = models[[model_name]]$method,
      trControl = ctrl,
      tuneGrid = models[[model_name]]$tuneGrid,
      preProcess = models[[model_name]]$preProcess,
      metric = "ROC"  # Optimize model performance using ROC/AUC
    )
    print(model_list[[model_name]])  # Print detailed model summary
  }, error = function(e) {
    message("Model ", model_name, " training failed: ", e$message)  # Capture and print errors
  })
}
# 5. Result Extraction (Add Null Check)-----------------------------------------------
# Extract valid models (exclude models that failed to train)
valid_models <- model_list[!sapply(model_list, is.null)]

# Generate cross-validation performance comparison table
if (length(valid_models) > 0) {
  results <- resamples(valid_models)
  print(summary(results, metric = "ROC"))  # Print cross-validation ROC metrics (mean ± SD)
   
  # Evaluate model performance on test set
  test_perf <- data.frame()  # Initialize data frame to store test set results
  for (model_name in names(valid_models)) {
    # Predict probabilities for positive class ("Yes") on test set
    pred <- predict(valid_models[[model_name]], test_X, type = "prob")
    # Calculate ROC curve and AUC
    roc_obj <- roc(test_Y, pred[, "Yes"])
    # Append results to data frame
    test_perf <- rbind(test_perf, data.frame(
      Model = model_name,
      AUC = round(roc_obj$auc, 3),  # Round AUC to 3 decimal places
      CI = paste(round(ci(roc_obj), 3), collapse = "-")  # 95% CI for AUC (rounded to 3 decimals)
    ))
  }
  print(test_perf)  # Print test set performance metrics
  
  ##### Ensemble Model Construction (Revised Version: Sep 5)
  # 1. Simple Average Ensemble
  # 1.1 Collect predicted probabilities from all valid base models (handle KNN-specific data)
  ensemble_preds <- lapply(names(valid_models), function(model_name) {
    if(model_name == "KNN") {
      pred_data <- test_X_knn  # Use KNN-specific preprocessed test data
    } else {
      pred_data <- test_X      # Use regular test data for other models
    }
    # Predict positive class probabilities
    predict(valid_models[[model_name]], pred_data, type = "prob")[, "Yes"]
  })
  names(ensemble_preds) <- names(valid_models)
  
  # 1.2 Convert list of predictions to data frame
  pred_df <- as.data.frame(ensemble_preds)
  
  # 1.3 Calculate simple average of predicted probabilities
  pred_df$Average <- rowMeans(pred_df)
  
  # 1.4 Evaluate simple average ensemble performance
  avg_roc <- roc(test_Y, pred_df$Average)
  avg_auc <- round(auc(avg_roc), 3)
  cat("Simple Average Ensemble Test Set AUC:", avg_auc, "\n")
  
  # 2. AUC-Weighted Average Ensemble
  # 2.1 Calculate weights based on base model test set AUC (higher AUC = higher weight)
  model_aucs <- sapply(pred_df[, names(valid_models)], 
                       function(x) auc(roc(test_Y, x)))
  weights <- model_aucs / sum(model_aucs)  # Normalize weights to sum to 1
  
  # 2.2 Calculate weighted average of predicted probabilities
  pred_df$Weighted <- rowSums(t(t(pred_df[, names(valid_models)]) * weights))
  
  # 2.3 Evaluate weighted average ensemble performance
  weighted_roc <- roc(test_Y, pred_df$Weighted)
  weighted_auc <- round(auc(weighted_roc), 3)
  cat("Weighted Average Ensemble Test Set AUC:", weighted_auc, "\n")
  
  # 2.4 Create weight summary table
  weight_summary <- data.frame(
    Model = names(weights),
    Weight = round(weights, 3),  # Round weights to 3 decimal places
    AUC = round(model_aucs, 3),  # Base model test set AUC
    stringsAsFactors = FALSE     # Avoid factor conversion for compatibility
  )
  # Print formatted weight summary
  cat("\n=== Model Weight Allocation Table ===\n")
  print(knitr::kable(weight_summary, align = c("l", "r", "r")))

  # 3. Stacking Ensemble (Meta-Model)
  # 3.1 Enhanced stacking data preparation function (fix KNN compatibility)
  prepare_stacking_data <- function(models, data) {
    preds <- sapply(names(models), function(name) {
      # Apply KNN-specific preprocessing if model is KNN
      if(models[[name]]$method == "knn") {
        data_processed <- predict(knn_preProc, data)
        pred <- predict(models[[name]], data_processed, type = "prob")
      } else {
        pred <- predict(models[[name]], data, type = "prob")
      }
      
      # Unified probability extraction (handles both "X1" and "Yes" column names)
      positive_col <- grep("X1|Yes", colnames(pred), value = TRUE)[1]
      if(is.na(positive_col)) positive_col <- 2  # Default to 2nd column if no match
      return(pred[, positive_col])
    })
    
    # Ensure column names match model names
    colnames(preds) <- names(models)
    return(preds)
  }
  
  # 3.2 Build stacking dataset (add data validation)
  train_stacking <- prepare_stacking_data(valid_models, train_X)
  
  # Check consistency between model names and prediction column names
  if(!identical(colnames(train_stacking), names(valid_models))) {
    warning("Mismatch between base model names and prediction feature names!")
  }
  
  # Create stacking dataset (predictors = base model predictions; target = original labels)
  stack_data <- data.frame(train_stacking, Outcome = factor(train_Y, levels = c("X0", "X1")))
  
  # 3.3 Configure enhanced glmnet meta-model (elastic net)
  stack_ctrl <- trainControl(
    method = "cv",               # 5-fold cross-validation for meta-model
    number = 5,
    classProbs = TRUE,           # Enable probability prediction
    summaryFunction = twoClassSummary,  # Use two-class metrics
    allowParallel = FALSE,       # Disable parallel to avoid conflicts
    verboseIter = TRUE           # Print meta-model training progress
  )
  
  # 3.4 Train elastic net meta-model (with cross-validation for lambda)
  lambda_seq <- 10^seq(-3, 3, length = 50)  # Fine-grained lambda sequence (regularization strength)
  
  # # Original meta-model training code (commented out for revision)
  # stack_model <- train(
  #   Outcome ~ .,
  #   data = stack_data,
  #   method = "glmnet",
  #   family = "binomial",
  #   tuneGrid = expand.grid(
  #     alpha = seq(0, 1, 0.1),  # Comprehensive alpha range (0=ridge, 1=lasso)
  #     lambda = lambda_seq
  #   ),
  #   metric = "ROC",
  #   trControl = stack_ctrl
  # )
  
  # 3.5 Revised stacking data preparation (use corrected labels)
  stack_data <- data.frame(train_stacking, Outcome = factor(train_Y))  # Use original factor labels directly
  
  # 3.6 Check for NA values or outliers in outcome column
  colSums(is.na(stack_data))    # Count NA values per column
  table(stack_data$Outcome)     # Verify outcome distribution (should show "No" and "Yes")
  
  # 3.7 Remove rows with NA values (if any)
  stack_data <- na.omit(stack_data)
  
  # 3.8 Retrain stacking meta-model with cleaned data
  stack_model <- train(
    Outcome ~ .,
    data = stack_data,
    method = "glmnet",
    family = "binomial",
    tuneGrid = expand.grid(alpha = seq(0, 1, 0.1), lambda = lambda_seq),
    metric = "ROC",
    trControl = stack_ctrl
  )  
  
  # 3.9 Evaluate stacking ensemble performance
  test_stacking <- prepare_stacking_data(valid_models, test_X)
  test_stacking_df <- data.frame(test_stacking)
  
  # Safe prediction function (handles column name inconsistencies)
  safe_predict <- function(model, newdata) {
    pred <- predict(model, newdata = newdata, type = "prob")
    # Extract positive class probability (handles "X1" or "Yes" column names)
    positive_col <- grep("X1|Yes", colnames(pred), value = TRUE)[1]
    return(pred[, positive_col])
  }
  
  # Predict probabilities using stacking meta-model
  stack_pred <- safe_predict(stack_model, test_stacking_df)
  
  # Calculate performance metrics for stacking ensemble
  stack_roc <- roc(test_Y, stack_pred)
  stack_auc <- round(auc(stack_roc), 3)
  ci_auc <- round(ci.auc(stack_roc), 3)  # 95% CI for AUC
  
  cat("Optimized Stacking Ensemble Performance:\n")
  cat("AUC:", stack_auc, "95% CI:", paste(ci_auc, collapse = "-"), "\n")
  
  # # Visualize ROC curve for stacking model (commented out if not needed)
  # plot(stack_roc, main = paste("Stacking Model ROC (AUC =", stack_auc, ")"), 
  #     print.auc = TRUE, print.thres = "best")
  
  # --------------------------
  # 7. Final Performance Comparison
  # --------------------------
  # Combine base model and ensemble model performance results
  final_perf <- rbind(
    test_perf,
    data.frame(
      Model = c("Average_Ensemble", "Weighted_Ensemble","Stack_Ensemble"),
      AUC = c(avg_auc, weighted_auc, stack_auc),
      CI = NA  # Add CI calculation here if needed
    )
  )
  
  # Sort results by AUC (descending order)
  final_perf <- final_perf %>% arrange(desc(AUC))
  print(final_perf)  # Print final performance ranking
 
  
  ###### Save Models as Independent RDS Files 
  # Create directory for model storage (skip if directory already exists)
  dir.create("model_onlyTYG_rds", showWarnings = FALSE)
  # Save valid base models (only models that trained successfully)
  for (model_name in names(valid_models)) {
    saveRDS(valid_models[[model_name]], 
            file = file.path("model_onlyTYG_rds", paste0(model_name, ".rds")))
  }
  # Save stacking ensemble model (only if it exists)
  if (!is.null(stack_model)) {
    saveRDS(stack_model, file = file.path("model_onlyTYG_rds", "stacked_ensemble.rds"))
  }
  # Save ensemble predictions (numeric vectors, not model objects)
  saveRDS(pred_df$Average, file = file.path("model_onlyTYG_rds", "simple_average.rds"))
  saveRDS(pred_df$Weighted, file = file.path("model_onlyTYG_rds", "weighted_average.rds"))  # Fix: Use pred_df$Weighted instead of undefined "weighted_avg"
  # Notes:
  # 1. Base models (max 4) are saved as individual RDS files
  # 2. Stacking ensemble model is saved as "stacked_ensemble.rds"
  # 3. Simple average and weighted average predictions are saved as numeric vectors
  # 4. All files are stored in the "model_onlyTYG_rds" directory (current working directory)
  
  # Calculate and save training set standardization parameters (for new data preprocessing)
  train_means <- colMeans(train_X)  # Mean of each feature in training set
  train_sds <- apply(train_X, 2, sd)# Standard deviation of each feature in training set
  saveRDS(list(train_means = train_means, train_sds = train_sds), 
          file = "model_onlyTYG_rds/standardization_params.rds")
  # Re-save to ensure parameters are stored in the correct directory
  saveRDS(list(train_means = train_means, train_sds = train_sds), 
          file = "model_onlyTYG_rds/standardization_params.rds")
  
 #### Visualization Code (Revised Version: Sep 5)
  library(pROC)      # For ROC curve analysis
  library(PRROC)     # For Precision-Recall (PR) curve analysis
  library(rmda)      # For Decision Curve Analysis (DCA)
  library(ggplot2)   # For advanced plotting
  library(patchwork) # For combining multiple plots
  
  # Initialize list to store plot data (ROC, PR, DCA curves)
  plot_data <- list(
    roc_curves = list(),
    pr_curves = list(),
    dca_curves = list()
  )
  
  # Define color scheme for models (consistent across all plots)
  model_colors <- c(
    RF = "#BFDFD2", 
    XGB = "#51999F",
    SVM = "#4198AC",
    GLM = "#7BC0CD",
    KNN = "#DBCB92",
    Average_Ensemble = "#ECB66C",
    Weighted_Ensemble = "#EA9E58",
    Stack_Ensemble = "#ED8D5A"
  ) 
  
  # Calculate ROC curves for base models
  for(model_name in names(valid_models)) {
    # Predict positive class probabilities on test set
    pred <- predict(valid_models[[model_name]], test_X, type = "prob")[, "Yes"]
    # Calculate ROC curve
    plot_data$roc_curves[[model_name]] <- roc(test_Y, pred)
  }
  
  # Calculate ROC curves for ensemble models
  plot_data$roc_curves$Average_Ensemble <- roc(test_Y, pred_df$Average)
  plot_data$roc_curves$Weighted_Ensemble <- roc(test_Y, pred_df$Weighted)
  plot_data$roc_curves$Stack_Ensemble <- stack_roc
  
  ## Plot ROC Curves
  roc_plot <- ggplot() +
    # Add diagonal line (random guess, AUC = 0.5)
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "ROC Curve Comparison",
         x = "1 - Specificity", 
         y = "Sensitivity",
         color = "Model") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
  
  # Add ROC curves for each model
  for(model_name in names(plot_data$roc_curves)) {
    roc_df <- data.frame(
      Sensitivity = plot_data$roc_curves[[model_name]]$sensitivities,
      Specificity = plot_data$roc_curves[[model_name]]$specificities,
      Model = model_name
    )
    
    roc_plot <- roc_plot +
      geom_line(data = roc_df, 
                aes(x = 1 - Specificity, y = Sensitivity, color = Model),
                linewidth = 0.8)
  }
  
  # Add AUC values to legend labels
  auc_labels <- sapply(names(plot_data$roc_curves), function(x) {
    paste0(x, " (AUC = ", round(auc(plot_data$roc_curves[[x]]), 3), ")")
  })
  
  roc_plot <- roc_plot +
    scale_color_manual(values = model_colors,
                       labels = auc_labels)
  
  print(roc_plot)
  
  ## Fix: If plot is not displayed, reset graphics device
  # Clear existing graphics devices
  try(dev.off(), silent = TRUE)
  
  # Create new graphics device (8x6 inches, 300 DPI)
  dev.new(width = 8, height = 6, unit = "in", res = 300)
  
  # Reconstruct ROC plot with data validation
  roc_plot <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "ROC Curve Comparison",
         x = "1 - Specificity", 
         y = "Sensitivity",
         color = "Model") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))  # Center plot title
  
  # Add valid ROC curves (skip models with missing data)
  for(model_name in names(plot_data$roc_curves)) {
    if(!is.null(plot_data$roc_curves[[model_name]])) {
      roc_df <- data.frame(
        Sensitivity = plot_data$roc_curves[[model_name]]$sensitivities,
        Specificity = plot_data$roc_curves[[model_name]]$specificities,
        Model = model_name
      )
      
      roc_plot <- roc_plot +
        geom_line(data = roc_df, 
                  aes(x = 1 - Specificity, y = Sensitivity, color = Model),
                  linewidth = 0.8)
    } else {
      warning(paste("Missing ROC data for model:", model_name))
    }
  }
  
  # Add legend labels with AUC values (exclude NA values)
  auc_values <- sapply(plot_data$roc_curves, function(x) {
    if(!is.null(x)) round(auc(x), 3) else NA
  })
  legend_labels <- paste0(names(auc_values), " (AUC = ", auc_values, ")")
  
  # Finalize ROC plot
  roc_plot <- roc_plot +
    scale_color_manual(values = model_colors,
                       labels = legend_labels[!is.na(auc_values)])
  
  print(roc_plot) 
  
# Plot Precision-Recall (PR) Curves  
  # Convert factor labels to numeric (PRROC package requires 0 = negative, 1 = positive)
  y_numeric <- as.numeric(test_Y) - 1
  
  # Calculate PR curves for base models
  for(model_name in names(valid_models)) {
    # Predict positive class probabilities on test set
    pred <- predict(valid_models[[model_name]], test_X, type = "prob")[, "Yes"]
    # Calculate PR curve
    plot_data$pr_curves[[model_name]] <- pr.curve(
      scores.class0 = pred,
      weights.class0 = y_numeric,
      curve = TRUE  # Return full curve data
    )
  }
  
  # Calculate PR curves for ensemble models
  plot_data$pr_curves$Average_Ensemble <- pr.curve(
    scores.class0 = pred_df$Average,
    weights.class0 = y_numeric,
    curve = TRUE
  )
  plot_data$pr_curves$Weighted_Ensemble <- pr.curve(
    scores.class0 = pred_df$Weighted,
    weights.class0 = y_numeric,
    curve = TRUE
  )
  plot_data$pr_curves$Stack_Ensemble <- pr.curve(
    scores.class0 = stack_pred,
    weights.class0 = y_numeric,
    curve = TRUE
  )  
  
  # Create PR curve plot
  pr_plot <- ggplot() +
    labs(title = "Precision-Recall Curve Comparison",
         x = "Recall", 
         y = "Precision",
         color = "Model") +
    theme_minimal(base_size = 12) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +  # Fix axis range to [0,1]
    theme(legend.position = "right")
  
  # Add PR curves for each model
  for(model_name in names(plot_data$pr_curves)) {
    pr_df <- data.frame(
      Recall = plot_data$pr_curves[[model_name]]$curve[, 1],
      Precision = plot_data$pr_curves[[model_name]]$curve[, 2],
      Model = model_name
    )
    
    pr_plot <- pr_plot +
      geom_line(data = pr_df, 
                aes(x = Recall, y = Precision, color = Model),
                linewidth = 0.8)
  }
  
  # Add AUC-PR values to legend labels
  prauc_labels <- sapply(names(plot_data$pr_curves), function(x) {
    paste0(x, " (AUC-PR = ", round(plot_data$pr_curves[[x]]$auc.integral, 3), ")")
  })
  
  pr_plot <- pr_plot +
    scale_color_manual(values = model_colors,
                       labels = prauc_labels)
  
  print(pr_plot)
  
# Plot Decision Curve Analysis (DCA) Curves  
  # Ensure rmda package is loaded (for DCA calculations)
  library(rmda)
  
  # Safe method: Create DCA data frame using original test set labels
  dca_data <- data.frame(
    Outcome = ifelse(testData$CVD_CAD == "Yes", 1, 0),  # Convert labels to 0/1 (1 = positive)
    # Add predicted probabilities for base models
    GLM = predict(valid_models$GLM, test_X, type = "prob")[, "Yes"],
    RF = predict(valid_models$RF, test_X, type = "prob")[, "Yes"],
    XGB = predict(valid_models$XGB, test_X, type = "prob")[, "Yes"],
    SVM = predict(valid_models$SVM, test_X, type = "prob")[, "Yes"],
    KNN = predict(valid_models$KNN, test_X_knn, type = "prob")[, "Yes"],
    # Add predicted probabilities for ensemble models
    Average = pred_df$Average,
    Weighted = pred_df$Weighted,
    Stacking = stack_pred
  )
  
  # Re-verify outcome distribution in DCA data
  cat("Outcome distribution in DCA data:\n")
  table(dca_data$Outcome)
  cat("Unique values in Outcome column:", unique(dca_data$Outcome), "\n")  # Should be 0 and 1
  
  # Check data summary (identify potential issues)
  summary(dca_data)
  
  # Check for NA values and remove rows with NA if present
  colSums(is.na(dca_data))
  dca_data <- na.omit(dca_data) 
  
  # Install and load dcurves package (for enhanced DCA plotting)
  if (!require(dcurves)) install.packages("dcurves")
  library(dcurves)
  
  # Run DCA (exclude potentially problematic models, e.g., KNN)
  dca_result_fixed <- dca(
    Outcome ~ GLM + RF + XGB + SVM + Average + Weighted + Stacking,
    data = dca_data,
    thresholds = seq(0, 0.4, by = 0.01)  # Threshold range relevant for clinical decision-making
  )
  
  # Plot DCA with custom colors
  final_plot <- plot(dca_result_fixed, smooth = TRUE) +
    scale_color_manual(
      values = c(
        'GLM' = '#7BC0CD',
        'RF' = '#BFDFD2', 
        'XGB' = '#51999F',
        'SVM' = '#4198AC',
        'Average' = '#ECB66C',
        'Weighted' = '#EA9E58',
        'Stacking' = '#ED8D5A'
      )
    ) +
    theme_minimal() +
    theme(legend.position = "topright",
          legend.title = element_blank()) +
    labs(title = "Decision Curve Analysis (Excluding Poor Performing Models)",
         x = "Threshold Probability", 
         y = "Net Benefit")
  
  # Display DCA plot
  print(final_plot)
  
  # Save high-resolution DCA plot
  ggsave("DCA_Curve_Fixed.png", final_plot, 
         width = 10, height = 6, dpi = 300)
  
# DCA Plot (Revised Version: Sep 5, 16:30)
  # Run full DCA including all models
  dca_result <- dca(
    Outcome ~ GLM + RF + XGB + SVM + KNN + Average + Weighted + Stacking,
    data = dca_data,
    thresholds = seq(0, 0.4, by = 0.01)
  )
  
  # Plot DCA with custom styling
  dca_plot <- plot(dca_result, smooth = TRUE) +
    # Manually specify color scheme (matches ROC/PR plots)
    ggplot2::scale_color_manual(
      name = "Prediction Models",  # Legend title
      values = c(
        'GLM' = '#7BC0CD',
        'RF' = '#BFDFD2',
        'XGB' = '#51999F', 
        'SVM' = '#4198AC',
        'KNN' = '#DBCB92',
        'Average' = '#ECB66C',
        'Weighted' = '#EA9E58',
        'Stacking' = '#ED8D5A',
        'All' = 'black',      # Color for "Treat all" strategy
        'None' = 'gray50'     # Color for "Treat none" strategy
      ),
      # Define legend display order (prioritize better-performing models)
      breaks = c('Stacking', 'Weighted', 'Average', 'GLM', 'RF', 'XGB', 'SVM', 'KNN', 'All', 'None')
    ) +
    # Improve plot theme and labels
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.8, "lines")
    ) +
    ggplot2::labs(
      title = "Decision Curve Analysis for CAD Prediction",
      x = "Threshold Probability", 
      y = "Net Benefit",
      color = "Models & Strategies"
    ) +
    # Adjust y-axis range to better display net benefit values
    ggplot2::coord_cartesian(ylim = c(-0.1, 1.0))
  
  # Display DCA plot
  print(dca_plot)
  
  # Uncomment to save high-resolution DCA plot
  #ggplot2::ggsave("DCA_Full_With_Legend.png", dca_plot, 
  #                width = 12, height = 8, dpi = 300, bg = "white")
  
## Combine ROC, PR, and DCA Plots into One Figure  
  combined_plot <- (roc_plot + theme(legend.position = "none")) + 
    (pr_plot + theme(legend.position = "none")) + 
    (dca_plot + theme(legend.position = "right")) +
    plot_layout(ncol = 2, guides = "collect") +  # Arrange in 2 columns, share legend
    plot_annotation(title = "Comprehensive Model Performance Evaluation",
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))  # Center main title
  
  print(combined_plot) 
  
  
  
## Plot Calibration Curves (Optional: Comment out if not needed)
  # Load required packages for calibration curves
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Fix test set label encoding to match training set
  test_Y <- factor(test_Y, levels = c("X0", "X1"))
  levels(test_Y) <- c("No", "Yes")  # Map "X0" to "No" and "X1" to "Yes"
  
  # Alternative: Re-extract labels directly from original test data
  test_Y <- testData$CVD_CAD
  
  # Verify corrected label distribution
  cat("Corrected test set label distribution:\n")
  print(table(test_Y))
  cat("Numeric conversion result:\n")
  test_Y_numeric <- ifelse(test_Y == "Yes", 1, 0)
  print(table(test_Y_numeric))
  
  # 1. Build calibration data frame ------------------------------------------------------------
  # Create data frame with predicted probabilities and actual labels
  calibration_data <- data.frame(
    Outcome = ifelse(test_Y == "Yes", 1, 0),  # Convert target to numeric (0 = No, 1 = Yes)
    
    # Predicted probabilities for base models
    GLM = predict(valid_models$GLM, test_X, type = "prob")[, "Yes"],
    RF = predict(valid_models$RF, test_X, type = "prob")[, "Yes"],
    XGB = predict(valid_models$XGB, test_X, type = "prob")[, "Yes"],
    SVM = predict(valid_models$SVM, test_X, type = "prob")[, "Yes"],
    KNN = predict(valid_models$KNN, test_X_knn, type = "prob")[, "Yes"],
    
    # Predicted probabilities for ensemble models
    Average = pred_df$Average,
    Weighted = pred_df$Weighted,
    Stacking = stack_pred
  )
  
  # Remove rows with NA values
  calibration_data <- na.omit(calibration_data)
  
  # Check calibration data frame structure
  cat("Calibration data frame dimensions:", dim(calibration_data), "\n")
  cat("Column summaries:\n")
  print(summary(calibration_data))
  
  # 2. Convert data to long format (for ggplot2 compatibility) --------------------------------------------------------
  cal_data_long <- calibration_data %>%
    pivot_longer(
      cols = -Outcome,  # Exclude outcome column
      names_to = "Model",
      values_to = "Predicted_Probability"
    )
  
  # 3. Calculate smoothed calibration curves --------------------------------------------------------
  # Use LOESS (Local Polynomial Regression) to smooth observed proportions
  cal_curve_data <- cal_data_long %>%
    group_by(Model) %>%
    mutate(Smoothed_Observed = predict(
      loess(Outcome ~ Predicted_Probability, span = 0.75)  # Span controls smoothness (0.75 = balanced)
    )) %>%
    ungroup()
  
  # 4. Plot calibration curves ----------------------------------------------------------
  cal_plot <- ggplot(cal_curve_data, aes(x = Predicted_Probability, y = Smoothed_Observed, color = Model)) +
    geom_line(linewidth = 1.2) +
    # Add ideal calibration line (predicted probability = observed proportion)
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    
    # Use custom color scheme (matches other plots)
    scale_color_manual(
      values = c(
        'GLM' = '#7BC0CD',
        'RF' = '#BFDFD2',
        'XGB' = '#51999F',
        'SVM' = '#4198AC',
        'KNN' = '#DBCB92',
        'Average' = '#ECB66C',
        'Weighted' = '#EA9E58',
        'Stacking' = '#ED8D5A'
      )
    ) +
    
    # Set equal axis scales and range [0,1]
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "Model Calibration Curves",
      x = "Predicted Probability",
      y = "Smoothed Observed Proportion",
      color = "Model",
      caption = "Dashed line = ideal calibration (predicted probability = observed proportion)"
    ) +
    
    # Improve plot theme
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  # 5. Display and save calibration plot --------------------------------------------------------
  print(cal_plot)
  
  # Save high-resolution calibration plot
  ggsave("Calibration_Curves.png", cal_plot, 
         width = 10, height = 6, dpi = 300, bg = "white")
  
  # 6. (Optional) Add binned observed proportions to calibration plot --------------------------------------
  # Calculate binned calibration data (10 bins, 0-0.1 to 0.9-1.0)
  binned_data <- cal_data_long %>%
    group_by(Model) %>%
    mutate(Bin = cut(Predicted_Probability, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
    group_by(Model, Bin) %>%
    summarize(
      Mean_Predicted = mean(Predicted_Probability),  # Mean predicted probability per bin
      Mean_Observed = mean(Outcome),                # Mean observed proportion per bin
      N = n(),                                      # Sample size per bin
      .groups = 'drop'
    ) %>%
    filter(N > 5)  # Keep only bins with ≥5 samples (reliable estimates)
  
  # Add binned points to calibration plot
  cal_plot_with_points <- cal_plot +
    geom_point(data = binned_data, 
               aes(x = Mean_Predicted, y = Mean_Observed, size = N),
               alpha = 0.6, shape = 16) +  # Semi-transparent circular points
    scale_size_continuous(range = c(2, 6), name = "Sample Size") +  # Adjust point size range
    labs(title = "Model Calibration Curves (With Binned Observations)")
  
  print(cal_plot_with_points)