# --------------------------
# 0. Environment Preparation and Data Loading
# --------------------------
library(tidyverse)
library(caret)
library(pROC)
library(SuperLearner)
library(doParallel)
library(caretEnsemble)
#install.packages("caretEnsemble") # Uncomment to install if not present
# --------------------------
# Corrected Machine Learning Modeling Code
# --------------------------
library(tidyverse)
library(caret)
library(pROC)
library(doParallel)
library(survey)
set.seed(42)  # Fix global random seed for reproducibility
cl <- makeCluster(detectCores()-1)  # Enable parallel computing
registerDoParallel(cl)
# Load data
#cleaned_data<-read.csv("cleaned_CAD_data.csv") # Uncomment to load your data
# Alternative: Load CSV file
#cleaned_data <- read.csv("cleaned_CAD_data.csv", stringsAsFactors = FALSE)
# Optional: Verify data loading success
# head(cleaned_data)  # View first 6 rows
# dim(cleaned_data)   # Check data dimensions (rows x columns)
# str(cleaned_data)   # Check data structure (variable types)
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
      error = function(e) NULL  # Catch errors and return NULL for failed models
    )
    if (!is.null(model)) {
      p_value <- summary(model)$coefficients[2, 4]
      if (p_value < 0.05) significant_vars <- c(significant_vars, var) # Select vars with p<0.05
    }
  }
  return(significant_vars)
}

# 1.2 Execute univariate selection
outcome <- "CVD_CAD"  # Target variable (Coronary Artery Disease)
excluded_vars <- c("SEQN", "TyG_quartile")  # Variables to exclude (e.g., ID, derived quartiles)
covariates <- setdiff(colnames(cleaned_data), c(outcome, excluded_vars)) # All predictors except excluded

# Split data into training (80%) and testing (20%) sets
trainIndex <- createDataPartition(cleaned_data[[outcome]], p = 0.8, list = FALSE)
trainData <- cleaned_data[trainIndex, ]
testData <- cleaned_data[-trainIndex, ]

# Create weighted survey design (for complex sampling data, e.g., NHANES)
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
# Confirm target variable format (critical fix)
cleaned_data$CVD_CAD <- factor(cleaned_data$CVD_CAD, 
                               levels = c("No", "Yes"), 
                               labels = c("No", "Yes"))

# Example output (expected class distribution):
#  No  Yes 
# 8000  132 
# Handle class imbalance with SMOTE if needed
#  Use SMOTE for oversampling minority class

# Load required packages
library(tidyverse)
library(caret)
library(pROC)
library(doParallel)
library(DMwR)

# 1. Data Splitting (Enforce 8:2 Ratio) ------------------------------------------------
set.seed(1010)  # Fix seed for consistent split
trainIndex <- createDataPartition(
  cleaned_data$CVD_CAD, 
  p = 0.8, 
  list = FALSE,
  times = 1  # Ensure single split
)
trainData <- cleaned_data[trainIndex, ]
testData <- cleaned_data[-trainIndex, ]
# colnames(trainData) # Uncomment to check variable names

# Extract features (X) and labels (Y)
train_X <- trainData %>% select(-CVD_CAD)  # Predictors for training
train_Y <- trainData$CVD_CAD              # Labels for training
test_X <- testData %>% select(-CVD_CAD)    # Predictors for testing
test_Y <- testData$CVD_CAD                # Labels for testing

# Define columns to retain (critical predictors)
columns_to_keep <- c("RIDAGEYR", "RIAGENDR", "LBXTC", "LBDLDL", "TyG", "LBXGLU", "BPXSY1", "BPXDI1", "SMQ020")
# Select columns for training predictors
train_X <- train_X %>% select(all_of(columns_to_keep))
# Select columns for testing predictors
test_X <- test_X %>% select(all_of(columns_to_keep))

# Check class distribution in training set (before balancing)
print(table(train_Y))

# Apply SMOTE for oversampling minority class
set.seed(42)
balanced_data <- SMOTE(CVD_CAD ~ ., data = cbind(train_X, CVD_CAD = train_Y))
train_X <- balanced_data[, -ncol(balanced_data)]  # Update training predictors
train_Y <- balanced_data$CVD_CAD                  # Update training labels
# Check class distribution after SMOTE (should be balanced)
print(table(train_Y))
# Convert all factor variables to numeric
train_X <- train_X %>% mutate_if(is.factor, as.numeric)
test_X <- test_X %>% mutate_if(is.factor, as.numeric)

# Add KNN-specific preprocessing (after SMOTE)
# KNN requires standardized features (mean=0, SD=1)
knn_preProc <- preProcess(train_X, method = c("center", "scale", "nzv")) # Center/scale + remove near-zero variance
train_X_knn <- predict(knn_preProc, train_X)  # Preprocessed predictors for KNN training
test_X_knn <- predict(knn_preProc, test_X)    # Preprocessed predictors for KNN testing
# Ensure factor variables are correctly converted (KNN needs numeric input)
train_X_knn <- train_X_knn %>% mutate(across(where(is.factor), as.numeric))
test_X_knn <- test_X_knn %>% mutate(across(where(is.factor), as.numeric))


# 2. Cross-Validation Configuration (Critical Fix) -----------------------------------------------
ctrl <- trainControl(
  method = "cv",               # 10-fold cross-validation
  number = 10,
  classProbs = TRUE,           # Enable probability prediction (required for ROC/AUC)
  summaryFunction = twoClassSummary,  # Use two-class metrics (sens, spec, ROC)
  savePredictions = "final",   # Save final CV predictions
  allowParallel = TRUE,        # Enable parallel computing
  verboseIter = TRUE           # Print training progress
)

# 3. Model Parameter Correction ---------------------------------------------------------
# Define list of base models with optimized parameters
models <- list(
  "RF" = list(
    method = "rf",             # Random Forest
    tuneGrid = data.frame(mtry = c(2, sqrt(ncol(train_X)))), # Tune mtry (2 vs sqrt(p))
    preProcess = c("center", "scale")  # Standardize features
  ),
  "XGB" = list(
    method = "xgbTree",        # XGBoost
    tuneGrid = expand.grid(
      nrounds = 100,           # Number of boosting rounds
      max_depth = c(4, 6),     # Tree depth (tune 4 vs 6)
      eta = c(0.1, 0.3),       # Learning rate (tune 0.1 vs 0.3)
      gamma = 0,               # Minimum loss reduction for splits
      colsample_bytree = 0.8,  # Fraction of features per tree
      min_child_weight = 1,    # Minimum sum of instance weights per leaf
      subsample = 0.8          # Fraction of samples per tree
    ),
    preProcess = c("center", "scale")  # Standardize features
  ),
  "SVM" = list(
    method = "svmRadial",      # Radial kernel SVM
    tuneLength = 5,            # Tune 5 hyperparameter combinations
    preProcess = c("center", "scale")  # SVM requires standardized features
  ),
  "GLM" = list(
    method = "glm",            # Logistic Regression (Generalized Linear Model)
    family = binomial(),       # Binomial family for binary classification
    preProcess = c("center", "scale")  # Standardize features
  ),
  "KNN" = list(
    method = "knn",            # K-Nearest Neighbors
    tuneGrid = expand.grid(
      k = seq(5, 25, by = 5)   # Tune k (5, 10, 15, 20, 25)
    ),
    preProcess = NULL,         # Already preprocessed separately (knn_preProc)
    metric = "ROC"             # Optimize for ROC
  )
)
# 4. Parallel Model Training (Add Error Handling) -------------------------------------------
model_list <- list()  # Initialize list to store trained models
for (model_name in names(models)) {
  cat("\n--- Training:", model_name, "---\n")
  
  # Use KNN-specific preprocessed data for KNN model
  if(model_name == "KNN") {
    train_data <- train_X_knn
    test_data <- test_X_knn
  } else {
    train_data <- train_X
    test_data <- test_X
  }
  
  tryCatch({
    # Train model using caret::train
    model_list[[model_name]] <- train(
      x = train_data,
      y = train_Y,
      method = models[[model_name]]$method,
      trControl = ctrl,
      tuneGrid = models[[model_name]]$tuneGrid,
      preProcess = models[[model_name]]$preProcess,
      metric = "ROC"  # Optimize model for ROC/AUC
    )
    print(model_list[[model_name]])  # Print model summary (best parameters + CV metrics)
  }, error = function(e) {
    # Catch and print errors for failed models
    message("Model ", model_name, " training failed: ", e$message)
  })
}
# 5. Result Extraction (Add Null Check) -----------------------------------------------
# Filter out null models (models that failed to train)
valid_models <- model_list[!sapply(model_list, is.null)]

# Compare cross-validation performance of valid models
if (length(valid_models) > 0) {
  results <- resamples(valid_models)
  print(summary(results, metric = "ROC"))  # Print CV ROC metrics (mean + SD)
   
  # Evaluate performance on test set
  test_perf <- data.frame()  # Initialize data frame to store test metrics
  for (model_name in names(valid_models)) {
    # Predict probabilities for test set (positive class = "Yes")
    pred <- predict(valid_models[[model_name]], test_X, type = "prob")
    # Calculate ROC curve and AUC
    roc_obj <- roc(test_Y, pred[, "Yes"])
    # Append results to data frame
    test_perf <- rbind(test_perf, data.frame(
      Model = model_name,
      AUC = round(roc_obj$auc, 3),  # Round AUC to 3 decimal places
      CI = paste(round(ci(roc_obj), 3), collapse = "-")  # 95% CI for AUC
    ))
  }
  print(test_perf)  # Print test set performance
  
  ##### Ensemble Model Construction (Revised: Sep 5)
  # 1. Simple Average Ensemble
  # Collect predicted probabilities from all valid base models
  ensemble_preds <- lapply(names(valid_models), function(model_name) {
    if(model_name == "KNN") {
      pred_data <- test_X_knn  # Use KNN-specific preprocessed data
    } else {
      pred_data <- test_X
    }
    # Predict probabilities for positive class ("Yes")
    predict(valid_models[[model_name]], pred_data, type = "prob")[, "Yes"]
  })
  names(ensemble_preds) <- names(valid_models)
  
  # Convert list of predictions to data frame
  pred_df <- as.data.frame(ensemble_preds)
  
  # Calculate simple average of probabilities
  pred_df$Average <- rowMeans(pred_df)
  
  # Evaluate simple average ensemble
  avg_roc <- roc(test_Y, pred_df$Average)
  avg_auc <- round(auc(avg_roc), 3)
  cat("Simple Average Ensemble Test Set AUC:", avg_auc, "\n")
  
  # 2. AUC-Weighted Average Ensemble
  # Calculate weights based on base model test set AUC (higher AUC = higher weight)
  model_aucs <- sapply(pred_df[, names(valid_models)], 
                       function(x) auc(roc(test_Y, x)))
  weights <- model_aucs / sum(model_aucs)  # Normalize weights to sum to 1
  
  # Calculate weighted average of probabilities
  pred_df$Weighted <- rowSums(t(t(pred_df[, names(valid_models)]) * weights))
  
  # Evaluate weighted average ensemble
  weighted_roc <- roc(test_Y, pred_df$Weighted)
  weighted_auc <- round(auc(weighted_roc), 3)
  cat("Weighted Average Ensemble Test Set AUC:", weighted_auc, "\n")
  
  # Create weight summary table
  weight_summary <- data.frame(
    Model = names(weights),
    Weight = round(weights, 3),  # Round weights to 3 decimal places
    AUC = round(model_aucs, 3),  # Base model AUCs
    stringsAsFactors = FALSE     # Avoid factor conversion for compatibility
  )
  # Print weight summary
  cat("\n=== Model Weight Allocation ===\n")
  print(knitr::kable(weight_summary, align = c("l", "r", "r")))

  # 3. Stacking Ensemble (Meta-Model)
  # 3.1 Enhanced stacking data preparation (fix KNN compatibility)
  prepare_stacking_data <- function(models, data) {
    preds <- sapply(names(models), function(name) {
      # Use KNN-specific preprocessing for KNN model
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
  
  # Check feature name consistency
  if(!identical(colnames(train_stacking), names(valid_models))) {
    warning("Mismatch between base model names and prediction feature names!")
  }
  
  # Create stacking dataset (predictors = base model predictions; target = original label)
  stack_data <- data.frame(train_stacking, Outcome = factor(train_Y, levels = c("X0", "X1")))
  
  # 3.3 Configure enhanced glmnet meta-model (elastic net)
  stack_ctrl <- trainControl(
    method = "cv",               # 5-fold CV for meta-model
    number = 5,
    classProbs = TRUE,           # Enable probability prediction
    summaryFunction = twoClassSummary,  # Two-class metrics
    allowParallel = FALSE,       # Disable parallel for meta-model (avoid conflicts)
    verboseIter = TRUE           # Print meta-model training progress
  )
  
  # 3.4 Train elastic net meta-model (with cross-validation for lambda)
  lambda_seq <- 10^seq(-3, 3, length = 50)  # Fine-grained lambda sequence (regularization strength)
  
  # Fix: Correct outcome label formatting
  # 1. Prepare stacking data (use corrected labels)
  stack_data <- data.frame(train_stacking, Outcome = factor(train_Y)) # Directly use factor labels
  
  # 2. Check for NA or outliers in outcome column
  colSums(is.na(stack_data))    # Check NA counts
  table(stack_data$Outcome)     # Should show "No" and "Yes" counts
  
  # 3. Remove rows with NA (if any)
  stack_data <- na.omit(stack_data)
  
  # 3.5 Retrain stacking meta-model
  stack_model <- train(
    Outcome ~ .,
    data = stack_data,
    method = "glmnet",           # Elastic net (combines L1 and L2 regularization)
    family = "binomial",         # Binomial family for binary classification
    tuneGrid = expand.grid(
      alpha = seq(0, 1, 0.1),    # Tune alpha (0 = ridge, 1 = lasso, 0.5 = elastic net)
      lambda = lambda_seq        # Tune regularization strength
    ),
    metric = "ROC",              # Optimize meta-model for ROC/AUC
    trControl = stack_ctrl
  ) 
  
  # 3.6 Evaluate stacking ensemble
  # Prepare test set predictions for stacking
  test_stacking <- prepare_stacking_data(valid_models, test_X)
  test_stacking_df <- data.frame(test_stacking)
  
  # Safe prediction function (handles column name inconsistencies)
  safe_predict <- function(model, newdata) {
    pred <- predict(model, newdata = newdata, type = "prob")
    # Extract positive class probability (handles "X1" or "Yes")
    positive_col <- grep("X1|Yes", colnames(pred), value = TRUE)[1]
    return(pred[, positive_col])
  }
  
  # Predict probabilities with stacking meta-model
  stack_pred <- safe_predict(stack_model, test_stacking_df)
  
  # Evaluate stacking performance
  stack_roc <- roc(test_Y, stack_pred)
  stack_auc <- round(auc(stack_roc), 3)
  ci_auc <- round(ci.auc(stack_roc), 3)  # 95% CI for AUC
  
  cat("Optimized Stacking Ensemble Performance:\n")
  cat("AUC:", stack_auc, "95% CI:", paste(ci_auc, collapse = "-"), "\n")
  
  # --------------------------
  # 7. Final Performance Comparison
  # --------------------------
  # Combine base model and ensemble performance
  final_perf <- rbind(
    test_perf,
    data.frame(
      Model = c("Average_Ensemble", "Weighted_Ensemble","Stack_Ensemble"),
      AUC = c(avg_auc, weighted_auc, stack_auc),
      CI = NA  # Add CI calculation if needed
    )
  )
  
  # Sort by AUC (descending order)
  final_perf <- final_perf %>% arrange(desc(AUC))
  print(final_perf)  # Print final performance ranking
 
  
  ###### Save Models as RDS Files 
  # Create directory for saving models (if not exists)
  dir.create("model_rds", showWarnings = FALSE)
  # Save base models (only valid/trained models)
  for (model_name in names(valid_models)) {
    saveRDS(valid_models[[model_name]], 
            file = file.path("model_rds", paste0(model_name, ".rds")))
  }
  # Save stacking ensemble model (meta-model object)
  if (!is.null(stack_model)) {
    saveRDS(stack_model, file = file.path("model_rds", "stacked_ensemble.rds"))
  }
  # Save ensemble predictions (numeric vectors, not model objects)
  saveRDS(pred_df$Average, file = file.path("model_rds", "simple_average.rds"))
  saveRDS(pred_df$Weighted, file = file.path("model_rds", "weighted_average.rds")) # Fix: Use pred_df$Weighted instead of weighted_avg
  # Notes:
  # 1. Base models (max 4) saved as individual RDS files
  # 2. Stacking model saved as "stacked_ensemble.rds"
  # 3. Simple/weighted average predictions saved as numeric vectors
  # 4. All files stored in "model_rds" directory (current working directory)
  
  # Save standardization parameters (for new data preprocessing)
  train_means <- colMeans(train_X)  # Mean of training features
  train_sds <- apply(train_X, 2, sd)# SD of training features
  saveRDS(list(train_means = train_means, train_sds = train_sds), 
          file = "model_rds/standardization_params.rds")
  
 #### Visualization Code (Revised: Sep 5)
  library(pROC)      # ROC curve analysis
  library(PRROC)     # Precision-Recall curve analysis
  library(rmda)      # Decision Curve Analysis (DCA)
  library(ggplot2)   # Advanced plotting
  library(patchwork) # Multi-plot arrangement
  
  # Initialize list to store plot data
  plot_data <- list(
    roc_curves = list(),
    pr_curves = list(),
    dca_curves = list()
  )
  
  # Define color scheme for models (consistent across plots)
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
    pred <- predict(valid_models[[model_name]], test_X, type = "prob")[, "Yes"]
    plot_data$roc_curves[[model_name]] <- roc(test_Y, pred)
  }
  
  # Calculate ROC curves for ensemble models
  plot_data$roc_curves$Average_Ensemble <- roc(test_Y, pred_df$Average)
  plot_data$roc_curves$Weighted_Ensemble <- roc(test_Y, pred_df$Weighted)
  plot_data$roc_curves$Stack_Ensemble <- stack_roc
  
  ## Plot ROC Curves
  roc_plot <- ggplot() +
    # Add diagonal line (random guess, AUC=0.5)
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
  
  # Add AUC labels to legend
  auc_labels <- sapply(names(plot_data$roc_curves), function(x) {
    paste0(x, " (AUC = ", round(auc(plot_data$roc_curves[[x]]), 3), ")")
  })
  
  roc_plot <- roc_plot +
    scale_color_manual(values = model_colors,
                       labels = auc_labels)
  
  print(roc_plot)
  
  ## Fix: If plot not displayed, reset graphics device
  # Clear graphics device
  try(dev.off(), silent = TRUE)
  
  # Create new graphics device (8x6 inches, 300 DPI)
  dev.new(width = 8, height = 6, unit = "in", res = 300)
  
  # Reconstruct ROC plot with validation
  roc_plot <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "ROC Curve Comparison",
         x = "1 - Specificity", 
         y = "Sensitivity",
         color = "Model") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))  # Center title
  
  # Validate ROC data and add curves
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
  
  # Add legend labels with AUC
  auc_values <- sapply(plot_data$roc_curves, function(x) {
    if(!is.null(x)) round(auc(x), 3) else NA
  })
  legend_labels <- paste0(names(auc_values), " (AUC = ", auc_values, ")")
  
  # Final ROC plot
  roc_plot <- roc_plot +
    scale_color_manual(values = model_colors,
                       labels = legend_labels[!is.na(auc_values)])
  
  print(roc_plot) 
  
# Plot Precision-Recall (PR) Curves  
  # Convert factor labels to numeric (required by PRROC package: 1 = positive, 0 = negative)
  y_numeric <- as.numeric(test_Y) - 1
  
  # Calculate PR curves for base models
  for(model_name in names(valid_models)) {
    pred <- predict(valid_models[[model_name]], test_X, type = "prob")[, "Yes"]
    plot_data$pr_curves[[model_name]] <- pr.curve(
      scores.class0 = pred,       # Scores for positive class
      weights.class0 = y_numeric, # Binary labels (0/1)
      curve = TRUE                # Return full curve data
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
  
  # Plot PR curves
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
  
  # Add AUC-PR labels to legend
  prauc_labels <- sapply(names(plot_data$pr_curves), function(x) {
    paste0(x, " (AUC-PR = ", round(plot_data$pr_curves[[x]]$auc.integral, 3), ")")
  })
  
  pr_plot <- pr_plot +
    scale_color_manual(values = model_colors,
                       labels = prauc_labels)
  
  print(pr_plot)
  
# Plot Decision Curve Analysis (DCA) Curves  
  # Ensure rmda package is loaded
  library(rmda)
  
  # Safe method: Create DCA data frame with numeric labels (0 = No, 1 = Yes)
  dca_data <- data.frame(
    Outcome = ifelse(testData$CVD_CAD == "Yes", 1, 0), # Use original test labels
    
    # Add base model predictions
    GLM = predict(valid_models$GLM, test_X, type = "prob")[, "Yes"],
    RF = predict(valid_models$RF, test_X, type = "prob")[, "Yes"],
    XGB = predict(valid_models$XGB, test_X, type = "prob")[, "Yes"],
    SVM = predict(valid_models$SVM, test_X, type = "prob")[, "Yes"],
    KNN = predict(valid_models$KNN, test_X_knn, type = "prob")[, "Yes"],
    
    # Add ensemble model predictions
    Average = pred_df$Average,
    Weighted = pred_df$Weighted,
    Stacking = stack_pred
  )
  
  # Validate outcome distribution
  cat("DCA data outcome distribution:\n")
  table(dca_data$Outcome)
  cat("Unique values in outcome:", unique(dca_data$Outcome), "\n") # Should be 0 and 1
  
  # Check data summary
  summary(dca_data)
  
  # Check for NA values (remove if present)
  colSums(is.na(dca_data))
  dca_data <- na.omit(dca_data) 
  
  # Install and load dcurves package (for enhanced DCA plotting)
  if (!require(dcurves)) install.packages("dcurves")
  library(dcurves)
  
  # Run DCA (exclude KNN if it causes issues; adjust thresholds to [0, 0.4])
  dca_result_fixed <- dca(
    Outcome ~ GLM + RF + XGB + SVM + Average + Weighted + Stacking,
    data = dca_data,
    thresholds = seq(0, 0.4, by = 0.01)  # Threshold range for clinical relevance
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

# Updated DCA (Sep 5, 16:30)
  # Run full DCA with all models
  dca_result <- dca(
    Outcome ~ GLM + RF + XGB + SVM + KNN + Average + Weighted + Stacking,
    data = dca_data,
    thresholds = seq(0, 0.4, by = 0.01)
  )
  
  # Plot DCA with custom styling
  dca_plot <- plot(dca_result, smooth = TRUE) +
    # Custom color palette (matches ROC/PR plots)
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
        'All' = 'black',      # "Treat all" strategy color
        'None' = 'gray50'     # "Treat none" strategy color
      ),
      # Define legend order (prioritize best models first)
      breaks = c('Stacking', 'Weighted', 'Average', 'GLM', 'RF', 'XGB', 'SVM', 'KNN', 'All', 'None')
    ) +
    # Improve theme and labels
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
    # Adjust y-axis range to better display net benefit
    ggplot2::coord_cartesian(ylim = c(-0.1, 1.0))
  
  # Display DCA plot
  print(dca_plot)
  
  # Uncomment to save high-resolution DCA plot
  # ggplot2::ggsave("DCA_Full_With_Legend.png", dca_plot, 
  #                width = 12, height = 8, dpi = 300, bg = "white")
  
## Combine ROC, PR, and DCA Plots into One Figure  
  combined_plot <- (roc_plot + theme(legend.position = "none")) + 
    (pr_plot + theme(legend.position = "none")) + 
    (dca_plot + theme(legend.position = "right")) +
    plot_layout(ncol = 2, guides = "collect") +  # Arrange in 2 columns, share legend
    plot_annotation(title = "Comprehensive Model Performance Evaluation",
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))) # Center main title
  
  print(combined_plot) 
  

## Plot Calibration Curves (Optional: Comment out if not needed)
  # Load required packages
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Fix test set label encoding (match training set)
  test_Y <- factor(test_Y, levels = c("X0", "X1"))
  levels(test_Y) <- c("No", "Yes")  # Map X0 -> No, X1 -> Yes
  
  # Alternative: Re-extract labels from original test data
  test_Y <- testData$CVD_CAD
  
  # Validate corrected label distribution
  cat("Corrected test set label distribution:\n")
  print(table(test_Y))
  cat("Numeric conversion result:\n")
  test_Y_numeric <- ifelse(test_Y == "Yes", 1, 0)
  print(table(test_Y_numeric))
  
  # 1. Build calibration data frame ------------------------------------------------------------
  # Create data frame with predicted probabilities and actual labels
  calibration_data <- data.frame(
    Outcome = ifelse(test_Y == "Yes", 1, 0),  # Convert target to numeric (0/1)
    
    # Base model predicted probabilities
    GLM = predict(valid_models$GLM, test_X, type = "prob")[, "Yes"],
    RF = predict(valid_models$RF, test_X, type = "prob")[, "Yes"],
    XGB = predict(valid_models$XGB, test_X, type = "prob")[, "Yes"],
    SVM = predict(valid_models$SVM, test_X, type = "prob")[, "Yes"],
    KNN = predict(valid_models$KNN, test_X_knn, type = "prob")[, "Yes"],
    
    # Ensemble model predicted probabilities
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
  
  # 2. Convert data to long format (for ggplot2) --------------------------------------------------------
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
      loess(Outcome ~ Predicted_Probability, span = 0.75)  # Span = 0.75 (balance smoothness/fidelity)
    )) %>%
    ungroup()
  
  # 4. Plot calibration curves ----------------------------------------------------------
  cal_plot <- ggplot(cal_curve_data, aes(x = Predicted_Probability, y = Smoothed_Observed, color = Model)) +
    geom_line(linewidth = 1.2) +
    # Add ideal calibration line (predicted probability = observed proportion)
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    
    # Use custom color scheme
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
    
    # Set equal axis scales and ranges
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "Model Calibration Curves",
      x = "Predicted Probability",
      y = "Smoothed Observed Proportion",
      color = "Model",
      caption = "Dashed line = ideal calibration (predicted = observed)"
    ) +
    
    # Improve theme
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
  
  # 6. (Optional) Add binned observed proportions ------------------------------------------------------
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
               alpha = 0.6, shape = 16) +  # Semi-transparent circles
    scale_size_continuous(range = c(2, 6), name = "Sample Size") +  # Adjust point size range
    labs(title = "Model Calibration Curves (With Binned Observations)")
  
  print(cal_plot_with_points)