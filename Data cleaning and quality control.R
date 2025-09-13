# Load required libraries
library(haven)
library(dplyr)
library(purrr)

############################
### 1. Data loading & initial integration ###
############################

# ----- 1.1 Define data-reading function -----
read_nhanes_xpt <- function(file_prefix, cycles = c("F", "G", "H", "I", "J"), select_vars = NULL) {
  # Purpose: Batch-read XPT files from multiple cycles and bind them
  # Args:
  #   file_prefix - file prefix (e.g. DEMO_)
  #   cycles      - cycle suffixes
  #   select_vars - variables to keep
  map_dfr(cycles, ~{
    file <- paste0(file_prefix, "_", .x, ".XPT")
    data <- read_xpt(file)
    if (!is.null(select_vars)) data <- data %>% select(any_of(select_vars))
    data
  })
}

# ----- 1.2 Read each module -----
# Demographics
demo_data <- read_nhanes_xpt("DEMO", 
                             select_vars = c("SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1", 
                                             "INDFMPIR", "WTMEC2YR", "SDMVSTRA", "SDMVPSU"))

# Coronary heart-disease outcome (MCQ160C)
mcq_data <- read_nhanes_xpt("MCQ", select_vars = c("SEQN", "MCQ160C"))

# Lipid & glucose data
lipid_data <- list(
  hdl = read_nhanes_xpt("HDL", select_vars = c("SEQN", "LBDHDD")),
  tchol = read_nhanes_xpt("TCHOL", select_vars = c("SEQN", "LBXTC")),
  trigly = read_nhanes_xpt("TRIGLY", select_vars = c("SEQN", "LBDLDL", "LBXTR")),
  glu = read_nhanes_xpt("GLU", select_vars = c("SEQN", "LBXGLU"))
) %>% reduce(full_join, by = "SEQN")

# Covariate data
covar_data <- list(
  bmx = read_nhanes_xpt("BMX", select_vars = c("SEQN", "BMXBMI")),
  bpx = read_nhanes_xpt("BPX", select_vars = c("SEQN", "BPXSY1", "BPXDI1")),
  smq = read_nhanes_xpt("SMQ", select_vars = c("SEQN", "SMQ020")),
  alq = read_nhanes_xpt("ALQ", select_vars = c("SEQN", "ALQ101"))
) %>% reduce(full_join, by = "SEQN")

# ----- 1.3 Merge all data -----
merged_data <- list(demo_data, mcq_data, lipid_data, covar_data) %>% 
  reduce(full_join, by = "SEQN")

################################
### 2. Systematic data-quality workflow ###
################################
# ----- 2.1 Initial sample filtering -----
step1_data <- merged_data %>%
  filter(
    RIDAGEYR >= 18,                  # keep adults
    MCQ160C %in% c(1, 2)             # definite CAD status
  )

# Exclude samples with SMQ020 == 9
step1_1_data <- step1_data %>%
  filter(SMQ020 != 9)

# Exclude samples with ALQ101 == 9
step1_2_data <- step1_1_data %>%
  filter(ALQ101 != 9)

# ----- 2.2 Missing-value handling -----
step2_data <- step1_2_data %>%
  filter(
    !is.na(RIDAGEYR),                # age
    !is.na(RIAGENDR),                # sex
    !is.na(BMXBMI),                  # BMI
    !is.na(BPXSY1) & !is.na(BPXDI1), # blood pressure
    !is.na(LBXGLU),                  # glucose
    !is.na(LBXTR),                   # triglycerides
    !is.na(SMQ020),                  # smoking
    !is.na(ALQ101)                   # alcohol
  )

# ----- 2.3 Outlier treatment -----
# (1) Hard clinical-range filters
step3_1_data <- step2_data %>%
  mutate(
    BMXBMI = ifelse(between(BMXBMI, 15, 50), BMXBMI, NA_real_),
    LBXGLU = ifelse(between(LBXGLU, 70, 300), LBXGLU, NA_real_),
    LBXTR = ifelse(between(LBXTR, 30, 1000), LBXTR, NA_real_),
    BPXSY1 = ifelse(between(BPXSY1, 90, 200), BPXSY1, NA_real_),
    BPXDI1 = ifelse(between(BPXDI1, 40, 120), BPXDI1, NA_real_),
  )

# (2) IQR filter (1.5 × IQR)
calculate_iqr_filter <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lower <- q[1] - 1.5 * iqr
  upper <- q[2] + 1.5 * iqr
  ifelse(x >= lower & x <= upper, x, NA)
}

step3_2_data <- step3_1_data %>%
  mutate(
    RIDAGEYR = calculate_iqr_filter(RIDAGEYR),
    BMXBMI = calculate_iqr_filter(BMXBMI),
    LBXGLU = calculate_iqr_filter(LBXGLU),
    LBXTR = calculate_iqr_filter(LBXTR),
    BPXSY1 = calculate_iqr_filter(BPXSY1),
    BPXDI1 = calculate_iqr_filter(BPXDI1)
  )

step3_data <- step3_2_data %>%
  filter(complete.cases(.))       # drop rows with NA introduced by outlier removal

# ----- 2.35 Compute TyG index (new step) -----
step3_3_data <- step3_data %>%
  mutate(
    TyG = log((LBXTR * LBXGLU) / 2)
  )

# ----- 2.4 Standardisation -----
step4_data <- step3_3_data %>%
  mutate(
    RIDAGEYR = as.numeric(scale(RIDAGEYR)),
    BMXBMI = as.numeric(scale(BMXBMI)),
    LBXGLU = as.numeric(scale(LBXGLU)),
    LBXTR = as.numeric(scale(LBXTR)),
    BPXSY1 = as.numeric(scale(BPXSY1)),
    BPXDI1 = as.numeric(scale(BPXDI1)),
    TyG = as.numeric(scale(TyG))
  )

# ----- 2.5 Encode categorical variables -----
step5_data <- step4_data %>%
  mutate(
    RIAGENDR = factor(RIAGENDR, labels = c("Male", "Female")),
    SMQ020 = factor(SMQ020, labels = c("Non-smoker", "Smoker")),
    ALQ101 = factor(ALQ101, labels = c("Non-drinker", "Drinker")),
    CVD_CAD = factor(MCQ160C, labels = c("No", "Yes")),
    TyG_quartile = ntile(TyG, 4) %>% factor()
  )
step5_1_data <- step5_data %>%
  select(-MCQ160C)

cleaned_data <- step5_1_data

# ----- 2.6 Save cleaned dataset -----
write.csv(cleaned_data, "cleaned_CAD_data.csv", row.names = FALSE)

################################
### 3. Quality-control report ###
################################

# Track sample-size changes
qc_report <- data.frame(
  Stage = c("Raw Data", "Age≥18 & Clear CAD", "Remove SMQ020 = 9", "Remove ALQ101 = 9", "Missing Values Removed", 
            "Outliers Removed", "Calculate TyG", "Standardization", "Categorical Encoding", "Remove MCQ160C", "Final Cleaned Data"),
  N = c(nrow(merged_data),
        nrow(step1_data),
        nrow(step1_1_data),
        nrow(step1_2_data),
        nrow(step2_data),
        nrow(step3_data),
        nrow(step3_3_data),
        nrow(step4_data),
        nrow(step5_data),
        nrow(step5_1_data),
        nrow(cleaned_data))
)

print(qc_report)

# Summary of final variables
summary(cleaned_data)

# ----- 2.7 Normality tests (large sample) -----
# Load required packages
suppressPackageStartupMessages({
  if (!require(ggplot2)) install.packages("ggplot2")
  if (!require(ggpubr)) install.packages("ggpubr")
  if (!require(gridExtra)) install.packages("gridExtra")
  if (!require(nortest)) install.packages("nortest") # multiple normality tests
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
  library(nortest)
})

# 1. Select continuous variables to test
continuous_vars <- c("RIDAGEYR", "BMXBMI", "LBXGLU", "LBXTR", "BPXSY1", "BPXDI1", "TyG")

# 2. Run large-sample normality tests and collect results
# Lilliefors (K-S corrected) and Anderson-Darling tests
normality_test_results <- lapply(continuous_vars, function(var) {
  lillie_result <- lillie.test(cleaned_data[[var]])
  ad_result <- ad.test(cleaned_data[[var]])
  
  data.frame(
    Variable = var,
    Lillie_p_value = round(lillie_result$p.value, 4),
    AD_p_value = round(ad_result$p.value, 4),
    Lillie_Normality = ifelse(lillie_result$p.value > 0.05, "Pass", "Fail"),
    AD_Normality = ifelse(ad_result$p.value > 0.05, "Pass", "Fail")
  )
}) %>% bind_rows()

cat("=== Large-sample normality tests (α = 0.05) ===\n")
print(normality_test_results)

# 3. Visual inspection (histograms + Q-Q plots)
create_hist_qq <- function(data, var_name) {
  p_hist <- ggplot(data, aes_string(x = var_name)) +
    geom_histogram(aes(y = ..density..), bins = 50, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(color = "red", linewidth = 1) +
    stat_function(fun = dnorm, args = list(mean = mean(data[[var_name]], na.rm = TRUE), 
                                           sd = sd(data[[var_name]], na.rm = TRUE)), 
                  color = "darkred", linetype = "dashed", linewidth = 1) +
    ggtitle(paste("Histogram of", var_name)) +
    theme_minimal()
  
  p_qq <- ggqqplot(data, x = var_name, color = "blue", alpha = 0.5) +
    ggtitle(paste("Q-Q Plot of", var_name)) +
    theme_minimal()
  
  list(hist = p_hist, qq = p_qq)
}

normality_plots <- lapply(continuous_vars, function(v) create_hist_qq(cleaned_data, v))

# 4. Arrange plots (2 variables per row)
for (i in seq(1, length(continuous_vars), by = 2)) {
  var1 <- continuous_vars[i]
  var2 <- ifelse(i+1 <= length(continuous_vars), continuous_vars[i+1], NULL)
  
  if (!is.null(var2)) {
    grid.arrange(
      normality_plots[[var1]]$hist, normality_plots[[var2]]$hist,
      normality_plots[[var1]]$qq,   normality_plots[[var2]]$qq,
      ncol = 2,
      top = paste("Normality checks:", var1, "&", var2)
    )
  } else {
    grid.arrange(
      normality_plots[[var1]]$hist, normality_plots[[var1]]$qq,
      ncol = 2,
      top = paste("Normality checks:", var1)
    )
  }
}

# 5. Save test results
write.csv(normality_test_results, "normality_test_large_sample_results.csv", row.names = FALSE)
cat("\nLarge-sample normality-test results saved to 'normality_test_large_sample_results.csv'\n")

# 6. Brief recommendations based on normality
cat("\n=== Preliminary recommendations based on normality tests ===\n")
non_normal_vars_lillie <- normality_test_results$Variable[normality_test_results$Lillie_Normality == "Fail"]
non_normal_vars_ad <- normality_test_results$Variable[normality_test_results$AD_Normality == "Fail"]
non_normal_vars <- union(non_normal_vars_lillie, non_normal_vars_ad)

if (length(non_normal_vars) > 0) {
  cat("The following variables may deviate from normality; consider:\n")
  cat(paste("-", non_normal_vars, collapse = "\n"), "\n")
  cat("  • Use non-parametric tests (e.g. Mann-Whitney U, Kruskal-Wallis)\n")
  cat("  • Apply transformations (log, square-root, etc.)\n")
  cat("  • With very large samples, parametric tests may still be robust; combine graphical and domain knowledge to decide\n")
} else {
  cat("All continuous variables appear normally distributed; parametric tests can be used.\n")
}

# 7. Extra: skewness & kurtosis for large samples
cat("\n=== Skewness & kurtosis (|skew| > 3 or kurtosis > 10 may indicate severe non-normality) ===\n")
library(moments)
skew_kurtosis <- lapply(continuous_vars, function(var) {
  data_vector <- cleaned_data[[var]]
  data.frame(
    Variable = var,
    Skewness = round(skewness(data_vector), 4),
    Kurtosis = round(kurtosis(data_vector), 4)
  )
}) %>% bind_rows()
print(skew_kurtosis)