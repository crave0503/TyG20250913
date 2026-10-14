# =============================================================
# Association of the triglyceride-glucose (TyG) index with coronary
# artery disease (CAD) -- complete analysis pipeline
#
# Cohorts
#   Derivation cohort : NHANES 2009-2018 (cycles F, G, H, I, J)
#   Validation cohort : NHANES 2005-2008 (cycles D, E), a temporally
#                       distinct sample drawn from the same national
#                       survey programme
#
# Sections (run in order 1 -> 13; Section 14 can be run at any time)
#    1. Derivation cohort: cleaning and quality control
#    2. Validation cohort: cleaning and quality control
#    3. Baseline characteristics and survey-weighted univariate regression
#    4. Normality tests
#    5. Machine-learning models (Full / TyG-Ablated / TyG-Only)
#    6. Temporally distinct validation
#    7. SHAP interpretability analysis
#    8. Derivation vs validation distribution comparison
#    9. Participant flow diagram
#   10. Supplementary Table 1: missingness proportions
#   11. Supplementary Table 2: AUC-PR with bootstrap 95% CI
#   12. Supplementary Table 4: Hosmer-Lemeshow calibration
#   13. Session information (software and package versions)
#
# Software: R 4.3.2. Package versions are recorded by Section 13
#           (see also sessionInfo.txt).
#
# Input data: download the required .XPT files from the CDC NHANES
#   website (https://wwwn.cdc.gov/nchs/nhanes/) and place them in the
#   folders configured in Section 0. File names follow the standard
#   NHANES convention, e.g. DEMO_F.XPT, MCQ_F.XPT, ... (derivation)
#   and DEMO_D.XPT, MCQ_D.XPT, ... (validation).
# =============================================================


# =============================================================
# 0. Global configuration  (EDIT THESE PATHS BEFORE RUNNING)
# =============================================================
# Folder holding the derivation-cohort XPT files (NHANES 2009-2018)
XPT_DIR_DERIVATION <- "data/nhanes_2009_2018"
# Folder holding the validation-cohort XPT files (NHANES 2005-2008)
XPT_DIR_VALIDATION <- "data/nhanes_2005_2008"
# Folder where every table, figure and model object will be written
OUT_DIR <- "output"

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Outcome-coding convention used throughout (NHANES raw codes):
#   MCQ160C : 1 = Yes (CAD), 2 = No
#   SMQ020  : 1 = Yes (ever smoker), 2 = No
#   ALQ101  : 1 = Yes (ever drinker), 2 = No
#   RIAGENDR: 1 = Male, 2 = Female
# All factor() calls below pass explicit levels = c(1, 2) so that labels
# are attached in the correct direction.


# =============================================================
# 1. Derivation cohort: data cleaning and quality control
#    Input : XPT files in XPT_DIR_DERIVATION
#    Output: cleaned_CAD_data_fixed.csv            (standardised, modelling)
#            cleaned_CAD_data_original_scale.csv   (raw scale, descriptive)
#            qc_report_train.csv                   (stepwise N, flow diagram)
# =============================================================
pkgs <- c("haven","dplyr","purrr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
}

# ---- 1.1 Read and merge raw data ----
read_nhanes_xpt <- function(file_prefix, cycles = c("F","G","H","I","J"), select_vars = NULL) {
  map_dfr(cycles, ~{
    data <- read_xpt(file.path(XPT_DIR_DERIVATION, paste0(file_prefix, "_", .x, ".XPT")))
    if (!is.null(select_vars)) data <- data %>% select(any_of(select_vars))
    data
  })
}

demo_data  <- read_nhanes_xpt("DEMO",
              select_vars = c("SEQN","RIDAGEYR","RIAGENDR","RIDRETH1","INDFMPIR",
                              "WTMEC2YR","SDMVSTRA","SDMVPSU"))
mcq_data   <- read_nhanes_xpt("MCQ", select_vars = c("SEQN","MCQ160C"))
lipid_data <- list(
  hdl    = read_nhanes_xpt("HDL",    select_vars = c("SEQN","LBDHDD")),
  tchol  = read_nhanes_xpt("TCHOL",  select_vars = c("SEQN","LBXTC")),
  trigly = read_nhanes_xpt("TRIGLY", select_vars = c("SEQN","LBDLDL","LBXTR")),
  glu    = read_nhanes_xpt("GLU",    select_vars = c("SEQN","LBXGLU"))
) %>% reduce(full_join, by = "SEQN")
covar_data <- list(
  bmx = read_nhanes_xpt("BMX", select_vars = c("SEQN","BMXBMI")),
  bpx = read_nhanes_xpt("BPX", select_vars = c("SEQN","BPXSY1","BPXDI1")),
  smq = read_nhanes_xpt("SMQ", select_vars = c("SEQN","SMQ020")),
  alq = read_nhanes_xpt("ALQ", select_vars = c("SEQN","ALQ101"))
) %>% reduce(full_join, by = "SEQN")
merged_data <- list(demo_data, mcq_data, lipid_data, covar_data) %>% reduce(full_join, by = "SEQN")
cat("=== 1. Merged raw sample size:", nrow(merged_data), "===\n")

# ---- 1.2 Quality control (stepwise) ----
qc <- list()
qc$raw <- nrow(merged_data)

# 1.2.1 Adults with an unambiguous CAD diagnosis (MCQ160C 1 = Yes, 2 = No)
step1_data <- merged_data %>%
  filter(RIDAGEYR >= 18, MCQ160C %in% c(1, 2))
qc$step1 <- nrow(step1_data)

# 1.2.2 Remove "don't know" (9) replies to the smoking and alcohol questions
step1_1_data <- step1_data %>% filter(SMQ020 != 9)
qc$step1_1 <- nrow(step1_1_data)
step1_2_data <- step1_1_data %>% filter(ALQ101 != 9)
qc$step1_2 <- nrow(step1_2_data)

# 1.2.3 Missing-value handling (complete cases)
step2_data <- step1_2_data %>%
  filter(!is.na(RIDAGEYR), !is.na(RIAGENDR), !is.na(BMXBMI),
         !is.na(BPXSY1) & !is.na(BPXDI1),
         !is.na(LBXGLU), !is.na(LBXTR), !is.na(SMQ020), !is.na(ALQ101))
qc$step2 <- nrow(step2_data)

# 1.2.4 Implausible values: clinical-threshold replacement followed by IQR
#       filtering, applied in a single complete.cases step. Splitting this
#       into two complete.cases steps would change the IQR quantile basis
#       and shift the final N away from 6062.
step3_1_data <- step2_data %>%
  mutate(
    BMXBMI = ifelse(between(BMXBMI, 15, 50), BMXBMI, NA_real_),
    LBXGLU = ifelse(between(LBXGLU, 70, 300), LBXGLU, NA_real_),
    LBXTR  = ifelse(between(LBXTR, 30, 1000), LBXTR, NA_real_),
    BPXSY1 = ifelse(between(BPXSY1, 90, 200), BPXSY1, NA_real_),
    BPXDI1 = ifelse(between(BPXDI1, 40, 120), BPXDI1, NA_real_)
  )

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
    BMXBMI   = calculate_iqr_filter(BMXBMI),
    LBXGLU   = calculate_iqr_filter(LBXGLU),
    LBXTR    = calculate_iqr_filter(LBXTR),
    BPXSY1   = calculate_iqr_filter(BPXSY1),
    BPXDI1   = calculate_iqr_filter(BPXDI1)
  )
step3_data <- step3_2_data %>% filter(complete.cases(.))
qc$step3 <- nrow(step3_data)

# Breakdown of the two exclusion types (used only for the flow diagram;
# the filtering itself remains a single combined step)
n_clin_excl <- nrow(step2_data) - sum(complete.cases(step3_1_data))
qc$step3_1b <- nrow(step2_data) - n_clin_excl
qc$step3_1b_excl <- n_clin_excl
qc$step3_iqr_excl <- (nrow(step2_data) - nrow(step3_data)) - n_clin_excl

# 1.2.5 TyG index
step3_3_data <- step3_data %>%
  mutate(TyG = log((LBXTR * LBXGLU) / 2))

# ---- 1.3 Raw-scale version (for descriptive statistics) ----
orig_data <- step3_3_data %>%
  mutate(
    RIAGENDR = factor(RIAGENDR, levels = c(1,2), labels = c("Male","Female")),
    SMQ020   = factor(SMQ020,   levels = c(1,2), labels = c("Smoker","Non-smoker")),
    ALQ101   = factor(ALQ101,   levels = c(1,2), labels = c("Drinker","Non-drinker")),
    CVD_CAD  = factor(MCQ160C,  levels = c(1,2), labels = c("Yes","No")),
    TyG_quartile = ntile(TyG, 4) %>% factor()
  ) %>%
  select(-MCQ160C)
write.csv(orig_data, file.path(OUT_DIR, "cleaned_CAD_data_original_scale.csv"), row.names = FALSE)

# ---- 1.4 Standardised version (for modelling) ----
std_data <- step3_3_data %>%
  mutate(
    RIDAGEYR = as.numeric(scale(RIDAGEYR)),
    BMXBMI   = as.numeric(scale(BMXBMI)),
    LBXGLU   = as.numeric(scale(LBXGLU)),
    LBXTR    = as.numeric(scale(LBXTR)),
    BPXSY1   = as.numeric(scale(BPXSY1)),
    BPXDI1   = as.numeric(scale(BPXDI1)),
    TyG      = as.numeric(scale(TyG))
  ) %>%
  mutate(
    RIAGENDR = factor(RIAGENDR, levels = c(1,2), labels = c("Male","Female")),
    SMQ020   = factor(SMQ020,   levels = c(1,2), labels = c("Smoker","Non-smoker")),
    ALQ101   = factor(ALQ101,   levels = c(1,2), labels = c("Drinker","Non-drinker")),
    CVD_CAD  = factor(MCQ160C,  levels = c(1,2), labels = c("Yes","No")),
    TyG_quartile = ntile(TyG, 4) %>% factor()
  ) %>%
  select(-MCQ160C)
write.csv(std_data, file.path(OUT_DIR, "cleaned_CAD_data_fixed.csv"), row.names = FALSE)

# ---- 1.5 QC report (drives the participant flow diagram) ----
qc_report <- data.frame(
  Stage = c("Raw merged data",
            "Age>=18 & MCQ160C in 1,2",
            "Remove SMQ020 = 9 and ALQ101 = 9",
            "Missing values removed",
            "Clinically implausible values removed",
            "IQR outlier filtering (final)"),
  N = c(qc$raw, qc$step1, qc$step1_2, qc$step2, qc$step3_1b, qc$step3),
  Excluded = c(NA, qc$raw - qc$step1, qc$step1 - qc$step1_2,
               qc$step1_2 - qc$step2, qc$step3_1b_excl, qc$step3_iqr_excl),
  Reason = c("Total NHANES 2009-2018",
             "Age < 18 or no clear CAD diagnosis (MCQ160C not 1/2)",
             "Smoking/alcohol status unknown (9)",
             "Missing values in any key variable",
             "Values outside clinically plausible thresholds",
             "Values outside 1.5 x IQR range")
)
write.csv(qc_report, file.path(OUT_DIR, "qc_report_train.csv"), row.names = FALSE)
print(qc_report)

# ---- 1.6 Outcome distribution check ----
cat("\n===== Derivation cohort: outcome distribution =====\n")
print(table(std_data$CVD_CAD))
cat("CAD prevalence:", round(mean(std_data$CVD_CAD == "Yes") * 100, 2), "%\n")
cat("Expected: n = 6062, Yes (CAD) = 206, No = 5856\n")

cat("\n[OK] Derivation cohort outputs written:\n")
cat("  cleaned_CAD_data_fixed.csv           (standardised, modelling)\n")
cat("  cleaned_CAD_data_original_scale.csv  (raw scale, descriptive)\n")
cat("  qc_report_train.csv                  (stepwise N, flow diagram)\n")


# =============================================================
# 2. Validation cohort: data cleaning and quality control
#    Input : XPT files in XPT_DIR_VALIDATION  (NHANES 2005-2008, D-E)
#    Output: cleaned_CAD_externalvalidation_data_fixed.csv
#            cleaned_CAD_externalvalidation_data_original_scale.csv
#            qc_report_valid.csv
#    Note  : no BMI and no alcohol variable, matching the original
#            validation-cohort definition.
# =============================================================
pkgs <- c("haven","dplyr","purrr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
}

read_nhanes_xpt_val <- function(file_prefix, cycles = c("D","E"), select_vars = NULL) {
  map_dfr(cycles, ~{
    data <- read_xpt(file.path(XPT_DIR_VALIDATION, paste0(file_prefix, "_", .x, ".XPT")))
    if (!is.null(select_vars)) data <- data %>% select(any_of(select_vars))
    data
  })
}

demo_data  <- read_nhanes_xpt_val("DEMO", select_vars = c("SEQN","RIDAGEYR","RIAGENDR"))
mcq_data   <- read_nhanes_xpt_val("MCQ",  select_vars = c("SEQN","MCQ160C"))
lipid_data <- list(
  tchol  = read_nhanes_xpt_val("TCHOL",  select_vars = c("SEQN","LBXTC")),
  trigly = read_nhanes_xpt_val("TRIGLY", select_vars = c("SEQN","LBDLDL","LBXTR")),
  glu    = read_nhanes_xpt_val("GLU",    select_vars = c("SEQN","LBXGLU"))
) %>% reduce(full_join, by = "SEQN")
covar_data <- list(
  bpx = read_nhanes_xpt_val("BPX", select_vars = c("SEQN","BPXSY1","BPXDI1")),
  smq = read_nhanes_xpt_val("SMQ", select_vars = c("SEQN","SMQ020"))
) %>% reduce(full_join, by = "SEQN")
merged_data <- list(demo_data, mcq_data, lipid_data, covar_data) %>% reduce(full_join, by = "SEQN")
cat("=== 2. Merged raw sample size:", nrow(merged_data), "===\n")

qc <- list()
qc$raw <- nrow(merged_data)

# 2.1 Adults with an unambiguous CAD diagnosis
step1_data <- merged_data %>%
  filter(RIDAGEYR >= 18, MCQ160C %in% c(1, 2))
qc$step1 <- nrow(step1_data)

# 2.2 Keep only valid smoking replies (1/2); refuse (7) and don't know (9) removed
step1_1_data <- step1_data %>% filter(SMQ020 %in% c(1, 2))
qc$step2 <- nrow(step1_1_data)

# 2.3 Missing values
step2_data <- step1_1_data %>%
  filter(!is.na(RIDAGEYR), !is.na(RIAGENDR),
         !is.na(BPXSY1) & !is.na(BPXDI1),
         !is.na(LBXGLU), !is.na(LBXTR), !is.na(SMQ020))
qc$step3 <- nrow(step2_data)

# 2.4 Implausible values: clinical thresholds + IQR, single complete.cases step
step3_1_data <- step2_data %>%
  mutate(
    LBXGLU = ifelse(between(LBXGLU, 70, 300), LBXGLU, NA_real_),
    LBXTR  = ifelse(between(LBXTR, 30, 1000), LBXTR, NA_real_),
    BPXSY1 = ifelse(between(BPXSY1, 90, 200), BPXSY1, NA_real_),
    BPXDI1 = ifelse(between(BPXDI1, 40, 120), BPXDI1, NA_real_)
  )
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
    LBXGLU   = calculate_iqr_filter(LBXGLU),
    LBXTR    = calculate_iqr_filter(LBXTR),
    BPXSY1   = calculate_iqr_filter(BPXSY1),
    BPXDI1   = calculate_iqr_filter(BPXDI1)
  )
step3_data <- step3_2_data %>% filter(complete.cases(.))
qc$step5 <- nrow(step3_data)

n_clin_excl <- nrow(step2_data) - sum(complete.cases(step3_1_data))
qc$step4 <- nrow(step2_data) - n_clin_excl
qc$step4_excl <- n_clin_excl
qc$step5_iqr_excl <- (nrow(step2_data) - nrow(step3_data)) - n_clin_excl

# 2.5 TyG index
step3_3_data <- step3_data %>%
  mutate(TyG = log((LBXTR * LBXGLU) / 2))

# ---- 2.6 Raw-scale version ----
orig_data <- step3_3_data %>%
  mutate(
    RIAGENDR = factor(RIAGENDR, levels = c(1,2), labels = c("Male","Female")),
    SMQ020   = factor(SMQ020,   levels = c(1,2), labels = c("Smoker","Non-smoker")),
    CVD_CAD  = factor(MCQ160C,  levels = c(1,2), labels = c("Yes","No")),
    TyG_quartile = ntile(TyG, 4) %>% factor()
  ) %>%
  select(-MCQ160C)
write.csv(orig_data, file.path(OUT_DIR, "cleaned_CAD_externalvalidation_data_original_scale.csv"), row.names = FALSE)

# ---- 2.7 Standardised version ----
std_data <- step3_3_data %>%
  mutate(
    RIDAGEYR = as.numeric(scale(RIDAGEYR)),
    LBXGLU   = as.numeric(scale(LBXGLU)),
    LBXTR    = as.numeric(scale(LBXTR)),
    BPXSY1   = as.numeric(scale(BPXSY1)),
    BPXDI1   = as.numeric(scale(BPXDI1)),
    TyG      = as.numeric(scale(TyG))
  ) %>%
  mutate(
    RIAGENDR = factor(RIAGENDR, levels = c(1,2), labels = c("Male","Female")),
    SMQ020   = factor(SMQ020,   levels = c(1,2), labels = c("Smoker","Non-smoker")),
    CVD_CAD  = factor(MCQ160C,  levels = c(1,2), labels = c("Yes","No")),
    TyG_quartile = ntile(TyG, 4) %>% factor()
  ) %>%
  select(-MCQ160C)
write.csv(std_data, file.path(OUT_DIR, "cleaned_CAD_externalvalidation_data_fixed.csv"), row.names = FALSE)

# ---- 2.8 QC report ----
qc_report <- data.frame(
  Stage = c("Raw merged data",
            "Age>=18 & MCQ160C in 1,2",
            "Remove SMQ020 = 7 and 9",
            "Missing values removed",
            "Clinically implausible values removed",
            "IQR outlier filtering (final)"),
  N = c(qc$raw, qc$step1, qc$step2, qc$step3, qc$step4, qc$step5),
  Excluded = c(NA, qc$raw - qc$step1, qc$step1 - qc$step2,
               qc$step2 - qc$step3, qc$step4_excl, qc$step5_iqr_excl),
  Reason = c("Total NHANES 2005-2008",
             "Age < 18 or no clear CAD diagnosis (MCQ160C not 1/2)",
             "Smoking status refused (7) or unknown (9)",
             "Missing values in any key variable",
             "Values outside clinically plausible thresholds",
             "Values outside 1.5 x IQR range")
)
write.csv(qc_report, file.path(OUT_DIR, "qc_report_valid.csv"), row.names = FALSE)
print(qc_report)

cat("\n===== Validation cohort: outcome distribution =====\n")
print(table(std_data$CVD_CAD))
cat("CAD prevalence:", round(mean(std_data$CVD_CAD == "Yes") * 100, 2), "%\n")
cat("Expected: Yes (CAD) = 128, No = 3366\n")

cat("\n[OK] Validation cohort outputs written:\n")
cat("  cleaned_CAD_externalvalidation_data_fixed.csv\n")
cat("  cleaned_CAD_externalvalidation_data_original_scale.csv\n")
cat("  qc_report_valid.csv\n")


# =============================================================
# 3. Baseline characteristics and survey-weighted univariate regression
#    Input : cleaned_CAD_data_fixed.csv (standardised)
#            cleaned_CAD_data_original_scale.csv (raw scale)
#    Output: Table1_univariate.csv            (weighted univariate OR, 95% CI, p)
#            baseline_characteristics.csv     (overall / CAD / no-CAD columns)
# =============================================================
pkgs <- c("survey","dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
}
std_data <- read.csv(file.path(OUT_DIR, "cleaned_CAD_data_fixed.csv"), stringsAsFactors = FALSE)
orig_data <- read.csv(file.path(OUT_DIR, "cleaned_CAD_data_original_scale.csv"), stringsAsFactors = FALSE)

# Set factor reference levels
std_data$CVD_CAD  <- factor(std_data$CVD_CAD,  levels = c("No","Yes"))
std_data$RIAGENDR <- factor(std_data$RIAGENDR, levels = c("Female","Male"))   # reference = Female
std_data$SMQ020   <- factor(std_data$SMQ020,   levels = c("Non-smoker","Smoker"))
std_data$ALQ101   <- factor(std_data$ALQ101,   levels = c("Non-drinker","Drinker"))
std_data$RIDRETH1 <- factor(std_data$RIDRETH1)

cat("Outcome distribution (full derivation cohort):\n"); print(table(std_data$CVD_CAD))

# ---- 3.1 Baseline characteristics table (raw scale) ----
orig_data$CVD_CAD <- factor(orig_data$CVD_CAD, levels = c("No","Yes"))
continuous_vars <- c("RIDAGEYR","BMXBMI","BPXSY1","BPXDI1","LBXTR","LBXGLU",
                     "LBDHDD","LBXTC","LBDLDL","TyG")
categorical_vars <- c("RIAGENDR","SMQ020","ALQ101","RIDRETH1","INDFMPIR")

summary_one <- function(x) {
  paste0(round(median(x, na.rm=TRUE), 2), " [", round(quantile(x, 0.25, na.rm=TRUE), 2),
         ", ", round(quantile(x, 0.75, na.rm=TRUE), 2), "]")
}
baseline <- data.frame(Variable = character(), Overall = character(),
                       CAD_Yes = character(), CAD_No = character(),
                       p_value = numeric(), stringsAsFactors = FALSE)
for (v in continuous_vars) {
  if (!(v %in% names(orig_data))) next
  pv <- tryCatch(wilcox.test(orig_data[[v]] ~ orig_data$CVD_CAD)$p.value, error = function(e) NA)
  baseline <- rbind(baseline, data.frame(
    Variable = v,
    Overall = summary_one(orig_data[[v]]),
    CAD_Yes = summary_one(orig_data[[v]][orig_data$CVD_CAD == "Yes"]),
    CAD_No  = summary_one(orig_data[[v]][orig_data$CVD_CAD == "No"]),
    p_value = pv, stringsAsFactors = FALSE))
}
for (v in categorical_vars) {
  if (!(v %in% names(orig_data))) next
  tab <- table(orig_data[[v]], orig_data$CVD_CAD)
  pv <- tryCatch(chisq.test(tab)$p.value, error = function(e) NA)
  levs <- rownames(tab)
  for (lv in levs) {
    n_yes <- tab[lv, "Yes"]; n_no <- tab[lv, "No"]
    pct_yes <- round(n_yes / sum(tab[, "Yes"]) * 100, 1)
    pct_no  <- round(n_no / sum(tab[, "No"]) * 100, 1)
    baseline <- rbind(baseline, data.frame(
      Variable = paste0(v, ": ", lv),
      Overall = paste0(n_yes + n_no, " (", round((n_yes+n_no)/nrow(orig_data)*100, 1), "%)"),
      CAD_Yes = paste0(n_yes, " (", pct_yes, "%)"),
      CAD_No  = paste0(n_no, " (", pct_no, "%)"),
      p_value = pv, stringsAsFactors = FALSE))
  }
}
write.csv(baseline, file.path(OUT_DIR, "baseline_characteristics.csv"), row.names = FALSE)
cat("\n===== Baseline characteristics (first 10 rows) =====\n"); print(head(baseline, 10))

# ---- 3.2 Survey-weighted univariate logistic regression (OR per 1 SD) ----
design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
                    data = std_data, nest = TRUE)
outcome <- "CVD_CAD"
excluded_vars <- c("SEQN","TyG_quartile","WTMEC2YR","SDMVSTRA","SDMVPSU")
covariates <- setdiff(colnames(std_data), c(outcome, excluded_vars))
cat("\nUnivariate candidate variables:", paste(covariates, collapse=", "), "\n")

uni <- data.frame(Variable = character(), OR = character(), p_value = numeric(),
                  Significant = character(), stringsAsFactors = FALSE)
for (var in covariates) {
  f <- as.formula(paste(outcome, "~", var))
  m <- tryCatch(svyglm(f, design = design, family = quasibinomial()), error = function(e) NULL)
  if (is.null(m)) { cat("  [skipped] ", var, "\n"); next }
  co <- summary(m)$coefficients
  for (i in 2:nrow(co)) {
    est <- co[i, 1]; se <- co[i, 2]; pv <- co[i, 4]
    or <- exp(est); lo <- exp(est - 1.96*se); hi <- exp(est + 1.96*se)
    nm <- if (nrow(co) > 2) paste0(var, " (", rownames(co)[i], ")") else var
    uni <- rbind(uni, data.frame(
      Variable = nm,
      OR = paste0(round(or, 2), " (", round(lo, 2), "-", round(hi, 2), ")"),
      p_value = pv,
      Significant = ifelse(pv < 0.05, "Significant", "Not significant"),
      stringsAsFactors = FALSE))
  }
}
uni$p_value <- round(uni$p_value, 6)
write.csv(uni, file.path(OUT_DIR, "Table1_univariate.csv"), row.names = FALSE)
cat("\n===== Weighted univariate results (all) =====\n")
print(uni, row.names = FALSE)

cat("\n[OK] baseline_characteristics.csv, Table1_univariate.csv written\n")


# =============================================================
# 4. Normality tests (Lilliefors and Anderson-Darling)
#    Input : cleaned_CAD_data_original_scale.csv
#    Output: normality_test_results.csv
#            normality_check.png (300 dpi) / .pdf
# =============================================================
pkgs <- c("nortest","ggplot2","patchwork","dplyr","tidyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
}

d <- read.csv(file.path(OUT_DIR, "cleaned_CAD_data_original_scale.csv"), stringsAsFactors = FALSE)

cont_vars <- c("RIDAGEYR","BMXBMI","BPXSY1","BPXDI1","LBXTR","LBXGLU",
               "LBDHDD","LBXTC","LBDLDL","TyG")
cont_vars <- cont_vars[cont_vars %in% names(d)]

res <- data.frame(Variable = character(), n = integer(),
                  Lilliefors_D = numeric(), Lilliefors_p = character(),
                  AD_A = numeric(), AD_p = character(),
                  Conclusion = character(), stringsAsFactors = FALSE)
for (v in cont_vars) {
  x <- d[[v]][!is.na(d[[v]])]
  ll <- lillie.test(x)
  ad <- ad.test(x)
  concl <- if (ll$p.value < 0.001 && ad$p.value < 0.001) "Reject normality (p<0.001)" else "Other"
  res <- rbind(res, data.frame(
    Variable = v, n = length(x),
    Lilliefors_D = round(ll$statistic, 4),
    Lilliefors_p = format.pval(ll$p.value, digits = 3),
    AD_A = round(ad$statistic, 3),
    AD_p = format.pval(ad$p.value, digits = 3),
    Conclusion = concl))
}
write.csv(res, file.path(OUT_DIR, "normality_test_results.csv"), row.names = FALSE)
cat("===== Normality test results =====\n")
print(res, row.names = FALSE)

# Histograms with fitted normal density
plot_var <- function(v) {
  x <- d[[v]][!is.na(d[[v]])]
  df <- data.frame(x = x)
  m <- mean(x); s <- sd(x)
  ggplot(df, aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30,
                   fill = "#7BC0CD", color = "white", alpha = 0.8) +
    stat_function(fun = dnorm, args = list(mean = m, sd = s),
                  color = "#ED8D5A", linewidth = 1) +
    labs(title = v, x = "", y = "") +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}
plist <- lapply(cont_vars, plot_var)
p_combined <- wrap_plots(plist, ncol = 2) +
  plot_annotation(title = "Normality Check (Lilliefors & Anderson-Darling: all p < 0.001)",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
ggsave(file.path(OUT_DIR, "normality_check.png"), p_combined, width = 12, height = 16, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "normality_check.pdf"), p_combined, width = 12, height = 16, bg = "white")
cat("\n[OK] normality_test_results.csv, normality_check.png/.pdf written\n")


# =============================================================
# 5. Machine-learning models for the three feature configurations
#    (Full / TyG-Ablated / TyG-Only)
#    Input : cleaned_CAD_data_fixed.csv
#    Pipeline: 80/20 split (seed = 1010) -> SMOTE (seed = 42) ->
#              five base learners (RF/XGB/GLM/SVM/KNN, 10-fold CV) ->
#              three ensembles (simple average / AUC-weighted / stacking)
#    Output per configuration <tag> in {full, notyg, onlytyg}:
#      internal_perf_<tag>.csv
#      internal_ROC_<tag>.png/.pdf, internal_PR_<tag>.png/.pdf,
#      internal_DCA_<tag>.png/.pdf, calibration_<tag>.png/.pdf
#      model_rds_<tag>/
# =============================================================
pkgs <- c("tidyverse","caret","pROC","doParallel","DMwR","PRROC","ggplot2","patchwork","glmnet")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
}
set.seed(42)
cl <- makeCluster(max(1, detectCores()-1)); registerDoParallel(cl)

run_configuration <- function(config_tag, config_label, features) {

  cat("\n=============================================\n")
  cat("Configuration:", config_label, "| n features:", length(features), "\n")
  cat("=============================================\n")

  # --- data and split ---
  cleaned_data <- read.csv(file.path(OUT_DIR, "cleaned_CAD_data_fixed.csv"), stringsAsFactors = FALSE)
  cleaned_data$CVD_CAD <- factor(cleaned_data$CVD_CAD, levels = c("No","Yes"))
  cat("Outcome distribution:\n"); print(table(cleaned_data$CVD_CAD))

  set.seed(1010)
  trainIndex <- createDataPartition(cleaned_data$CVD_CAD, p = 0.8, list = FALSE, times = 1)
  trainData <- cleaned_data[trainIndex, ]; testData <- cleaned_data[-trainIndex, ]

  train_X <- trainData %>% select(-CVD_CAD); train_Y <- trainData$CVD_CAD
  test_X  <- testData %>% select(-CVD_CAD);  test_Y  <- testData$CVD_CAD

  train_X <- train_X %>% select(all_of(features))
  test_X  <- test_X  %>% select(all_of(features))

  # Explicit numeric coding of categorical variables, identical across all
  # configurations and the validation step (Male = 1, Female = 2;
  # Smoker = 1, Non-smoker = 2). Numeric coding is also required so that
  # SMOTE can operate on all columns.
  encode_cat <- function(d) d %>% mutate(
    RIAGENDR = ifelse(RIAGENDR == "Male", 1, 2),
    SMQ020   = ifelse(SMQ020 == "Smoker", 1, 2))
  train_X <- encode_cat(train_X); test_X <- encode_cat(test_X)
  cat("Train/test distribution:\n"); print(table(train_Y)); print(table(test_Y))

  # --- SMOTE oversampling ---
  set.seed(42)
  balanced_data <- SMOTE(CVD_CAD ~ ., data = cbind(train_X, CVD_CAD = train_Y))
  train_X <- balanced_data[, -ncol(balanced_data)]
  train_Y <- balanced_data$CVD_CAD
  cat("Train distribution after SMOTE:\n"); print(table(train_Y))

  # KNN-specific preprocessing
  knn_preProc <- preProcess(train_X, method = c("center","scale","nzv"))
  train_X_knn <- predict(knn_preProc, train_X) %>% mutate(across(where(is.factor), as.numeric))
  test_X_knn  <- predict(knn_preProc, test_X)  %>% mutate(across(where(is.factor), as.numeric))

  # --- cross-validation settings and base learners ---
  ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE,
                       summaryFunction = twoClassSummary, savePredictions = "final",
                       allowParallel = TRUE, verboseIter = FALSE)
  models <- list(
    RF  = list(method="rf",  tuneGrid=data.frame(mtry=c(2, round(sqrt(ncol(train_X))))), preProcess=c("center","scale")),
    XGB = list(method="xgbTree",
               tuneGrid=expand.grid(nrounds=100, max_depth=c(4,6), eta=c(0.1,0.3), gamma=0,
                                    colsample_bytree=0.8, min_child_weight=1, subsample=0.8),
               preProcess=c("center","scale")),
    SVM = list(method="svmRadial", tuneLength=5, preProcess=c("center","scale")),
    GLM = list(method="glm", family=binomial(), preProcess=c("center","scale")),
    KNN = list(method="knn", tuneGrid=expand.grid(k=seq(5,25,by=5)), preProcess=NULL)
  )

  # --- train base learners ---
  model_list <- list()
  for (nm in names(models)) {
    cat("\n--- Training:", nm, "---\n")
    tr_X <- if (nm == "KNN") train_X_knn else train_X
    model_list[[nm]] <- tryCatch(
      train(x = tr_X, y = train_Y, method = models[[nm]]$method,
            trControl = ctrl, tuneGrid = models[[nm]]$tuneGrid,
            preProcess = models[[nm]]$preProcess, metric = "ROC"),
      error = function(e) { message("Model ", nm, " failed: ", e$message); NULL })
  }
  valid_models <- model_list[!sapply(model_list, is.null)]
  cat("\nBase models trained:", names(valid_models), "\n")

  # --- test-set evaluation and ensembles ---
  test_perf <- data.frame()
  ensemble_preds <- lapply(names(valid_models), function(nm) {
    pred_data <- if (nm == "KNN") test_X_knn else test_X
    predict(valid_models[[nm]], pred_data, type = "prob")[, "Yes"]
  })
  names(ensemble_preds) <- names(valid_models)
  pred_df <- as.data.frame(ensemble_preds)

  for (nm in names(valid_models)) {
    roc_obj <- roc(test_Y, pred_df[[nm]], levels = c("No","Yes"), direction = "<")
    test_perf <- rbind(test_perf, data.frame(
      Model = nm, AUC = round(auc(roc_obj), 3),
      CI = paste(round(ci(roc_obj), 3), collapse = "-")))
  }

  # Simple average
  pred_df$Average <- rowMeans(pred_df)
  avg_roc <- roc(test_Y, pred_df$Average, levels = c("No","Yes"), direction = "<")
  avg_auc <- round(auc(avg_roc), 3)

  # AUC-weighted average
  model_aucs <- sapply(names(valid_models), function(x) auc(roc(test_Y, pred_df[[x]])))
  weights <- model_aucs / sum(model_aucs)
  pred_df$Weighted <- as.numeric(rowSums(t(t(as.matrix(pred_df[, names(valid_models)])) * weights)))
  weighted_roc <- roc(test_Y, pred_df$Weighted, levels = c("No","Yes"), direction = "<")
  weighted_auc <- round(auc(weighted_roc), 3)
  weight_summary <- data.frame(Model = names(weights), Weight = round(weights, 3),
                               AUC = round(model_aucs, 3))
  cat("\n=== Model weights ===\n"); print(weight_summary)

  # Stacking (glmnet meta-learner)
  prepare_stacking_data <- function(models, data) {
    sapply(names(models), function(nm) {
      if (models[[nm]]$method == "knn") {
        d <- predict(knn_preProc, data); pred <- predict(models[[nm]], d, type = "prob")
      } else pred <- predict(models[[nm]], data, type = "prob")
      positive_col <- grep("X1|Yes", colnames(pred), value = TRUE)[1]
      if (is.na(positive_col)) positive_col <- 2
      pred[, positive_col]
    })
  }
  train_stacking <- prepare_stacking_data(valid_models, train_X)
  stack_data <- data.frame(train_stacking, Outcome = factor(train_Y, levels = c("No","Yes")))
  stack_data <- na.omit(stack_data)
  stack_ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                             summaryFunction = twoClassSummary, allowParallel = FALSE)
  lambda_seq <- 10^seq(-3, 3, length = 50)
  set.seed(42)
  stack_model <- train(Outcome ~ ., data = stack_data, method = "glmnet", family = "binomial",
                       tuneGrid = expand.grid(alpha = seq(0, 1, 0.1), lambda = lambda_seq),
                       metric = "ROC", trControl = stack_ctrl)
  test_stacking <- prepare_stacking_data(valid_models, test_X)
  stack_pred <- predict(stack_model, newdata = data.frame(test_stacking), type = "prob")[, "Yes"]
  stack_roc <- roc(test_Y, stack_pred, levels = c("No","Yes"), direction = "<")
  stack_auc <- round(auc(stack_roc), 3)
  cat("\nStacked ensemble test AUC:", stack_auc, "\n")

  # --- performance table ---
  final_perf <- rbind(test_perf, data.frame(
    Model = c("Average_Ensemble","Weighted_Ensemble","Stack_Ensemble"),
    AUC = c(avg_auc, weighted_auc, stack_auc), CI = NA))
  final_perf <- final_perf %>% arrange(desc(AUC))
  cat("\n===== Internal validation (20% test set), ranked =====\n"); print(final_perf)
  write.csv(final_perf, file.path(OUT_DIR, paste0("internal_perf_", config_tag, ".csv")), row.names = FALSE)

  # --- save models ---
  dir.create(file.path(OUT_DIR, paste0("model_rds_", config_tag)), showWarnings = FALSE)
  for (nm in names(valid_models))
    saveRDS(valid_models[[nm]], file.path(OUT_DIR, paste0("model_rds_", config_tag), paste0(nm, ".rds")))
  saveRDS(stack_model, file.path(OUT_DIR, paste0("model_rds_", config_tag), "stacked_ensemble.rds"))
  saveRDS(pred_df$Average,  file.path(OUT_DIR, paste0("model_rds_", config_tag), "simple_average.rds"))
  saveRDS(pred_df$Weighted, file.path(OUT_DIR, paste0("model_rds_", config_tag), "weighted_average.rds"))
  saveRDS(pred_df,          file.path(OUT_DIR, paste0("model_rds_", config_tag), "ensemble_preds.rds"))
  cat("\nModels saved to model_rds_", config_tag, "/\n", sep = "")

  # --- figures (PNG 300 dpi + PDF vector) ---
  save_fig <- function(p, name, w = 10, h = 6) {
    ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
    ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  }
  model_colors <- c(RF="#BFDFD2", XGB="#51999F", SVM="#4198AC", GLM="#7BC0CD", KNN="#DBCB92",
                    Average_Ensemble="#ECB66C", Weighted_Ensemble="#EA9E58", Stack_Ensemble="#ED8D5A")

  all_rocs <- c(lapply(names(valid_models), function(nm) roc(test_Y, pred_df[[nm]], levels=c("No","Yes"), direction="<")),
                list(Average_Ensemble=avg_roc, Weighted_Ensemble=weighted_roc, Stack_Ensemble=stack_roc))
  names(all_rocs) <- c(names(valid_models), "Average_Ensemble","Weighted_Ensemble","Stack_Ensemble")

  # ROC
  roc_df_all <- do.call(rbind, lapply(names(all_rocs), function(nm) {
    r <- all_rocs[[nm]]
    data.frame(Sensitivity = r$sensitivities, FPR = 1 - r$specificities,
               Model = nm, AUC = round(auc(r), 3))
  }))
  p_roc <- ggplot(roc_df_all, aes(x = FPR, y = Sensitivity, color = Model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 0.8) +
    labs(title = paste0("ROC Curves - ", config_label, " (Internal Validation)"),
         x = "1 - Specificity", y = "Sensitivity", color = "Model (AUC)") +
    scale_color_manual(values = model_colors,
                       labels = paste0(unique(roc_df_all$Model), " (AUC=", unique(roc_df_all$AUC), ")")) +
    theme_minimal(base_size = 12) + theme(legend.position = "right")
  print(p_roc); save_fig(p_roc, paste0("internal_ROC_", config_tag), 10, 6)

  # PR
  y_num <- as.numeric(test_Y) - 1
  pr_list <- lapply(names(all_rocs), function(nm) {
    pr.curve(scores.class0 = all_rocs[[nm]]$predictor, weights.class0 = y_num, curve = TRUE)
  })
  names(pr_list) <- names(all_rocs)
  pr_df_all <- do.call(rbind, lapply(names(pr_list), function(nm) {
    data.frame(Recall = pr_list[[nm]]$curve[,1], Precision = pr_list[[nm]]$curve[,2],
               Model = nm, AUCPR = round(pr_list[[nm]]$auc.integral, 3))
  }))
  p_pr <- ggplot(pr_df_all, aes(x = Recall, y = Precision, color = Model)) +
    geom_line(linewidth = 0.8) +
    labs(title = paste0("PR Curves - ", config_label, " (Internal Validation)"),
         x = "Recall", y = "Precision", color = "Model (AUC-PR)") +
    scale_color_manual(values = model_colors,
                       labels = paste0(unique(pr_df_all$Model), " (AUC-PR=", unique(pr_df_all$AUCPR), ")")) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    theme_minimal(base_size = 12) + theme(legend.position = "right")
  print(p_pr); save_fig(p_pr, paste0("internal_PR_", config_tag), 10, 6)

  # DCA (net benefit computed manually, avoids external dependencies)
  outcome01 <- as.numeric(test_Y == "Yes")
  calc_nb <- function(prob, t) {
    pred <- ifelse(prob >= t, 1, 0)
    tp <- sum(outcome01 == 1 & pred == 1); fp <- sum(outcome01 == 0 & pred == 1)
    tp/length(outcome01) - fp/length(outcome01) * (t/(1-t))
  }
  thr <- seq(0.01, 0.4, by = 0.01)
  dca_list <- lapply(names(all_rocs), function(nm) {
    data.frame(threshold = thr, NB = sapply(thr, function(t) calc_nb(all_rocs[[nm]]$predictor, t)),
               Model = nm)
  })
  dca_df <- do.call(rbind, dca_list)
  prevalence <- mean(outcome01)
  dca_ref <- data.frame(threshold = thr,
                        TreatAll = prevalence - (1-prevalence)*thr/(1-thr),
                        TreatNone = 0)
  p_dca <- ggplot(dca_df, aes(x = threshold, y = NB, color = Model)) +
    geom_line(linewidth = 0.8) +
    geom_line(data = dca_ref, aes(x = threshold, y = TreatAll), linetype = "dashed", color = "gray50") +
    geom_line(data = dca_ref, aes(x = threshold, y = TreatNone), linetype = "dotted", color = "black") +
    labs(title = paste0("Decision Curve Analysis - ", config_label, " (Internal)"),
         x = "Threshold Probability", y = "Net Benefit", color = "Model") +
    coord_cartesian(ylim = c(-0.05, 0.35)) +
    theme_minimal(base_size = 12) + theme(legend.position = "right")
  print(p_dca); save_fig(p_dca, paste0("internal_DCA_", config_tag), 10, 6)

  # Calibration (test set, equal-frequency deciles; rank-based binning avoids
  # non-unique quantile breaks caused by tied RF probabilities)
  cal_data <- data.frame(Outcome = outcome01,
                         lapply(names(all_rocs), function(nm) all_rocs[[nm]]$predictor))
  names(cal_data)[-1] <- names(all_rocs)
  cal_long <- cal_data %>% pivot_longer(-Outcome, names_to = "Model", values_to = "Pred")
  cal_bin <- cal_long %>% group_by(Model) %>%
    mutate(Bin = cut(rank(Pred, ties.method = "first"), breaks = 10, labels = FALSE)) %>%
    group_by(Model, Bin) %>%
    summarise(MeanPred = mean(Pred), MeanObs = mean(Outcome), N = n(), .groups = "drop") %>%
    filter(N >= 10)
  p_cal <- ggplot(cal_bin, aes(x = MeanPred, y = MeanObs, color = Model)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_point(size = 2, alpha = 0.8) + geom_line(linewidth = 0.5) +
    labs(title = paste0("Calibration - ", config_label, " (Internal)"),
         x = "Predicted Probability", y = "Observed Rate") +
    coord_equal(xlim = c(0,1), ylim = c(0,1)) +
    theme_minimal(base_size = 12) + theme(legend.position = "right")
  print(p_cal); save_fig(p_cal, paste0("calibration_", config_tag), 9, 6)

  cat("\n[OK]", config_label, "completed: internal_perf_", config_tag, ".csv and 4 figure sets (PNG+PDF).\n", sep = "")
}

# Full configuration (all nine variables)
run_configuration("full", "Full Model",
                  c("RIDAGEYR","RIAGENDR","LBXTC","LBDLDL","TyG","LBXGLU","BPXSY1","BPXDI1","SMQ020"))
# TyG-Ablated configuration (TyG removed, its components retained)
run_configuration("notyg", "No-TyG Model",
                  c("RIDAGEYR","RIAGENDR","LBXTC","LBDLDL","LBXGLU","BPXSY1","BPXDI1","SMQ020"))
# TyG-Only configuration (TyG retained, its components removed)
run_configuration("onlytyg", "TyG-Only Model",
                  c("RIDAGEYR","RIAGENDR","LBDLDL","TyG","BPXSY1","BPXDI1","SMQ020"))


# =============================================================
# 6. Temporally distinct validation (TyG-Only random forest)
#    Input : cleaned_CAD_externalvalidation_data_fixed.csv
#            model_rds_onlytyg/RF.rds  (produced by Section 5)
#    Note  : validation cohort is drawn from the same survey programme
#            as the derivation cohort, so this assesses temporal
#            generalizability rather than transportability.
#    Output: external_validation_main_metrics.csv
#            external_validation_confusion_matrix.csv
#            external_validation_metrics_by_threshold.csv
#            validation_ROC/PR/DCA/calibration .png + .pdf
# =============================================================
pkgs <- c("dplyr","pROC","caret","PRROC","ggplot2","ResourceSelection")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
}

external_data <- read.csv(file.path(OUT_DIR, "cleaned_CAD_externalvalidation_data_fixed.csv"), stringsAsFactors = FALSE)
external_data$CVD_CAD <- factor(external_data$CVD_CAD, levels = c("No","Yes"))
cat("Validation outcome distribution:\n"); print(table(external_data$CVD_CAD))

RF_model <- readRDS(file.path(OUT_DIR, "model_rds_onlytyg", "RF.rds"))

columns_to_keep <- c("RIDAGEYR","RIAGENDR","LBDLDL","TyG","BPXSY1","BPXDI1","SMQ020")
external_X <- external_data %>% select(all_of(columns_to_keep))
external_Y <- external_data$CVD_CAD

# Categorical coding aligned with training (Male = 1, Female = 2; Smoker = 1)
external_X <- external_X %>%
  mutate(
    RIAGENDR = as.numeric(factor(RIAGENDR, levels = c("Male","Female"))),
    SMQ020   = as.numeric(factor(SMQ020,   levels = c("Smoker","Non-smoker")))
  )
if (anyNA(external_X)) { keep <- complete.cases(external_X); external_X <- external_X[keep, ]; external_Y <- external_Y[keep] }

pred <- predict(RF_model, external_X, type = "prob")
positive_col <- grep("Yes", colnames(pred), value = TRUE)[1]
pred_prob <- pred[, positive_col]

# --- discrimination ---
roc_obj <- roc(external_Y, pred_prob, levels = c("No","Yes"), direction = "<")
auc_val <- as.numeric(auc(roc_obj))
cat("AUC:", round(auc_val, 4), " 95%CI:", paste(round(ci(roc_obj), 4), collapse="-"), "\n")

# --- threshold-based metrics ---
pred_class <- factor(ifelse(pred_prob > 0.5, "Yes", "No"), levels = c("No","Yes"))
cm <- confusionMatrix(pred_class, external_Y, positive = "Yes")
acc  <- cm$overall["Accuracy"]; sens <- cm$byClass["Sensitivity"]
spec <- cm$byClass["Specificity"]; ppv <- cm$byClass["Precision"]
npv  <- cm$byClass["Neg Pred Value"]; f1 <- cm$byClass["F1"]
bal_acc <- cm$byClass["Balanced Accuracy"]; kappa <- cm$overall["Kappa"]

# Youden-optimal threshold
youden <- roc_obj$sensitivities + roc_obj$specificities - 1
best_idx <- which.max(youden)
best_thr <- roc_obj$thresholds[best_idx]
best_sens <- roc_obj$sensitivities[best_idx]; best_spec <- roc_obj$specificities[best_idx]

# --- calibration metrics ---
outcome01 <- as.numeric(external_Y == "Yes")
brier <- mean((outcome01 - pred_prob)^2)
hl <- hoslem.test(outcome01, pred_prob, g = 10)
pred_c <- pmin(pmax(pred_prob, 0.001), 0.999)
logit_p <- qlogis(pred_c)
cal_fit <- glm(outcome01 ~ logit_p, family = binomial)
cal_intercept <- unname(coef(cal_fit)[1]); cal_slope <- unname(coef(cal_fit)[2])
cal_fit_off <- glm(outcome01 ~ offset(logit_p), family = binomial)
cal_int_off <- unname(coef(cal_fit_off)[1])

cat("\n===== Temporally distinct validation (TyG-Only RF) =====\n")
cat("AUC:", round(auc_val, 4), "\n")
cat("Accuracy:", round(acc, 4), " Sensitivity:", round(sens, 4), " Specificity:", round(spec, 4), "\n")
cat("PPV:", round(ppv, 4), " NPV:", round(npv, 4), " F1:", round(f1, 4), "\n")
cat("Balanced Accuracy:", round(bal_acc, 4), " Kappa:", round(kappa, 4), "\n")
cat("Brier:", round(brier, 4), "\n")
cat("Hosmer-Lemeshow: X2 =", round(hl$statistic, 2), "df =", hl$parameter, "p =", format(hl$p.value, digits=3), "\n")
cat("Calibration intercept (2-param):", round(cal_intercept, 3), " slope:", round(cal_slope, 3), "\n")
cat("Calibration intercept (offset-only):", round(cal_int_off, 3), "\n")
cat("Youden-optimal threshold:", round(best_thr, 3), " Sens:", round(best_sens, 4), " Spec:", round(best_spec, 4), "\n")
print(cm$table)

# --- tables ---
main_metrics <- data.frame(
  Metric = c("AUC","Accuracy","Sensitivity","Specificity","PPV","NPV","F1",
             "Balanced_Accuracy","Kappa","Brier_Score","HL_X2","HL_df","HL_p",
             "Calibration_intercept","Calibration_slope","Calibration_intercept_offset",
             "Best_Threshold","Best_Sensitivity","Best_Specificity","Prevalence"),
  Value = c(auc_val, acc, sens, spec, ppv, npv, f1, bal_acc, kappa, brier,
            as.numeric(hl$statistic), hl$parameter, hl$p.value,
            cal_intercept, cal_slope, cal_int_off, best_thr, best_sens, best_spec,
            mean(outcome01)))
write.csv(main_metrics, file.path(OUT_DIR, "external_validation_main_metrics.csv"), row.names = FALSE)

thresholds <- seq(0.1, 0.9, by = 0.1)
metrics_by_thr <- do.call(rbind, lapply(thresholds, function(t) {
  pc <- factor(ifelse(pred_prob > t, "Yes","No"), levels = c("No","Yes"))
  cmt <- confusionMatrix(pc, external_Y, positive = "Yes")
  data.frame(Threshold = t, Accuracy = cmt$overall["Accuracy"],
             Sensitivity = cmt$byClass["Sensitivity"], Specificity = cmt$byClass["Specificity"],
             PPV = cmt$byClass["Precision"], NPV = cmt$byClass["Neg Pred Value"],
             F1 = cmt$byClass["F1"])
}))
write.csv(metrics_by_thr, file.path(OUT_DIR, "external_validation_metrics_by_threshold.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(cm$table), file.path(OUT_DIR, "external_validation_confusion_matrix.csv"))

# --- figures ---
save_fig <- function(p, name, w = 8, h = 6) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}
roc_df <- data.frame(Sensitivity = roc_obj$sensitivities, FPR = 1 - roc_obj$specificities)
p_roc <- ggplot(roc_df, aes(FPR, Sensitivity)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_line(color = "#BFDFD2", linewidth = 1.2) +
  annotate("text", x = 0.6, y = 0.25, label = paste0("AUC = ", round(auc_val, 3)),
           color = "#4198AC", size = 5) +
  labs(title = "ROC - Temporally Distinct Validation (TyG-Only RF)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 12)
print(p_roc); save_fig(p_roc, "validation_ROC")

prc <- pr.curve(scores.class0 = pred_prob, weights.class0 = outcome01, curve = TRUE)
pr_df <- data.frame(Recall = prc$curve[,1], Precision = prc$curve[,2])
p_pr <- ggplot(pr_df, aes(Recall, Precision)) +
  geom_line(color = "#51999F", linewidth = 1.2) +
  annotate("text", x = 0.6, y = 0.3, label = paste0("AUC-PR = ", round(prc$auc.integral, 3)),
           color = "#4198AC", size = 5) +
  labs(title = "PR - Temporally Distinct Validation (TyG-Only RF)",
       x = "Recall", y = "Precision") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_minimal(base_size = 12)
print(p_pr); save_fig(p_pr, "validation_PR")

# DCA
thr <- seq(0.01, 0.4, by = 0.01)
nb <- sapply(thr, function(t) {
  pc <- ifelse(pred_prob >= t, 1, 0)
  tp <- sum(outcome01 == 1 & pc == 1); fp <- sum(outcome01 == 0 & pc == 1)
  tp/length(outcome01) - fp/length(outcome01) * (t/(1-t))
})
prevalence <- mean(outcome01)
dca_df <- data.frame(threshold = thr,
                     Model = nb,
                     TreatAll = prevalence - (1-prevalence)*thr/(1-thr),
                     TreatNone = 0)
dca_long <- data.frame(
  threshold = rep(thr, 3),
  Strategy  = rep(c("Model","TreatAll","TreatNone"), each = length(thr)),
  NB        = c(nb, dca_df$TreatAll, dca_df$TreatNone))
p_dca <- ggplot(dca_long, aes(threshold, NB, color = Strategy, linetype = Strategy)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(Model = "#ED8D5A", TreatAll = "gray50", TreatNone = "black")) +
  labs(title = "DCA - Temporally Distinct Validation (TyG-Only RF)",
       x = "Threshold Probability", y = "Net Benefit") +
  coord_cartesian(ylim = c(-0.05, 0.4)) +
  theme_minimal(base_size = 12)
print(p_dca); save_fig(p_dca, "validation_DCA")

# Calibration plot (equal-frequency deciles; rank-based binning)
grp <- cut(rank(pred_prob, ties.method = "first"), breaks = 10, labels = FALSE)
cal_tab <- data.frame(
  mean_pred = tapply(pred_prob, grp, mean),
  obs_rate  = tapply(outcome01, grp, mean))
p_cal <- ggplot(cal_tab, aes(mean_pred, obs_rate)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, color = "#1F4EFF") +
  geom_smooth(method = "loess", se = FALSE, color = "#B45309") +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = "Calibration - Temporally Distinct Validation (TyG-Only RF)",
       x = "Predicted Probability", y = "Observed Rate",
       caption = paste0("H-L p=", format(hl$p.value, digits=2),
                        "; intercept=", round(cal_int_off, 2), "; slope=", round(cal_slope, 2))) +
  theme_minimal(base_size = 12)
print(p_cal); save_fig(p_cal, "validation_calibration")

cat("\n[OK] Validation outputs written (tables + ROC/PR/DCA/calibration figures).\n")


# =============================================================
# 7. SHAP interpretability analysis (TyG-Only random forest)
#    Input : cleaned_CAD_data_fixed.csv, model_rds_onlytyg/RF.rds
#    Method: fastshap (nsim = 100, seed = 123), explaining the
#            SMOTE-balanced training set
#    Output: SHAP_importance.csv
#            SHAP_core_analysis.png/.pdf, SHAP_beeswarm.png/.pdf
# =============================================================
pkgs <- c("tidyverse","caret","DMwR","fastshap","ggplot2","patchwork","ggbeeswarm","reshape2","randomForest")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
}
set.seed(42)

cleaned_data <- read.csv(file.path(OUT_DIR, "cleaned_CAD_data_fixed.csv"), stringsAsFactors = FALSE)
cleaned_data$CVD_CAD <- factor(cleaned_data$CVD_CAD, levels = c("No","Yes"))

set.seed(1010)
trainIndex <- createDataPartition(cleaned_data$CVD_CAD, p = 0.8, list = FALSE, times = 1)
trainData <- cleaned_data[trainIndex, ]
train_X <- trainData %>% select(-CVD_CAD); train_Y <- trainData$CVD_CAD
columns_to_keep <- c("RIDAGEYR","RIAGENDR","LBDLDL","TyG","BPXSY1","BPXDI1","SMQ020")
train_X <- train_X %>% select(all_of(columns_to_keep))
# Categorical coding consistent with Sections 5 and 6
train_X <- train_X %>% mutate(
  RIAGENDR = ifelse(RIAGENDR == "Male", 1, 2),
  SMQ020   = ifelse(SMQ020 == "Smoker", 1, 2))

set.seed(42)
balanced_data <- SMOTE(CVD_CAD ~ ., data = cbind(train_X, CVD_CAD = train_Y))
train_X <- balanced_data[, -ncol(balanced_data)]
train_Y <- balanced_data$CVD_CAD
cat("SHAP input dimensions:", nrow(train_X), "x", ncol(train_X), "\n")

RF_model <- readRDS(file.path(OUT_DIR, "model_rds_onlytyg", "RF.rds"))
rf_final <- RF_model$finalModel
shap_data_matrix <- as.matrix(train_X)

# Explicit predict.randomForest wrapper (avoids S3 dispatch issues)
pfun <- function(object, newdata) predict(object, newdata, type = "prob")[, "Yes"]
set.seed(123)
shap_values <- fastshap::explain(rf_final, X = shap_data_matrix, nsim = 100, pred_wrapper = pfun, parallel = FALSE)

# --- feature-importance table ---
shap_importance <- data.frame(
  feature = colnames(shap_values),
  mean_abs_shap = colMeans(abs(shap_values))
) %>% arrange(desc(mean_abs_shap))
write.csv(shap_importance, file.path(OUT_DIR, "SHAP_importance.csv"), row.names = FALSE)
cat("\n===== SHAP feature importance (descending) =====\n"); print(shap_importance, row.names = FALSE)

my_colors <- c("#BFDFD2","#51999F","#4198AC","#7BC0CD","#DBCB92","#ECB66C","#EA9E58","#ED8D5A")

# --- importance bar chart ---
p_importance <- ggplot(shap_importance, aes(x = mean_abs_shap, y = reorder(feature, mean_abs_shap), fill = mean_abs_shap)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = round(mean_abs_shap, 4)), hjust = -0.1, size = 3.5) +
  scale_fill_gradientn(colors = my_colors[c(2,4,5,7)]) +
  labs(x = "Mean |SHAP value|", y = "", title = "Feature Importance (SHAP)") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 12) + theme(legend.position = "none", plot.title = element_text(face = "bold"))

# --- beeswarm summary ---
shap_long <- reshape2::melt(as.matrix(shap_values))
names(shap_long) <- c("row_id","feature","shap_value")
# melt expands columns first: all rows of one column, then the next column
row_idx <- rep(seq_len(nrow(shap_data_matrix)), ncol(shap_data_matrix))
col_idx <- rep(seq_len(ncol(shap_data_matrix)), each = nrow(shap_data_matrix))
shap_long$feature_value <- shap_data_matrix[cbind(row_idx, col_idx)]
p_beeswarm <- ggplot(shap_long, aes(x = shap_value, y = reorder(feature, abs(shap_value), FUN = median), color = feature_value)) +
  geom_quasirandom(alpha = 0.6, size = 1.2) +
  scale_color_gradient(low = my_colors[1], high = my_colors[5], name = "Feature value") +
  labs(x = "SHAP value", y = "", title = "SHAP Beeswarm Summary") +
  theme_minimal(base_size = 12) + theme(legend.position = "right")

# --- TyG dependence plot ---
tyg_data <- data.frame(shap_value = shap_values[, "TyG"], feature_value = shap_data_matrix[, "TyG"])
p_dependence_tyg <- ggplot(tyg_data, aes(x = feature_value, y = shap_value, color = feature_value)) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1) +
  scale_color_gradientn(colors = my_colors[c(1,3,5,7)], name = "TyG (z)") +
  labs(x = "TyG index (standardized)", y = "SHAP value", title = "SHAP Dependence - TyG Index") +
  theme_minimal(base_size = 12) + theme(legend.position = "right")

# --- TyG x age interaction plot ---
interaction_data <- data.frame(
  shap_value = shap_values[, "TyG"],
  tyg_value = shap_data_matrix[, "TyG"],
  age_value = shap_data_matrix[, "RIDAGEYR"])
interaction_data$age_group <- cut(interaction_data$age_value,
                                  breaks = quantile(interaction_data$age_value, c(0, 0.33, 0.66, 1), na.rm = TRUE),
                                  labels = c("Young","Middle","Old"), include.lowest = TRUE)
p_interaction <- ggplot(interaction_data, aes(x = tyg_value, y = shap_value, color = age_group)) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_color_manual(values = c(Young = my_colors[1], Middle = my_colors[5], Old = my_colors[7]), name = "Age group") +
  labs(x = "TyG index (standardized)", y = "SHAP value", title = "SHAP Interaction - TyG x Age") +
  theme_minimal(base_size = 12)

# --- combined figure ---
save_fig <- function(p, name, w = 14, h = 10) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}
combined_plot <- (p_importance + p_beeswarm) / (p_dependence_tyg + p_interaction) +
  plot_annotation(title = "SHAP Analysis - TyG-Only Random Forest",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))
print(combined_plot)
save_fig(combined_plot, "SHAP_core_analysis")

ggsave(file.path(OUT_DIR, "SHAP_beeswarm.png"), p_beeswarm, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "SHAP_beeswarm.pdf"), p_beeswarm, width = 10, height = 7, bg = "white")

cat("\n[OK] SHAP outputs written: SHAP_importance.csv, SHAP_core_analysis.png/.pdf, SHAP_beeswarm.png/.pdf\n")


# =============================================================
# 8. Derivation vs validation distribution comparison
#    Input : cleaned_CAD_data_original_scale.csv
#            cleaned_CAD_externalvalidation_data_original_scale.csv
#    Method: continuous variables -> Mann-Whitney U; categorical -> chi-square
#    Output: distribution_comparison_results.csv
#            distribution_boxplot.png/.pdf
# =============================================================
pkgs <- c("ggplot2","dplyr","tidyr","patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
}

train_data <- read.csv(file.path(OUT_DIR, "cleaned_CAD_data_original_scale.csv"), stringsAsFactors = FALSE)
validation_data <- read.csv(file.path(OUT_DIR, "cleaned_CAD_externalvalidation_data_original_scale.csv"), stringsAsFactors = FALSE)
train_data$Dataset <- "Training (2009-2018)"
validation_data$Dataset <- "Temporally distinct (2005-2008)"

continuous_vars <- c("RIDAGEYR","BMXBMI","BPXSY1","BPXDI1","LBXTR","LBXGLU",
                     "LBDHDD","LBXTC","LBDLDL","TyG")
results <- data.frame(Variable = character(), Train_Median = numeric(),
                      Validation_Median = numeric(), p_value = numeric(),
                      stringsAsFactors = FALSE)
for (v in continuous_vars) {
  if (v %in% names(train_data) && v %in% names(validation_data)) {
    pv <- wilcox.test(train_data[[v]], validation_data[[v]])$p.value
    results <- rbind(results, data.frame(
      Variable = v,
      Train_Median = round(median(train_data[[v]], na.rm = TRUE), 2),
      Validation_Median = round(median(validation_data[[v]], na.rm = TRUE), 2),
      p_value = pv))
  }
}
categorical_vars <- c("RIAGENDR","SMQ020","ALQ101")
for (v in categorical_vars) {
  if (v %in% names(train_data) && v %in% names(validation_data)) {
    tt <- table(train_data[[v]]); vt <- table(validation_data[[v]])
    lv <- union(names(tt), names(vt))
    m <- rbind(Train = tt[lv], Validation = vt[lv]); m[is.na(m)] <- 0
    pv <- tryCatch(chisq.test(m)$p.value, error = function(e) NA)
    results <- rbind(results, data.frame(
      Variable = v, Train_Median = NA, Validation_Median = NA, p_value = pv))
  }
}
results$Significant <- ifelse(!is.na(results$p_value) & results$p_value < 0.05, "Yes", "No")
results$p_value <- round(results$p_value, 6)
write.csv(results, file.path(OUT_DIR, "distribution_comparison_results.csv"), row.names = FALSE)
cat("\n===== Distribution comparison results =====\n"); print(results, row.names = FALSE)

# Side-by-side boxplots
common <- intersect(intersect(continuous_vars, names(train_data)), names(validation_data))
plot_list <- lapply(common, function(v) {
  d <- rbind(train_data[, c("Dataset", v)], validation_data[, c("Dataset", v)])
  ggplot(d, aes(x = Dataset, y = .data[[v]], fill = Dataset)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    labs(title = v, x = "", y = v) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 11))
})
combined <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(OUT_DIR, "distribution_boxplot.png"), combined, width = 12, height = 14, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "distribution_boxplot.pdf"), combined, width = 12, height = 14, bg = "white")
print(combined)

cat("\n[OK] distribution_comparison_results.csv, distribution_boxplot.png/.pdf written\n")


# =============================================================
# 9. Participant flow diagram (CONSORT-style)
#    Input : qc_report_train.csv, qc_report_valid.csv
#    Output: participant_flow.png (300 dpi) / .pdf / .svg (if svglite)
# =============================================================
pkgs <- c("ggplot2","dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
}

qc_t <- read.csv(file.path(OUT_DIR, "qc_report_train.csv"), stringsAsFactors = FALSE)
qc_v <- read.csv(file.path(OUT_DIR, "qc_report_valid.csv"), stringsAsFactors = FALSE)

make_flow <- function(qc, cohort_title) {
  short_stage <- c(
    "Raw merged data" = "Total NHANES participants",
    "Age>=18 & MCQ160C in 1,2" = "Adults with clear CAD status",
    "Remove SMQ020 = 9 and ALQ101 = 9" = "Remove unknown (9) smoking/alcohol",
    "Remove SMQ020 = 7 and 9" = "Remove refused (7)/unknown (9) smoking",
    "Missing values removed" = "Complete data on key variables",
    "Clinically implausible values removed" = "Clinically plausible values",
    "IQR outlier filtering (final)" = "After 1.5 x IQR filtering (final)"
  )
  short_reason <- c(
    "Age < 18 or no clear CAD diagnosis (MCQ160C not 1/2)" = "Age < 18 or no clear CAD diagnosis",
    "Smoking/alcohol status unknown (9)" = "Smoking/alcohol unknown (9)",
    "Smoking status refused (7) or unknown (9)" = "Smoking refused (7) or unknown (9)",
    "Missing values in any key variable" = "Missing values in key variables",
    "Values outside clinically plausible thresholds" = "Outside clinically plausible thresholds",
    "Values outside 1.5 x IQR range" = "Outside 1.5 x IQR range"
  )
  n <- nrow(qc)
  boxes <- data.frame(
    y = seq(1, n, 1), Stage = qc$Stage, N = qc$N, Excluded = qc$Excluded,
    Reason = qc$Reason, stringsAsFactors = FALSE)
  boxes$label <- paste0(short_stage[boxes$Stage], "\n(n = ", boxes$N, ")")
  boxes$excl_label <- ifelse(is.na(boxes$Excluded), "",
    paste0("Excluded ", boxes$Excluded, ":\n", short_reason[boxes$Reason]))
  boxes$cohort <- cohort_title
  boxes
}

flow_t <- make_flow(qc_t, "Derivation cohort: NHANES 2009-2018 (cycles F-J)")
flow_v <- make_flow(qc_v, "Temporally distinct validation cohort: NHANES 2005-2008 (cycles D-E)")

# Layout parameters
box_w <- 3.0; box_h <- 0.6
x_left <- 2.0
x_right <- 7.2
x_excl_t <- 5.2
x_excl_v <- 9.8
gap <- 1.15

flow_t$y_center <- seq(14.5 - box_h/2, by = -gap, length.out = nrow(flow_t))
flow_v$y_center <- seq(14.5 - box_h/2, by = -gap, length.out = nrow(flow_v))
y_span_v <- range(flow_v$y_center)
flow_v$y_center <- flow_v$y_center - (y_span_v[1] - 0.8)

p <- ggplot() + theme_void() +
  coord_cartesian(xlim = c(0, 10.8), ylim = c(0, 15.5))

p <- p + annotate("text", x = x_left, y = 15.3, label = "Derivation cohort: NHANES 2009-2018\n(cycles F-J)",
                  size = 3.6, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = x_right, y = 15.3, label = "Temporally distinct validation:\nNHANES 2005-2008 (cycles D-E)",
           size = 3.6, fontface = "bold", lineheight = 0.9)

for (i in 1:nrow(flow_t)) {
  yc <- flow_t$y_center[i]
  p <- p +
    geom_rect(aes(xmin = x_left - box_w/2, xmax = x_left + box_w/2,
                  ymin = yc - box_h/2, ymax = yc + box_h/2),
              fill = "white", color = "black", linewidth = 0.4) +
    annotate("text", x = x_left, y = yc, label = flow_t$label[i], size = 2.9, lineheight = 0.9)
  if (flow_t$excl_label[i] != "") {
    p <- p +
      geom_rect(aes(xmin = x_excl_t - 1.35, xmax = x_excl_t + 1.35,
                    ymin = yc - box_h/2, ymax = yc + box_h/2),
                fill = "grey92", color = "grey40", linewidth = 0.3) +
      annotate("text", x = x_excl_t, y = yc, label = flow_t$excl_label[i], size = 2.4, lineheight = 0.85) +
      geom_segment(aes(x = x_left + box_w/2, xend = x_excl_t - 1.35,
                       y = yc, yend = yc), linewidth = 0.3)
  }
  if (i < nrow(flow_t)) {
    p <- p + geom_segment(aes(x = x_left, xend = x_left,
                              y = yc - box_h/2, yend = flow_t$y_center[i+1] + box_h/2),
                          linewidth = 0.3)
  }
}
for (i in 1:nrow(flow_v)) {
  yc <- flow_v$y_center[i]
  p <- p +
    geom_rect(aes(xmin = x_right - box_w/2, xmax = x_right + box_w/2,
                  ymin = yc - box_h/2, ymax = yc + box_h/2),
              fill = "white", color = "black", linewidth = 0.4) +
    annotate("text", x = x_right, y = yc, label = flow_v$label[i], size = 2.9, lineheight = 0.9)
  if (flow_v$excl_label[i] != "") {
    p <- p +
      geom_rect(aes(xmin = x_excl_v - 1.35, xmax = x_excl_v + 1.35,
                    ymin = yc - box_h/2, ymax = yc + box_h/2),
                fill = "grey92", color = "grey40", linewidth = 0.3) +
      annotate("text", x = x_excl_v, y = yc, label = flow_v$excl_label[i], size = 2.4, lineheight = 0.85) +
      geom_segment(aes(x = x_right + box_w/2, xend = x_excl_v - 1.35,
                       y = yc, yend = yc), linewidth = 0.3)
  }
  if (i < nrow(flow_v)) {
    p <- p + geom_segment(aes(x = x_right, xend = x_right,
                              y = yc - box_h/2, yend = flow_v$y_center[i+1] + box_h/2),
                          linewidth = 0.3)
  }
}

p <- p +
  annotate("rect", xmin = x_left - box_w/2, xmax = x_left + box_w/2,
           ymin = 0.05, ymax = 0.35, fill = "#EAF2EE", color = "#1F4EFF", linewidth = 0.5) +
  annotate("text", x = x_left, y = 0.20, label = "Final training cohort",
           size = 2.6, fontface = "bold") +
  annotate("rect", xmin = x_right - box_w/2, xmax = x_right + box_w/2,
           ymin = 0.05, ymax = 0.35, fill = "#EAF2EE", color = "#1F4EFF", linewidth = 0.5) +
  annotate("text", x = x_right, y = 0.20, label = "Final validation cohort",
           size = 2.6, fontface = "bold")

print(p)

ggsave(file.path(OUT_DIR, "participant_flow.png"), p, width = 10, height = 12, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "participant_flow.pdf"), p, width = 10, height = 12, bg = "white")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(OUT_DIR, "participant_flow.svg"), p, width = 10, height = 12, bg = "white")
  cat("SVG (vector) version written\n")
} else {
  cat("(svglite not installed; SVG skipped, PDF is already vector)\n")
}

cat("\n[OK] participant_flow.png (300 dpi) / .pdf / .svg written\n")
cat("Counts are taken from qc_report_train.csv / qc_report_valid.csv.\n")


# =============================================================
# 10. Supplementary Table 1: missingness proportions
#     Input : XPT files in XPT_DIR_DERIVATION
#     Scope : candidate population = adults with clear CAD status
#             (MCQ160C in {1,2}), i.e. the population just before
#             complete-case exclusion in the derivation cohort.
#     Output: SuppTable1_missing_proportions.csv
# =============================================================
pkgs <- c("haven","dplyr","purrr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
}

read_nhanes_xpt <- function(file_prefix, cycles = c("F","G","H","I","J"), select_vars = NULL) {
  map_dfr(cycles, ~{
    data <- read_xpt(file.path(XPT_DIR_DERIVATION, paste0(file_prefix, "_", .x, ".XPT")))
    if (!is.null(select_vars)) data <- data %>% select(any_of(select_vars))
    data
  })
}
demo_data  <- read_nhanes_xpt("DEMO", select_vars = c("SEQN","RIDAGEYR","RIAGENDR","RIDRETH1","INDFMPIR"))
mcq_data   <- read_nhanes_xpt("MCQ",  select_vars = c("SEQN","MCQ160C"))
lipid_data <- list(
  hdl    = read_nhanes_xpt("HDL",    select_vars = c("SEQN","LBDHDD")),
  tchol  = read_nhanes_xpt("TCHOL",  select_vars = c("SEQN","LBXTC")),
  trigly = read_nhanes_xpt("TRIGLY", select_vars = c("SEQN","LBDLDL","LBXTR")),
  glu    = read_nhanes_xpt("GLU",    select_vars = c("SEQN","LBXGLU"))
) %>% reduce(full_join, by = "SEQN")
covar_data <- list(
  bmx = read_nhanes_xpt("BMX", select_vars = c("SEQN","BMXBMI")),
  bpx = read_nhanes_xpt("BPX", select_vars = c("SEQN","BPXSY1","BPXDI1")),
  smq = read_nhanes_xpt("SMQ", select_vars = c("SEQN","SMQ020")),
  alq = read_nhanes_xpt("ALQ", select_vars = c("SEQN","ALQ101"))
) %>% reduce(full_join, by = "SEQN")
merged_data <- list(demo_data, mcq_data, lipid_data, covar_data) %>% reduce(full_join, by = "SEQN")

step1_data <- merged_data %>% filter(RIDAGEYR >= 18, MCQ160C %in% c(1, 2))
cat("Candidate population (adults with clear CAD status) n =", nrow(step1_data), "\n")

vars <- c("RIDAGEYR","RIAGENDR","BMXBMI","BPXSY1","BPXDI1",
          "LBXGLU","LBXTR","SMQ020","ALQ101",
          "LBDLDL","LBXTC","LBDHDD","INDFMPIR","RIDRETH1")
var_labels <- c(RIDAGEYR="Age",
                RIAGENDR="Gender",
                BMXBMI="BMI",
                BPXSY1="Systolic blood pressure",
                BPXDI1="Diastolic blood pressure",
                LBXGLU="Fasting glucose",
                LBXTR="Triglycerides",
                SMQ020="Smoking status",
                ALQ101="Alcohol consumption",
                LBDLDL="LDL cholesterol",
                LBXTC="Total cholesterol",
                LBDHDD="HDL cholesterol",
                INDFMPIR="Family income (PIR)",
                RIDRETH1="Race/ethnicity")
miss_tab <- do.call(rbind, lapply(vars, function(v) {
  n_miss <- sum(is.na(step1_data[[v]]))
  data.frame(Variable = v,
             Label = var_labels[[v]],
             N_total = nrow(step1_data),
             N_missing = n_miss,
             Percent_missing = round(100 * n_miss / nrow(step1_data), 2))
}))
miss_tab <- miss_tab %>% arrange(desc(Percent_missing))
print(miss_tab)
write.csv(miss_tab, file.path(OUT_DIR, "SuppTable1_missing_proportions.csv"), row.names = FALSE)
cat("\n[OK] SuppTable1_missing_proportions.csv written\n")


# =============================================================
# 11. Supplementary Table 2: AUC-PR with bootstrap 95% CI
#     Input : cleaned_CAD_data_fixed.csv and model_rds_*/ objects
#             (produced by Section 5)
#     Method: reproduce the seed-1010 80/20 split to recover the exact
#             internal validation labels; read stored prediction
#             probabilities; compute AUC-PR with PRROC::pr.curve;
#             bootstrap B = 1000 row resamples for the 2.5/97.5 percentiles.
#     Output: SuppTable2_AUCPR_bootstrapCI.csv
# =============================================================
pkgs <- c("dplyr","caret","pROC","PRROC")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
}
B <- 1000

cleaned_data <- read.csv(file.path(OUT_DIR, "cleaned_CAD_data_fixed.csv"), stringsAsFactors = FALSE)
cleaned_data$CVD_CAD <- factor(cleaned_data$CVD_CAD, levels = c("No","Yes"))
set.seed(1010)
trainIndex <- createDataPartition(cleaned_data$CVD_CAD, p = 0.8, list = FALSE, times = 1)
test_Y  <- cleaned_data[-trainIndex, ]$CVD_CAD
y_num   <- as.numeric(test_Y) - 1          # No = 0, Yes = 1
cat("Internal validation set n =", length(y_num), " positives =", sum(y_num), "\n")

groups <- c(full = "Full", notyg = "TyG-Ablated", onlytyg = "TyG-Only")
results <- list()

for (g in names(groups)) {
  dir_g <- file.path(OUT_DIR, paste0("model_rds_", g))
  pred_df <- readRDS(file.path(dir_g, "ensemble_preds.rds"))
  stack_model <- readRDS(file.path(dir_g, "stacked_ensemble.rds"))

  base_cols <- intersect(c("RF","XGB","SVM","GLM","KNN"), names(pred_df))
  stack_pred <- predict(stack_model,
                        newdata = data.frame(pred_df[, base_cols, drop = FALSE]),
                        type = "prob")[, "Yes"]

  probs <- cbind(pred_df[, c(base_cols, "Average", "Weighted"), drop = FALSE],
                 Stack_Ensemble = stack_pred)
  cat("\n===", groups[[g]], "group, model columns:", paste(colnames(probs), collapse = ", "), "===\n")

  for (nm in colnames(probs)) {
    ap <- pr.curve(scores.class0 = probs[[nm]], weights.class0 = y_num, curve = FALSE)$auc.integral
    set.seed(2026)                                   # fixed seed: same bootstrap draws for all models
    bs <- replicate(B, {
      idx <- sample.int(length(y_num), replace = TRUE)
      pr.curve(scores.class0 = probs[[nm]][idx], weights.class0 = y_num[idx], curve = FALSE)$auc.integral
    })
    ci <- quantile(bs, probs = c(0.025, 0.975))
    cat(sprintf("  %-16s AUC-PR = %.3f  95%%CI: %.3f-%.3f\n", nm, ap, ci[1], ci[2]))
    results[[paste(g, nm)]] <- data.frame(
      Group = groups[[g]], Model = nm,
      AUC_PR = round(ap, 3),
      CI_lower = round(unname(ci[1]), 3),
      CI_upper = round(unname(ci[2]), 3))
  }
}

out <- do.call(rbind, results); rownames(out) <- NULL
cat("\n===== Summary table =====\n"); print(out)
write.csv(out, file.path(OUT_DIR, "SuppTable2_AUCPR_bootstrapCI.csv"), row.names = FALSE)
cat("\n[OK] SuppTable2_AUCPR_bootstrapCI.csv written\n")


# =============================================================
# 12. Supplementary Table 4: Hosmer-Lemeshow calibration table
#     Input : cleaned_CAD_externalvalidation_data_fixed.csv
#             model_rds_onlytyg/RF.rds
#     Output: SuppTable4_hl_groups.csv   (decile groups)
#             SuppTable4_hl_summary.csv  (chi-square / df / p / intercept / slope)
# =============================================================
pkgs <- c("dplyr","pROC","ResourceSelection")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
}

external_data <- read.csv(file.path(OUT_DIR, "cleaned_CAD_externalvalidation_data_fixed.csv"), stringsAsFactors = FALSE)
external_data$CVD_CAD <- factor(external_data$CVD_CAD, levels = c("No","Yes"))
RF_model <- readRDS(file.path(OUT_DIR, "model_rds_onlytyg", "RF.rds"))

columns_to_keep <- c("RIDAGEYR","RIAGENDR","LBDLDL","TyG","BPXSY1","BPXDI1","SMQ020")
external_X <- external_data %>% select(all_of(columns_to_keep))
external_Y <- external_data$CVD_CAD
external_X <- external_X %>%
  mutate(
    RIAGENDR = as.numeric(factor(RIAGENDR, levels = c("Male","Female"))),
    SMQ020   = as.numeric(factor(SMQ020,   levels = c("Smoker","Non-smoker")))
  )
if (anyNA(external_X)) { keep <- complete.cases(external_X); external_X <- external_X[keep, ]; external_Y <- external_Y[keep] }

pred <- predict(RF_model, external_X, type = "prob")
pred_prob <- pred[, grep("Yes", colnames(pred), value = TRUE)[1]]
outcome01 <- as.numeric(external_Y == "Yes")
cat("Validation set n =", length(outcome01), " positives =", sum(outcome01), "\n")

hl <- hoslem.test(outcome01, pred_prob, g = 10)
# Column names differ across package versions; use column 2 (the positive class)
cat("H-L observed column names:", paste(dimnames(hl$observed)[[2]], collapse = "|"), "\n")
cat("H-L expected column names:", paste(dimnames(hl$expected)[[2]], collapse = "|"), "\n")
pred_c <- pmin(pmax(pred_prob, 0.001), 0.999)
logit_p <- qlogis(pred_c)
cal_fit    <- glm(outcome01 ~ logit_p, family = binomial)
cal_fit_off<- glm(outcome01 ~ offset(logit_p), family = binomial)

grp_tab <- data.frame(
  Bin            = 1:nrow(hl$observed),
  N              = rowSums(hl$observed),
  Observed_Yes   = hl$observed[, 2],
  Expected_Yes   = round(hl$expected[, 2], 2),
  Observed_Rate  = round(hl$observed[, 2] / rowSums(hl$observed), 4),
  Expected_Rate  = round(hl$expected[, 2] / rowSums(hl$observed), 4)
)
print(grp_tab)
write.csv(grp_tab, file.path(OUT_DIR, "SuppTable4_hl_groups.csv"), row.names = FALSE)

summary_tab <- data.frame(
  Metric = c("Hosmer-Lemeshow X-squared","Degrees of freedom","p-value",
             "Calibration intercept (2-parameter model)","Calibration slope",
             "Calibration intercept (offset-only)"),
  Value = c(round(as.numeric(hl$statistic), 2),
            hl$parameter,
            format(hl$p.value, digits = 3, scientific = TRUE),
            round(unname(coef(cal_fit)[1]), 3),
            round(unname(coef(cal_fit)[2]), 3),
            round(unname(coef(cal_fit_off)[1]), 3))
)
print(summary_tab)
write.csv(summary_tab, file.path(OUT_DIR, "SuppTable4_hl_summary.csv"), row.names = FALSE)

cat("\n[OK] SuppTable4_hl_groups.csv and SuppTable4_hl_summary.csv written\n")


# =============================================================
# 13. Session information (software and package versions)
#     Output: sessionInfo.txt
# =============================================================
pkgs <- c("randomForest","xgboost","glmnet","e1071","caret",
          "pROC","PRROC","rmda","DMwR","ggplot2","dplyr",
          "purrr","haven","gridExtra","RColorBrewer")
loaded <- character()
for (p in pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    suppressMessages(library(p, character.only = TRUE))
    loaded <- c(loaded, p)
  }
}
cat("Packages loaded successfully:", paste(loaded, collapse = ", "), "\n")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("[OK] sessionInfo.txt written (R version and all loaded package versions)\n")

cat("\n===== Key package versions =====")
cat("\nR version:", R.version.string, "\n")
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
  cat(sprintf("%-15s: %s\n", p, ifelse(is.na(v), "NOT INSTALLED", v)))
}
