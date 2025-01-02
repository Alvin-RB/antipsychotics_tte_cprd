# -----------------------
# ANTIPSYCHOTICS TTE STUDY
# -----------------------
# Last run: 08/05/24

# EFFECTIVENESS OUTCOMES: MORTALITY AND DISCONTINUATION ####

# Clear memory
rm(list = ls())

# Packages
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(gtsummary)
library(gt)
library(forcats)
library(lubridate)
library(tidyr)
library(stringr)
library(reshape2)
library(webshot)
library(broom)
library(mice)
library(mitools)
library(miceadds)
library(reshape2)
library(grid)
library(gridExtra)
library(survival)
library(survminer)
library(cmprsk)
library(tidycmprsk)
library(ggsurvfit)
library(forestplot)
library(rlang)
library(patchwork)
library(RColorBrewer)
library(parallel)
library(cowplot)
library(scales)
library(mitml)
library(janitor)

# Set file path
path <- anonymised

# Set working directory
setwd(paste0(path))

#Load files
load(file = "tte_prevalent_cohort.Rdata")
load(file = "tte_imputed_long.Rdata")
tte_imputed <- as.mids(tte_imputed_long)

#Define adjustment variables
adjustment_set <- c("age_atprevcohortentry", "gender", "ethnicity_cat_cprdhes", "pat_2019imd_quintile", "diag_prev", "region", "prevcohortentry_year", "gpconsults_last6m", "smoking_status_cat", 
                    "prioralcoholabuse", "priorsubstanceabuse", "priordyslipidaemia", "priordiabetes", "priorhypertension", "priorcerebrovasculardisease", "priormyocardialinfarction", "priorrenaldisease", "priorliverdisease", 
                    "apuse_prior2years", "lipiddrugs_prior2years", "hypertensiondrugs_prior2years", "antidiabetics_prior2years", "antidepressant_prior2years", "moodstab_prior2years",
                    "baseline_totalcholesterol", "baseline_ldl", "baseline_hdl", "baseline_triglycerides", "baseline_systolicbp", "baseline_diastolicbp", "baseline_glucose", "baseline_hba1c", "baseline_weightkg", "baseline_bmi_cat")

# MORTALITY

# Formulas
mortality_itt_formula <- as.formula("Surv(time = mortality_fuptime_days_cen, event = diedby2y) ~ trt_group") # define formula
mortality_itt_adjusted_formula <- paste0(deparse1(mortality_itt_formula), " + ", paste(adjustment_set, collapse = " + "))

mortality_pp_formula <- as.formula("Surv(time = mortality_pp_fuptime_days_cen, event = diedby2yr_pp) ~ trt_group") # define formula
mortality_pp_adjusted_formula <- paste0(deparse1(mortality_pp_formula), " + ", paste(adjustment_set, collapse = " + "))

# Models

# ITT
mortality_itt_cox_unadjusted <-  coxph(mortality_itt_formula, data = tte_prevalent_cohort) # Unadjusted cox regression
mortality_itt_cox_adjusted <- with(tte_imputed, coxph(as.formula(mortality_itt_adjusted_formula))) # adjusted, full MI dataset

# PP
mortality_pp_cox_unadjusted <-  coxph(mortality_pp_formula,  data = tte_prevalent_cohort) # Unadjusted cox regression
mortality_pp_cox_adjusted <- with(tte_imputed, coxph(as.formula(mortality_pp_adjusted_formula))) # adjusted, full MI dataset

# TIME TO DISCONTINUATION ####

# Formulas
discont_formula <- as.formula("Surv(time = discontinuation_fuptime_days_cen, event = discontinued) ~ trt_group") # define formula
discont_adjusted_formula <- paste0(deparse1(discont_formula), " + ", paste(adjustment_set, collapse = " + "))

# Models

discont_cox_unadjusted <-  coxph(discont_formula, data = tte_prevalent_cohort)
discont_cox_adjusted <- with(tte_imputed, coxph(as.formula(discont_adjusted_formula))) # adjusted, full MI dataset

# COMPETING RISKS REGRESSION ####
# Primary competing risks analyses codes people who died prior to discontinuation as 2 for the event variables and uses the mortality fup time for them.
# This ensures those people are appropriately treated as having a competing risk of death (precluding hospitalisation), whereas the others remain the same (i.e. discontinued without dying, and had neither events)
discont_crr_unadjusted_formula <- as.formula("Surv(time = discdeath_crr_fuptime, event = discdeath_crr) ~ trt_group") # define formula
discont_crr_adjusted_formula <- paste0(deparse1(discont_crr_unadjusted_formula), " + ", paste(adjustment_set, collapse = " + "))

# Cumulative incidence
discont_cuminc <- cuminc(discont_crr_unadjusted_formula, data = tte_prevalent_cohort) # used for table
discont_crr_model <- survfit2(discont_crr_unadjusted_formula, data = tte_prevalent_cohort) # used for plot

# competing risks regression, unadjusted
discont_crr_unadjusted <- crr(discont_crr_unadjusted_formula, data = tte_prevalent_cohort) # unadjusted crr regression

# Competing risks regression, adjusted, using MI data
# the crr function does not seem compatible with MICE functions, therefore analysis does separately only each imputed dataset and combined at end

# Subset the dataframe by ".imp" variable and create separate data frames (exclude raw data)
imp_values <- unique(tte_imputed_long$.imp[tte_imputed_long$.imp != 0])

# Create separate data frames for each value except 0 (raw data)
for (value in imp_values) {
  assign(paste0("discont_data_", value), subset(tte_imputed_long, .imp == value))}

mc.cores <- detectCores() - 1 # define number of cores to use

# Run the models over each data frame in parallel (using mclapply) and store the results in a list
discont_crr_adjusted_results_list <- mclapply(imp_values, function(value) {
  crr(as.formula(discont_crr_adjusted_formula), data = get(paste0("discont_data_", value)))
}, mc.cores = mc.cores)

# Remove all the individual data frames that are no longer needed
rm(list = ls(pattern = "^discont_data_\\d+$"))

# Remove all formulae
rm(list = ls(pattern = "formula"))

# Tables

# Generate regression results table for cox models
reg_table <- function(model) {
  tab <- as.data.frame(tbl_regression(model,
                                      include = trt_group,
                                      estimate_fun = purrr::partial(style_ratio, digits = 8),
                                      pvalue_fun = function(x) style_pvalue(x, digits = 3))) %>%
    clean_names() %>%
    separate(x95_percent_ci, into = c("ci_high", "ci_low"), 
             sep = ", ") %>%
    mutate_at(vars(2:5), as.numeric) %>% # convert to numeric
    mutate_at(vars(2:4), ~round(exp(. * -1), 2)) %>% # inverse the log hazards, exponentiate it, round to 2dp
    select(group = characteristic, hr = log_hr, ci_low, ci_high, p_value) %>%
    filter(row_number() > 2) %>%
    mutate(result = paste0(format(hr, nsmall = 2), # nsmall ensures 2 digits are used in all cells (otherwise 1.00 might show as just 1)
                           " (", paste0(format(ci_low, nsmall = 2), ", ", format(ci_high, nsmall = 2), ")")),
           result = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", result),
           label = case_when(group == "Olanzapine" ~ "Aripiprazole vs. Olanzapine",
                             group == "Quetiapine" ~ "Aripiprazole vs. Quetiapine",
                             group == "Risperidone" ~ "Aripiprazole vs. Risperidone",
                             TRUE ~ "Unknown"))
  
  return(tab)}

# Function to extract pooled crr results
extract_pooled_crr_results <- function(results_list) {
  pooled_results <- MIcombine(results_list, method = "rubin")
  pooled_coefs <- coef(pooled_results) # extract coefficient
  pooled_ci <- confint(pooled_results) # extract CIs
  pooled_vcov <- vcov(pooled_results) # extract variance
  pooled_se <- sqrt(diag(pooled_vcov)) # calculate standard error
  pooled_p_values <- 2 * (1 - pnorm(abs(pooled_coefs / pooled_se))) # Calculate p-values
  
  pooled_df <- data.frame( # add to df
    rowname = rownames(pooled_ci),
    hr = round(exp(pooled_coefs * -1), 2), # inverse the log hazards, exponentiate it, round to 2dp
    ci_high = round(exp(pooled_ci[, 1] * -1), 2),
    ci_low = round(exp(pooled_ci[, 2] * -1), 2),
    p_value = round(pooled_p_values, 3)) %>%
    mutate(p_value = ifelse(p_value < 0.001, "<0.001", format(p_value, nsmall = 3))) %>% # use "<0.001" for really small p values
    filter(grepl("trt_group", rowname)) %>% # keep only the main coefficients
    rename(trt_group = rowname) %>%
    mutate(trt_group = gsub("trt_group", "", trt_group),
           result = paste0(
             format(hr, nsmall = 2), # nsmall ensures 2 digits are used in all cells (otherwise 1.00 might show as just 1)
             " (", paste0(format(ci_low, nsmall = 2), ", ", format(ci_high, nsmall = 2), ")")),
           result = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", result),
           label = case_when(trt_group == "Olanzapine" ~ "Aripiprazole vs. Olanzapine",
                             trt_group == "Quetiapine" ~ "Aripiprazole vs. Quetiapine",
                             trt_group == "Risperidone" ~ "Aripiprazole vs. Risperidone",
                             TRUE ~ "Unknown"))
  return(pooled_df)}

# Mortality

# Unadjusted
mortality_itt_cox_unadjusted_tab <- reg_table(mortality_itt_cox_unadjusted)
mortality_pp_cox_unadjusted_tab <- reg_table(mortality_pp_cox_unadjusted)

# Adjusted
mortality_adjustedresults_itt <- reg_table(mortality_itt_cox_adjusted)
mortality_adjustedresults_pp <- reg_table(mortality_pp_cox_adjusted)

# Discontinuation

# Unadjusted
discont_unadjusted_cox_tab <- reg_table(discont_cox_unadjusted)
discont_unadjustedmodel_crr_tab <- reg_table(discont_crr_unadjusted)

# Adjusted
discont_adjustedresults_cox_tab <- reg_table(discont_cox_adjusted)
discont_adjustedresults_crr_tab <- extract_pooled_crr_results(discont_crr_adjusted_results_list)

mortality_table <- mortality_itt_cox_unadjusted_tab %>%
  select(label, cox_itt_unadjusted = result) %>%
  left_join(select(mortality_adjustedresults_itt, label, cox_itt_adjusted = result), by = "label") %>%
  left_join(select(mortality_pp_cox_unadjusted_tab, label, cox_pp_unadjusted = result), by = "label") %>%
  left_join(select(mortality_adjustedresults_pp, label, cox_pp_adjusted = result), by = "label") %>%
  mutate(Outcome = "Mortality") %>%
  rename(Comparison = label)
 
discontinuation_table <- discont_unadjusted_cox_tab %>%
  select(label, cox_itt_unadjusted = result) %>%
  left_join(select(discont_adjustedresults_cox_tab, label, cox_itt_adjusted = result), by = "label") %>%
  left_join(select(discont_unadjustedmodel_crr_tab, label, crr_itt_unadjusted = result), by = "label") %>%
  left_join(select(discont_adjustedresults_crr_tab, label, crr_itt_adjusted = result), by = "label") %>%
  mutate(Outcome = "Discontinuation") %>%
  rename(Comparison = label)

discontinuation_mortality_table <- discontinuation_table %>%
  bind_rows(mortality_table) %>%
  gt(
    groupname_col = "Outcome") %>%
  tab_style(
    style = cell_text(style = "italic"), # italicize group names
    locations = cells_row_groups()) %>%
  sub_missing(missing_text = "-") %>% # use - for empty cells
  tab_spanner(label = "Cause-specific hazards model",  
              columns = c(cox_itt_unadjusted, cox_itt_adjusted), id = "itt_cox") %>%
  tab_spanner(label = "Subdistribution hazards model",  # add spanner for ITT SHRs column
              columns = c(crr_itt_unadjusted, crr_itt_adjusted), id = "itt_crr") %>%
  tab_spanner(label = "Cause-specific hazards model",  
              columns = c(cox_pp_unadjusted, cox_pp_adjusted), id = "pp_cox") %>%
  tab_spanner(label = "Intention-to-treat", 
              columns = c(cox_itt_unadjusted, cox_itt_adjusted, crr_itt_unadjusted, crr_itt_adjusted), id = "itt") %>%
  tab_spanner(label = "Per-protocol",  
              columns = c(cox_pp_unadjusted, cox_pp_adjusted), id = "pp") %>%
  tab_style(style = cell_text(weight = "bold"), # format spanners bold
            locations = cells_column_spanners(spanners = c("itt", "pp", "itt_cox", "itt_crr", "pp_cox"))) %>%
  cols_label( # tidy column labels
    cox_itt_unadjusted = "Unadjusted",
    cox_pp_unadjusted = "Unadjusted",
    crr_itt_unadjusted = "Unadjusted",
    cox_itt_adjusted = "Adjusted",
      crr_itt_adjusted = "Adjusted",
      cox_pp_adjusted = "Adjusted") %>%
  cols_align(
    align = "center", # center allign columns
    columns = everything()) %>%
  tab_footnote("Estimates are hazard ratios (95% confidence interval) reported after 2y of follow-up.") %>%
  tab_footnote(footnote = "Death was accounted for as a potential competing risk using the Fine-Gray subdistribution hazards model. Estimates are subdistribution hazard ratios.",
               locations = cells_column_spanners(spanners = "itt_crr")) %>%
  tab_footnote(footnote = "For discontinatuon, death was treated as a censoring event. Estimates are cause-specific hazard ratios.",
               locations = cells_column_spanners(spanners = c("itt_cox", "pp_cox"))) %>%
  tab_footnote("Adjusted for pre-specified baseline covariates: age, sex, ethnicity, SMI diagnosis category, prior use of antipsychotics, level of deprivation (quintile), geographic region, 
               calendar year of index date, number of primary care consultations in prior six months, smoking status, comorbidities (dyslipidaemia, diabetes, hypertension, cerebrovascular disease, 
               myocardial infarction, renal failure, liver disease, alcohol misuse, substance misuse), concomitant medications (lipid-regulating medications, antihypertensives, antidiabetics, 
               antidepressants, mood stabilisers) and cardiometabolic values (total cholesterol, LDL-C, HDL-C, triglycerides, systolic blood pressure, diastolic blood pressure, glucose, HbA1c, weight, 
               BMI category).", locations = cells_column_labels(columns = c(cox_itt_adjusted, crr_itt_adjusted, cox_pp_adjusted))) %>%
  tab_footnote(footnote = "Patients were censored at discontinuation or switch to one of the other study medications.",
             locations = cells_column_spanners(spanners = "pp"))

# Package versions
packageVersion("dplyr") # 1.1.4
packageVersion("tibble") # 3.2.1
packageVersion("purrr") # 1.0.2
packageVersion("ggplot2") # 3.5.1
packageVersion("gtsummary") # 2.0.1
packageVersion("gt") # 0.11.0
packageVersion("forcats") # 1.0.0
packageVersion("lubridate") # 1.9.3
packageVersion("tidyr") # 1.3.1
packageVersion("stringr") # 1.5.1
packageVersion("reshape2") # 1.4.4
packageVersion("webshot") # 0.5.5
packageVersion("broom") # 1.0.6
packageVersion("mice") # 3.16.0
packageVersion("mitools") # 2.4
packageVersion("miceadds") # 3.17.44
packageVersion("grid") # 4.4.1
packageVersion("gridExtra") # 2.3
packageVersion("survival") # 3.7.0
packageVersion("survminer") # 0.4.9
packageVersion("cmprsk") # 2.2.12
packageVersion("tidycmprsk") # 1.0.0
packageVersion("ggsurvfit") # 1.1.0
packageVersion("forestplot") # 3.1.3
packageVersion("rlang") # 1.1.4
packageVersion("patchwork") # 1.3.0
packageVersion("RColorBrewer") # 1.1.3
packageVersion("parallel") # 4.4.1
packageVersion("cowplot") # 1.1.3
packageVersion("scales") # 1.3.0
packageVersion("mitml") # 0.4.5
packageVersion("janitor") # 2.2.0
