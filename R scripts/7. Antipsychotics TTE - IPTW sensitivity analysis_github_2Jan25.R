# -----------------------
# ANTIPSYCHOTICS TTE STUDY
# -----------------------

# IPTW SENSITIVITY ANALYSIS

# Conducted for total cholesterol at one year and for psychiatric hospitalisation
# Last run 05/11/2024

# Clear memory
rm(list = ls())

library(dplyr)
library(readr)
library(forcats)
library(tidyr)
library(survival)
library(survminer)
library(webshot)
library(mice)
library(mitools)
library(WeightIt)
library(lubridate)
library(cobalt)
library(gtsummary)
library(gt)
library(ggplot2)
library(grid)
library(gridExtra)
library(MatchThem)
library(janitor)
library(patchwork)

# Set file path
path <- anonymised 

# Set working directory
setwd(paste0(path))

#Load files
load(file = "tte_imputed_long.Rdata")
tte_imputed <- as.mids(tte_imputed_long)

adjustment_set <- c("age_atprevcohortentry", "gender", "ethnicity_cat_cprdhes", "pat_2019imd_quintile", "diag_prev", "region", "prevcohortentry_year", "gpconsults_last6m", "smoking_status_cat", 
                    "prioralcoholabuse", "priorsubstanceabuse", "priordyslipidaemia", "priordiabetes", "priorhypertension", "priorcerebrovasculardisease", "priormyocardialinfarction", "priorrenaldisease", "priorliverdisease", 
                    "apuse_prior2years", "lipiddrugs_prior2years", "hypertensiondrugs_prior2years", "antidiabetics_prior2years", "antidepressant_prior2years", "moodstab_prior2years",
                    "baseline_totalcholesterol", "baseline_ldl", "baseline_hdl", "baseline_triglycerides", "baseline_systolicbp", "baseline_diastolicbp", "baseline_glucose", "baseline_hba1c", "baseline_weightkg", "baseline_bmi_cat")

# Calculate inverse probability weights
# Use MatchThem package (which is designed to calculate weights from imputed data)
tte_imputed_weighted <- MatchThem::weightthem(trt_group ~ age_atprevcohortentry + gender + ethnicity_cat_cprdhes + pat_2019imd_quintile + diag_prev + region + prevcohortentry_year + gpconsults_last6m +
                                  smoking_status_cat + prioralcoholabuse + priorsubstanceabuse + priordyslipidaemia + priordiabetes + priorhypertension + 
                                  priorcerebrovasculardisease + priormyocardialinfarction + apuse_prior2years + lipiddrugs_prior2years + hypertensiondrugs_prior2years +
                                  antidiabetics_prior2years + antidepressant_prior2years + moodstab_prior2years + baseline_totalcholesterol + baseline_ldl + baseline_hdl + 
                                  baseline_triglycerides + baseline_systolicbp + baseline_diastolicbp + baseline_glucose + baseline_hba1c + baseline_weightkg + baseline_bmi_cat,
                                datasets = tte_imputed, # the mids object created by mice
                                approach = "within",  # use a within methods rather than across (within methods have been shown to be more consistent)
                                method = "ps", # default propensity score method which uses logistic regression
                                stabilize = TRUE, # stabilize weights to reduce the impact of extreme weights
                                estimand = "ATE") # specify estimand of interest

# Observe success of covariate balance
loveplot <- love.plot(tte_imputed_weighted,
          binary = "std",
          abs = TRUE,
          threshold = c(m = .1),
          colors = c("darkorange", "darkgreen"),
          shapes = c("square", "circle"),
          sample.names = c("Unweighted", "IPTW"),
          var.names = c(age_atprevcohortentry = "Age at cohort entry",
                        gender = "Sex",
                        ethnicity_cat_cprdhes = "Ethnicity",
                        diag_prev = "SMI diagnosis",
                        pat_2019imd_quintile = "Deprivation quintile",
                        region = "Region",
                        prevcohortentry_year = "Index year",
                        gpconsults_last6m = "Number of GP consultations",
                        smoking_status_cat = "Smoking status",
                        prioralcoholabuse = "Alcohol misuse",
                        priorsubstanceabuse = "Substance misuse",
                        priordyslipidaemia = "Dyslipidaemia",
                        priordiabetes = "Diabetes",
                        priorhypertension = "Hypertension",
                        priorrenaldisease = "Renal disease",
                        priorcerebrovasculardisease = "Cerebrovascular disease",
                        priormyocardialinfarction = "Myocardial infarction",
                        priorliverdisease = "Liver disease",
                        apuse_prior2years = "Prior antipsychotic use",
                        lipiddrugs_prior2years = "Lipid-regulating medications",
                        hypertensiondrugs_prior2years = "Antihypertensives",
                        antidiabetics_prior2years = "Antidiabetics",
                        antidepressant_prior2years = "Antidepressants",
                        moodstab_prior2years = "Mood stabilisers",
                        baseline_totalcholesterol = "Total cholesterol",
                        baseline_ldl = "LDL-C",
                        baseline_hdl = "HDL-C",
                        baseline_glucose = "Glucose",
                        baseline_hba1c = "HbA1c",
                        baseline_weightkg = "Body weight",
                        baseline_triglycerides = "Triglycerides",
                        baseline_systolicbp = "Systolic BP",
                        baseline_diastolicbp = "Diastolic BP",
                        baseline_bmi_cat = "BMI category"),
          limits = c(0, .65),
          position = c(.91, .10),
          which.treat = "Aripiprazole") +
  theme(legend.box.background = element_rect(), 
        legend.box.margin = margin(1, 1, 1, 1))

baltable <- bal.tab(tte_imputed_weighted,
                    stats = c("m", "v"), 
                    threshold = c(m = .1), 
                    which.treat = "Aripiprazole")

# Calculate mean weights for each patient-treatment group combination
mean <- MatchThem::complete(tte_imputed_weighted, action = "long", include = FALSE) %>%
  select(patid, trt_group, weights) %>%
  group_by(patid, trt_group) %>%
  summarise(weights = mean(weights), .groups = 'drop')

# Summary statistics table with Mean (SD) format, ensuring 2 decimal places
group_stats <- mean %>%
  group_by(trt_group) %>%
  summarise(
    `Mean (SD)` = sprintf("%.2f (%.2f)", mean(weights), sd(weights)),
    Min = sprintf("%.2f", min(weights)),
    Max = sprintf("%.2f", max(weights))
  ) %>%
  gt() %>%
  cols_label(trt_group = "Treatment Group") %>%
  tab_options(
    table.font.size = 14,
    column_labels.font.size = 15,
    data_row.padding = px(5)
  ) %>%
  tab_style(
    style = list(cell_text(align = "center")),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "grey80"),
    locations = cells_column_labels()
  )

# Density plot for weight distribution without grid lines
mean_weight_distribution <- ggplot(mean, aes(x = weights, fill = trt_group)) +
  geom_density(alpha = 0.5) +
  labs(x = "Mean inverse probability weight", y = "Density", fill = "Treatment group") +
  scale_x_continuous(limits = c(0, 4)) +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#CC79A7", "#F0E442")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    panel.grid = element_blank()
  )

# Combine plot and table with patchwork
mean_weight_distribution_combined <- mean_weight_distribution / patchwork::wrap_table(group_stats, panel = "cols", space = "free_y") +
  plot_layout(ncol = 1, heights = c(3, 1))

# ANALYSIS

# Total cholesterol at one year

tte_imputed_weighted_long <- MatchThem::complete(tte_imputed_weighted, action = "long")

# ITT: Limit to those alive at 1 year
tte_imputed_weighted_long_tc <- tte_imputed_weighted_long %>%
  filter(diedby1y == 0)

#  Per protocol: Limit to those alive and active at 1 year
tte_imputed_weighted_long_tc_pp <- tte_imputed_weighted_long %>%
  filter(active_1y == 0)

# Create an empty list to store the regression results
regression_results <- list()
regression_results_pp <- list()

# Iterate over unique .imp values
for (imp_val in unique(tte_imputed_weighted_long_tc$.imp)) {
  # Subset the data for the current imp value
  subset_data <- subset(tte_imputed_weighted_long_tc, .imp == imp_val)
  subset_data_pp <- subset(tte_imputed_weighted_long_tc_pp, .imp == imp_val)
  
  # linear regression with weights
  model <- lm(oney_totalcholesterol ~ trt_group, data = subset_data, weights = subset_data$weights)
  model_pp <- lm(oney_totalcholesterol ~ trt_group, data = subset_data_pp, weights = subset_data_pp$weights)
  
  # Store the regression results in the list
  regression_results[[as.character(imp_val)]] <- summary(model)
  regression_results_pp[[as.character(imp_val)]] <- summary(model_pp)}

# Combine the regression results using pool
combined_results_itt <- MIcombine(regression_results)
combined_results_pp <- MIcombine(regression_results_pp)


# Function to process combined regression results
process_combined_results <- function(combined_results, result_type) {
  # Extract coefficients
  coefficients <- combined_results$coefficients
  
  # Extract estimates and standard errors
  estimates <- coefficients[, "Estimate"]
  std_errors <- coefficients[, "Std. Error"]
  
  # Calculate 95% confidence intervals
  z_value <- qnorm(0.975)  # For 95% CI
  lower_ci <- estimates - z_value * std_errors
  upper_ci <- estimates + z_value * std_errors
  
  # Combine results into a data frame
  pooled_df <- data.frame(
    term = rownames(coefficients),  # Coerce row labels into a column
    Estimate = estimates,
    Std.Error = std_errors,
    CI.Lower = lower_ci,
    CI.Upper = upper_ci
  )
  
  # Add a 'type' column based on the result_type argument
  pooled_df <- pooled_df %>%
    mutate(type = ifelse(grepl("itt", result_type, ignore.case = TRUE), "ITT", 
                         ifelse(grepl("pp", result_type, ignore.case = TRUE), "PP", NA)))
  
  return(pooled_df)}

pooled_df_itt <- process_combined_results(combined_results_itt, "itt")
pooled_df_pp <- process_combined_results(combined_results_pp, "pp")

tc_weighted <- pooled_df_itt %>% 
  clean_names() %>%
  select(term, estimate, ci_lower, ci_upper) %>% 
  rename(itt_estimate = estimate, itt_ci_lower = ci_lower, itt_ci_upper = ci_upper) %>%  # Rename estimate and CI columns for ITT
  left_join(
    pooled_df_pp %>% 
      clean_names() %>%
      select(term, estimate, ci_lower, ci_upper) %>% 
      rename(pp_estimate = estimate, pp_ci_lower = ci_lower, pp_ci_upper = ci_upper),  # Rename estimate and CI columns for PP
    by = "term"
  ) %>%
  filter(grepl("trt_group", term)) %>% # Filter to main coefficients of interest
  mutate(term = gsub("trt_group", "", term)) %>% # Remove trt_group variable name
  mutate(
    # Create result strings for ITT and PP with CIs, padded to 2 decimal places
    result_itt = paste0(format(round(itt_estimate * -1, 2), nsmall = 2), " (", 
                        format(round(itt_ci_lower * -1, 2), nsmall = 2), ", ", 
                        format(round(itt_ci_upper * -1, 2), nsmall = 2), ")"),
    result_pp = paste0(format(round(pp_estimate * -1, 2), nsmall = 2), " (", 
                       format(round(pp_ci_lower * -1, 2), nsmall = 2), ", ", 
                       format(round(pp_ci_upper * -1, 2), nsmall = 2), ")")
  ) %>%
  mutate(
    label = case_when(
      term == "Olanzapine" ~ "Aripiprazole vs. Olanzapine",
      term == "Quetiapine" ~ "Aripiprazole vs. Quetiapine",
      term == "Risperidone" ~ "Aripiprazole vs. Risperidone",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(label, result_itt, result_pp)  # Select relevant columns

# Create the gt table
tc_weighted_gt <- tc_weighted %>%
  gt() %>%
  tab_header(title = "Pooled IPTW estimates for ITT and PP") %>%
  cols_label(
    label = "Comparison",
    result_itt = "ITT",
    result_pp = "PP")

# Psychiatric hospitalisation

cox_model <- with(tte_imputed_weighted, coxph(Surv(psychhosp_main_fuptime_days_cen, hospadmit_2yr_any) ~ trt_group))
cox_model <- tbl_regression(cox_model, exp = FALSE,
               pvalue_fun = function(x) style_pvalue(x, digits = 3))

cox_model_tab <- as.data.frame(cox_model) %>%
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
                           TRUE ~ "Unknown")) %>%
  select(label, result, p_value) %>%
  mutate(type = "IPTW ITT")

cox_model_pp <- with(tte_imputed_weighted, coxph(Surv(psychhosp_pp_fuptime_days_cen, hospadmit_2yr_pp) ~ trt_group))
cox_model_pp <- tbl_regression(cox_model_pp, exp = FALSE,
               pvalue_fun = function(x) style_pvalue(x, digits = 3))

cox_model_pp_tab <- as.data.frame(cox_model_pp) %>%
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
                           TRUE ~ "Unknown")) %>%
  select(label, result, p_value) %>%
  mutate(type = "IPTW PP")

psychhosp_weighted <- cox_model_tab %>%
  bind_rows(cox_model_pp_tab) %>%
  gt(groupname_col = "type")

# IPTW for psychiatric hospitalisation, which includes prior hospitalisation in the calculation of the weights

#Load files
load(file = "tte_imputed_long.Rdata")

tte_imputed_priorhosp <- tte_imputed_long %>%
  filter(!is.na(priorhosp_any)) %>%
  as.mids()

adjustment_set <- c("age_atprevcohortentry", "gender", "ethnicity_cat_cprdhes", "pat_2019imd_quintile", "diag_prev", "region", "prevcohortentry_year", "gpconsults_last6m", "smoking_status_cat", 
                    "prioralcoholabuse", "priorsubstanceabuse", "priordyslipidaemia", "priordiabetes", "priorhypertension", "priorcerebrovasculardisease", "priormyocardialinfarction", "priorrenaldisease", "priorliverdisease", 
                    "apuse_prior2years", "lipiddrugs_prior2years", "hypertensiondrugs_prior2years", "antidiabetics_prior2years", "antidepressant_prior2years", "moodstab_prior2years",
                    "baseline_totalcholesterol", "baseline_ldl", "baseline_hdl", "baseline_triglycerides", "baseline_systolicbp", "baseline_diastolicbp", "baseline_glucose", "baseline_hba1c", "baseline_weightkg", "baseline_bmi_cat",
                    "priorhosp_any")

# Calculate inverse probability weights
# Use MatchThem package (which is designed to calculate weights from imputed data)
tte_imputed_weighted_priorhosp <- MatchThem::weightthem(trt_group ~ age_atprevcohortentry + gender + ethnicity_cat_cprdhes + pat_2019imd_quintile + diag_prev + region + prevcohortentry_year + gpconsults_last6m +
                                                smoking_status_cat + prioralcoholabuse + priorsubstanceabuse + priordyslipidaemia + priordiabetes + priorhypertension + 
                                                priorcerebrovasculardisease + priormyocardialinfarction + apuse_prior2years + lipiddrugs_prior2years + hypertensiondrugs_prior2years +
                                                antidiabetics_prior2years + antidepressant_prior2years + moodstab_prior2years + baseline_totalcholesterol + baseline_ldl + baseline_hdl + 
                                                baseline_triglycerides + baseline_systolicbp + baseline_diastolicbp + baseline_glucose + baseline_hba1c + baseline_weightkg + baseline_bmi_cat + priorhosp_any,
                                              datasets = tte_imputed_priorhosp, # the mids object created by mice
                                              approach = "within",  # use a within methods rather than across (within methods have been shown to be more consistent)
                                              method = "ps", # default propensity score method which uses logistic regression
                                              stabilize = TRUE, # stabilize weights to reduce the impact of extreme weights
                                              estimand = "ATE") # specify estimand of interest

# Observe success of covariate balance
loveplot <- love.plot(tte_imputed_weighted_priorhosp,
                      binary = "std",
                      abs = TRUE,
                      threshold = c(m = .1),
                      colors = c("darkorange", "darkgreen"),
                      shapes = c("square", "circle"),
                      sample.names = c("Unweighted", "IPTW"),
                      var.names = c(age_atprevcohortentry = "Age at cohort entry",
                                    gender = "Sex",
                                    ethnicity_cat_cprdhes = "Ethnicity",
                                    diag_prev = "SMI diagnosis",
                                    pat_2019imd_quintile = "Deprivation quintile",
                                    region = "Region",
                                    prevcohortentry_year = "Index year",
                                    gpconsults_last6m = "Number of GP consultations",
                                    smoking_status_cat = "Smoking status",
                                    prioralcoholabuse = "Alcohol misuse",
                                    priorsubstanceabuse = "Substance misuse",
                                    priordyslipidaemia = "Dyslipidaemia",
                                    priordiabetes = "Diabetes",
                                    priorhypertension = "Hypertension",
                                    priorrenaldisease = "Renal disease",
                                    priorcerebrovasculardisease = "Cerebrovascular disease",
                                    priormyocardialinfarction = "Myocardial infarction",
                                    priorliverdisease = "Liver disease",
                                    apuse_prior2years = "Prior antipsychotic use",
                                    lipiddrugs_prior2years = "Lipid-regulating medications",
                                    hypertensiondrugs_prior2years = "Antihypertensives",
                                    antidiabetics_prior2years = "Antidiabetics",
                                    antidepressant_prior2years = "Antidepressants",
                                    moodstab_prior2years = "Mood stabilisers",
                                    baseline_totalcholesterol = "Total cholesterol",
                                    baseline_ldl = "LDL-C",
                                    baseline_hdl = "HDL-C",
                                    baseline_glucose = "Glucose",
                                    baseline_hba1c = "HbA1c",
                                    baseline_weightkg = "Body weight",
                                    baseline_triglycerides = "Triglycerides",
                                    baseline_systolicbp = "Systolic BP",
                                    baseline_diastolicbp = "Diastolic BP",
                                    baseline_bmi_cat = "BMI category",
                                    priorhosp_any = "Prior hospitalisation"),
                      limits = c(0, .65),
                      position = c(.91, .10),
                      which.treat = "Aripiprazole") +
  theme(legend.box.background = element_rect(), 
        legend.box.margin = margin(1, 1, 1, 1))

baltable <- bal.tab(tte_imputed_weighted_priorhosp,
                    stats = c("m", "v"), 
                    threshold = c(m = .1), 
                    which.treat = "Aripiprazole")

# Calculate mean weights for each patient-treatment group combination
mean <- MatchThem::complete(tte_imputed_weighted_priorhosp, action = "long", include = FALSE) %>%
  select(patid, trt_group, weights) %>%
  group_by(patid, trt_group) %>%
  summarise(weights = mean(weights), .groups = 'drop')

# Summary statistics table with Mean (SD) format, ensuring 2 decimal places
group_stats <- mean %>%
  group_by(trt_group) %>%
  summarise(
    `Mean (SD)` = sprintf("%.2f (%.2f)", mean(weights), sd(weights)),
    Min = sprintf("%.2f", min(weights)),
    Max = sprintf("%.2f", max(weights))) %>%
  gt() %>%
  cols_label(trt_group = "Treatment Group") %>%
  tab_options(
    table.font.size = 14,
    column_labels.font.size = 15,
    data_row.padding = px(5)) %>%
  tab_style(
    style = list(cell_text(align = "center")),
    locations = cells_body(columns = everything())) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "grey80"),
    locations = cells_column_labels())

# Density plot for weight distribution without grid lines
mean_weight_distribution <- ggplot(mean, aes(x = weights, fill = trt_group)) +
  geom_density(alpha = 0.5) +
  labs(x = "Mean inverse probability weight", y = "Density", fill = "Treatment group") +
  scale_x_continuous(limits = c(0, 4)) +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#CC79A7", "#F0E442")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    panel.grid = element_blank())

# Combine plot and table with patchwork
mean_weight_distribution_combined <- mean_weight_distribution / patchwork::wrap_table(group_stats, panel = "cols", space = "free_y") +
  plot_layout(ncol = 1, heights = c(3, 1))

# ANALYSIS

# Psychiatric hospitalisation

cox_model_priorhosp <- with(tte_imputed_weighted_priorhosp, coxph(Surv(psychhosp_main_fuptime_days_cen, hospadmit_2yr_any) ~ trt_group))
cox_model_priorhosp <- tbl_regression(cox_model_priorhosp, exp = FALSE,
                            pvalue_fun = function(x) style_pvalue(x, digits = 3))

cox_model_tab_priorhosp <- as.data.frame(cox_model_priorhosp) %>%
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
                           TRUE ~ "Unknown")) %>%
  select(label, result, p_value) %>%
  mutate(type = "IPTW ITT")

cox_model_pp_priorhosp <- with(tte_imputed_weighted_priorhosp, coxph(Surv(psychhosp_pp_fuptime_days_cen, hospadmit_2yr_pp) ~ trt_group))
cox_model_pp_priorhosp <- tbl_regression(cox_model_pp_priorhosp, exp = FALSE,
                               pvalue_fun = function(x) style_pvalue(x, digits = 3))

cox_model_pp_tab_priorhosp <- as.data.frame(cox_model_pp_priorhosp) %>%
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
                           TRUE ~ "Unknown")) %>%
  select(label, result, p_value) %>%
  mutate(type = "IPTW PP")

psychhosp_weighted_priorhosp <- cox_model_tab_priorhosp %>%
  bind_rows(cox_model_pp_tab_priorhosp) %>%
  gt(groupname_col = "type") %>%
  tab_footnote("Variables in the weights were the pre-specified covariates plus prior psychiatric hospitalisation.")

# Package versions
packageVersion("dplyr") # 1.1.4
packageVersion("readr") # 2.1.5
packageVersion("forcats") # 1.0.0
packageVersion("tidyr") # 1.3.1
packageVersion("survival") # 3.7.0
packageVersion("survminer") # 0.4.9
packageVersion("webshot") # 0.5.5
packageVersion("mice") # 3.16.0
packageVersion("mitools") # 2.4
packageVersion("WeightIt") # 1.3.0
packageVersion("lubridate") # 1.9.3
packageVersion("cobalt") # 4.5.5
packageVersion("gtsummary") # 2.0.1
packageVersion("gt") # 0.11.0
packageVersion("ggplot2") # 3.5.1
packageVersion("grid") # 4.4.1
packageVersion("gridExtra") # 2.3
packageVersion("MatchThem") # 1.2.1
packageVersion("janitor") # 2.2.0
packageVersion("patchwork") # 1.3.0
