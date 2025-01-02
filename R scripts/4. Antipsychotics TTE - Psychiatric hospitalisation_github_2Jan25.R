# -----------------------
# ANTIPSYCHOTICS TTE STUDY
# -----------------------
# Last run: 01/11/24

# EFFECTIVENESS OUTCOME: PSYCHIATRIC HOSPITALISATION ####

# Clear memory
rm(list = ls())

# Packages
invisible(sapply(c("dplyr", "tibble", "purrr", "ggplot2", "gtsummary", "gt", "forcats", "lubridate", "tidyr",
                   "mice", "stringr", "reshape2", "webshot", "broom", "mitools", "miceadds", "grid", "gridExtra",
                   "survival", "survminer", "cmprsk", "tidycmprsk", "ggsurvfit", "forestplot", "rlang", "patchwork",
                   "RColorBrewer", "parallel", "cowplot", "scales", "mitml", "janitor", "tidylog"), 
                 library, character.only = TRUE))

# Set file path
path <- anonymised

# Set working directory
setwd(paste0(path))

# Load files
load(file = "tte_prevalent_cohort.Rdata")
load(file = "tte_imputed_long.Rdata")
tte_imputed <- as.mids(tte_imputed_long)

# Define adjustment variables
adjustment_set <- c("age_atprevcohortentry", "gender", "ethnicity_cat_cprdhes", "pat_2019imd_quintile", "diag_prev", "region", "prevcohortentry_year", "gpconsults_last6m", "smoking_status_cat", 
                    "prioralcoholabuse", "priorsubstanceabuse", "priordyslipidaemia", "priordiabetes", "priorhypertension", "priorcerebrovasculardisease", "priormyocardialinfarction", "priorrenaldisease", "priorliverdisease", 
                    "apuse_prior2years", "lipiddrugs_prior2years", "hypertensiondrugs_prior2years", "antidiabetics_prior2years", "antidepressant_prior2years", "moodstab_prior2years",
                    "baseline_totalcholesterol", "baseline_ldl", "baseline_hdl", "baseline_triglycerides", "baseline_systolicbp", "baseline_diastolicbp", "baseline_glucose", "baseline_hba1c", "baseline_weightkg", "baseline_bmi_cat")


# Model formulas

# ITT - cox regression (any hospitalisation record with a primary psychiatric ICD code)
psychhosp_itt_formula <- as.formula("Surv(time = psychhosp_main_fuptime_days_cen, event = hospadmit_2yr_any) ~ trt_group") # define ITT formula
psychhosp_itt_adjusted_formula <- paste0(deparse1(psychhosp_itt_formula), " + ", paste(adjustment_set, collapse = " + "))

# PP - cox regression (patient censored at switch or discontinuation)
psychhosp_pp_formula <- as.formula("Surv(time = psychhosp_pp_fuptime_days_cen, event = hospadmit_2yr_pp) ~ trt_group") # PP formula
psychhosp_pp_adjusted_formula <- paste0(deparse1(psychhosp_pp_formula), " + ", paste(adjustment_set, collapse = " + "))

# ITT - cox regression - Post-hoc sensitivity analysis: additional adjustment for prior hospitalisation
adjustment_set_priorhosp <- c(adjustment_set, "priorhosp_any")
psychhosp_itt_adjusted_formula_model2 <- paste0(deparse1(psychhosp_itt_formula), " + ", paste(adjustment_set_priorhosp, collapse = " + "))

# PP - cox regression - Post-hoc sensitivity analysis: additional adjustment for prior hospitalisation
psychhosp_pp_adjusted_formula_model2 <- paste0(deparse1(psychhosp_pp_formula), " + ", paste(adjustment_set_priorhosp, collapse = " + "))

# ITT - cox regression - sensitivity analysis: psychiatric hospitalisation definition change (any hospitalisation record with a primary psychiatric ICD code & corresponding psychiatric consultant codes)
psychhosp_definitionsens_formula <- as.formula("Surv(time = psychhosp_sens_fuptime_days_cen, event = hospadmit_psychspecialty_2yr_any) ~ trt_group") # define formula
psychhosp_definitionsens_adjusted_formula <- paste0(deparse1(psychhosp_definitionsens_formula), " + ", paste(adjustment_set, collapse = " + "))

# Competing risks - ITT
psychhosp_itt_crr_formula <- as.formula("Surv(time = hospdeath_crr_fuptime, event = hospdeath_crr) ~ trt_group") # ITT formula
psychhosp_itt_crr_adjusted_formula <- paste0(deparse1(psychhosp_itt_crr_formula), " + ", paste(adjustment_set, collapse = " + "))

# Competing risks - PP
psychhosp_pp_crr_formula <- as.formula("Surv(time = hospdeath_pp_crr_fuptime, event = hospdeath_pp_crr) ~ trt_group") # PP formula
psychhosp_pp_crr_adjusted_formula <- paste0(deparse1(psychhosp_pp_crr_formula), " + ", paste(adjustment_set, collapse = " + "))

# Competing risks - ITT - Post-hoc sensitivity analysis: additional adjustment for prior hospitalisation
psychhosp_itt_crr_adjusted_formula_model2 <- paste0(deparse1(psychhosp_itt_crr_formula), " + ", paste(adjustment_set_priorhosp, collapse = " + "))

# Competing risks - PP - Post-hoc sensitivity analysis: additional adjustment for prior hospitalisation
psychhosp_pp_crr_adjusted_formula_model2 <- paste0(deparse1(psychhosp_pp_crr_formula), " + ", paste(adjustment_set_priorhosp, collapse = " + "))

# Models

# Unadjusted
psychhosp_itt_cox_unadj <-  coxph(psychhosp_itt_formula, data = tte_prevalent_cohort) # ITT - Cox regression, unadjusted
psychhosp_pp_cox_unadj <-  coxph(psychhosp_pp_formula, data = tte_prevalent_cohort) # PP - Cox regression, unadjusted

psychhosp_crr <-  crr(psychhosp_itt_crr_formula, data = tte_prevalent_cohort) # ITT unadjusted crr regression
psychhosp_pp_crr <-  crr(psychhosp_pp_crr_formula, data = tte_prevalent_cohort) # PP unadjusted crr regression

# Adjusted
psychhosp_itt_cox_adjusted <- with(tte_imputed, coxph(as.formula(psychhosp_itt_adjusted_formula))) # ITT - adjusted, using MI dataset
psychhosp_itt_cox_adjusted_model2 <- with(tte_imputed, coxph(as.formula(psychhosp_itt_adjusted_formula_model2))) # ITT - adjusted, full MI dataset

psychhosp_pp_cox_adjusted <- with(tte_imputed, coxph(as.formula(psychhosp_pp_adjusted_formula))) # PP - adjusted, using MI dataset
psychhosp_pp_cox_adjusted_model2 <- with(tte_imputed, coxph(as.formula(psychhosp_pp_adjusted_formula_model2))) # PP - adjusted, full MI dataset

# Competing risks regression ####

# adjusted, using MI data
# the crr function does not seem compatible with MICE functions, therefore analysis does separately on each imputed dataset and combined at end

# Subset the dataframe by ".imp" variable and create separate data frames (exclude raw data)
imp_values <- unique(tte_imputed_long$.imp[tte_imputed_long$.imp != 0])

# Create separate data frames for each imp value except 0 (raw data)
for (value in imp_values) {
  assign(paste0("psychhosp_data_", value), subset(tte_imputed_long, .imp == value & hes_apc_e == 1))}

# Run the models over each data frame in parallel (using mclapply) and store the results in a list
mc.cores <- detectCores() - 1 # define number of cores to use

# ITT
psychhosp_crr_results_list <- mclapply(imp_values, function(value) {
  crr(as.formula(psychhosp_itt_crr_adjusted_formula), data = get(paste0("psychhosp_data_", value)))}, mc.cores = mc.cores)

# ITT - Post-hoc sensitivity analysis: additional adjustment for prior hospitalisation
psychhosp_crr_hospsens_results_list <- mclapply(imp_values, function(value) {
  crr(as.formula(psychhosp_itt_crr_adjusted_formula_model2), data = get(paste0("psychhosp_data_", value)))}, mc.cores = mc.cores)

# PP
psychhosp_crr_pp_results_list <- mclapply(imp_values, function(value) {
  crr(as.formula(psychhosp_pp_crr_adjusted_formula), data = get(paste0("psychhosp_data_", value)))}, mc.cores = mc.cores)

# PP - Post-hoc sensitivity analysis: additional adjustment for prior hospitalisation
psychhosp_crr_pp_hospsens_results_list <- mclapply(imp_values, function(value) {
  crr(as.formula(psychhosp_pp_crr_adjusted_formula_model2), data = get(paste0("psychhosp_data_", value)))}, mc.cores = mc.cores)

# Remove all the individual data frames
rm(list = ls(pattern = "^psychhosp_data_\\d+$"))

# TABLES ####

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

# Unadjusted tables
itt_unadjustedmodel <- reg_table(psychhosp_itt_cox_unadj)
itt_unadjustedmodel_crr <- reg_table(psychhosp_crr)

pp_unadjustedmodel <- reg_table(psychhosp_pp_cox_unadj)
pp_unadjustedmodel_crr <- reg_table(psychhosp_pp_crr)

# Adjusted tables
itt_model1 <- reg_table(psychhosp_itt_cox_adjusted)
itt_model2 <- reg_table(psychhosp_itt_cox_adjusted_model2)

psychhosp_itt_crr_model1_tab <- extract_pooled_crr_results(psychhosp_crr_results_list)
psychhosp_itt_crr_model2_tab <- extract_pooled_crr_results(psychhosp_crr_hospsens_results_list)

pp_model1 <- reg_table(psychhosp_pp_cox_adjusted)
pp_model2 <- reg_table(psychhosp_pp_cox_adjusted_model2)

psychhosp_pp_crr_model1_tab <- extract_pooled_crr_results(psychhosp_crr_pp_results_list)
psychhosp_pp_crr_model2_tab <- extract_pooled_crr_results(psychhosp_crr_pp_hospsens_results_list)

# Number of events / denominator
psychhosp_events <- tte_prevalent_cohort %>%
  filter(hes_apc_e == 1) %>%
  group_by(trt_group) %>%
  summarise(itt_events = paste(sum(hospadmit_2yr_any)),
            pp_events = paste(sum(hospadmit_2yr_pp)),
            n = n())

# Main results table
psych_hosp_main_results <- itt_model1 %>%
  left_join(itt_model2, by = "label") %>%
  select(label, model1 = result.x, model2 = result.y) %>%
  rename(Comparison = label) %>%
  gt() %>%
  tab_footnote(paste0("Analyses based on ", psychhosp_events[[1, 2]], " events among ", psychhosp_events[[1, 4]], " patients in the ", psychhosp_events[[1, 1]], " group, ",
                      psychhosp_events[[2, 2]], " events among ", psychhosp_events[[2, 4]], " patients in the ", psychhosp_events[[2, 1]], " group, ",
                      psychhosp_events[[3, 2]], " events among ", psychhosp_events[[3, 4]], " patients in the ", psychhosp_events[[3, 1]], " group, and ",
                      psychhosp_events[[4, 2]], " events among ", psychhosp_events[[4, 4]], " patients in the ", psychhosp_events[[4, 1]], " group.")) %>%
  tab_footnote("Estimates are cause-specific hazard ratios (95% confidence interval) reported after 2y of follow-up, where death was treated as a censoring event.") %>%
  tab_footnote("Estimates are cause-specific hazard ratios (95% confidence interval) reported after 2y of follow-up, where death was treated as a censoring event.") %>%
  tab_footnote("Adjusted for pre-specified baseline covariates: age, sex, ethnicity, SMI diagnosis category, prior use of antipsychotics, level of deprivation (quintile), geographic region, 
               calendar year of index date, number of primary care consultations in prior six months, smoking status, comorbidities (dyslipidaemia, diabetes, hypertension, cerebrovascular disease, 
               myocardial infarction, renal failure, liver disease, alcohol misuse, substance misuse), concomitant medications (lipid-regulating medications, antihypertensives, antidiabetics, 
               antidepressants, mood stabilisers) and cardiometabolic values (total cholesterol, LDL-C, HDL-C, triglycerides, systolic blood pressure, diastolic blood pressure, glucose, HbA1c, weight, 
               BMI category).", locations = cells_column_labels(columns = model1)) %>%
  tab_footnote("As above, but additionally adjusted for psychiatric hospitalization in the prior two years.", locations = cells_column_labels(columns = model2))

# Table of unadjusted and adjusted ITT and PP results
psych_hosp_results_itt <- itt_unadjustedmodel %>%
  select(label, unadjusted = result) %>%
  left_join(select(itt_model1, label, model1 = result), by = "label") %>%
  left_join(select(itt_model2, label, model2 = result), by = "label") %>%
  left_join(select(itt_unadjustedmodel_crr, label, unadjusted_crr = result), by = "label") %>%
  left_join(select(psychhosp_itt_crr_model1_tab, label, model1_crr = result), by = "label") %>%
  left_join(select(psychhosp_itt_crr_model2_tab, label, model2_crr = result), by = "label") %>%
  mutate(type = "Intention-to-treat") %>%
  rename(Comparison = label)

psych_hosp_results_pp <- pp_unadjustedmodel %>%
  select(label, unadjusted = result) %>%
  left_join(select(pp_model1, label, model1 = result), by = "label") %>%
  left_join(select(pp_model2, label, model2 = result), by = "label") %>%
  left_join(select(pp_unadjustedmodel_crr, label, unadjusted_crr = result), by = "label") %>%
  left_join(select(psychhosp_pp_crr_model1_tab, label, model1_crr = result), by = "label") %>%
  left_join(select(psychhosp_pp_crr_model2_tab, label, model2_crr = result), by = "label") %>%
  mutate(type = "Per-protocol")  %>%
  rename(Comparison = label)

psych_hosp_results <- psych_hosp_results_itt %>%
  bind_rows(psych_hosp_results_pp) %>%
  gt(
    groupname_col = "type") %>%
  tab_style(
    style = cell_text(style = "italic"), # italicize group names
    locations = cells_row_groups()) %>%
  sub_missing(missing_text = "-") %>% # use - for empty cells
  tab_spanner(label = "Cause-specific hazards model", # add spanner for ITT CSHRs column
              columns = c(unadjusted, model1, model2), id = "cox") %>%
  tab_spanner(label = "Subdistribution hazards model",  # add spanner for ITT SHRs column
              columns = c(unadjusted_crr, model1_crr, model2_crr), id = "crr") %>%
  tab_style(style = cell_text(weight = "bold"), # format spanners bold
            locations = cells_column_spanners(spanners = c("cox", "crr"))) %>%
  cols_label( # tidy column labels
    unadjusted = "Unadjusted",
    unadjusted_crr = "Unadjusted",
    model1 = "Model 1",
    model1_crr = "Model 1",
    model2 = "Model 2",
    model2_crr = "Model 2") %>%
  cols_align(
    align = "center", # center allign columns
    columns = everything()) %>%
  tab_footnote("Estimates are hazard ratios (95% confidence interval) reported after 2y of follow-up.") %>%
  tab_footnote(footnote = "Death was accounted for as a potential competing risk using the Fine-Gray subdistribution hazards model. Estimates are subdistribution hazard ratios.",
               locations = cells_column_spanners(spanners = "crr")) %>%
  tab_footnote(footnote = "Death was treated as a censoring event. Estimates are cause-specific hazard ratios.",
               locations = cells_column_spanners(spanners = "cox")) %>%
  tab_footnote("Adjusted for pre-specified baseline covariates: age, sex, ethnicity, SMI diagnosis category, prior use of antipsychotics, level of deprivation (quintile), geographic region, 
               calendar year of index date, number of primary care consultations in prior six months, smoking status, comorbidities (dyslipidaemia, diabetes, hypertension, cerebrovascular disease, 
               myocardial infarction, renal failure, liver disease, alcohol misuse, substance misuse), concomitant medications (lipid-regulating medications, antihypertensives, antidiabetics, 
               antidepressants, mood stabilisers) and cardiometabolic values (total cholesterol, LDL-C, HDL-C, triglycerides, systolic blood pressure, diastolic blood pressure, glucose, HbA1c, weight, 
               BMI category).", locations = cells_column_labels(columns = c(model1, model1_crr))) %>%
  tab_footnote("As above, but additionally adjusted for psychiatric hospitalization in the prior two years.", locations = cells_column_labels(columns = c(model2, model2_crr))) %>%
  tab_footnote(footnote = "Patients were censored at discontinuation or switch to one of the other study medications.",
             locations = cells_row_groups(2))

# Forest plots
setwd(paste0(path))

forestplot_combined <- function(data, xlab) {
  
  data %>%
    group_by(model) %>%
    forestplot::forestplot(fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                           labeltext = c(label),
                           clip = c(0.90, 1.15),
                           xlog = TRUE,
                           mean = hr,
                           lower = ci_low,
                           upper = ci_high,
                           vertices = TRUE,
                           grid = FALSE,
                           xlab = "Adjusted HR",
                           ci.vertices.height = .05,
                           boxsize = 0.15,
                           legend_args = fpLegend(pos = "top")) %>%
    fp_add_lines("#555555") %>%
    fp_set_style(box = c("#377EB8", "#E69F00") %>%
                   lapply(function(x) gpar(fill = x, col = "#555555")),
                 default = gpar(vertices = TRUE))}


tiff("psychhosp_model1and2_forestplot.tiff", units = "in", width = 7, height = 4, res = 300)

forestplot_combined(combined_forest_model, "Adjusted HR")

dev.off()

tiff("psychhosp_model1and2_forestplot_pp.tiff", units = "in", width = 7, height = 4, res = 300)

forestplot_combined(combined_forest_model_pp, "Adjusted HR")

dev.off()

# INTERACTIONS

# Function to perform interaction model fitting and testing
fit_interaction_model <- function(interaction_var, adjustment_set, model, subgroup_name) {
  interact_formula <- paste("Surv(psychhosp_main_fuptime_days_cen, hospadmit_2yr_any) ~ trt_group *", interaction_var, "+", adjustment_set)
  
  interaction_result <- as.mitml.result(with(tte_imputed, coxph(as.formula(interact_formula)))) # adjusted analysis, using MI dataset
  
  test_result <- testModels(interaction_result, model, method = c("D1"), use = "wald")
  interactions_df <- as.data.frame(test_result$test) %>%
    mutate(subgroup = subgroup_name)
  
  return(interactions_df)
}

# Standard adjustment set
adjustment_set_str <- paste(adjustment_set, collapse = " + ")

# Adjustment set with prior hosp
adjustment_set_priorhosp <- paste(c(adjustment_set, "priorhosp_any"), collapse = " + ")

# Main models
main_model <- as.mitml.result(psychhosp_itt_cox_adjusted) # itt model
phosp_model <- as.mitml.result(psychhosp_itt_cox_adjusted_model2) # itt model with prior hosp adjustment

# Define interaction variables and subgroup names
interactions_list <- list(
  list(var = "diag_prev", name = "Diagnosis"),
  list(var = "ethnicity_cat_cprdhes", name = "Ethnicity"),
  list(var = "gender", name = "Gender"),
  list(var = "age_atprevcohortentry", name = "Age"),
  list(var = "apuse_prior2years", name = "Prior non-study AP use"),
  list(var = "prevcohortentry_year", name = "Time-period"))

# Initialize interaction results
interactions <- do.call(rbind, lapply(interactions_list, function(x) {
  fit_interaction_model(x$var, adjustment_set_str, main_model, x$name)
}))

# Interactions using adjustment set including prior hosp
interactions_phosp <- do.call(rbind, lapply(interactions_list, function(x) {
  fit_interaction_model(x$var, adjustment_set_priorhosp, phosp_model, x$name)
}))

# Finalize results for standard adjustment set
interactions <- interactions %>%
  rename(pvalue = `P(>F)`) %>%
  mutate(pvalue = round(pvalue, 3),
         pvalue = ifelse(pvalue < 0.001, "<0.001",
                         ifelse(pvalue == 1.000, ">0.999", format(pvalue, nsmall = 3)))) # use "<0.001" or ">0.999" for really small/large p values

# Finalize results for prior hosp adjustment set
interactions_phosp <- interactions_phosp %>%
  rename(pvalue_phosp = `P(>F)`) %>%
  mutate(pvalue_phosp = round(pvalue_phosp, 3),
         pvalue_phosp = ifelse(pvalue_phosp < 0.001, "<0.001",
                               ifelse(pvalue_phosp == 1.000, ">0.999", format(pvalue_phosp, nsmall = 3)))) # use "<0.001" or ">0.999" for really small/large p values

psychhosp_interaction_table <- interactions %>%
  select(subgroup, pvalue) %>%
  left_join(select(interactions_phosp, subgroup, pvalue_phosp), by = "subgroup") %>%
  gt() %>%
  tab_spanner(label = "Wald test P value",  # add spanner for N column
              columns = c(pvalue, pvalue_phosp), id = "pp_span") %>%
  cols_label( # tidy column labels
    subgroup = "Subgroup",
    pvalue = "Model 1",
    pvalue_phosp = "Model 2 (with prior hospitalisation)") %>%
  cols_align(
    align = "center", # center allign columns
    columns = everything())

# Robust standard errors ####

# Add pracid to imputed data
robust_data <- complete(tte_imputed, action = "long", include = TRUE) %>%
  left_join(tte_prevalent_cohort %>% select(patid, pracid), by = "patid") %>%
  as.mids()

# Model 1

psychhosp_itt_cox_adjusted_robust <- with(robust_data, {
  # Create the formula with clustering by pracid
  adjusted_formula_with_cluster <- update(as.formula(psychhosp_itt_adjusted_formula), . ~ . + cluster(pracid))
  
  # Fit the Cox proportional hazards model
  coxph(adjusted_formula_with_cluster)})

# Model 2
psychhosp_itt_cox_adjusted_model2_robust <- with(robust_data, {
  # Create the formula with clustering by pracid
  adjusted_formula_with_cluster_model2 <- update(as.formula(psychhosp_itt_adjusted_formula_model2), . ~ . + cluster(pracid))
  
  # Fit the Cox proportional hazards model
  coxph(adjusted_formula_with_cluster_model2)})

# Define a function to handle the Cox model fitting and result processing
create_table_robust <- function(model, model_name) {
  # Process the results
  result_table <- as.data.frame(model %>%
                                  tbl_regression(include = trt_group, exp = FALSE)) %>%
    clean_names() %>%
    separate(x95_percent_ci, into = c("ci_high", "ci_low"), sep = ", ") %>%
    mutate_at(vars(2:5), as.numeric) %>% # convert to numeric
    mutate_at(vars(2:4), ~round(exp(. * -1), 2)) %>% # inverse the log hazards, exponentiate it, round to 2dp
    select(group = characteristic, hr = log_hr, ci_low, ci_high, p_value) %>%
    filter(row_number() > 2) %>%
    mutate(model = model_name)
  
  return(result_table)}

# Apply the function to both models
robust_results_table_model1 <- create_table_robust(psychhosp_itt_cox_adjusted_robust, "Model 1")
robust_results_table_model2 <- create_table_robust(psychhosp_itt_cox_adjusted_model2_robust, "Model 2")

psychhosp_pp_cox_adjusted_robust <- with(robust_data, {
  # Create the formula with clustering by pracid
  adjusted_formula_with_cluster <- update(as.formula(psychhosp_pp_adjusted_formula), . ~ . + cluster(pracid))
  
  # Fit the Cox proportional hazards model
  coxph(adjusted_formula_with_cluster)})

# Model 2
psychhosp_pp_cox_adjusted_model2_robust <- with(robust_data, {
  # Create the formula with clustering by pracid
  adjusted_formula_with_cluster_model2 <- update(as.formula(psychhosp_pp_adjusted_formula_model2), . ~ . + cluster(pracid))
  
  # Fit the Cox proportional hazards model
  coxph(adjusted_formula_with_cluster_model2)})

# Apply the function to both models
robust_results_table_pp_model1 <- create_table_robust(psychhosp_pp_cox_adjusted_robust, "Model 1")
robust_results_table__pp_model2 <- create_table_robust(psychhosp_pp_cox_adjusted_model2_robust, "Model 2")

psych_hosp_robust_itt <- robust_results_table_model1 %>%
  bind_rows(robust_results_table_model2)

psych_hosp_robust_pp <- robust_results_table_pp_model1 %>%
  bind_rows(robust_results_table__pp_model2)

psych_hosp_robust <- psych_hosp_robust_itt %>%
  left_join(psych_hosp_robust_pp, by = c("group", "model")) %>%
  gt(groupname_col = "model") %>%
  
  # Format the ITT columns with HR (CI low, CI high) [p-value] format
  fmt_number(columns = c(hr.x, ci_low.x, ci_high.x), decimals = 2) %>%
  fmt_number(columns = p_value.x, decimals = 3) %>%
  cols_merge(
    columns = c(hr.x, ci_low.x, ci_high.x, p_value.x),
    pattern = "{1} ({2}, {3}) [{4}]"
  ) %>%
  
  # Format the PP columns with HR (CI low, CI high) [p-value] format
  fmt_number(columns = c(hr.y, ci_low.y, ci_high.y), decimals = 2) %>%
  fmt_number(columns = p_value.y, decimals = 3) %>%
  cols_merge(
    columns = c(hr.y, ci_low.y, ci_high.y, p_value.y),
    pattern = "{1} ({2}, {3}) [{4}]"
  ) %>%
  
  # Rename merged columns for display
  cols_label(
    hr.x = "ITT", 
    hr.y = "PP"
  ) %>%
  
  # Add spanners for ITT and PP sections
  tab_spanner(label = "ITT", columns = hr.x) %>%
  tab_spanner(label = "PP", columns = hr.y) %>%
  
  # Italicize group names
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_row_groups())

rm(list = ls(pattern = "robust"))

#SUBGROUP ANALYSIS - PSYCHIATRIC HOSPITALISATION ####
# These are done using standard cox regression

# Subset data based on subgroups
# The below only currently works on subgroup defining variables that do not contain missing data (so ethnicity not included)

# Define subgroup defining fields
fields <- c("diag_prev", "diag_prev", "diag_prev", "gender", "gender", "age_atcohortentry_cat_two", "age_atcohortentry_cat_two", 
            "apuse_prior2years", "apuse_prior2years", "cohortentry_timeperiod", "cohortentry_timeperiod", "cohortentry_timeperiod")

# Define subgroup defining values (in the same order as above)
values <- c("schizophrenia", "bipolar", "other psychosis", "Male", "Female", "Under 50", "50 or over", 
            "Yes", "No", "2005-2009", "2010-2014", "2015+")

# Define function to filter mice imputed dataset based on subgroups (converts mids object to long format, filters dataset, converts back to mids)
extract_subgroups <- function(field, value) {
  filtered_data <- complete(tte_imputed, action = "long", include = TRUE) %>%
    filter(.data[[field]] == value & hes_apc_e == 1)
  mids_object <- as.mids(filtered_data)
  return(mids_object)}

# Initialize an empty vector to store the mids names (referred to in later analysis)
mids_names <- c()

# Loop over the fields and values
for (i in 1:length(fields)) {
  field <- fields[i]
  value <- values[i]
  
  # Convert data to mids object
  mids_object <- extract_subgroups(field, value)
  
  # Create the variable name for the mids object, remove any non alphanumeric characters (e.g. in the cohort entry year variable)
  mids_name <- paste0("tte_imputed_", field, "_", gsub("[^[:alnum:]]", "", tolower(value)))
  
  # Assign the filtered mids object to the environment
  assign(mids_name, mids_object, envir = .GlobalEnv)
  
  # Add the mids name to a vector
  mids_names <- c(mids_names, mids_name)}

# Run subgroup analyses ####
# This essentially repeats the standard analysis as above (cox for ITT and PP)

# Initialize an empty list to store the results
subgroup_results_list <- list()

# Iterate over the mids names
for (mids_name in mids_names) {
  # Fit ITT cox model
  cox_model <- with(get(mids_name), 
                    coxph(as.formula(paste("Surv(psychhosp_main_fuptime_days_cen, hospadmit_2yr_any) ~ trt_group +", paste(adjustment_set, collapse = " + ")))))
  
  cox_model_hospsens <- with(get(mids_name), 
                             coxph(as.formula(paste("Surv(psychhosp_main_fuptime_days_cen, hospadmit_2yr_any) ~ trt_group + priorhosp_any +", paste(adjustment_set, collapse = " + ")))))
  
  # Fit PP cox model
  cox_model_pp <- with(get(mids_name), 
                       coxph(as.formula(paste("Surv(psychhosp_pp_fuptime_days_cen, hospadmit_2yr_pp) ~ trt_group +", paste(adjustment_set, collapse = " + ")))))
  
  # Fit PP cox model
  cox_model_pp_hospsens <- with(get(mids_name), 
                                coxph(as.formula(paste("Surv(psychhosp_pp_fuptime_days_cen, hospadmit_2yr_pp) ~ trt_group + priorhosp_any +", paste(adjustment_set, collapse = " + ")))))
  
  # Store all results in the list
  subgroup_results_list[[mids_name]] <- list(cox_model = cox_model, 
                                             cox_model_pp = cox_model_pp,
                                             cox_model_itt_hospsens = cox_model_hospsens,
                                             cox_model_pp_hospsens = cox_model_pp_hospsens)}

# Remove the interim mids objects
rm(list = mids_names)

# Ethnicity subgroup analysis

# Create sepatate data frames for each ethnicity
tte_imputed_ethnicity_black <- complete(tte_imputed, action = "long", include = TRUE) %>%
  filter(ethnicity_cat_cprdhes == "black" & hes_apc_e == 1 & .imp > 0)

tte_imputed_ethnicity_asian <- complete(tte_imputed, action = "long", include = TRUE) %>%
  filter(ethnicity_cat_cprdhes == "asian" & hes_apc_e == 1 & .imp > 0)

tte_imputed_ethnicity_white <- complete(tte_imputed, action = "long", include = TRUE) %>%
  filter(ethnicity_cat_cprdhes == "white" & hes_apc_e == 1 & .imp > 0)

tte_imputed_ethnicity_mixedother <- complete(tte_imputed, action = "long", include = TRUE) %>%
  filter((ethnicity_cat_cprdhes == "other" | ethnicity_cat_cprdhes == "mixed") & hes_apc_e == 1 & .imp > 0) %>% # too few to consider mixed and other separately
  mutate(ethnicity_cat_cprdhes = "Mixed/Other")

labels <- c("black", "asian", "white", "mixedother")

# Update adjustment set to remove ethnicity (did not seem necessary to do this in the other subgroups? but this one did not run without doing this)
variable_to_remove <- "ethnicity_cat_cprdhes"
adjustment_set_withoutethnicity <- adjustment_set[adjustment_set != variable_to_remove]

# Create an empty list to store the analysis results
eth_results_list <- list()
eth_results_list_pp <- list()

eth_results_list_hospsens <- list()
eth_results_list_pp_hospsens <- list()

# Loop over each data frame and its corresponding label
for (current_label in labels) {
  # Get the current data frame
  current_df <- get(paste("tte_imputed_ethnicity_", current_label, sep = ""))
  
  # Get unique values of .imp column
  unique_imps <- unique(current_df$.imp)
  
  # Loop over each unique .imp value
  for (imp_value in unique_imps) {
    # Subset the imputed dataset based on the .imp value
    subset_data <- current_df[current_df$.imp == imp_value, ]
    
    # Perform survival analysis on the subsetted data (ITT and PP)
    survival_result <- coxph(as.formula(paste("Surv(psychhosp_main_fuptime_days_cen, hospadmit_2yr_any) ~ trt_group +", paste(adjustment_set_withoutethnicity, collapse = " + "))), data = subset_data)
    survival_result_pp <- coxph(as.formula(paste("Surv(psychhosp_pp_fuptime_days_cen, hospadmit_2yr_pp) ~ trt_group +", paste(adjustment_set_withoutethnicity, collapse = " + "))), data = subset_data)
    
    survival_result_hospsens <- coxph(as.formula(paste("Surv(psychhosp_main_fuptime_days_cen, hospadmit_2yr_any) ~ trt_group + priorhosp_any +", paste(adjustment_set_withoutethnicity, collapse = " + "))), data = subset_data)
    survival_result_pp_hospsens <- coxph(as.formula(paste("Surv(psychhosp_pp_fuptime_days_cen, hospadmit_2yr_pp) ~ trt_group + priorhosp_any + ", paste(adjustment_set_withoutethnicity, collapse = " + "))), data = subset_data)
    
    # Store the results in the list
    eth_results_list[[paste(current_label, imp_value, sep = "_")]] <- survival_result
    eth_results_list_pp[[paste(current_label, imp_value, sep = "_")]] <- survival_result_pp
    eth_results_list_hospsens[[paste(current_label, imp_value, sep = "_")]] <- survival_result_hospsens
    eth_results_list_pp_hospsens[[paste(current_label, imp_value, sep = "_")]] <- survival_result_pp_hospsens}
  
  # Pool the results for the current label using Rubin's rules
  pooled_results <- MIcombine(eth_results_list[paste(current_label, unique_imps, sep = "_")])
  pooled_results_pp <- MIcombine(eth_results_list_pp[paste(current_label, unique_imps, sep = "_")])
  pooled_results_hospsens <- MIcombine(eth_results_list_hospsens[paste(current_label, unique_imps, sep = "_")])
  pooled_results_pp_hospsens <- MIcombine(eth_results_list_pp_hospsens[paste(current_label, unique_imps, sep = "_")])
  
  # Extract pooled coefficients and CIs
  pooled_coefs <- coef(pooled_results)
  pooled_ci <- confint(pooled_results)
  
  pooled_coefs_pp <- coef(pooled_results_pp)
  pooled_ci_pp <- confint(pooled_results_pp)
  
  pooled_coefs_itt_hospsens <- coef(pooled_results_hospsens)
  pooled_ci_itt_hospsens <- confint(pooled_results_hospsens)
  
  pooled_coefs_pp_hospsens <- coef(pooled_results_pp_hospsens)
  pooled_ci_pp_hospsens <- confint(pooled_results_pp_hospsens)
  
  # Calculate variance-covariance matrix for pooled estimates
  pooled_vcov <- vcov(pooled_results)
  pooled_vcov_pp <- vcov(pooled_results_pp)
  pooled_vcov_itt_hospsens <- vcov(pooled_results_hospsens)
  pooled_vcov_pp_hospsens <- vcov(pooled_results_pp_hospsens)
  
  # Calculate standard errors for pooled estimates
  pooled_se <- sqrt(diag(pooled_vcov))
  pooled_se_pp <- sqrt(diag(pooled_vcov_pp))
  pooled_se_itt_hospsens <- sqrt(diag(pooled_vcov_itt_hospsens))
  pooled_se_pp_hospsens <- sqrt(diag(pooled_vcov_pp_hospsens))
  
  # Calculate pooled p-values
  pooled_p_values <- 2 * (1 - pnorm(abs(pooled_coefs / pooled_se)))
  pooled_p_values_pp <- 2 * (1 - pnorm(abs(pooled_coefs_pp / pooled_se_pp)))
  pooled_p_values_itt_hospsens <- 2 * (1 - pnorm(abs(pooled_coefs_itt_hospsens / pooled_se_itt_hospsens)))
  pooled_p_values_pp_hospsens <- 2 * (1 - pnorm(abs(pooled_coefs_pp_hospsens / pooled_se_pp_hospsens)))
  
  # Create a data frame for pooled results (exponentiated values only) for the current label
  
  pooled_df <- data.frame(
    rowname = rownames(pooled_ci),
    exp_coefficient = round(exp(pooled_coefs * -1), 2),
    exp_CI_upper = round(exp(pooled_ci[, 1] * -1), 2),  # Round the exp_CI_upper column
    exp_CI_lower = round(exp(pooled_ci[, 2] * -1), 2),  # Round the exp_CI_lower column
    pvalue = round(pooled_p_values, 3),
    Label = rep(current_label, nrow(pooled_ci)))
  
  pooled_df_pp <- data.frame(
    rowname = rownames(pooled_ci_pp),
    exp_coefficient = round(exp(pooled_coefs_pp * -1), 2),
    exp_CI_upper = round(exp(pooled_ci_pp[, 1] * -1), 2),  # Round the exp_CI_upper column
    exp_CI_lower = round(exp(pooled_ci_pp[, 2] * -1), 2),  # Round the exp_CI_lower column
    pvalue = round(pooled_p_values_pp, 3),
    Label = rep(current_label, nrow(pooled_ci_pp)))
  
  pooled_df_itt_hospsens <- data.frame(
    rowname = rownames(pooled_ci_itt_hospsens),
    exp_coefficient = round(exp(pooled_coefs_itt_hospsens * -1), 2),
    exp_CI_upper = round(exp(pooled_ci_itt_hospsens[, 1] * -1), 2),  # Round the exp_CI_upper column
    exp_CI_lower = round(exp(pooled_ci_itt_hospsens[, 2] * -1), 2),  # Round the exp_CI_lower column
    pvalue = round(pooled_p_values_itt_hospsens, 3),
    Label = rep(current_label, nrow(pooled_ci_itt_hospsens)))
  
  pooled_df_pp_hospsens <- data.frame(
    rowname = rownames(pooled_ci_pp_hospsens),
    exp_coefficient = round(exp(pooled_coefs_pp_hospsens * -1), 2),
    exp_CI_upper = round(exp(pooled_ci_pp_hospsens[, 1] * -1), 2),  # Round the exp_CI_upper column
    exp_CI_lower = round(exp(pooled_ci_pp_hospsens[, 2] * -1), 2),  # Round the exp_CI_lower column
    pvalue = round(pooled_p_values_pp_hospsens, 3),
    Label = rep(current_label, nrow(pooled_ci_pp_hospsens)))
  
  # Filter rows where rowname contains "trt_group"
  pooled_df <- pooled_df[grep("trt_group", pooled_df$rowname), ]
  pooled_df_pp <- pooled_df_pp[grep("trt_group", pooled_df_pp$rowname), ]
  pooled_df_itt_hospsens <- pooled_df_itt_hospsens[grep("trt_group", pooled_df_itt_hospsens$rowname), ]
  pooled_df_pp_hospsens <- pooled_df_pp_hospsens[grep("trt_group", pooled_df_pp_hospsens$rowname), ]
  
  # Assign the pooled_df with the label as a separate data frame
  assign(paste(current_label, "_pooled_df", sep = ""), pooled_df)
  assign(paste(current_label, "_pooled_df_pp", sep = ""), pooled_df_pp)
  assign(paste(current_label, "_pooled_df_itt_hospsens", sep = ""), pooled_df_itt_hospsens)
  assign(paste(current_label, "_pooled_df_pp_hospsens", sep = ""), pooled_df_pp_hospsens)}

# Tables

# ITT

age_50plus <- reg_table(subgroup_results_list$tte_imputed_age_atcohortentry_cat_two_50orover$cox_model) %>%
  mutate(category = "50+", subgroup = "Age")

age_50less <- reg_table(subgroup_results_list$tte_imputed_age_atcohortentry_cat_two_under50$cox_model) %>%
  mutate(category = "<50", subgroup = "Age")

schiz <- reg_table(subgroup_results_list$tte_imputed_diag_prev_schizophrenia$cox_model) %>%
  mutate(category = "Schizophrenia", subgroup = "Diagnosis")

bipolar <- reg_table(subgroup_results_list$tte_imputed_diag_prev_bipolar$cox_model) %>%
  mutate(category = "Bipolar disorder", subgroup = "Diagnosis")

otherpsychosis <- reg_table(subgroup_results_list$tte_imputed_diag_prev_otherpsychosis$cox_model) %>%
  mutate(category = "Other psychosis", subgroup = "Diagnosis")

female <- reg_table(subgroup_results_list$tte_imputed_gender_female$cox_model) %>%
  mutate(category = "Female", subgroup = "Gender")

male <- reg_table(subgroup_results_list$tte_imputed_gender_male$cox_model) %>%
  mutate(category = "Male", subgroup = "Gender")

priorapuse_yes <- reg_table(subgroup_results_list$tte_imputed_apuse_prior2years_yes$cox_model) %>%
  mutate(category = "Yes", subgroup = "Prior non-study AP use")

priorapuse_no <- reg_table(subgroup_results_list$tte_imputed_apuse_prior2years_no$cox_model) %>%
  mutate(category = "No", subgroup = "Prior non-study AP use")

timeperiod_2005 <- reg_table(subgroup_results_list$tte_imputed_cohortentry_timeperiod_20052009$cox_model) %>%
  mutate(category = "2005-2009", subgroup = "Time-period")

timeperiod_2010 <- reg_table(subgroup_results_list$tte_imputed_cohortentry_timeperiod_20102014$cox_model) %>%
  mutate(category = "2010-2014", subgroup = "Time-period")

timeperiod_2015 <- reg_table(subgroup_results_list$tte_imputed_cohortentry_timeperiod_2015$cox_model) %>%
  mutate(category = "2015+", subgroup = "Time-period")

extract_subgroup_results_ethnicity <- function(data) {
  data %>%
    rename(hr = exp_coefficient,
           ci_low = exp_CI_lower,
           ci_high = exp_CI_upper,
           p_value = pvalue,
           category = Label) %>%
    mutate(subgroup = "Ethnicity",
           group = gsub("trt_group", "", rowname),
           label = case_when(group == "Olanzapine" ~ "Aripiprazole vs. Olanzapine",
                             group == "Quetiapine" ~ "Aripiprazole vs. Quetiapine",
                             group == "Risperidone" ~ "Aripiprazole vs. Risperidone",
                             TRUE ~ "Unknown"))}

white <- extract_subgroup_results_ethnicity(white_pooled_df)
asian <- extract_subgroup_results_ethnicity(asian_pooled_df)
black <- extract_subgroup_results_ethnicity(black_pooled_df)
mixedother <- extract_subgroup_results_ethnicity(mixedother_pooled_df)

psych_subgroups <- bind_rows(age_50plus, age_50less, schiz, bipolar, otherpsychosis, 
                             priorapuse_yes, priorapuse_no, female, male,
                             timeperiod_2005, timeperiod_2010, timeperiod_2015) %>%
  bind_rows(asian, black, mixedother, white) %>%
  select(subgroup, label, category, term = group, hr, ci_low, ci_high, p_value) %>%
  mutate(category = str_to_sentence(category),
         category = case_when(category == "Mixedother" ~ "Mixed/Other",
                              category == "Other psychosis" ~ "Other psychoses",
                              TRUE ~ category)) %>%
  arrange(subgroup, category) %>%
  left_join(select(interactions, subgroup, pvalue), by = "subgroup") %>%
  group_by(label, subgroup) %>%
  mutate(
    # Only show subgroup and p-value for first row in group
    subgroup_display = if_else(row_number() == 1, subgroup, ""),
    p_value_display = if_else(row_number() == 1, as.character(pvalue), "")) %>%
  ungroup()


# ITT - additional adjustment for prior hospitalisation

age_50plus_model2 <- reg_table(subgroup_results_list$tte_imputed_age_atcohortentry_cat_two_50orover$cox_model_itt_hospsens) %>%
  mutate(category = "50+", subgroup = "Age")

age_50less_model2 <- reg_table(subgroup_results_list$tte_imputed_age_atcohortentry_cat_two_under50$cox_model_itt_hospsens) %>%
  mutate(category = "<50", subgroup = "Age")

schiz_model2 <- reg_table(subgroup_results_list$tte_imputed_diag_prev_schizophrenia$cox_model_itt_hospsens) %>%
  mutate(category = "Schizophrenia", subgroup = "Diagnosis")

bipolar_model2 <- reg_table(subgroup_results_list$tte_imputed_diag_prev_bipolar$cox_model_itt_hospsens) %>%
  mutate(category = "Bipolar disorder", subgroup = "Diagnosis")

otherpsychosis_model2 <- reg_table(subgroup_results_list$tte_imputed_diag_prev_otherpsychosis$cox_model_itt_hospsens) %>%
  mutate(category = "Other psychosis", subgroup = "Diagnosis")

female_model2 <- reg_table(subgroup_results_list$tte_imputed_gender_female$cox_model_itt_hospsens) %>%
  mutate(category = "Female", subgroup = "Gender")

male_model2 <- reg_table(subgroup_results_list$tte_imputed_gender_male$cox_model_itt_hospsens) %>%
  mutate(category = "Male", subgroup = "Gender")

priorapuse_yes_model2 <- reg_table(subgroup_results_list$tte_imputed_apuse_prior2years_yes$cox_model_itt_hospsens) %>%
  mutate(category = "Yes", subgroup = "Prior non-study AP use")

priorapuse_no_model2 <- reg_table(subgroup_results_list$tte_imputed_apuse_prior2years_no$cox_model_itt_hospsens) %>%
  mutate(category = "No", subgroup = "Prior non-study AP use")

timeperiod_2005_model2 <- reg_table(subgroup_results_list$tte_imputed_cohortentry_timeperiod_20052009$cox_model_itt_hospsens) %>%
  mutate(category = "2005-2009", subgroup = "Time-period")

timeperiod_2010_model2 <- reg_table(subgroup_results_list$tte_imputed_cohortentry_timeperiod_20102014$cox_model_itt_hospsens) %>%
  mutate(category = "2010-2014", subgroup = "Time-period")

timeperiod_2015_model2 <- reg_table(subgroup_results_list$tte_imputed_cohortentry_timeperiod_2015$cox_model_itt_hospsens) %>%
  mutate(category = "2015+", subgroup = "Time-period")

white_model2 <- extract_subgroup_results_ethnicity(white_pooled_df_itt_hospsens)
asian_model2 <- extract_subgroup_results_ethnicity(asian_pooled_df_itt_hospsens)
black_model2 <- extract_subgroup_results_ethnicity(black_pooled_df_itt_hospsens)
mixedother_model2 <- extract_subgroup_results_ethnicity(mixedother_pooled_df_itt_hospsens)

psych_subgroups_model2 <- bind_rows(age_50plus_model2, age_50less_model2, schiz_model2, bipolar_model2, otherpsychosis_model2, 
                                      priorapuse_yes_model2, priorapuse_no_model2, female_model2, male_model2,
                                      timeperiod_2005_model2, timeperiod_2010_model2, timeperiod_2015_model2) %>%
  bind_rows(asian_model2, black_model2, mixedother_model2, white_model2) %>%
  select(subgroup, category, label, term = group, hr, ci_low, ci_high, p_value) %>%
  mutate(category = str_to_sentence(category),
         category = case_when(category == "Mixedother" ~ "Mixed/Other",
                              category == "Other psychosis" ~ "Other psychoses",
                              TRUE ~ category)) %>%
  arrange(subgroup, category) %>%
  left_join(select(interactions_phosp, subgroup, pvalue_phosp), by = "subgroup") %>%
  group_by(label, subgroup) %>%
  mutate(
    # Only show subgroup and p-value for first row in group
    subgroup_display = if_else(row_number() == 1, subgroup, ""),
    p_value_display = if_else(row_number() == 1, as.character(pvalue_phosp), "")) %>%
  ungroup()

# Forest plots

# Set wd to save plots
setwd(paste0(path))

forestplot_subgroup <- function(data, xlab) {
  
  # Add bold formatting to subgroup values
  
  # Create the forest plot with clean labels
  data %>%
    group_by(label) %>%
    forestplot::forestplot(
      fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI, fpDrawDiamondCI),
      labeltext = c(subgroup_display, category, p_value_display),
      clip = c(0.40, 2.00),
      xlog = TRUE,
      mean = hr,
      lower = ci_low,
      upper = ci_high,
      vertices = TRUE,
      grid = FALSE,
      xlab = "Adjusted HR",
      ci.vertices.height = .05,
      boxsize = 0.15,
      # Format text appearance, including increasing legend text size
      txt_gp = fpTxtGp(
        label = gpar(cex = 1),
        xlab = gpar(cex = 1.25),
        title = gpar(cex = 1, fontface = "bold"),  # Bold headers
        ticks = gpar(cex = 0.9),
        legend = gpar(cex = 1.25)  # Increase legend text size
      )
    ) %>%
    fp_add_lines("#555555") %>%
    fp_add_header(
      "Subgroup", "Category", 
      expression(bold(paste("Wald ", bolditalic("P"))))) %>%
    fp_set_style(
      box = c("#377EB8", "#E69F00", "#D55E00") %>%
        lapply(function(x) gpar(fill = x, col = "#555555")),
      default = gpar(vertices = TRUE)
    ) %>%
    fp_set_zebra_style("#F5F9F9")}

tiff("psychhosp_subgroup_forestplot.tiff", units = "in", width = 10, height = 13, res = 300)

forestplot_subgroup(psych_subgroups, "Adjusted HR")

dev.off()

tiff("psychhosp_subgroup_forestplot_model2.tiff", units = "in", width = 10, height = 13, res = 300)

forestplot_subgroup(psych_subgroups_model2, "Adjusted HR")

dev.off()

# Definition change sensitivity analysis ####

psychhosp_definitionsens_cox_unadj <-  coxph(psychhosp_definitionsens_formula, data = tte_prevalent_cohort) # Cox regression, unadjusted
psychhosp_definitionsens_cox_adjusted_model1 <- with(tte_imputed, coxph(as.formula(psychhosp_definitionsens_adjusted_formula))) # adjusted, full MI dataset

psychhosp_definitionsens_adjusted_formula_model2 <- paste0(deparse1(psychhosp_definitionsens_formula), " + ", paste(adjustment_set, "+ priorhosp_any", collapse = " + "))
psychhosp_definitionsens_cox_adjusted_model2 <- with(tte_imputed, coxph(as.formula(psychhosp_definitionsens_adjusted_formula_model2))) # adjusted, full MI dataset

definitionchange_unadj <- reg_table(psychhosp_definitionsens_cox_unadj)
definitionchange_adj_model1 <- reg_table(psychhosp_definitionsens_cox_adjusted_model1)
definitionchange_adj_model2 <- reg_table(psychhosp_definitionsens_cox_adjusted_model2)

definitionchange <- definitionchange_unadj %>%
  select(label, unadjusted = result) %>%
  left_join(select(definitionchange_adj_model1, label, model1 = result), by = "label") %>%
  left_join(select(definitionchange_adj_model2, label, model2 = result), by = "label") %>%
  rename(Comparison = label) %>%
  gt() %>%
  tab_footnote("Mode931 adjusts for pre-specified covariates. Model 2 adjusts for pre-specified covariates plus prior psychiatric hospitalisation")

# IPCW
load(file = "tte_prevalent_cohort.Rdata")
load(file = "tte_imputed_long.Rdata")

tte_imputed_long <- tte_imputed_long %>%
  left_join(select(tte_prevalent_cohort, patid, pp_censor, sixm_totalcholesterol_flag, oney_totalcholesterol_flag), by = "patid") %>%
  filter(hes_apc_e == 1) %>%
  select(-hes_apc_e) %>%
  mutate(adherence = case_when(hospadmit_2yr_pp == 1 ~ 1,
                               psychhosp_pp_fuptime_days_cen >= 730 ~ 1,
                               TRUE ~ 0))

adjustment_set_2y <- c("age_atprevcohortentry", "gender", "ethnicity_cat_cprdhes", "pat_2019imd_quintile", "diag_prev", "region", "prevcohortentry_year", "gpconsults_last6m", "smoking_status_cat", 
                       "prioralcoholabuse", "priorsubstanceabuse", "priordyslipidaemia", "priordiabetes", "priorhypertension", "priorcerebrovasculardisease", "priormyocardialinfarction", "priorrenaldisease", "priorliverdisease", 
                       "apuse_prior2years", "lipiddrugs_prior2years", "hypertensiondrugs_prior2years", "antidiabetics_prior2years", "antidepressant_prior2years", "moodstab_prior2years",
                       "baseline_totalcholesterol", "baseline_ldl", "baseline_hdl", "baseline_triglycerides", "baseline_systolicbp", "baseline_diastolicbp", "baseline_glucose", "baseline_hba1c", "baseline_weightkg", "baseline_bmi_cat",
                       "sixm_totalcholesterol", "sixm_ldl", "sixm_hdl", "sixm_triglycerides", "sixm_systolicbp", "sixm_diastolicbp", "sixm_glucose", "sixm_hba1c", "sixm_weightkg",
                       "oney_totalcholesterol", "oney_ldl", "oney_hdl", "oney_triglycerides", "oney_systolicbp", "oney_diastolicbp", "oney_glucose", "oney_hba1c", "oney_weightkg",
                       "sixm_totalcholesterol_flag", "oney_totalcholesterol_flag")

calculate_censor_weights <- function(data, adjustment_vars) {
  results_list <- lapply(1:25, function(imp_value) {
    # Prepare data for current imputation
    test_data <- data %>% 
      filter(.imp == imp_value) %>%
      select(patid, trt_group, adherence, all_of(adjustment_vars)) %>%
      drop_na(adherence, trt_group, all_of(adjustment_vars))
    
    # Fit both censoring models
    censor_models <- list(
      full = glm(
        as.formula(paste("adherence == 1 ~ trt_group +", paste(adjustment_vars, collapse = " + "))),
        data = test_data,
        family = binomial(link = "logit")),
      simple = glm(
        adherence == 1 ~ trt_group,
        data = test_data,
        family = binomial(link = "logit")))
    
    # Calculate weights with trimming
    test_data %>%
      mutate(
        censor_weight_denom = predict(censor_models$full, type = "response"),
        censor_weight_num = predict(censor_models$simple, type = "response"),
        stabilized_censor_weight = (1 - censor_weight_num) / (1 - censor_weight_denom),
        trimmed_stabilized_censor_weight = ifelse(
          is.na(stabilized_censor_weight), 
          NA, 
          pmin(
            stabilized_censor_weight, 
            quantile(stabilized_censor_weight, 0.995, na.rm = TRUE)
          )
        ),
        .imp = imp_value
      ) %>%
      select(patid, .imp, trt_group, adherence, 
             censor_weight_denom, censor_weight_num, 
             stabilized_censor_weight, trimmed_stabilized_censor_weight)
  })
  
  bind_rows(results_list)}

results_2y <- calculate_censor_weights(tte_imputed_long, adjustment_set_2y)

# Functions

# Function to calculate pooled statistics
calculate_group_stats <- function(data, weight_var) {
  data %>%
    group_by(.imp, trt_group) %>%
    summarise(mean_weight = mean(!!sym(weight_var)),
              within_var = var(!!sym(weight_var)),
              min_weight = min(!!sym(weight_var)),
              max_weight = max(!!sym(weight_var)),
              .groups = "drop") %>%
    group_by(trt_group) %>%
    summarise(pooled_mean = mean(mean_weight),
              pooled_sd = sqrt(mean(within_var) + (1 + (1 / 25)) * var(mean_weight)),
              min = min(min_weight),
              max = max(max_weight),
              .groups = "drop") %>%
    mutate(`Mean (SD)` = sprintf("%.2f (%.2f)", pooled_mean, pooled_sd),
           Min = sprintf("%.2f", min),
           Max = sprintf("%.2f", max)) %>%
    select(trt_group, `Mean (SD)`, Min, Max) %>%
    gt() %>%
    cols_label(trt_group = "Treatment Group") %>%
    tab_options(table.font.size = 14,
                column_labels.font.size = 15,
                data_row.padding = px(5)) %>%
    tab_style(style = list(cell_text(align = "center")),
              locations = cells_body(columns = everything())) %>%
    tab_style(style = cell_borders(sides = "bottom", color = "grey80"),
              locations = cells_column_labels())}

group_stats_censor_2y <- calculate_group_stats(results_2y, "trimmed_stabilized_censor_weight")

# Remove objects matching specific naming patterns
rm(list = ls(pattern = "group_stats_censor|censor_weights_meansd_plot|censor_weights_dist_plot|
             censor_weights_dist_plot|unweightedcensoring_distribution_plot|group_stats_denom"))

# Cox REGRESSION, with covariate adjustment and inverse probability of censoring weights

# Add weights to mids object
tte_imputed_pp_psychhosp <- tte_imputed_long %>%
  left_join(results_2y, by = c("patid", "trt_group", ".imp", "adherence")) %>%
  as.mids()

# PP - cox regression - Post-hoc sensitivity analysis: additional adjustment for prior hospitalisation and IPCW
psychhosp_pp_formula <- as.formula("Surv(time = psychhosp_pp_fuptime_days_cen, event = hospadmit_2yr_pp) ~ trt_group") # PP formula

adjustment_set <- c("age_atprevcohortentry", "gender", "ethnicity_cat_cprdhes", "pat_2019imd_quintile", "diag_prev", "region", "prevcohortentry_year", "gpconsults_last6m", "smoking_status_cat", 
                    "prioralcoholabuse", "priorsubstanceabuse", "priordyslipidaemia", "priordiabetes", "priorhypertension", "priorcerebrovasculardisease", "priormyocardialinfarction", "priorrenaldisease", "priorliverdisease", 
                    "apuse_prior2years", "lipiddrugs_prior2years", "hypertensiondrugs_prior2years", "antidiabetics_prior2years", "antidepressant_prior2years", "moodstab_prior2years",
                    "baseline_totalcholesterol", "baseline_ldl", "baseline_hdl", "baseline_triglycerides", "baseline_systolicbp", "baseline_diastolicbp", "baseline_glucose", "baseline_hba1c", "baseline_weightkg", "baseline_bmi_cat")
adjustment_set_priorhosp <- c(adjustment_set, "priorhosp_any")

psychhosp_pp_adjusted_formula_model1 <- paste0(deparse1(psychhosp_pp_formula), " + ", paste(adjustment_set, collapse = " + "))
psychhosp_pp_adjusted_formula_model2 <- paste0(deparse1(psychhosp_pp_formula), " + ", paste(adjustment_set_priorhosp, collapse = " + "))

psychhosp_pp_cox_adjusted_model1_ipcw <- with(tte_imputed_pp_psychhosp, 
                                              coxph(as.formula(psychhosp_pp_adjusted_formula_model1), 
                                                    weights = trimmed_stabilized_censor_weight))

psychhosp_pp_cox_adjusted_model2_ipcw <- with(tte_imputed_pp_psychhosp, 
                                              coxph(as.formula(psychhosp_pp_adjusted_formula_model2), 
                                                    weights = trimmed_stabilized_censor_weight))

pp_model1_ipcw <- reg_table(psychhosp_pp_cox_adjusted_model1_ipcw) %>%
  select(label, Model1 = result, Model1_p = p_value)

pp_model2_ipcw <- reg_table(psychhosp_pp_cox_adjusted_model2_ipcw) %>%
  select(label, Model2 = result, Model2_p = p_value)

pp_psychhosp_ipcw <- pp_model1_ipcw %>%
  left_join(pp_model2_ipcw, by = "label") %>%
  mutate(Analysis = "IPCW") %>%
  select(Analysis, everything()) %>%
  gt() %>%
  tab_footnote(footnote = "IPCW with trimming at 99.5th percentile.")

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
packageVersion("mice") # 3.16.0
packageVersion("stringr") # 1.5.1
packageVersion("reshape2") # 1.4.4
packageVersion("webshot") # 0.5.5
packageVersion("broom") # 1.0.6
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
packageVersion("tidylog") # 1.1.0
