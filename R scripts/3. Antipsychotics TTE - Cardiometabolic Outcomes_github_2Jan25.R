# -----------------------
# ANTIPSYCHOTICS TTE STUDY
# -----------------------
# Last run: 05/11/2024

# CARDIOMETABOLIC OUTCOMES ####

# Clear memory
rm(list = ls())

# Packages
invisible(sapply(c("dplyr", "tibble", "purrr", "ggplot2", "gtsummary", "gt", "forcats", "lubridate", "tidyr", "mice", "stringr", 
                   "reshape2", "webshot", "broom", "mitools", "miceadds", "grid", "gridExtra", "car", "forestplot", "rlang", "patchwork",
                   "RColorBrewer", "parallel", "cowplot", "scales", "mitml", "janitor", "tidylog"), 
                 library, character.only = TRUE))

# Set file path
path <- anonymised

# Set working directory
setwd(paste0(path))

# Load files ####
load(file = "tte_prevalent_cohort.Rdata")
load(file = "tte_imputed_long.Rdata")
tte_imputed <- as.mids(tte_imputed_long)

# DEFINE VECTORS AND FUNCTIONS ####

# Define adjustment variables
adjustment_set <- c("age_atprevcohortentry", "gender", "ethnicity_cat_cprdhes", "pat_2019imd_quintile", "diag_prev", "region", "prevcohortentry_year", "gpconsults_last6m", "smoking_status_cat", 
                    "prioralcoholabuse", "priorsubstanceabuse", "priordyslipidaemia", "priordiabetes", "priorhypertension", "priorcerebrovasculardisease", "priormyocardialinfarction", "priorrenaldisease", "priorliverdisease", 
                    "apuse_prior2years", "lipiddrugs_prior2years", "hypertensiondrugs_prior2years", "antidiabetics_prior2years", "antidepressant_prior2years", "moodstab_prior2years",
                    "baseline_totalcholesterol", "baseline_ldl", "baseline_hdl", "baseline_triglycerides", "baseline_systolicbp", "baseline_diastolicbp", "baseline_glucose", "baseline_hba1c", "baseline_weightkg", "baseline_bmi_cat")

# Define cardiometabolic outcome variables
outcome_names <- c("totalcholesterol", "ldl", "hdl", "triglycerides", "tchdlratio", "weightkg", "systolicbp", "diastolicbp", "hba1c", "glucose")
prefixes <- c("sixm_", "oney_", "twoy_")
full_names <- lapply(prefixes, function(prefix) paste0(prefix, outcome_names))
study_outcomes <- unlist(full_names)

# Define cardiometabolic baseline variables
baseline_prefix <- c("baseline_")
baseline_names <- lapply(baseline_prefix, function(prefix) paste0(baseline_prefix, outcome_names))
baseline_names <- unlist(baseline_names)

#Specify order for eventual tables, and include formatted outcome name
order <- as.data.frame(outcome_names) %>%
  tibble::rownames_to_column() %>%
  mutate(order = as.numeric(rowname)) %>%
  select(-rowname) %>%
  rename(outcome = outcome_names) %>%
  mutate(Outcome = c("Total cholesterol (mmol/L)", "LDL-C (mmol/L)", "HDL-C (mmol/L)", "Triglycerides (mmol/L)", "TC:HDL ratio", 
                     "Weight (kg)", "Systolic blood pressure (mm Hg)", "Diastolic blood pressure (mm Hg)", "HbA1c (mmol/mol)", "Glucose (mmol/L)"))

remove(prefixes, full_names, baseline_prefix)

# Function for tidying data frame
clean_df <- function(df) {
  df %>%
    mutate(outcome = as.factor(outcome)) %>% 
    filter(grepl("trt_group", term)) %>% # filter to main coefficients of interest
    mutate(term = gsub("trt_group", "", term), # remove trt_group variable name so it just states the group name alone
           estimate = round(estimate * -1, 2), # round estimates to 2 decimals
           ci_high = round(conf.low * -1, 2),
           ci_low = round(conf.high * -1, 2),
           pvalue = round(p.value, 3),
           pvalue = ifelse(pvalue < 0.001, "<0.001", format(pvalue, nsmall = 3)), # use "<0.001" for really small p values
           result = paste(format(estimate, nsmall = 2), " (", format(ci_low, nsmall = 2), ", ", format(ci_high, nsmall = 2), ")"), # paste all into one cell
           result = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", result),
           result = gsub("\\s+", " ", result), # remove unnecessary spaces in pasted result
           result = gsub(" ,", ",", result),
           term = case_when(term == "Olanzapine" ~ "Aripiprazole vs. Olanzapine",
                                    term == "Quetiapine" ~ "Aripiprazole vs. Quetiapine",
                                    term == "Risperidone" ~ "Aripiprazole vs. Risperidone",
                                    TRUE ~ "Unknown")) %>%
    select(outcome, term, estimate, ci_low, ci_high, pvalue, result) %>%
    separate(outcome, into = c("prefix", "outcome"), sep = "_", remove = FALSE) %>% # separate outcome into time point and outcome name
    mutate(prefix = ifelse(prefix == "sixm", "6m",
                           ifelse(prefix == "oney", "1y",
                                  ifelse(prefix == "twoy", "2y", prefix)))) %>%
    left_join(order, by = "outcome") %>% # merge with cleaned outcome name
    arrange(order) %>% # arrange by specificed order
    select(-outcome, -order) %>%
    select(prefix, Outcome, everything())}
  
# Function for creating a presentation ready gt table
clean_gt_table <- function(data) {
  data %>%
    select(prefix, Outcome, term, result) %>%
    rename('Time-point' = prefix) %>%
    pivot_wider(names_from = term, values_from = result, names_sep = "_") %>% # pivot long to wide 
    gt(
      groupname_col = "Outcome"
    ) %>%
    tab_style(
      style = cell_text(style = "italic"), # italicize group names
      locations = cells_row_groups()) %>%
    sub_missing(missing_text = "-") %>% # use - for empty cells
    cols_align(
      align = "center", # center align columns
      columns = everything()) %>%
    tab_footnote(
      footnote = "Values are adjusted mean difference (95% confidence interval)")}

# Extract subpopulations from imputed datasets
extract_population <- function(field, value) {
  filtered_data <- complete(tte_imputed, action = "long", include = TRUE) %>%
    filter(.data[[field]] == value)
  mids_object <- as.mids(filtered_data)
  return(mids_object)}

# Parallel
mc.cores <- detectCores() - 2 # define number of cores to use for parallel processes

# ADJUSTED, USING MI DATASET ####

# Dataset filtering (remove people who had died by the time-point)
tte_imputed_itt_6m <- extract_population("diedby6m", 0)
tte_imputed_itt_1y <- extract_population("diedby1y", 0)
tte_imputed_itt_2y <- extract_population("diedby2y", 0)

# Outcome names
sixm_itt_outcomes <- study_outcomes[grepl("sixm", study_outcomes)]
oney_itt_outcomes <- study_outcomes[grepl("oney", study_outcomes)]
twoy_itt_outcomes <- study_outcomes[grepl("twoy", study_outcomes)]

# Initialise lists to store results
results_itt_6m <- list()
results_itt_1y <- list()
results_itt_2y <- list()

# 6 months
for (outcome in sixm_itt_outcomes) { # loop over each outcome variable
  
  # fit adjusted model for each imputed dataset
  adj_models <- with(tte_imputed_itt_6m, lm(as.formula(paste(outcome, "~ trt_group + ", paste(adjustment_set, collapse = " + ")))))
  
  # extract results, combine across imputed datasets, add to list
  adj_results <- pool(adj_models)
  results_itt_6m[[outcome]] <- adj_results}

# One year
for (outcome in oney_itt_outcomes) { # loop over each outcome variable
  
  # fit adjusted model for each imputed dataset
  adj_models <- with(tte_imputed_itt_1y, lm(as.formula(paste(outcome, "~ trt_group + ", paste(adjustment_set, collapse = " + ")))))
  
  # extract results, combine across imputed datasets, add to list
  adj_results <- pool(adj_models)
    results_itt_1y[[outcome]] <- adj_results}

# Two years
for (outcome in twoy_itt_outcomes) { # loop over each outcome variable
  
  # fit adjusted model for each imputed dataset
  adj_models <- with(tte_imputed_itt_2y, lm(as.formula(paste(outcome, "~ trt_group + ", paste(adjustment_set, collapse = " + ")))))
  
  # extract results, combine across imputed datasets, add to list
  adj_results <- pool(adj_models)
    results_itt_2y[[outcome]] <- adj_results}

# Combine all time-point results into data frame
results_itt_6m_df <- results_itt_6m %>% 
  map_dfr(~ tidy(.x, conf.int = TRUE), .id = "outcome")

results_itt_1y_df <- results_itt_1y %>% 
  map_dfr(~ tidy(.x, conf.int = TRUE), .id = "outcome")

results_itt_2y_df <- results_itt_2y %>% 
  map_dfr(~ tidy(.x, conf.int = TRUE), .id = "outcome")

adjusted_itt_results_df <- rbind(results_itt_6m_df, results_itt_1y_df, results_itt_2y_df) %>%
  clean_df()

# Convert to gt table
cardiometabolic_itt_adjustedresults_MI_table <- clean_gt_table(adjusted_itt_results_df)

# PER PROTOCOL - ADJUSTED, USING MI DATASET ####

# Dataset filtering (remove people who had died, switched or discontinued by the time point)
tte_imputed_pp_6m <- extract_population("active_6m", 1)
tte_imputed_pp_1y <- extract_population("active_1y", 1)
tte_imputed_pp_2y <- extract_population("active_2y", 1)

# Outcome names
sixm_outcomes <- study_outcomes[grepl("sixm", study_outcomes)]
oney_outcomes <- study_outcomes[grepl("oney", study_outcomes)]
twoy_outcomes <- study_outcomes[grepl("twoy", study_outcomes)]

#Initialise lists to store results
results_pp_6m <- list()
results_pp_1y <- list()
results_pp_2y <- list()

# 6 months
for (outcome in sixm_outcomes) { # loop over each outcome variable
  
  # fit adjusted model for each imputed dataset
  adj_models <- with(tte_imputed_pp_6m, lm(as.formula(paste(outcome, "~ trt_group + ", paste(adjustment_set, collapse = " + ")))))
  
  # extract results, combine across imputed datasets, add to list
  adj_results <- pool(adj_models)
  results_pp_6m[[outcome]] <- adj_results}

# One year
for (outcome in oney_outcomes) { # loop over each outcome variable

  # fit adjusted model for each imputed dataset
  adj_models <- with(tte_imputed_pp_1y, lm(as.formula(paste(outcome, "~ trt_group + ", paste(adjustment_set, collapse = " + ")))))
  
  # extract results, combine across imputed datasets, add to list
  adj_results <- pool(adj_models)
  results_pp_1y[[outcome]] <- adj_results}

# Two years
for (outcome in twoy_outcomes) { # loop over each outcome variable
  
  # fit adjusted model for each imputed dataset
  adj_models <- with(tte_imputed_pp_2y, lm(as.formula(paste(outcome, "~ trt_group + ", paste(adjustment_set, collapse = " + ")))))
  
  # extract results, combine across imputed datasets, add to list
  adj_results <- pool(adj_models)
  results_pp_2y[[outcome]] <- adj_results}

# Extract and combine results into data frame
results_pp_6m_df <- results_pp_6m %>% 
  map_dfr(~ tidy(.x, conf.int = TRUE), .id = "outcome")

results_pp_1y_df <- results_pp_1y %>% 
  map_dfr(~ tidy(.x, conf.int = TRUE), .id = "outcome")

results_pp_2y_df <- results_pp_2y %>% 
  map_dfr(~ tidy(.x, conf.int = TRUE), .id = "outcome")

adjusted_pp_results_df <- rbind(results_pp_6m_df, results_pp_1y_df, results_pp_2y_df) %>%
  clean_df()

#Convert to gt table
cardiometabolic_pp_adjustedresults_MI_table <- clean_gt_table(adjusted_pp_results_df)

# Combined ITT and PP results table

# calculate ITT denominators for footnote
sixm_itt_denom <- as.data.frame(table(tte_prevalent_cohort$trt_group, tte_prevalent_cohort$diedby6m == 0))
oney_itt_denom <- as.data.frame(table(tte_prevalent_cohort$trt_group, tte_prevalent_cohort$diedby1y == 0))
twoy_itt_denom <- as.data.frame(table(tte_prevalent_cohort$trt_group, tte_prevalent_cohort$diedby2y == 0))

itt_denom <- sixm_itt_denom %>%
  rename(sixm_n = Freq) %>%
  filter(Var2 == TRUE) %>%
  left_join(oney_itt_denom, by = c("Var1", "Var2")) %>%
  rename(oney_n = Freq) %>%
  left_join(twoy_itt_denom, by = c("Var1", "Var2")) %>%
  rename(twoy_n = Freq) %>%
  select(-Var2)

itt_denominator_footnote <- paste("Intention-to-treat results reported for patients alive at each time point. At 6m, 1y and 2y, denominators for Aripiprazole were:", 
                  paste(itt_denom[1, c("sixm_n", "oney_n", "twoy_n")], collapse = ", "),
                  paste("; Olanzapine:", 
                        paste(itt_denom[2, c("sixm_n", "oney_n", "twoy_n")], collapse = ", "),
                        paste("; Quetiapine:", 
                              paste(itt_denom[3, c("sixm_n", "oney_n", "twoy_n")], collapse = ", "),
                              paste("; and Risperidone:", 
                                    paste(itt_denom[4, c("sixm_n", "oney_n", "twoy_n")], collapse = ", "),
                                    paste(", respectively.")))))

itt_denominator_footnote <- gsub("\\s+([,;])", "\\1", itt_denominator_footnote)

# calculate PP denominators for footnote
sixm_pp_denom <- as.data.frame(table(tte_prevalent_cohort$trt_group, tte_prevalent_cohort$active_6m == 1))
oney_pp_denom <- as.data.frame(table(tte_prevalent_cohort$trt_group, tte_prevalent_cohort$active_1y == 1))
twoy_pp_denom <- as.data.frame(table(tte_prevalent_cohort$trt_group, tte_prevalent_cohort$active_2y == 1))

pp_denom <- sixm_pp_denom %>%
  rename(sixm_n = Freq) %>%
  filter(Var2 == TRUE) %>%
  left_join(oney_pp_denom, by = c("Var1", "Var2")) %>%
  rename(oney_n = Freq) %>%
  left_join(twoy_pp_denom, by = c("Var1", "Var2")) %>%
  rename(twoy_n = Freq) %>%
  select(-Var2)

pp_denominator_footnote <- paste("Per-protocol results reported for patients alive and who had not switched/discontinued at each time point. Denominators for Aripiprazole were:", 
                                  paste(pp_denom[1, c("sixm_n", "oney_n", "twoy_n")], collapse = ", "),
                                  paste("; Olanzapine:", 
                                        paste(pp_denom[2, c("sixm_n", "oney_n", "twoy_n")], collapse = ", "),
                                        paste("; Quetiapine:", 
                                              paste(pp_denom[3, c("sixm_n", "oney_n", "twoy_n")], collapse = ", "),
                                              paste("; and Risperidone:", 
                                                    paste(pp_denom[4, c("sixm_n", "oney_n", "twoy_n")], collapse = ", "),
                                                    paste(", respectively.")))))

pp_denominator_footnote <- gsub("\\s+([,;])", "\\1", pp_denominator_footnote)
remove(sixm_itt_denom, oney_itt_denom, twoy_itt_denom, itt_denom, sixm_pp_denom, oney_pp_denom, twoy_pp_denom, pp_denom)

combined_adjusted_results_1 <- adjusted_itt_results_df %>%
  select(prefix, Outcome, term, result) %>%
  pivot_wider(names_from = term, values_from = result, names_sep = "_") # pivot long to wide

combined_adjusted_results_2 <- adjusted_pp_results_df %>%
  select(prefix, Outcome, term, result) %>%
  pivot_wider(names_from = term, values_from = result, names_sep = "_") # pivot long to wide

combined_adjusted_results <- combined_adjusted_results_1 %>%
  left_join(combined_adjusted_results_2, by = c("prefix", "Outcome")) %>%
  select(prefix, Outcome, contains("x"), contains("y")) %>% # re-order
  rename('Time-point' = prefix) %>%
  gt(
    groupname_col = "Outcome" # group by outcome
  ) %>%
  tab_style(
    style = cell_text(style = "italic"), # italicize group names
    locations = cells_row_groups()) %>%
  tab_spanner(label = "Intention-to-treat",  # add spanner for ITT column
              columns = contains(".x"), id = "itt") %>%
  tab_spanner(label = "Per-protocol",  # add spanner for ITT column
              columns = contains(".y"), id = "pp") %>%
  tab_style(style = cell_text(weight = "bold"), # format spanners bold
            locations = cells_column_spanners(spanners = c("itt",  "pp"))) %>%
  tab_style(
    style = cell_borders(
      sides = c("right"),
      color = "lightgrey",
      weight = px(1.5)),
    locations = cells_body(
      columns = c("Aripiprazole vs. Risperidone.x"))) %>% # add line to separate ITT and PP results
  cols_label( # tidy column labels
    'Aripiprazole vs. Risperidone.x' = "Aripiprazole vs. Risperidone",
    'Aripiprazole vs. Risperidone.y' = "Aripiprazole vs. Risperidone",
    'Aripiprazole vs. Olanzapine.x' = "Aripiprazole vs. Olanzapine",
    'Aripiprazole vs. Olanzapine.y' = "Aripiprazole vs. Olanzapine",
    'Aripiprazole vs. Quetiapine.x' = "Aripiprazole vs. Quetiapine",
    'Aripiprazole vs. Quetiapine.y' = "Aripiprazole vs. Quetiapine") %>%
  sub_missing(missing_text = "-") %>% # use - for empty cells
  cols_align(
    align = "center", # center allign columns
    columns = everything()
  ) %>%
  tab_header(title = md("**Cardiometabolic outcomes**")) %>% # overall title # add footnotes
  tab_footnote(footnote = "mmol/L, millimoles per litre; LDL-C, low-density lipoprotein cholesterol; HDL-C, high-density lipoprotein cholesterol; TC:HDL, total cholesterol to high-density lipoprotein; 
               kg, kilogram; mm Hg, millimetres of mercury; mmol/mol, millimoles per mole.") %>%
  tab_footnote(footnote = "Values are adjusted mean difference (95% confidence interval). Missing values were replaced using multiple imputation.") %>%
  tab_footnote(footnote = paste(itt_denominator_footnote),
               locations = cells_column_spanners(spanners = "itt")) %>%
    tab_footnote(footnote = paste(pp_denominator_footnote),
               locations = cells_column_spanners(spanners = "pp"))

# Robust standard errors sensitivity analysis
# ITT
# Add pracid to imputed dataset
robust_data <- complete(tte_imputed, action = "long", include = TRUE) %>%
  filter(diedby1y == 0) %>%
  left_join(tte_prevalent_cohort %>% select(patid, pracid), by = "patid") %>%
  as.mids()

datlist <- miceadds::mids2datlist(robust_data)

# Apply the model to each dataset in datlist
robust_model_pracid <- lapply(datlist, function(data) {
  formula <- as.formula(paste("oney_totalcholesterol", "~", "trt_group", "+", paste(adjustment_set, collapse = " + ")))
  miceadds::lm.cluster(
    data = data, 
    formula = formula, 
    cluster = "pracid")})

# extract parameters and covariance matrix
betas <- lapply(robust_model_pracid, FUN=function(rr){coef(rr) })
vars <- lapply(robust_model_pracid, FUN=function(rr){vcov(rr) })

# Extract pooled results
robust_results <- summary(miceadds::pool_mi(qhat=betas, u=vars))
rm(datlist)

# Per protocol
# Robust standard errors sensitivity analysis
# Add pracid to imputed dataset
robust_data_pp <- complete(tte_imputed, action = "long", include = TRUE) %>%
  filter(active_1y == 1) %>%
  left_join(tte_prevalent_cohort %>% select(patid, pracid), by = "patid") %>%
  as.mids()

datlist <- miceadds::mids2datlist(robust_data_pp)

# Apply the model to each dataset in datlist
robust_model_pracid <- lapply(datlist, function(data) {
  formula <- as.formula(paste("oney_totalcholesterol", "~", "trt_group", "+", paste(adjustment_set, collapse = " + ")))
  miceadds::lm.cluster(
    data = data, 
    formula = formula, 
    cluster = "pracid")})

# extract parameters and covariance matrix
betas <- lapply(robust_model_pracid, FUN=function(rr){coef(rr) })
vars <- lapply(robust_model_pracid, FUN=function(rr){vcov(rr) })

# Extract pooled results
robust_results_pp <- summary(miceadds::pool_mi(qhat=betas, u=vars))
rm(datlist)

# Create table
# Function to process results
format_results <- function(results_df) {
  as.data.frame(results_df) %>%
    rownames_to_column() %>%
    clean_names() %>%
    filter(grepl("trt_group", rowname)) %>%
    select(rowname, estimate = results, conf.low = upper, conf.high = lower) %>%
    mutate(rowname = case_when(grepl("Olanzapine", rowname) ~ "Aripiprazole vs. Olanzapine",
                               grepl("Quetiapine", rowname) ~ "Aripiprazole vs. Quetiapine",
                               grepl("Risperidone", rowname) ~ "Aripiprazole vs. Risperidone",
                               TRUE ~ "Unknown")) %>%
    mutate(across(c(estimate, conf.low, conf.high), ~ round(. * -1, 2))) %>%
    mutate(estimate_ci = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)) %>%
    select(rowname, estimate_ci)}

# Create table
results_combined <- format_results(robust_results) %>%
  left_join(
    format_results(robust_results_pp),
    by = "rowname",
    suffix = c("_itt", "_pp"))

# Create GT table
robust_results_tab <- results_combined %>%
  gt() %>%
  tab_style(style = cell_text(style = "italic"),
            locations = cells_row_groups()) %>%
  tab_spanner(label = md("*Accounting for clustering*"),
              columns = everything(),
              level = 1) %>%
  tab_spanner(label = "Intention-to-treat",
              columns = contains("_itt")) %>%
  tab_spanner(label = "Per-protocol",
              columns = contains("_pp")) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_spanners()) %>%
  cols_label(estimate_ci_itt = "Estimate (95% CI)",
             estimate_ci_pp = "Estimate (95% CI)") %>%
  tab_footnote("Estimated using robust standard errors accounting for clustering by general practice.")

rm(list = ls(pattern = "robust"))

#SUBGROUP ANALYSIS - TOTAL CHOLESTEROL AT ONE YEAR ####

# Subset data based on subgroups ####
# The below only currently works on subgroup defining variables that do not contain missing data (so ethnicity not included)

# Define subgroup defining fields
fields <- c("diag_prev", "diag_prev", "diag_prev", "gender", "gender", "age_atcohortentry_cat_two", "age_atcohortentry_cat_two", 
            "apuse_prior2years", "apuse_prior2years", "cohortentry_timeperiod", "cohortentry_timeperiod", "cohortentry_timeperiod")

# Define subgroup defining values
values <- c("schizophrenia", "bipolar", "other psychosis", "Male", "Female", "Under 50", "50 or over", 
            "Yes", "No", "2005-2009", "2010-2014", "2015+")

# Initialize an empty vector to store the mids names (referred to in later analysis)
mids_names <- c()

# Loop over the fields and values
for (i in 1:length(fields)) {
  field <- fields[i]
  value <- values[i]
  
  # Convert data to mids object
  mids_object <- extract_population(field, value)
  
  # Create name for each mids object using the subgroup defining name, remove any non alphanumeric characters
  mids_name <- paste0("tte_imputed_", field, "_", gsub("[^[:alnum:]]", "", tolower(value)))
  
  # Assign the mids object to the evnrionment, using mids_name
  assign(mids_name, mids_object, envir = .GlobalEnv)
  
  # Add the mids name to a vector
  mids_names <- c(mids_names, mids_name)}

# Run subgroup analyses ####

# Define variables to be excluded per iteration (e.g. diagnosis when looking at diagnostic subgroups)
# Continous variables were not excluded
exclude <- c("diag_prev", "diag_prev", "diag_prev", "gender", "gender", NA, NA, "apuse_prior2years", "apuse_prior2years", NA, NA, NA)

adj_results_subgroup_itt_list <- list()  # Initialize an empty list

# Run the model in each subgroup (excluding the subgroup defining variable of the current subgroup from the adjustment set)
for (i in seq_along(mids_names)) {
  mids_name <- mids_names[i] 
  current_exclude <- exclude[i] # variable to exclude
  
  if (!is.na(current_exclude)) {
    current_adjustment_set <- setdiff(adjustment_set, current_exclude)
  } else {
    current_adjustment_set <- adjustment_set
  }
  
  adj_formula <- parse(text = paste("oney_totalcholesterol ~ trt_group + ", paste(current_adjustment_set, collapse = " + ")))
  adj_models <- with(get(mids_name), lm(eval(adj_formula)))
  adj_results <- pool(adj_models)
  
  adj_results_subgroup_itt_list[[mids_name]] <- adj_results  # Save the pooled results to the list
}

# Create data frame with adjustment formula to check that the right variables were included/excluded
adjustment_df <- data.frame(mids_name = character(), adjustment_set = character(), stringsAsFactors = FALSE)

for (i in seq_along(mids_names)) {
  mids_name <- mids_names[i]
  current_exclude <- exclude[i]
  
  if (!is.na(current_exclude)) {
    current_adjustment_set <- setdiff(adjustment_set, current_exclude)
  } else {
    current_adjustment_set <- adjustment_set
  }
  
  current_results <- data.frame(
    mids_name = mids_name,
    adjustment_set = paste(current_adjustment_set, collapse = ", "),
    stringsAsFactors = FALSE
  )
  
  adjustment_df <- rbind(adjustment_df, current_results)
}

# Put results into DF
subgroup_results_df <- NULL

for (i in seq_along(adj_results_subgroup_itt_list)) {
  mids_name <- names(adj_results_subgroup_itt_list)[i]
  adj_result <- adj_results_subgroup_itt_list[[mids_name]]
  
  # Extract coefficients, p-values, and CIs using tidy()
  tidy_results <- tidy(adj_result, conf.int = TRUE)
  
  # Add mids_name column
  tidy_results$mids_name <- mids_name
  
  # Append to subgroup_results_df
  subgroup_results_df <- rbind(subgroup_results_df, tidy_results)}

subgroup_results_df <- subgroup_results_df %>%
  select(mids_name, term, estimate, conf.low, conf.high, pvalue = p.value) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ . * -1)) %>%
  filter(grepl("trt_group", term)) %>%
  mutate(subgroup = case_when(grepl("diag", mids_name) ~ "Diagnosis",
                              grepl("timeperiod", mids_name) ~ "Time-period",
                              grepl("apuse", mids_name) ~ "Prior non-study AP use",
                              grepl("age", mids_name) ~ "Age",
                              grepl("gender", mids_name) ~ "Gender",
                               TRUE ~ "None")) %>%
  mutate(category = str_extract(mids_name, "(?<=_)[^_]+$")) %>%
  mutate(category = str_to_title(category)) %>%
  mutate(category = case_when(category == "Otherpsychosis" ~ "Other psychosis",
                              category == "Under50" ~ "<50",
                              category == "50orover" ~ "50+",
                              category == "20052009" ~ "2005-2009",
                              category == "20102014" ~ "2010-2014",
                              category == "2015" ~ "2015+",
                              category == "Bipolar" ~ "Bipolar disorder",
                              TRUE ~ category)) %>%
  select(-mids_name)

# Ethnicity
tte_imputed_ethnicity_black <- complete(tte_imputed, action = "long", include = TRUE) %>%
  filter(ethnicity_cat_cprdhes == "black" & diedby1y == 0 & .imp > 0)

tte_imputed_ethnicity_asian <- complete(tte_imputed, action = "long", include = TRUE) %>%
  filter(ethnicity_cat_cprdhes == "asian" & diedby1y == 0 & .imp > 0)

tte_imputed_ethnicity_white <- complete(tte_imputed, action = "long", include = TRUE) %>%
  filter(ethnicity_cat_cprdhes == "white" & diedby1y == 0 & .imp > 0)

tte_imputed_ethnicity_mixedother <- complete(tte_imputed, action = "long", include = TRUE) %>%
  filter((ethnicity_cat_cprdhes == "other" | ethnicity_cat_cprdhes == "mixed") & diedby1y == 0 & .imp > 0) %>% # too few to consider mixed and other separately
  mutate(ethnicity_cat_cprdhes = "Mixed/Other")

labels <- c("black", "asian", "white", "mixedother")
exclude_var <- c("ethnicity_cat_cprdhes")
adjusted_vars <- setdiff(adjustment_set, exclude_var)

# Create empty lists to store the analysis results
ethnicity_subgroup_results_list <- list()

# Loop over each data frame and its corresponding label
for (current_label in labels) {
  # Get the current data frame
  current_df <- get(paste("tte_imputed_ethnicity_", current_label, sep = ""))
  
  # Get unique values of .imp column
  unique_imps <- unique(current_df$.imp)
  
  # Initialize an empty list for storing individual results
  results_list <- list()
  
  # Loop over each unique .imp value
  for (imp_value in unique_imps) {
    # Subset the imputed dataset based on the .imp value
    subset_data <- current_df[current_df$.imp == imp_value, ]
    
    # Perform linear regression
    eth_formula <- as.formula(paste("oney_totalcholesterol ~ trt_group + ", paste(adjusted_vars, collapse = " + ")))
    eth_model <- lm(eval(eth_formula), data = subset_data)
    
    # Save the individual model results to the list
    results_list[[as.character(imp_value)]] <- eth_model}
  
  # Save the results list for the current data frame
  ethnicity_subgroup_results_list[[current_label]] <- results_list}

pooled_results_black <- pool(ethnicity_subgroup_results_list$black)
pooled_results_black_df <- as.data.frame(summary(pooled_results_black, conf.int = TRUE)) %>%
  filter(grepl("trt_group", term)) %>%
  mutate(category = "Black")

pooled_results_asian <- pool(ethnicity_subgroup_results_list$asian)
pooled_results_asian_df <- as.data.frame(summary(pooled_results_asian, conf.int = TRUE)) %>%
  filter(grepl("trt_group", term)) %>%
  mutate(category = "Asian")

pooled_results_white <- pool(ethnicity_subgroup_results_list$white)
pooled_results_white_df <- as.data.frame(summary(pooled_results_white, conf.int = TRUE)) %>%
  filter(grepl("trt_group", term)) %>%
  mutate(category = "White")

pooled_results_mixedother <- pool(ethnicity_subgroup_results_list$mixedother)
pooled_results_mixedother_df <- as.data.frame(summary(pooled_results_mixedother, conf.int = TRUE)) %>%
  filter(grepl("trt_group", term)) %>%
  mutate(category = "Mixed/Other")

pooled_results_ethnicity <- pooled_results_white_df %>%
  bind_rows(pooled_results_black_df) %>%
  bind_rows(pooled_results_asian_df) %>%
  bind_rows(pooled_results_mixedother_df) %>%
  select(category, term, estimate, conf.low = `2.5 %`, conf.high = `97.5 %`, pvalue = p.value) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ . * -1)) %>%
  mutate(subgroup = "Ethnicity")

all_subgroups <- subgroup_results_df %>%
  bind_rows(pooled_results_ethnicity) %>%
  mutate(term = gsub("trt_group", "", term)) %>% # remove trt_group variable name so it just states the group name alone
  select(subgroup, category, everything()) %>%
  arrange(subgroup, category) %>%
  mutate(term = case_when(term == "Olanzapine" ~ "Aripiprazole vs. Olanzapine",
                   term == "Quetiapine" ~ "Aripiprazole vs. Quetiapine",
                   term == "Risperidone" ~ "Aripiprazole vs. Risperidone",
                    TRUE ~ "Unknown"))

subgroup_forestplot <- all_subgroups

all_subgroups <- all_subgroups %>%
  mutate(across(c(estimate, conf.low, conf.high, pvalue), \(x) round(x, 2))) %>%
  mutate(result = paste(format(estimate, nsmall = 2), " (", format(conf.low, nsmall = 2), ", ", format(conf.high, nsmall = 2), ")")) %>% # paste all into one cell
  mutate(result = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", result)) %>% # remove unnecessary spaces in pasted result
  select(-estimate, -conf.low, -conf.high, -pvalue) %>%
  pivot_wider(names_from = term, values_from = result, names_sep = "_") %>% # pivot long to wide
  select(subgroup, category, everything())

all_subgroups_gt <- all_subgroups %>%
  gt(groupname_col = "subgroup") %>%
  tab_style(style = cell_text(style = "italic"), # italicise group names
            locations = cells_row_groups()) %>%
  tab_style(style = cell_text(weight = "bold"), # bold group names
            locations = cells_row_groups()) %>%
  sub_missing(missing_text = "-") %>% # use - for empty cells
  cols_align(align = "center", # center allign columns
             columns = everything()) %>%
  cols_label(category = "Subgroup") %>%
  tab_header(title = md("**Subgroup analysis of total cholesterol at one-year.**")) %>%
  tab_footnote(footnote = "Values are adjusted mean difference (95% confidence interval). Missing values were replaced using multiple imputation.")

# INTERACTIONS ####

# Test for interactions
adjustment_set_str <- paste(adjustment_set, collapse = " + ")

# Diagnosis
diag_formula <- paste("oney_totalcholesterol ~ trt_group * diag_prev +", adjustment_set_str)
diagnosis <- as.data.frame(mi.anova(tte_imputed, diag_formula, type=3))

# Ethnicity
eth_formula <- paste("oney_totalcholesterol ~ trt_group * ethnicity_cat_cprdhes +", adjustment_set_str)
ethnicity <- as.data.frame(mi.anova(tte_imputed, eth_formula, type=3))

# Gender
gender_formula <- paste("oney_totalcholesterol ~ trt_group * gender +", adjustment_set_str)
gender <- as.data.frame(mi.anova(tte_imputed, gender_formula, type=3))

# Age
age_formula <- paste("oney_totalcholesterol ~ trt_group * age_atprevcohortentry +", adjustment_set_str)
age <- as.data.frame(mi.anova(tte_imputed, age_formula, type=3))

# Prior use of non-study AP
apuse_formula <- paste("oney_totalcholesterol ~ trt_group * apuse_prior2years +", adjustment_set_str)
priorapuse <- as.data.frame(mi.anova(tte_imputed, apuse_formula, type=3))

# Time period
timeperiod_formula <- paste("oney_totalcholesterol ~ trt_group * prevcohortentry_year +", adjustment_set_str)
timeperiod <- as.data.frame(mi.anova(tte_imputed, timeperiod_formula, type=3))

interactions <- rbind(diagnosis, ethnicity, gender, age, timeperiod, priorapuse)  %>%
  rownames_to_column(var = "name") %>%
  filter(grepl(":", name)) %>% # filter to just the interaction variables
  rename(interact_pvalue = anova.table.Pr..F.) %>%
  mutate(interact_pvalue = round(interact_pvalue, 3),
         interact_pvalue = ifelse(interact_pvalue < 0.001, "<0.001", format(interact_pvalue, nsmall = 3))) %>% # use "<0.001" for really small p values
  mutate(subgroup = case_when(grepl("diag", name) ~ "Diagnosis",
                              grepl("ethnicity", name) ~ "Ethnicity",
                              grepl("gender", name) ~ "Gender",
                              grepl("age_at", name) ~ "Age",
                              grepl("prevcohortentry_year", name) ~ "Time-period",
                              grepl("apuse_prior2years", name) ~ "Prior non-study AP use",
                              TRUE ~ NA))

remove(diagnosis, ethnicity, gender, age, timeperiod, priorapuse, diag_formula, eth_formula, age_formula, timeperiod_formula, apuse_formula, gender_formula)

# Forest plots

setwd(paste0(path))

# Forest plot of all cardiometabolic outcomes
tiff("tc_oneyear_subgroup_forestplot.tiff", units = "in", width = 11, height = 13, res = 300)

subgroup_forestplot %>%
  left_join(select(interactions, subgroup, interact_pvalue), by = "subgroup") %>%
  group_by(term, subgroup) %>%
  mutate(
    # Only show subgroup and p-value for first row in group
    subgroup_display = if_else(row_number() == 1, subgroup, ""),
    p_value_display = if_else(row_number() == 1, as.character(interact_pvalue), "")) %>%
  ungroup() %>%
  mutate(category = case_when(category == "Other psychosis" ~ "Other psychoses",
                              TRUE ~ category)) %>%
  group_by(term) %>%
  forestplot(clip = c(-0.15, 0.15),
             fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             ci.vertices.height = 0.2,
             boxsize = 0.25,
             xlab = "Adjusted mean difference (mmol/L)",
             legend_args = fpLegend(pos = "top", gp = NULL, title = NULL),
             labeltext = c(subgroup_display, category, p_value_display),
             line.margin = .1,
             mean = estimate,
             lower = conf.high,
             upper = conf.low,
             txt_gp = fpTxtGp(
               label = gpar(cex = 1),
               xlab = gpar(cex = 1.25),
               title = gpar(cex = 1, fontface = "bold"),  # Bold headers
               ticks = gpar(cex = 0.9),
               legend = gpar(cex = 1.25)  # Increase legend text size
             )) %>%
  fp_add_lines("#555555") %>%
  fp_add_header(
    "Subgroup", "Category", 
    expression(bold(paste("ANOVA ", bolditalic("P"))))) %>%
  fp_set_style(
    box = c("#377EB8", "#E69F00", "#D55E00") %>%
      lapply(function(x) gpar(fill = x, col = "#555555")),
    default = gpar(vertices = TRUE)) %>%
  fp_set_zebra_style("#F5F9F9")
dev.off()

# Forest plot of all cardiometabolic outcomes
tiff("cardiometabolic_overall_forestplot.tiff", units = "in", width = 10, height = 13, res = 300)

adjusted_itt_results_df %>%
  group_by(term) %>%
  forestplot(clip = c(-1.75, 0.80),
             fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             ci.vertices.height = 0.2,
             boxsize = 0.25,
             xlab = "Adjusted mean difference (mmol/L)",
             legend_args = fpLegend(pos = "top", gp = NULL, title = NULL),
             labeltext = c(Outcome, prefix),
             line.margin = .1,
             mean = estimate,
             lower = ci_low,
             upper = ci_high) %>%
  fp_add_lines("#555555") %>%
  fp_add_header("Outcome") %>%
  fp_set_style(
    box = c("#377EB8", "#E69F00", "#D55E00") %>%
      lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))

dev.off()

setwd(paste0(path))

# ITT descriptives ####
# Table of means and SDs at each time point

#Means and SD based on the imputed dataset
cardiometabolic_descriptives_itt <- tte_imputed %>%
  mice::complete("long") %>%
  select(trt_group, .imp, all_of(baseline_names), all_of(study_outcomes)) %>%
  group_by(trt_group, .imp) %>%
  summarise_all(list(mean = ~ mean(., na.rm = TRUE),
                     var = ~ var(., na.rm = TRUE))) %>%
  ungroup()

# Calculate pooled means by averaging the means across imputations
pooled_means <- cardiometabolic_descriptives_itt %>%
  group_by(trt_group) %>%
  summarise_at(vars(ends_with("_mean")), ~ mean(.))

# Calculate within-imputation variance (W) by averaging the variances across imputations
within_variance <- cardiometabolic_descriptives_itt %>%
  group_by(trt_group) %>%
  summarise_at(
    vars(ends_with("_var")),
    ~ mean(., na.rm = TRUE)
  ) %>%
  rename_with(~ str_replace(., "_var", "_within_var"))

# Calculate between-imputation variance (B) by calculating the variance of the means across imputations
between_variance <- cardiometabolic_descriptives_itt %>%
  group_by(trt_group) %>%
  summarise_at(
    vars(ends_with("_mean")),
    ~ var(., na.rm = TRUE)
  ) %>%
  rename_with(~ str_replace(., "_mean", "_between_var"))

# Calculate Rubin's Rules SDs using within and between variances
pooled_sd <- within_variance %>%
  inner_join(between_variance, by = "trt_group") %>%
  mutate(across(ends_with("_within_var"), 
                ~ sqrt(. + (1 + 1 / max(cardiometabolic_descriptives_itt$.imp)) * get(str_replace(cur_column(), "_within_var", "_between_var"))),
                .names = "{str_replace(.col, '_within_var', '_sd')}")
  ) %>%
  select(trt_group, ends_with("_sd"))

# Convert means from wide to long format
means_long <- pooled_means %>%
  melt(id.vars = "trt_group", variable.name = "variable", value.name = "mean") %>%
  dcast(variable ~ trt_group, value.var = "mean")

# Convert SDs from wide to long format
sd_long <- pooled_sd %>%
  melt(id.vars = "trt_group", variable.name = "variable", value.name = "sd") %>%
  dcast(variable ~ trt_group, value.var = "sd")

# Rename treatment group columns with _mean and _sd suffixes
names(means_long)[-1] <- paste0(names(means_long)[-1], "_mean")
names(sd_long)[-1] <- paste0(names(sd_long)[-1], "_sd")

# Remove suffixes from the variable column in both data frames
means_long$variable <- gsub("_mean", "", means_long$variable)
sd_long$variable <- gsub("_sd", "", sd_long$variable)

# Merge means and SDs based on the variable column
merged_df <- merge(means_long, sd_long, by = "variable") %>%
  select(sort(names(.))) %>%
  select(variable, everything()) %>%
  mutate_at(vars(2:9), ~round(., 2)) %>%
  mutate(Aripiprazole = paste0(format(Aripiprazole_mean, nsmall = 2), " (", format(Aripiprazole_sd, nsmall = 2), ")"),
         Olanzapine = paste0(format(Olanzapine_mean, nsmall = 2), " (", format(Olanzapine_sd, nsmall = 2), ")"),
         Quetiapine = paste0(format(Quetiapine_mean, nsmall = 2), " (", format(Quetiapine_sd, nsmall = 2), ")"),
         Risperidone = paste0(format(Risperidone_mean, nsmall = 2), " (", format(Risperidone_sd, nsmall = 2), ")")) %>%
  select(variable, 10:13) %>%
  mutate(Aripiprazole = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", Aripiprazole)) %>%  # remove unnecessary spaces in pasted result
  mutate(Olanzapine = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", Olanzapine)) %>%
  mutate(Quetiapine = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", Quetiapine)) %>%
  mutate(Risperidone = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", Risperidone)) %>%
  separate(variable, into = c("prefix", "outcome"), sep = "_", remove = FALSE) %>%
  left_join(order, by = "outcome") %>%
  mutate(prefix = if_else(prefix == "baseline", "Baseline",
                          if_else(prefix == "sixm", "6m",
                                  if_else(prefix == "oney", "1y",
                                          if_else(prefix == "twoy", "2y", "")))),
         time_order = if_else(prefix == "6m", 1,
                              if_else(prefix == "1y", 2,
                                      if_else(prefix == "2y", 3, 0)))) %>%
  arrange(order, time_order) %>%
  select(prefix, Outcome, Aripiprazole, Olanzapine, Quetiapine, Risperidone) %>%
  rename('Time-point' = prefix) 

# Convert to gt table
cardiometabolic_values <- merged_df %>%
  gt(groupname_col = "Outcome") %>%
  tab_style(style = cell_text(style = "italic"), # italicise group names
            locations = cells_row_groups()) %>%
  sub_missing(missing_text = "-") %>% # use - for empty cells
  cols_align(
    align = "center", # center allign columns
    columns = everything()) %>%
  tab_header(title = md("**Cardiometabolic values (ITT)**")) %>% # overall title # add footnotes
  tab_footnote(footnote = "Values are mean (SD). Missing values were replaced using multiple imputation.")

# PP descriptives
cardiometabolic_descriptives_pp <- tte_imputed %>%
  mice::complete("long") %>%
  select(patid, .imp, trt_group, active_6m, active_1y, active_2y, all_of(baseline_names), all_of(study_outcomes)) %>%
  mutate(across(contains("sixm"), ~ if_else(active_6m != 1, NA_real_, .))) %>% # set values to NA if not active at time point
  mutate(across(contains("oney"), ~ if_else(active_1y != 1, NA_real_, .))) %>%
  mutate(across(contains("twoy"), ~ if_else(active_2y != 1, NA_real_, .))) %>%
  select(-contains("active"), -patid) %>%
  group_by(trt_group, .imp) %>%
  summarise_all(list(mean = ~ mean(., na.rm = TRUE),
                     var = ~ var(., na.rm = TRUE))) %>%
  ungroup()

# Calculate pooled means by averaging the means across imputations
pooled_means <- cardiometabolic_descriptives_pp %>%
  group_by(trt_group) %>%
  summarise_at(vars(ends_with("_mean")), ~ mean(.))

# Calculate within-imputation variance (W) by averaging the variances across imputations
within_variance <- cardiometabolic_descriptives_pp %>%
  group_by(trt_group) %>%
  summarise_at(
    vars(ends_with("_var")),
    ~ mean(., na.rm = TRUE)
  ) %>%
  rename_with(~ str_replace(., "_var", "_within_var"))

# Calculate between-imputation variance (B) by calculating the variance of the means across imputations
between_variance <- cardiometabolic_descriptives_pp %>%
  group_by(trt_group) %>%
  summarise_at(
    vars(ends_with("_mean")),
    ~ var(., na.rm = TRUE)
  ) %>%
  rename_with(~ str_replace(., "_mean", "_between_var"))

# Calculate Rubin's Rules SDs using within and between variances
pooled_sd <- within_variance %>%
  inner_join(between_variance, by = "trt_group") %>%
  mutate(across(ends_with("_within_var"), 
                ~ sqrt(. + (1 + 1 / max(cardiometabolic_descriptives_pp$.imp)) * get(str_replace(cur_column(), "_within_var", "_between_var"))),
                .names = "{str_replace(.col, '_within_var', '_sd')}")
  ) %>%
  select(trt_group, ends_with("_sd"))

# Convert means from wide to long format
means_long <- pooled_means %>%
  melt(id.vars = "trt_group", variable.name = "variable", value.name = "mean") %>%
  dcast(variable ~ trt_group, value.var = "mean")

# Convert SDs from wide to long format
sd_long <- pooled_sd %>%
  melt(id.vars = "trt_group", variable.name = "variable", value.name = "sd") %>%
  dcast(variable ~ trt_group, value.var = "sd")

# Rename treatment group columns with _mean and _sd suffixes
names(means_long)[-1] <- paste0(names(means_long)[-1], "_mean")
names(sd_long)[-1] <- paste0(names(sd_long)[-1], "_sd")

# Remove suffixes from the variable column in both data frames
means_long$variable <- gsub("_mean", "", means_long$variable)
sd_long$variable <- gsub("_sd", "", sd_long$variable)

# Merge means and SDs based on the variable column
merged_df_pp <- merge(means_long, sd_long, by = "variable") %>%
  select(sort(names(.))) %>%
  select(variable, everything()) %>%
  mutate_at(vars(2:9), ~round(., 2)) %>%
  mutate(Aripiprazole = paste0(format(Aripiprazole_mean, nsmall = 2), " (", format(Aripiprazole_sd, nsmall = 2), ")"),
         Olanzapine = paste0(format(Olanzapine_mean, nsmall = 2), " (", format(Olanzapine_sd, nsmall = 2), ")"),
         Quetiapine = paste0(format(Quetiapine_mean, nsmall = 2), " (", format(Quetiapine_sd, nsmall = 2), ")"),
         Risperidone = paste0(format(Risperidone_mean, nsmall = 2), " (", format(Risperidone_sd, nsmall = 2), ")")) %>%
  select(variable, 10:13) %>%
  mutate(Aripiprazole = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", Aripiprazole)) %>%  # remove unnecessary spaces in pasted result
  mutate(Olanzapine = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", Olanzapine)) %>%
  mutate(Quetiapine = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", Quetiapine)) %>%
  mutate(Risperidone = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", Risperidone)) %>%
  separate(variable, into = c("prefix", "outcome"), sep = "_", remove = FALSE) %>%
  left_join(order, by = "outcome") %>%
  mutate(prefix = if_else(prefix == "baseline", "Baseline",
                          if_else(prefix == "sixm", "6m",
                                  if_else(prefix == "oney", "1y",
                                          if_else(prefix == "twoy", "2y", "")))),
         time_order = if_else(prefix == "6m", 1,
                              if_else(prefix == "1y", 2,
                                      if_else(prefix == "2y", 3, 0)))) %>%
  arrange(order, time_order) %>%
  select(prefix, Outcome, Aripiprazole, Olanzapine, Quetiapine, Risperidone) %>%
  rename('Time-point' = prefix) 

# Convert to gt table
cardiometabolic_values_pp <- merged_df_pp %>%
  gt(groupname_col = "Outcome") %>%
  tab_style(style = cell_text(style = "italic"), # italicise group names
            locations = cells_row_groups()) %>%
  sub_missing(missing_text = "-") %>% # use - for empty cells
  cols_align(
    align = "center", # center allign columns
    columns = everything()) %>%
  tab_header(title = md("**Cardiometabolic values (PP)**")) %>% # overall title # add footnotes
  tab_footnote(footnote = "Values are mean (SD). Missing values were replaced using multiple imputation.")

# Combined ITT and PP descriptives
cardiometabolic_values_ittandpp <- merged_df %>%
  left_join(merged_df_pp, by =c("Outcome", "Time-point")) %>%
  gt(
    groupname_col = "Outcome"
  ) %>%
  tab_style(
    style = cell_text(style = "italic"), # italicise group names
    locations = cells_row_groups()) %>%
  tab_spanner(label = "Intention-to-treat",  # add spanner for ITT SHRs column
              columns = c(Aripiprazole.x, Olanzapine.x, Quetiapine.x, Risperidone.x), id = "itt") %>%
  tab_spanner(label = "Per-protocol",  # add spanner for PP SHRs column
              columns = c(Aripiprazole.y, Olanzapine.y, Quetiapine.y, Risperidone.y), id = "pp") %>%
  tab_style(style = cell_text(weight = "bold"), # format spanners bold
            locations = cells_column_spanners(spanners = c("itt",  "pp"))) %>%
  tab_style(
    style = cell_borders(
      sides = c("right"),
      color = "lightgrey",
      weight = px(1.5)),
    locations = cells_body(
      columns = c(Risperidone.x))) %>%
      cols_label( # tidy column labels
    Aripiprazole.x = "Aripiprazole",
    Olanzapine.x = "Olanzapine",
    Quetiapine.x = "Quetiapine",
    Risperidone.x = "Risperidone",
    Aripiprazole.y = "Aripiprazole",
      Olanzapine.y = "Olanzapine",
      Quetiapine.y = "Quetiapine",
      Risperidone.y = "Risperidone") %>%
  sub_missing(missing_text = "-") %>% # use - for empty cells
  cols_align(
    align = "center", # center allign columns
    columns = everything()
  ) %>%
  tab_header(title = md("**Cardiometabolic values**")) %>% # overall title # add footnotes
  tab_footnote(footnote = "mmol/L, millimoles per litre; LDL-C, low-density lipoprotein cholesterol; HDL-C, high-density lipoprotein cholesterol; TC:HDL, total cholesterol to high-density lipoprotein; 
               kg, kilogram; mm Hg, millimetres of mercury; mmol/mol, millimoles per mole.") %>%
  tab_footnote(footnote = "Values are mean (SD). Missing values were replaced using multiple imputation.") %>%
  tab_footnote(footnote = "Values were only included if patient had not switched or discontinued by the time-point.",
               locations = cells_column_spanners(spanners = "pp"))

# IPCW sensitivity analysis ####

# Load files ####
tte_imputed_long_ipcw <- tte_imputed_long %>%
  left_join(select(tte_prevalent_cohort, patid, baseline_totalcholesterol_flag, sixm_totalcholesterol_flag, oney_totalcholesterol_flag), by = "patid")

# Define adjustment sets
adjustment_set <- c("age_atprevcohortentry", "gender", "ethnicity_cat_cprdhes", "pat_2019imd_quintile", "diag_prev", "region", "prevcohortentry_year", "gpconsults_last6m", "smoking_status_cat", 
                    "prioralcoholabuse", "priorsubstanceabuse", "priordyslipidaemia", "priordiabetes", "priorhypertension", "priorcerebrovasculardisease", "priormyocardialinfarction", "priorrenaldisease", "priorliverdisease", 
                    "apuse_prior2years", "lipiddrugs_prior2years", "hypertensiondrugs_prior2years", "antidiabetics_prior2years", "antidepressant_prior2years", "moodstab_prior2years",
                    "baseline_totalcholesterol", "baseline_ldl", "baseline_hdl", "baseline_triglycerides", "baseline_systolicbp", "baseline_diastolicbp", "baseline_glucose", "baseline_hba1c", "baseline_weightkg", "baseline_bmi_cat")

adjustment_set_6m <- c("age_atprevcohortentry", "gender", "ethnicity_cat_cprdhes", "pat_2019imd_quintile", "diag_prev", "region", "prevcohortentry_year", "gpconsults_last6m", "smoking_status_cat", 
                       "prioralcoholabuse", "priorsubstanceabuse", "priordyslipidaemia", "priordiabetes", "priorhypertension", "priorcerebrovasculardisease", "priormyocardialinfarction", "priorrenaldisease", "priorliverdisease", 
                       "apuse_prior2years", "lipiddrugs_prior2years", "hypertensiondrugs_prior2years", "antidiabetics_prior2years", "antidepressant_prior2years", "moodstab_prior2years",
                       "baseline_totalcholesterol", "baseline_ldl", "baseline_hdl", "baseline_triglycerides", "baseline_systolicbp", "baseline_diastolicbp", "baseline_glucose", "baseline_hba1c", "baseline_weightkg", "baseline_bmi_cat",
                       "baseline_totalcholesterol_flag")

adjustment_set_1y <- c("age_atprevcohortentry", "gender", "ethnicity_cat_cprdhes", "pat_2019imd_quintile", "diag_prev", "region", "prevcohortentry_year", "gpconsults_last6m", "smoking_status_cat", 
                       "prioralcoholabuse", "priorsubstanceabuse", "priordyslipidaemia", "priordiabetes", "priorhypertension", "priorcerebrovasculardisease", "priormyocardialinfarction", "priorrenaldisease", "priorliverdisease", 
                       "apuse_prior2years", "lipiddrugs_prior2years", "hypertensiondrugs_prior2years", "antidiabetics_prior2years", "antidepressant_prior2years", "moodstab_prior2years",
                       "baseline_totalcholesterol", "baseline_ldl", "baseline_hdl", "baseline_triglycerides", "baseline_systolicbp", "baseline_diastolicbp", "baseline_glucose", "baseline_hba1c", "baseline_weightkg", "baseline_bmi_cat",
                       "sixm_totalcholesterol", "sixm_ldl", "sixm_hdl", "sixm_triglycerides", "sixm_systolicbp", "sixm_diastolicbp", "sixm_glucose", "sixm_hba1c", "sixm_weightkg",
                       "baseline_totalcholesterol", "sixm_totalcholesterol_flag")

adjustment_set_2y <- c("age_atprevcohortentry", "gender", "ethnicity_cat_cprdhes", "pat_2019imd_quintile", "diag_prev", "region", "prevcohortentry_year", "gpconsults_last6m", "smoking_status_cat", 
                       "prioralcoholabuse", "priorsubstanceabuse", "priordyslipidaemia", "priordiabetes", "priorhypertension", "priorcerebrovasculardisease", "priormyocardialinfarction", "priorrenaldisease", "priorliverdisease", 
                       "apuse_prior2years", "lipiddrugs_prior2years", "hypertensiondrugs_prior2years", "antidiabetics_prior2years", "antidepressant_prior2years", "moodstab_prior2years",
                       "baseline_totalcholesterol", "baseline_ldl", "baseline_hdl", "baseline_triglycerides", "baseline_systolicbp", "baseline_diastolicbp", "baseline_glucose", "baseline_hba1c", "baseline_weightkg", "baseline_bmi_cat",
                       "sixm_totalcholesterol", "sixm_ldl", "sixm_hdl", "sixm_triglycerides", "sixm_systolicbp", "sixm_diastolicbp", "sixm_glucose", "sixm_hba1c", "sixm_weightkg",
                       "oney_totalcholesterol", "oney_ldl", "oney_hdl", "oney_triglycerides", "oney_systolicbp", "oney_diastolicbp", "oney_glucose", "oney_hba1c", "oney_weightkg",
                       "baseline_totalcholesterol", "sixm_totalcholesterol_flag", "oney_totalcholesterol_flag")

calculate_censor_weights <- function(data, adjustment_vars, time_point) {
  results_list <- lapply(1:25, function(imp_value) {
    # Define the outcome variable and died variable based on the time point
    outcome_var <- paste0("active_", time_point)
    died_var <- paste0("diedby", time_point)
    
    # Filter data for the current imputation value
    filtered_data <- data %>%
      filter(.imp == imp_value, !!sym(died_var) == 0) %>%
      select(patid, trt_group, !!sym(outcome_var), all_of(adjustment_vars))
    
    # Censor models
    censor_model <- glm(as.formula(paste(outcome_var, "== 1 ~ trt_group +", paste(adjustment_vars, collapse = " + "))),
                        data = filtered_data,
                        family = binomial(link = "logit"))
    
    censor_model_simple <- glm(as.formula(paste(outcome_var, "== 0 ~ trt_group")),
                               data = filtered_data,
                               family = binomial(link = "logit"))
    
    # Add predictions and stabilized weights
    filtered_data <- filtered_data %>%
      mutate(censor_weight_denom = predict(censor_model, type = "response"),
             censor_weight_num = predict(censor_model_simple, type = "response"),
             stabilized_censor_weight = (1 - censor_weight_num) / (1 - censor_weight_denom),
             trimmed_stabilized_censor_weight = ifelse(stabilized_censor_weight > quantile(stabilized_censor_weight, 0.995), 
                                                       quantile(stabilized_censor_weight, 0.995), stabilized_censor_weight),
             .imp = imp_value) %>%
      select(patid, .imp, trt_group, !!sym(outcome_var), censor_weight_denom, censor_weight_num, stabilized_censor_weight, 
             trimmed_stabilized_censor_weight)    
    return(filtered_data)})
  
  # Combine results into a single data frame
  final_results <- bind_rows(results_list)
  return(final_results)}

results_6m <- calculate_censor_weights(tte_imputed_long_ipcw, adjustment_set_6m, "6m")
results_1y <- calculate_censor_weights(tte_imputed_long_ipcw, adjustment_set_1y, "1y")
results_2y <- calculate_censor_weights(tte_imputed_long_ipcw, adjustment_set_2y, "2y")

# Review the distributions of weights 

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

# Calculate group statistics for both weight types
group_stats_censor_1y <- calculate_group_stats(results_1y, "trimmed_stabilized_censor_weight")

# LINEAR REGRESSION, with covariate adjustment and inverse probability of censoring weights

# Add weights to mids object
tte_imputed_pp_6m_cen <- tte_imputed_long_ipcw %>%
  left_join(results_6m, by = c("patid", "trt_group", ".imp", "active_6m")) %>%
  filter(active_6m == 1) %>%
  as.mids()

tte_imputed_pp_1y_cen <- tte_imputed_long_ipcw %>%
  left_join(results_1y, by = c("patid", "trt_group", ".imp", "active_1y")) %>%
  filter(active_1y == 1) %>%
  as.mids()

tte_imputed_pp_2y_cen <- tte_imputed_long_ipcw %>%
  left_join(results_2y, by = c("patid", "trt_group", ".imp", "active_2y")) %>%
  filter(active_2y == 1) %>%
  as.mids()

rm(results_6m, results_1y, results_2y, tte_imputed_long_ipcw)

# Outcomes

# Outcome names
outcome_names <- c("totalcholesterol", "ldl", "hdl", "triglycerides", "tchdlratio", "weightkg", "systolicbp", "diastolicbp", "hba1c", "glucose")
prefixes <- c("sixm_", "oney_", "twoy_")
study_outcomes <- lapply(prefixes, function(prefix) paste0(prefix, outcome_names))
study_outcomes <- unlist(study_outcomes)

sixm_outcomes <- study_outcomes[grepl("sixm", study_outcomes)]
oney_outcomes <- study_outcomes[grepl("oney", study_outcomes)]
twoy_outcomes <- study_outcomes[grepl("twoy", study_outcomes)]

# Define cardiometabolic baseline variables
baseline_prefix <- c("baseline_")
baseline_names <- lapply(baseline_prefix, function(prefix) paste0(baseline_prefix, outcome_names))
baseline_names <- unlist(baseline_names)

#Initialise lists to store results
results_pp_6m <- list()
results_pp_1y <- list()
results_pp_2y <- list()

# One year
for (outcome in oney_outcomes) {
  
  # fit adjusted model for each imputed dataset
  adj_models <- with(tte_imputed_pp_1y_cen, 
                     lm(as.formula(paste(outcome, "~ trt_group + ", paste(adjustment_set, collapse = " + "))),
                        weights = trimmed_stabilized_censor_weight))
  
  # extract results, combine across imputed datasets, add to list
  adj_results <- pool(adj_models)
  results_pp_1y[[outcome]] <- adj_results}

rm(tte_imputed_pp_1y_cen, adj_models, adj_results)
gc()

results_pp_1y_df <- results_pp_1y %>% 
  map_dfr(~ tidy(.x, conf.int = TRUE), .id = "outcome")

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
packageVersion("car") # 3.1.2
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
