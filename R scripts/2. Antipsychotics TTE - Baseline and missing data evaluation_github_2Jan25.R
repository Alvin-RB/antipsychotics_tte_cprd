# -----------------------
# ANTIPSYCHOTICS TTE STUDY
# -----------------------

# BASELINE CHARACTERISTICS AND MISSING DATA EVALUATION

# Last updated: 15/05/2024

#Clear memory
rm(list = ls())

#Packages
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
library(ggplot2)
library(gtsummary)
library(forcats)
library(gridExtra)
library(lubridate)
library(readxl)
library(reshape2)
library(gt)
library(webshot)
library(ftExtra)
library(tidylog)
library(scales)
library(janitor)
library(RColorBrewer)

# Set file path
path <- anonymised
# Set working directory
setwd(paste0(path))

#Load data
load(file = "tte_prevalent_cohort.Rdata")

# Calculate the number of patients per practice
patients_perpractice_mean <- (n_distinct(tte_prevalent_cohort$patid)) / (n_distinct(tte_prevalent_cohort$pracid))

# Accrual plot

# Calculate number of patients per month/year per group
accrual <- tte_prevalent_cohort %>%
  mutate(month_year = lubridate::floor_date(cohortentrydate, unit = "month")) %>%
  group_by(month_year, trt_group) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  arrange(month_year) %>%
  group_by(trt_group) %>%
  mutate(cumulative_count = cumsum(count))

accrual_plot <- ggplot(accrual, aes(x = month_year, y = cumulative_count, color = trt_group)) +
  geom_line(linewidth = .85) +
  labs(x = "Year", y = "Cumulative patient accrual", color = "Study group") +
  scale_x_date(date_labels = "%Y", date_breaks = "2 year", limits = as.Date(c("2005-01-01", "2017-12-31"))) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  guides(colour = guide_legend(override.aes = list(size=12))) +
  theme(axis.line = element_line(),
        legend.text = element_text(size = 12),  # Adjust legend text size
        axis.title = element_text(size = 12))  # Adjust axis label size

# Baseline characteristics tables
tte_prevalent_cohort <- tte_prevalent_cohort %>%
  mutate(ethnicity_cat_cprdhes = str_to_title(ethnicity_cat_cprdhes),
         diag_prev = case_when(diag_prev == "bipolar" ~ "Bipolar disorder",
                               diag_prev == "other psychosis" ~ "Other non-organic psychoses",
                               diag_prev == "schizophrenia" ~ "Schizophrenia",
                               TRUE ~ diag_prev),
         smoking_status_cat = case_when(smoking_status_cat == "Never" ~ "Never smoked",
                                        smoking_status_cat == "Ex" ~ "Ex-smoker",
                                        smoking_status_cat == "Current" ~ "Current smoker",
                                        TRUE ~ smoking_status_cat))

# Specify variable types for baseline tables
type_list <- list(
  age_atprevcohortentry ~ "continuous",
  age_atcohortentry_cat_two ~ "categorical",
  gender ~ "categorical",
  age_atfirstdiag ~ "continuous",
  years_diagtoindex ~ "continuous",
  cohortentry_timeperiod ~ "categorical",
  patprac_2019imd_quintile ~ "categorical",
  ethnicity_cat_cprdhes ~ "categorical",
  diag_prev ~ "categorical",
  baselinedose_olanz ~ "continuous",
  apuse_prior2years ~ "dichotomous",
  priorhosp_any ~ "dichotomous",
  antidepressant_prior2years ~ "dichotomous",
  moodstab_prior2years ~ "dichotomous",
  lipiddrugs_prior2years ~ "dichotomous",
  antidiabetics_prior2years ~ "dichotomous",
  hypertensiondrugs_prior2years ~ "dichotomous",
  baseline_totalcholesterol ~ "continuous",
  baseline_bmi ~ "continuous",
  baseline_bmi_cat ~ "categorical",
  smoking_status_cat ~ "categorical",
  priordyslipidaemia ~ "dichotomous",
  priorcerebrovasculardisease ~ "dichotomous", 
  priormyocardialinfarction ~ "dichotomous", 
  priordiabetes ~ "dichotomous",
  priorhypertension ~ "dichotomous",
  prioralcoholabuse ~ "dichotomous", 
  priorsubstanceabuse ~ "dichotomous",
  priorliverdisease ~ "dichotomous",
  priorrenaldisease ~ "dichotomous",
  gpconsults_last6m ~ "continuous")

# Specify variable labels for baseline tables
label_list <- list(
  age_atprevcohortentry = "Age at baseline (years), median (IQR)",
  age_atcohortentry_cat_two = "Age at baseline (years), n (%)",
  age_atfirstdiag = "Age at diagnosis (years), median (IQR)",
  years_diagtoindex ~ "Years from first diagnosis to index date, median (IQR)",
  cohortentry_timeperiod = "Index date time period, n (%)",
  patprac_2019imd_quintile = "IMD quintile, n (%)",
  gender = "Sex, n (%)",
  ethnicity_cat_cprdhes = "Ethnicity, n (%)",
  diag_prev = "Diagnosis, n (%)",
  apuse_prior2years = "Prescribed an antipsychotic in last 2y, n (%)",
  priorhosp_any = "Psychiatric hospitalisation in last 2y, n (%)",
  baselinedose_olanz = "Baseline dose (Olanzapine equivalent, mg), mean (SD)",
  antidepressant_prior2years = "Prescribed an antidepressant",
  moodstab_prior2years = "Prescribed a mood stabiliser",
  lipiddrugs_prior2years = "Prescribed a lipid regulating medication",
  antidiabetics_prior2years = "Prescribed an antidiabetic",
  hypertensiondrugs_prior2years = "Prescribed an antihypertensive",
  baseline_totalcholesterol ~ "Total cholesterol (mmol/L), mean (SD)",
  baseline_bmi = "BMI (mg/m2), mean (SD)",
  baseline_bmi_cat = "BMI category, n (%)",
  smoking_status_cat ~ "Smoking status, n (%)",
  priordyslipidaemia = "Prior dyslipidaemia",
  priorcerebrovasculardisease = "Prior cerebrovascular disease", 
  priormyocardialinfarction = "Prior myocardial infarction", 
  priordiabetes = "Prior diabetes",
  priorhypertension = "Prior hypertension",
  prioralcoholabuse = "Prior alcohol misuse", 
  priorsubstanceabuse = "Prior substance misuse",
  priorliverdisease = "Prior liver disease",
  priorrenaldisease = "Prior renal disease",
  gpconsults_last6m = "Number of primary care consults in last 6m, median (IQR)")

# Specify continuous variable statistics for baseline tables
statistic_list <- list(
  age_atprevcohortentry ~ "{median} ({p25}, {p75})",
  age_atfirstdiag ~ "{median} ({p25}, {p75})",
  years_diagtoindex ~ "{median} ({p25}, {p75})",
  baseline_bmi ~ "{mean} ({sd})",
  baseline_totalcholesterol ~ "{mean} ({sd})",
  gpconsults_last6m ~ "{median} ({p25}, {p75})",
  baselinedose_olanz ~ "{mean} ({sd})")

# Specify continuous variable digits for baseline tables
digits_list <- list(
  all_categorical() ~ c(0, 1),
  age_atprevcohortentry ~ c(0, 0),
  age_atfirstdiag ~ c(0, 0),
  years_diagtoindex ~ c(1, 1),
  baseline_bmi ~ c(1, 1),
  baseline_totalcholesterol ~ c(2, 2),
  gpconsults_last6m ~ c(0, 0),
  baselinedose_olanz ~ c(1, 1))

# Overall baseline characteristics, stratified by group
characteristics_baseline <- tte_prevalent_cohort %>%
  select(trt_group, age_atprevcohortentry, gender, ethnicity_cat_cprdhes, diag_prev, age_atfirstdiag, years_diagtoindex, apuse_prior2years, priorhosp_any,
         cohortentry_timeperiod, patprac_2019imd_quintile, priordyslipidaemia, priordiabetes, priorhypertension, 
         priorcerebrovasculardisease, priormyocardialinfarction, prioralcoholabuse, priorsubstanceabuse, priorliverdisease, priorrenaldisease,
         antidepressant_prior2years, moodstab_prior2years, lipiddrugs_prior2years, antidiabetics_prior2years, hypertensiondrugs_prior2years,
         baseline_totalcholesterol, baseline_bmi, baseline_bmi_cat, smoking_status_cat, gpconsults_last6m, baselinedose_olanz)

characteristics_baseline <- as.data.frame(characteristics_baseline %>%
                         tbl_summary(
                           by = trt_group,
                           statistic = statistic_list,
                           digits = digits_list,
                           label = label_list,
                           sort = everything() ~ "alphanumeric")) %>%
  rename_all(~gsub("\\*", "", .)) %>%
  mutate(Characteristic = str_replace_all(Characteristic, "_", ""),
         row_num = row_number()) %>%
  add_row(row_num = 31.5, Characteristic = "Comorbidities, n (%)") %>%
  add_row(row_num = 40.5, Characteristic = "Concomitant medications, n (%)") %>%
  arrange(row_num) %>%
  mutate(row_num = row_number(),
         Characteristic = if_else(row_num >= 32 & row_num <= 47,
                                  str_to_sentence(gsub("Prior|Prescribed an|Prescribed a", "", Characteristic)),
                                  Characteristic)) %>%
  gt() %>%
  tab_style(style = cell_text(weight = "bold"), # format spanners bold
            locations = cells_column_labels(everything())) %>%
  tab_style(style = cell_text(weight = "bold"), # format spanners bold
            locations = cells_body(columns = 1, rows = c(1, 2, 5, 12, 16, 17, 18, 19, 21, 25, 32, 42, 48, 50, 52, 58, 63, 64))) %>%
  tab_footnote("IMD, index of multiple deprivation; BMI, body mass index; mg, miligrams") %>%
  tab_footnote("IMD defined accroding to the patient's postcode or where this was missing, the GP practice postcode.", locations = cells_body(columns = 1, rows = 25)) %>%
  tab_footnote("Comorbidities determined using the patient's available full medical history.", locations = cells_body(columns = 1, rows = 32)) %>%
  tab_footnote("Concomitant medications defined according to prescriptions in the prior two years.", locations = cells_body(columns = 1, rows = 42)) %>%
  tab_footnote("The most frequently prescribed were citalopram, mirtazapine, sertraline, amitriptyline, fluoxetine and venlafaxine.", locations = cells_body(columns = 1, rows = 43)) %>%
  tab_footnote("The most frequently prescribed were lithium and sodium valproate.", locations = cells_body(columns = 1, rows = 44)) %>%
  tab_footnote("The most frequently prescribed were simvastatin and atorvastatin.", locations = cells_body(columns = 1, rows = 45)) %>%
  tab_footnote("The most frequently prescribed were metformin, gliclazide and insulin.", locations = cells_body(columns = 1, rows = 46)) %>%
  tab_footnote("The most frequently prescribed were amlodipine, ramipril and bendroflumethiazide.", locations = cells_body(columns = 1, rows = 47)) %>%
  tab_footnote("Calculated according to the Defined Daily Dose method and expressed as an olanzapine equivalent dose.", locations = cells_body(columns = 1, rows = 64)) %>%
  
  cols_hide(row_num) %>%
  sub_missing(missing_text = "") %>% # use - for empty cells
  cols_align(
    align = "center", # center allign columns
    columns = 2:5)

# Overall cohort characteristics (not stratified)
overall_characteristics_baseline <- tte_prevalent_cohort %>%
  select(trt_group, age_atprevcohortentry, age_atcohortentry_cat_two, gender, ethnicity_cat_cprdhes, diag_prev, age_atfirstdiag, apuse_prior2years, priorhosp_any,
         cohortentry_timeperiod, patprac_2019imd_quintile, priordyslipidaemia, priordiabetes, priorhypertension, 
         priorcerebrovasculardisease, priormyocardialinfarction, prioralcoholabuse, priorsubstanceabuse, priorliverdisease, priorrenaldisease,
         antidepressant_prior2years, moodstab_prior2years, lipiddrugs_prior2years, antidiabetics_prior2years, hypertensiondrugs_prior2years,
         baseline_totalcholesterol, baseline_bmi, baseline_bmi_cat, smoking_status_cat, gpconsults_last6m, baselinedose_olanz) %>% 
  tbl_summary(
    type = type_list,
    statistic = statistic_list,
    digits = digits_list,
    label = label_list,
    sort = everything() ~ "alphanumeric") %>%
  bold_labels() %>%
  modify_caption("**Table 1. Sample characterisitcs**")

# Baseline characteristics according to whether or not total cholesterol was observed at one year, in those alive at 1y
characteristics_missingprimary <- tte_prevalent_cohort %>%
  filter(diedby1y != 1) %>%
  select(trt_group, age_atprevcohortentry, age_atcohortentry_cat_two, gender, ethnicity_cat_cprdhes, diag_prev, age_atfirstdiag, apuse_prior2years, priorhosp_any,
         cohortentry_timeperiod, patprac_2019imd_quintile, priordyslipidaemia, priordiabetes, priorhypertension, 
         priorcerebrovasculardisease, priormyocardialinfarction, prioralcoholabuse, priorsubstanceabuse, priorliverdisease, priorrenaldisease,
         antidepressant_prior2years, moodstab_prior2years, lipiddrugs_prior2years, antidiabetics_prior2years, hypertensiondrugs_prior2years,
         baseline_totalcholesterol, baseline_bmi, baseline_bmi_cat, smoking_status_cat, gpconsults_last6m, baselinedose_olanz, oney_totalcholesterol) %>% 
  mutate(missing_primary = case_when(is.na(oney_totalcholesterol) ~ "Missing",
                                     TRUE ~ "Observed")) %>%
  select(-oney_totalcholesterol) %>%
  tbl_summary(
    by = missing_primary,
    type = type_list,
    statistic = statistic_list,
    digits = digits_list,
    label = label_list,
    sort = everything() ~ "alphanumeric") %>%
  bold_labels() %>%
  modify_caption("**Table 1. Characteristics of patients with missing and observed primary outcome data**")

# Cardiometabolic data availability summary
# Create a table for each time point, calculate percentages with patients alive at the time point as the denominator
baseline_availability <- tte_prevalent_cohort %>%
  select(trt_group, baseline_totalcholesterol_flag, baseline_hdl_flag, baseline_ldl_flag, baseline_triglycerides_flag, baseline_tchdlratio_flag,
         baseline_weightkg_flag, baseline_systolicbp_flag, baseline_diastolicbp_flag, baseline_hba1c_flag, baseline_glucose_flag) %>%
  tbl_summary(by = trt_group,
              label = list(
                baseline_totalcholesterol_flag ~ "Total cholesterol at baseline",
                baseline_hdl_flag ~ "HDL-C at baseline",
                baseline_ldl_flag ~ "LDL-C at baseline",
                baseline_triglycerides_flag ~ "Triglycerides at baseline",
                baseline_tchdlratio_flag ~ "TC:HDL ratio at baseline",
                baseline_weightkg_flag ~ "Weight at baseline",
                baseline_systolicbp_flag ~ "Systolic blood pressure at baseline",
                baseline_diastolicbp_flag ~ "Diastolic blood pressure at baseline",
                baseline_hba1c_flag ~ "HbA1c at baseline",
                baseline_glucose_flag ~ "Glucose at baseline"))

sixm_availability <- tte_prevalent_cohort %>%
  filter(diedby6m == 0) %>%
  select(trt_group, sixm_totalcholesterol_flag, sixm_hdl_flag, sixm_ldl_flag, sixm_triglycerides_flag, sixm_tchdlratio_flag,
         sixm_weightkg_flag, sixm_systolicbp_flag, sixm_diastolicbp_flag, sixm_hba1c_flag, sixm_glucose_flag) %>%
  tbl_summary(by = trt_group,
              label = list(
                sixm_totalcholesterol_flag ~ "Total cholesterol at 6m",
                sixm_hdl_flag ~ "HDL-C at 6m",
                sixm_ldl_flag ~ "LDL-C at 6m",
                sixm_triglycerides_flag ~ "Triglycerides at 6m",
                sixm_tchdlratio_flag ~ "TC:HDL ratio at 6m",
                sixm_weightkg_flag ~ "Weight at 6m",
                sixm_systolicbp_flag ~ "Systolic blood pressure at 6m",
                sixm_diastolicbp_flag ~ "Diastolic blood pressure at 6m",
                sixm_glucose_flag ~ "Glucose at 6m", 
                sixm_hba1c_flag ~ "HbA1c at 6m"))

oney_availability <- tte_prevalent_cohort %>%
  filter(diedby1y == 0) %>%
  select(trt_group, oney_totalcholesterol_flag, oney_hdl_flag, oney_ldl_flag, oney_triglycerides_flag, oney_tchdlratio_flag,
         oney_weightkg_flag, oney_systolicbp_flag, oney_diastolicbp_flag, oney_hba1c_flag, oney_glucose_flag) %>%
  tbl_summary(by = trt_group,
              label = list(
                oney_totalcholesterol_flag ~ "Total cholesterol at 1y",
                oney_hdl_flag ~ "HDL-C at 1y",
                oney_ldl_flag ~ "LDL-C at 1y",
                oney_tchdlratio_flag ~ "TC:HDL ratio at 1y",
                oney_triglycerides_flag ~ "Triglycerides at 1y",
                oney_weightkg_flag ~ "Weight at 1y",
                oney_systolicbp_flag ~ "Systolic blood pressure at 1y",
                oney_diastolicbp_flag ~ "Diastolic blood pressure at 1y",
                oney_hba1c_flag ~ "HbA1c at 1y",
                oney_glucose_flag ~ "Glucose at 1y"))

twoy_availability <- tte_prevalent_cohort %>%
  filter(diedby6m == 0) %>%
  select(trt_group, twoy_totalcholesterol_flag, twoy_hdl_flag, twoy_ldl_flag, twoy_triglycerides_flag, twoy_tchdlratio_flag,
         twoy_weightkg_flag, twoy_systolicbp_flag, twoy_diastolicbp_flag, twoy_hba1c_flag, twoy_glucose_flag) %>%
  tbl_summary(by = trt_group,
              label = list(
                twoy_totalcholesterol_flag ~ "Total cholesterol at 2y",
                twoy_hdl_flag ~ "HDL-C at 2y",
                twoy_ldl_flag ~ "LDL-C at 2y",
                twoy_triglycerides_flag ~ "Triglycerides at 2y",
                twoy_tchdlratio_flag ~ "TC:HDL ratio at 2y",
                twoy_weightkg_flag ~ "Weight at 2y",
                twoy_systolicbp_flag ~ "Systolic blood pressure at 2y",
                twoy_diastolicbp_flag ~ "Diastolic blood pressure at 2y",
                twoy_glucose_flag ~ "Glucose at 2y", 
                twoy_hba1c_flag ~ "HbA1c at 2y"))

any_availability <- tte_prevalent_cohort %>%
  filter(diedby6m == 0) %>%
  select(trt_group, totalcholesterol_bin, hdl_bin, ldl_bin, triglycerides_bin, tchdlratio_bin, weightkg_bin, systolicbp_bin, diastolicbp_bin, hba1c_bin, glucose_bin) %>%
  tbl_summary(by = trt_group,
              label = list(
                totalcholesterol_bin ~ "Total cholesterol - any time-point",
                hdl_bin ~ "HDL-C - any time-point",
                ldl_bin ~ "LDL-C - any time-point",
                triglycerides_bin ~ "Triglycerides - any time-point",
                tchdlratio_bin ~ "TC:HDL ratio - any time-point",
                weightkg_bin ~ "Weight - any time-point",
                systolicbp_bin ~ "Systolic blood pressure - any time-point",
                diastolicbp_bin ~ "Diastolic blood pressure - any time-point",
                hba1c_bin ~ "HbA1c - any time-point",
                glucose_bin ~ " Glucose - any time-point"))

data_availability <- tbl_stack(tbls = list(baseline_availability, sixm_availability, oney_availability, twoy_availability, any_availability))

dataavailability <- as.data.frame(data_availability) %>%
  clean_names() %>%
  rename_at(vars(contains("_")), ~gsub("_.*$", "", .)) %>% # remove the group numbers from the field name
  mutate(outcome_num = case_when(grepl("(?i)total cholesterol", characteristic) ~ 1,
                                 grepl("(?i)HDL-C", characteristic) ~ 2,
                                 grepl("(?i)LDL-C", characteristic) ~ 3,
                                 grepl("(?i)Triglycerides", characteristic) ~ 4,
                                 grepl("(?i)TC:HDL ratio", characteristic) ~ 5,
                                 grepl("(?i)Weight", characteristic) ~ 6,
                                 grepl("(?i)Systolic", characteristic) ~ 7,
                                 grepl("(?i)Diastolic", characteristic) ~ 8,
                                 grepl("(?i)Hba1c", characteristic) ~ 9,
                                 grepl("(?i)Glucose", characteristic) ~ 10,
                                 TRUE ~ 0)) %>%
  arrange(outcome_num) %>%
  mutate(outcome_type = case_when(grepl("(?i)total cholesterol", characteristic) ~ "Total cholesterol",
                                  grepl("(?i)HDL-C", characteristic) ~ "HDL-C",
                                  grepl("(?i)LDL-C", characteristic) ~ "LDL-C",
                                  grepl("(?i)Triglycerides", characteristic) ~ "Triglycerides",
                                  grepl("(?i)TC:HDL ratio", characteristic) ~ "TC:HDL ratio",
                                  grepl("(?i)Weight", characteristic) ~ "Weight",
                                  grepl("(?i)Systolic", characteristic) ~ "Systolic blood pressure",
                                  grepl("(?i)Diastolic", characteristic) ~ "Diastolic blood pressure",
                                  grepl("(?i)Hba1c", characteristic) ~ "Hba1c",
                                  grepl("(?i)Glucose", characteristic) ~ "Glucose",
                                  TRUE ~ "unk")) %>%
  mutate(outcome_time = case_when(grepl("(?i)baseline", characteristic) ~ "At baseline",
                                  grepl("(?i)6m", characteristic) ~ "At 6m",
                                  grepl("(?i)1y", characteristic) ~ "At 1y",
                                  grepl("(?i)2y", characteristic) ~ "At 2y",
                                  grepl("(?i)time-point", characteristic) ~ "Any time-point")) %>%
  select(outcome_type, outcome_time, aripiprazole, olanzapine, quetiapine, risperidone) %>%
  gt(
    groupname_col = "outcome_type"
  ) %>%
  tab_style(
    style = cell_text(style = "italic"), # italicize group names
    locations = cells_row_groups()) %>%
  cols_label( # tidy column labels
    aripiprazole = "Aripiprazole",
    olanzapine = "Olanzapine",
    quetiapine = "Quetiapine",
    risperidone = "Risperidone",
    outcome_time = "Outcome, n (%)") %>%
  sub_missing(missing_text = "-") %>% # use - for empty cells
  cols_align(
    align = "center", # center align columns
    columns = everything()) %>%
  tab_footnote(
    footnote = "Values were considered baseline if measured on the index date or within the two years prior (using the most recent value, where multiple values were available).
               Time-windows for cardiometabolic outcomes included a period before and after the study outcome time-point (e.g., at 6 ± 3 months, 12 ± 3 months, and 24 ± 6 months), with the measurement closest to the time-point used.
    Percentages are of patients alive at each time-point. Patients are counted in the denominator for the 'any time-point' measure if alive at six months.")

#Most frequent co-medications

# Mood stabilisers
load("SMIMoodStabilisers_C.Rdata")
moodstabilisers <- tte_prevalent_cohort %>%
  select(patid, cohortentrydate, trt_group) %>%
  inner_join(select(SMImoodstabilisers_C, patid, MoodStabiliser, issuedate),  by = "patid") %>%
  filter(issuedate >= cohortentrydate - days(731) & issuedate <= cohortentrydate) %>%
  select(patid, MoodStabiliser) %>%
  distinct() %>%
  select(-patid) %>%
  tbl_summary(sort = all_categorical() ~ "frequency")

# Lipid-regulating drugs
load("SMIlipiddrugs_C.Rdata")
lipiddrugs <- tte_prevalent_cohort %>%
  select(patid, cohortentrydate, trt_group) %>%
  inner_join(select(SMIlipiddrugs_C, patid, LipidRegulatingMed, issuedate), by = "patid") %>%
  filter(issuedate >= cohortentrydate - days(731) & issuedate <= cohortentrydate) %>%
  select(patid, LipidRegulatingMed) %>%
  distinct() %>%
  select(-patid) %>%
  tbl_summary(sort = all_categorical() ~ "frequency")

# Antidiabetics
load("SMIantidiabetics_C.Rdata")
antidiabetics <- tte_prevalent_cohort %>%
  select(patid, cohortentrydate, trt_group) %>%
  inner_join(SMIantidiabetics_C, by = "patid") %>%
  select(patid, cohortentrydate, trt_group, Antidiabetic, issuedate) %>%
  filter(issuedate >= cohortentrydate - days(731) & issuedate <= cohortentrydate) %>%
  select(patid, Antidiabetic) %>%
  distinct() %>%
  select(-patid) %>%
  tbl_summary(sort = all_categorical() ~ "frequency")

# Antihypertensives
load("SMIhypertensionmeds_C.Rdata")
antihypertensive <- tte_prevalent_cohort %>%
  select(patid, cohortentrydate, trt_group) %>%
  inner_join(SMIhypertensionmeds_C, by = "patid") %>%
  select(patid, cohortentrydate, trt_group, antihypertensive, issuedate) %>%
  filter(issuedate >= cohortentrydate - days(731) & issuedate <= cohortentrydate) %>%
  select(patid, antihypertensive) %>%
  distinct() %>%
  select(-patid) %>%
  tbl_summary(sort = all_categorical() ~ "frequency")

# Antidepressants
load("SMIantidepressants.Rdata")
antidepressants <- tte_prevalent_cohort %>%
  select(patid, cohortentrydate, trt_group) %>%
  inner_join(SMIantidepressants, by = "patid") %>%
  select(patid, cohortentrydate, trt_group, Antidepressant, issuedate) %>%
  filter(issuedate >= cohortentrydate - days(731) & issuedate <= cohortentrydate) %>%
  select(patid, Antidepressant) %>%
  distinct() %>%
  select(-patid) %>%
  tbl_summary(sort = all_categorical() ~ "frequency")

# Antipsychotics
load("antipsychotics_combined.Rdata")
# issuedate < cohortentrydate used as only wanting to count the non-study APs
antipsychotics_combined <- antipsychotics_combined %>%
  select(patid,AP, issuedate) %>%
  filter(AP != "Prochlorperazine")

aps <- tte_prevalent_cohort %>%
  select(patid, cohortentrydate, trt_group) %>%
  inner_join(antipsychotics_combined, by = "patid") %>%
  select(patid, cohortentrydate, trt_group, AP, issuedate) %>%
  filter(issuedate >= cohortentrydate - days(731) & issuedate < cohortentrydate) %>%
  select(patid, AP) %>%
  distinct() %>%
  select(-patid) %>%
  tbl_summary(sort = all_categorical() ~ "frequency")

# Was there differential rates of starting statins after index date?
load("SMIlipiddrugs_C.Rdata")

statins_afterindex <- tte_prevalent_cohort %>%
  select(patid, cohortentrydate) %>%
  left_join(SMIlipiddrugs_C, by = "patid") %>%
  distinct() %>%
  group_by(patid) %>%
  summarise(
    prior = ifelse(any(issuedate >= cohortentrydate - days(731) & issuedate <= cohortentrydate), 1, 0), # received statins in prior 2y
    oneyear = ifelse(any(issuedate > cohortentrydate & issuedate <= cohortentrydate + days(365)), 1, 0), # received within 1y of index
    twoyears = ifelse(any(issuedate > cohortentrydate + days(365) & issuedate <= cohortentrydate + days(730)), 1, 0)) %>% # received between 1-2y after index
  ungroup() %>%
  mutate(across(2:4, ~ ifelse(is.na(.), 0, .)),
         within2y = oneyear + twoyears,
         within2y = case_when(within2y == 2 ~ 1, # received within 2y of index
                              TRUE ~ within2y))

statinspatients <- tte_prevalent_cohort %>%
  select(patid, trt_group, lipiddrugs_prior2years, diedby1y, diedby2y) %>%
  left_join(statins_afterindex, by = "patid")

statinspatients_tab <- statinspatients %>%
  filter(diedby1y == 0) %>%
#  filter(lipiddrugs_prior2years == "No") %>% # restrict table to those alive at 1y, to match main ITT outcome
  select(-patid, -diedby1y, -diedby2y) %>%
  tbl_summary(by = trt_group,
              label = list(
                prior = "Statins at baseline",
                oneyear = "Prescribed a lipid regulating medication within 1 year of index date",
                twoyears = "Prescribed a lipid regulating medication between 1-2 years after index date",
                within2y = "Prescribed a lipid regulating medication within 2 years of index date"),
              digits = list(
                all_categorical() ~ c(0, 1)))

# Package versions
packageVersion("dplyr") # 1.1.4
packageVersion("tidyr") # 1.3.1
packageVersion("stringr") # 1.5.1
packageVersion("readr") # 2.1.5
packageVersion("purrr") # 1.0.2
packageVersion("ggplot2") # 3.5.1
packageVersion("gtsummary") # 2.0.1
packageVersion("forcats") # 1.0.0
packageVersion("gridExtra") # 2.3
packageVersion("lubridate") # 1.9.3
packageVersion("readxl") # 1.4.3
packageVersion("reshape2") # 1.4.4
packageVersion("gt") # 0.11.0
packageVersion("webshot") # 0.5.5
packageVersion("ftExtra") # 0.6.4
packageVersion("tidylog") # 1.1.0
packageVersion("scales") # 1.3.0
packageVersion("janitor") # 2.2.0
packageVersion("RColorBrewer") # 1.1.3
