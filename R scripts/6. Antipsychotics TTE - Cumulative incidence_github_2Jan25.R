# -----------------------
# ANTIPSYCHOTICS TTE STUDY
# -----------------------
# Last run: 14/05/2024

# CUMULATIVE INCIDENCE OF EFFECTIVENESS OUTCOMES ####

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

# Set file path
path <- anonymised
# Set working directory
setwd(paste0(path))

#Load files
load(file = "tte_prevalent_cohort.Rdata")
load(file = "tte_imputed_long.Rdata")
tte_imputed <- as.mids(tte_imputed_long)

# PSYCHIATRIC HOSPITALISATION ####

# FORMULAS

# Competing risks - ITT
psychhosp_itt_crr_formula <- as.formula("Surv(time = hospdeath_crr_fuptime, event = hospdeath_crr) ~ trt_group") # ITT formula

# Competing risks - PP
psychhosp_pp_crr_formula <- as.formula("Surv(time = hospdeath_pp_crr_fuptime, event = hospdeath_pp_crr) ~ trt_group") # PP formula

# COMPETING RISKS CUMULATIVE INCIDENCE

# ITT
psychhosp_cuminc <- cuminc(psychhosp_itt_crr_formula, data = tte_prevalent_cohort) # used for table

# PP
psychhosp_pp_cuminc <- cuminc(psychhosp_pp_crr_formula, data = tte_prevalent_cohort) # used for table

# MORTALITY ####

mortality_itt_formula <- as.formula("Surv(time = mortality_fuptime_days_cen, event = diedby2y) ~ trt_group") # define formula

# KAPLAN-MEIER ANALYSIS
mortality_itt_fit <- surv_fit(mortality_itt_formula, data = tte_prevalent_cohort) #Fit KM survival model
mortality_itt_logrank <- survdiff(mortality_itt_formula, data = tte_prevalent_cohort)
mortality_itt_logrank_pval <- na.omit(as.numeric(sub(".*p= ([0-9.e+-]+).*", "\\1", capture.output(print(mortality_itt_logrank))))) # Extract p value
mortality_itt_logrank_pval <- sprintf("%.8f", mortality_itt_logrank_pval[!is.na(mortality_itt_logrank_pval)]) # convert from scientific notation

mortality_pp_fit <- surv_fit(mortality_pp_formula, data = tte_prevalent_cohort)
mortality_pp_logrank <- survdiff(mortality_pp_formula, data = tte_prevalent_cohort)
mortality_pp_logrank_pval <- na.omit(as.numeric(sub(".*p= ([0-9.e+-]+).*", "\\1", capture.output(print(mortality_pp_logrank))))) # Extract p value
mortality_pp_logrank_pval <- sprintf("%.8f", mortality_pp_logrank_pval[!is.na(mortality_pp_logrank_pval)]) # convert from scientific notation

# TIME TO DISCONTINUATION ####

# COMPETING RISKS REGRESSION ####
# Primary competing risks analyses codes people who died prior to discontinuation as 2 for the event variables and uses the mortality fup time for them.
# This ensures those people are appropriately treated as having a competing risk of death (precluding hospitalisation), whereas the others remain the same (i.e. discontinued without dying, and had neither events)
discont_crr_formula <- as.formula("Surv(time = discdeath_crr_fuptime, event = discdeath_crr) ~ trt_group") # define formula

# Cumulative incidence
discont_cuminc <- cuminc(discont_crr_formula, data = tte_prevalent_cohort) # used for table

# TABLES/FIGURES ####

# Functions

# Generate survival tables
generate_KM_survival_table <- function(model) {
  tbl_survfit(model, times = c(180, 365, 730), label_header = "**{time} days**", 
              estimate_fun = function(x) style_percent(x, symbol = TRUE, digits = 1, columns = "survival"), reverse = TRUE)}

generate_crr_survival_table <- function(model) {
  model %>%
    tbl_cuminc(times = c(180, 365, 730), label_header = "**{time} days**",
               estimate_fun = function(x) style_percent(x, symbol = FALSE, digits = 1))}

# Convert gtsummary regression tables into a data frame
convert_gtsummary_to_df <- function(gtsummarytable, resulttype) {
  as.data.frame(gtsummarytable) %>%
    rename(!!resulttype := `**HR**`, treatment = `**Characteristic**`)}

# Hospitalisation

psychhosp_crr_incidencetable <- generate_crr_survival_table(psychhosp_cuminc)
psychhosp_crr_pp_incidencetable <- generate_crr_survival_table(psychhosp_pp_cuminc)
psychhosp_crr_combined_incidencetable <- tbl_stack(tbls = list(psychhosp_crr_incidencetable, psychhosp_crr_pp_incidencetable),
                                                   group_header = c("Intention-to-treat", "Per-protocol"))

# Plots

#Survival tables
mortality_itt_survtable <- generate_KM_survival_table(mortality_itt_fit)
mortality_pp_survtable <- generate_KM_survival_table(mortality_pp_fit)
mortality_combined_survtable <- tbl_stack(tbls = list(mortality_itt_survtable, mortality_pp_survtable),
          group_header = c("Intention-to-treat", "Per-protocol")) 

# Discontinuation

#Survival table
discont_cuminc_tab <- generate_crr_survival_table(discont_cuminc)

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
