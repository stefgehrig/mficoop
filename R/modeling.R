# load libraries
library(tidyverse)    # v2.0.0
library(brms)         # v2.19.0
library(ordbetareg)   # v0.7.2
library(cmdstanr)     # v0.5.3
library(countrycode)  # v1.5.0
library(cluster)      # v2.1.4
library(RColorBrewer) # v1.1-3

# load data and fcts
df <- read_csv2("data/data_mix_integr.csv")
source("R/functions.R")

# placeholder for cultural cluster variable (created right before model fitting)
df$ling_cluster <- "cluster_NA"

# set reference level for institutional type to a common category
df$inst_type <- relevel(factor(df$inst_type), ref = "NGO")

# mcmc sampler config
options(mc.cores = parallel::detectCores(logical = TRUE))
cores <- parallel::detectCores(logical = TRUE)-1
warmup <- 1e3
iter   <- 5e3
chains <- 2

####################
#### main model ####
####################
# define variables
outcome_var  <- "par30"
intact_vars  <- c("cooperation", "wdi_agricul")
fin_inst_env <- c("wdi_privcredit", 
                  "wgi_reg",
                  "wgi_law",
                  "wgi_sta",
                  "wgi_acc",
                  "wgi_eff",
                  "wgi_cor")
macroeco_env <- c("wdi_gdp", 
                  "wdi_gdp_sq", 
                  "wdi_gdpgrowth",  
                  "wdi_remittances", 
                  "wdi_laborforce", 
                  "wdi_manufac", 
                  "wdi_agricul",
                  "swiid_gini") 
nominal_vars <- c("year", 
                  "inst_type")
orgnz_charct <- c("age",
                  "log_avg_loan", 
                  "log_numborrs", 
                  "log_ass_loan")
random_itcpt <- c("ling_cluster", 
                  "country.name")
# run model
run_model(
  data                = df,
  outcome_var         = outcome_var,
  intact_vars         = intact_vars,
  fin_inst_env        = fin_inst_env,
  macroeco_env        = macroeco_env,
  orgnz_charct        = orgnz_charct,
  nominal_vars        = nominal_vars,
  random_itcpt        = random_itcpt,
  priorset            = 1,
  modelname           = "main",
  warmup              = warmup,
  iter                = iter,
  chains              = chains, 
  withincores         = floor(cores / chains),
  cores               = cores,
  seed                = 12345,
  get_prior_sims_too  = TRUE
)


###############################
# varying precision parameter #
###############################
run_model(
  data                = df,
  outcome_var         = outcome_var,
  intact_vars         = intact_vars,
  fin_inst_env        = fin_inst_env,
  macroeco_env        = macroeco_env,
  orgnz_charct        = orgnz_charct,
  nominal_vars        = nominal_vars,
  random_itcpt        = random_itcpt,
  model_phi_group     = TRUE,
  priorset            = 1,
  modelname           = "phi",
  warmup              = warmup,
  iter                = iter,
  chains              = chains, 
  withincores         = floor(cores / chains),
  cores               = cores,
  seed                = 12345
)

################################
#### mfi varying intercepts ####
################################
run_model(
  data                = df,
  outcome_var         = outcome_var,
  intact_vars         = intact_vars,
  fin_inst_env        = fin_inst_env,
  macroeco_env        = macroeco_env,
  orgnz_charct        = orgnz_charct,
  nominal_vars        = nominal_vars,
  random_itcpt        = c("ling_cluster", "country.name", "mfiname"),
  priorset            = 1,
  modelname           = "mfi",
  warmup              = warmup,
  iter                = iter,
  chains              = chains, 
  withincores         = floor(cores / chains),
  cores               = cores,
  seed                = 12345
)

###############################
#### no varying intercepts ####
###############################
run_model(
  data                = df,
  outcome_var         = outcome_var,
  intact_vars         = intact_vars,
  fin_inst_env        = fin_inst_env,
  macroeco_env        = macroeco_env,
  orgnz_charct        = orgnz_charct,
  nominal_vars        = nominal_vars,
  random_itcpt        = c(),
  priorset            = 1,
  modelname           = "noitcpts",
  warmup              = warmup,
  iter                = iter,
  chains              = chains,
  withincores         = floor(cores / chains),
  cores               = cores,
  seed                = 12345
)

################################
#### country varying slopes ####
################################
run_model(
  data                = df,
  outcome_var         = outcome_var,
  intact_vars         = intact_vars,
  fin_inst_env        = fin_inst_env,
  macroeco_env        = macroeco_env,
  orgnz_charct        = orgnz_charct,
  nominal_vars        = nominal_vars,
  random_itcpt        = random_itcpt,
  priorset            = 1,
  modelname           = "ranslopes",
  ranslopes           = TRUE,
  warmup              = warmup,
  iter                = iter,
  chains              = 2, 
  withincores         = 7,
  cores               = cores,
  seed                = 12345
)

#######################
#### exclude Haiti ####
#######################
run_model(
  data                = df %>% filter(country.name != "Haiti"),
  outcome_var         = outcome_var,
  intact_vars         = intact_vars,
  fin_inst_env        = fin_inst_env,
  macroeco_env        = macroeco_env,
  orgnz_charct        = orgnz_charct,
  nominal_vars        = nominal_vars,
  random_itcpt        = random_itcpt,
  priorset            = 1,
  modelname           = "nohaiti",
  warmup              = warmup,
  iter                = iter,
  chains              = chains, 
  withincores         = floor(cores / chains),
  cores               = cores,
  seed                = 12345
)

#######################
#### par90 outcome ####
#######################
# run modelling
run_model(
  data                = df,
  outcome_var         = "par90",
  intact_vars         = intact_vars,
  fin_inst_env        = fin_inst_env,
  macroeco_env        = macroeco_env,
  orgnz_charct        = orgnz_charct,
  nominal_vars        = nominal_vars,
  random_itcpt        = random_itcpt,
  priorset            = 1,
  modelname           = "par90",
  warmup              = warmup,
  iter                = iter,
  chains              = chains, 
  withincores         = floor(cores / chains),
  cores               = cores,
  seed                = 12345
)

####################################
#### alternative adjustment set ####
####################################
run_model(
  data                = df,
  outcome_var         = outcome_var,
  intact_vars         = c("cooperation", "wdi_rural"),
  fin_inst_env        = c("wdi_atmdens",     
                          "wdi_bankdens",    
                          "efi_score",       
                          "db_getcredit",    
                          "db_enfcontract",  
                          "db_registprop" ), 
  macroeco_env        = c("hdi",             
                          "hdi_sq",          
                          "wdi_gdpgrowth",   
                          "wdi_inflation",   
                          "wdi_fdi",         
                          "wdi_rural",       
                          "fgi_ineq"),       
  orgnz_charct        = orgnz_charct,
  nominal_vars        = nominal_vars,
  random_itcpt        = random_itcpt,
  priorset            = 1,
  modelname           = "otheradj",
  warmup              = warmup,
  iter                = iter,
  chains              = chains, 
  withincores         = floor(cores / chains),
  cores               = cores,
  seed                = 12345
)

#################################################
#### all country-level variable interactions ####
#################################################
run_model(
  data                = df,
  outcome_var         = outcome_var,
  intact_vars         = c("cooperation",
                          fin_inst_env,
                          macroeco_env),
  fin_inst_env        = fin_inst_env,
  macroeco_env        = macroeco_env,
  orgnz_charct        = orgnz_charct,
  nominal_vars        = nominal_vars,
  random_itcpt        = random_itcpt,
  priorset            = 1, 
  modelname           = "allintact",
  warmup              = warmup,
  iter                = iter,
  chains              = chains, 
  withincores         = floor(cores / chains),
  cores               = cores,
  seed                = 12345
)

#######################
#### weaker priors ####
#######################
run_model(
  data                = df,
  outcome_var         = outcome_var,
  intact_vars         = intact_vars,
  fin_inst_env        = fin_inst_env,
  macroeco_env        = macroeco_env,
  orgnz_charct        = orgnz_charct,
  nominal_vars        = nominal_vars,
  random_itcpt        = random_itcpt,
  priorset            = 2,
  modelname           = "weaker",
  warmup              = warmup,
  iter                = iter,
  chains              = chains, 
  withincores         = floor(cores / chains),
  cores               = cores,
  seed                = 12345
)
