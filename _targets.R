#### Workflow Orchestration by Targets =======================================================================
library(targets)
library(tarchetypes)
library(tidyverse)

#### Source Function ####
tar_source("R/clean_data_raw.R")
tar_source("dag/dgp_anthropometry.R")
tar_source("R/anthropometry_models.R")
tar_source("R/sbc.R")

# Pipeline
list(
  # Load the raw data
  tar_target(
    name = data_raw,
    command = read.csv("data/raw/Data_raw.csv")
  ),
  # Clean the data
  tar_target(
    name = data_clean,
    command = clean_data_raw(data_raw)
  ),
  
  ## Exploratory Data Analysis =================================================
  
  # Prepare the data for EDA 
  tar_target(
    name = data_eda,
    command = prepare_data_eda(data_clean)
  ),
  # Explore the Antropometrics 
  tar_quarto(
    name = antropometrics_report,
    path = "analysis/01_explore_anthropometry.qmd",
    quiet = TRUE
  ),
  
  ## Total effect of PIP on post-intervention BMI ==============================
  
  ## Prepare the data 
  tar_target(
    name = bmi_total_effect_data,
    command = prepare_data_bmi_total(data_clean)
  ),
  # Generative Model for Param Recovery
  tar_target(
    name = gen_bmi_total,
    command = dgp_total_bmi(n = 200) 
  ),
  # Recovery Model 
  tar_target(
    name = recovery_bmi_total,
    command = bmi_total_model(gen_bmi_total),
    memory = "transient",
    garbage_collection = TRUE
  ),
  # BMI Total Effect Model 
  tar_target(
    name = bmi_total_effect,
    command = bmi_total_model(bmi_total_effect_data),
    memory = "transient",
    garbage_collection = TRUE
  ),
  # Prior BMI Total effect model
  tar_target(
    name = bmi_total_effect_prior,
    command = update(bmi_total_effect, sample_prior = "only"),
    memory = "transient",
    garbage_collection = TRUE
  ),
  # Simulation Based Calibration 
  tar_target(
    name = bmi_total_effect_sbc,
    command = sbc_bmi_total(bmi_total_effect_data),
    memory = "transient",
    garbage_collection = TRUE
  ),
  # Final Report on the BMI Total Effect 
  tar_quarto(
    name = bmi_total_effect_report,
    path = "reports/treatment_bmi_effect.qmd",
    quiet = TRUE
  ),
  #### Anthropometric outcomes, joint multivariate model ####
  tar_target(
    name = anthropometric_data,
    command = prepare_anthropometric_data(data_clean)
  ),
  tar_target(
    name = gen_anthropometric,
    command = dgp_anthropometric(n = 200)
  ),
  tar_target(
    name = anthropometric_recovery,
    command = anthropometric_model(gen_anthropometric)
  ),
  tar_target(
    name = anthropometric_effect,
    command = anthropometric_model(anthropometric_data)
  ),
  tar_target(
    name = anthropometric_priors,
    command = update(anthropometric_effect, sample_prior = "only")
  ),
  tar_target(
    name = anthropometric_sbc,
    command = sbc_anthropometric_model(anthropometric_data)
  ),
  tar_quarto(
    name = anthropometric_report,
    path = "reports/anthropometric_report.qmd"
  ),
  ####  Mediated effect through behavioral latent ####
  ## CFA analysis on the behaviour indicators
  tar_target(
    name = gen_cfa_behaviour,
    command = dgp_cfa_behaviour(n = 200)
  )
)
