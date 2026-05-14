#### Workflow Orchestration by Targets =======================================================================
library(targets)
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

  ## Total effect of PIP on post-intervention BMI ==============================
  
  # Generative Model for Param Recovery
  tar_target(
    name = gen_bmi_total,
    command = dgp_total_bmi(n = 200) 
  ),
  # Recovery Model 
  tar_target(
    name = recovery_bmi_total,
    command = bmi_total_model(gen_bmi_total) 
  ),
  # BMI Total Effect Model 
  tar_target(
    name = bmi_total_effect,
    command = bmi_total_model(data_clean) 
  ),
  # Prior BMI Total effect model
  tar_target(
    name = bmi_total_effect_prior,
    command = update(bmi_total_effect, sample_prior = "only") 
  ),
  # Simulation Based Calibration 
  tar_target(
    name = bmi_total_effect_sbc,
    command = sbc_bmi_total(data_clean)
  )
)
