#### Workflow Orchestration by Targets =======================================================================
library(targets)
library(tidyverse)

#### Source Function ####
tar_source("R/clean_data_raw.R")
tar_source("dag/dgp_anthropometry.R")
tar_source("R/anthropometry_models.R")

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
  tar_target(
    name = gen_bmi_total,
    command = dgp_total_bmi(n = 200)
  ),
  tar_target(
    name = recovery_bmi_total,
    command = bmi_total_model(gen_bmi_total)
  ),
  tar_target(
    name = bmi_total_effect,
    command = bmi_total_model(data_clean)
  )
)
