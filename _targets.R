#### Workflow Orchestration by Targets =======================================================================
library(targets)
library(tidyverse)

#### Source Function ####
tar_source("R/clean_data_raw.R")

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
  )
)
