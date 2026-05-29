#### Workflow Orchestration by Targets =======================================================================
library(targets)
library(tarchetypes)
library(tidyverse)
library(cmdstanr)

#### Source Function ####
tar_source("R/clean_data_raw.R")
tar_source("dag/dgp_anthropometry.R")
tar_source("R/anthropometry_models.R")
tar_source("R/sbc.R")
tar_source("R/bayesian_path_model_staircase.R")
tar_source("targets_factories.R")

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

  #### Total effect of PIP on post-intervention BMI #### 
  total_treatment_effect_factory(data_clean  = data_clean),

  #### Anthropometric outcomes, joint multivariate model #### 
  anthropometric_outcomes_factory(data_clean = data_clean),

  #### Composite Model =============================================================
  
  #### Stage 1 (anchor): Model A:  Mass_post ~ Treatment + Mass_pre 
  composite_models_stage_1_factory(data_clean = data_clean),

  ####  Stage 2 (a-paths — does treatment move each mediator?)
  composite_models_stage_2_factory(data_clean = data_clean),

  #### Stage 3 (b-paths — does each mediator predict outcome?)
  composite_models_stage_3_factory(data_clean = data_clean),

  ####  Stage 4 (joint mediators) 
  composite_models_stage_4_factory(data_clean = data_clean),

  #### Stage 5 (moderation — if a barrier matters, single-item)
  composite_models_stage_5_factory(data_clean = data_clean),

  #### Stage 6 (full path model — mediation decomposition)
  composite_models_stage_6_factory(data_clean = data_clean)

)

