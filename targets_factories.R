
#### Total effect of PIP on post-intervention BMI factory ####
total_treatment_effect_factory <- function(data_clean){
  list(
    ## Prepare the data 
    tar_target(
      name = bmi_total_effect_data,
      command = prepare_data_bmi_total(data_clean)
    ),
    # Generative Model for Param Recover
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
    )
  )
}

#### Anthropometric outcomes, joint multivariate model factory ####
anthropometric_outcomes_factory <- function(data_clean){
  list(
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
    )
  )
}

#### Composite models Stage 1 factory #####
composite_models_stage_1_factory <- function(data_clean){
  list(
    ## Normal Likelihood
    tar_target(
      name = stage_1_model_A_normal_pp,
      command = model_A(data_clean, prior = 1, stan_file_A = "Stan/Stage_1_stan_models/ Model_A.stan")
    ),
    tar_target(
      name = stage_1_model_A_normal,
      command = model_A(data_clean, prior = 0, stan_file_A = "Stan/Stage_1_stan_models/ Model_A.stan")
    ),
    ## Skew Normal Likelihood
    tar_target(
      name = stage_1_model_A_skewnormal_pp,
      command = model_A(data_clean, prior = 1, stan_file_A = "Stan/Stage_1_stan_models/Model_Skew_A.stan")
    ),
    tar_target(
      name = stage_1_model_A_skewnormal,
      command = model_A(data_clean, prior = 0, stan_file_A = "Stan/Stage_1_stan_models/Model_Skew_A.stan")
    ),
    ## Student T Likelihood
    tar_target(
      name = stage_1_model_A_tnormal_pp,
      command = model_A(data_clean, prior = 1, stan_file_A = "Stan/Stage_1_stan_models/Model_StudentT_A.stan")
    ),
    tar_target(
      name = stage_1_model_A_tnormal,
      command = model_A(data_clean, prior = 0, stan_file_A = "Stan/Stage_1_stan_models/Model_StudentT_A.stan")
    )
  )
}

#### Composite models Stage 2 factory #####
composite_models_stage_2_factory <- function(data_clean){
  list(
    # Model B1: Diet ~ Treatment + Diet_pre
    tar_target(
      name = stage_2_model_B1,
      command = model_B1(data_clean)
    ),
    tar_target(
      name = stage_2_model_B1_prior,
      command = model_B1(data_clean, prior = 1)
    ),
    # Model B2: Activity  ~ Treatment + Activity_pre
    tar_target(
      name = stage_2_model_B2,
      command = model_B2(data_clean)
    ),
    tar_target(
      name = stage_2_model_B2_prior,
      command = model_B2(data_clean, prior = 1)
    ),
    # Model B3: SelfCare  ~ Treatment + SelfCare_pre
    tar_target(
      name = stage_2_model_B3,
      command = model_B3(data_clean)
    ),
    tar_target(
      name = stage_2_model_B3_prior,
      command = model_B3(data_clean, prior = 1)
    )
  )
}

#### Composite models Stage 3 factory #####
composite_models_stage_3_factory <- function(data_clean){
  list(
    # Model C1: Mass_post ~ Treatment + Mass_pre + Diet
    tar_target(
      name = stage_3_model_C1,
      command = model_C1(data_clean)
    ),
    tar_target(
      name = stage_3_model_C1_prior,
      command = model_C1(data_clean, prior = 1)
    ),
    # Model C2: Mass_post ~ Treatment + Mass_pre + Activity
    tar_target(
      name = stage_3_model_C2,
      command = model_C2(data_clean) 
    ),
    tar_target(
      name = stage_3_model_C2_prior,
      command = model_C2(data_clean, prior = 1) 
    ),
    # Model C3: Mass_post ~ Treatment + Mass_pre + SelfCare
    tar_target(
      name = stage_3_model_C3,
      command = model_C3(data_clean)
    ),
    tar_target(
      name = stage_3_model_C3_prior,
      command = model_C3(data_clean, prior = 1)
    )
  )
}

#### Composite models Stage 4 factory #####
composite_models_stage_4_factory <- function(data_clean){
  list(
    # Model D:  Mass_post ~ Treatment + Mass_pre + Diet + Activity + SelfCare
    tar_target(
      name = stage_4_model_D,
      command = model_D(data_clean)
    ),
    tar_target(
      name = stage_4_model_D_prior,
      command = model_D(data_clean, prior = 1)
    )
  )
}

#### Composite models Stage 5 factory #####
composite_models_stage_5_factory <- function(data_clean){
  list(
    # Model E:  Mass_post ~ Treatment * [chosen_barrier] + Mass_pre + [mediators]
    tar_target(
      name = stage_5_model_E1_motivation,
      command = model_E1(data_clean) 
    ),
    tar_target(
      name = stage_5_model_E1_motivation_prior,
      command = model_E1(data_clean, prior = 1) 
    )
  )
}

#### Composite models Stage 6 factory #####
composite_models_stage_6_factory <- function(data_clean){
  list(
    # Diet, Activity, SelfCare ~ Treatment (+ baselines)
    # # Mass_post ~ Treatment + Mass_pre + Diet + Activity + SelfCare
    tar_target(
      name = stage_6_model_F1,
      command = model_F1(data_clean)
    ),
    tar_target(
      name = stage_6_model_F1_prior,
      command = model_F1(data_clean, prior = 1)
    )
  )
} 
