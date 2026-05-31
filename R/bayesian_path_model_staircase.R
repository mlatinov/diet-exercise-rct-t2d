
#### Functions to run the Stan path models staircase ####

#### Model A:  Mass_post ~ Treatment + Mass_pre ####
model_A <- function(
  data_clean, 
  prior = 0, 
  stan_file_A,
  iter_sampling = 1000,
  chains = 4,
  adapt_delta = 0.95
){
  # Get the model 
  model_a <- cmdstanr::cmdstan_model(stan_file = stan_file_A)
  
  # Fit the model 
  fit_model_a <- model_a$sample(
    data = list(
      N         = nrow(data_clean),
      treatment = data_clean$treatment,

      # Mass post Indice
      J_mass_post     = 2,
      mass_post_items = as.matrix(data_clean[,c("bmi_post","weight_post")]),
      mass_post_sign  = c(1, 1),
      
      # Mass pre Indice
      J_mass_pre     = 2,
      mass_pre_items = as.matrix(data_clean[,c("bmi_pre", "weight_pre")]), 
      mass_pre_sign  = c(1, 1),
      prior_only     = prior
    ),
    iter_sampling = 1000,
    chains        = 4,
    output_dir = "stan_results/",
    seed       = 42,
    adapt_delta = adapt_delta
  )
  return(fit_model_a)
}

#### Model B1: Diet ~ Treatment + Diet_pre ####
model_B1 <- function(
  data_clean,
  prior = 0,
  stan_file,
  iter_sampling = 1000,
  chains = 4
){
  
  # Get the model 
  b1 <- cmdstan_model(stan_file = stan_file)
  
  # Fit the model B1
  fit_model_b1 <- b1$sample(
    data = list(
      N          = nrow(data_clean),
      treatment  = data_clean$treatment,
      diet_pre   = data_clean$diet_score_pre,
      diet_post  = data_clean$diet_score_post,
      prior_only = prior
    ),
    iter_sampling = iter_sampling,
    chains        = chains,
    output_dir    = "stan_results/",
    seed          = 42 
  )
  return(fit_model_b1)
}

#### Model B2: Activity  ~ Treatment + Activity_pre
model_B2 <- function(
  data_clean,
  prior = 0,
  stan_file,
  iter_sampling = 1000,
  chains        = 4
){
  
  # Get the model 
  b2 <- cmdstan_model(stan_file = stan_file)
  
  # Fit the model 
  fit_b2 <- b2$sample(
    data = list(
      N = nrow(data_clean),
      treatment = data_clean$treatment,
      intensity_pre = data_clean$exercise_intensity_pre,
      duration_pre  = data_clean$exercise_duration_pre,
      frequency_pre = data_clean$exercise_type_1_freq_pre,
      intensity_post = data_clean$exercise_intensity_post,
      duration_post  = data_clean$exercise_duration_post,
      frequency_post = data_clean$exercise_type_1_freq_post,
      prior_only     = prior
    ),
    iter_sampling = iter_sampling,
    chains        = chains,
    output_dir = "stan_results/",
    seed = 42
  )
  return(fit_b2)
}

#### Model B3: SelfCare  ~ Treatment + SelfCare_pre ####
model_B3 <- function(
  data_clean,
  prior = 0,
  stan_file,
  iter_sampling = 1000,
  chains        = 4
){
  
  # Get the model 
  b3 <- cmdstan_model(stan_file = stan_file)
  
  # Fit the model 
  fit_b3 <- b3$sample(
    data = list(
      N         = nrow(data_clean),
      treatment = data_clean$treatment,
      self_care_pre    = data_clean$self_care_score_pre,
      J_self_care_post     = 2,
      self_care_post_items = as.matrix(data_clean[,c("self_care_score_post","exercise_adherence")]),
      self_care_post_sing  = c(1, 1),
      prior_only     = prior
    ),
    iter_sampling = iter_sampling,
    chains        = chains,
    output_dir = "stan_results/",
    seed = 42
  )
  return(fit_b3)
}

#### Model C1: Mass_post ~ Treatment + Mass_pre + Diet ####
model_C1 <- function(data_clean, prior = 0){

  # Get the model 
  c1 <- cmdstanr::cmdstan_model(stan_file = "Stan/Model_C1.stan")

  # Fit the model 
  fit_c1 <- c1$sample(
    data = list(
      N         = nrow(data_clean),
      treatment = data_clean$treatment, 
      # Mass post Indice
      J_mass_post     = 2,
      mass_post_items = as.matrix(data_clean[,c("bmi_post","weight_post")]),
      mass_post_sign  = c(1, 1),
      
      # Mass pre Indice
      J_mass_pre     = 2,
      mass_pre_items = as.matrix(data_clean[,c("bmi_pre", "weight_pre")]), 
      mass_pre_sign  = c(1, 1),

      # Diet post Indice
      J_diet_post     = 2,
      diet_post_items = as.matrix(data_clean[,c("diet_score_post","diet_adherence")]),
      diet_post_sign  = c(1, 1),
      prior_only     = prior
    ),
    # Sampler Arguments 
    iter_sampling =  1000,
    chains        = 4,
    output_dir    = "stan_results/",
    seed          = 42 
  )
  return(fit_c1)
}

#### Model C2: Mass_post ~ Treatment + Mass_pre + Activity ####
model_C2 <- function(data_clean, prior = 0){

  # Get the model 
  c2 <- cmdstanr::cmdstan_model(stan_file = "Stan/ Model_C2.stan")

  # Fit the C2 model 
  fit_c2 <- c2$sample(
    data = list(
      N         = nrow(data_clean),
      treatment = data_clean$treatment,

      # Mass post Indice
      J_mass_post     = 2,
      mass_post_items = as.matrix(data_clean[,c("bmi_post","weight_post")]),
      mass_post_sign  = c(1, 1),
      
      # Mass pre Indice
      J_mass_pre     = 2,
      mass_pre_items = as.matrix(data_clean[,c("bmi_pre", "weight_pre")]), 
      mass_pre_sign  = c(1, 1),

      # Activity Post
      intensity_post = data_clean$exercise_intensity_post,
      duration_post  = data_clean$exercise_duration_post,
      frequency_post = data_clean$exercise_type_1_freq_post,
      prior_only     = prior
    ),
    # Sampler Arguments
    iter_sampling = 1000,
    chains        = 4,
    output_dir    = "stan_results/",
    seed          = 42   
  )
  return(fit_c2)

}

#### Model C3: Mass_post ~ Treatment + Mass_pre + SelfCare ####
model_C3 <- function(data_clean, prior = 0){

  # Get the model 
  c3 <- cmdstanr::cmdstan_model(stan_file = "Stan/Model_C3.stan")

  # Fit the C3 model 
  fit_c3 <- c3$sample(
    data = list(
      N         = nrow(data_clean),
      treatment = data_clean$treatment,

      # Mass post Indice
      J_mass_post     = 2,
      mass_post_items = as.matrix(data_clean[,c("bmi_post","weight_post")]),
      mass_post_sign  = c(1, 1),
      
      # Mass pre Indice
      J_mass_pre     = 2,
      mass_pre_items = as.matrix(data_clean[,c("bmi_pre", "weight_pre")]), 
      mass_pre_sign  = c(1, 1),

      # Self Care
      J_self_care_post     = 2,
      self_care_post_items = as.matrix(data_clean[,c("self_care_score_post","exercise_adherence")]),
      self_care_post_sing  = c(1, 1),
      prior_only     = prior
    ),
    # Sampler Arguments
    iter_sampling = 1000,
    chains        = 4,
    output_dir    = "stan_results/",
    seed          = 42   
  )
  return(fit_c3)
}

#### Model D:  Mass_post ~ Treatment + Mass_pre + Diet + Activity + SelfCare ####
model_D <- function(data_clean, prior = 0){

  # Get the model
  d <- cmdstanr::cmdstan_model(stan_file = "Stan/ Model_D.stan")

  # Fit model D 
  fit_d <- d$sample(
    data = list(
      N         = nrow(data_clean),
      treatment = data_clean$treatment,

      # Mass post Indice
      J_mass_post     = 2,
      mass_post_items = as.matrix(data_clean[,c("bmi_post","weight_post")]),
      mass_post_sign  = c(1, 1),
      
      # Mass pre Indice
      J_mass_pre     = 2,
      mass_pre_items = as.matrix(data_clean[,c("bmi_pre", "weight_pre")]), 
      mass_pre_sign  = c(1, 1),

      # Diet post Indice
      J_diet_post     = 2,
      diet_post_items = as.matrix(data_clean[,c("diet_score_post","diet_adherence")]),
      diet_post_sign  = c(1, 1),

      # Self Care
      J_self_care_post     = 2,
      self_care_post_items = as.matrix(data_clean[,c("self_care_score_post","exercise_adherence")]),
      self_care_post_sing  = c(1, 1),

      # Activity Post
      intensity_post = data_clean$exercise_intensity_post,
      duration_post  = data_clean$exercise_duration_post,
      frequency_post = data_clean$exercise_type_1_freq_post,
      prior_only     = prior
    ),
    # Sampler Arguments
    iter_sampling = 1000,
    chains        = 4,
    output_dir    = "stan_results/",
    seed          = 42   
  )
  return(fit_d)
}

#### Model E:  Mass_post ~ Treatment * motivation + Mass_pre + [mediators] ####
model_E1 <- function(data_clean, prior = 0){

  # Get the model 
  e1 <- cmdstanr::cmdstan_model(stan_file = "Stan/Model_E1_motivation.stan")

  # Fit Model E
  fit_e1 <- e1$sample(
      data = list(
      N                  = nrow(data_clean),
      treatment          = data_clean$treatment,
      barrier_motivation = data_clean$lack_motivation, 

      # Mass post Indice
      J_mass_post     = 2,
      mass_post_items = as.matrix(data_clean[,c("bmi_post","weight_post")]),
      mass_post_sign  = c(1, 1),
      
      # Mass pre Indice
      J_mass_pre     = 2,
      mass_pre_items = as.matrix(data_clean[,c("bmi_pre", "weight_pre")]), 
      mass_pre_sign  = c(1, 1),

      # Diet post Indice
      J_diet_post     = 2,
      diet_post_items = as.matrix(data_clean[,c("diet_score_post","diet_adherence")]),
      diet_post_sign  = c(1, 1),

      # Self Care
      J_self_care_post     = 2,
      self_care_post_items = as.matrix(data_clean[,c("self_care_score_post","exercise_adherence")]),
      self_care_post_sing  = c(1, 1),

      # Activity Post
      intensity_post = data_clean$exercise_intensity_post,
      duration_post  = data_clean$exercise_duration_post,
      frequency_post = data_clean$exercise_type_1_freq_post,
      prior_only     = prior
    ),
    # Sampler Arguments
    iter_sampling = 1000,
    chains        = 4,
    output_dir    = "stan_results/",
    seed          = 42   
  )
  return(fit_e1)
}

####  Model F1 (full path model — mediation decomposition) ####
# Diet, Activity, SelfCare ~ Treatment (+ baselines)
 # Mass_post ~ Treatment + Mass_pre + Diet + Activity + SelfCare
model_F1 <- function(data_clean, prior = 0){

  # Get the model 
  f1 <- cmdstanr::cmdstan_model(stan_file = "Stan/Model_F1.stan")

  # Fit the F1 model
  fit_f1 <- f1$sample(
    data = list(
      N                  = nrow(data_clean),
      treatment          = data_clean$treatment,

      # Mass post Indice
      J_mass_post     = 2,
      mass_post_items = as.matrix(data_clean[,c("bmi_post","weight_post")]),
      mass_post_sign  = c(1, 1),
      
      # Mass pre Indice
      J_mass_pre     = 2,
      mass_pre_items = as.matrix(data_clean[,c("bmi_pre", "weight_pre")]), 
      mass_pre_sign  = c(1, 1),

      # Diet Indice
      J_diet_post     = 2,
      diet_post_items = as.matrix(data_clean[,c("diet_score_post","diet_adherence")]),
      diet_post_sign  = c(1, 1),
      diet_pre        = data_clean$diet_score_pre,

      # Self Care
      J_self_care_post     = 2,
      self_care_post_items = as.matrix(data_clean[,c("self_care_score_post","exercise_adherence")]),
      self_care_post_sing  = c(1, 1),
      self_care_pre        = data_clean$self_care_score_pre,

      # Activity 
      intensity_post = data_clean$exercise_intensity_post,
      duration_post  = data_clean$exercise_duration_post,
      frequency_post = data_clean$exercise_type_1_freq_post,
      intensity_pre  = data_clean$exercise_intensity_pre,
      duration_pre   = data_clean$exercise_duration_pre,
      frequency_pre  = data_clean$exercise_type_1_freq_pre,
      prior_only     = prior
    ),
    # Sampler Arguments
    iter_sampling = 1000,
    chains        = 4,
    output_dir    = "stan_results/",
    seed          = 42   
  )
  return(fit_f1)
}




