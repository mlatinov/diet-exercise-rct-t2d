#### Data Generative Processes to Support the Validation of all Anthropometry models ==================================
library(mvtnorm)
source("R/clean_data_raw.R")

## Total effect of PIP on post-intervention BMI
dgp_total_bmi <- function(
    n,
    treatment_p          = 0.5,
    bmi_pre_mu           = 30,
    bmi_pre_sd           = 5,
    beta_0               = 0,
    beta_treatment       = -2,
    beta_bmi_pre         = 0.95,
    treatment_bmi_pre    = 0.3,
    sigma_intercept      = log(3),   
    beta_sigma_treatment = 0         
) {
  # Simulate Treatment ,BMI Pre, 
  treatment <- rbinom(n = n, size = 1, prob = treatment_p)
  bmi_pre <- rnorm(n = n, mean = bmi_pre_mu, sd = bmi_pre_sd)
  
  # Center the Bmi pre
  bmi_pre_c <- bmi_pre - mean(bmi_pre)
  
  # Mean Submodel 
  bmi_post_mu <- beta_0 + beta_treatment * treatment + beta_bmi_pre * bmi_pre_c +
    treatment_bmi_pre * treatment * bmi_pre_c
  
  # Sigma Submodel 
  log_sigma <- sigma_intercept + beta_sigma_treatment * treatment
  sigma     <- exp(log_sigma)
  
  # Sample
  bmi_post <- rnorm(n = n, mean = bmi_post_mu, sd = sigma)

  # Combine and return a dataframe
  data <- data.frame(
    treatment = treatment,
    bmi_pre = bmi_pre_c,
    bmi_post = bmi_post
  )
  return(data)
}

#### Direct Effect of PIP on BMI post intervention ####
dgp_direct_bmi <- function(
    n,
    treatment_p          = 0.5,
    bmi_pre_mu           = 30,
    bmi_pre_sd           = 5,
    height_mu            = 165,
    height_sd            = 15,
    weight_post_mu       = 80,
    weight_post_sd       = 15,
    beta_0               = 0,
    beta_treatment       = -2,
    beta_bmi_pre         = 0.95,
    beta_height          = -0.4,
    beta_weight_post     = 0.3,
    treatment_bmi_pre    = 0.3,
    sigma_intercept      = log(3),   
    beta_sigma_treatment = 0    
  ){
  # Adjustment set contains BMIPre, Height, WtPost 
  
  # Simulate Treatment ,BMI Pre, Height , Wt Post 
  treatment <- rbinom(n = n, size = 1, prob = treatment_p)
  bmi_pre <- rnorm(n = n, mean = bmi_pre_mu, sd = bmi_pre_sd)
  height  <- rnorm(n = n, mean = height_mu, sd = height_sd )
  weight_post <- rnorm(n = n, mean = weight_post_mu, sd = weight_post_sd)
  
  # Center all the numerical predictors 
  bmi_pre_c <- bmi_pre - mean(bmi_pre)
  height_c  <- height  - mean(height)
  weight_c  <- weight_post  - mean(weight_post)
  
  # Mean Submodel 
  bmi_post_mu <- beta_0 
    + beta_treatment * treatment 
    + beta_height  * height_c
    + beta_weight_post * weight_c
    + beta_bmi_pre * bmi_pre_c 
    + treatment_bmi_pre * treatment * bmi_pre_c
  
  # Sigma Submodel 
  log_sigma <- sigma_intercept + beta_sigma_treatment * treatment
  sigma     <- exp(log_sigma)
  
  # Sample
  bmi_post <- rnorm(n = n, mean = bmi_post_mu, sd = sigma)
  
  # Combine and return a dataframe
  data <- data.frame(
    treatment = treatment,
    bmi_pre = bmi_pre_c,
    height  = height_c,
    weight_post = weight_c,
    bmi_post = bmi_post
  )
  return(data)
}

#### Function to Simulate the Anthropometric outcomes, joint multivariate model
dgp_anthropometric <- function(
    n = 200,
    treatment_p = 0.5,
    
    # Predictor distributions
    bmi_pre_mu   = 30, bmi_pre_sd   = 5,
    waist_pre_mu = 95, waist_pre_sd = 10,
    hip_pre_mu   = 100, hip_pre_sd  = 10,
    
    # BMI submodel
    bmi_beta_0           = 0,
    bmi_beta_treatment   = -2,
    bmi_beta_pre         = 0.95,
    bmi_treatment_pre    = 0.3,
    
    # Waist submodel
    waist_beta_0         = 0,
    waist_beta_treatment = -3,
    waist_beta_pre       = 0.95,
    waist_treatment_pre  = 0.2,
    
    # Hip submodel
    hip_beta_0           = 0,
    hip_beta_treatment   = -1.5,
    hip_beta_pre         = 0.95,
    hip_treatment_pre    = 0.1,
    
    # Residual SDs per outcome
    sigma_bmi   = 2.5,
    sigma_waist = 4,
    sigma_hip   = 4,
    
    # Residual correlations Rho between outcomes
    rho_bmi_waist = 0.6,
    rho_bmi_hip   = 0.5,
    rho_waist_hip = 0.7
  ){
  
  # Simulate Treatment ,BMI Pre, Weist Pre, Hip Pre
  treatment <- rbinom(n = n, size = 1, prob = treatment_p)
  bmi_pre   <- rnorm(n = n, mean = bmi_pre_mu, sd = bmi_pre_sd)
  waist_pre <- rnorm(n = n, mean = waist_pre_mu, sd = waist_pre_sd)
  hip_pre   <- rnorm(n = n, mean = hip_pre_mu,   sd = hip_pre_sd)
  
  # Center all the numerical predictors 
  bmi_pre_c <- bmi_pre - mean(bmi_pre)
  waist_pre_c  <- waist_pre  - mean(waist_pre)
  hip_pre_c  <- hip_pre  - mean(hip_pre)
  
  # BMI Post Submodel 
  bmi_post_mu <- (
    bmi_beta_0 
    + bmi_beta_treatment * treatment 
    + bmi_beta_pre * bmi_pre_c
  )
  
  # Waist Post Submodel 
  waist_post_mu <- (
    waist_beta_0 
    + waist_beta_treatment * treatment 
    + waist_beta_pre * waist_pre_c
  ) 
  
  # Hip Post Submodel 
  hip_post_mu <-(
    hip_beta_0 
    + hip_beta_treatment * treatment 
    + hip_beta_pre * hip_pre_c
  )
  
  # Build mu matrix 
  mu_matrix <- cbind(
    bmi_post_mu,
    waist_post_mu,
    hip_post_mu
  )
  
  # Build Sum Matrix 
  R <- matrix(c(
    1,             rho_bmi_waist, rho_bmi_hip,
    rho_bmi_waist, 1,             rho_waist_hip,
    rho_bmi_hip,   rho_waist_hip, 1
  ), nrow = 3, byrow = TRUE)
  
  # Diagonalize it 
  D <- diag(c(sigma_bmi, sigma_waist, sigma_hip))
  Sigma <- D %*% R %*% D
  
  # Draw n Obervations 
  noise <- rmvnorm(n, mean = c(0, 0, 0), sigma = Sigma)
  y     <- mu_matrix + noise
  
  # Return a Data frame 
  data <-   data.frame(
    treatment  = treatment,
    bmi_pre    = bmi_pre_c,
    waist_pre  = waist_pre_c,
    hip_pre    = hip_pre_c,
    bmi_post   = y[, 1],
    waist_post = y[, 2],
    hip_post   = y[, 3]
  )
  return(data)
}

####  Confirmatory factor analysis on the behaviour indicators ####
dgp_cfa_behaviour <- function(
    n = 200,
    # Prior SDs for the loadings lambda_k ~ Normal(0, sd)
    sd_lambda_self_care_score_post = 1,
    sd_lambda_q_score_post         = 1,
    sd_lambda_diet_adherence     = 1,
    sd_lambda_exercise_adherence = 1,
    # Prior SDs for the intercepts nu_k ~ Normal(0, sd)
    sd_nu_self_care_score_post  = 5,
    sd_nu_q_score_post    = 5,
    sd_nu_diet_adherence = 1,    
    sd_nu_exercise_adherence = 1,    
    # Prior SDs for the residuals psi_k ~ HalfNormal(0, sd)
    sd_psi_self_care_score_post  = 3,
    sd_psi_q_score_post          = 3
  ){
  # Draw factor loadings from priors 
  # SDSCA is the marker: 
  lambda_self_care_score_post <- abs(rnorm(1, 0, sd_lambda_self_care_score_post))
  lambda_q_score_post    <- rnorm(1, 0, sd_lambda_q_score_post)
  lambda_diet_adherence  <- rnorm(1, 0, sd_lambda_diet_adherence)
  lambda_exercise_adherence <- rnorm(1, 0, sd_lambda_exercise_adherence)
  
  # Draw indicator intercepts from priors 
  nu_self_care_score_post   <- rnorm(1, 0, sd_nu_self_care_score_post)
  nu_q_score_post    <- rnorm(1, 0, sd_nu_q_score_post)
  nu_diet_adherence  <- rnorm(1, 0, sd_nu_diet_adherence)
  nu_exercise_adherence <- rnorm(1, 0, sd_nu_exercise_adherence)
  
  # Draw residual SDs from priors 
  psi_self_care_score_post <- abs(rnorm(1, 0, sd_psi_self_care_score_post))
  psi_q_score_post         <- abs(rnorm(1, 0, sd_psi_q_score_post))
  
  # Sample latent variable per observation 
  eta <- rnorm(n, 0, 1)
  
  # Continuous indicators (SDSCA, Qpost) 
  self_care_score_post <- ( 
    nu_self_care_score_post 
    + lambda_self_care_score_post 
    * eta 
    + rnorm(n, 0, psi_self_care_score_post)
  )
  q_score_post <- ( 
    nu_q_score_post 
    + lambda_q_score_post  
    * eta 
    + rnorm(n, 0, psi_q_score_post)
  )
  # Binary indicators (Adhere, EXR_ADH) via probit 
  p_diet_adherence     <- pnorm(nu_diet_adherence  + lambda_diet_adherence  * eta)
  p_exercise_adherence <- pnorm(nu_exercise_adherence + lambda_exercise_adherence * eta)
  diet_adherence     <- rbinom(n, 1, p_diet_adherence)
  exercise_adherence <- rbinom(n, 1, p_exercise_adherence)
  
  # Prepare the data to be used with Stan 
  raw_dgp_data <- data_frame(
    self_care_score_post  = self_care_score_post,
    q_score_post   = q_score_post,
    diet_adherence = diet_adherence,
    exercise_adherence = exercise_adherence
  )
  stan_data <- prepare_cfa_behaviour_data(raw_dgp_data)
  
  # Return Stan data + true parameters for SBC and recovery checking
  list(
    data = stan_data,
    true_params = list(
      # Lambdas
      lambda_self_care_score_post = lambda_self_care_score_post,
      lambda_q_score_post         = lambda_q_score_post ,
      lambda_diet_adherence     = lambda_diet_adherence,
      lambda_exercise_adherence = lambda_exercise_adherence,
      # Nu
      nu_self_care_score_post = nu_self_care_score_post,
      nu_q_score_post         = nu_q_score_post ,
      nu_diet_adherence      = nu_diet_adherence,
      nu_exercise_adherence  = nu_exercise_adherence,
      # PSI
      psi_self_care_score_post = psi_self_care_score_post,
      psi_q_score_post         = psi_q_score_post,
      # Eta 
      eta  = eta            
    )
  )
}
####  CFA on the behaviour indicators — 3 continuous + 2 binary ####
dgp_cfa_behaviour2 <- function(
    n = 200,
    # Prior SDs for the loadings  lambda_k ~ Normal(0, sd)
    sd_lambda_self_care_score_post = 1,
    sd_lambda_diet_score_post      = 1,
    sd_lambda_exercise_post_total  = 1,
    sd_lambda_diet_adherence       = 1,
    sd_lambda_exercise_adherence   = 1,
    # Prior SDs for the intercepts  nu_k ~ Normal(0, sd)
    sd_nu_self_care_score_post     = 0.5,
    sd_nu_diet_score_post          = 0.5,
    sd_nu_exercise_post_total      = 0.5,
    sd_nu_diet_adherence           = 0.5,
    sd_nu_exercise_adherence       = 0.5,
    # Residual SDs for continuous indicators  psi_k ~ LogNormal(0, sd)
    meanlog_psi = 0,
    sdlog_psi   = 0.5
){
  
  # Draw factor loadings from priors 
  lambda_self_care_score_post <- rnorm(1, 0, sd_lambda_self_care_score_post)
  lambda_diet_score_post      <- rnorm(1, 0, sd_lambda_diet_score_post)
  lambda_exercise_post_total  <- rnorm(1, 0, sd_lambda_exercise_post_total)
  lambda_diet_adherence       <- rnorm(1, 0, sd_lambda_diet_adherence)
  lambda_exercise_adherence   <- rnorm(1, 0, sd_lambda_exercise_adherence)
  
  # Draw indicator intercepts from priors 
  nu_self_care_score_post <- rnorm(1, 0, sd_nu_self_care_score_post)
  nu_diet_score_post      <- rnorm(1, 0, sd_nu_diet_score_post)
  nu_exercise_post_total  <- rnorm(1, 0, sd_nu_exercise_post_total)
  nu_diet_adherence       <- rnorm(1, 0, sd_nu_diet_adherence)
  nu_exercise_adherence   <- rnorm(1, 0, sd_nu_exercise_adherence)
  
  # Draw residual SDs for continuous indicators 
  psi_self_care_score_post <- rlnorm(1, meanlog_psi, sdlog_psi)
  psi_diet_score_post      <- rlnorm(1, meanlog_psi, sdlog_psi)
  psi_exercise_post_total  <- rlnorm(1, meanlog_psi, sdlog_psi)
  
  # Sample latent variable per observation (unit variance) 
  eta <- rnorm(n, 0, 1)
  
  # Continuous indicators 
  self_care_score_post <- nu_self_care_score_post +
    lambda_self_care_score_post * eta +
    rnorm(n, 0, psi_self_care_score_post)
  
  diet_score_post <- nu_diet_score_post +
    lambda_diet_score_post * eta +
    rnorm(n, 0, psi_diet_score_post)
  
  exercise_post_total <- nu_exercise_post_total +
    lambda_exercise_post_total * eta +
    rnorm(n, 0, psi_exercise_post_total)
  
  # Binary indicators via probit
  p_diet_adherence     <- pnorm(nu_diet_adherence     + lambda_diet_adherence     * eta)
  p_exercise_adherence <- pnorm(nu_exercise_adherence + lambda_exercise_adherence * eta)
  diet_adherence       <- rbinom(n, 1, p_diet_adherence)
  exercise_adherence   <- rbinom(n, 1, p_exercise_adherence)
  
  # Assemble raw data 
  raw_dgp_data <- data.frame(
    self_care_score_post = self_care_score_post,
    diet_score_post      = diet_score_post,
    exercise_post_total  = exercise_post_total,
    diet_adherence       = diet_adherence,
    exercise_adherence   = exercise_adherence
  )
  
  stan_data <- prepare_cfa_behaviour_data(raw_dgp_data)
  
  # Return Stan data + true parameters
  list(
    data = stan_data,
    true_params = list(
      # Loadings 
      lambda_self_care_score_post = lambda_self_care_score_post,
      lambda_diet_score_post      = lambda_diet_score_post,
      lambda_exercise_post_total  = lambda_exercise_post_total,
      lambda_diet_adherence       = lambda_diet_adherence,
      lambda_exercise_adherence   = lambda_exercise_adherence,
      # Intercepts
      nu_self_care_score_post = nu_self_care_score_post,
      nu_diet_score_post      = nu_diet_score_post,
      nu_exercise_post_total  = nu_exercise_post_total,
      nu_diet_adherence       = nu_diet_adherence,
      nu_exercise_adherence   = nu_exercise_adherence,
      # Residual SDs 
      psi_self_care_score_post = psi_self_care_score_post,
      psi_diet_score_post      = psi_diet_score_post,
      psi_exercise_post_total  = psi_exercise_post_total,
      # Latent values
      eta = eta
    )
  )
}

