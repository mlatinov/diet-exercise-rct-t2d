#### Data Generative Processes to Support the Validation of all Anthropometry models ==================================
library(mvtnorm)

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
    sd_lambda_SDSCA   = 1,
    sd_lambda_Qpost   = 1,
    sd_lambda_Adhere  = 1,
    sd_lambda_EXR_ADH = 1,
    # Prior SDs for the intercepts nu_k ~ Normal(0, sd)
    sd_nu_SDSCA       = 5,
    sd_nu_Qpost       = 5,
    sd_nu_Adhere      = 1,    
    sd_nu_EXR_ADH     = 1,    
    # Prior SDs for the residuals psi_k ~ HalfNormal(0, sd)
    sd_psi_SDSCA      = 3,
    sd_psi_Qpost      = 3
  ){
  # Draw factor loadings from priors 
  # SDSCA is the marker: 
  lambda_SDSCA   <- abs(rnorm(1, 0, sd_lambda_SDSCA))
  lambda_Qpost   <- rnorm(1, 0, sd_lambda_Qpost)
  lambda_Adhere  <- rnorm(1, 0, sd_lambda_Adhere)
  lambda_EXR_ADH <- rnorm(1, 0, sd_lambda_EXR_ADH)
  
  # Draw indicator intercepts from priors 
  nu_SDSCA   <- rnorm(1, 0, sd_nu_SDSCA)
  nu_Qpost   <- rnorm(1, 0, sd_nu_Qpost)
  nu_Adhere  <- rnorm(1, 0, sd_nu_Adhere)
  nu_EXR_ADH <- rnorm(1, 0, sd_nu_EXR_ADH)
  
  # Draw residual SDs from priors 
  psi_SDSCA <- abs(rnorm(1, 0, sd_psi_SDSCA))
  psi_Qpost <- abs(rnorm(1, 0, sd_psi_Qpost))
  
  # Sample latent variable per observation 
  eta <- rnorm(n, 0, 1)
  
  # Continuous indicators (SDSCA, Qpost) 
  SDSCA <- nu_SDSCA + lambda_SDSCA * eta + rnorm(n, 0, psi_SDSCA)
  Qpost <- nu_Qpost + lambda_Qpost * eta + rnorm(n, 0, psi_Qpost)
  
  # Binary indicators (Adhere, EXR_ADH) via probit 
  p_Adhere  <- pnorm(nu_Adhere  + lambda_Adhere  * eta)
  p_EXR_ADH <- pnorm(nu_EXR_ADH + lambda_EXR_ADH * eta)
  Adhere    <- rbinom(n, 1, p_Adhere)
  EXR_ADH   <- rbinom(n, 1, p_EXR_ADH)
  
  # Return data + true parameters for SBC and recovery checking
  list(
    data = data.frame(
      SDSCA   = SDSCA,
      Qpost   = Qpost,
      Adhere  = Adhere,
      EXR_ADH = EXR_ADH
    ),
    true_params = list(
      lambda_SDSCA   = lambda_SDSCA,
      lambda_Qpost   = lambda_Qpost,
      lambda_Adhere  = lambda_Adhere,
      lambda_EXR_ADH = lambda_EXR_ADH,
      nu_SDSCA       = nu_SDSCA,
      nu_Qpost       = nu_Qpost,
      nu_Adhere      = nu_Adhere,
      nu_EXR_ADH     = nu_EXR_ADH,
      psi_SDSCA      = psi_SDSCA,
      psi_Qpost      = psi_Qpost,
      eta            = eta            
    )
  )
}

