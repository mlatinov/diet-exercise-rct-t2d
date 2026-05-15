#### Data Generative Processes to Support the Validation of all Anthropometry models ==================================


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
