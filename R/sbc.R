
#### Functions to support Simulation Based Calibrations ####
library(SBC)

# SBC for the BMI total effect
sbc_bmi_total <- function(template_data){
  
  # Subset the data to pass it in generator 
  n <- nrow(template_data)
  x_treatment <- as.numeric(template_data$treatment)
  x_bmi_pre   <- template_data$bmi_pre
  
  # Model priors 
  priors <- c(
    brms::set_prior("normal(30, 5)",          class = "Intercept"),
    brms::set_prior("normal(0, 2)",           class = "b", coef = "treatment"),
    brms::set_prior("normal(1, 0.2)",         class = "b", coef = "bmi_pre"),
    brms::set_prior("normal(0, 0.3)",         class = "b", coef = "treatment:bmi_pre"),
    brms::set_prior("normal(log(2.5), 0.5)",  class = "Intercept", dpar = "sigma"),
    brms::set_prior("normal(0, 0.5)",         class = "b", coef = "treatment", dpar = "sigma")
  )
  
  # Generator function 
  generator_fn <- function() {
    # Draw from priors directly (match the model )
    b_Intercept                     <- rnorm(1, 30, 5)
    b_treatment                     <- rnorm(1, 0, 2)
    b_bmi_pre                       <- rnorm(1, 1, 0.2)
    
    # Interaction term
    b_treatment_bmi_pre             <- rnorm(1, 0, 0.3)

    # Draw the Sigma and treatment submodel 
    b_sigma_Intercept               <- rnorm(1, log(2.5), 0.5)
    b_sigma_treatment               <- rnorm(1, 0, 0.5)
    
    # Sub sigma model 
    log_sigma <- b_sigma_Intercept + b_sigma_treatment * x_treatment
    sigma     <- exp(log_sigma)
      
    # Mean model with interaction 
    mu <- b_Intercept + b_treatment * x_treatment + b_bmi_pre * x_bmi_pre +
      b_treatment_bmi_pre * x_treatment * x_bmi_pre
    
    # Simulate response from likelihood
    y  <- rnorm(n, mu, sigma)
    
    # Return list with the requared format variables and generated 
    list(
      variables = list(
        b_Intercept = b_Intercept,
        b_treatment = b_treatment,                             
        b_bmi_pre   = b_bmi_pre,
        "b_treatmentintervention:bmi_pre" = b_treatment_bmi_pre,
        b_sigma_Intercept = b_sigma_Intercept,
        b_sigma_treatment = b_sigma_treatment
      ),
      generated = data.frame(
        bmi_post  = y,
        treatment = x_treatment,
        bmi_pre   = x_bmi_pre
      )
    )
  }
  
  # Generate datasets 
  generator <- SBC_generator_function(generator_fn)
  datasets  <- generate_datasets(generator, n_sims = 100)
  
  # Backend brms SBC fit 
  backend <- SBC_backend_brms(
    bf(
      bmi_post ~ treatment * bmi_pre,
      sigma ~ treatment
    ),
    template_data = template_data,
    prior         = priors,
    family        = gaussian(),
    chains = 2, 
    warmup = 500,
    iter   = 2000,
    init   = 0.1
  )
  
  # Compute SBC 
  results <- compute_SBC(
    datasets, backend,
    keep_fits = FALSE,
    cache_mode     = "results",
    cache_location = "sbc_cache"
  )
  
  # Return the SBC object for ploting  
  return(results)
}
