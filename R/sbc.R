
#### Functions to support Simulation Based Calibrations ####
library(SBC)

# SBC for the BMI total effect
sbc_bmi_total <- function(data){
  
  # Template data
  template_data <- data %>% 
    select(treatment, bmi_post, bmi_pre) %>%
    mutate(treatment = as.numeric(treatment))
  
  # Subset the data to pass it in generator 
  n <- nrow(template_data)
  x_treatment <- as.numeric(template_data$treatment)
  x_bmi_pre   <- template_data$bmi_pre
  
  # Model priors 
  priors <- c(
    brms::prior(normal(30, 5), class = "Intercept"),
    brms::prior(normal(0, 3), class = "b", coef = "treatment"),
    brms::prior(normal(1, 0.5), class = "b", coef = "bmi_pre"),
    brms::prior(normal(0, 2), class = "sigma", lb = 0)
  )
  
  # Generator function 
  generator_fn <- function() {
    # Draw from priors directly (mnatch the model )
    b_Intercept <- rnorm(1, 30, 5)
    b_treatment <- rnorm(1, 0,  3)
    b_bmi_pre   <- rnorm(1, 1,  0.5)
    sigma       <- abs(rnorm(1, 0, 2))            # half-normal(0, 2)
    
    # Simulate response from likelihood
    mu <- b_Intercept + b_treatment * x_treatment + b_bmi_pre * x_bmi_pre
    y  <- rnorm(n, mu, sigma)
    
    # Return list with the requared format variables and generated 
    list(
      variables = list(
        b_Intercept = b_Intercept,
        b_treatment = b_treatment,
        b_bmi_pre   = b_bmi_pre,
        sigma       = sigma
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
    bmi_post ~ treatment + bmi_pre,
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
    cache_mode     = "results",
    cache_location = "sbc_cache"
  )
  
  # Return the SBC object for ploting  
  return(results)
}
