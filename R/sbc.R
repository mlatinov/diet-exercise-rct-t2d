
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

#### SBC on Anthropometric outcomes, joint multivariate model ####
sbc_anthropometric_model <- function(anthropometric_data){
  
  # Subset the data to pass it to the generator 
  n <- nrow(anthropometric_data)
  x_treatment <- anthropometric_data$treatment
  x_bmi_pre   <- anthropometric_data$bmi_pre
  x_waist_pre <- anthropometric_data$waist_pre
  x_hip_pre   <- anthropometric_data$hip_pre

  # Model priors 
  priors <- c(
    # BMI submodel ==============
    brms::set_prior("normal(30, 5)",   class = "Intercept", resp = "bmipost"),
    brms::set_prior("normal(0, 2)",    class = "b", coef = "treatment", resp = "bmipost"),
    brms::set_prior("normal(1, 0.2)",  class = "b", coef = "bmi_pre", resp = "bmipost"),
    brms::set_prior("normal(0, 0.3)",  class = "b", coef = "treatment:bmi_pre", resp = "bmipost"),
    brms::set_prior("normal(0, 3)",    class = "sigma", resp = "bmipost"),
    
    # Waist Submodel ===========
    brms::set_prior("normal(95, 10)",  class = "Intercept", resp = "waistpost"),
    brms::set_prior("normal(0, 5)",    class = "b", coef = "treatment", resp = "waistpost"),
    brms::set_prior("normal(1, 0.2)",  class = "b", coef = "waist_pre", resp = "waistpost"),
    brms::set_prior("normal(0, 0.3)",  class = "b", coef = "treatment:waist_pre", resp = "waistpost"),
    brms::set_prior("normal(0, 5)",    class = "sigma", resp = "waistpost"),
    
    # Hip Submodel =============
    brms::set_prior("normal(100, 10)", class = "Intercept", resp = "hippost"),
    brms::set_prior("normal(0, 3)",    class = "b", coef = "treatment", resp = "hippost"),
    brms::set_prior("normal(1, 0.2)",  class = "b", coef = "hip_pre", resp = "hippost"),
    brms::set_prior("normal(0, 0.3)",  class = "b", coef = "treatment:hip_pre", resp = "hippost"),
    brms::set_prior("normal(0, 5)",    class = "sigma", resp = "hippost"),
    
    # Residual correlation matrix 
    brms::set_prior("lkj(2)",          class = "rescor")
  )
  
  # Generator Function 
  generator_f <- function() {
    # BMI submodel parameters
    bmi_Intercept       <- rnorm(1, 30,  5)
    bmi_treatment       <- rnorm(1, 0,   2)
    b_bmi_pre           <- rnorm(1, 1,   0.2)
    bmi_treatment_pre   <- rnorm(1, 0,   0.3)
    sigma_bmi           <- abs(rnorm(1, 0, 3))
    
    # Waist submodel parameters 
    waist_Intercept     <- rnorm(1, 95,  10)
    waist_treatment     <- rnorm(1, 0,   5)
    b_waist_pre         <- rnorm(1, 1,   0.2)
    waist_treatment_pre <- rnorm(1, 0,   0.3)
    sigma_waist         <- abs(rnorm(1, 0, 5))
    
    # Hip submodel parameters 
    hip_Intercept       <- rnorm(1, 100, 10)
    hip_treatment       <- rnorm(1, 0,   3)
    b_hip_pre           <- rnorm(1, 1,   0.2)
    hip_treatment_pre   <- rnorm(1, 0,   0.3)
    sigma_hip           <- abs(rnorm(1, 0, 5))
    
    # Residual correlation matrix from LKJ 
    R <- clusterGeneration::rcorrmatrix(d = 3, alpha = 2)
    rho_bmi_waist <- R[1, 2]
    rho_bmi_hip   <- R[1, 3]
    rho_waist_hip <- R[2, 3]
    
    # Build Sigma = D R D
    D     <- diag(c(sigma_bmi, sigma_waist, sigma_hip))
    Sigma <- D %*% R %*% D
    
    # Linear predictors 
    mu_bmi <- bmi_Intercept +
      bmi_treatment     * x_treatment +
      b_bmi_pre         * x_bmi_pre +
      bmi_treatment_pre * x_treatment * x_bmi_pre
    
    mu_waist <- waist_Intercept +
      waist_treatment     * x_treatment +
      b_waist_pre         * x_waist_pre +
      waist_treatment_pre * x_treatment * x_waist_pre
    
    mu_hip <- hip_Intercept +
      hip_treatment     * x_treatment +
      b_hip_pre         * x_hip_pre +
      hip_treatment_pre * x_treatment * x_hip_pre
    
    # Multivariate sampling 
    mu_matrix <- cbind(mu_bmi, mu_waist, mu_hip)
    noise     <- mvtnorm::rmvnorm(n, mean = c(0, 0, 0), sigma = Sigma)
    y         <- mu_matrix + noise
    
    # Return 
    list(
      variables = list(
        b_bmipost_Intercept              = bmi_Intercept,
        b_bmipost_treatment              = bmi_treatment,
        b_bmipost_bmi_pre                = b_bmi_pre,
        "b_bmipost_treatment:bmi_pre"    = bmi_treatment_pre,
        sigma_bmipost                    = sigma_bmi,
        
        b_waistpost_Intercept            = waist_Intercept,
        b_waistpost_treatment            = waist_treatment,
        b_waistpost_waist_pre            = b_waist_pre,
        "b_waistpost_treatment:waist_pre" = waist_treatment_pre,
        sigma_waistpost                  = sigma_waist,
        
        b_hippost_Intercept              = hip_Intercept,
        b_hippost_treatment              = hip_treatment,
        b_hippost_hip_pre                = b_hip_pre,
        "b_hippost_treatment:hip_pre"    = hip_treatment_pre,
        sigma_hippost                    = sigma_hip,
        
        "rescor__bmipost__waistpost"  = rho_bmi_waist,
        "rescor__bmipost__hippost"    = rho_bmi_hip,
        "rescor__waistpost__hippost"  = rho_waist_hip
      ),
      generated = data.frame(
        bmi_post   = y[, 1],
        waist_post = y[, 2],
        hip_post   = y[, 3],
        treatment  = x_treatment,
        bmi_pre    = x_bmi_pre,
        waist_pre  = x_waist_pre,
        hip_pre    = x_hip_pre
      )
    )
  }
  # Generate datasets 
  generator <- SBC_generator_function(generator_f)
  datasets  <- generate_datasets(generator, n_sims = 100)
  
  # Backend brms SBC fit 
  backend <- SBC_backend_brms(
    formula = 
      bf(bmi_post   ~ treatment * bmi_pre) +
      bf(waist_post ~ treatment * waist_pre) +
      bf(hip_post   ~ treatment * hip_pre) +
      set_rescor(TRUE),
    template_data = anthropometric_data,
    prior         = priors,
    family        = gaussian(),
    chains = 2, 
    warmup = 500,
    iter   = 3000,
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
