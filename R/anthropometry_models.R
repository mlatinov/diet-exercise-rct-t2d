#### Functions to support Anthropometry Bayesian models #####
library(brms)

bmi_total_model <- function(data) {
  
  # Define formula
  formula <- brmsformula(
    bmi_post ~ treatment * bmi_pre, # set the target 
    sigma    ~ treatment            # model sigma as f(treatment)
  )

  # Define priors
  priors <- c(
    set_prior("normal(30, 5)",          class = "Intercept"),
    set_prior("normal(0, 2)",           class = "b", coef = "treatment"),
    set_prior("normal(1, 0.2)",         class = "b", coef = "bmi_pre"),
    set_prior("normal(0, 0.3)",         class = "b", coef = "treatment:bmi_pre"),
    set_prior("normal(log(2.5), 0.5)",  class = "Intercept", dpar = "sigma"),
    set_prior("normal(0, 0.5)",         class = "b", coef = "treatment", dpar = "sigma")
  )

  # Fit the model
  model <- brm(
    formula = formula,
    data = data,
    family = gaussian(),
    prior = priors,
    chains = 4,
    sample_prior = "yes",
    iter = 2000,
    seed = 42
  )

  # Return the fitted model
  return(model)
}

####  Anthropometric outcomes, joint multivariate model ####
anthropometric_model <- function(data) {
  
  # Priors 
  priors <- c(
    # BMI submodel ==============
    set_prior("normal(30, 5)",        class = "Intercept", resp = "bmipost"),
    set_prior("normal(0, 2)",         class = "b", coef = "treatment", resp = "bmipost"),
    set_prior("normal(1, 0.2)",       class = "b", coef = "bmi_pre", resp = "bmipost"),
    set_prior("normal(0, 0.3)",       class = "b", coef = "treatment:bmi_pre", resp = "bmipost"),
    set_prior("normal(0, 3)",         class = "sigma", resp = "bmipost"),
    
    # Waist Submodel ===========
    set_prior("normal(95, 10)",       class = "Intercept", resp = "waistpost"),
    set_prior("normal(0, 5)",         class = "b", coef = "treatment", resp = "waistpost"),
    set_prior("normal(1, 0.2)",       class = "b", coef = "waist_pre", resp = "waistpost"),
    set_prior("normal(0, 0.3)",       class = "b", coef = "treatment:waist_pre", resp = "waistpost"),
    set_prior("normal(0, 5)",         class = "sigma", resp = "waistpost"),
    
    # Hip Submodel =============
    set_prior("normal(100, 10)",      class = "Intercept", resp = "hippost"),
    set_prior("normal(0, 3)",         class = "b", coef = "treatment", resp = "hippost"),
    set_prior("normal(1, 0.2)",       class = "b", coef = "hip_pre", resp = "hippost"),
    set_prior("normal(0, 0.3)",       class = "b", coef = "treatment:hip_pre", resp = "hippost"),
    set_prior("normal(0, 5)",         class = "sigma", resp = "hippost"),
    
    # Residual correlation matrix 
    set_prior("lkj(2)",               class = "rescor")
  )
  # Joint Model 
  joint_model <- brm(
    formula = 
      bf(bmi_post   ~ treatment * bmi_pre) +
      bf(waist_post ~ treatment * waist_pre) +
      bf(hip_post   ~ treatment * hip_pre) +
      set_rescor(TRUE),
    data = data,
    family = gaussian(),
    prior = priors,
    sample_prior = "yes",
    iter = 3000,
    seed = 42
  )
  return(joint_model)
}
