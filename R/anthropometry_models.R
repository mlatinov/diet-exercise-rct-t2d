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
