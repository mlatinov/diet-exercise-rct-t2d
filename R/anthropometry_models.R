#### Functions to support Anthropometry Bayesian models #####
library(brms)

bmi_total_model <- function(data) {
  # Convert the factors to numeric
  model_data <- data %>%
    mutate(
      treatment = as.numeric(treatment)
    )

  # Define formula
  formula <- brmsformula(bmi_post ~ treatment + bmi_pre)

  # Define priors
  priors <- c(
    prior(normal(30, 5), class = "Intercept"),
    prior(normal(0, 3), class = "b", coef = "treatment"),
    prior(normal(1, 0.5), class = "b", coef = "bmi_pre"),
    prior(normal(0, 2), class = "sigma", lb = 0)
  )

  # Fit the model
  model <- brm(
    formula = formula,
    data = model_data,
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
