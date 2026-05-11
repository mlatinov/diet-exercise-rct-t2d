#### Data Generative Processes to Support the Validation of all Anthropometry models ==================================


## Total effect of PIP on post-intervention BMI
dgp_total_bmi <- function(
  n,
  treatment_p = 0.5,
  bmi_pre_mu = 30, bmi_pre_sd = 5, bmi_post_sd = 2,
  beta_treatment = -2, beta_bmi_pre = 0.95, beta_0 = 0
) {
  # Simulate Treatment
  treatment <- rbinom(n = n, size = 1, prob = treatment_p)

  # BMI Pre
  bmi_pre <- rnorm(n = n, mean = bmi_pre_mu, sd = bmi_pre_sd)

  # BMI Post
  bmi_post_mu <- beta_0 + beta_treatment * treatment + beta_bmi_pre * bmi_pre

  # Sample
  bmi_post <- rnorm(n = n, mean = bmi_post_mu, sd = bmi_post_sd)

  # Combine and return a dataframe
  data <- data.frame(
    treatment = treatment,
    bmi_pre = bmi_pre,
    bmi_post = bmi_post
  )
  return(data)
}
