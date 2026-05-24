
#### Functions to run the Stan path models staircase ####

## Model A:  Mass_post ~ Treatment + Mass_pre
model_A <- function(data_clean){
  
  # Get the model 
  model_a <- cmdstanr::cmdstan_model(stan_file = "Stan/ Model_A.stan")

  # Fit the model 
  fit_model_a <- model_a$sample(
    data = list(
      N         = nrow(data_clean),
      treatment = data_clean$treatment,

      # Mass post Indice
      J_mass_post     = 2,
      mass_post_items = as.matrix(data_clean[,c("bmi_post","weight_post")]),
      mass_post_sign  = c(1, 1),
      
      # Mass pre Indice
      J_mass_pre     = 2,
      mass_pre_items = as.matrix(data_clean[,c("bmi_pre", "weight_pre")]), 
      mass_pre_sign  = c(1, 1)
    ),
    iter_sampling = 1000,
    chains        = 4,
    output_dir = "stan_results/",
    seed       = 42 
  )


}