// Include Stan lib Function 
functions {
    #include "../lib/utils.stanfunctions"
    #include "../lib/composites.stanfunctions"
    #include "../lib/diagnostics.stanfunctions"
}

// Data Block 
data{
    int<lower=1> N;
    vector[N] treatment;
    int<lower=0,upper=1> prior_only;

    // Item Block per contruct
    int<lower=1> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[J_mass_post] mass_post_sign;
    int<lower=1> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[J_mass_pre] mass_pre_sign;
}

// Data Transform Block 
transformed data{
    // Build the compostites 
    vector[N] mass_post = composite(mass_post_items, mass_post_sign);
    vector[N] mass_pre  = composite(mass_pre_items,  mass_pre_sign);

    // Z score the composites 
    vector[N] mass_post_stand = zscore(mass_post);
    vector[N] mass_pre_stand  = zscore(mass_pre);
}

// Model Paramters Block 
parameters{
    real alpha;
    real beta_treatment;
    real beta_mass_pre;
    real<lower=0.001> omega;
    real alpha_skew;
}

// Model Block 
model{
    // Priors
    alpha          ~ normal(0, 0.5);
    beta_treatment ~ normal(0, 0.5);
    beta_mass_pre  ~ normal(0.8, 0.3);
    omega          ~ exponential(1);
    alpha_skew     ~ normal(0, 3); 

    // Model Likelihood
    if(prior_only == 0){
        vector[N] xi = 
            alpha 
            + beta_treatment * treatment 
            + beta_mass_pre *  mass_pre_stand;

        mass_post_stand ~ skew_normal(xi, omega, alpha_skew);
    } 
}
generated quantities {

  // LOCATION  =================================================================
  vector[N] xi =
      alpha
      + beta_treatment * treatment
      + beta_mass_pre  * mass_pre_stand;

  // TRUE conditional mean and residual SD  ===================================
  vector[N] mu       = sn_mean(xi, omega, alpha_skew);
  real sd_resid      = sn_sd(omega, alpha_skew);

  // POSTERIOR PREDICTIVE GENERATION ==========================================
  vector[N] mass_post_rep = skew_normal_predictive_rng(xi, omega, alpha_skew);

  // MODEL FIT / INFORMATION CRITERIA =========================================
  vector[N] log_lik = skew_normal_pointwise_loglik(mass_post_stand, xi, omega, alpha_skew);
  real R2           = bayes_R2_skew_normal(mu, omega, alpha_skew);

  // TREATMENT EFFECT ESTIMATION ==============================================
  real treatment_effect     = beta_treatment;
  real treatment_effect_std = standardized_effect(beta_treatment, sd_resid);

  // Adjusted means on the MEAN scale 
  real adjusted_mean_control = sn_mean(
        [alpha + beta_mass_pre * mean(mass_pre_stand)]', omega, alpha_skew)[1];

  real adjusted_mean_treated = sn_mean(
        [alpha + beta_treatment + beta_mass_pre * mean(mass_pre_stand)]', omega, alpha_skew)[1];

  real ATE = adjusted_mean_treated - adjusted_mean_control;   

  // DIRECTIONAL / PROBABILITY STATEMENTS =====================================
  int treatment_reduces_mass    = effect_lt(beta_treatment, 0);

  int treatment_increases_mass  = effect_gt(beta_treatment, 0);
  
  int treatment_in_rope         = in_rope(beta_treatment, -0.10, 0.10);
  
  int treatment_large_reduction = effect_lt(beta_treatment, -0.30);

  // RESIDUAL DIAGNOSTICS ====================================================
  vector[N] raw_resid     = raw_residuals(mass_post_stand, mu);
  vector[N] pearson_resid = pearson_residuals(mass_post_stand, mu, sd_resid);

  // CALIBRATION =============================================================
  vector[N] pit = skew_normal_pit(mass_post_stand, xi, omega, alpha_skew);

  // POSTERIOR PREDICTIVE CHECKS ==============================================
  int p_mean = ppc_indicator_mean(mass_post_stand, mass_post_rep);
  int p_sd   = ppc_indicator_sd(mass_post_stand, mass_post_rep);
  int p_max  = ppc_indicator_max(mass_post_stand, mass_post_rep);

}
