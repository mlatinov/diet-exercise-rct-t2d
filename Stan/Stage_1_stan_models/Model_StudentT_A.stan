// Include Stan lib Function 
functions {
    #include "../lib/utils.stanfunctions"
    #include "../lib/composites.stanfunctions"
    #include "../lib/diagnostics.stanfunctions"
}

// Data Input Block 
data{
    int<lower=1> N;
    vector[N] treatment;
    int<lower=0,upper=1> prior_only;

    // Item Block per contruct
    int<lower=1> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[J_mass_post] mass_post_sign;
    int<lower=1> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[J_mass_pre] mass_pre_sign;
}

// Data Transformation block 
transformed data{
    // Build the Composites
    vector[N] mass_post = composite(mass_post_items, mass_post_sign); 
    vector[N] mass_pre  = composite(mass_pre_items, mass_pre_sign);

    // Z transformation 
    vector[N] mass_post_stand = zscore(mass_post);
    vector[N] mass_pre_stand  = zscore(mass_pre);
}

// Model Parameters
parameters{
    real alpha;
    real beta_treatment;
    real beta_mass_pre;
    real<lower=0.001> sigma;
    real<lower=1> nu;
}

// Model Block 
model{
    // Priors
    alpha          ~ normal(0, 0.5);
    beta_treatment ~ normal(0, 0.5);
    beta_mass_pre  ~ normal(0.8, 0.3);
    sigma          ~ exponential(1);
    nu             ~ gamma(2, 0.1);

    // Model Likelihood 
    if(prior_only == 0){
        mass_post_stand ~ student_t(nu,
            alpha
            + beta_treatment * treatment
            + beta_mass_pre  * mass_pre_stand
            ,sigma
        );
    }
}
// Additional Calculations
generated quantities {

  // LINEAR PREDICTOR ==========================================================
  vector[N] mu =
      alpha
      + beta_treatment * treatment
      + beta_mass_pre  * mass_pre_stand;

  // TRUE residual ============================================================
  real sd_resid = student_t_sd(sigma, nu);

  // POSTERIOR PREDICTIVE GENERATION ==========================================
  vector[N] mass_post_rep = student_t_predictive_rng(mu, sigma, nu);

  // MODEL FIT / INFORMATION CRITERIA =========================================
  vector[N] log_lik = student_t_pointwise_loglik(mass_post_stand, mu, sigma, nu);
  real R2           = bayes_R2_student_t(mu, sigma, nu);

  // TREATMENT EFFECT ESTIMATION ==============================================
  real treatment_effect     = beta_treatment;
  real treatment_effect_std = standardized_effect(beta_treatment, sd_resid);

  real adjusted_mean_control = alpha + beta_mass_pre * mean(mass_pre_stand);
  real adjusted_mean_treated = alpha + beta_treatment + beta_mass_pre * mean(mass_pre_stand);
  
  real ATE = adjusted_mean_treated - adjusted_mean_control; 

  // DIRECTIONAL / PROBABILITY STATEMENTS =====================================
  int treatment_reduces_mass    = effect_lt(beta_treatment, 0);
  int treatment_increases_mass  = effect_gt(beta_treatment, 0);
  int treatment_in_rope         = in_rope(beta_treatment, -0.10, 0.10);
  int treatment_large_reduction = effect_lt(beta_treatment, -0.30);

  // RESIDUAL DIAGNOSTICS (use heavy-tail-corrected SD) =======================
  vector[N] raw_resid     = raw_residuals(mass_post_stand, mu);
  vector[N] pearson_resid = pearson_residuals(mass_post_stand, mu, sd_resid);

  // CALIBRATION (Student-t PIT) ==============================================
  vector[N] pit = student_t_pit(mass_post_stand, mu, sigma, nu);

  // POSTERIOR PREDICTIVE CHECKS ==============================================
  int p_mean = ppc_indicator_mean(mass_post_stand, mass_post_rep);
  int p_sd   = ppc_indicator_sd(mass_post_stand, mass_post_rep);
  int p_max  = ppc_indicator_max(mass_post_stand, mass_post_rep);

}
