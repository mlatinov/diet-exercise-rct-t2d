// Include Stan lib Function 
functions {
    #include "lib/utils.stanfunctions"
    #include "lib/composites.stanfunctions"
    #include "lib/diagnostics.stanfunctions"
}

// Input Data Block 
data{
    int<lower=1> N;
    vector[N] treatment;
    int<lower=0,upper=1> prior_only;

    // Composite Building Blocks =========================================================================

    // Mass Blocks 
    int<lower=0> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[J_mass_pre] mass_pre_sign;
    int<lower=0> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[J_mass_post] mass_post_sign;

    // Diet Block
    int<lower=0> J_diet_post; matrix[N, J_diet_post] diet_post_items; vector[J_diet_post] diet_post_sign;
}
// Data Transformation 
transformed data {
   // Build the Composites 
   vector[N] mass_pre  = composite(mass_pre_items, mass_pre_sign);
   vector[N] mass_post = composite(mass_post_items, mass_post_sign);
   vector[N] diet_post = composite(diet_post_items, diet_post_sign);

   // Standartize all Composits
   vector[N] mass_pre_stand   = zscore(mass_pre);
   vector[N] mass_post_stand  = zscore(mass_post);
   vector[N] diet_post_stand  = zscore(diet_post);

    // For later back transformation 
   real mass_post_composite_sd = sd(mass_post);
}
// Model Paramters 
parameters{
    real alpha;
    real beta_treatment;
    real beta_mass_pre;
    real beta_diet_post;
    real<lower=0.001> sigma;
}
// Model Block
model{
    // Priors
    alpha ~ normal(0, 0.5);
    beta_treatment ~ normal(0, 0.5);
    beta_mass_pre  ~ normal(0.8, 0.3);
    beta_diet_post ~ normal(0, 0.5);
    sigma ~ exponential(1);

    // Model Likelihood 
    if(prior_only == 0){
        mass_post_stand ~ normal(
            alpha 
            + beta_treatment * treatment 
            + beta_mass_pre  * mass_pre_stand 
            + beta_diet_post * diet_post_stand
            ,sigma
        );
    }
} 
// Aditional Calculations 
generated quantities {

    // EXPECTED VALUES / LINEAR PREDICTOR ====================================
    vector[N] mu =
        alpha
        + beta_treatment * treatment
        + beta_mass_pre  * mass_pre_stand
        + beta_diet_post * diet_post_stand;

    // POSTERIOR PREDICTIVE GENERATION =========================================
    vector[N] mass_post_rep = normal_predictive_rng(mu, sigma);

    // MODEL FIT / INFORMATION CRITERIA ========================================
    vector[N] log_lik = normal_pointwise_loglik(mass_post_stand, mu, sigma);
    real R2           = bayes_R2_gaussian(mu, sigma);

    // EFFECT ESTIMATION ========================================================

    // Direct treatment effect after conditioning on diet 
    real direct_treatment_effect = beta_treatment;
    real direct_treatment_effect_std = standardized_effect(beta_treatment, sigma);

    // Diet effect on metabolic burden 
    real diet_effect = beta_diet_post;
    real diet_effect_std = standardized_effect(beta_diet_post, sigma);

    // ADJUSTED EXPECTED OUTCOMES =================================================

    // Expected metabolic burden at average baseline mass + average diet
    real adjusted_mean_control = 
        alpha
        + beta_mass_pre * mean(mass_pre_stand)
        + beta_diet_post * mean(diet_post_stand);

    real adjusted_mean_treated =
        alpha
        + beta_treatment
        + beta_mass_pre * mean(mass_pre_stand)
        + beta_diet_post * mean(diet_post_stand);

    // Residual treatment effect after adjusting for diet
    real residual_ATE = adjusted_mean_treated- adjusted_mean_control;

    // DIET EFFECT INTERPRETATION ================================================

    // Better diet reduces metabolic burden
    int diet_reduces_mass = effect_lt(beta_diet_post, 0);

    // Better diet worsens metabolic burden
    int diet_increases_mass = effect_gt(beta_diet_post, 0);

    // Diet effect practically negligible
    int diet_effect_in_rope = in_rope(beta_diet_post, -0.10, 0.10);

    // Strong clinically meaningful diet effect
    int diet_large_effect = effect_lt(beta_diet_post, -0.30);

    // DIRECT TREATMENT EFFECT INTERPRETATION =======================================

    // Remaining treatment effect after conditioning on diet
    int treatment_reduces_mass = effect_lt(beta_treatment, 0);
    int treatment_effect_in_rope = in_rope(beta_treatment, -0.10, 0.10);

    // RESIDUAL DIAGNOSTICS =========================================================
    vector[N] raw_resid     = raw_residuals(mass_post_stand, mu);
    vector[N] pearson_resid = pearson_residuals(mass_post_stand, mu, sigma);

    // CALIBRATION DIAGNOSTICS ======================================================
    vector[N] pit = normal_pit(mass_post_stand, mu, sigma);

    // POSTERIOR PREDICTIVE CHECKS ==================================================
    int p_mean = ppc_indicator_mean(mass_post_stand, mass_post_rep);
    int p_sd   = ppc_indicator_sd(mass_post_stand, mass_post_rep);
    int p_max  = ppc_indicator_max(mass_post_stand, mass_post_rep);

}
