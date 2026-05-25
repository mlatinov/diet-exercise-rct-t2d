// Include Stan lib Function 
functions {
    #include "lib/utils.stanfunctions"
    #include "lib/composites.stanfunctions"
    #include "lib/diagnostics.stanfunctions"
}
// Data input Block 
data{
    int<lower=1> N;
    vector[N] treatment;

    // Composite Building block 
    int<lower=1> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[J_mass_pre]  mass_pre_sign;
    int<lower=1> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[J_mass_post] mass_post_sign;
    int<lower=1> J_self_care_post; matrix[N, J_self_care_post] self_care_post_items; vector[J_self_care_post] self_care_post_sing;
}
// Transform Data Block 
transformed data {
   // Build the indices 
   vector[N] mass_pre  = composite(mass_pre_items, mass_pre_sign);
   vector[N] mass_post = composite(mass_post_items, mass_post_sign);
   vector[N] self_care_post = composite(self_care_post_items, self_care_post_sing);

   // Standatize all indices 
   vector[N] mass_pre_stand  = zscore(mass_pre);
   vector[N] mass_post_stand = zscore(mass_post);
   vector[N] self_care_post_stand = zscore(self_care_post); 
}
// Model Paramters block 
parameters{
    real alpha;
    real beta_treatment;
    real beta_mass_pre;
    real beta_self_care_post;
    real<lower=0.001> sigma;
}
// Model block
model{
    // Priors 
    alpha ~ normal(0, 0.5);
    beta_treatment      ~ normal(0, 0.5);
    beta_mass_pre       ~ normal(0.8, 0.3);
    beta_self_care_post ~ normal(0, 0.5);
    sigma ~ exponential(1);
    
    // Model Likelihood 
    mass_post_stand ~ normal(
        alpha 
        + beta_treatment      * treatment 
        + beta_mass_pre       * mass_pre_stand
        + beta_self_care_post * self_care_post_stand
        , sigma
    );
}
// Additional Calculatations 
generated quantities {
   // Linear Predictor 
   vector[N] mu = alpha 
        + beta_treatment      * treatment 
        + beta_mass_pre       * mass_pre_stand
        + beta_self_care_post * self_care_post_stand;

    // POSTERIOR PREDICTIVE GENERATION =======================================
    vector[N] mass_post_rep = normal_predictive_rng(mu, sigma);

    // MODEL FIT / INFORMATION CRITERIA ======================================
    vector[N] log_lik = normal_pointwise_loglik(mass_post_stand, mu, sigma);
    real R2 = bayes_R2_gaussian(mu, sigma);

    // EFFECT ESTIMATION =======================================================

    // Residual direct treatment effect
    real direct_treatment_effect     = beta_treatment;
    real direct_treatment_effect_std = standardized_effect(beta_treatment, sigma);

    // Self Care effect on metabolic burden
    real self_care_effect     = beta_self_care_post;
    real self_care_effect_std = standardized_effect(beta_self_care_post, sigma);

     // ADJUSTED EXPECTED OUTCOMES ===============================================

    // Expected metabolic burden at average self care + average baseline severity
    real adjusted_mean_control =
        alpha
        + beta_mass_pre * mean(mass_pre_stand)
        + beta_self_care_post * mean(self_care_post_stand);

    real adjusted_mean_treated =
        alpha
        + beta_treatment
        + beta_mass_pre * mean(mass_pre_stand)
        + beta_self_care_post * mean(self_care_post_stand);

    // Residual treatment effect after conditioning on self care 
    real residual_ATE = adjusted_mean_treated - adjusted_mean_control;

    // SELF CARE EFFECT INTERPRETATION ===========================================

    // More self care reduces metabolic burden
    int self_care_reduces_mass = effect_lt(beta_self_care_post, 0);

    // More self care worsens metabolic burden
    int self_care_increases_mass = effect_gt(beta_self_care_post, 0);

    // Self Care effect practically negligible
    int activity_effect_in_rope = in_rope(beta_self_care_post, -0.10, 0.10);

    // Strong clinically meaningful self care effect
    int activity_large_effect = effect_lt(beta_self_care_post, -0.30);

    // DIRECT TREATMENT EFFECT INTERPRETATION ====================================

    // Remaining treatment effect after conditioning on activity
    int treatment_reduces_mass = effect_lt(beta_treatment, 0);
    int treatment_effect_in_rope = in_rope(beta_treatment, -0.10, 0.10);

    // Strong residual treatment effect
    int treatment_large_reduction = effect_lt(beta_treatment, -0.30);

    // RESIDUAL DIAGNOSTICS =======================================================
    vector[N] raw_resid     = raw_residuals(mass_post_stand, mu);
    vector[N] pearson_resid = pearson_residuals(mass_post_stand, mu, sigma);

    // CALIBRATION DIAGNOSTICS ===================================================
    vector[N] pit = normal_pit(mass_post_stand, mu, sigma);

    // POSTERIOR PREDICTIVE CHECKS ===============================================
    int p_mean = ppc_indicator_mean(mass_post_stand, mass_post_rep);
    int p_sd   = ppc_indicator_sd(mass_post_stand, mass_post_rep);
    int p_max  = ppc_indicator_max(mass_post_stand, mass_post_rep);

}
