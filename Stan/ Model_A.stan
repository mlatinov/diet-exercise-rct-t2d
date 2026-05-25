
// Include Stan lib Function 
functions {
    #include "lib/utils.stanfunctions"
    #include "lib/composites.stanfunctions"
    #include "lib/diagnostics.stanfunctions"
}
// Input Data Block 
data{
    int<lower=0> N; // N Observations 
    vector[N] treatment;

    // Item Block per contruct
    int<lower=1> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[J_mass_post] mass_post_sign;
    int<lower=1> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[J_mass_pre] mass_pre_sign;
}
// Transformation Data Block 
transformed data {
    // Build the Indices
   vector[N] mass_post = composite(mass_post_items, mass_post_sign);
   vector[N] mass_pre  = composite(mass_pre_items, mass_pre_sign);

   // Standatize the Indices
   vector[N] mass_post_stand = zscore(mass_post);
   vector[N] mass_pre_stand  = zscore(mass_pre);

   // For later back transformation 
   real mass_post_composite_sd = sd(mass_post);
}
// Model Parameters Block
parameters{
    real beta_treatment;
    real beta_mass_pre;
    real alpha;
    real<lower=0> sigma;
}
// Model Block 
model{
    // Priors
    alpha          ~ normal(0, 0.5);
    beta_treatment ~ normal(0, 0.5);
    beta_mass_pre  ~ normal(0.8, 0.3);
    sigma          ~ exponential(1);

    // Model Likelihood 
    mass_post_stand ~ normal(
        alpha 
        + beta_treatment * treatment 
        + beta_mass_pre  * mass_pre_stand
        ,sigma
    );
}
// Additional Calculations 
generated quantities {

    // EXPECTED VALUES / LINEAR PREDICTOR ====================
    vector[N] mu = 
        alpha
        + beta_treatment * treatment
        + beta_mass_pre  * mass_pre_stand;

    // POSTERIOR PREDICTIVE GENERATION 
    // Simulated replicated outcomes for PPCs
    vector[N] mass_post_rep = normal_predictive_rng(mu, sigma);

    // MODEL FIT / INFORMATION CRITERIA ========================

    // Pointwise log-likelihood for:
    vector[N] log_lik = normal_pointwise_loglik(mass_post_stand, mu, sigma);

    // Bayesian R²
    real R2 = bayes_R2_gaussian(mu, sigma);

    // TREATMENT EFFECT ESTIMATION =================================

    // Raw treatment effect
    real treatment_effect = beta_treatment;

    // Standardized effect size
    real treatment_effect_std = standardized_effect(beta_treatment, sigma);

    // Adjusted means (population-average predictions)
    real adjusted_mean_control = alpha + beta_mass_pre * mean(mass_pre_stand);

    real adjusted_mean_treated =
        alpha
        + beta_treatment
        + beta_mass_pre * mean(mass_pre_stand);

    // Average Treatment Effect
    real ATE = adjusted_mean_treated - adjusted_mean_control;

    // DIRECTIONAL / PROBABILITY STATEMENTS =========================

    // Probability treatment reduces metabolic burden
    int treatment_reduces_mass = effect_lt(beta_treatment, 0);

    // Probability treatment increases metabolic burden
    int treatment_increases_mass = effect_gt(beta_treatment, 0);

    // Probability effect practically negligible
    int treatment_in_rope = in_rope(beta_treatment, -0.10, 0.10);

    // Strong clinically meaningful reduction
    int treatment_large_reduction = effect_lt(beta_treatment, -0.30);

    // RESIDUAL DIAGNOSTICS ===========================================

    // Raw residuals
    vector[N] raw_resid = raw_residuals(mass_post_stand, mu);

    // Standardized / Pearson residuals
    vector[N] pearson_resid = pearson_residuals(mass_post_stand, mu, sigma);
    
    // CALIBRATION DIAGNOSTICS

    // PIT values should be Uniform(0,1)
    vector[N] pit = normal_pit(mass_post_stand, mu, sigma);

    // POSTERIOR PREDICTIVE CHECKS =====================================

    // Bayesian posterior predictive p-values
    int p_mean = ppc_indicator_mean(mass_post_stand, mass_post_rep);
    int p_sd   = ppc_indicator_sd(mass_post_stand, mass_post_rep);
    int p_max  = ppc_indicator_max(mass_post_stand, mass_post_rep);

}
