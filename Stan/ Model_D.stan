// Include Stan lib Function 
functions {
    #include "lib/utils.stanfunctions"
    #include "lib/composites.stanfunctions"
    #include "lib/diagnostics.stanfunctions"
}
// Data Input block 
data{
    int<lower=1> N;
    vector[N] treatment;

    // Composite Building Block 

    // Mass
    int<lower=1> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[J_mass_post] mass_post_sign;
    int<lower=1> J_mass_pre ; matrix[N, J_mass_pre ] mass_pre_items ; vector[J_mass_pre ] mass_pre_sign;
    
    // Diet 
    int<lower=1> J_diet_post; matrix[N, J_diet_post] diet_post_items; vector[J_diet_post] diet_post_sign;
    
    // SelfCare
    int<lower=2> J_self_care_post; matrix[N, J_self_care_post] self_care_post_items; vector[J_self_care_post] self_care_post_sing;
    
    // Activity
    vector[N] intensity_post;
    vector[N] duration_post ;
    vector[N] frequency_post;
}
// Data Transformation block 
transformed data {
   // Build the Indices
   vector[N] mass_post  = composite(mass_post_items, mass_post_sign);
   vector[N] mass_pre   = composite(mass_pre_items , mass_pre_sign);
   vector[N] diet_post  = composite(diet_post_items, diet_post_sign);
   vector[N] self_care_post = composite(self_care_post_items, self_care_post_sing);
   vector[N] activity_post  = activity_index(intensity_post, duration_post, frequency_post);

   // Standartize indices
   vector[N] mass_post_stand = zscore(mass_post);
   vector[N] mass_pre_stand  = zscore(mass_pre);
   vector[N] diet_post_stand = zscore(diet_post);
   vector[N] self_care_post_stand  = zscore(self_care_post);
}
// Model Paramters 
parameters{
    real alpha;
    real beta_treatment;
    real beta_mass_pre;
    real beta_diet_post;
    real beta_self_care_post;
    real beta_activity_post;
    real<lower=0.001> sigma;
}
// Model Block 
model{
    alpha ~ normal(0, 0.5);
    beta_treatment ~ normal(0, 0.5);
    beta_mass_pre  ~ normal(0.8, 0.3);
    beta_diet_post ~ normal(0, 0.5);
    beta_self_care_post ~ normal(0, 0.5);
    beta_activity_post  ~ normal(0, 0.5);
    sigma ~ exponential(1);
    
    // Model Likelihood 
    mass_post_stand ~ normal(
        alpha 
        + beta_treatment      * treatment
        + beta_self_care_post * self_care_post_stand
        + beta_diet_post      * diet_post_stand
        + beta_self_care_post * self_care_post_stand
        + beta_activity_post  * activity_post
        ,sigma
    );
}
// Additional Calculations 
generated quantities {

    // EXPECTED VALUES / LINEAR PREDICTOR =======================================
    vector[N] mu =
        alpha
        + beta_treatment      * treatment
        + beta_mass_pre       * mass_pre_stand
        + beta_diet_post      * diet_post_stand
        + beta_self_care_post * self_care_post_stand
        + beta_activity_post  * activity_post;

    // POSTERIOR PREDICTIVE GENERATION =========================================
    vector[N] mass_post_rep = normal_predictive_rng(mu, sigma);

    // MODEL FIT / INFORMATION CRITERIA =======================================
    vector[N] log_lik = normal_pointwise_loglik(mass_post_stand, mu, sigma);
    real R2           = bayes_R2_gaussian(mu, sigma);

    // DIRECT TREATMENT EFFECT ================================================

    // Remaining treatment effect after adjusting for ALL mediators
    real direct_treatment_effect     = beta_treatment;
    real direct_treatment_effect_std = standardized_effect(beta_treatment, sigma);

    // MEDIATOR EFFECTS ========================================================
    real diet_effect     = beta_diet_post;
    real diet_effect_std = standardized_effect(beta_diet_post, sigma);
    
    real selfcare_effect     = beta_self_care_post;
    real selfcare_effect_std = standardized_effect(beta_self_care_post, sigma);

    real activity_effect     = beta_activity_post;
    real activity_effect_std = standardized_effect(beta_activity_post, sigma);

    // ADJUSTED EXPECTED OUTCOMES ==============================================
    real adjusted_mean_control =
        alpha
        + beta_mass_pre       * mean(mass_pre_stand)
        + beta_diet_post      * mean(diet_post_stand)
        + beta_self_care_post * mean(self_care_post_stand)
        + beta_activity_post  * mean(activity_post);

    real adjusted_mean_treated =
        alpha
        + beta_treatment
        + beta_mass_pre       * mean(mass_pre_stand)
        + beta_diet_post      * mean(diet_post_stand)
        + beta_self_care_post * mean(self_care_post_stand)
        + beta_activity_post  * mean(activity_post);

    // Residual treatment effect after all mediators
    real residual_ATE = adjusted_mean_treated - adjusted_mean_control;

    // DIRECTIONAL / PROBABILITY STATEMENTS ====================================
    
    // Treatment effects
    int treatment_reduces_mass   = effect_lt(beta_treatment, 0);
    int treatment_increases_mass = effect_gt(beta_treatment, 0);
    int treatment_effect_in_rope = in_rope(beta_treatment, -0.10, 0.10);

    // Diet effect 
    int diet_reduces_mass   = effect_lt(beta_diet_post, 0);
    int diet_effect_in_rope = in_rope(beta_diet_post, -0.10, 0.10);

    // Self-care effect 
    int selfcare_reduces_mass   = effect_lt(beta_self_care_post, 0);
    int selfcare_effect_in_rope = in_rope(beta_self_care_post, -0.10, 0.10);

    // Activity effect 
    int activity_reduces_mass   = effect_lt(beta_activity_post, 0);
    int activity_effect_in_rope = in_rope(beta_activity_post, -0.10, 0.10);

    // RESIDUAL DIAGNOSTICS ====================================================
    vector[N] raw_resid     = raw_residuals(mass_post_stand, mu);
    vector[N] pearson_resid = pearson_residuals(mass_post_stand, mu, sigma);

    // CALIBRATION DIAGNOSTICS =================================================
    vector[N] pit = normal_pit(mass_post_stand, mu, sigma);
    
    // POSTERIOR PREDICTIVE CHECKS
    int p_mean = ppc_indicator_mean(mass_post_stand, mass_post_rep);
    int p_sd   = ppc_indicator_sd(mass_post_stand, mass_post_rep);
    int p_max  = ppc_indicator_max(mass_post_stand, mass_post_rep);

}

