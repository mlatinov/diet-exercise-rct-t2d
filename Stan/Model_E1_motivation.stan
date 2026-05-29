// Include Stan lib Function 
functions {
    #include "lib/utils.stanfunctions"
    #include "lib/composites.stanfunctions"
    #include "lib/diagnostics.stanfunctions"
}

// Data Input Block 
data{
    int<lower=1> N;
    vector[N] treatment;
    int<lower=0,upper=1> prior_only;
    vector[N] barrier_motivation;

    // Composite Building blocks ============================================================= 
    
    // Mass
    int<lower=1> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[J_mass_pre] mass_pre_sign;
    int<lower=1> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[J_mass_post] mass_post_sign;
    
    // Diet
    int<lower=1> J_diet_post; matrix[N, J_diet_post] diet_post_items; vector[J_diet_post] diet_post_sign;
    
    // Self Care
    int<lower=1> J_self_care_post; matrix[N, J_self_care_post] self_care_post_items; vector[J_self_care_post] self_care_post_sing;
    
    // Activity
    vector[N] intensity_post;
    vector[N] duration_post;
    vector[N] frequency_post;
}
// Transform Data Block 
transformed data {
   // Build Composite Indices
   vector[N] mass_pre  = composite(mass_pre_items, mass_pre_sign);
   vector[N] mass_post = composite(mass_post_items, mass_post_sign);
   vector[N] diet_post = composite(diet_post_items, diet_post_sign);
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
    real beta_barrier_motivation;
    real beta_treatment_barrier;
    real beta_diet_post;
    real beta_self_care_post;
    real beta_activity_post;
    real<lower=0.001> sigma;
}
// Model Block 
model{
    // Priors 
    alpha ~ normal(0, 0.5);
    beta_treatment ~ normal(0, 0.5);
    beta_mass_pre  ~ normal(0.8, 0.5);
    beta_diet_post ~ normal(0, 0.5);
    beta_barrier_motivation ~ normal(0, 0.5);
    beta_treatment_barrier  ~ normal(0, 0.3);
    beta_self_care_post     ~ normal(0, 0.5);
    beta_activity_post      ~ normal(0, 0.5);
    sigma ~ exponential(1);

    // Model Likelihood 
    if(prior_only == 0){
        mass_post_stand ~ normal(
            alpha
            + beta_treatment * treatment
            + beta_barrier_motivation * barrier_motivation
            + beta_treatment_barrier  * treatment .* barrier_motivation
            + beta_mass_pre  * mass_pre_stand
            + beta_diet_post * diet_post_stand
            + beta_self_care_post * self_care_post_stand
            + beta_activity_post  * activity_post
            ,sigma
        );
    }
}
// Additonal Calculations 
generated quantities {

    // EXPECTED VALUES / LINEAR PREDICTOR ========================================
    vector[N] mu =
        alpha
        + beta_treatment * treatment
        + beta_barrier_motivation * barrier_motivation
        + beta_treatment_barrier  * (treatment .* barrier_motivation)
        + beta_mass_pre * mass_pre_stand
        + beta_diet_post * diet_post_stand
        + beta_self_care_post * self_care_post_stand
        + beta_activity_post * activity_post;

    // POSTERIOR PREDICTIVE GENERATION ============================================
    vector[N] mass_post_rep = normal_predictive_rng(mu, sigma);

    // MODEL FIT / INFORMATION CRITERIA ===========================================
    vector[N] log_lik = normal_pointwise_loglik(mass_post_stand, mu, sigma);
    real R2           = bayes_R2_gaussian(mu, sigma);

    // MAIN EFFECTS ===============================================================
    real treatment_effect     = beta_treatment;
    real treatment_effect_std = standardized_effect(beta_treatment, sigma);

    real barrier_effect     = beta_barrier_motivation;
    real barrier_effect_std = standardized_effect(beta_barrier_motivation, sigma);

    // MODERATION EFFECT ==========================================================
    real interaction_effect     = beta_treatment_barrier;
    real interaction_effect_std = standardized_effect(beta_treatment_barrier, sigma);

    // SIMPLE INTERPRETATION OF INTERACTION ========================================

    // Treatment effect when barrier = 0
    real treatment_effect_low_barrier = beta_treatment;

    // Treatment effect when barrier = 1
    real treatment_effect_high_barrier = beta_treatment + beta_treatment_barrier;

    // Difference between high vs low barrier treatment effects
    real moderation_difference = treatment_effect_high_barrier - treatment_effect_low_barrier;

    // DIRECTIONAL / PROBABILITY STATEMENTS =============================================

    // Interaction effect 

    // Higher barriers REDUCE treatment effectiveness
    int interaction_reduces_treatment = effect_gt(beta_treatment_barrier, 0);

    // Higher barriers IMPROVE treatment effectiveness
    int interaction_improves_treatment = effect_lt(beta_treatment_barrier, 0);

    // Practically negligible interaction
    int interaction_in_rope = in_rope(beta_treatment_barrier, -0.10, 0.10);

    // Strong moderation effect
    int strong_interaction = effect_gt(beta_treatment_barrier, 0.30);
    
    // MEDIATOR EFFECTS ================================================================
    real diet_effect     = beta_diet_post;
    real selfcare_effect = beta_self_care_post;
    real activity_effect = beta_activity_post;

    // ADJUSTED EXPECTED OUTCOMES ======================================================

    // Expected outcome for average participant under control
    real adjusted_mean_control =
        alpha
        + beta_barrier_motivation * mean(barrier_motivation)
        + beta_mass_pre  * mean(mass_pre_stand)
        + beta_diet_post * mean(diet_post_stand)
        + beta_self_care_post * mean(self_care_post_stand)
        + beta_activity_post  * mean(activity_post);

    // Expected outcome for average participant under treatment
    real adjusted_mean_treated =
        alpha
        + beta_treatment
        + beta_barrier_motivation * mean(barrier_motivation)
        + beta_treatment_barrier  * mean(barrier_motivation)
        + beta_mass_pre  * mean(mass_pre_stand)
        + beta_diet_post * mean(diet_post_stand)
        + beta_self_care_post * mean(self_care_post_stand)
        + beta_activity_post  * mean(activity_post);

    real conditional_ATE = adjusted_mean_treated - adjusted_mean_control;

    // RESIDUAL DIAGNOSTICS ====================================================
    vector[N] raw_resid     = raw_residuals(mass_post_stand, mu);
    vector[N] pearson_resid = pearson_residuals(mass_post_stand, mu, sigma);

    // CALIBRATION DIAGNOSTICS ==================================================
    vector[N] pit = normal_pit(mass_post_stand, mu, sigma);

    // POSTERIOR PREDICTIVE CHECKS
    int p_mean = ppc_indicator_mean(mass_post_stand, mass_post_rep);
    int p_sd   = ppc_indicator_sd(mass_post_stand, mass_post_rep);
    int p_max  = ppc_indicator_max(mass_post_stand, mass_post_rep);

}



