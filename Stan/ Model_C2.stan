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

    // Composit Building blocks ==============================================================================

    // Mass 
    int<lower=1> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[J_mass_pre] mass_pre_sign;
    int<lower=1> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[J_mass_post] mass_post_sign;

    // Activity Post 
    vector[N] intensity_post;
    vector[N] duration_post;
    vector[N] frequency_post;
}
// Data Transformation Block 
transformed data {
   // Build each indices
   vector[N] mass_pre  = composite(mass_pre_items, mass_pre_sign);
   vector[N] mass_post = composite(mass_post_items, mass_post_sign);
   vector[N] activity_post = activity_index(intensity_post, duration_post, frequency_post);

   // Standartize the Mess Composites 
   vector[N] mass_pre_stand  = zscore(mass_pre);
   vector[N] mass_post_stand = zscore(mass_post);
   
   // For later back transformation 
   real mass_post_composite_sd = sd(mass_post);
}
// Model Parameters 
parameters{
    real alpha;
    real beta_treatment;
    real beta_mass_pre;
    real beta_activity_post;
    real<lower=0.001> sigma;
}
// Model Block 
model{
    // Priors 
    alpha ~ normal(0, 0.5);
    beta_treatment     ~ normal(0, 0.5);
    beta_activity_post ~ normal(0, 0.5);
    beta_mass_pre      ~ normal(0.8, 0.3);
    sigma ~ exponential(1); 

    // Model Likelihood 
    if(prior_only == 0){
        mass_post_stand ~ normal(  
            alpha 
            + beta_treatment     * treatment 
            + beta_activity_post * activity_post 
            + beta_mass_pre      * mass_pre_stand
            ,sigma
        );
    }
}
// Additional Calculations 
generated quantities {

    // EXPECTED VALUES / LINEAR PREDICTOR ======================================
    vector[N] mu = 
        alpha
        + beta_treatment     * treatment
        + beta_activity_post * activity_post
        + beta_mass_pre      * mass_pre_stand;

    // POSTERIOR PREDICTIVE GENERATION =======================================
    vector[N] mass_post_rep = normal_predictive_rng(mu, sigma);

    // MODEL FIT / INFORMATION CRITERIA ======================================
    vector[N] log_lik = normal_pointwise_loglik(mass_post_stand, mu, sigma);
    real R2 = bayes_R2_gaussian(mu, sigma);

    // EFFECT ESTIMATION =======================================================

    // Residual direct treatment effect
    real direct_treatment_effect = beta_treatment;
    real direct_treatment_effect_std = standardized_effect(beta_treatment, sigma);

    // Activity effect on metabolic burden
    real activity_effect = beta_activity_post;
    real activity_effect_std = standardized_effect(beta_activity_post, sigma);

    // ADJUSTED EXPECTED OUTCOMES ===============================================

    // Expected metabolic burden at average activity + average baseline severity
    real adjusted_mean_control =
        alpha
        + beta_mass_pre * mean(mass_pre_stand)
        + beta_activity_post * mean(activity_post);

    real adjusted_mean_treated =
        alpha
        + beta_treatment
        + beta_mass_pre * mean(mass_pre_stand)
        + beta_activity_post * mean(activity_post);

    // Residual treatment effect after conditioning on activity
    real residual_ATE = adjusted_mean_treated - adjusted_mean_control;
    
    // ACTIVITY EFFECT INTERPRETATION ===========================================

    // More activity reduces metabolic burden
    int activity_reduces_mass = effect_lt(beta_activity_post, 0);

    // More activity worsens metabolic burden
    int activity_increases_mass = effect_gt(beta_activity_post, 0);

    // Activity effect practically negligible
    int activity_effect_in_rope = in_rope(beta_activity_post, -0.10, 0.10);

    // Strong clinically meaningful activity effect
    int activity_large_effect = effect_lt(beta_activity_post, -0.30);

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


