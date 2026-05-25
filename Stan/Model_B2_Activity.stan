// Include Stan lib Function 
functions {
    #include "lib/utils.stanfunctions"
    #include "lib/composites.stanfunctions"
    #include "lib/diagnostics.stanfunctions"
}

// Input data block 
data{
    int<lower=1> N;
    vector[N] treatment;

    // Indices Building blocks 
    vector<lower=0>[N] intensity_pre;
    vector<lower=0>[N] duration_pre;
    vector<lower=0>[N] frequency_pre;

    vector<lower=0>[N] intensity_post;
    vector<lower=0>[N] duration_post;
    vector<lower=0>[N] frequency_post;
}
// Input data Trasnformation Block 
transformed data {
   // Build the Activity Mets Indices
   vector[N] activity_pre  = activity_index(intensity_pre, duration_pre,frequency_pre);
   vector[N] activity_post = activity_index(intensity_post, duration_post, frequency_post); 
}
// Model paramters 
parameters{
    real alpha;
    real beta_treatment;
    real beta_activity_pre;
    real<lower=0.001> sigma;
}
// Model Block 
model{
    // Priors 
    alpha             ~ normal(0, 0.5);
    beta_treatment    ~ normal(0, 0.5);
    beta_activity_pre ~ normal(0.8, 0.3);
    sigma             ~ exponential(1); 

    // Model Likelihood 
    activity_post ~ normal(
        alpha 
        + beta_treatment    * treatment 
        + beta_activity_pre * activity_pre
        ,sigma
    );
}
// Aditional Calculations 
generated quantities {
    
    // EXPECTED VALUES / LINEAR PREDICTOR ===========================
    vector[N] mu = alpha + beta_treatment * treatment + beta_activity_pre * activity_pre;

    // POSTERIOR PREDICTIVE GENERATION ==============================
    // Simulated replicated outcomes for PPCs
    vector[N] activity_post_rep = normal_predictive_rng(mu, sigma);
    
    // MODEL FIT / INFORMATION CRITERIA =============================
    // Pointwise log-likelihood for:
    vector[N] log_lik = normal_pointwise_loglik(activity_post, mu, sigma);
    
    // Bayesian R2
    real R2 = bayes_R2_gaussian(mu, sigma);
    
    // TREATMENT EFFECT ESTIMATION =================================
    
    // Raw treatment effect
    real treatment_effect = beta_treatment;
    
    // Standardized effect size
    real treatment_effect_std = standardized_effect(beta_treatment, sigma);
    
    // Adjusted means (population-average predictions)
    real adjusted_mean_control = alpha + beta_activity_pre * mean(activity_post);
    
    real adjusted_mean_treated =
        alpha
        + beta_treatment
        + beta_activity_pre * mean(activity_post);
        
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
    vector[N] raw_resid = raw_residuals(activity_post, mu);

    // Standardized / Pearson residuals
    vector[N] pearson_resid = pearson_residuals(activity_post, mu, sigma);
    
    // CALIBRATION DIAGNOSTICS ==========================================

    // PIT values should be Uniform(0,1)
    vector[N] pit = normal_pit(activity_post, mu, sigma);

    // POSTERIOR PREDICTIVE CHECKS =====================================

    // Bayesian posterior predictive p-values
    int p_mean = ppc_indicator_mean(activity_post, activity_post_rep);
    int p_sd   = ppc_indicator_sd(activity_post, activity_post_rep);
    int p_max  = ppc_indicator_max(activity_post, activity_post_rep);

}
