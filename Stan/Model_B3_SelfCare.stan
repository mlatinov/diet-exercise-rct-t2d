// Include Stan lib Function 
functions {
    #include "lib/utils.stanfunctions"
    #include "lib/composites.stanfunctions"
    #include "lib/diagnostics.stanfunctions"
}

// Input data Block 
data{
    int<lower=1> N;
    vector[N] treatment;
    vector[N] self_care_pre;

    // Indices Blocks 
    int<lower=1> J_self_care_post; 
    matrix[N, J_self_care_post] self_care_post_items; 
    vector[J_self_care_post] self_care_post_sing;
}

// Data Transformation block 
transformed data {
   // Build the Indice 
   vector[N] self_care = composite(self_care_post_items, self_care_post_sing);

   // Standartize the Self Care Composite and self care pre
   vector[N] self_care_stand     = zscore(self_care);
   vector[N] self_care_pre_stand = zscore(self_care_pre);
   
   // For later back transformation 
   real self_care_post_composite_sd = sd(self_care);  
}

// Model Paramters 
parameters{
    real alpha;
    real beta_treatment;
    real beta_self_care_pre;
    real<lower=0.001> sigma;
}

// Model Block 
model{
    // Priors
    alpha              ~ normal(0, 0.5);
    beta_treatment     ~ normal(0, 0.5);
    beta_self_care_pre ~ normal(0.8, 0.3);
    sigma ~ exponential(1);

    // Model Likelihood
    self_care_stand ~ normal(
        alpha 
        + beta_treatment     * treatment 
        + beta_self_care_pre * self_care_pre_stand
        , sigma
    );
}

// Additional Calculations 
generated quantities {

    // EXPECTED VALUES / LINEAR PREDICTOR =========================
    vector[N] mu = alpha + beta_treatment * treatment + beta_self_care_pre * self_care_pre_stand;

    // POSTERIOR PREDICTIVE GENERATION ============================
    // Simulated replicated outcomes for PPCs
    vector[N] self_care_rep = normal_predictive_rng(mu, sigma);
    
    // MODEL FIT / INFORMATION CRITERIA ===========================
    // Pointwise log-likelihood for:
    vector[N] log_lik = normal_pointwise_loglik(self_care_stand, mu, sigma);
    
    // Bayesian R2
    real R2 = bayes_R2_gaussian(mu, sigma);
    
    // TREATMENT EFFECT ESTIMATION =================================
    
    // Raw treatment effect
    real treatment_effect = beta_treatment;
    
    // Standardized effect size
    real treatment_effect_std = standardized_effect(beta_treatment, sigma);
    
    // Adjusted means (population-average predictions)
    real adjusted_mean_control = alpha + beta_self_care_pre * mean(self_care_stand);
    
    real adjusted_mean_treated =
        alpha
        + beta_treatment
        + beta_self_care_pre * mean(self_care_pre_stand);
        
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
    vector[N] raw_resid = raw_residuals(self_care_stand, mu);

    // Standardized / Pearson residuals
    vector[N] pearson_resid = pearson_residuals(self_care_stand, mu, sigma);
    
    // CALIBRATION DIAGNOSTICS ===========================================

    // PIT values should be Uniform(0,1)
    vector[N] pit = normal_pit(self_care_stand, mu, sigma);

    // POSTERIOR PREDICTIVE CHECKS =====================================

    // Bayesian posterior predictive p-values
    int p_mean = ppc_indicator_mean(self_care_stand, self_care_rep);
    int p_sd   = ppc_indicator_sd(self_care_stand, self_care_rep);
    int p_max  = ppc_indicator_max(self_care_stand, self_care_rep);
}
