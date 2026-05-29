// Include Stan lib Function 
functions {
    #include "lib/utils.stanfunctions"
    #include "lib/composites.stanfunctions"
    #include "lib/diagnostics.stanfunctions"
}

// Input Data block 
data{

    int<lower=1> N; // Number of observations 
    vector[N] treatment;
    int<lower=0,upper=1> prior_only; // Prior Switch 
    vector[N] diet_pre;

    // Item blocks 
    int<lower=1> J_diet_post;  matrix[N, J_diet_post]  diet_post_items; vector[J_diet_post] diet_post_sign;
}
// Data Transformation Block 
transformed data {

   // Build the Indices 
   vector[N] diet_post  = composite(diet_post_items, diet_post_sign);
   
   // Standartize the Indices and the Diet post  
   vector[N] diet_pre_stand  = zscore(diet_pre);
   vector[N] diet_post_stand = zscore(diet_post);

    // For later back transformation 
   real diet_post_sd = sd(diet_post);
}
// Model Paramters 
parameters{
    real alpha;
    real beta_treatment;
    real beta_diet_pre;
    real<lower=0.001> sigma;
}
// Model Block 
model{
    // Priors 
    alpha          ~ normal(0, 0.5);
    beta_treatment ~ normal(0, 0.5);
    beta_diet_pre  ~ normal(0.8, 0.3);
    sigma          ~ exponential(1);
    
    // Model Likelihood
    if(prior_only == 0){
        diet_post_stand ~ normal(
            alpha 
            + beta_treatment * treatment 
            + beta_diet_pre  * diet_pre_stand
            ,sigma
        );
    }
}
// Aditional calculations 
generated quantities {
    
    // EXPECTED VALUES / LINEAR PREDICTOR ====================
    vector[N] mu =  alpha + beta_treatment * treatment + beta_diet_pre * diet_pre_stand;
    
    // POSTERIOR PREDICTIVE GENERATION 
    // Simulated replicated outcomes for PPCs
    vector[N] diet_post_rep = normal_predictive_rng(mu, sigma);
    
    // MODEL FIT / INFORMATION CRITERIA ========================
    // Pointwise log-likelihood for:
    vector[N] log_lik = normal_pointwise_loglik(diet_post_stand, mu, sigma);
    
    // Bayesian R2
    real R2 = bayes_R2_gaussian(mu, sigma);
    
    // TREATMENT EFFECT ESTIMATION =================================
    
    // Raw treatment effect
    real treatment_effect = beta_treatment;
    
    // Standardized effect size
    real treatment_effect_std = standardized_effect(beta_treatment, sigma);
    
    // Adjusted means (population-average predictions)
    real adjusted_mean_control = alpha + beta_diet_pre * mean(diet_post_stand);
    
    real adjusted_mean_treated =
        alpha
        + beta_treatment
        + beta_diet_pre * mean(diet_post_stand);
        
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
    vector[N] raw_resid = raw_residuals(diet_post_stand, mu);

    // Standardized / Pearson residuals
    vector[N] pearson_resid = pearson_residuals(diet_post_stand, mu, sigma);
    
    // CALIBRATION DIAGNOSTICS

    // PIT values should be Uniform(0,1)
    vector[N] pit = normal_pit(diet_post_stand, mu, sigma);

    // POSTERIOR PREDICTIVE CHECKS =====================================

    // Bayesian posterior predictive p-values
    int p_mean = ppc_indicator_mean(diet_post_stand, diet_post_rep);
    int p_sd   = ppc_indicator_sd(diet_post_stand, diet_post_rep);
    int p_max  = ppc_indicator_max(diet_post_stand, diet_post_rep);

}


