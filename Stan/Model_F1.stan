
// Include Stan lib Function 
functions {
    #include "lib/utils.stanfunctions"
    #include "lib/composites.stanfunctions"
    #include "lib/diagnostics.stanfunctions"
    #include "lib/mediation.stanfunctions"
}

// Data Input Block 
data{
    int<lower=1> N;
    vector[N] treatment;

    // =============== Composites Building Block ======================
    // Mass
    int<lower=1> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[J_mass_pre] mass_pre_sign;
    int<lower=1> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[J_mass_post] mass_post_sign;

    // Diet 
    int<lower=1> J_diet_post; matrix[N, J_diet_post] diet_post_items; vector[J_diet_post] diet_post_sign;
    vector[N] diet_pre;

    // Self Care 
    int<lower=1> J_self_care_post; matrix[N, J_self_care_post] self_care_post_items; vector[J_self_care_post] self_care_post_sing;
    vector[N] self_care_pre;

    // Activity
    vector[N] intensity_pre;
    vector[N] duration_pre;
    vector[N] frequency_pre;
    vector[N] intensity_post;
    vector[N] duration_post;
    vector[N] frequency_post;
}
// Data Transformation Block 
transformed data {
   // ============== Build the Composite Indices ====================
   
   // Mass 
   vector[N] mass_pre  = composite(mass_pre_items, mass_pre_sign);
   vector[N] mass_post = composite(mass_post_items, mass_post_sign);

   // Diet 
   vector[N] diet_post = composite(diet_post_items, diet_post_sign);

   // Self Care 
   vector[N] self_care_post = composite(self_care_post_items, self_care_post_sing);

   // Activity 
   vector[N] activity_pre  = activity_index(intensity_pre, duration_pre, frequency_pre);
   vector[N] activity_post = activity_index(intensity_post, duration_post, frequency_post);

   // ====================  Z score Standartize the Indices =================
   // Mass
   vector[N] mass_pre_stand  = zscore(mass_pre);
   vector[N] mass_post_stand = zscore(mass_post);
   
   // Diet 
   vector[N] diet_post_stand = zscore(diet_post);
   vector[N] diet_pre_stand  = zscore(diet_pre);

   // Self Care 
   vector[N] self_care_post_stand  = zscore(self_care_post);
   vector[N] self_care_pre_stand   = zscore(self_care_pre);
}
// Models Paramters 
parameters{
    // Diet Submodel Parameters 
    real alpha_diet;
    real beta_treatment_diet;
    real beta_diet_pre;
    real<lower=0.001> sigma_diet;

    // Self Care Submodel Paramters 
    real alpha_self_care;
    real beta_treatment_self_care;
    real beta_self_care_pre;
    real<lower=0.001> sigma_self_care;

    // Activity Submodel Paramters 
    real alpha_activity;
    real beta_treatment_activity;
    real beta_activity_pre;
    real<lower=0.001> sigma_activity;

    // Mass Post final Equation Paramters
    real alpha_mass;
    real beta_treatment_mass;
    real beta_mass_pre;
    real beta_diet_mass;
    real beta_activity_mass;
    real beta_selfcare_mass;
    real<lower=0.001>sigma_mass;  
}
// Model 
model{
    // ========  Diet Submodel Priors and Likelihood ==============
    alpha_diet ~ normal(0, 0.5);
    beta_treatment_diet ~ normal(0, 0.5);
    beta_diet_pre       ~ normal(0, 0.5);
    sigma_diet ~ exponential(1);

    diet_post_stand ~ normal(
        alpha_diet 
        + beta_treatment_diet * treatment
        + beta_diet_pre       * diet_pre_stand
        , sigma_diet
    );

    // ========== Self Care Submodel Priors and Likelihood ========================
    alpha_self_care ~ normal(0, 0.5);
    beta_treatment_self_care ~ normal(0, 0.5);
    beta_self_care_pre       ~ normal(0, 0.5);
    sigma_self_care ~ exponential(1);
    
    self_care_post_stand ~ normal(
        alpha_self_care 
        + beta_treatment_self_care * treatment
        + beta_self_care_pre       * self_care_pre_stand
        ,sigma_self_care
    );

    // ========== Activity Submodel Priors and Likelihood ==========================
    alpha_activity ~ normal(0, 0.5);
    beta_treatment_activity ~ normal(0, 0.5);
    beta_activity_pre       ~ normal(0, 0.5);
    sigma_activity ~ exponential(1);

    activity_post ~ normal(
        alpha_activity 
        + beta_treatment_activity * treatment
        + beta_activity_pre       * activity_pre
        ,sigma_activity
    );

    // ==============  Mass Model Final Equation Priors and Likelihood ================
    alpha_mass ~ normal(0, 0.5); 
    beta_treatment_mass ~ normal(0, 0.5);
    beta_mass_pre ~  normal(0.8, 0.3);
    beta_diet_mass ~ normal(0, 0.5);
    beta_activity_mass ~ normal(0, 0.5);
    beta_selfcare_mass ~ normal(0, 0.5); 
    sigma_mass ~ exponential(1);

    mass_post_stand ~ normal(
        alpha_mass 
        + beta_treatment_mass * treatment
        + beta_mass_pre       * mass_pre_stand
        + beta_diet_mass      * diet_post_stand
        + beta_activity_mass  * activity_post
        + beta_selfcare_mass  * self_care_post
        ,sigma_mass
    );
}
// Additional Calculations 
generated quantities {

    // LINEAR PREDICTOR (FINAL MASS EQUATION) =====================================
    vector[N] mu =
        alpha_mass
        + beta_treatment_mass * treatment
        + beta_mass_pre       * mass_pre_stand
        + beta_diet_mass      * diet_post_stand
        + beta_activity_mass  * activity_post
        + beta_selfcare_mass  * self_care_post_stand;

    // POSTERIOR PREDICTIVE GENERATION ===========================================
    vector[N] mass_post_rep = normal_predictive_rng(mu, sigma_mass);

    // MODEL FIT / INFORMATION CRITERIA ==========================================
    vector[N] log_lik = normal_pointwise_loglik(mass_post_stand, mu, sigma_mass);
    real R2           = bayes_R2_gaussian(mu, sigma_mass);

    // DIRECT EFFECT ============================================================

    // Remaining treatment effect after accounting for mediators
    real direct_effect     = beta_treatment_mass;
    real direct_effect_std = standardized_effect(beta_treatment_mass, sigma_mass);

    // INDIVIDUAL MEDIATOR INDIRECT EFFECTS (NIEs) ==============================
    
    // Treatment -> Diet -> Mass
    real NIE_diet = beta_treatment_diet* beta_diet_mass;

    // Treatment -> Activity -> Mass
    real NIE_activity = beta_treatment_activity * beta_activity_mass;

    // Treatment -> SelfCare -> Mass
    real NIE_selfcare = beta_treatment_self_care * beta_selfcare_mass;

    // TOTAL INDIRECT EFFECT ===================================================
    real total_indirect_effect = NIE_diet+ NIE_activity + NIE_selfcare;

    // TOTAL EFFECT ============================================================
    real total_effect = direct_effect + total_indirect_effect;

    // PROPORTION MEDIATED =====================================================
    real proportion_mediated_total    = proportion_mediated(total_effect, total_indirect_effect);
    real proportion_mediated_diet     = proportion_mediated(total_effect, NIE_diet);
    real proportion_mediated_activity = proportion_mediated(total_effect, NIE_activity);
    real proportion_mediated_selfcare = proportion_mediated(total_effect, NIE_selfcare);
    
    // MEDIATION DECOMPOSITION ========================================================
    vector[4] mediation_diet     = mediation_linear(beta_treatment_diet, beta_diet_mass, beta_treatment_mass);
    vector[4] mediation_activity = mediation_linear(beta_treatment_activity, beta_activity_mass, beta_treatment_mass);
    vector[4] mediation_selfcare = mediation_linear(beta_treatment_self_care, beta_selfcare_mass, beta_treatment_mass);

    // DIRECTIONAL / PROBABILITY STATEMENTS ============================================

    // Total treatment effect 
    int total_effect_reduces_mass = effect_lt(total_effect, 0);
    int total_effect_in_rope      = in_rope(total_effect, -0.10, 0.10);

    // Direct effect 
    int direct_effect_reduces_mass = effect_lt(direct_effect, 0);
    int direct_effect_in_rope      = in_rope(direct_effect, -0.10, 0.10);

    // Diet mediation 
    int diet_mediation_present = effect_lt(NIE_diet, 0);
    int diet_mediation_in_rope = in_rope(NIE_diet, -0.05, 0.05);

    // Activity mediation 
    int activity_mediation_present = effect_lt(NIE_activity, 0);
    int activity_mediation_in_rope = in_rope(NIE_activity, -0.05, 0.05);

    // Self-care mediation 
    int selfcare_mediation_present = effect_lt(NIE_selfcare, 0);
    int selfcare_mediation_in_rope = in_rope(NIE_selfcare, -0.05, 0.05);

    // RESIDUAL DIAGNOSTICS =============================================
    vector[N] raw_resid     = raw_residuals(mass_post_stand, mu);
    vector[N] pearson_resid = pearson_residuals(mass_post_stand, mu, sigma_mass);
    
    // CALIBRATION DIAGNOSTICS ==============================================
    vector[N] pit = normal_pit(mass_post_stand, mu, sigma_mass);

    // POSTERIOR PREDICTIVE CHECKS ==========================================
    int p_mean = ppc_indicator_mean(mass_post_stand, mass_post_rep);
    int p_sd   = ppc_indicator_sd(mass_post_stand, mass_post_rep);
    int p_max  = ppc_indicator_max(mass_post_stand, mass_post_rep);

}

