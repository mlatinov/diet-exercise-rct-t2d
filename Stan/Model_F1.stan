functions {
    // Z Score Standartization 
    vector zscore(vector x) {
        return (x - mean(x)) / sd(x);
    }
    // Build an equal-weight composite from a block of items with per-item signs.
    vector composite(matrix items, vector sign) {
        int N = rows(items);
        int J = cols(items);
        vector[N] acc = rep_vector(0, N);
        for (j in 1:J) acc += sign[j] * zscore(items[, j]);
        return acc / J;
  }
      // MET-style activity composite
    vector activity_index(
        vector intensity,
        vector duration,
        vector frequency
    ) {
        int N = rows(intensity);
        
        // Log-transform to stabilize scale
        vector[N] log_intensity = log1p(intensity);
        vector[N] log_duration  = log1p(duration);
        vector[N] log_frequency = log1p(frequency);
        
        // Standardize components
        vector[N] z_intensity = zscore(log_intensity);
        vector[N] z_duration  = zscore(log_duration);
        vector[N] z_frequency = zscore(log_frequency);

        // Equal-weight additive MET proxy
        vector[N] activity;
        
        // Calculate the MEts 
        activity =
        (z_intensity + z_duration + z_frequency) / 3;
        
        return activity;
    }
}
// Data Input Block 
data{
    int<lower=1> N;
    vector[N] treatment;

    // =============== Composites Building Block ======================
    // Mass
    int<lower=1> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[N] mass_pre_sign;
    int<lower=1> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[N] mass_post_sign;

    // Diet 
    int<lower=1> J_diet_post; matrix[N, J_diet_post] diet_post_items; vector[N] diet_post_sign;
    vector[N] diet_pre;

    // Self Care 
    int<lower=1> J_self_care_post; matrix[N, J_self_care_post] self_care_post_items; vector[N] self_care_post_sign;
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
   vector[N] self_care_post = composite(self_care_post_items, self_care_post_sign);

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
    real<lower=0> sigma_diet;

    // Self Care Submodel Paramters 
    real alpha_self_care;
    real beta_treatment_self_care;
    real beta_self_care_pre;
    real<lower=0> sigma_self_care;

    // Activity Submodel Paramters 
    real alpha_activity;
    real beta_treatment_activity;
    real beta_activity_pre;
    real<lower=0> sigma_activity;

    // Mass Post final Equation Paramters
    real alpha_mass;
    real beta_treatment_mass;
    real beta_mass_pre;
    real beta_diet_mass;
    real beta_activity_mass;
    real beta_selfcare_mass;
    real<lower=0>sigma_mass;  
}
// Model 
model{
    // ========  Diet Submodel Priors and Likelihood ==============
    alpha_diet ~ normal(0, 1);
    beta_treatment_diet ~ normal(0 ,1);
    beta_diet_pre       ~ normal(0, 1);
    sigma_diet ~ exponential(1);

    diet_post_stand ~ normal(
        alpha_diet 
        + beta_treatment_diet * treatment
        + beta_diet_pre       * diet_pre_stand
        , sigma_diet
    );

    // ========== Self Care Submodel Priors and Likelihood ========================
    alpha_self_care ~ normal(0, 1);
    beta_treatment_self_care ~ normal(0, 1);
    beta_self_care_pre       ~ normal(0, 1);
    sigma_self_care ~ exponential(1);
    
    self_care_post_stand ~ normal(
        alpha_self_care 
        + beta_treatment_self_care * treatment
        + beta_self_care_pre       * self_care_pre_stand
        ,sigma_self_care
    );

    // ========== Activity Submodel Priors and Likelihood ==========================
    alpha_activity ~ normal(0, 1);
    beta_treatment_activity ~ normal(0 ,1);
    beta_activity_pre       ~ normal(0 ,1);
    sigma_activity ~ exponential(1);

    activity_post ~ normal(
        alpha_activity 
        + beta_treatment_activity * treatment
        + beta_activity_pre       * activity_pre
        ,sigma_activity
    );

    // ==============  Mass Model Final Equation Priors and Likelihood ================
    alpha_mass ~ normal(0, 1); 
    beta_treatment_mass ~ normal(0, 1);
    beta_mass_pre ~  normal(0, 1);
    beta_diet_mass ~ normal(0 ,1);
    beta_activity_mass ~ normal(0, 1);
    beta_selfcare_mass ~ normal(0, 1); 
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
   // Linear Predictors 
   vector[N] mu = alpha_mass 
        + beta_treatment_mass * treatment
        + beta_mass_pre       * mass_pre_stand
        + beta_diet_mass      * diet_post_stand
        + beta_activity_mass  * activity_post
        + beta_selfcare_mass  * self_care_post;

    // Posterior Predictive Draws & Pointwise log Likehood 
    vector[N] mass_post_stand_rep;
    vector[N] log_lik;
    for(i in 1:N){
        mass_post_stand_rep[i] = normal_rng(mu[i], sigma_mass);
        log_lik[i]             = normal_lpdf(mass_post_stand[i] | mu[i], sigma_mass);
    }

    // Bayesian R2
    real R2;
    {
        real fit_var = variance(mu);
        R2           = fit_var / (fit_var + square(sigma_mass)); 
    }
}
