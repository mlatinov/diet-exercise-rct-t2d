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
    int<lower=2> J_self_care_post; matrix[N, J_self_care_post] self_care_post_items; vector[J_self_care_post] self_care_post_sign;
    
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
   vector[N] self_care_post = composite(self_care_post_items, self_care_post_sign);
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
    real<lower=0> sigma;
}
// Model Block 
model{
    alpha ~ normal(0, 1);
    beta_treatment ~ normal(0, 1);
    beta_mass_pre  ~ normal(0, 1);
    beta_diet_post ~ normal(0, 1);
    beta_self_care_post ~ normal(0, 1);
    beta_activity_post  ~ normal(0, 1);
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
   // Linear Predictor 
   vector[N] mu = 
     alpha 
        + beta_treatment      * treatment
        + beta_self_care_post * self_care_post_stand
        + beta_diet_post      * diet_post_stand
        + beta_self_care_post * self_care_post_stand
        + beta_activity_post  * activity_post;
    
    // Posterior Predictive Draws & Poitwise log Likeluhood 
    vector[N] mass_post_stand_rep;
    vector[N] log_lik;
    for(i in 1:N){
        mass_post_stand_rep[i] = normal_rng(mu[i], sigma);
        log_lik[i]             = normal_lpdf(mass_post_stand[i] | mu[i], sigma);
    }

    // Bayesian R2
    real R2;
    {
        real fit_var = variance(mu);
        R2           = fit_var / (fit_var + square(sigma)); 
    }
}