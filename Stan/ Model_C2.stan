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

    // Composit Building blocks ==============================================================================

    // Mass 
    int<lower=1> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[N] mass_pre_sign;
    int<lower=1> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[N] mass_post_sign;

    // Activity Post 
    vector[N] intensity_post;
    vector[N] duration_post ;
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
    real sigma;
}
// Model Block 
model{
    // Priors 
    alpha ~ normal(0, 1);
    beta_treatment     ~ normal(0, 1);
    beta_activity_post ~ normal(0, 1);
    beta_mass_pre      ~ normal(0, 1);
    sigma ~ exponential(1); 

    // Model Likelihood 
    mass_post_stand ~ normal(
        alpha 
        + beta_treatment     * treatment 
        + beta_activity_post * activity_post 
        + beta_mass_pre      * mass_pre_stand
        ,sigma
    );
}
// Additional Calculations 
generated quantities {
   // Linear predictor 
   vector[N] mu = 
        alpha      
        + beta_treatment     * treatment 
        + beta_activity_post * activity_post 
        + beta_mass_pre      * mass_pre_stand;

    // Posterior Predictive Reps & Pointwise log-likelihood (for LOO / WAIC) 
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

