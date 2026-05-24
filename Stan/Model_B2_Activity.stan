functions {
    // Zscore function 
    vector zscore(vector x) {
        return (x - mean(x)) / sd(x);
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
    real<lower=0> sigma;
}
// Model Block 
model{
    // Priors 
    alpha             ~ normal(0, 1);
    beta_treatment    ~ normal(0, 1);
    beta_activity_pre ~ normal(0, 1);
    sigma             ~ exponential(1); 

    // Model Likelihood 
    activity_post ~ normal(alpha + beta_treatment * treatment + beta_activity_pre * activity_pre, sigma);
}
// Aditional Calculations 
generated quantities {
   // Linear Predictions 
   vector[N] mu = alpha + beta_treatment * treatment + beta_activity_pre * activity_pre;

   // Posterior predictive draws & Pointwise log-likelihood (for LOO / WAIC)
   vector[N] activity_post_rep;
   vector[N] log_lik;
   for(i in 1:N){
    activity_post_rep[i] = normal_rng(mu[i], sigma);
    log_lik[i]           = normal_lpdf(activity_post[i] | mu[i], sigma);
   }

   // Bayesian R2
   real R2;
   {
    real var_fit = variance(mu);
    R2           = var_fit / (var_fit + square(sigma)); 
   }
   
   // Counterfactual predictions
   vector[N] mu_treated;
   vector[N] mu_control;
   for(i in 1:N){
    mu_treated[i] =
    alpha
    + beta_treatment
    + beta_activity_pre * activity_pre[i];
    
    mu_control[i] =
    alpha
    + beta_activity_pre * activity_pre[i]; 
    }

    // ITE ATE & Adjusted Group Means 
    vector[N] ITE = mu_treated - mu_control;
    real ATE      = mean(ITE);
    real mean_treated =  mean(mu_treated);
    real mean_control =  mean(mu_control);

}
