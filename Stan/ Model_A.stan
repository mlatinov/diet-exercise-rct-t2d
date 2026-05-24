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
}
// Input Data Block 
data{
    int<lower=0> N; // N Observations 
    vector[N] treatment;

    // Item Block per contruct
    int<lower=1> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[J_mass_post] mass_post_sign;
    int<lower=1> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[J_mass_pre] mass_pre_sign;
}
// Transformation Data Block 
transformed data {
    // Build the Indices
   vector[N] mass_post = composite(mass_post_items, mass_post_sign);
   vector[N] mass_pre  = composite(mass_pre_items, mass_pre_sign);

   // Standatize the Indices
   vector[N] mass_post_stand = zscore(mass_post);
   vector[N] mass_pre_stand  = zscore(mass_pre);

   // For later back transformation 
   real mass_post_composite_sd = sd(mass_post);
}
// Model Parameters Block
parameters{
    real beta_treatment;
    real beta_mass_pre;
    real alpha;
    real<lower=0> sigma;
}
// Model Block 
model{
    // Priors
    beta_treatment ~ normal(0,   0.5);
    beta_mass_pre  ~ normal(0.8, 0.3);
    alpha          ~ normal(0,   0.5);
    sigma          ~ exponential(1);

    // Model Likelihood 
    mass_post_stand ~ normal(
        alpha + beta_treatment * treatment + beta_mass_pre * mass_pre_stand, 
        sigma
    );
}
// Additional Calculations 
generated quantities {
    // Linear predictor  per observation 
    vector[N] mu = alpha + beta_treatment * treatment + beta_mass_pre * mass_pre_stand;

    // Posterior predictive draws
    vector[N] mass_post_rep;
    for (i in 1:N)
        mass_post_rep[i] = normal_rng(mu[i], sigma);

    // Pointwise log-likelihood (for LOO / WAIC) 
    vector[N] log_lik;
    for (i in 1:N)
        log_lik[i] = normal_lpdf(mass_post_stand[i] | mu[i], sigma);

    // Bayesian R^2 
    real R2;
    {
        real var_fit = variance(mu);
        R2 = var_fit / (var_fit + square(sigma));
    }

    //  Average Treatment Effect 
    real ATE_std       = beta_treatment;                          // standardized scale
    real ATE_composite = beta_treatment * mass_post_composite_sd; // composite's own scale
}