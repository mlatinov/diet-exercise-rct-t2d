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
// Data input Block 
data{
    int<lower=1> N;
    vector[N] treatment;

    // Composite Building block 
    int<lower=1> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[J_mass_pre]  mass_pre_sign;
    int<lower=1> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[J_mass_post] mass_post_sign;
    int<lower=1> J_self_care_post; matrix[N, J_self_care_post] self_care_post_items; vector[J_self_care_post] self_care_post_sign;
}
// Transform Data Block 
transformed data {
   // Build the indices 
   vector[N] mass_pre  = composite(mass_pre_items, mass_pre_sign);
   vector[N] mass_post = composite(mass_post_items, mass_post_sign);
   vector[N] self_care_post = composite(self_care_post_items, self_care_post_sign);

   // Standatize all indices 
   vector[N] mass_pre_stand  = zscore(mass_pre);
   vector[N] mass_post_stand = zscore(mass_post);
   vector[N] self_care_post_stand = zscore(self_care_post); 
}
// Model Paramters block 
parameters{
    real alpha;
    real beta_treatment;
    real beta_mass_pre;
    real beta_self_care_post;
    real<lower=0> sigma;
}
// Model block
model{
    // Priors 
    alpha ~ normal(0, 1);
    beta_treatment      ~ normal(0, 1);
    beta_mass_pre       ~ normal(0, 1);
    beta_self_care_post ~ normal(0, 1);
    sigma ~ exponential(1);
    
    // Model Likelihood 
    mass_post_stand ~ normal(
        alpha 
        + beta_treatment      * treatment 
        + beta_mass_pre       * mass_pre_stand
        + beta_self_care_post * self_care_post_stand
        , sigma
    );
}
// Additional Calculatations 
generated quantities {
   // Linear Predictor 
   vector[N] mu = alpha 
        + beta_treatment      * treatment 
        + beta_mass_pre       * mass_pre_stand
        + beta_self_care_post * self_care_post_stand;

    // Posterior Predictive Draws & Pointwise log Likehood 
    vector[N] mass_post_stand_rep;
    vector[N] log_lik;
    for(i in 1:N){
        mass_post_stand_rep[i] = normal_rng(mu[i], sigma);
        log_lik[i]             = normal_lpdf(mass_post_stand[i] | mu[i], sigma); 
    }

    // Bayesian R2
    real R2;
    {
        real var_fit = variance(mu);
        R2           = var_fit / (var_fit + square(sigma));
    }
     
}