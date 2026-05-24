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
    int<lower=1> N;
    vector[N] treatment;

    // Composite Building Blocks =========================================================================

    // Mass Blocks 
    int<lower=0> J_mass_pre;  matrix[N, J_mass_pre]  mass_pre_items;  vector[J_mass_pre] mass_pre_sign;
    int<lower=0> J_mass_post; matrix[N, J_mass_post] mass_post_items; vector[J_mass_post] mass_post_sign;

    // Diet Block
    int<lower=0> J_diet_post; matrix[N, J_diet_post] diet_post_items; vector[J_diet_post] diet_post_sign;
}
// Data Transformation 
transformed data {
   // Build the Composites 
   vector[N] mass_pre  = composite(mass_pre_items, mass_pre_sign);
   vector[N] mass_post = composite(mass_post_items, mass_post_sign);
   vector[N] diet_post = composite(diet_post_items, diet_post_sign);

   // Standartize all Composits
   vector[N] mass_pre_stand   = zscore(mass_pre);
   vector[N] mass_post_stand  = zscore(mass_post);
   vector[N] diet_post_stand  = zscore(diet_post);

    // For later back transformation 
   real mass_post_composite_sd = sd(mass_post);
}
// Model Paramters 
parameters{
    real alpha;
    real beta_treatment;
    real beta_mass_pre;
    real beta_diet_post;
    real<lower=0> sigma;
}
// Model Block
model{
    // Priors
    alpha ~ normal(0, 1);
    beta_treatment ~ normal(0, 1);
    beta_mass_pre  ~ normal(0, 1);
    beta_diet_post ~ normal(0, 1);
    sigma ~ exponential(1);

    // Model Likelihood 
    mass_post_stand ~ normal(
        alpha 
        + beta_treatment * treatment 
        + beta_mass_pre  * mass_pre_stand 
        + beta_diet_post * diet_post_stand
        ,sigma
    );
} 
// Aditional Calculations 
generated quantities {
   // Linear Predictor
   vector[N] mu =  alpha 
        + beta_treatment * treatment 
        + beta_mass_pre  * mass_pre_stand 
        + beta_diet_post * diet_post_stand;

    // Posterior Predictions & Pointwise log-likelihood (for LOO / WAIC)
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
    
    //  Average Treatment Effect 
    real ATE_std       = beta_treatment;                          // standardized scale
    real ATE_composite = beta_treatment * mass_post_composite_sd; // composite's own scale
}