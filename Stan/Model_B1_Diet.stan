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
// Input Data block 
data{

    int<lower=0> N; // Number of observations 
    vector[N] treatment;
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
    real<lower=0> sigma;
}
// Model Block 
model{
    // Priors 
    alpha          ~ normal(0, 1);
    beta_treatment ~ normal(0, 1);
    beta_diet_pre  ~ normal(0, 1);
    sigma          ~ exponential(1);
    
    // Model Likelihood
    diet_post_stand ~ normal(
        alpha + beta_treatment * treatment + beta_diet_pre * diet_pre_stand,
        sigma
    );
}
// Aditional calculations 
generated quantities {
   // Linear predictors
   vector[N] mu =  alpha + beta_treatment * treatment + beta_diet_pre * diet_pre_stand;

   // Posterior predictive draws 
   vector[N] diet_post_rep;
   for(i in 1:N){
    diet_post_rep[i] = normal_rng(mu[i], sigma);
   }
   
   // Pointwise log-likelihood (for LOO / WAIC) 
   vector[N] log_link;
   for(i in 1:N){
    log_link[i] = normal_lpdf(diet_post_stand[i] | mu[i], sigma);
   }

   // Bayesian R^2 
   real R2;
   {
    real var_fit = variance(mu);
    R2           = var_fit / (var_fit + square(sigma)); 
   }

   // Avetage Treatment Effect ATE 
   real ATE_std       = beta_treatment;
   real ATE_composite = beta_treatment * diet_post_sd;

}

