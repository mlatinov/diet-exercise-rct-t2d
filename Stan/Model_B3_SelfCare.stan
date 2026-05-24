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

// Input data Block 
data{
    int<lower=1> N;
    vector[N] treatment;
    vector[N] self_care_pre;

    // Indices Blocks 
    int<lower=1> J_self_care_post; matrix[N, J_self_care_post] self_cate_post_items; vector[J_self_care_post] self_cate_post_sing;
}

// Data Transformation block 
transformed data {
   // Build the Indice 
   vector[N] self_care = composite(self_cate_post_items, self_cate_post_sing);

   // Standartize the Self Care Composite and self care pre
   vector[N] self_care_stand     = zscore(self_care);
   vector[N] self_care_pre_stand = zscore(self_care_pre);
   
   // For later back transformation 
   real self_care_post_composite_sd = sd(self_care);  
}

// Model Paramters 
parameters{
    real alpha;
    real beta_treatment;
    real beta_self_care_pre;
    real<lower=0> sigma;
}

// Model Block 
model{
    // Priors
    alpha              ~ normal(0, 1);
    beta_treatment     ~ normal(0, 1);
    beta_self_care_pre ~ normal(0, 1);
    sigma ~ exponential(1);

    // Model Likelihood
    self_care_stand ~ normal(alpha + beta_treatment * treatment + beta_self_care_pre * self_care_pre_stand, sigma);
}

// Additional Calculations 
generated quantities {
    // Linear predictor
    vector[N] mu = alpha + beta_treatment * treatment + beta_self_care_pre * self_care_pre_stand;

    // Posterior predictive draws & Pointwise log-likelihood (for LOO / WAIC)
    vector[N] self_care_stand_rep;
    vector[N] log_lik;
    for(i in 1:N){
        self_care_stand_rep[i] = normal_rng(mu[i], sigma);
        log_lik[i]             = normal_lpdf(self_care_stand[i] | mu[i], sigma); 
    }

    // Bayesian R2
    real R2;
    {
        real var_fit = variance(mu);
        R2           = var_fit / (var_fit + square(sigma)); 
    }

    //  Average Treatment Effect 
    real ATE_std       = beta_treatment;                               // standardized scale
    real ATE_composite = beta_treatment * self_care_post_composite_sd; // composite's own scale
}