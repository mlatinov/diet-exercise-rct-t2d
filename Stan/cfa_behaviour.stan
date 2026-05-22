// Data Block
data{
  int<lower=1> N;
  int<lower=1> J_cont;
  matrix[N, J_cont] Y_cont;
  int<lower=1> J_bin;
  array[N, J_bin] int<lower=0, upper=1> Y_bin;
}

// Transformation Block
transformed data {
  // Standartize all numerical indicators 
   vector[J_cont] Y_mean;
   vector[J_cont] Y_sd;
   matrix[N, J_cont] Z;
   for (j in 1:J_cont){
    Y_mean[j] = mean(Y_cont[, j]);
    Y_sd[j]   = sd(Y_cont[, j]);
    Z[, j]    = (Y_cont[, j] - Y_mean[j]) / Y_sd[j];
   }
}

// Parameters Block
parameters{
  // Non-centred latent variable
  vector[N] eta_raw;

  // Continuous Parameters
  vector[J_cont] nu_cont;
  vector<lower=0.05>[J_cont] sigma_cont;
  vector[J_cont] lambda_cont;         

  // Binary Parameters
  vector[J_bin] nu_bin;
  vector[J_bin] lambda_bin;
}

transformed parameters {
  // Non-centred: eta has unit variance, but parameterised through eta_raw
  vector[N] eta = eta_raw;
}

// Model
model{
  // Latent variable prior — unit variance
  eta_raw ~ std_normal();

  // Continuous Priors
  nu_cont      ~ normal(0, 0.5);
  lambda_cont  ~ normal(0, 1);          
  sigma_cont   ~ lognormal(0, 0.5);

  // Binary Priors
  nu_bin     ~ normal(0, 0.5);
  lambda_bin ~ normal(0, 1);

  // Continuous Normal Likelihood
  for (j in 1:J_cont)
    Z[, j] ~ normal(nu_cont[j] + lambda_cont[j] * eta, sigma_cont[j]);

  // Binary Bernoulli probit Likelihood
  for (j in 1:J_bin)
    Y_bin[, j] ~ bernoulli(Phi(nu_bin[j] + lambda_bin[j] * eta));
}

generated quantities {
  // Anchor orientation on the first continuous loading if it's negative, flip everything
  real sign_flip = lambda_cont[1] < 0 ? -1.0 : 1.0;
  vector[J_cont] lambda_cont_oriented = sign_flip * lambda_cont;
  vector[J_bin]  lambda_bin_oriented  = sign_flip * lambda_bin;
  vector[N]      eta_oriented         = sign_flip * eta;

  // R2 per indicator (sign-invariant, but use oriented for consistency)
  vector[J_cont] R2_cont;
  vector[J_bin]  R2_bin;
  for (j in 1:J_cont)
    R2_cont[j] = square(lambda_cont[j]) / (square(lambda_cont[j]) + square(sigma_cont[j]));
  for (j in 1:J_bin)
    R2_bin[j] = square(lambda_bin[j]) / (square(lambda_bin[j]) + 1.0);

  // Back-transform continuous loadings/intercepts/sigmas to original Y scale (oriented)
  vector[J_cont] lambda_raw = lambda_cont_oriented .* Y_sd;
  vector[J_cont] sigma_raw  = sigma_cont .* Y_sd;
  vector[J_cont] nu_raw     = Y_mean + nu_cont .* Y_sd;

  // Posterior predictive draws — continuous, original scale
  matrix[N, J_cont] Y_cont_rep;
  for (j in 1:J_cont)
    for (i in 1:N)
      Y_cont_rep[i, j] = Y_mean[j]
        + Y_sd[j] * normal_rng(nu_cont[j] + lambda_cont[j] * eta[i], sigma_cont[j]);

  // Posterior predictive draws — binary
  array[N, J_bin] int Y_bin_rep;
  for (j in 1:J_bin)
    for (i in 1:N)
      Y_bin_rep[i, j] = bernoulli_rng(Phi(nu_bin[j] + lambda_bin[j] * eta[i]));

  // Pointwise log-likelihood
  vector[N] log_lik;
  for (i in 1:N) {
    real ll = 0;
    for (j in 1:J_cont)
      ll += normal_lpdf(Z[i, j] | nu_cont[j] + lambda_cont[j] * eta[i], sigma_cont[j]);
    for (j in 1:J_bin)
      ll += bernoulli_lpmf(Y_bin[i, j] | Phi(nu_bin[j] + lambda_bin[j] * eta[i]));
    log_lik[i] = ll;
  }
}