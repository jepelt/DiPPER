data {
  int<lower=0> N;                 // Number of samples
  int<lower=0> K;                 // Number of taxa/features
  int<lower=0> P;                 // Number of predictors (excl. intercept)
  array[K, N] int<lower=0, upper=1> y; // Presence/absence matrix
  matrix[N, P] X;                 // Design matrix

  // User-defined hyperparameters for priors
  real prior_alpha_mean;          // Prior mean for fixed intercepts
  real<lower=0> prior_alpha_sd;   // Prior SD for fixed intercepts
  real<lower=0> prior_tau_sd;     // Prior SD for global scale tau
  vector[P - 1] prior_cov_mean;   // Prior means for covariates
  vector<lower=0>[P - 1] prior_cov_sd; // Prior SDs for covariates
}

parameters {
  vector[K] alpha;                // Taxon-specific fixed intercepts
  matrix[P - 1, K] beta_cov;      // Coefficients for covariates

  // Parameters for the Symmetric Laplace prior
  vector[K] z_norm;               // Standard normal auxiliary variables
  vector<lower=0>[K] z_exp;       // Standard exponential auxiliary variables
  real<lower=0> tau;              // Global scale parameter
}

transformed parameters {
  vector[K] z;                    // Unscaled symmetric Laplace variables
  vector[K] beta;                 // Scaled differential prevalence parameters
  matrix[P, K] B;                 // Full coefficient matrix

  for (k in 1:K) {
    // Construct Symmetric Laplace variables using stochastic representation
    // (nu = 0.5 -> theta = 0, tau_sq = 8.0)
    z[k] = sqrt(8.0 * z_exp[k] + 1e-8) * z_norm[k];
  }

  // Final differential prevalence parameters
  beta = z * tau;

  // Combine into the full coefficient matrix B (first row is beta)
  B[1, ] = beta';
  if (P > 1) {
    B[2:P, ] = beta_cov;
  }
}

model {
  // --- Priors ---
  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);

  if (P > 1) {
    for (p in 1:(P - 1)) {
      beta_cov[p] ~ normal(prior_cov_mean[p], prior_cov_sd[p]);
    }
  }

  // Laplace Hyperparameters
  tau ~ normal(0, prior_tau_sd);

  // Latent variables for the Symmetric Laplace prior
  z_norm ~ std_normal();
  z_exp ~ exponential(1.0);

  // --- Likelihood ---
  for (k in 1:K) {
    y[k] ~ bernoulli_logit_glm(X, alpha[k], B[, k]);
  }
}
