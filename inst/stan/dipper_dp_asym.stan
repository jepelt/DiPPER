data {
  int<lower=0> N;                                // Number of samples
  int<lower=0> K;                                // Number of taxa
  int<lower=0> P;                                // Number of predictors
  array[K, N] int<lower=0, upper=1> y;           // Presence/absence matrix
  matrix[N, P] X;                                // Design matrix

  // User-defined hyperparameters for priors
  real prior_alpha_mean;                         // Prior mean for alpha
  real<lower=0> prior_alpha_sd;                  // Prior SD for alpha
  real<lower=0> prior_nu_sd;                     // Standard deviation for nu
  real<lower=0> prior_tau_sd;                    // Standard deviation for tau
  vector[P - 1] prior_cov_mean;                  // Prior means for covariates
  vector<lower=0>[P - 1] prior_cov_sd;           // Prior SDs for covariates
}

parameters {
  vector[K] alpha;                               // Taxon-specific intercepts
  matrix[P - 1, K] beta_cov;                     // Coefficients for covariates

  // Parameters for the Asymmetric Laplace prior
  vector[K] z_norm;
  vector<lower=0>[K] z_exp;
  real<lower=0, upper=1> nu;
  real<lower=0> tau;
}

transformed parameters {
  vector[K] z;
  vector[K] beta;                                // Differential prevalence estimates
  matrix[P, K] B;

  real theta = (1.0 - 2.0 * nu) / (nu * (1.0 - nu));
  real tau_sq = 2.0 / (nu * (1.0 - nu));

  for (k in 1:K) {
    z[k] = theta * z_exp[k] + sqrt(tau_sq * z_exp[k] + 1e-8) * z_norm[k];
  }

  beta = z * tau;

  B[1, ] = beta';
  if (P > 1) {
    B[2:P, ] = beta_cov;
  }
}

model {
  // --- Priors ---
  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);

  // Covariate priors applied dynamically based on user/R-wrapper input
  if (P > 1) {
    for (p in 1:(P - 1)) {
      beta_cov[p] ~ normal(prior_cov_mean[p], prior_cov_sd[p]);
    }
  }

  // Hyperparameters
  nu ~ double_exponential(0.5, prior_nu_sd);
  tau ~ normal(0, prior_tau_sd);

  // Latent variables for the Asymmetric Laplace
  z_norm ~ std_normal();
  z_exp ~ exponential(1.0);

  // --- Likelihood ---
  for (k in 1:K) {
    y[k] ~ bernoulli_logit_glm(X, alpha[k], B[, k]);
  }
}
