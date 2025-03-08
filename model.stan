functions {
    real linear_link(real input) {
        real output;
        output = input * (1 / (1 + exp(-3 * input)));
        return output;
    }
}


data {
  int<lower=1> N;  // Total number of observations
  int<lower=1> J;  // Number of entities
  int<lower=1> temps;  // Number of time points (max observed)
  int<lower=1> id[N];   // Entity index for each observation
  int <lower=0, upper=1> statut[N];
  vector[N] y;    // Response variable
  vector[N] delaidiag;    // Predictor variable
}

parameters {
  vector[J] beta;                  // Fixed effects coefficients
  vector[J] alpha;                 // Random intercepts for each group
  vector[J] gamma; 
  vector[J] tau;                    // Random slopes for each group (time effect)
  real<lower=0> sigma;             // Standard deviation of residuals
  real<lower=0> sigma_alpha;         // SD of the group intercepts
  real<lower=0> sigma_gamma;         // SD of the group slopes
  real<lower=0> sigma_beta;
  real<lower=0> sigma_tau;
  real eta; 
}
model {
  // Priors
  beta ~ normal(0, sigma_beta);
  alpha ~ normal(0, sigma_alpha);
  gamma ~ normal(0, sigma_gamma);
  tau ~ normal(-10, sigma_tau);
  sigma ~ cauchy(0, 2);
  sigma_alpha ~ cauchy(0, 2);
  sigma_gamma ~ cauchy(0, 2);
  sigma_beta ~ cauchy(0,2);
  sigma_tau ~ cauchy(0,2);
  eta ~ normal(-2.6, 0.2);

  // Likelihood
  for (i in 1:N) {
    real spline = linear_link(delaidiag[i] -(tau[id[i]]));
    y[i] ~ normal(alpha[id[i]] + gamma[id[i]] * delaidiag[i] + beta[id[i]]*pow((exp(eta)/(1+exp(eta))), 1 - statut[i])*spline, sigma);
  }
}

