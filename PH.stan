//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int <lower=0> K; //num covariates
  int <lower=0> N; // num uncensored obs
  
  matrix[N, K] x; // covariates for uncensored obs
  int N_cens; // num censored obs
  matrix[N_cens, K] x_cens; // covariates for censored obs
}


parameters {
  vector[K] beta; // slopes without intercept
}

model {
  beta ~ normal(0, 2); // prior
  vector[N] log_theta = x * beta;
  vector [N_cens] log_theta_c = x_cens * beta;
  real log_denom = log_sum_exp(log_theta_c); //log_sum_exp is defined as the logarithm of the sum of exponentials of the input values
  target += log_theta - log_denom;
}

