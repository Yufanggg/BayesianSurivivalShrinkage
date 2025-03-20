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
  int<lower=0> p; // num main effect
  int<lower=0> q; // num interaction effect
  int <lower=0> N; // num uncensored obs
  matrix[N, p] x; // covariates for uncensored obs
  matrix[N, q] x_int;
  
  int g1[q];
  int g2[q];
  
  int <lower=0> N_cens; // num censored obs
  matrix[N_cens, p] x_cens; // covariates for censored obs
  matrix[N_cens, q] x_int_cens;
}


parameters {
  vector[p] Beta; // coefficients for design matrix; without intercept
  vector[q] Beta_int;
  
  real<lower=0.01, upper=1> tau2int;
  real<lower=0>  tau2[p];
  real<lower=0>  gam2[p];
}

model {
  //Prior
  coefficients ~ normal(0, 1);
  tau2int ~ uniform(0.01,1);
  for (k in 1:p){
    gam2[k] ~ inv_gamma(0.5, 1);
    tau2[k] ~ inv_gamma(0.5, 1/gam2[k]);
    }
    
  for (i in 1:p){
    Beta[i] ~ normal(0, sqrt(tau2[i]));
    }
  for (i in 1:q){
    Beta_int[i] ~ normal(0, sqrt(sqrt(tau2[g1[i]]*tau2[g2[i]])*tau2int));
    }
    
  vector[N] log_theta = x * Beta + x_int * Beta_int;
  vector [N_cens] log_theta_c = x_cens * Beta + + x_int_cens * Beta_int;
  real log_denom = log_sum_exp(log_theta_c); //log_sum_exp is defined as the logarithm of the sum of exponentials of the input values
  target += log_theta - log_denom;
}

