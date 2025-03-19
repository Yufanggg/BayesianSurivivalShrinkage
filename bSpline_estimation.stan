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


functions{
  vector log_hi(matrix X, vector Beta, matrix  X_int, vector Beta_int, matrix bSpline_basis, vector coefficients){
    vector[num_elements(bSpline_basis)] log_hi = bSpline_basis * coefficients + X * Beta + X_int * Beta;
    
    return log_hi;
  }
  
  vector approxfun(vector x, vector x_pred, vector y_pred) {
    // x_pred, y_pred are the known points 
    int N = num_elements(x);
    int K = num_elements(x_pred);
    vector[N] y;
    
    for (n in 1:N) {
      int i = 1;
      
      // Check if x[n] exactly matches any x_pred value
      for (k in 1:K) {
        if (x[n] == x_pred[k]) {
          y[n] = y_pred[k];
          continue;
        }
      }
      
      // Find the interval for x[n]
      while (i < K && x_pred[i + 1] < x[n]) {
        i += 1;
      }
      
      // Linear interpolation
      if (i < K) {
        y[n] = y_pred[i] + (y_pred[i + 1] - y_pred[i]) * (x[n] - x_pred[i]) / (x_pred[i + 1] - x_pred[i]);
      } else {
        y[n] = y_pred[K]; // Extrapolate using the last interval
      }
    }
    
    return y;
  }

  
  real H_i(real Ti, vector T_known, vector locates, vector weights, matrix X, vector Beta, matrix X_int, vector Beta_int, matrix bSpline_basis, vector coefficients ){
    
    vector[num_elements(locates)] hi_seq = exp(log_hi(X, Beta, X_int, Beta_int, bSpline_basis, coefficients));
    vector[num_elements(locates)] hi_values = approxfun(Ti * (1+locates)/2, T_known, hi_seq);
    real integral_ = Ti * sum(weights .* hi_values);
    return integral_/2;
  }
    
  }

data {
  int<lower=0> p; // num main effect
  int<lower=0> q; // num interaction effect
  int <lower=0> N; // num uncensored obs
  int <lower=0> M; // num unique time points
  

  
  int g1[q];
  int g2[q];
  
  vector[N] t; // event time (non-strict decreasing)
  matrix[N, p] x; // covariates for uncensored obs
  matrix[N, q] x_int;
  
  int <lower=0> N_cens; // num censored obs
  vector[N_cens] t_cens; // censoring time
  matrix[N_cens, p] x_cens; // covariates for censored obs
  matrix[N_cens, q] x_int_cens;
  
  matrix[M, 6] bSpline_basis;
  vector[M] uniqueT; //montously increasing unique time points
  
  int<lower=0> N_new;
  matrix[N_new, p] x_new;
  matrix[N_new, q] x_new_int;
  
  vector[15] locates;
  vector[15] weights;
  
}

// The parameters accepted by the model. Our model
parameters {
  vector[6] coefficients;
  vector[p] Beta; // coefficients for design matrix;
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
  
  //log-likelihood, represented by [target]
  //log(f(t)) = -H(t) + log(h(t))for uncensored data
  for (n in 1:N){
    target += (-H_i(t[n], uniqueT, locates, weights, x, Beta, x_int, Beta_int, bSpline_basis, coefficients) + log_hi(x, Beta, x_int, Beta_int, bSpline_basis, coefficients));
  }
   
  
  // log(S(t)) = - H(t)for censored data
  for (n in 1:N_cens){
    target += (-H_i(t_cens[n], uniqueT, locates, weights, x_cens, Beta, x_int_cens, Beta_int, bSpline_basis, coefficients));
  }
  
}

generated quantities{// predicting the survival time on the new/test dataset by Monte Carlo sampling (needs to be verfied)
  matrix[N_new, M] survival_prob; // Nobs * unqiue time points 
  for (m in 1:M){
    real Hi = H_i(uniqueT[m], uniqueT, locates, weights, x_new, Beta, x_new_int, Beta_int, bSpline_basis, coefficients);
    survival_prob[, m] = exp(-Hi);
  }
}

