//
// This Stan program defines a survival model of cox regression.
// The baseline hazard in this stan program is assumed to be weibull distributed.
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  // response and time variables
  int<lower=0> nobs; // number of observations, including events, and censored
  vector[nobs] t_points;  // observed time points (non-strict decreasing)
  vector[nobs] event_flag;  // whether or not an event happened at the conresponding time point
  
  int<lower=0> unique_nevent; // # number of events
  vector[unique_nevent] unique_event_times;
  
  // predictor matrices (time-fixed)
  int<lower=0> p; // num main effect
  int<lower=0> q; // num interaction effect
  
  matrix[nobs,p] x; // for rows with both events and cencored
  matrix[nobs,q] x_int; 
  
  
  
  // link the interaction effect with the corresponding main effects
  int g1[q];
  int g2[q];
  
}

parameters {
  vector[p] Beta; // coefficients for design matrix;
  vector[q] Beta_int;
  
  real<lower=0.01, upper=1> tau2int;
  real<lower=0>  tau2[p];
  real<lower=0>  gam2[p];
  }

model {
       // pre-allocated variables
    vector[nobs] eta; // for events & right censored
    
    real log_denom;


    // prior
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
      
      // Partial log-likelihood contribution
      if (nobs > 0) {
        eta = x * Beta + x_int * Beta_int;
        }

    
    // Assumes: 
    // - event_indices is sorted
    // - eta is the linear predictor vector
    // - nevent is the number of events
    // - t_points is also sorted
    
    for (j in 1:unique_nevent) {
      unique_event_time = unique_event_times[j]
      real log_num; real log_denom;
      for (t in 1:nobs){
        if(t_points[t] == unique_event_time && event_flag[i] == 1){
          log_num = log_sum_exp(log_num, eta[i]);
        }
        if (t_points[i] >= unique_event_time) {
          log_denom = log_sum_exp(log_denom, eta[i]);
      }
    }
    target += log_num - log_denom
      
}


