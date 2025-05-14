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
  int<lower=0> nevent; // number of events
  vector[nobs] t_points;  // observed time points (non-strict decreasing)
  vector[nobs] event_flag;  // whether or not an event happened at the conresponding time point
  
  real last_event_time; // # time point where the last event occurs
  int event_indices[nevent]; # the row ID where events happened
  
  
  // predictor matrices (time-fixed)
  int<lower=0> p; // num main effect
  int<lower=0> q; // num interaction effect
  
  matrix[nobs,p] x; // for rows with both events and cencored
  matrix[nobs,q] x_int; 
  
  
  
  // link the interaction effect with the corresponding main effects
  int g1[q];
  int g2[q];
  
  
  // for prediction
  int<lower=0> nnew;
  vector[nnew] t_new;
  matrix[nnew,p] x_new; // for rows with events
  matrix[nnew,q] x_int_new;
  
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
    for (i in 1:p){
      Beta[i] ~ normal(0, 100);
      }
      
    for (i in 1:q){
      Beta_int[i] ~ normal(0, 100);
      }

      
    // partical log-likelihood, represented by [target]
    if (nobs > 0) {
      eta = x * Beta + x_int * Beta_int;
      }
    
    // Assumes: 
    // - event_indices is sorted
    // - eta is the linear predictor vector
    // - nevent is the number of events


    for (n in 1:nevent) {
      int current_idx = event_indices[n];
      if (n == 1){
        if (current_idx == 1){
          log_denom = eta[current_idx];  // First event
        } else {
          log_denom = log_sum_exp(eta[1:current_idx]); 
        }
      } else {
         int prev_idx = event_indices[n - 1];
         if (current_idx != prev_idx + 1) {
          int a = prev_idx + 1;
          int b = current_idx;
          log_denom = log_sum_exp(log_denom, log_sum_exp(eta[a:b]));
        } else {
            log_denom = log_sum_exp(log_denom, eta[current_idx]);
        }
      }
      // Add log-likelihood contribution
      target += eta[current_idx] - log_denom;
    }
}

generated quantities{
      // Predicting the survival time on the new/test dataset
      vector[nnew] eta_new;
      eta_new = x_new * Beta + x_int_new * Beta_int;
}

