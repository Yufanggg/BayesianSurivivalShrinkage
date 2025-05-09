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
  int<lower=0> nobs;
  int<lower=0> nevent;
  vector[nobs] t_point;  // time of events (non-strict decreasing)
  int<lower=1> event_indices[nevent];  // indices of events // time of right censoring (non-strict decreasing)
  
  
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

      
    // partical log-likelihood, represented by [target]
    if (nobs > 0) {
      eta = x * Beta + x_int * Beta_int;
      }
    
    // Ensure event_indices is declared as int[] in the data block
    int last_event_index = event_indices[num_elements(event_indices)];  // last event time point
    if (last_event_index == nobs) {  // if the last time point is the event time point
      log_denom = 0;
      } else {  // if the last time point is a censored time point
        int start_idx = last_event_index;
        int end_idx = num_elements(eta);
        vector[end_idx - start_idx + 1] eta_censored_tail = eta[start_idx:end_idx];
        log_denom = log_sum_exp(eta_censored_tail);  // contribution from censored cases
        }
        
    for (n in 1:nevent) {
      if (n != 1 && event_indices[n] != event_indices[n - 1] + 1) {
        int a = event_indices[n - 1] -1;
        int b = event_indices[n];
        log_denom = log_sum_exp(log_denom, log_sum_exp(eta[a:b]));
        } else {
          log_denom = log_sum_exp(log_denom, eta[event_indices[n]]);
          }
        target += eta[event_indices[n]] - log_denom;  // log-likelihood contribution
        }
          
}

generated quantities{
      // Predicting the survival time on the new/test dataset
      vector[nnew] eta_new;
      eta_new = x_new * Beta + x_int_new * Beta_int;
}

