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
  
  int<lower=0> unique_nevent;
  int first_indices_events[unique_nevent];
  int last_indices_events[unique_nevent];

  
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
  }

model {
       // pre-allocated variables
    vector[nobs] eta; // for events & right censored
    real log_denom_lhs;
    real numerator;


    // prior
    for (i in 1:p){
      Beta[i] ~ normal(0, 100);
      }
      
    for (i in 1:q){
      Beta_int[i] ~ normal(0, 100);
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
      int start = first_indices_events[j];
      int end = last_indices_events[j];
      int len = end - start + 1;
      real log_len = log(len);
      if (len == 1){
        numerator = eta[start];
      } else {
        numerator = sum(eta[start:end]);
      }
      
      if (j == 1){ // initalized
        if (start != 1){
          log_denom_lhs = log_sum_exp(eta[1:end]);
        } else {
            log_denom_lhs = log_sum_exp(eta[start:end]);
        }
      } else {
        int last_end = last_indices_events[j-1];
        log_denom_lhs = log_sum_exp(log_denom_lhs, log_sum_exp(eta[last_end+1:end]));
      }
      
      vector[len] diff;
      for (ell in 1:len){
        if (ell == 1){
          diff[ell] = log_denom_lhs;// or some other appropriate value
        } else {
          diff[ell] = log_diff_exp(log_denom_lhs, log_sum_exp(eta[start:end]) - log_len + log(ell -1 ));
        }
      }
      
      target += numerator - sum(diff);
      }
}

generated quantities{
  // Predicting the survival time on the new/test dataset
  vector[nnew] eta_new;
  eta_new = x_new * Beta + x_int_new * Beta_int;
}
