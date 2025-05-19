//
// This Stan program defines a survival model of cox regression.
// The baseline hazard in this stan program is assumed to be weibull distributed.
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
functions {
  int num_unique_starts(vector t) { //The first calculates how many different survival times occur in the data.
    if (size(t) == 0) return 0;
    int us = 1;
    for (n in 2:size(t)) {
      if (t[n] != t[n - 1]) us += 1;
    }
    return us;
  }
  
  array[] int unique_starts(vector t, int J) { // compute the value J to send into the function that computes the position in the array of failure times where each new failure time starts, plus an end point that goes one past the target
    array[J + 1] int starts;
    if (J == 0) return starts;
    starts[1] = 1;
    int pos = 2;
    for (n in 2:size(t)) {
      if (t[n] != t[n - 1]) {
    starts[pos] = n;
    pos += 1;
      }
    }
    starts[J + 1] = size(t) + 1;
    return starts;
  }
}


data {
  // response and time variables
  int<lower=0> nobs; // number of observations, including events, and censored
  vector[nobs] t_points;  // observed time points (non-strict decreasing)
  vector[nobs] event_flag;  // whether or not an event happened at the conresponding time point
  
  real last_event_time; // # time point where the last event occurs
  
  
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

transformed data {
  int<lower=0> J = num_unique_starts(t_points); // number of unique survival time
  array[J + 1] int<lower=0> starts = unique_starts(t_points, J);
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
    
    for (j in 1:J) {
      int start = starts[j];
      int end = starts[j + 1] - 1;
      int len = end - start + 1;
      real log_len = log(len);
      
      if (j == 1){ // first unqiue time point, initializ the log_denom
        if (len == 1){
          log_denom = eta[start]; // single event
        } else {
          log_denom = log_sum_exp(eta[start:end]); // multiple events at first time point
        }
      } else {
        log_denom = log_sum_exp(append_row(rep_vector(log_denom, 1), eta[start:end]));
        }
      // Add log-likelihood contribution
      if (event_flag[j] == 1){
        target += log_sum_exp(eta[start:end]) - log_denom;
      }
    } 
      
}

generated quantities{
      // Predicting the survival time on the new/test dataset
      vector[nnew] eta_new;
      eta_new = x_new * Beta + x_int_new * Beta_int;
}

