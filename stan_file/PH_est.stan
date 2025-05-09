//
// This Stan program defines a survival model of cox regression.
// The baseline hazard in this stan program is assumed to be weibull distributed.
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
  int find_largest_index(vector t_rcens, real t_event_point) {
    int max_index = -1;// Initialize to -1 or another indicator for no valid index
    for (j in 1:num_elements(t_rcens)) {
      if (t_rcens[j] > t_event_point) {
        max_index = j; // Update max_index to the current index
        }
        }
        return max_index;
        }
        
}

data {
  // response and time variables
  int<lower=0> nevent;
  int<lower=0> nrcens;
  vector[nevent] t_event;  // time of events (non-strict decreasing)
  vector[nrcens] t_rcens;  // time of right censoring (non-strict decreasing)
  
  
  // predictor matrices (time-fixed)
  int<lower=0> p; // num main effect
  int<lower=0> q; // num interaction effect
  
  matrix[nevent,p] x_event; // for rows with events
  matrix[nevent,q] x_int_event; 
  
  matrix[nrcens,p] x_rcens; // for rows with right censoring
  matrix[nrcens,q] x_int_rcens; // for rows with right censoring
  
  
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
    vector[nevent] eta_event; // for events
    vector[nrcens] eta_rcens; // for right censored
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
    if (nevent > 0) {
      eta_event = x_event * Beta + x_int_event * Beta_int;
      }
      
    // real log_denom = log_sum_exp(eta_rcens);
    if (nrcens > 0){
       eta_rcens = x_rcens * Beta + x_int_rcens * Beta_int;
    }
    
    for (n in 1:nevent) {
      if  (nrcens > 0) {
        // Find indices in t_rcens that are larger than t_event[i]
        int largest_indices = find_largest_index(t_rcens, t_event[n]);
        if(largest_indices != -1){
          log_denom = log_sum_exp(eta_rcens[1:largest_indices]);
        } else {
          log_denom = log_sum_exp(0, 0);
        }
      }
      // log_denom = log_sum_exp(eta_rcens[1:])
      log_denom = log_sum_exp(log_denom, eta_event[n]);
      target += eta_event[n] - log_denom;   // log likelihood
  }
      
    
}
