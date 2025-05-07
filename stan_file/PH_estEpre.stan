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
  int<lower=0> nevent;
  int<lower=0> nrcens;
  vector[nevent] t_event;  // time of events
  vector[nrcens] t_rcens;  // time of right censoring
  
  
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
  real<lower=0>  shape;
  real<lower=0>  lambda;
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
      
    shape ~ lognormal(0, 100);
    lambda ~ lognormal(0, 100);
      

    // partical log-likelihood, represented by [target]
    if (nrcens > 0) {
      eta_rcens = x_rcens * Beta + x_int_rcens * Beta_int;
      log_denom = log_sum_exp(eta_rcens); //log_sum_exp is defined as the logarithm of the sum of exponentials of the input values
      }
      
    if (nevent > 0) {
      eta_event = x_event * Beta + x_int_event * Beta_int;
      
      for (n in 1:nevent){
        log_denom = log_sum_exp(log_denom, eta_event[n]);
        target += eta_event[n] - log_denom; // log likelihood
      }
      }
      
    
}

generated quantities{
      // Predicting the survival time on the new/test dataset
      vector[nnew] risk_ratio_pred;  // 
      vector[nnew] eta_new;
      eta_new = x_new * Beta + x_int_new * Beta_int;
        
      risk_ratio_pred = exp(eta_new);
}

