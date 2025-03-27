//
// This Stan program defines a survival model of cox regression.
// The baseline hazard in this stan program is assumed to be exponential distributed.
//

// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//


functions{
   /**
  * Log hazard for exponential distribution
  *
  * @param eta Vector, linear predictor
  * @return A vector
  */
  vector exponential_log_haz(vector eta, real lambda) {
    return eta + log(lambda);
  }
  
  /**
  * Log survival for exponential distribution
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @return A vector
  */
  vector exponential_log_surv(vector eta, vector t, real lambda) {
    vector[rows(eta)] res;
    res = (- t .* exp(eta)) * lambda;
    return res;
  }
  
  
}


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
  real<lower=0>  lambda;
  }

model {
    // pre-allocated variables
    vector[nevent] eta_event; // for events
    vector[nrcens] eta_rcens; // for right censored

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
    
    lambda ~ lognormal(0, 1);

    // log-likelihood, represented by [target]
    if (nevent > 0) {
      eta_event = x_event * Beta + x_int_event * Beta_int;
      target +=  exponential_log_haz (eta_event, lambda);
      target +=  exponential_log_surv(eta_event, t_event, lambda); // uncensored data log(f(t)) = log(h(t)) + log(S(t))
      }
      
    if (nrcens > 0) {
      eta_rcens = x_rcens * Beta + x_int_rcens * Beta_int;
      target +=  exponential_log_surv(eta_rcens, t_rcens, lambda); // right censored data log(S(t))
      }
}

generated quantities{
      // Predicting the survival time on the new/test dataset
      vector[nnew] survival_prob;  // 
      survival_prob = exp(exponential_log_surv(x_new * Beta + x_int_new * Beta_int, t_new, lambda));
}

