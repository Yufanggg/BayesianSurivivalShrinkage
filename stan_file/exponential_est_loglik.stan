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
  
  // for log_lik
  int<lower=0> N; //number of Obs = nevent + nrcens
  vector[N] t; // time of Obs
  vector[N] status; // status of Obs
  matrix[N, p] X; // for rows with Obs
  matrix[N, q] X_int; // for rows with Obs
  
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
    
    lambda ~ lognormal(0, 100);

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

// save the log_lik explicitly
generated quantities{
  vector [N] log_lik;
  for (n in 1:N){
    if (status[n] == 1){
      log_lik[n] = exponential_lpdf(t[n] | exp(X[n] * Beta + X_int[n] * Beta_int) * lambda); // uncensored data log(f(t)) = log(h(t)) + log(S(t))
    } else {
      log_lik[n] = exponential_lccdf(t[n] | exp(X[n] * Beta + X_int[n] * Beta_int) * lambda); // right censored data log(S(t))
    }
  }
}



