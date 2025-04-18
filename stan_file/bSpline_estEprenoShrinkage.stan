//
// This Stan program defines a survival model of cox regression.
// The baseline hazard in this stan program is approximated by bSplines.
// The integral for the bSpline baseline hazard is obatined by the Gauss-Kronrod quadrature.
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//


functions{
  /**
  * Log hazard for B-spline model
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param coefs Vector, B-spline coefficients
  * @return A vector
  */
  vector bspline_log_haz(vector eta, matrix bSpline_basis, vector coefs) {
    return bSpline_basis * coefs + eta;
  }
  
   /**
  * Evaluate log survival from the log hazard evaluated at
  * quadrature points and a corresponding vector of quadrature weights
  *
  * @param qwts Vector, the quadrature weights
  * @param log_hazard Vector, log hazard at the quadrature points
  * @param qnodes Integer, the number of quadrature points for each individual
  * @param N Integer, the number of individuals (ie. rows(log_hazard) / qnodes)
  * @return A vector
  */
  real quadrature_log_surv(vector qwts, vector log_hazard) {
    return - dot_product(qwts, exp(log_hazard)); // sum across all individuals
  }
  }


data {
  
  //-----for model fitting-----/
  // dimensions
  int<lower=0> p;          // num. cols in main effects (time-fixed)
  int<lower=0> q;          // num. cols in interaction effects (time-fixed)
  
  
/*  int<lower=0> nevent;     // num. rows w/ an event (ie. not censored)
  int<lower=0> nrcens;     // num. rows w/ right censoring*/
  
  int<lower=0> Nevent;     // num. rows w/ an event;      used only w/ quadrature
  int<lower=0> Nrcens;     // num. rows w/ right cens;    used only w/ quadrature

  int<lower=0> qnodes;     // num. nodes for GK quadrature
  int<lower=0> qevent;     // num. quadrature points for rows w/ an event
  int<lower=0> qrcens;     // num. quadrature points for rows w/ right censoring

  int<lower=0> nvars;      // num. aux parameters for baseline hazard



  // response and time variables
/*  vector[nevent] t_event;  // time of events
  vector[nrcens] t_rcens;  // time of right censoring*/


  vector[Nevent] epts_event;  // time of events
  vector[qevent] qpts_event;  // qpts for time of events
  vector[qrcens] qpts_rcens;  // qpts for time of right censoring



  // predictor matrices (time-fixed), with quadrature
  matrix[Nevent,p] x_epts_event; // for rows with events
  matrix[Nevent,q] x_int_epts_event; // for rows with events
  
  matrix[qevent,p] x_qpts_event; // for rows with events
  matrix[qevent,q] x_int_qpts_event; // for rows with events
  
  matrix[qrcens,p] x_qpts_rcens; // for rows with right censoring
  matrix[qrcens,q] x_int_qpts_rcens; // for rows with right censoring



  // basis matrices for B-splines, with quadrature
  matrix[Nevent,nvars] basis_epts_event; // at event time
  matrix[qevent,nvars] basis_qpts_event; // at qpts for event time
  matrix[qrcens,nvars] basis_qpts_rcens; // at qpts for right censoring time



  // GK quadrature weights, with (b-a)/2 scaling already incorporated
  vector[qevent] qwts_event;
  vector[qrcens] qwts_rcens;
  
  //-----for model prediction-----/
  int<lower=0> nnew;
  int<lower=0> qevent_new;
  vector[nnew] t_new;
  matrix[qevent_new, p] x_new_qpts_event;
  matrix[qevent_new, q] x_new_int_qpts_event;
  matrix[qevent_new,nvars] basis_qpts_event_new;
  vector[qevent_new] eta_qpts_event_new;

}

// The parameters accepted by the model. Our model
parameters {
  vector[nvars] coefs;
  vector[p] Beta; // coefficients for design matrix;
  vector[q] Beta_int;
  
  real<lower=0.01, upper=1> tau2int;
  real<lower=0>  tau2[p];
  real<lower=0>  gam2[p];
  
}


model {
  //pre-allocate variables
  vector[Nevent] eta_epts_event; // for event times
  vector[qevent] eta_qpts_event; // for qpts for event time
  vector[qrcens] eta_qpts_rcens; // for qpts for right censoring time

  
  //log-scale hazard models
  vector[Nevent] lhaz_epts_event;
  vector[qevent] lhaz_qpts_event;
  vector[qrcens] lhaz_qpts_rcens;
  
  
  //Prior
  coefs ~ normal(0, 100);
  tau2int ~ uniform(0.01,1);
  for (k in 1:p){
    gam2[k] ~ inv_gamma(0.5, 1);
    tau2[k] ~ inv_gamma(0.5, 1/gam2[k]);
    }
    
  for (i in 1:p){
    Beta[i] ~ normal(0, sqrt(tau2[i]));
    }
  for (i in 1:q){
    Beta_int[i] ~ normal(0, sqrt(tau2int));
    }
  
  
  //log-likelihood, represented by [target]
  if (Nevent > 0) {
    eta_epts_event = x_epts_event * Beta + x_int_epts_event * Beta_int;
    lhaz_epts_event = bspline_log_haz(eta_epts_event, basis_epts_event, coefs);// B-splines, on log haz scale 
    target +=  lhaz_epts_event;
    }
    
  if (qevent > 0){
    eta_qpts_event = x_qpts_event * Beta + x_int_qpts_event * Beta_int;;
    lhaz_qpts_event = bspline_log_haz(eta_qpts_event, basis_qpts_event, coefs);
    target +=  quadrature_log_surv(qwts_event, lhaz_qpts_event); // log(f(t)) = -H(t) + log(h(t))for uncensored data
  }
  
  
  if (qrcens > 0) {
    eta_qpts_rcens = x_qpts_rcens * Beta + x_int_qpts_rcens * Beta_int;
    lhaz_qpts_rcens = bspline_log_haz(eta_qpts_rcens, basis_qpts_rcens, coefs);
    target +=  quadrature_log_surv(qwts_rcens, lhaz_qpts_rcens); // log(S(t)) = - H(t)for right censored data
    }

}

generated quantities{
  // Predicting the survival time on the new/test dataset
      vector[nnew] survival_prob;  // 
      
      vector[qevent_new] eta_epts_event_new = x_new_qpts_event * Beta + x_new_int_qpts_event * Beta_int;
      vector[qevent_new] lhaz_epts_event_new = bspline_log_haz(eta_qpts_event_new, basis_qpts_event_new, coefs);
      vector[qevent_new] quadrature_log_surv_qwtsindiv = - (eta_qpts_event_new .* exp(lhaz_epts_event_new));
      matrix[qnodes, nnew] quadrature_log_surv_indiv = to_matrix(quadrature_log_surv_qwtsindiv, qnodes, nnew);
      
      for (n in 1:nnew){
        survival_prob[n] = exp(sum(quadrature_log_surv_indiv[,n]));
      }
  }

