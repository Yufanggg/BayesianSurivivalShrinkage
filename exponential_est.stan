data {
      int<lower=0> p; // num main effect
      int<lower=0> q; // num interaction effect
      int <lower=0> N; // num uncensored obs
      int <lower=0> M; // num unique time points
      
      int g1[q];
      int g2[q];

      vector[N] t; // event time (non-strict decreasing)
      matrix[N, p] x; // covariates for uncensored obs
      matrix[N, q] x_int;

      int <lower=0> N_cens; // num censored obs
      vector[N_cens] t_cens; // censoring time
      matrix[N_cens, p] x_cens; // covariates for censored obs
      matrix[N_cens, q] x_int_cens;
      
      vector[M] uniqueT; //montously increasing unique time points
      int<lower=0> N_new;
      matrix[N_new, p] x_new;
      matrix[N_new, q] x_int_new;
      vector[N_new] t_new;

    }

parameters {
      real<lower = 0> lambda; // parameters in expoential assumption
      vector[p] Beta; // coefficients for design matrix;
      vector[q] Beta_int;
      
      real<lower=0.01, upper=1> tau2int;
      real<lower=0>  tau2[p];
      real<lower=0>  gam2[p];
    }
    
model {
    // prior
    lambda ~ lognormal(0, 1);
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

    // log-likelihood, represented by [target]

    target += exponential_lpdf(t | lambda * exp(x*Beta + x_int * Beta_int)); // uncensored data log(f(t))

    target += exponential_lccdf(t_cens | lambda * exp(x_cens*Beta + x_int_cens * Beta_int)); // right censored data log(S(t))
    }
    
generated quantites{
      // Predicting the survival time on the new/test dataset
      matrix[N_new, N_new] survival_prob;  // Unique time points * rows of new data
      for (n in 1:N_new){
        survival_prob[m,] = exp(-lambda*t_new[n] + x_new*Beta) //vectorize over the x_new.
      }
    }