data {
      int <lower=0> K; //num covariates
      int <lower=0> N; // num uncensored obs
      int <lower=0> M; // num unique time points

      vector[N] t; // event time (non-strict decreasing)
      matrix[N, K] x; // covariates for uncensored obs

      int N_cens; // num censored obs
      vector[N_cens] t_cens; // censoring time
      matrix[N_cens, K] x_cens; // covariates for censored obs
      
      vector[M] uniqueT; //montously increasing unique time points
      int<lower=0> N_new;
      matrix[N_new, K] x_new;

    }

parameters {
      real<lower = 0> lambda; // parameters in expoential assumption
      vector[K] Beta; // coefficients
    }
    
model {
    // prior
    lambda ~ lognormal(0, 1);
    Beta ~ normal (0, 2);

    // log-likelihood, represented by [target]

    target += exponential_lpdf(t | lambda * exp(x*Beta)); // uncensored data log(f(t))

    target += exponential_lccdf(t_cens | lambda * exp(x*Beta)); // right censored data log(S(t))
    }
    
generated quantities{
      // Predicting the survival time on the new/test dataset
      matrix[M, N_new] survival_prob;  // Unique time points * rows of new data
      for (m in 1:M){
        survival_prob[m,] = exp(-lambda*uniqueT[m] + x_new*Beta) //vectorize over the x_new.
      }
    }