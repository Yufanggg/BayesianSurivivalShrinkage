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

        
      }
parameters {
        real<lower = 0> alpha; // parameters in weibull assumption
        real<lower = 0> sigma; // parameters in weibull assumption
        vector[K] Beta; // coefficients
       }
model {
       // Priors
       alpha ~ lognormal(0, 1);
       sigma ~ lognormal(0, 1);
       Beta ~ normal (0, 2);

       // log-likelihood, represented by [target]
       target += weibull_lpdf(t | alpha, sigma/exp(x*Beta/alpha)); // weibull_lpdf() is at log scale // uncensored data log(f(t))
       target += weibull_lccdf(t_cens | alpha, sigma/exp(x_cens*Beta/alpha)); // weibull_lccdf() is at log scale  // right censored data log(S(t)) //log complementary cumulative distribution function (log-CCDF) for an weibull distribution
     }
generated quantities{
      // Predicting the survival time on the new/test dataset
      matrix[M, N_new] survival_prob;  // Unique time points * rows of new data
      for (m in 1:M){
        survival_prob[m,] = exp(-(uniqueT[m]/sigma)^alpha) + x_new*Beta) //vectorize over the x_new.
      }
    }