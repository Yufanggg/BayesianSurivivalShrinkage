//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
  real gauss_kronrod_quad_Hazard_basline(real Ti, vector T, matrix bSpline_basis, vector coefficients){
    # Original nodes and weights for Gauss-Kronrod quadrature on [-1, 1]
    real gauss_nodes[7] = c(-0.949107912342759, -0.741531185599394, -0.405845151377397, 0.000000000000000, 0.405845151377397, 0.741531185599394, 0.949107912342759);
    real gauss_weights[7] = c(0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469, 0.381830050505119, 0.279705391489277, 0.129484966168870);
    
    #Scaling nodes and weights to the interval [0, 10]
    real locat <- (gauss_nodes + 1) / 2 * 1.3859386 + 0.1299756
    real weights <-gauss_weights * 0.6929693
    
    // vectorize all computation
    vector[num_elements(locat)] u_vec = Ti * (rep_vector(1, num_elements(locat))/2);
    vector [num_elements(locat)] logh_0i_vec = vectorized_logh_0i(u_vec, T, bSpline_basis, coefficients);
    real integral_ = sum(to_vector(weights) .* exp(logh_0i_vec));
/*    print("Ti:", Ti);
    // // Print intermediate values for debugging
    // print("u_vec: ", u_vec);
    // print("logh_0i_vec: ", logh_0i_vec);
    // // Print the integral value
    print("integral_: ", integral_);*/
    
    
    return integral_/2;
    }
    
  vector vectorized_logh_0i(vector u_vec, vector T, matrix bSpline_basis, vector coefficients){
    // define the log baseline hazard function in a vectorized way
    vector[num_elements(u_vec)] logh_0i_vec;
    vector[num_elements(T)] logh_0t = bSpline_basis * coefficients;
    for (i in 1:num_elements(u_vec)){
      int index = find_index(u_vec[i], T);
      if (index != -1){
        logh_0i_vec[i] = logh_0t[index];
        }
      else {
        logh_0i_vec[i] = interpolation(u_vec[i], T, logh_0t);
        }
        }
      return logh_0i_vec;
      }
 real logh_0i(real u, vector T, matrix bSpline_basis, vector coefficients){
   // define the log baseline hazard function in a vectorized way
   vector[num_elements(T)] logh_0t = bSpline_basis * coefficients;
   real logh_0i;
   
   int index = find_index(u, T);
   if (index != -1){
     logh_0i = logh_0t[index];
     }
   else{
     logh_0i = interpolation(u, T, logh_0t);
     }
   return logh_0i;
   }
   
 real interpolation (real u, vector T, vector logh_0t){
   int lower_index = 1;
   real logh_0i_interp;
   int upper_index = num_elements(T);
   for (j in 1:(num_elements(T)-1)){
     if (T[j] <= u && u < T[j + 1]) {
       lower_index = j;
       upper_index = j + 1;
       break;
       }
       }
     real t1 = T[lower_index];
     real t2 = T[upper_index];
     real h1 = logh_0t[lower_index];
     real h2 = logh_0t[upper_index];
     logh_0i_interp = h1 + (u - t1) * (h2 - h1) / (t2 - t1);
     
     return logh_0i_interp;
     }
     
 int find_index(real u, vector t) {
   for (i in 1:num_elements(t)) {
     if (u == t[i]) {
       return i; // Return the index if u is found in t
       }
   }
   return -1; //Return -1 if u is not found in t
   }
   
}
      
data {
  int <lower=0> K; //num covariates
  int <lower=0> N; // num uncensored obs
  int <lower=0> M; // num unique time points
  int <lower=0> O; // num unqiue time points in t
  int <lower=0> Q; // num unqiue time points in t_cens
  
  vector[N] t; // event time (non-strict decreasing)
  matrix[N, K] x; // covariates for uncensored obs
  
  int N_cens; // num censored obs
  vector[N_cens] t_cens; // censoring time
  matrix[N_cens, K] x_cens; // covariates for censored obs
  
  matrix[N + N_cens, 6] bSpline_basis;
  vector[M] uniqueT; //montously increasing unique time points
  vector[O] uniquet; //montously increasing unique time points for t
  vector[Q] uniquet_cens; //montously increasing unique time points for t_cens
  
  int<lower=0> N_new;
  matrix[N_new, K] x_new;
  }
      
parameters {
  vector[6] coefficients;
  vector[K] Beta; // coefficients for design matrix;
  } 
  
model {
  // Priors
  coefficients ~ normal(0, 1);
  Beta ~ normal (0, 2);
  // print("uniqueT:", uniqueT);
  // log-likelihood, represented by [target]
  for (Ti in uniqueT){
    real H_0i_ = gauss_kronrod_quad_Hazard_basline(Ti, uniqueT, bSpline_basis, coefficients);
    real logh_0i_ = logh_0i(Ti, uniqueT, bSpline_basis, coefficients);
    
    if (find_index(Ti, uniquet) != -1) {
      vector[N] H_i = H_0i_ * exp(x*Beta);
      vector[N] logh_i = logh_0i_ + x*Beta;
      target += -H_i + logh_i; // uncensored data log(f(t)) = -H(t) + log(h(t))
      }
    if (find_index(Ti, uniquet_cens) != 1){
      vector[N_cens] H_i = H_0i_ * exp(x_cens*Beta);
      target += -H_i;  // at log scale,  right censored data log(S(t)),  log(S(t)) = -H(t)
      }
    }
  }

generated quantities{ // predicting the survival time on the new/test dataset by Monte Carlo sampling (needs to be verfied)
  matrix[M, N_new] survival_prob; // unqiue time points * rows of new data
  for (m in 1:M){
    real H_0i_ = gauss_kronrod_quad_Hazard_basline(m, uniqueT, bSpline_basis, coefficients);
    print("The time points: ", M);
    print("The value of H_0i_ is: ", H_0i_);
    survival_prob[m,] = to_row_vector(exp(-H_0i_ + x_new*Beta)); //vectorize over the x_new, output is the time points * obs
      }
  }
  