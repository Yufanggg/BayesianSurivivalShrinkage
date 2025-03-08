set.seed(50234)
library(KMsurv)
library(survival)

## Reference paper: Generating survival times to simulate Cox proportional hazards models with time-varying covariates
SurvivalDataGenerator = function(n_samples, Betas, Design_matrix, BaselinePara, BaselineDistr = "exponential"){
  u = runif(n_samples, 0, 1)
  scale = BaselinePara[1]
  if (scale <= 0) stop("Error: The scale must be greater than 0")
  survival_time = -log(u) / (scale * exp(Design_matrix %*% Betas)) # for exponential baseline distribution
  
  if (BaselineDistr == "weibull") {
    shape = BaselinePara[2]
    if (shape <= 0) stop("Error: The scale must be greater than 0 for weibull baseline distribution")
    survival_time = survival_time^(1/shape)
  }
  
  else if (BaselineDistr == "Gompertz") {
    shape = BaselinePara[2]
    survival_time = (1 / shape) * log(1 + shape * survival_time)
  }
  
  # simulate censored data from a uniform distribution
  censtime <- runif(n_samples, 3, 6)
  
  # indicator variable
  status <- as.numeric(survival_time <= censtime) # the sum of status is the number of event
  
  # observed time
  obstime <- survival_time * status + censtime * (1 - status) 
  
  stan_data <- list(
    N = sum(status == 1),
    t = obstime[status == 1],
    N_cens = sum(status == 0),
    t_cens = obstime[status == 0], # Censoring time used in the simulation
    K = dim(Design_matrix)[2], # number of covariates,
    x = as.matrix(Design_matrix[status == 1,]),
    x_cens = as.matrix(Design_matrix[status == 0,]),
    covariate = Design_matrix,
    obstime = obstime,
    status = status,
    Beta = Betas
  )
  return (stan_data)
}


