
find_main_effect_indices <- function(interaction_effects, main_effects) {
  interaction_indices <- list()
  
  for (interaction in interaction_effects) {
    features <- unlist(strsplit(interaction, ":"))
    indices <- match(features, main_effects)
    interaction_indices[[interaction]] <- indices
  }
  
  return(interaction_indices)
}


g <- function(main_indices_for_int, index){
  g_val = c()
  for (i in 1:length(main_indices_for_int)){
    #print(main_indices_for_int[[i]][index])
    g_val = c(g_val, main_indices_for_int[[i]][index])
  }
  return(g_val)
}


Bayesian_Survival_includingbaseline <- function(stan_data, baseline_assumption = "exponential", school = "Bayesian", 
                                                niter = 10000, 
                                                nwarmup = 1000,
                                                thin = 10,
                                                chains = 1) {
  if (school == "Bayesian") {
    if (baseline_assumption == "exponential") {
      # compile the model
      bayesian_model <- stan_model("./exponential_est.stan")
    }
    
    else if (baseline_assumption == "weibull") {
      # compile the model
      bayesian_model <- stan_model("./weibull_est.stan")
    }
    
    else if (baseline_assumption == "bSplines") {
      message("We utilized B-splines to estimate the baseline cumulative hazard function.")
      
      
      # compile the model
      bayesian_model <- stan_model("./bSpline_est.stan")
    }
    
    # Model fitting and summary
    bayesian_model_fit <- suppressWarnings(
      sampling(
        bayesian_model,
        data = stan_data,
        iter = niter,
        warmup = nwarmup,
        thin = 10,
        chain = 1
      )
    )
    
    return(bayesian_model_fit)
  }
}