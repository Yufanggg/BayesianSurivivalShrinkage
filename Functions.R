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
      bayesian_model <- stan_model("./exponential.stan")
    }
    
    else if (baseline_assumption == "weibull") {
      # compile the model
      bayesian_model <- stan_model("./weibull.stan")
    }
    
    else if (baseline_assumption == "bSplines") {
      message("We utilized B-splines to estimate the baseline cumulative hazard function.")
      time_combined <- sort(unique(c(stan_data$t, stan_data$t_cens)))
      notes = sort(runif(5, min(time_combined), max(time_combined)))
      bSpline_basis <- bSpline(time_combined, knots = notes, degree = 1, intercept = FALSE) # The B-spline basis is calculated using the method implemented in the splines2 package
      
      
      # Out the corresponding information in stan data
      stan_data$bSpline_basis <- bSpline_basis
      stan_data$M = length(time_combined)
      stan_data$uniqueT = time_combined
      
      library(mvQuad)
      
      # Create a Gauss-Kronrod grid
      grid <- createNIGrid(dim = 1, type = "GHe", level = 15)
      
      # Get nodes and weights
      locates <- getNodes(grid)
      weights <- getWeights(grid)
      stan_data$locates <- as.vector(locates)
      stan_data$weights <- as.vector(weights)
      
      # compile the model
      bayesian_model <- stan_model("./bSpline_estimation.stan")
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