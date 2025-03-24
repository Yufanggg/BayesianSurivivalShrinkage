# Function to find indices of main effects for interaction effects
# @param: interaction_effects, a vector with the names of the interaction effects
# @param: main_effects, a vector with the names of the main effects
# @return: a list of the indices of main effects for given interaction effects:
#          interaction_indices[interaction][main effect index] 
find_main_effect_indices <- function(interaction_effects, main_effects) {
  interaction_indices <- list()
  
  for (interaction in interaction_effects) {
    features <- unlist(strsplit(interaction, ":"))
    indices <- match(features, main_effects)
    interaction_indices[[interaction]] <- indices
  }
  
  return(interaction_indices)
}

# Function to extract specific indices from a list of indices
# @param: main_indices_for_int, a list of the indices of main effects for given interaction effects 
# @param: index, the index of first or second main effect to be extract for a given interaction effect
# @return: the index of the first or second (suggested by index) main effect for interaction effects in the list.
g <- function(main_indices_for_int, index){
  g_val = c()
  for (i in 1:length(main_indices_for_int)){
    #print(main_indices_for_int[[i]][index])
    g_val = c(g_val, main_indices_for_int[[i]][index])
  }
  return(g_val)
}

# Function to fit a Bayesian survival model with different baseline assumptions
# @param: stan_data, a list including items for the corresponding stan model
Bayesian_Survival_model <- function(stan_data, baseline_assumption = "exponential", school = "Bayesian", 
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
    
    # Model fitting
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

# Function to check the Bayesian survival model using various criteria
# @param: model, a model fitting by the bayesian method.
# @param: the metric for model fitting criteria
# @return: not return any things but give a few plots
Bayesian_Survival_model_check <- function(model, criteria = c("Bayesian R-square", "MCMC chain", "MCMC autocorrelations")) {
  library(bayesplot)
  
  if ("Bayesian R-square" %in% criteria) {
    Bayes_r2 <- bayes_R2(model)  # Calculate Bayesian R-square
    mcmc_hist(Bayes_r2)  # Plot Bayesian R-square
  }
  
  if ("MCMC chain" %in% criteria) {
    mcmc_trace(as.array(model))  # Plot MCMC chains
  }
  
  if ("MCMC autocorrelations" %in% criteria) {
    mcmc_acf(as.array(model))  # Plot MCMC autocorrelations
  }
}

# Function to extract the results from the bayesian model that are useful for model comparsion.
# @param: model, a model fitting by the bayesian method.
# @param: the results to be extract from the fitted model
# @return: a list including (1) the estimated coefficients for the design matrix (i.e., Beta_bayesian_est);
#         (2) predictive survival probability for the testing dataset (i.e., sp); and (3) the information 
#         of selected/non-selected variables(i.e., variableSelection): TRUE if selected, otherwise FALSE. 
Bayesian_Survival_result_Extract<- function(model, criteria = c("DesignCoefficients", "Prediction_SurvivalProb", "variableSelection")) {
  # Ensure the Output object is available
  Output <- summary(model)$summary
  
  if ("DesignCoefficients" %in% criteria) {
    Beta_bayesian_est <- Output[grep("^Beta", rownames(Output)), "mean", drop = FALSE]
    print(Beta_bayesian_est)
  }
  
  if ("Prediction_SurvivalProb" %in% criteria) {
    sp <- Output[grep("^survival_prob", rownames(Output)), "mean", drop = FALSE]
    print(sp)
  }
  
  if ("variableSelection" %in% criteria) {
    Beta_bayesian_est_LB <- Output[grep("^Beta", rownames(Output)), "2.5%", drop = FALSE]
    Beta_bayesian_est_UB <- Output[grep("^Beta", rownames(Output)), "97.5%", drop = FALSE]
    variableSelection <- (Beta_bayesian_est_LB > 0) & (Beta_bayesian_est_UB < 0)
    print(variableSelection)
  }
  
  model_result <- list(
    Beta_bayesian_est = Beta_bayesian_est,
    sp = sp,
    variableSelection = variableSelection
  )
  return(model_result)
}


# Function to evaluate the  bayesian survival model with certain criteria
# @model_result: a list of model results extracted from the bayesian survival model via Bayesian_Survival_result_Extract function.
# @ComparsionValues: a list of ground-truth/model results extract from other models but organized in the same way with `model_result`.
# @return: a list of model performance metric including 
Model_comparsion <- function(model_result, ComparsionValues, ComparsionCriteria = c("MSE_est", "IBS", "CI", "Selection")){
  
  
  
  model_metric = list(
    MSE_est = MSE_est,
    IBS = IBS,
    CI = CI,
    Selection = Selection,
  )
}