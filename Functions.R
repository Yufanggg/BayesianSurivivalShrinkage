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
Bayesian_Survival_model_check <- function(bayesian_model_fit, criteria = c("MCMC chain", "MCMC autocorrelations")) {
  # Extract posterior samples
  posterior_samples <- extract(bayesian_model_fit)
  
  for (i in 1:length(posterior_samples)){
    # Convert posterior samples to mcmc object
    mcmc_samples <- mcmc(as.matrix(posterior_samples[[i]]))
    # Plot MCMC chains
    plot(mcmc_samples)
    # Plot MCMC autocorrelations
    autocorr.plot(mcmc_samples)
  }
}

# Function to extract the results from the Bayesian model that are useful for model comparison.
# @param bayesian_model_fit: a model fitted by the Bayesian method.
# @param criteria: the results to be extracted from the fitted model.
# @return: a list including (1) the estimated coefficients for the design matrix (i.e., Beta_bayesian_est);
#         (2) predictive survival probability for the testing dataset (i.e., sp); and (3) the information 
#         of selected/non-selected variables (i.e., variableSelection): TRUE if selected, otherwise FALSE.

Bayesian_Survival_result_Extract <- function(bayesian_model_fit, criteria = c("DesignCoefficients", "Prediction_SurvivalProb", "variableSelection")) {
  # Ensure the Output object is available
  Output <- summary(bayesian_model_fit)$summary
  
  # Initialize the result list
  model_result <- list()
  
  if ("DesignCoefficients" %in% criteria) {
    Beta_bayesian_est <- Output[grep("^Beta", rownames(Output)), "mean", drop = FALSE]
    model_result$Beta_bayesian_est <- Beta_bayesian_est
  }
  
  if ("Prediction_SurvivalProb" %in% criteria) {
    sp <- Output[grep("^survival_prob", rownames(Output)), "mean", drop = FALSE]
    model_result$sp <- sp
  }
  
  if ("variableSelection" %in% criteria) {
    Beta_bayesian_est_LB <- Output[grep("^Beta", rownames(Output)), "2.5%", drop = FALSE]
    Beta_bayesian_est_UB <- Output[grep("^Beta", rownames(Output)), "97.5%", drop = FALSE]
    includes_zero <- (Beta_bayesian_est_LB <= 0) & (Beta_bayesian_est_UB >= 0)
    model_result$variableSelection <- !includes_zero 
  }
  
  return(model_result)
}


# Function to evaluate the Bayesian survival model with certain criteria
# @model_result: a list of model results extracted from the Bayesian survival model via Bayesian_Survival_result_Extract function.
# @ComparsionValues: a list of ground-truth/model results extracted from other models, including the following components: Beta_reference, test_t, test_status, and Selection_reference.
# @return: a list of model performance metrics including MSE_est, IBS, C_index, and FDR.
Model_performance_eval <- function(model_result, ComparsionValues, Criteria = c("MSE_est", "IBS", "C_index", "Selection")) {
  
  model_metric <- list()
  
  if ("MSE_est" %in% Criteria) {
    Beta_reference <- ComparsionValues$Beta_reference
    # EVAL 1: Evaluating the model performance at the level of Beta's.
    Beta_bayesian_est <- model_result$Beta_bayesian_est # Extract the Betas
    MSE_est <- (Beta_bayesian_est - Beta_reference) ** 2 # MSE_est is a vector with the length being the number of Betas
    model_metric$rMSE <- MSE_est
  }
  
  if ("IBS" %in% Criteria) {
    # EVAL 2: Evaluating the model performance at the level of prediction
    test_t <- ComparsionValues$test_t
    test_status <- ComparsionValues$test_status
    
    sp <- model_result$sp
    Brier_score <- (sp - test_status)**2 # IBrier-score is a scalar
    model_metric$Brier_score <- mean(Brier_score)
  }
  
  if ("C_index" %in% Criteria) {
    # EVAL 2: Evaluating the model performance at the level of prediction
    sp <- model_result$sp
    C_index <- Cindex(Surv(test_t, test_status), sp, test_t)
    model_metric$C_index <- C_index
  }
  
  if ("Selection" %in% Criteria) {
    # EVAL 3: Evaluating the model variable selection
    variableSelection_reference <- ComparsionValues$Selection_reference
    selection <- model_result$variableSelection
    FDR <- sum(selection & !variableSelection_reference) / sum(selection) # Only selected in the estimation but not for the ground-truth / selected in the estimation
    model_metric$FDR <- FDR
  }
  
  if (length(model_metric) == 0) {
    stop("Sorry, we cannot evaluate the model without the criteria.")
  }
  
  return(model_metric)
}
