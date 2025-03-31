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
      message("We assume that the baseline hazard function is exponentially distributed.")
      bayesian_model <- stan_model("./stan_file/exponential_est.stan")
    }
    
    else if (baseline_assumption == "weibull") {
      # compile the model
      message("We assume that the baseline hazard function is weibully distributed.")
      bayesian_model <- stan_model("./stan_file/weibull_est.stan")
    }
    
    else if (baseline_assumption == "bSplines") {
      message("We utilized B-splines to estimate the log baseline hazard function.")
      
      # compile the model
      bayesian_model <- stan_model("./stan_file/bSpline_est.stan")
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


# Generate a simulate dateset
DataGenerator <- function(n_samples = 210, n_features = 10) {
  sim_data <- list()
  
  # generate design matrix in the formate of dataframe
  # Generate covariance matrix
  Sigma <- matrix(0.3, nrow = n_features, ncol = n_features)
  diag(Sigma) <- 1
  
  design_matrix <- mvrnorm(n_samples, mu = rep(0, n_features), Sigma = Sigma)
  design_matrix <- as.data.frame(design_matrix)
  names(design_matrix) <- paste0("Var", 1:n_features)
  
  # Generate interaction effects
  for (i in 1:(n_features - 1)) {
    for (j in (i + 1):n_features) {
      interaction_term <- design_matrix[, i] * design_matrix[, j]
      colname <- paste0("Var", i, ":Var", j)
      design_matrix[[colname]] <- interaction_term
    }
  }

  print(names(design_matrix))
  
  # generate the cofficents of design matrix
  # Main effects
  beta_main <- c(0.5, -0.4, 0.3, 0.2, 0.1, rep(0, 5))
  
  # Interaction effects
  beta_inter <- matrix(0, nrow = n_features, ncol = n_features)
  beta_inter[1, 2] <- 0.3
  beta_inter[1, 3] <- 0.2
  beta_inter[1, 5] <- -0.1
  beta_inter[2, 3] <- -0.1
  beta_inter[2, 4] <- 0.3
  beta_inter[3, 4] <- 0.2
  beta_inter[1, 6] <- 0.3
  beta_inter[1, 9] <- 0.2
  beta_inter[3, 6] <- -0.2
  beta_inter[4, 7] <- 0.1
  beta_inter[6, 7] <- 0.3
  beta_inter[6, 8] <- 0.1
  beta_inter[9, 10] <- -0.2
  
  beta_inter_t = t(beta_inter)
  Beta = c(beta_main, beta_inter_t[lower.tri(beta_inter_t, diag = FALSE)])
  names(Beta) = names(design_matrix)
  sim_data$Beta = Beta
  
  
  # Generate simulated survival data
  survival_data <- simsurv(
    dist = "weibull",
    lambdas = 2,
    # scale
    gammas = 5,
    # shape for weibull
    betas = Beta,
    x = design_matrix,
    mixture = FALSE,
    maxt = 5  # Maximum follow-up time
  )
  
  survival_data <- survival_data |>
    mutate(
      censtime = runif(n_samples, 0.5, 2),
      status = as.numeric(eventtime <= censtime),
      obstime = pmin(eventtime, censtime)
    ) 
  survival_data <- survival_data[, c("status", "obstime")]
  
  dataset <- cbind(survival_data, design_matrix)
  sim_data$dataset <- dataset
  
  return(sim_data)
}


cross_simulation <- function(whole_dataset, baseline_modelling = "weibull", num_folds = 5){
  # create the cross-validation folds
  folds <- createFolds(whole_dataset$status, k = num_folds, list = TRUE, returnTrain = FALSE)
  cross_val_metric <- list()
  rMSE_s = matrix(nrow = 5, ncol = 55); Brier_scores =  rep(NA, 5); C_indices =  rep(NA, 5); FDR_vals =  rep(NA, 5)
  for (i in 1:num_folds){
    fold_indices <- folds[[i]]
    training_data = whole_dataset[-fold_indices, ]
    testing_data = whole_dataset[fold_indices, ]
    stan_data <- stan_data_Constructer(training_dataset = training_data, testing_dataset = testing_data, baseline_modelling = baseline_modelling, obs_window = 5)
    model_fit <- Bayesian_Survival_model(stan_data = stan_data, baseline_assumption = baseline_modelling)
    
    #extract the info from the bayesian model fit
    model_result <- Bayesian_Survival_result_Extract(bayesian_model_fit = model_fit)
    
    ComparsionValues <- list(
      Beta_reference = c(beta, rep(0, 45)),
      test_t = testing_data$obstime,
      test_status = testing_data$status,
      Selection_reference = c(rep(TRUE, 10), rep(FALSE, 45))
    )
    model_metric <- Model_performance_eval(model_result = model_result, ComparsionValues)
    rMSE_s[i, ] <- model_metric$rMSE
    Brier_scores[i] <- model_metric$Brier_score
    C_indices[i] <- model_metric$C_index
    FDR_vals[i] <- model_metric$FDR
  }
  cross_val_metric$rMSE_s <- rMSE_s
  cross_val_metric$Brier_scores <- Brier_scores
  cross_val_metric$C_indices <- C_indices
  cross_val_metric$FDR_vals <- FDR_vals
  
  return(cross_val_metric)
}



# Function to construct stan_data for model fitting
# @para: training_dataset: 
stan_data_Constructer <- function(training_dataset, withPrediction = FALSE, testing_dataset = NULL, baseline_modelling = "bSplines", obs_window = 5){
  if (withPrediction) {
    if (is.null(testing_dataset)) {
      stop("For predictions, the testing_dataset cannot be empty!")
    }
    if (baseline_modelling == "exponential") {
      source("./function_file/exponential_stan_constructorwithPred.R")
      stan_data = stan_exponential_data_Constructer(training_dataset = training_dataset, testing_dataset = testing_dataset)
    }
    
    if (baseline_modelling == "weibull") {
      source("./function_file/weibull_stan_constructorwithPred.R")
      stan_data = stan_weibull_data_Constructer(training_dataset = training_dataset, testing_dataset = testing_dataset)
    }
    
    if (baseline_modelling == "bSplines") {
      source("./function_file/bSpline_stan_constructorwithPred.R")
      stan_data = stan_bSpline_data_Constructer(training_dataset = training_dataset, testing_dataset = testing_dataset,
        obs_window = 5)
    }
  }
  else{
    if (baseline_modelling == "exponential") {
      source("./function_file/exponential_stan_constructornoPred.R")
      stan_data = stan_exponential_data_Constructer1(training_dataset = training_dataset)
    }
    
    if (baseline_modelling == "weibull") {
      source("./function_file/weibull_stan_constructornoPred.R")
      stan_data = stan_weibull_data_Constructer1(training_dataset = training_dataset)
    }
    
    if (baseline_modelling == "bSplines") {
      source("./function_file/bSpline_stan_constructornoPred.R")
      stan_data = stan_bSpline_data_Constructer1(training_dataset = training_dataset,
                                                obs_window = 5)
    }
  }
  return(stan_data)
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

Bayesian_Survival_result_Extract <- function(bayesian_model_fit,
                                             criteria = c("DesignCoefficients",
                                                          "Prediction_SurvivalProb",
                                                          "variableSelection",
                                                          "baseline")) {
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
    includes_zero <- (Beta_bayesian_est_LB <= 0) &
      (Beta_bayesian_est_UB >= 0)
    model_result$variableSelection <- !includes_zero
  }
  
  if ("baseline" %in% criteria) {
    if (model == "exponential") {
      lambda <- Output[grep("^lambda", rownames(Output)), "50%", drop = FALSE]
      model_result$baselinePara <- lambda
      # Extract posterior samples
      posterior_samples <- extract(bayesian_model_fit,
                                   pars = "lambda",
                                   permuted = TRUE)
      
      # Posterior samples for beta
      lambda_samples <- posterior_samples$lambda
      # Plot the posterior distribution
      hist(
        lambda_samples,
        breaks = 30,
        main = "Posterior Distribution of lambda",
        xlab = "Lambda"
      )
    }
    if (model == "weibull") {
      lambda <- Output[grep("^lambda", rownames(Output)), "50%", drop = FALSE]
      shape <- Output[grep("^shape", rownames(Output)), "50%", drop = FALSE]
      model_result$baselinePara <- c(lambda, shape)
      
      # Extract posterior samples
      posterior_samples <- extract(bayesian_model_fit,
                                   pars = "lambda",
                                   permuted = TRUE)
      
      # Posterior samples for beta
      lambda_samples <- posterior_samples$lambda
      # Plot the posterior distribution
      hist(
        lambda_samples,
        breaks = 30,
        main = "Posterior Distribution of lambda",
        xlab = "lambda"
      )
      
      # Extract posterior samples
      posterior_samples <- extract(bayesian_model_fit,
                                   pars = "shape",
                                   permuted = TRUE)
      
      # Posterior samples for beta
      shape_samples <- posterior_samples$shape
      # Plot the posterior distribution
      hist(
        shape_samples,
        breaks = 30,
        main = "Posterior Distribution of shape",
        xlab = "shape"
      )
    }
    if (model == "bSplines") {
      coefs <- Output[grep("^coefs", rownames(Output)), "mean", drop = FALSE]
      model_result$baselinePara = coefs
      
      # Extract posterior samples
      posterior_samples <- extract(bayesian_model_fit,
                                   pars = "coefs",
                                   permuted = TRUE)
      
      # Posterior samples for beta
      coefs_samples <- posterior_samples$coefs
      # Plot the posterior distribution
      hist(
        coefs_samples[, 1],
        breaks = 30,
        main = "Posterior Distribution of lambda",
        xlab = "coefs[1]"
      )
      hist(
        coefs_samples[, 2],
        breaks = 30,
        main = "Posterior Distribution of lambda",
        xlab = "coefs[2]"
      )
      hist(
        coefs_samples[, 3],
        breaks = 30,
        main = "Posterior Distribution of lambda",
        xlab = "coefs[3]"
      )
    }
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
    C_index <- rcorr.cens(sp, Surv(test_t, test_status))[["C Index"]] # rcorr.cens
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


# Function for model evaluation visualization
# @model_metric: a list, including rMSE_s, Brier_scores, C_indices, C_indices
#   rMSE_s is a matrix with [nFold, nParams], Brier_scores, C_indices, and FDR_vals
#   are vectors with the length of nFold.
metric_visualization <- function(cross_val_metric){
  rMSE_s <- cross_val_metric$rMSE_s
  Brier_scores <- cross_val_metric$Brier_scores
  C_indices <- cross_val_metric$C_indices
  FDR_vals <- cross_val_metric$FDR_vals

  # Adjust the margins
  par(mfrow = c(2, 2))
  #-----------------------------
  # visual the result of rMSE_s
  #-----------------------------
  mean_values <- colMeans(rMSE_s)
  sd_values <- apply(rMSE_s, 2, sd)
  # Plot the mean values
  plot(mean_values, type = "o", col = "blue", xlab = "estimated Betas", ylab = "estimated Value of Betas", ylim = c(0, 2), main = "MSE for Beta estimation", pch = 16)
  # Add error bars
  arrows(1:55, mean_values - sd_values, 1:55, mean_values + sd_values, angle = 90, code = 3, length = 0.05, col = "red")

  # Add a legend
  legend("topright", legend = c("Mean Values", "Error Bars"), col = c("blue", "red"), pch = 16, lty = 1)

  #-----------------------------
  # visual the result of Brier_scores
  #-----------------------------
  boxplot(Brier_scores, main = "Brier_scores for the model prediction", ylab = "Brier_scores", col = "blue")


  #-----------------------------
  # visual the result of C_indices
  #-----------------------------
  boxplot(C_indices, main = "C_indices for the model prediction", ylab = "C_indices", col = "green")

  #-----------------------------
  # visual the result of FDR_vals
  #-----------------------------
  boxplot(FDR_vals, main = "FDR for the variable selection", ylab = "FDR", col = "orange")


}
