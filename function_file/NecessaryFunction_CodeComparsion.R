# Necessary function for code comparsion
# Generate a simulate dateset
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
    gammas = 10,
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


# Function to construct stan_data for model fitting
# @para: training_dataset: 
stan_data_Constructer <- function(training_dataset, withPrediction = FALSE, testing_dataset = NULL, baseline_modelling = "bSplines", obs_window = 5){
  if (withPrediction) {
    if (is.null(testing_dataset)) {
      stop("For predictions, the testing_dataset cannot be empty!")
    }
    if (baseline_modelling == "bSplines") {
      source("./function_file/bSpline_stan_constructorwithPred.R")
      stan_data = stan_bSpline_data_Constructer(training_dataset = training_dataset, testing_dataset = testing_dataset,
                                                obs_window = 5)
    }
  }
  else{
    if (baseline_modelling == "bSplines") {
      source("./function_file/bSpline_stan_constructornoPred.R")
      stan_data = stan_bSpline_data_Constructer1(training_dataset = training_dataset,
                                                 obs_window = 5)
    }
  }
  return(stan_data)
}

# Function to fit a Bayesian survival model with different baseline assumptions
# @param: stan_data, a list including items for the corresponding stan model
Bayesian_Survival_model <- function(stan_data,
                                    baseline_assumption = "exponential",
                                    school = "Bayesian",
                                    niter = 10000,
                                    nwarmup = 1000,
                                    thin = 10,
                                    chains = 1) {
  if (school == "Bayesian") {
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


Bayesian_Survival_result_Extract <- function(bayesian_model_fit, model_type,
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
    if (model_type == "bSplines") {
      coefs <- Output[grep("^coefs", rownames(Output)), "mean", drop = FALSE]
      model_result$baselinePara = coefs
      
      # # Extract posterior samples
      # posterior_samples <- extract(bayesian_model_fit,
      #                              pars = "coefs",
      #                              permuted = TRUE)
      
      # # Posterior samples for beta
      # coefs_samples <- posterior_samples$coefs
      # # Plot the posterior distribution
      # hist(
      #   coefs_samples[, 1],
      #   breaks = 30,
      #   main = "Posterior Distribution of lambda",
      #   xlab = "coefs[1]"
      # )
      # hist(
      #   coefs_samples[, 2],
      #   breaks = 30,
      #   main = "Posterior Distribution of lambda",
      #   xlab = "coefs[2]"
      # )
      # hist(
      #   coefs_samples[, 3],
      #   breaks = 30,
      #   main = "Posterior Distribution of lambda",
      #   xlab = "coefs[3]"
      # )
    }
  }
  
  return(model_result)
}



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