#Schoenfeld residuals 
compute_baseline_hazard <- function(beta, X, time, status) {
  # beta: numeric vector of Cox coefficients
  # X: matrix or data frame of covariates
  # time: numeric vector of observed times
  # status: numeric or logical vector (1 = event, 0 = censored)
  
  # Compute risk scores
  risk_scores <- exp(as.matrix(X) %*% beta)
  
  # All unique times (including censored)
  all_times <- sort(unique(time))
  
  baseline_hazard <- numeric(length(all_times))
  last_hazard <- 0
  
  
  for (i in seq_along(all_times)) {
    t <- all_times[i]
    
    # Number of events at time t
    d_i <- sum(time == t & status == 1)
    
    if (d_i > 0) {
      # Risk set: subjects with time >= t
      risk_set <- sum(risk_scores[time >= t])
      
      # Breslow estimator for hazard at time t
      last_hazard <- d_i / risk_set
    }
    
    
    
    # If no event, carry forward previous hazard
    baseline_hazard[i] <- last_hazard
  }
  # Return as data frame
  return(baseline_hazard)
}


Schoenfeld_resid <- function(beta, X, times, status) {
  # Extract event times and their original indices
  event_times <- times[status == 1]
  event_indices <- which(status == 1)
  sorted_event_indices <- order(event_times)
  sorted_event_times <- event_times[sorted_event_indices]
  sorted_event_indices <- event_indices[sorted_event_indices]
  #
  # event_timPoSorted <- sort(times[status == 2])
  # event_timPoSortedIndex <- sort(times[status == 2], index.return = TRUE)$ix
  # indices_list <- lapply(event_timPoSorted, function(t_i) which(times >= t_i))
  # R_ij <- matrix(data = NA, nrow = length(indices_list), ncol = ncol(X))
  
  Output <- list()
  # Build risk sets: indices where times >= t_i
  indices_list <- lapply(sorted_event_times, function(t_i)
    which(times >= t_i))
  
  # Initialize residual matrix
  R_ij <- matrix(NA, nrow = length(indices_list), ncol = ncol(X))
  
  # Loop over covariates and event times
  
  # if (length(unique(event_times)) == length(event_times)){
  for (i in 1:length(indices_list)) {
    # print(sorted_event_times[i])
    for (j in 1:ncol(X)) {
      R_i <- indices_list[[i]]
      X_k_all <- X[R_i, ]
      exp_part <- exp(X_k_all %*% beta)
      deno <- sum(exp_part)
      
      X_kj_all <- X[R_i, j]
      nemo <- sum(X_kj_all * exp_part)
      E_ij <- nemo / deno
      
      t_iIndex <- sorted_event_indices[i]
      R_ij[i, j] <- X[t_iIndex, j] - E_ij
    }
  }
  
  Output$R_ij <- R_ij
  Output$sorted_event_times <- t

return(Output)
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

# Function to fit a Bayesian survival model with different baseline assumptions
# @param: stan_data, a list including items for the corresponding stan model
Bayesian_Survival_model <- function(stan_data, baseline_modelling = "exponential", withPrediction = TRUE, shrinkage = TRUE, log_likSaving = FALSE, havingInt = TRUE,
                                                niter = 10000, 
                                                nwarmup = 1000,
                                                thin = 10,
                                                chains = 1) {
  if (!shrinkage){
    message("we only applied the non shrinage code for bSpline baseline hazard or PH models.")
  }
  
  
  if (!withPrediction) {
    if (baseline_modelling == "exponential") {
      # compile the model
      message("We assume that the baseline hazard function is exponentially distributed.")
      if (log_likSaving == TRUE){
        bayesian_model <- rstan::stan_model("./stan_file/exponential_est_loglik.stan")
      } else {
        bayesian_model <- rstan::stan_model("./stan_file/exponential_est.stan")
        }
    }
    
    else if (baseline_modelling == "weibull") {
      # compile the model
      message("We assume that the baseline hazard function is weibully distributed.")
      if (log_likSaving == TRUE){
        bayesian_model <- rstan::stan_model("./stan_file/weibull_est_loglik.stan")
      } else {
        bayesian_model <- rstan::stan_model("./stan_file/weibull_est.stan")
      }
     
    }
    
    else if (baseline_modelling == "bSplines") {
      message("We utilized B-splines to estimate the log baseline hazard function.")
      
      # compile the model
      if (log_likSaving == TRUE) {
        if (havingInt == TRUE) {
          bayesian_model <- rstan::stan_model("./stan_file/bSpline_est_loglik.stan")
        } else {
          bayesian_model <- rstan::stan_model("./stan_file/bSpline_est_loglik_noInt.stan")
        }
        
      } else {
        if (havingInt == TRUE){
          bayesian_model <- rstan::stan_model("./stan_file/bSpline_est.stan")
        } else {
          bayesian_model <- rstan::stan_model("./stan_file/bSpline_est_noInt.stan")
        }
        
      }
    }
    
    else if (baseline_modelling == "none"){
      message("we used the partical likelihood to estimate the cofficients of covariates")
      if (havingInt == TRUE){
        bayesian_model <- rstan::stan_model("./stan_file/PH_est.stan")
      } else {
        bayesian_model <- rstan::stan_model("./stan_file/PH_est_noInt.stan")
      }
      
    }
    
  } else {
    if (baseline_modelling == "exponential") {
      # compile the model
      message("We assume that the baseline hazard function is exponentially distributed.")
      bayesian_model <-
        rstan::stan_model("./stan_file/exponential_estEpre.stan")
    }
    
    else if (baseline_modelling == "weibull") {
      # compile the model
      message("We assume that the baseline hazard function is weibully distributed.")
      bayesian_model <- rstan::stan_model("./stan_file/weibull_estEpre.stan")
    }
    
    else if (baseline_modelling == "bSplines") {
      message("We utilized B-splines to estimate the log baseline hazard function.")
      
      # compile the model
      if (shrinkage){
        bayesian_model <- rstan::stan_model("./stan_file/bSpline_estEpre.stan")
      } else {
        message("No shrinkaged bSplines modelling!")
        bayesian_model <- rstan::stan_model("./stan_file/bSpline_estEprenoShrinkage.stan")
      }
      
    }
    
    else if (baseline_modelling == "none"){
      message("we used the partical likelihood to estimate the cofficients of covariates")
      # compile the model
      if (shrinkage){
        bayesian_model <- rstan::stan_model("./stan_file/PH_estEpre.stan")
      } else {
        message("No shrinkaged partical likelihood modelling!")
        bayesian_model <- rstan::stan_model("./stan_file/PH_estEprenoShrinkage.stan")
      }
    }
  }
  
    
    # Model fitting
    bayesian_model_fit <- suppressWarnings(
      rstan::sampling(
        bayesian_model,
        data = stan_data,
        iter = niter,
        warmup = nwarmup,
        thin = thin,
        chain = chains
      )
    )
    
    return(bayesian_model_fit)
}



# Generate a simulate dateset
DataGenerator <- function(n_samples, n_features, multicolinearity = 0.3, TimeVarying = FALSE) {
  sim_data <- list()
  
  # generate design matrix in the format of dataframe
  # Generate covariance matrix
  Sigma <- matrix(multicolinearity, nrow = n_features, ncol = n_features)
  diag(Sigma) <- 1
  
  design_matrix <- MASS::mvrnorm(n_samples, mu = rep(0, n_features), Sigma = Sigma)
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

  # print(names(design_matrix))
  
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
  
  if (TimeVarying == FALSE) {
    # Generate simulated survival data
    survival_data <- simsurv::simsurv(
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
    
    # survival_data <- survival_data |>
    #   mutate(
    #     censtime = runif(n_samples, 0.5, 5),
    #     status = as.numeric(eventtime <= censtime),
    #     obstime = pmin(round(eventtime, 4), censtime)
    #   )
    # survival_data <- survival_data[, c("status", "obstime")]
    
    dataset <- cbind(survival_data, design_matrix)
    sim_data$dataset <- dataset
  } else {
    # define a time-varying effec function
    Beta_timevarying <- function(t) {
      -0.5 + 0.3 * log(t + 1e-5)  # avoid log(0)
    }
    
    survival_data <- simsurv::simsurv(
      dist = "weibull",
      lambdas = 2,
      # scale
      gammas = 10,
      # shape for weibull
      betas = Beta,
      tde = c(Var1 = -0.5),
      tdefunction = function(t) log(t + 1e-5),  # time transformation,
      x = design_matrix,
      mixture = FALSE,
      maxt = 5  # Maximum follow-up time
    )
    dataset <- cbind(survival_data, design_matrix)
    sim_data$dataset <- dataset
    
  }
  
  
  
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

Bayesian_Survival_result_Extract <- function(bayesian_model_fit, model_type,
                                             criteria = c("DesignCoefficients",
                                                          "Prediction_SurvivalProb",
                                                          "variableSelection",
                                                          "baseline"), uncertainty = FALSE) {
  # Ensure the Output object is available
  Output <- summary(bayesian_model_fit)$summary
  
  # Initialize the result list
  model_result <- list()
  
  if ("variableSelection" %in% criteria) {
    Beta_bayesian_est_LB <- Output[grep("^Beta", rownames(Output)), "2.5%", drop = FALSE]
    Beta_bayesian_est_UB <- Output[grep("^Beta", rownames(Output)), "97.5%", drop = FALSE]
    includes_zero <- (Beta_bayesian_est_LB <= 0) &
      (Beta_bayesian_est_UB >= 0)
    model_result$variableSelection <- !includes_zero
  }
  
  if (uncertainty == FALSE){
    if ("DesignCoefficients" %in% criteria) {
      Beta_bayesian_est <- Output[grep("^Beta", rownames(Output)), "mean", drop = FALSE]
      model_result$Beta_bayesian_est <- Beta_bayesian_est
    }
    
    if ("Prediction_SurvivalProb" %in% criteria) {
      if (model_type != "none"){
        eta_pred <- Output[grep("^eta_new", rownames(Output)), "50%", drop = FALSE]
        model_result$eta_pred <- eta_pred
        sp <- Output[grep("^survival_prob", rownames(Output)), "50%", drop = FALSE]
        model_result$sp <- sp
      } else {
        eta_pred <- Output[grep("^eta_new", rownames(Output)), "50%", drop = FALSE]
        model_result$eta_pred <- eta_pred
      }
      
    }
    
    if ("baseline" %in% criteria) {
      if (model_type == "exponential") {
        lambda <- Output[grep("^lambda", rownames(Output)), "50%", drop = FALSE]
        model_result$baselinePara <- lambda
      } else if (model_type == "weibull") {
        lambda <- Output[grep("^lambda", rownames(Output)), "50%", drop = FALSE]
        shape <- Output[grep("^shape", rownames(Output)), "50%", drop = FALSE]
        model_result$baselinePara <- c(lambda, shape)
      } else if (model_type == "bSplines") {
        coefs <- Output[grep("^coefs", rownames(Output)), "mean", drop = FALSE]
        model_result$baselinePara = coefs
      } else {
        # No actions for 
      }
    }
  } else {
    
    if ("DesignCoefficients" %in% criteria) {
      Beta_bayesian_est <- Output[grep("^Beta", rownames(Output)), "mean", drop = FALSE]
      model_result$Beta_bayesian_est <- Beta_bayesian_est
      
      model_result$Beta_bayesian_estLB <- Output[grep("^Beta", rownames(Output)), "2.5%", drop = FALSE]
      model_result$Beta_bayesian_estUB <- Output[grep("^Beta", rownames(Output)), "97.5%", drop = FALSE]
    }
    
    if ("Prediction_SurvivalProb" %in% criteria) {
      eta_pred <- Output[grep("^eta_new", rownames(Output)), "50%", drop = FALSE]
      model_result$eta_pred <- eta_pred
      
      model_result$eta_predLB <- Output[grep("^eta_new", rownames(Output)), "2.5%", drop = FALSE]
      model_result$eta_predUB <- Output[grep("^eta_new", rownames(Output)), "97.5%", drop = FALSE]
      
      if (model_type != "none"){
        sp <- Output[grep("^survival_prob", rownames(Output)), "50%", drop = FALSE]
        model_result$sp <- sp
        
        model_result$spLB <- Output[grep("^survival_prob", rownames(Output)), "2.5%", drop = FALSE]
        model_result$spUB <- Output[grep("^survival_prob", rownames(Output)), "97.5%", drop = FALSE]
      }
      
    }
    
    if ("baseline" %in% criteria) {
      if (model_type == "exponential") {
        lambda <- Output[grep("^lambda", rownames(Output)), "50%", drop = FALSE]
        lambdaLB <- Output[grep("^lambda", rownames(Output)), "2.5%", drop = FALSE]
        lambdaUB <- Output[grep("^lambda", rownames(Output)), "97.5%", drop = FALSE]
        
        model_result$baselinePara <- list(lambda, lambdaLB, lambdaUB)
      } else if (model_type == "weibull") {
        lambda <- Output[grep("^lambda", rownames(Output)), "50%", drop = FALSE]
        lambdaLB <- Output[grep("^lambda", rownames(Output)), "2.5%", drop = FALSE]
        lambdaUB <- Output[grep("^lambda", rownames(Output)), "97.5%", drop = FALSE]
        
        shape <- Output[grep("^shape", rownames(Output)), "50%", drop = FALSE]
        shapeLB <- Output[grep("^shape", rownames(Output)), "2.5%", drop = FALSE]
        shapeUB <- Output[grep("^shape", rownames(Output)), "97.5%", drop = FALSE]
        
        model_result$baselinePara <- list(lambda, shape, lambdaLB, lambdaUB, shapeLB, shapeUB)
      } else if (model_type == "bSplines") {
        coefs <- Output[grep("^coefs", rownames(Output)), "mean", drop = FALSE]
        coefsLB <- Output[grep("^coefs", rownames(Output)), "2.5%", drop = FALSE]
        coefsUB <- Output[grep("^coefs", rownames(Output)), "97.5%", drop = FALSE]
        
        model_result$baselinePara = list(coefs, coefsLB, coefsUB)
      } else {
        # No actions for 
      }
    }
  }
  
  return(model_result)
}


# Function for cross-validation


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


# Cross valiating
cross_validation <- function(whole_dataset, baseline_modelling = "bSplines", num_folds = 5, obs_window = 3415){
  # create the cross-validation folds
  folds <- createFolds(whole_dataset$status, k = num_folds, list = TRUE, returnTrain = FALSE)
  cross_val_metric <- list()
  Brier_scores =  rep(NA, num_folds); C_indices =  rep(NA, num_folds)
  for (i in 1:num_folds){
    fold_indices <- folds[[i]]
    training_data = whole_dataset[-fold_indices, ]
    testing_data = whole_dataset[fold_indices, ]
    stan_data_bSplines_cv <- stan_data_Constructer(training_dataset = training_data, withPrediction = TRUE, testing_dataset = testing_data, baseline_modelling = "bSplines", obs_window = obs_window)
    model_fit_bSplines <- Bayesian_Survival_model(stan_data = stan_data_bSplines_cv, withPrediction = TRUE, baseline_assumption = "bSplines")
   
    #model diagnotics
    diagnostics_bSplines <- summary(model_fit_bSplines)$summary[, c("Rhat", "n_eff")]
    cat("The covergence for", i, "\n")
    print(diagnostics_bSplines)
    
    #extract the info from the bayesian model fit
    model_result_bSplines <- Bayesian_Survival_result_Extract(bayesian_model_fit = model_fit_bSplines, model_type = "bSplines", criteria = "Prediction_SurvivalProb")
    
    predicted_survP <- model_result_bSplines$sp
    predicted_eta <- model_result_bSplines$eta_pre
    
    real_status <- testing_data$status
    time <- testing_data$obstime
    Brier_scores[i] = mean((predicted_survP - real_status)^2)
    C_indices[i] = rcorr.cens(-predicted_eta, Surv(time, real_status))[['C Index']]# rcorr.cens
    
  }
  cross_val_metric$Brier_scores <- Brier_scores
  cross_val_metric$C_indices <- C_indices
  return(cross_val_metric)
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

# Obtaining the cvl
cross_validated_log_likelihood <- function(model_fit_bSpline, nObs) {
  cvl_s = rep(NA, nObs)
  log_lik_matrix <- extract_log_lik(model_fit_bSpline,
                                    parameter_name = "log_lik",
                                    merge_chains = TRUE) # S*N
  
  
  for (i in 1:ncol(log_lik_matrix)) {
    cvl_s[i] <- max(apply(log_lik_matrix[, -i], 1, sum))
  }
  return(sum(cvl_s))
}

max_log_likelihood <- function(model_fit_bSpline, nObs) {
  cvl_s = rep(NA, nObs)
  log_lik_matrix <- extract_log_lik(model_fit_bSpline,
                                    parameter_name = "log_lik",
                                    merge_chains = TRUE) # S*N
  log_lik <- max(apply(log_lik_matrix, 1, sum))
  return(log_lik)
}

# include the interaction effect for the dataframe, tailed for the real kindney data.
includingInter <- function(dataframe_) {
  #sum-coded all factor variables
  factor_var <- sapply(dataframe_, is.factor)
  for (var in names(dataframe_)[factor_var]) {
    if (var != "status") {
      contrasts(dataframe_[[var]]) <- contr.sum(levels(dataframe_[[var]]))
    }
  }
  numeric_var <- sapply(dataframe_, is.numeric)
  for (var in names(dataframe_)[numeric_var]){
    if (var != "time"){
      dataframe_[[var]] <- (dataframe_[[var]] - mean(dataframe_[[var]], na.rm = TRUE)) / sd(dataframe_[[var]], na.rm = TRUE)
    }
    
  }
  
  # print(colnames(dataframe_)[numeric_var])
  
  # build up interaction effect
  predictors <- subset(dataframe_, select = -c(status, time))
  
  interaction_formula <- as.formula(paste("~ (", paste(colnames(predictors), collapse = " + "), ")^2"))
  interaction_matrix <- model.matrix(interaction_formula, data = dataframe_)
  
  interaction_df <- as.data.frame(interaction_matrix)
  interaction_df <- interaction_df[, -1]
  
  # Exclude interactions involving categorical variables
  all_terms <- colnames(interaction_df)
  interaction_terms <- all_terms[grepl(":", all_terms)]
  # print(interaction_terms)
  factor_vars <- c("Recipientsex", "Donorsex", "Smoking", "InitialOnmachineindicator", "status", "InitialPrimaryDiseaseET_regroup", "Donorcauseofdeath_group")#colnames(dataframe_)[factor_var]
  # print(factor_vars)
  
  interaction_df$status = ifelse(dataframe_$status == "graftloss", 1, 0)
  interaction_df$obstime = dataframe_$time
  
  
  # print(interaction_terms[exclude_terms])
  
  
  return(interaction_df)
}
