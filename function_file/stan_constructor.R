# assumptions: all data are right censored data with the observed window being [0, obs_window]
# Construct stan data structure for baseline hazard of exponential and weibull.
stan_data_Constructer_PL <- function (training_dataset, testing_dataset = NULL, havingInt = TRUE){
  
  #----------------------------
  # Prepare data for model fitting
  #-----------------------------
  #----- organize the data regarding predictor
  Sorted_training_dataset <- training_dataset[order(-training_dataset$obstime, -training_dataset$status), ] # observed time points (non-strict decreasing)
  t_points = Sorted_training_dataset[, "obstime"]
  
  design__matrix = Sorted_training_dataset[,!(names(Sorted_training_dataset) %in% c("id", "obstime", "status"))]
  
  column_names = colnames(design__matrix)
  main_names =  column_names[!grepl(":", column_names)] #whether or not having intercept needs to be verified
  X_main = design__matrix[, main_names]
  
  int_names =  column_names[grepl(":", column_names)]
  X_int = design__matrix[, int_names]
  
  p <- dim(X_main)[2]
  q <- dim(X_int)[2]
  
  main_indices_for_int = find_main_effect_indices(int_names, main_names)
  g1 <- g(main_indices_for_int, 1)
  g2 <- g(main_indices_for_int, 2)
  

  
  event_flag <- Sorted_training_dataset$status == 1
  event_values <- t_points[event_flag]
  unique_event_times <- unique(event_values)
  
  # Find the first index in t_points for each unique event time where event_flag is TRUE
  first_indices_events <- sapply(unique_event_times, function(time) {
    which(t_points == time & event_flag)[1]
  })
  
  
  # Find the last index in t_points for each unique event time where event_flag is TRUE
  last_indices_events <- sapply(unique_event_times, function(time) {
    tail(which(t_points == time & event_flag), 1)
  })
  
  
  
  #----------------
  # Construct data
  #----------------
  if (havingInt == TRUE){
    stan_data = list(
      #----- for model fitting --------
      nobs = nrow(Sorted_training_dataset),
      first_indices_events = first_indices_events,
      last_indices_events = last_indices_events,
      unique_nevent = length(unique_event_times),
      p = p,
      q = q,
      
      x = X_main,
      x_int =  X_int,
      
      # link the interaction effect with the corresponding main effects
      g1 = g1,
      g2 = g2
    )
  } else {
    stan_data = list(
      #----- for model fitting --------
      nobs = nrow(Sorted_training_dataset),
      first_indices_events = first_indices_events,
      last_indices_events = last_indices_events,
      unique_nevent = length(unique_event_times),
      p = p,
      x = X_main
    )
  }
  
  
  
  #----------------------------
  # Prepare data for prediction
  #-----------------------------
  if(! is.null(testing_dataset)){
    nnew <- nrow(testing_dataset)
    t_new <- testing_dataset$obstime
    design__matrix_2 = testing_dataset[,!(names(testing_dataset) %in% c("id", "obstime", "status"))]
    
    column_names = colnames(design__matrix_2)
    main_names =  column_names[!grepl(":", column_names) &
                                 column_names != "(Intercept)"] #whether or not having intercept needs to be verified
    X_new_main = design__matrix_2[, main_names]
    
    int_names =  column_names[grepl(":", column_names)]
    X_new_int = design__matrix_2[, int_names]
    
    
    #----------------
    # add data to Constructed data for prediction
    #----------------
    
    stan_data = c(stan_data, list(
      nnew = nrow(testing_dataset),
      x_new = X_new_main,
      x_int_new = X_new_int,
      t_new = testing_dataset[, "obstime"]
    ))
  }
  
  return(stan_data)
}



# assumptions: all data are right censored data with the observed window being [0, obs_window]
# Construct stan data structure for baseline hazard of exponential and weibull.
stan_data_Constructer_noBSpline <- function (training_dataset, testing_dataset = NULL, log_lik_saving = TRUE, havingInt = TRUE) {
  
  #----------------------------
  # Prepare data for model fitting
  #-----------------------------
  #----- organize the data regarding predictor
  
  design__matrix = training_dataset[, !(names(training_dataset) %in% c("id", "obstime", "status"))]
  
  column_names = colnames(design__matrix)
  
  if (havingInt == TRUE) {
    main_names =  column_names[!grepl(":", column_names)] #whether or not having intercept needs to be verified
    X_main = design__matrix[, main_names]
    
    int_names =  column_names[grepl(":", column_names)]
    X_int = design__matrix[, int_names]
    
    p <- dim(X_main)[2]
    q <- dim(X_int)[2]
    
    main_indices_for_int = find_main_effect_indices(int_names, main_names)
    g1 <- g(main_indices_for_int, 1)
    g2 <- g(main_indices_for_int, 2)
    
    # x_event = X_main[training_dataset$status == 1,]
    # x_int_event =  X_int[training_dataset$status == 1,]
    
    
    
    #----------------
    # Construct data
    #----------------
    stan_data = list(
      #----- for model fitting --------
      nevent = nrow(training_dataset[training_dataset$status == 1, ]),
      nrcens = nrow(training_dataset[training_dataset$status == 0, ]),
      t_event = training_dataset[training_dataset$status == 1, "obstime"],
      t_rcens = training_dataset[training_dataset$status == 0, "obstime"],
      
      # predictor matrices (time-fixed)
      p = p,
      q = q,
      
      x_event = X_main[training_dataset$status == 1, ],
      x_int_event =  X_int[training_dataset$status == 1, ],
      
      x_rcens = X_main[training_dataset$status == 0, ],
      x_int_rcens = X_int[training_dataset$status == 0, ],
      
      # link the interaction effect with the corresponding main effects
      g1 = g1,
      g2 = g2
    )
    
    #----------------------------
    # Prepare data for log_lik
    #-----------------------------
    
    if (log_lik_saving == TRUE) {
      #----- for generate log_lik
      stan_data = c(
        stan_data,
        list(
          #--- for log_lik ----
          N = nrow(training_dataset),
          status = training_dataset$status,
          t = training_dataset$obstime,
          X = X_main,
          X_int = X_int
        )
      )
      
    }
    
    #----------------------------
    # Prepare data for prediction
    #-----------------------------
    if (!is.null(testing_dataset)) {
      nnew <- nrow(testing_dataset)
      t_new <- testing_dataset$obstime
      design__matrix_2 = testing_dataset[, !(names(testing_dataset) %in% c("id", "obstime", "status"))]
      
      column_names = colnames(design__matrix_2)
      main_names =  column_names[!grepl(":", column_names) &
                                   column_names != "(Intercept)"] #whether or not having intercept needs to be verified
      X_new_main = design__matrix_2[, main_names]
      
      int_names =  column_names[grepl(":", column_names)]
      X_new_int = design__matrix_2[, int_names]
      
      
      #----------------
      # add data to Constructed data for prediction
      #----------------
      
      stan_data = c(
        stan_data,
        list(
          nnew = nrow(testing_dataset),
          x_new = X_new_main,
          x_int_new = X_new_int,
          t_new = testing_dataset[, "obstime"]
        )
      )
    }
  } else {
    main_names =  column_names[!grepl(":", column_names)] #whether or not having intercept needs to be verified
    X_main = design__matrix[, main_names]
    
    # x_event = X_main[training_dataset$status == 1,]
    # x_int_event =  X_int[training_dataset$status == 1,]
    
    
    
    #----------------
    # Construct data
    #----------------
    stan_data = list(
      #----- for model fitting --------
      nevent = nrow(training_dataset[training_dataset$status == 1, ]),
      nrcens = nrow(training_dataset[training_dataset$status == 0, ]),
      t_event = training_dataset[training_dataset$status == 1, "obstime"],
      t_rcens = training_dataset[training_dataset$status == 0, "obstime"],
      
      # predictor matrices (time-fixed)
      
      x_event = X_main[training_dataset$status == 1, ],
      x_rcens = X_main[training_dataset$status == 0, ],
    )
    
    #----------------------------
    # Prepare data for log_lik
    #-----------------------------
    
    if (log_lik_saving == TRUE) {
      #----- for generate log_lik
      stan_data = c(
        stan_data,
        list(
          #--- for log_lik ----
          N = nrow(training_dataset),
          status = training_dataset$status,
          t = training_dataset$obstime,
          X = X_main
        )
      )
      
    }
    
    #----------------------------
    # Prepare data for prediction
    #-----------------------------
    if (!is.null(testing_dataset)) {
      nnew <- nrow(testing_dataset)
      t_new <- testing_dataset$obstime
      design__matrix_2 = testing_dataset[, !(names(testing_dataset) %in% c("id", "obstime", "status"))]
      
      column_names = colnames(design__matrix_2)
      main_names =  column_names[!grepl(":", column_names) &
                                   column_names != "(Intercept)"] #whether or not having intercept needs to be verified
      X_new_main = design__matrix_2[, main_names]
      
      
      #----------------
      # add data to Constructed data for prediction
      #----------------
      
      stan_data = c(stan_data,
                    list(
                      nnew = nrow(testing_dataset),
                      x_new = X_new_main,
                      t_new = testing_dataset[, "obstime"]
                    ))
    }
  }
  
  return(stan_data)
}


#' @param qnodes The number of nodes to use for the Gauss-Kronrod quadrature
#'   that is used to evaluate the cumulative hazard when \code{basehaz = "bs"}
#'   Options are 15 (the default), 11 or 7.

# assumptions: all data are right censored data with the observed window being [0, obs_window]
stan_data_Constructer_BSpline <- function (training_dataset, testing_dataset, obs_window, qnodes = 15, log_lik_saving = FALSE, havingInt = TRUE){
  
  #----------------------------
  # Prepare data for model fitting
  #-----------------------------
  
  #----- organize the data regarding predictor
  
  design__matrix = training_dataset[,!(names(training_dataset) %in% c("id", "obstime", "status"))]
  
  column_names = colnames(design__matrix)
  main_names =  column_names[!grepl(":", column_names)] #whether or not having intercept needs to be verified
  X_main = design__matrix[, main_names]
  
  if (havingInt == TRUE) {
    int_names =  column_names[grepl(":", column_names)]
    X_int = design__matrix[, int_names]
    
    p <- dim(X_main)[2]
    q <- dim(X_int)[2]
    
    main_indices_for_int = find_main_effect_indices(int_names, main_names)
    g1 <- g(main_indices_for_int, 1)
    g2 <- g(main_indices_for_int, 2)
    
    t_event = training_dataset[training_dataset$status == 1, "obstime"]
    t_rcens = training_dataset[training_dataset$status == 0, "obstime"]
    
    
    
    #----- baseline hazard
    basehaz <- handle_basehaz_surv(
      times          = training_dataset$obstime,
      status         = training_dataset$status,
      min_t          = 0,
      max_t          = max(c(training_dataset$obstime, obs_window), na.rm = TRUE)
    )
    nvars <- basehaz$nvars # number of basehaz aux parameters
    
    
    
    # model uses quadrature
    
    # standardised nodes and weights for quadrature
    qq <- get_quadpoints(nodes = qnodes)
    qp <- qq$points
    qw <- qq$weights
    
    # quadrature points, evaluated for each row of data
    qpts_event <- uapply(qp, unstandardise_qpts, 0, t_event)
    qpts_rcens <- uapply(qp, unstandardise_qpts, 0, t_rcens)
    
    
    # quadrature weights, evaluated for each row of data
    qwts_event <- uapply(qw, unstandardise_qwts, 0, t_event)
    qwts_rcens <- uapply(qw, unstandardise_qwts, 0, t_rcens)
    
    
    # times at events and all quadrature points
    cpts_list <- list(t_event, qpts_event, qpts_rcens)
    
    idx_cpts <- get_idx_array(sapply(cpts_list, length))
    # cpts     <- unlist(cpts_list) # as vector
    
    # number of quadrature points
    qevent <- length(qwts_event)
    qrcens <- length(qwts_rcens)
    
    
    #----- basis terms for baseline hazard
    
    basis_epts_event <- make_basis(t_event, basehaz)
    print("dim of basis_epts_event is:")
    print(dim(basis_epts_event))
    basis_qpts_event <- make_basis(qpts_event, basehaz)
    print("dim of qpts_event is:")
    print(dim(basis_qpts_event))
    basis_qpts_rcens <- make_basis(qpts_rcens, basehaz)
    print("dim of qpts_rcens is:")
    print(dim(basis_qpts_rcens))
    
    
    #----- model frames for generating predictor matrices
    
    id_event <- which(training_dataset$status == 1)
    id_rcens <-  which(training_dataset$status == 0)
    
    # combined model frame, with quadrature
    X_main_cpts <- rbind(X_main[id_event, ],
                         rep_rows(X_main[id_event, ], times = qnodes),
                         rep_rows(X_main[id_rcens, ], times = qnodes))
    X_int_cpts <- rbind(X_int[id_event, ],
                        rep_rows(X_int[id_event, ], times = qnodes),
                        rep_rows(X_int[id_rcens, ], times = qnodes))
    
    
    
    # time-fixed predictor matrices, with quadrature
    x_epts_event <- X_main_cpts[idx_cpts[1, 1]:idx_cpts[1, 2], , drop = FALSE]
    x_qpts_event <- X_main_cpts[idx_cpts[2, 1]:idx_cpts[2, 2], , drop = FALSE]
    x_qpts_rcens <- X_main_cpts[idx_cpts[3, 1]:idx_cpts[3, 2], , drop = FALSE]
    
    x_int_epts_event <- X_int_cpts[idx_cpts[1, 1]:idx_cpts[1, 2], , drop = FALSE]
    x_int_qpts_event <- X_int_cpts[idx_cpts[2, 1]:idx_cpts[2, 2], , drop = FALSE]
    x_int_qpts_rcens <- X_int_cpts[idx_cpts[3, 1]:idx_cpts[3, 2], , drop = FALSE]
    
    
    
    
    #----------------
    # Construct data
    #----------------
    stan_data = list(
      basis = basehaz$bs_basis,
      #----- for model fitting --------
      p = p,
      q = q,
      
      nvars = nvars,
      
      qnodes  = qnodes,
      
      
      Nevent       = sum(training_dataset$status == 1),
      Nrcens       = sum(training_dataset$status == 0),
      
      qevent = qevent,
      qrcens = qrcens,
      
      
      epts_event   = t_event,
      qpts_event = qpts_event,
      qpts_rcens = qpts_rcens,
      
      
      qwts_event = qwts_event,
      qwts_rcens = qwts_rcens,
      
      
      x_epts_event = x_epts_event,
      x_qpts_event = x_qpts_event,
      x_qpts_rcens = x_qpts_rcens,
      
      x_int_epts_event = x_int_epts_event,
      x_int_qpts_event = x_int_qpts_event,
      x_int_qpts_rcens = x_int_qpts_rcens,
      
      
      basis_epts_event = basis_epts_event,
      basis_qpts_event = basis_qpts_event,
      basis_qpts_rcens = basis_qpts_rcens,
      
      # link the interaction effect with the corresponding main effects
      g1 = g1,
      g2 = g2
    )
    
    #----------------------------
    # Prepare data for saving log_lik
    #-----------------------------
    if (log_lik_saving == TRUE) {

      t_all = training_dataset$obstime
      qpts_all <- uapply(qp, unstandardise_qpts, 0, t_all)
      
      qwts_all <- uapply(qw, unstandardise_qwts, 0, t_all)
      
      basis_epts_all <- make_basis(t_all, basehaz)
      basis_qpts_all <- make_basis(qpts_all, basehaz)
      
      
      
      # combined model frame, with quadrature
      X_main_cpts_all <- rbind(X_main, X_main[rep(1:nrow(X_main), each = qnodes), ])
      X_int_cpts_all <- rbind(X_int, X_int[rep(1:nrow(X_int), each = qnodes), ])
      
      
      
      # time-fixed predictor matrices, with quadrature
      Nobs <- nrow(training_dataset)
      x_epts_all <- X_main_cpts_all[1:Nobs, , drop = FALSE]
      x_qpts_all <- X_main_cpts_all[Nobs:(Nobs * (qnodes + 1) - 1), , drop = FALSE]
      
      x_int_epts_all <- X_int_cpts_all[1:Nobs, , drop = FALSE]
      x_int_qpts_all <- X_int_cpts_all[Nobs:(Nobs * (qnodes + 1) - 1), , drop = FALSE]
      
      
      stan_data = c(
        stan_data,
        list(
          Nobs = nrow(training_dataset),
          Nobs_qnode = nrow(training_dataset) * 15,
          status = training_dataset$status,
          x_epts_all = x_epts_all,
          x_int_epts_all = x_int_epts_all,
          x_qpts_all = x_qpts_all,
          x_int_qpts_all = x_int_qpts_all,
          
          basis_epts_all = basis_epts_all,
          basis_qpts_all = basis_qpts_all,
          qwts_all = qwts_all
        )
      )
    }
    
    
    #----------------------------
    # Prepare data for prediction
    #-----------------------------
    
    if (!is.null(testing_dataset)) {
      nnew <- nrow(testing_dataset)
      t_new <- testing_dataset$obstime
      design__matrix_2 = testing_dataset[, !(names(testing_dataset) %in% c("id", "obstime", "status"))]
      
      column_names = colnames(design__matrix_2)
      main_names =  column_names[!grepl(":", column_names) &
                                   column_names != "(Intercept)"] #whether or not having intercept needs to be verified
      X_new_main = design__matrix_2[, main_names]
      
      int_names =  column_names[grepl(":", column_names)]
      X_new_int = design__matrix_2[, int_names]
      
      
      #----------------------------
      # prediction the survival probability
      # for each participant at a specific
      # time point
      #-----------------------------
      
      # quadrature points, evaluated for each row of data
      qpts_new_event <- uapply(qp, unstandardise_qpts, 0, t_new) #500 first node for all t_new, 500 second node for all t_new
      
      
      # quadrature weights, evaluated for each row of data
      qwts_new_event <- uapply(qw, unstandardise_qwts, 0, t_new)
      
      
      
      # number of quadrature points
      qevent_new <- length(qwts_new_event)
      
      
      #----- basis terms for baseline hazard
      basis_new_qpts_event <- make_basis(qpts_new_event, basehaz) #500 first node for all t_new, 500 second node for all t_new
      
      
      
      #----- model frames for generating predictor matrices
      
      
      # combined model frame, with quadrature
      X_main_cpts_new <- rep_rows(X_new_main, times = qnodes) #X_main first for all t_new; X_main second for all t_new
      X_int_cpts_new <- rep_rows(X_new_int, times = qnodes)
      
      
      
      # time-fixed predictor matrices, with quadrature
      x_new_qpts_event <- X_main_cpts_new[, , drop = FALSE] #X_main first for all t_new; X_main second for all t_new
      
      x_new_int_qpts_event <- X_int_cpts_new[, , drop = FALSE]
      
      #----------------
      # add data to Constructed data for prediction
      #----------------
      
      stan_data = c(
        stan_data,
        list(
          #--- for prediction ----
          nnew = nrow(testing_dataset),
          qevent_new = qevent_new,
          t_new = testing_dataset$obstime,
          x_new_qpts_event = x_new_qpts_event,
          x_new_int_qpts_event = x_new_int_qpts_event,
          basis_qpts_event_new = basis_new_qpts_event,
          qwts_event_new = qwts_new_event,
          x_new = X_new_main,
          x_int_new = X_new_int
        )
      )
    }
  } else{
    p <- dim(X_main)[2]
    
    t_event = training_dataset[training_dataset$status == 1, "obstime"]
    t_rcens = training_dataset[training_dataset$status == 0, "obstime"]
    
    
    
    #----- baseline hazard
    basehaz <- handle_basehaz_surv(
      times          = training_dataset$obstime,
      status         = training_dataset$status,
      min_t          = 0,
      max_t          = max(c(training_dataset$obstime, obs_window), na.rm = TRUE)
    )
    nvars <- basehaz$nvars # number of basehaz aux parameters
    
    
    
    # model uses quadrature
    
    # standardised nodes and weights for quadrature
    qq <- get_quadpoints(nodes = qnodes)
    qp <- qq$points
    qw <- qq$weights
    
    # quadrature points, evaluated for each row of data
    qpts_event <- uapply(qp, unstandardise_qpts, 0, t_event)
    qpts_rcens <- uapply(qp, unstandardise_qpts, 0, t_rcens)
    
    
    # quadrature weights, evaluated for each row of data
    qwts_event <- uapply(qw, unstandardise_qwts, 0, t_event)
    qwts_rcens <- uapply(qw, unstandardise_qwts, 0, t_rcens)
    
    
    # times at events and all quadrature points
    cpts_list <- list(t_event, qpts_event, qpts_rcens)
    
    idx_cpts <- get_idx_array(sapply(cpts_list, length))
    # cpts     <- unlist(cpts_list) # as vector
    
    # number of quadrature points
    qevent <- length(qwts_event)
    qrcens <- length(qwts_rcens)
    
    
    #----- basis terms for baseline hazard
    
    basis_epts_event <- make_basis(t_event, basehaz)
    basis_qpts_event <- make_basis(qpts_event, basehaz)
    basis_qpts_rcens <- make_basis(qpts_rcens, basehaz)
    
    
    #----- model frames for generating predictor matrices
    
    id_event <- which(training_dataset$status == 1)
    id_rcens <-  which(training_dataset$status == 0)
    
    # combined model frame, with quadrature
    X_main_cpts <- rbind(X_main[id_event, ],
                         rep_rows(X_main[id_event, ], times = qnodes),
                         rep_rows(X_main[id_rcens, ], times = qnodes))
    
    
    
    # time-fixed predictor matrices, with quadrature
    x_epts_event <- X_main_cpts[idx_cpts[1, 1]:idx_cpts[1, 2], , drop = FALSE]
    x_qpts_event <- X_main_cpts[idx_cpts[2, 1]:idx_cpts[2, 2], , drop = FALSE]
    x_qpts_rcens <- X_main_cpts[idx_cpts[3, 1]:idx_cpts[3, 2], , drop = FALSE]
    
    
    
    
    #----------------
    # Construct data
    #----------------
    stan_data = list(
      p = p,
      basis = basehaz$bs_basis,
      #----- for model fitting --------
      nvars = nvars,
      
      qnodes  = qnodes,
      
      
      Nevent       = sum(training_dataset$status == 1),
      Nrcens       = sum(training_dataset$status == 0),
      
      qevent = qevent,
      qrcens = qrcens,
      
      
      epts_event   = t_event,
      qpts_event = qpts_event,
      qpts_rcens = qpts_rcens,
      
      
      qwts_event = qwts_event,
      qwts_rcens = qwts_rcens,
      
      
      x_epts_event = x_epts_event,
      x_qpts_event = x_qpts_event,
      x_qpts_rcens = x_qpts_rcens,
      
      
      basis_epts_event = basis_epts_event,
      basis_qpts_event = basis_qpts_event,
      basis_qpts_rcens = basis_qpts_rcens
    )
    
    #----------------------------
    # Prepare data for saving log_lik
    #-----------------------------
    if (log_lik_saving == TRUE) {
      
      t_all = training_dataset$obstime
      qpts_all <- uapply(qp, unstandardise_qpts, 0, t_all)
      
      qwts_all <- uapply(qw, unstandardise_qwts, 0, t_all)
      
      basis_epts_all <- make_basis(t_all, basehaz)
      basis_qpts_all <- make_basis(qpts_all, basehaz)
      
      
      
      # combined model frame, with quadrature
      X_main_cpts_all <- rbind(X_main, X_main[rep(1:nrow(X_main), each = qnodes), ])
      
      
      # time-fixed predictor matrices, with quadrature
      Nobs <- nrow(training_dataset)
      x_epts_all <- X_main_cpts_all[1:Nobs, , drop = FALSE]
      x_qpts_all <- X_main_cpts_all[Nobs:(Nobs * (qnodes + 1) - 1), , drop = FALSE]
      
      
      stan_data = c(
        stan_data,
        list(
          Nobs = nrow(training_dataset),
          Nobs_qnode = nrow(training_dataset) * 15,
          status = training_dataset$status,
          x_epts_all = x_epts_all,
          x_qpts_all = x_qpts_all,
          
          basis_epts_all = basis_epts_all,
          basis_qpts_all = basis_qpts_all,
          qwts_all = qwts_all
        )
      )
    }
    
    
    #----------------------------
    # Prepare data for prediction
    #-----------------------------
    
    if (!is.null(testing_dataset)) {
      nnew <- nrow(testing_dataset)
      t_new <- testing_dataset$obstime
      design__matrix_2 = testing_dataset[, !(names(testing_dataset) %in% c("id", "obstime", "status"))]
      
      column_names = colnames(design__matrix_2)
      main_names =  column_names[!grepl(":", column_names) &
                                   column_names != "(Intercept)"] #whether or not having intercept needs to be verified
      X_new_main = design__matrix_2[, main_names]
      
      
      #----------------------------
      # prediction the survival probability
      # for each participant at a specific
      # time point
      #-----------------------------
      
      # quadrature points, evaluated for each row of data
      qpts_new_event <- uapply(qp, unstandardise_qpts, 0, t_new) #500 first node for all t_new, 500 second node for all t_new
      
      
      # quadrature weights, evaluated for each row of data
      qwts_new_event <- uapply(qw, unstandardise_qwts, 0, t_new)
      
      
      
      # number of quadrature points
      qevent_new <- length(qwts_new_event)
      
      
      #----- basis terms for baseline hazard
      basis_new_qpts_event <- make_basis(qpts_new_event, basehaz) #500 first node for all t_new, 500 second node for all t_new
      
      
      
      #----- model frames for generating predictor matrices
      
      
      # combined model frame, with quadrature
      X_main_cpts_new <- rep_rows(X_new_main, times = qnodes) #X_main first for all t_new; X_main second for all t_new
      
      
      # time-fixed predictor matrices, with quadrature
      x_new_qpts_event <- X_main_cpts_new[, , drop = FALSE] #X_main first for all t_new; X_main second for all t_new
      
      #----------------
      # add data to Constructed data for prediction
      #----------------
      
      stan_data = c(
        stan_data,
        list(
          #--- for prediction ----
          nnew = nrow(testing_dataset),
          qevent_new = qevent_new,
          t_new = testing_dataset$obstime,
          x_new_qpts_event = x_new_qpts_event,
          basis_qpts_event_new = basis_new_qpts_event,
          qwts_event_new = qwts_new_event,
          x_new = X_new_main
          )
      )
    }
  }
  return(stan_data)
  
}  



#---------- internal
# Select rows of a matrix
#
# @param x A matrix.
# @param rows Logical or numeric vector stating which rows of 'x' to retain.
keep_rows <- function(x, rows = 1:nrow(x)) {
  x[rows, , drop = FALSE]
}  


# Replicate rows of a matrix or data frame
#
# @param x A matrix or data frame.
# @param ... Arguments passed to 'rep', namely 'each' or 'times'.
rep_rows <- function(x, ...) {
  if (is.null(x) || !nrow(x)) {
    return(x)
  } else if (is.matrix(x) || is.data.frame(x)) {
    x <- x[rep(1:nrow(x), ...), , drop = FALSE]
  } else {
    stop2("'x' must be a matrix or data frame.")
  }
  x
}

# From a vector of length M giving the number of elements (for example number
# of parameters or observations) for each submodel, create an indexing array 
# of dimension M * 2, where column 1 is the beginning index and 2 is the end index
#
# @param x A numeric vector
# @return A length(x) * 2 array
get_idx_array <- function(x) {
  as.array(do.call("rbind", lapply(1:length(x), function(i) {
    idx_beg <- ifelse(x[i] > 0L, sum(x[0:(i-1)]) + 1, 0L)
    idx_end <- ifelse(x[i] > 0L, sum(x[0:i]),         0L)
    c(idx_beg, idx_end)
  })))
}  


# Convert a standardised quadrature weight to an unstandardised value based on 
# the specified integral limits
#
# @param x An unstandardised quadrature weight
# @param a The lower limit(s) of the integral, possibly a vector
# @param b The upper limit(s) of the integral, possibly a vector
unstandardise_qwts <- function(x, a, b, na.ok = TRUE) {
  if (!identical(length(x), 1L) || !is.numeric(x))
    stop2("'x' should be a single numeric value.")
  if (!length(a) %in% c(1L, length(b)))
    stop2("'a' and 'b' should be vectors of length 1, or, be the same length.")
  if (!na.ok) {
    if (!all(is.numeric(a), is.numeric(b)))
      stop2("'a' and 'b' should be numeric.")
    if (any((b - a) < 0))
      stop2("The upper limits for the integral ('b' values) should be greater than ",
            "the corresponding lower limits for the integral ('a' values).")
  }
  ((b - a) / 2) * x
}


# Convert a standardised quadrature node to an unstandardised value based on 
# the specified integral limits
#
# @param x An unstandardised quadrature node
# @param a The lower limit(s) of the integral, possibly a vector
# @param b The upper limit(s) of the integral, possibly a vector
unstandardise_qpts <- function(x, a, b, na.ok = TRUE) {
  if (!identical(length(x), 1L) || !is.numeric(x))
    stop2("'x' should be a single numeric value.")
  if (!length(a) %in% c(1L, length(b)))
    stop2("'a' and 'b' should be vectors of length 1, or, be the same length.")
  if (!na.ok) {
    if (!all(is.numeric(a), is.numeric(b)))
      stop2("'a' and 'b' should be numeric.")
    if (any((b - a) < 0))
      stop2("The upper limits for the integral ('b' values) should be greater than ",
            "the corresponding lower limits for the integral ('a' values).")
  }
  ((b - a) / 2) * x + ((b + a) / 2)
}



# Return the cutpoints for a specified number of quantiles of 'x'
#
# @param x A numeric vector.
# @param nq Integer specifying the number of quantiles.
# @return A vector of percentiles corresponding to percentages 100*k/m for 
#   k=1,2,...,nq-1.
qtile <- function(x, nq = 2) {
  if (nq > 1) {
    probs <- seq(1, nq - 1) / nq
    return(quantile(x, probs = probs))
  } else if (nq == 1) {
    return(NULL)
  } else {
    stop("'nq' must be >= 1.")
  }
}




# Return a vector with internal knots for 'x', based on evenly spaced quantiles
#
# @param x A numeric vector.
# @param df The degrees of freedom. If specified, then 'df - degree - intercept'.
#   knots are placed at evenly spaced percentiles of 'x'. If 'iknots' is 
#   specified then 'df' is ignored.
# @param degree Non-negative integer. The degree for the spline basis.
# @param iknots Optional vector of internal knots.
# @return A numeric vector of internal knot locations, or NULL if there are
#   no internal knots.
get_iknots <- function(x, df = 5L, degree = 3L, iknots = NULL, intercept = FALSE) {
  
  # obtain number of internal knots
  if (is.null(iknots)) {
    nk <- df - degree - intercept
  } else {
    nk <- length(iknots)
  }
  
  # validate number of internal knots
  if (nk < 0) {
    stop2("Number of internal knots cannot be negative.")
  }
  
  # if no internal knots then return empty vector
  if (nk == 0) {
    return(numeric(0))
  }
  
  # obtain default knot locations if necessary
  if (is.null(iknots)) {
    iknots <- qtile(x, nq = nk + 1) # evenly spaced percentiles
  }
  
  # return internal knot locations, ensuring they are positive
  validate_positive_scalar(iknots)
  
  return(iknots)
}


# Throw error if parameter isn't a positive scalar
#
# @param x The object to test.
validate_positive_scalar <- function(x, not_greater_than = NULL) {
  nm <- deparse(substitute(x))
  if (is.null(x))
    stop(nm, " cannot be NULL", call. = FALSE)
  if (!is.numeric(x))
    stop(nm, " should be numeric", call. = FALSE)
  if (any(x <= 0)) 
    stop(nm, " should be postive", call. = FALSE)
  if (!is.null(not_greater_than)) {
    if (!is.numeric(not_greater_than) || (not_greater_than <= 0))
      stop("'not_greater_than' should be numeric and postive")
    if (!all(x <= not_greater_than))
      stop(nm, " should less than or equal to ", not_greater_than, call. = FALSE)
  }
}



# Function to return standardised GK quadrature points and weights
#
# @param nodes The required number of quadrature nodes
# @return A list with two named elements (points and weights) each
#   of which is a numeric vector with length equal to the number of
#   quadrature nodes
get_quadpoints <- function(nodes = 15) {
  if (!is.numeric(nodes) || (length(nodes) > 1L)) {
    stop("'qnodes' should be a numeric vector of length 1.")
  } else if (nodes == 21) {
    list(
      points = c(
        -0.9956571630258081,
        -0.9739065285171717,
        -0.9301574913557082,
        -0.8650633666889845,
        -0.7808177265864169,
        -0.6794095682990244,
        -0.5627571346686047,
        -0.4333953941292472,
        -0.2943928627014602,
        -0.1488743389816312,
        0,
        0.1488743389816312,
        0.2943928627014602,
        0.4333953941292472,
        0.5627571346686047,
        0.6794095682990244,
        0.7808177265864169,
        0.8650633666889845,
        0.9301574913557082,
        0.9739065285171717,
        0.9956571630258081
      ),
      weights = c(
        0.0116946388673719,
        0.0325581623079647,
        0.0547558965743519,
        0.0750396748109199,
        0.0931254545836976,
        0.1093871588022976,
        0.1234919762620659,
        0.1347092173114733,
        0.1427759385770601,
        0.1477391049013385,
        0.1494455540029169,
        0.1477391049013385,
        0.1427759385770601,
        0.1347092173114733,
        0.1234919762620659,
        0.1093871588022976,
        0.0931254545836976,
        0.0750396748109199,
        0.0547558965743519,
        0.0325581623079647,
        0.0116946388673719
      )
    )
    
  } else if (nodes == 15) {
    list(
      points = c(
        -0.991455371120812639207,
        -0.949107912342758524526,
        -0.86486442335976907279,
        -0.7415311855993944398639,
        -0.5860872354676911302941,
        -0.4058451513773971669066,
        -0.2077849550078984676007,
        0,
        0.2077849550078984676007,
        0.405845151377397166907,
        0.5860872354676911302941,
        0.741531185599394439864,
        0.86486442335976907279,
        0.9491079123427585245262,
        0.991455371120812639207),
      weights = c(
        0.0229353220105292249637,
        0.063092092629978553291,
        0.10479001032225018384,
        0.140653259715525918745,
        0.1690047266392679028266,
        0.1903505780647854099133,
        0.204432940075298892414,
        0.209482141084727828013,
        0.204432940075298892414,
        0.1903505780647854099133,
        0.169004726639267902827,
        0.140653259715525918745,
        0.1047900103222501838399,
        0.063092092629978553291,
        0.0229353220105292249637))      
  } else if (nodes == 11) {
    list(
      points = c(
        -0.984085360094842464496,
        -0.906179845938663992798,
        -0.754166726570849220441,
        -0.5384693101056830910363,
        -0.2796304131617831934135,
        0,
        0.2796304131617831934135,
        0.5384693101056830910363,
        0.754166726570849220441,
        0.906179845938663992798,
        0.984085360094842464496),
      weights = c(
        0.042582036751081832865,
        0.1152333166224733940246,
        0.186800796556492657468,
        0.2410403392286475866999,
        0.272849801912558922341,
        0.2829874178574912132043,
        0.272849801912558922341,
        0.241040339228647586701,
        0.186800796556492657467,
        0.115233316622473394025,
        0.042582036751081832865))     
  } else if (nodes == 7) {
    list(
      points = c(
        -0.9604912687080202834235,
        -0.7745966692414833770359,
        -0.4342437493468025580021,
        0,
        0.4342437493468025580021,
        0.7745966692414833770359,
        0.9604912687080202834235),
      weights = c(
        0.1046562260264672651938,
        0.268488089868333440729,
        0.401397414775962222905,
        0.450916538658474142345,
        0.401397414775962222905,
        0.268488089868333440729,
        0.104656226026467265194))      
  } else stop("'qnodes' must be either 7, 11 or 15.")  
}


# Unlist the result from an lapply call
#
# @param X,FUN,... Same as lapply
uapply <- function(X, FUN, ...) {
  unlist(lapply(X, FUN, ...))
}

# Construct a list with information about the baseline hazard
#

# @param times A numeric vector with eventtimes for each individual
# @param status A numeric vector with event indicators for each individual
# @param min_t Scalar, the minimum entry time across all individuals
# @param max_t Scalar, the maximum event or censoring time across all individuals
# @return A named list with the following elements:

#   user_df: integer specifying the input to the df argument
#   df: integer specifying the number of parameters to use for the 
#     baseline hazard.
#   knots: the knot locations for the baseline hazard.
#   bs_basis: The basis terms for the B-splines. This is passed to Stan
#     as the "model matrix" for the baseline hazard. It is also used in
#     post-estimation when evaluating the baseline hazard for posterior
#     predictions since it contains information about the knot locations
#     for the baseline hazard (this is implemented via splines::predict.bs). 
handle_basehaz_surv <- function(times, 
                                status,
                                min_t, max_t) {
  
  df     <- 5
  degree <- 3
  
  
  # tt <- times#sort(unique(times))#[status == 1]) # uncensored event times
  
  bknots <- c(min_t, max_t)
  # iknots <- get_iknots(tt, df = df, degree = degree)
  # basis  <- splines2::bSpline(tt, iknots = iknots, Boundary.knots = bknots, degree = 3, intercept = FALSE)
  iknots <- get_iknots(sort(times[status == 1]), df = df, degree = degree)
  basis  <- splines2::bSpline(sort(unique(times)), iknots = iknots, Boundary.knots = bknots, degree = 3, intercept = FALSE)
  
  
  nvars  <- ncol(basis)  # number of aux parameters, basis terms
  
  
  rstan::nlist(
    nvars, 
    iknots, 
    bknots,
    degree,
    basis,
    df = nvars,
    user_df = nvars,
    knots = iknots,
    bs_basis = basis)
}

# Return a vector with valid names for elements in the list passed to the
# 'basehaz_ops' argument of a 'stan_jm' or 'stan_surv' call
#
# @param basehaz_name A character string, the type of baseline hazard.
# @return A character vector, or NA if unmatched.
get_ok_basehaz_ops <- function(basehaz_name) {
  switch(basehaz_name,
         "bs"        = c("df", "knots", "degree"),
         NA)
}

# Return the integer representation for the baseline hazard, used by Stan
#
# @param basehaz_name A character string, the type of baseline hazard.
# @return An integer, or NA if unmatched.
basehaz_for_stan <- function(basehaz_name) {
  switch(basehaz_name, 
         "bs"          = 2L,
         NA)
}

# Return a vector with internal knots for 'x', based on evenly spaced quantiles
#
# @param x A numeric vector.
# @param df The degrees of freedom. If specified, then 'df - degree - intercept'.
#   knots are placed at evenly spaced percentiles of 'x'. If 'iknots' is 
#   specified then 'df' is ignored.
# @param degree Non-negative integer. The degree for the spline basis.
# @param iknots Optional vector of internal knots.
# @return A numeric vector of internal knot locations, or NULL if there are
#   no internal knots.
get_iknots <- function(x, df = 5L, degree = 3L, iknots = NULL, intercept = FALSE) {
  
  # obtain number of internal knots
  if (is.null(iknots)) {
    nk <- df - degree - intercept
  } else {
    nk <- length(iknots)
  }
  
  # validate number of internal knots
  if (nk < 0) {
    stop2("Number of internal knots cannot be negative.")
  }
  
  # if no internal knots then return empty vector
  if (nk == 0) {
    return(numeric(0))
  }
  
  # obtain default knot locations if necessary
  if (is.null(iknots)) {
    iknots <- qtile(x, nq = nk + 1) # evenly spaced percentiles
  }
  
  # return internal knot locations, ensuring they are positive
  validate_positive_scalar(iknots)
  
  return(iknots)
}

# Identify whether the type of baseline hazard requires an intercept in
# the linear predictor (NB splines incorporate the intercept into the basis).
#
# @param basehaz A list with info about the baseline hazard; see 'handle_basehaz'.
# @return A Logical.
has_intercept <- function(basehaz) {
  nm <- get_basehaz_name(basehaz)
  (nm %in% c("bs"))
}



# Return the spline basis for the given type of baseline hazard.
# 
# @param times A numeric vector of times at which to evaluate the basis.
# @param basehaz A list with info about the baseline hazard, returned by a 
#   call to 'handle_basehaz'.
# @param integrate A logical, specifying whether to calculate the integral of
#   the specified basis.
# @return A matrix.
make_basis <- function(times, basehaz, integrate = FALSE) {
  N <- length(times)
  K <- basehaz$nvars
  if (!N) { # times is NULL or empty vector
    return(matrix(0, 0, K))
  } 
  basis_matrix(times, basis = basehaz$basis)
}


aa <- function(x, ...) as.array  (x, ...)

# Evaluate a spline basis matrix at the specified times
#
# @param time A numeric vector.
# @param basis Info on the spline basis.
# @param integrate A logical, should the integral of the basis be returned?
# @return A two-dimensional array.
basis_matrix <- function(times, basis, integrate = FALSE) {
  out <- predict(basis, times)
  if (integrate) {
    class(basis) <- c("splines2")
    out <- predict(basis, times)
  }
  aa(out)
}

