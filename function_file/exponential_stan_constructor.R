
# assumptions: all data are right censored data with the observed window being [0, obs_window]
stan_exponential_data_Constructer <- function (training_dataset, testing_dataset){
  
  #----------------------------
  # Prepare data for model fitting
  #-----------------------------
  
  #----- organize the data regarding predictor
  
  design__matrix = training_dataset[,!(names(training_dataset) %in% c("id", "obstime", "status"))]
  
  column_names = colnames(design__matrix)
  main_names =  column_names[!grepl(":", column_names)]
  X_main = design__matrix[, main_names]
  
  int_names =  column_names[grepl(":", column_names)]
  X_int = design__matrix[, int_names]
  
  p <- dim(X_main)[2]
  q <- dim(X_int)[2]
  
  main_indices_for_int = find_main_effect_indices(int_names, main_names)

  g1 <- g(main_indices_for_int, 1)
  g2 <- g(main_indices_for_int, 2)

  
  
  #----------------------------
  # Prepare data for prediction
  #-----------------------------
  
    
  nnew <- nrow(testing_dataset)
  design__matrix_2 = testing_dataset[,!(names(testing_dataset) %in% c("id", "obstime", "status"))]
  
  column_names = colnames(design__matrix_2)
  main_names =  column_names[!grepl(":", column_names)] #whether or not having intercept needs to be verified
  X_new_main = design__matrix_2[, main_names]
  
  int_names =  column_names[grepl(":", column_names)]
  X_new_int = design__matrix_2[, int_names]

  
 
  
  
  # Prepares data and parameter input for Stan.
  #----------------
  # Construct data
  #----------------
  stan_data <- list(
    #----- for model fitting --------
    nevent = nrow(training_dataset[training_dataset$status == 1, ]),
    nrcens = nrow(training_dataset[training_dataset$status == 0, ]),
    t_event = training_dataset[training_dataset$status == 1, "obstime"],
    t_rcens = training_dataset[training_dataset$status == 0, "obstime"],
    
    # predictor matrices (time-fixed)
    p = p,
    q = q,
    
    x_event = X_main[training_dataset$status == 1,],
    x_int_event =  X_int[training_dataset$status == 1,],
    
    x_rcens = X_main[training_dataset$status == 0,],
    x_int_rcens = X_int[training_dataset$status == 0,],
    
    # link the interaction effect with the corresponding main effects
    g1 = g1,
    g2 = g2,
    
    #----- for model prediction --------
    nnew = nrow(testing_dataset),
    x_new = X_new_main,
    x_int_new = X_new_int,
    t_new = testing_dataset[, "obstime"]
  )

  return(stan_data)
  
}  

