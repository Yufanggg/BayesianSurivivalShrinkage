# Part of the rstanarm package for estimating model parameters

#' @param basehaz A character string indicating which baseline hazard or
#'   baseline survival distribution to use for the event submodel. 
#'   
#'   The following are available under a hazard scale formulation: 
#'   \code{"bs"}: A flexible parametric model using cubic B-splines to 
#'     model the \emph{log} baseline hazard. The default locations for the  
#'     internal knots, as well as the basis terms for the splines, are calculated 
#'     with respect to time. A closed form solution for the cumulative hazard 
#'     is \strong{not} available regardless of whether or not the model includes
#'     time-varying effects; instead, quadrature is used to evaluate 
#'     the cumulative hazard at each MCMC iteration. Therefore, if your model
#'     does not include any time-varying effects, then estimation using the 
#'   
#' @param basehaz_ops A named list specifying options related to the baseline
#'   hazard. Currently this can include: \cr
#'   \itemize{
#'     \item \code{degree}: A positive integer specifying the degree for the 
#'     M-splines or B-splines. The default is \code{degree = 3}, which
#'     corresponds to cubic splines. Note that specifying \code{degree = 0}
#'     is also allowed and corresponds to piecewise constant.
#'     \item \code{df}: A positive integer specifying the degrees of freedom 
#'     for the M-splines or B-splines. For M-splines (i.e. when 
#'     \code{basehaz = "ms"}), two boundary knots and \code{df - degree - 1} 
#'     internal knots are used to generate the spline basis. For B-splines 
#'     (i.e. when \code{basehaz = "bs"}), two boundary knots and 
#'     \code{df - degree} internal knots are used to generate the spline 
#'     basis. The difference is due to the fact that the M-spline basis
#'     includes an intercept, whereas the B-spline basis does not. The 
#'     default is \code{df = 6} for M-splines and \code{df = 5} for
#'     B-splines (i.e. two boundary knots and two internal knots when the
#'     default cubic splines are being used). The internal knots are placed 
#'     at equally spaced percentiles of the distribution of uncensored event 
#'     times.
#'     \item \code{knots}: A numeric vector explicitly specifying internal 
#'     knot locations for the M-splines or B-splines. Note that \code{knots} 
#'     cannot be specified if \code{df} is specified.
#'   }
#'   Note that for the M-splines and B-splines -- in addition to any internal
#'   \code{knots} -- a lower boundary knot is placed at the earliest entry time
#'   and an upper boundary knot is placed at the latest event or censoring time.
#'   These boundary knot locations are the default and cannot be changed by the
#'   user.
#'   
#' @param qnodes The number of nodes to use for the Gauss-Kronrod quadrature
#'   that is used to evaluate the cumulative hazard when \code{basehaz = "bs"}
#'   Options are 15 (the default), 11 or 7.
 
# assumptions: all data are right censored data with the observed window being [0, obs_window]
stan_bSpline_data_Constructer(dataset, obs_window, ){
  
  basehaz = "bs"
  qnodes = 15
}  
  
  #----------------
  # Construct data
  #----------------
  
  
  #----- dimensions and response vectors
  
  status = dataset$status # event indicator for each row of data
  
  
  
  # time variables for stan
  t_event # exact event time
  t_rcens # right censoring time
  
  
  # dimensions
  nevent <- sum(status == 1)
  nrcens <- sum(status == 0)
  

  
  #----- baseline hazard
  basehaz <- handle_basehaz_surv(basehaz        = "bs", 
                                 ok_basehaz     = "bs",
                                 times          = dataset$obstime, 
                                 status         = status,
                                 min_t          = 0,
                                 max_t          = max(c(dataset$obstime, obs_window), na.rm = TRUE)){
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
    cpts_list <- list(t_event,
                      qpts_event,
                      qpts_rcens)
    
    idx_cpts <- get_idx_array(sapply(cpts_list, length))
    cpts     <- unlist(cpts_list) # as vector 
    
    # number of quadrature points
    qevent <- length(qwts_event)
    qlcens <- length(qwts_lcens)
    qrcens <- length(qwts_rcens)
    qicens <- length(qwts_icenl)
    qdelay <- length(qwts_delay)
    
 
  #----- basis terms for baseline hazard
  
    basis_epts_event <- make_basis(t_event,    basehaz)
    basis_qpts_event <- make_basis(qpts_event, basehaz)
    basis_qpts_lcens <- make_basis(qpts_lcens, basehaz)
    basis_qpts_rcens <- make_basis(qpts_rcens, basehaz)
    basis_qpts_icenl <- make_basis(qpts_icenl, basehaz)
    basis_qpts_icenu <- make_basis(qpts_icenu, basehaz)
    basis_qpts_delay <- make_basis(qpts_delay, basehaz)
    

  
  #----- model frames for generating predictor matrices
  
  mf_event <- keep_rows(mf, status == 1)
  mf_rcens <- keep_rows(mf, status == 0)
 
  

    
    # combined model frame, with quadrature
    mf_cpts <- rbind(mf_event,
                     rep_rows(mf_event, times = qnodes),
                     rep_rows(mf_rcens, times = qnodes))
    

  
  #----- time-fixed predictor matrices
  ff        <- formula$fe_form
  x         <- make_x(ff, mf     )$x
  x_cpts    <- make_x(ff, mf_cpts)$x
  x_centred <- sweep(x_cpts, 2, colMeans(x), FUN = "-")
  K         <- ncol(x_cpts)
  
 
    
    # time-fixed predictor matrices, with quadrature
    # NB skip index 6 on purpose, since time fixed predictor matrix is 
    # identical for lower and upper limits of interval censoring time
    x_epts_event <- x_centred[idx_cpts[1,1]:idx_cpts[1,2], , drop = FALSE]
    x_qpts_event <- x_centred[idx_cpts[2,1]:idx_cpts[2,2], , drop = FALSE]
    x_qpts_rcens <- x_centred[idx_cpts[4,1]:idx_cpts[4,2], , drop = FALSE]
    
  
  
  #----- stan data
  
  standata <- nlist(
    K, 
    nvars,


    qnodes       = qnodes,
    
    Nevent       = nevent,
    Nrcens       = nrcens, 

    
    qevent       = qevent,
    qrcens       = qrcens,

    
    epts_event   = t_event,
    qpts_event   = qpts_event,
    qpts_rcens   = qpts_rcens,

    
    qwts_event   = qwts_event,
    qwts_rcens   =qwts_rcens,

    
    x_epts_event = x_epts_event,
    x_qpts_event = x_qpts_event,
    x_qpts_rcens = ix_qpts_rcens,

    
    basis_epts_event = basis_epts_event,
    basis_qpts_event = basis_qpts_event,
    basis_qpts_rcens = basis_qpts_rcens,
  )
}



#---------- internal
  
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
  
  
# Function to return standardised GK quadrature points and weights
#
# @param nodes The required number of quadrature nodes
# @return A list with two named elements (points and weights) each
#   of which is a numeric vector with length equal to the number of
#   quadrature nodes
get_quadpoints <- function(nodes = 15) {
    if (!is.numeric(nodes) || (length(nodes) > 1L)) {
      stop("'qnodes' should be a numeric vector of length 1.")
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
  
  

# Construct a list with information about the baseline hazard
#
# @param basehaz A string specifying the type of baseline hazard
# @param basehaz_ops A named list with elements: df, knots, degree
# @param ok_basehaz A list of admissible baseline hazards
# @param times A numeric vector with eventtimes for each individual
# @param status A numeric vector with event indicators for each individual
# @param min_t Scalar, the minimum entry time across all individuals
# @param max_t Scalar, the maximum event or censoring time across all individuals
# @return A named list with the following elements:
#   type: integer specifying the type of baseline hazard, 1L = weibull,
#     2L = b-splines, 3L = piecewise.
#   type_name: character string specifying the type of baseline hazard.
#   user_df: integer specifying the input to the df argument
#   df: integer specifying the number of parameters to use for the 
#     baseline hazard.
#   knots: the knot locations for the baseline hazard.
#   bs_basis: The basis terms for the B-splines. This is passed to Stan
#     as the "model matrix" for the baseline hazard. It is also used in
#     post-estimation when evaluating the baseline hazard for posterior
#     predictions since it contains information about the knot locations
#     for the baseline hazard (this is implemented via splines::predict.bs). 
handle_basehaz_surv <- function(basehaz, 
                                basehaz_ops,
                                ok_basehaz,
                                times, 
                                status,
                                min_t, max_t) {
  
  if (!basehaz %in% ok_basehaz)
    stop2("'basehaz' should be one of: ", comma(ok_basehaz))
  
  ok_basehaz_ops <- get_ok_basehaz_ops(basehaz)
  if (!all(names(basehaz_ops) %in% ok_basehaz_ops))
    stop2("'basehaz_ops' can only include: ", comma(ok_basehaz_ops))
  
  if (basehaz %in% c("ms", "bs", "piecewise")) {
    
    df     <- basehaz_ops$df
    knots  <- basehaz_ops$knots
    degree <- basehaz_ops$degree
    
    if (!is.null(df) && !is.null(knots))
      stop2("Cannot specify both 'df' and 'knots' for the baseline hazard.")
    
    if (is.null(df))
      df <- switch(basehaz,
                   "bs"        = 5L, # assumes no intercept
                   df) # NB this is ignored if the user specified knots
    
    if (is.null(degree))
      degree <- 3L # cubic splines
    
    tt <- times[status == 1] # uncensored event times
    if (is.null(knots) && !length(tt)) {
      warning2("No observed events found in the data. Censoring times will ",
               "be used to evaluate default knot locations for splines.")
      tt <- times
    }
    
    if (!is.null(knots)) {
      if (any(knots < min_t))
        stop2("'knots' cannot be placed before the earliest entry time.")
      if (any(knots > max_t))
        stop2("'knots' cannot be placed beyond the latest event time.")
    }
    
  }
  
  if (basehaz == "bs") {
    
    bknots <- c(min_t, max_t)
    iknots <- get_iknots(tt, df = df, iknots = knots, degree = degree)
    basis  <- get_basis(tt, iknots = iknots, bknots = bknots, degree = degree, type = "bs")      
    nvars  <- ncol(basis)  # number of aux parameters, basis terms
    
  } 
  
  nlist(type_name = basehaz, 
        type = basehaz_for_stan(basehaz),
        nvars, 
        iknots, 
        bknots,
        degree,
        basis,
        df = nvars,
        user_df = nvars,
        knots = if (basehaz == "bs") iknots else c(bknots[1], iknots, bknots[2]),
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
  switch(basehaz$type_name,
         "bs"          = basis_matrix(times, basis = basehaz$basis),
         stop2("Bug found: unknown type of baseline hazard."))
}

# Evaluate a spline basis matrix at the specified times
#
# @param time A numeric vector.
# @param basis Info on the spline basis.
# @param integrate A logical, should the integral of the basis be returned?
# @return A two-dimensional array.
basis_matrix <- function(times, basis, integrate = FALSE) {
  out <- predict(basis, times)
  if (integrate) {
    stopifnot(inherits(basis, "MSpline"))
    class(basis) <- c("ISpline", "splines2", "matrix")
    out <- predict(basis, times)
  }
  aa(out)
}
