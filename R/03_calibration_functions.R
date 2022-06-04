#' Generate model outputs for calibration from a parameter set
#'
#' \code{calibration_out} computes model outputs to be used for calibration 
#' routines.
#'
#' @param v_params_calib Vector of parameters that need to be calibrated.
#' @param l_params_all List with all parameters of the decision model.
#' @return 
#' A list with disease-free survival (DFS), disease-specific survival (DSS), and
#' overall survival (OS)
#' @export
calibration_out <- function(v_params_calib, l_params_all){
  # Substitute values of calibrated parameters in base-case with 
  # calibrated values
  l_params_all <- update_param_list(l_params_all = l_params_all, 
                                    params_updated = v_params_calib)

  # Run the Decision Model
  # When all the initial members of the cohort are BM-positive
  l_out_stm <- run_drinkmc(l_params_all = l_params_all, 
                                    Trt = FALSE)

  # Epidemiological Output ----------------------------------------------------
  # TODO: Summary of MC outcome (calculate prevalence)
  
  # output list
  l_out <- list(v_dfs_BMpos = v_dfs_BMpos[surv_time], 
                v_dfs_BMneg = v_dfs_BMneg[surv_time],
                v_os_BMpos  = v_os_BMpos[surv_time],
                v_os_BMneg  = v_os_BMneg[surv_time],
                v_dss_BMpos = v_dss_BMpos[surv_time],
                v_dss_BMneg = v_dss_BMneg[surv_time]) #TODO: Update this
  return(l_out)
}

#' Sample from prior distributions of calibrated parameters
#'
#' \code{sample.prior} generates a sample of parameter sets from their prior 
#' distribution.
#' @param n_samp Number of samples.
#' @param v_param_names Vector with parameter names.
#' @param v_ub Vector with lower bounds for each parameter.
#' @param v_lb Vector with upper bounds for each parameter.
#' @return 
#' A matrix with 3 rows and \code{n_samp} rows. Each row corresponds to a 
#' parameter set sampled from their prior distributions
#' @examples 
#' v_param_names  <- c("r_Mets_DC",
#'                     "r_BMneg_Mets",
#'                     "hr_BMpos_Mets")
#' v_lb <- c(r_Mets_DC     = 0.037,               # lower bound
#'           r_BMneg_Mets  = 0.001,
#'           hr_BMpos_Mets = 1.58))  
#' v_ub <- c(r_Mets_DC     = -log(1-(1-0.03))/60, # upper bound 
#'           r_BMneg_Mets  = 0.03,
#'           hr_BMpos_Mets = 4.72) 
#' log_prior(v_params = sample.prior(n_samp = 5))
#' sample.prior(2)
#' @export
sample.prior <- function(n_samp,
                         v_param_names = c("r_Mets_DC", 
                                           "r_BMneg_Mets", 
                                           "hr_BMpos_Mets"),
                         v_lb = c(r_Mets_DC     = (0.037)*12, # O'Connell 2004 JNCI Stg IV Fig1 & Fig2;
                                  r_BMneg_Mets  = (0.001)*12, 
                                  hr_BMpos_Mets = 1.58),
                         v_ub = c(r_Mets_DC     = (-log(1-(1-0.03))/60)*12, # Rutter 2013 JNCI Table 4 5yr RS Colon cancer Stage IV 80+ Lower bound
                                  r_BMneg_Mets  = (0.03)*12,  
                                  hr_BMpos_Mets = 4.72)){
  # Number of parameters
  n_param <- length(v_param_names)
  
  ### Transformed design 
  ## Transformed bounds
  v_lb_transf <- c(log(v_lb))
  v_ub_transf <- c(log(v_ub))
  # Find means and SDs of log-normal based on bounds 
  # assuming bounds are represent the 95% equal tailed interval for these 
  # distributions
  v_mu_transf <- (v_ub_transf + v_lb_transf)/2
  v_sd_transf <- (v_ub_transf - v_lb_transf)/(2*2) # *1.96
  
  ### Draw LHS from Uniform[0,1] distributions
  m_lhs_unif   <- lhs::randomLHS(n = n_samp, k = n_param)
  colnames(m_lhs_unif) <- v_param_names
  
  ### Transformed LHS
  ## Get values in Normal scale
  m_lhs_normal <- m_lhs_unif
  for (i in 1:n_param){
    m_lhs_normal[, i] <- qnorm(m_lhs_unif[,i], v_mu_transf[i], v_sd_transf[i])
  }
  
  m_param_samp <- m_lhs_normal
  colnames(m_param_samp) <- v_param_names
  ## Get values in Original scale
  m_param_samp[, 1:n_param] <- exp(m_lhs_normal[, 1:n_param])
  
  return(m_param_samp)
}

#' Evaluate log-prior of calibrated parameters
#'
#' \code{log_prior} computes a log-prior value for one (or multiple) parameter 
#' set(s) based on their prior distributions.
#' @param v_params Vector (or matrix) of model parameters.
#' @param v_param_names Vector with parameter names.
#' @param v_ub Vector with lower bounds for each parameter.
#' @param v_lb Vector with upper bounds for each parameter.
#' @return 
#' A scalar (or vector) with log-prior values.
#' @examples
#' v_param_names  <- c("r_Mets_DC",
#'                     "r_BMneg_Mets",
#'                     "hr_BMpos_Mets")
#' v_lb <- c(r_Mets_DC     = 0.037,               # lower bound
#'           r_BMneg_Mets  = 0.001,
#'           hr_BMpos_Mets = 1.58))  
#' v_ub <- c(r_Mets_DC     = -log(1-(1-0.03))/60, # upper bound 
#'           r_BMneg_Mets  = 0.03,
#'           hr_BMpos_Mets = 4.72) 
#' log_prior(v_params = sample.prior(n_samp = 5))
#' @export
log_prior <- function(v_params, 
                      v_param_names = c("r_Mets_DC", 
                                        "r_BMneg_Mets", 
                                        "hr_BMpos_Mets"),
                      v_lb = c(r_Mets_DC     = (0.037)*12, 
                               r_BMneg_Mets  = (0.001)*12, 
                               hr_BMpos_Mets = 1.58),
                      v_ub = c(r_Mets_DC     = (-log(1-(1-0.03))/60)*12, 
                               r_BMneg_Mets  = (0.03)*12, 
                               hr_BMpos_Mets = 4.72)){
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  colnames(v_params) <- v_param_names
  ## Number of parameters
  n_param <- length(v_param_names)
  ## Number of samples
  n_samp <- nrow(v_params)
  
  ### Transformed design 
  ## Transformed bounds
  v_lb_transf <- c(log(v_lb[1:n_param]))
  v_ub_transf <- c(log(v_ub[1:n_param]))
  # Find means and SDs normal based on bounds assuming bounds are represent 
  # the 95% equal tailed interval for these distributions
  v_mu_transf <- (v_ub_transf + v_lb_transf)/2
  v_sd_transf <- (v_ub_transf - v_lb_transf)/(2*2) # *1.96
  
  lprior <- rep(0, n_samp)
  
  for (i in seq(n_param)) { # i <- 1
    lprior <- lprior + dlnorm(v_params[, i], v_mu_transf[i], v_sd_transf[i], log = T)
  }
  # lprior <- lprior + dlnorm(v_params[, 1], v_mu_transf[1], v_sd_transf[1], log = T)
  # lprior <- lprior + dlnorm(v_params[, 2], v_mu_transf[2], v_sd_transf[2], log = T)
  # lprior <- lprior + dlnorm(v_params[, 3], v_mu_transf[3], v_sd_transf[3], log = T)
  # lprior <- lprior + logti(v_params[, 3], v_mu_transf[3], v_sd_transf[3], log = T)
  # lprior 
  
  return(lprior)
}

#' Evaluate prior of calibrated parameters
#'
#' \code{prior} computes a prior value for one (or multiple) parameter set(s).
#' @param v_params Vector (or matrix) of model parameters 
#' @return 
#' A scalar (or vector) with prior values.
#' @examples
#' v_param_names  <- c("r_Mets_DC",
#'                     "r_BMneg_Mets",
#'                     "hr_BMpos_Mets")
#' v_lb <- c(r_Mets_DC     = 0.037,               # lower bound
#'           r_BMneg_Mets  = 0.001,
#'           hr_BMpos_Mets = 1.58))  
#' v_ub <- c(r_Mets_DC     = -log(1-(1-0.03))/60, # upper bound 
#'           r_BMneg_Mets  = 0.03,
#'           hr_BMpos_Mets = 4.72) 
#' prior(v_params = sample.prior(n_samp = 5))
#' @export
prior <- function(v_params) { 
  v_prior <- exp(log_prior(v_params)) 
  return(v_prior)
}

#' Log-likelihood function for a parameter set
#'
#' \code{log_lik} computes a log-likelihood value for one (or multiple) 
#' parameter set(s).
#'
#' @param v_params Vector (or matrix) of model parameters.
#' @param l_params_all List with all parameters of the decision model. 
#' @return 
#' A scalar (or vector) with log-likelihood values.
#' @importFrom stats dnorm dunif quantile qunif rbeta rgamma sd
#' @examples 
#' \dontrun{
#' v_param_names  <- c("r_Mets_DC",
#'                     "r_BMneg_Mets",
#'                     "hr_BMpos_Mets")
#' v_lb <- c(r_Mets_DC     = 0.037,               # lower bound
#'           r_BMneg_Mets  = 0.001,
#'           hr_BMpos_Mets = 1.58))  
#' v_ub <- c(r_Mets_DC     = -log(1-(1-0.03))/60, # upper bound 
#'           r_BMneg_Mets  = 0.03,
#'           hr_BMpos_Mets = 4.72) 
#' v_target_names <- c("DFS", "OS", "DSS")
#' n_target       <- length(v_target_names)
#' log_lik(v_params = sample.prior(n_samp = 2))
#' }
#' @export
log_lik <- function(v_params,
                    l_params_all = load_all_params()){ # User defined
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  
  n_samp         <- nrow(v_params)
  v_target_names <- c("DFS", "OS", "DSS")
  n_target       <- length(v_target_names)
  v_llik <- matrix(0, nrow = n_samp, ncol = n_target) 
  colnames(v_llik) <- v_target_names
  v_llik_overall <- numeric(n_samp)
  for(j in 1:n_samp) { # j=1
    jj <- tryCatch( { 
      ###   Run model for parameter set "v_params" ###
      l_model_res <- calibration_out(v_params_calib = v_params[j, ], 
                                     l_params_all = l_params_all)
      
      ###  Calculate log-likelihood of model outputs to targets  ###
      ## TARGET 1: Disease-free survival ("DFS")
      ## Normal log-likelihood  
      v_llik[j, "DFS"] <- sum(dnorm(x = df_calibration_targets$S[1:2],
                                     mean = c(l_model_res$v_dfs_BMpos, 
                                              l_model_res$v_dfs_BMneg),
                                     sd = df_calibration_targets$se[1:2],
                                     log = T))
      
      ## TARGET 2: Overall survival ("OS")
      ## Normal log-likelihood
      v_llik[j, "OS"] <- sum(dnorm(x = df_calibration_targets$S[3:4],
                                     mean = c(l_model_res$v_os_BMpos, 
                                              l_model_res$v_os_BMneg),
                                     sd = df_calibration_targets$se[3:4],
                                     log = T))

      ## TARGET 3: Disease-specific survival ("DSS")
      ## Normal log-likelihood
      v_llik[j, "DSS"] <- sum(dnorm(x = df_calibration_targets$S[5:6],
                                         mean = c(l_model_res$v_dss_BMpos,
                                                  l_model_res$v_dss_BMneg),
                                         sd = df_calibration_targets$se[5:6],
                                         log = T))
      
      ## OVERALL
      ## can give different targets different weights (user must change this)
      v_weights <- rep(1, n_target)
      ## weighted sum
      v_llik_overall[j] <- v_llik[j, ] %*% v_weights
    }, error = function(e) NA) 
    if(is.na(jj)) { v_llik_overall <- -Inf }
  } ## End loop over sampled parameter sets
  
  ## return GOF
  return(v_llik_overall)
}

#' Parallel evaluation of log-likelihood function for a sets of parameters
#'
#' \code{log_lik_par} computes a log-likelihood value for one (or multiple) 
#' parameter set(s) using parallel computation.
#'
#' @param v_params Vector (or matrix) of model parameters.
#' @param l_params_all List with all parameters of the decision model. 
#' @return 
#' A scalar (or vector) with log-likelihood values.
#' @importFrom stats dnorm dunif quantile qunif rbeta rgamma sd
#' @import doParallel
#' @examples 
#' \dontrun{
#' #' v_param_names  <- c("r_Mets_DC",
#'                     "r_BMneg_Mets",
#'                     "hr_BMpos_Mets")
#' v_lb <- c(r_Mets_DC     = 0.037,               # lower bound
#'           r_BMneg_Mets  = 0.001,
#'           hr_BMpos_Mets = 1.58))  
#' v_ub <- c(r_Mets_DC     = -log(1-(1-0.03))/60, # upper bound 
#'           r_BMneg_Mets  = 0.03,
#'           hr_BMpos_Mets = 4.72) 
#' v_target_names <- c("DFS", "OS", "DSS")
#' n_target       <- length(v_target_names)
#' log_lik_par(v_params = sample.prior(n_samp = 2))
#' }
#' @export
log_lik_par <- function(v_params, 
                        l_params_all = load_all_params(),
                        ...) { 
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  
  n_samp <- nrow(v_params)
  
  ### Get OS
  os <- darthtools::get_os()
  
  no_cores <- parallel::detectCores() - 1
  
  print(paste0("Parallelized Likelihood calculations on ", os, " using ", no_cores, " cores"))
  
  n_time_init_likpar <- Sys.time()
  
  if(os == "macosx"){
    # Initialize cluster object
    cl <- parallel::makeForkCluster(no_cores) 
    doParallel::registerDoParallel(cl)
    v_llk <- foreach::foreach(i = 1:n_samp, .combine = c) %dopar% {
      log_lik(v_params[i, ], l_params_all) # i = 1
    }
    n_time_end_likpar <- Sys.time()
  }
  
  if(os == "windows"){
    # Initialize cluster object
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)
    opts <- list(attachExportEnv = TRUE)
    v_llk <- foreach::foreach(i = 1:n_samp, .combine = c,
                              .export = ls(globalenv()),
                              .packages=c(),
                              .options.snow = opts) %dopar% {
                                CBMmod::log_lik(v_params[i, ], ...)
                              }
    n_time_end_likpar <- Sys.time()
  }
  if(os == "linux"){
    # Initialize cluster object
    cl <- parallel::makeCluster(no_cores)
    doMC::registerDoMC(cl)
    v_llk <- foreach::foreach(i = 1:n_samp, .combine = c) %dopar% {
      CBMmod::log_lik(v_params[i, ], ...)
    }
    n_time_end_likpar <- Sys.time()
  }
  
  parallel::stopCluster(cl)
  n_time_total_likpar <- difftime(n_time_end_likpar, n_time_init_likpar, 
                                  units = "hours")
  print(paste0("Runtime: ", round(n_time_total_likpar, 2), " hrs."))
  #-# Try this: # PO
  rm(cl)        # PO
  gc()          # PO
  #-#           # PO
  return(v_llk)
}


#' Likelihood
#'
#' \code{likelihood} computes a likelihood value for one (or multiple) 
#' parameter set(s).
#'
#' @param v_params Vector (or matrix) of model parameters. 
#' @return 
#' A scalar (or vector) with likelihood values.
#' @examples
#' v_param_names  <- c("r_Mets_DC",
#'                     "r_BMneg_Mets",
#'                     "hr_BMpos_Mets")
#' v_lb <- c(r_Mets_DC     = 0.037,               # lower bound
#'           r_BMneg_Mets  = 0.001,
#'           hr_BMpos_Mets = 1.58)  
#' v_ub <- c(r_Mets_DC     = -log(1-(1-0.03))/60, # upper bound 
#'           r_BMneg_Mets  = 0.03,
#'           hr_BMpos_Mets = 4.72) 
#' v_target_names <- c("DFS", "OS", "DSS")
#' likelihood(v_params = sample.prior(n_samp = 2))
#' @export
likelihood <- function(v_params = v_params){ 
  v_like <- exp(log_lik_par(v_params)) 
  return(v_like)
}

#' Evaluate log-posterior of calibrated parameters
#'
#' \code{log_post} Computes a log-posterior value for one (or multiple) 
#' parameter set(s) based on the simulation model, likelihood functions and 
#' prior distributions.
#' @param v_params Vector (or matrix) of model parameters 
#' @return 
#' A scalar (or vector) with log-posterior values.
#' @examples 
#' v_param_names  <- c("r_Mets_DC",
#'                     "r_BMneg_Mets",
#'                     "hr_BMpos_Mets")
#' v_lb <- c(r_Mets_DC     = 0.037,               # lower bound
#'           r_BMneg_Mets  = 0.001,
#'           hr_BMpos_Mets = 1.58))  
#' v_ub <- c(r_Mets_DC     = -log(1-(1-0.03))/60, # upper bound 
#'           r_BMneg_Mets  = 0.03,
#'           hr_BMpos_Mets = 4.72) 
#' v_target_names <- c("DFS", "OS", "DSS")
#' n_target       <- length(v_target_names)
#' log_post(v_params = sample.prior(n_samp = 2))
#' @export
log_post <- function(v_params) { 
  v_lpost <- log_prior(v_params) + log_lik(v_params)
  return(v_lpost) 
}

#' Evaluate posterior of calibrated parameters
#'
#' \code{posterior} computes a posterior value for one (or multiple) parameter 
#' set(s).
#' @param v_params Vector (or matrix) of model parameters 
#' @return 
#' A scalar (or vector) with posterior values.
#' @examples
#' \dontrun{
#' v_param_names  <- c("r_Mets_DC",
#'                     "r_BMneg_Mets",
#'                     "hr_BMpos_Mets")
#' v_lb <- c(r_Mets_DC     = 0.037,               # lower bound
#'           r_BMneg_Mets  = 0.001,
#'           hr_BMpos_Mets = 1.58))  
#' v_ub <- c(r_Mets_DC     = -log(1-(1-0.03))/60, # upper bound 
#'           r_BMneg_Mets  = 0.03,
#'           hr_BMpos_Mets = 4.72) 
#' v_target_names <- c("DFS", "OS", "DSS")
#' n_target       <- length(v_target_names)
#' posterior (v_params = sample.prior(n_samp = 5))
#' }
#' @export
posterior <- function(v_params) { 
  v_posterior <- exp(log_post(v_params)) 
  return(v_posterior)
}