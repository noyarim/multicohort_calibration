#' Transition probability (TP) array for alcohol model
#' 
#' \code{create_tpm} generates the transition probability array for the CBM model.
#'
#' @param l_params_all List with all parameters of decision model
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE
#' @param Trt is this the Treat All strategy? (default is FALSE)
#' @return 
#' The transition probability array for the alcohol model.
#' @export
create_tpm <- function(l_params_all, err_stop = FALSE, verbose = FALSE, 
                   Trt = FALSE
){
  ### Depends on:
  ###   - Package `msm` to compute the matrix exponential
  with(as.list(l_params_all), {
    # Sanity Check
    if (n_cycles_year <= 0) {
      stop("Each year must have at least one cycle")
    }
    # Discrete-time model -----------------------------------------------------
    ## Time scale of the model
    t_scale <- 1/n_cycles_year
    
    # Continuous-time model ---------------------------------------------------
    ## Create transition rates
    
    # From ND ----
    ## to CD (adjusted by treatment status) ----
    v_r_ND_CD <- c(rep(r_ND_CD_5, 5), )
    v_r_ND_CD_Rx <- v_r_ND_CD*(1-Trt) + # No treatment
      v_r_ND_CD*hr_ND_CD_Rx*(Trt)     # Treatment
    ## To Death ----
    v_r_ND_Death <- 
      
    # From CD  ----
    ## to FD ----
    v_r_CD_FD <- 
    ## To Death ----
    
    # From FD  ----
    ## to CD ----
    v_r_FD_CD <- 
    ## To Death ----
    v_r_FD_Death <- 
  
    ## Transition Intensity Array -------------------------------------------
    a_Q <- array(0, dim = list(n_states, n_states, n_cycles), 
                 dimnames = list(v_names_states, v_names_states, v_age_names))
    
    ## Fill in array
    # From ND
    a_Q["ND", "CD", ]   <- 
    a_Q["ND", "Dead", ]  <- 
    # From CD
    a_Q["CD", "FD", ]   <- 
    a_Q["CD", "Dead", ]    <- 
    # From FD
    a_Q["FD", "CD", ]   <- 
    a_Q["FD", "Dead", ] <- 
    
    ### Compute transition probability array
    a_P <- sapply(seq_along(v_age_names), 
                  function(x) msm::MatrixExp(mat = a_G[,,x], t = t_scale), 
                  simplify="array")
    
      
    ## Check if Transition Probability array is valid
    # darthtools::check_transition_probability(a_P, err_stop = err_stop, 
    #                                          verbose = verbose)
    # darthtools::check_sum_of_transition_array(a_P, n_states = n_states, 
    #                                           n_cycles = n_cycles,
    #                                           err_stop = err_stop, 
    #                                           verbose = verbose)
    
    # Return transition probability array -------------------------------------
    return(a_P) 
  }
  )
}

