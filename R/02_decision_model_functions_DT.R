#' Transition probability (TP) array for alcohol model
#' 
#' \code{create_tpm} generates the transition probability array for the alcohol 
#' model.This code calculate discrete-time transition probability by using simple
#' exponential transformation of rates.
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
    v_r_ND_CD <- c(rep(r_ND_CD_16, 4),rep(r_ND_CD_20, 4),rep(r_ND_CD_24, 4),
                   rep(r_ND_CD_28, 4),rep(r_ND_CD_32, 4),rep(r_ND_CD_36, 4),
                   rep(r_ND_CD_40, 4),rep(r_ND_CD_44, 4),rep(r_ND_CD_48, 4),
                   rep(r_ND_CD_52, 4),rep(r_ND_CD_56, 4),rep(r_ND_CD_60, 4),
                   rep(r_ND_CD_64, 4),rep(r_ND_CD_68, 4),rep(r_ND_CD_72, 28))
    v_r_ND_CD_Rx <- v_r_ND_CD*(1-Trt) + # No treatment
      v_r_ND_CD*hr_ND_CD_Rx*(Trt)     # Treatment
    ## To Death ----
    v_r_ND_Death <- v_r_mort_year_age_cycle
      
    # From CD  ----
    ## to FD ----
    v_r_CD_FD <- c(rep(r_CD_FD_16, 4),rep(r_CD_FD_20, 4),rep(r_CD_FD_24, 4),
                   rep(r_CD_FD_28, 4),rep(r_CD_FD_32, 4),rep(r_CD_FD_36, 4),
                   rep(r_CD_FD_40, 4),rep(r_CD_FD_44, 4),rep(r_CD_FD_48, 4),
                   rep(r_CD_FD_52, 4),rep(r_CD_FD_56, 4),rep(r_CD_FD_60, 4),
                   rep(r_CD_FD_64, 4),rep(r_CD_FD_68, 4),rep(r_CD_FD_72, 28))
    ## To Death ----
    v_r_CD_Death <- v_r_mort_year_age_cycle * v_hr_mort_age
    
    # From FD  ----
    ## to CD ----
    v_r_FD_CD <- c(rep(r_FD_CD_16, 4),rep(r_FD_CD_20, 4),rep(r_FD_CD_24, 4),
                   rep(r_FD_CD_28, 4),rep(r_FD_CD_32, 4),rep(r_FD_CD_36, 4),
                   rep(r_FD_CD_40, 4),rep(r_FD_CD_44, 4),rep(r_FD_CD_48, 4),
                   rep(r_FD_CD_52, 4),rep(r_FD_CD_56, 4),rep(r_FD_CD_60, 4),
                   rep(r_FD_CD_64, 4),rep(r_FD_CD_68, 4),rep(r_FD_CD_72, 28))
    ## To Death ----
    v_r_FD_Death <- v_r_mort_year_age_cycle
  
    ## Transition Intensity Array -------------------------------------------
    a_Q <- array(0, dim = list(n_states, n_states, n_cycles), 
                 dimnames = list(v_names_states, v_names_states, v_age_names))
    
    ## Fill in array
    # From ND
    #a_Q["ND", "ND", ]   <- -(v_r_ND_CD+v_r_ND_Death)
    a_Q["ND", "CD", ]   <- v_r_ND_CD
    a_Q["ND", "Dead", ]  <- v_r_ND_Death
    # From CD
    #a_Q["CD", "CD", ]   <- -(v_r_CD_FD+v_r_CD_Death) 
    a_Q["CD", "FD", ]   <- v_r_CD_FD
    a_Q["CD", "Dead", ]    <- v_r_CD_Death
    # From FD
    #a_Q["FD", "FD", ]   <- -(v_r_FD_CD+v_r_FD_Death) 
    a_Q["FD", "CD", ]   <- v_r_FD_CD
    a_Q["FD", "Dead", ] <- v_r_FD_Death
    
    ### Compute transition probability array
    a_P <- 1-exp(-a_Q * t_scale)
    
    a_P["ND", "ND", ]   <- 1-a_P["ND", "CD", ]-a_P["ND","Dead",]
    a_P["CD", "CD", ]   <- 1-a_P["CD", "FD",]-a_P["CD","Dead",]
    a_P["FD", "FD", ]   <- 1-a_P["FD", "CD", ]-a_P["FD","Dead",]
    a_P["Dead","Dead",] <- 1
    
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

