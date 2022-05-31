#' 4-state alcohol model
#' 
#' \code{create_tpm} generates the transition probability array for the CBM model.
#'
#' @param a_P Array of age-specific transition matrices. 
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE
#' @param Trt is this the Treat All strategy? (default is FALSE)
#' @return 
#' The markov trace for the alcohol model. 
#' @export

run_drinkmc <- function(l_param_all,  err_stop = FALSE, verbose = FALSE, 
                       Trt = FALSE
){
  with(as.list(l_params_all), {
    # Sanity Check
    if (n_cycles_year <= 0) {
      stop("Each year must have at least one cycle")
    }
    # 
    a_P = create_tpm(l_params_all, err_stop = err_stop, verbose = verbose, 
                     Trt = Trt)
    # Discrete-time alcohol model -----------------------------------------------------
    mc_trace <- matrix(0,nrow = n_cycles_year, ncol=n_states) # Markov trace matrix
    
    for (t in n_cycles_year){
      mc_trace[t+1,] <- mc_trace[t,] %*% a_P[,,t]
    }

    # Return markov trace matrix -------------------------------------
    return(mc_trace) 
  }
  )
}