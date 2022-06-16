#' Prepare R data for target
#'
#' \code{prep_targ} prepare calibration targets data in one list 
#' 
#' @return R data (list) for calibration targets
#' @export
#' 
prep_target <- function(){
  
  cd_prev <- readxl:read_xlsx("data/target_beta_data_evenyr.xlsx", sheet='p')
  cd_prev_n <- readxl:read_xlsx("data/target_beta_data_evenyr.xlsx", sheet='n')
  #cd_prev_a <- readxl:read_xlsx("data/target_beta_data_evenyr.xlsx", sheet='a')
  #cd_prev_b <- readxl:read_xlsx("data/target_beta_data_evenyr.xlsx", sheet='b')
  nd_prev <- readxl::read_xlsx("data/target_beta_data_nonhvy_evenyr.xlsx", sheet="p_never")
  nd_prev_n <- readxl:read_xlsx("data/target_beta_data_evenyr.xlsx", sheet='n_never')
  #nd_prev_a <- readxl::read_xlsx("data/target_beta_data_nonhvy_evenyr.xlsx", sheet="a_never")
  #nd_prev_b <- readxl::read_xlsx("data/target_beta_data_nonhvy_evenyr.xlsx", sheet="b_never")
  fd_prev <- readxl::read_xlsx("data/target_beta_data_nonhvy_evenyr.xlsx", sheet="p_nonhvy")
  fd_prev_n <- readxl::read_xlsx("data/target_beta_data_nonhvy_evenyr.xlsx", sheet="n_nonhvy")
  #fd_prev_a <- readxl::read_xlsx("data/target_beta_data_nonhvy_evenyr.xlsx", sheet="a_nonhvy")
  #fd_prev_b <- readxl::read_xlsx("data/target_beta_data_nonhvy_evenyr.xlsx", sheet="b_nonhvy")  
  
  alc_targets <- list(Prev_CD = cd_prev, 
                      N_CD = cd_prev_n, 
                      Prev_ND = nd_prev, 
                      N_ND = nd_prev_n,
                      Prev_FD = fd_prev, 
                      N_FD = fd_prev_n)
  
  return(alc_targets)
}

  
