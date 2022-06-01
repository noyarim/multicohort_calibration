#' Base-case initial parameter set
#'
#' \code{load_params_init} generates the initial values of the Cancer Biomarker 
#' Model (CBMmod)
#' 
#' @param n_age_init Initial age of the cohort.
#' @param n_age_max Oldest age of the cohort.
#' @param n_cycles_year Number of cycles per year.
#' @return 
#' List of all parameters 
#' @export
load_params_init <- function(
  # Initial age of the cohort
  n_age_init = 16,
  # Oldest age of the cohort.
  n_age_max  = 100,
  # Number of cycles per year
  n_cycles_year = 1,
  ## Disease parameters
  # From never drinks to current drinkers
  r_ND_CD_16 = 0.005319545,
  r_ND_CD_20 = 0.017262072,
  r_ND_CD_24 = 0.029033276,
  r_ND_CD_28 = 0.031688322,
  r_ND_CD_32 = 0.029593226,
  r_ND_CD_36 = 0.027498131,
  r_ND_CD_40 = 0.030385763,
  r_ND_CD_44 = 0.033207841,
  r_ND_CD_48 = 0.030643089,
  r_ND_CD_52 = 0.028565132,
  r_ND_CD_56 = 0.024959867,
  r_ND_CD_60 = 0.017192874,
  r_ND_CD_64 = 0.009861861,
  r_ND_CD_68 = 0.002936822,
  r_ND_CD_72 = 0.002936822,
  # From current drinkers to former drinkers
  r_CD_FD_16 = 0.486965245,
  r_CD_FD_20 = 0.297701075,
  r_CD_FD_24 = 0.272752696,
  r_CD_FD_28 = 0.276080088,
  r_CD_FD_32 = 0.279407481,
  r_CD_FD_36 = 0.282734873,
  r_CD_FD_40 = 0.591892876,
  r_CD_FD_44 = 0.53503715,
  r_CD_FD_48 = 0.709829619,
  r_CD_FD_52 = 0.475097113,
  r_CD_FD_56 = 0.324259987,
  r_CD_FD_60 = 0.283815614,
  r_CD_FD_64 = 0.255740758,
  r_CD_FD_68 = 0.094477091,
  r_CD_FD_72 = 0.094477091,
  # From former drinkers to current drinkers
  r_FD_CD_16 = 0.028825709,
  r_FD_CD_20 = 0.070759726,
  r_FD_CD_24 = 0.078689093,
  r_FD_CD_28 = 0.048541515,
  r_FD_CD_32 = 0.047456544,
  r_FD_CD_36 = 0.053302603,
  r_FD_CD_40 = 0.145325323,
  r_FD_CD_44 = 0.08545241,
  r_FD_CD_48 = 0.06100452,
  r_FD_CD_52 = 0.028945683,
  r_FD_CD_56 = 0.017866895,
  r_FD_CD_60 = 0.012490663,
  r_FD_CD_64 = 0.009642906,
  r_FD_CD_68 = 0.004194026,
  r_FD_CD_72 = 0.004194026,
  
  # Hazard ratio for starting drinking among never drinkers under treatment 
  # versus BM-negative patients without treatment.
  hr_ND_CD_Rx = 1.00

  ){
  # Error checking ---------------------------------------------------------
  # Check if the selected distribution is in the available set

  if (!(n_age_init %% 1 == 0)) {
    stop("n_age_init should be an integer")
  }
  if (!(n_age_max %% 1 == 0)) {
    stop("n_age_max should be an integer")
  }
  if (n_age_init < 0 | n_age_init > 100) {
    stop("n_age_init should be greater than 0 and lower than 100")
  }
  if (n_age_max < 0 | n_age_max > 100) {
    stop("n_age_max should be greater than 0 and lower than 100")
  }
  if (!(all(c(r_ND_CD_16,r_ND_CD_20,r_ND_CD_24,r_ND_CD_28,r_ND_CD_32,r_ND_CD_36,
              r_ND_CD_40,r_ND_CD_44,r_ND_CD_48,r_ND_CD_52,r_ND_CD_56,r_ND_CD_60,
              r_ND_CD_64,r_ND_CD_68,r_ND_CD_72) >= 0))) { ## TODO: ADD ONE LINE OF CHECKING ALL ELEMENTS IN THE VECTOR
    stop("r_CD_FD should be a value greater than or euqal to zero")
  }
  if (!(all(c(r_CD_FD_16,r_CD_FD_20,r_CD_FD_24,r_CD_FD_28,r_CD_FD_32,r_CD_FD_36,
              r_CD_FD_40,r_CD_FD_44,r_CD_FD_48,r_CD_FD_52,r_CD_FD_56,r_CD_FD_60,
              r_CD_FD_64,r_CD_FD_68,r_CD_FD_72) >= 0))) { ## TODO: ADD ONE LINE OF CHECKING ALL ELEMENTS IN THE VECTOR
    stop("r_CD_FD should be a value greater than or euqal to zero")
  }
  if (!(all(c(r_FD_CD_16,r_FD_CD_20,r_FD_CD_24,r_FD_CD_28,r_FD_CD_32,r_FD_CD_36,
              r_FD_CD_40,r_FD_CD_44,r_FD_CD_48,r_FD_CD_52,r_FD_CD_56,r_FD_CD_60,
              r_FD_CD_64,r_FD_CD_68,r_FD_CD_72) >= 0))) { ## TODO: ADD ONE LINE OF CHECKING ALL ELEMENTS IN THE VECTOR
    stop("r_CD_FD should be a value greater than or euqal to zero")
  }
  if (n_cycles_year <= 0 | !(n_cycles_year%%1 == 0)) {
    stop("n_cycles_year must be an integer greater than zero")
  }
  
  # Number of cycles
  n_cycles <- (n_age_max - n_age_init)*n_cycles_year

  ### Create list of initial parameters
  l_params_init <- list(
    # Initial and final ages
    n_age_init = n_age_init,
    # Oldest age of the cohort.
    n_age_max  = n_age_max,
    # Number of cycles
    n_cycles = n_cycles,
    # Number of cycles per year
    n_cycles_year = n_cycles_year,
    ## Disease parameters
    r_ND_CD_16 = r_ND_CD_16,
    r_ND_CD_20 = r_ND_CD_20,
    r_ND_CD_24 = r_ND_CD_24,
    r_ND_CD_28 = r_ND_CD_28,
    r_ND_CD_32 = r_ND_CD_32,
    r_ND_CD_36 = r_ND_CD_36,
    r_ND_CD_40 = r_ND_CD_40,
    r_ND_CD_44 = r_ND_CD_44,
    r_ND_CD_48 = r_ND_CD_48,
    r_ND_CD_52 = r_ND_CD_52,
    r_ND_CD_56 = r_ND_CD_56,
    r_ND_CD_60 = r_ND_CD_60,
    r_ND_CD_64 = r_ND_CD_64,
    r_ND_CD_68 = r_ND_CD_68,
    r_ND_CD_72 = r_ND_CD_72,
    # From current drinkers to former drinkers
    r_CD_FD_16 = r_CD_FD_16,
    r_CD_FD_20 = r_CD_FD_20,
    r_CD_FD_24 = r_CD_FD_24,
    r_CD_FD_28 = r_CD_FD_28,
    r_CD_FD_32 = r_CD_FD_32,
    r_CD_FD_36 = r_CD_FD_36,
    r_CD_FD_40 = r_CD_FD_40,
    r_CD_FD_44 = r_CD_FD_44,
    r_CD_FD_48 = r_CD_FD_48,
    r_CD_FD_52 = r_CD_FD_52,
    r_CD_FD_56 = r_CD_FD_56,
    r_CD_FD_60 = r_CD_FD_60,
    r_CD_FD_64 = r_CD_FD_64,
    r_CD_FD_68 = r_CD_FD_68,
    r_CD_FD_72 = r_CD_FD_72,  
    # From former drinkers to current drinkers
    r_FD_CD_16 = r_FD_CD_16,
    r_FD_CD_20 = r_FD_CD_20,
    r_FD_CD_24 = r_FD_CD_24,
    r_FD_CD_28 = r_FD_CD_28,
    r_FD_CD_32 = r_FD_CD_32,
    r_FD_CD_36 = r_FD_CD_36,
    r_FD_CD_40 = r_FD_CD_40,
    r_FD_CD_44 = r_FD_CD_44,
    r_FD_CD_48 = r_FD_CD_48,
    r_FD_CD_52 = r_FD_CD_52,
    r_FD_CD_56 = r_FD_CD_56,
    r_FD_CD_60 = r_FD_CD_60,
    r_FD_CD_64 = r_FD_CD_64,
    r_FD_CD_68 = r_FD_CD_68,
    r_FD_CD_72 = r_FD_CD_72,
    # Hazard ratio for starting drinking among never drinkers under treatment 
    # versus BM-negative patients without treatment.
    hr_ND_CD_Rx = hr_ND_CD_Rx
  )
  return(l_params_init)
}

#' Load mortality data
#'
#' \code{load_mort_data} is used to load age-specific mortality from .csv file 
#' into vector.
#'
#' @param file String with the location and name of the file with mortality 
#' data. If \code{NULL}, \code{v_r_mort_year_sex_age} will be used as default
#' @return 
#' A vector with mortality by age in 2010
#' A vector with hazard rate ratio with heavy drinking by age
#' @export
load_mort_data <- function(file = NULL){
  # Load mortality data from file
  if(!is.null(file)) {
    df_r_mort_age <- read.csv(file = file)}
  else{
    # df_r_mort_by_age <- all_cause_mortality
    df_r_mort_age <- all_cause_mortality
  }
  # Vector with mortality rates
  v_r_mort_year_age  <- dplyr::select(df_r_mort_age, 
                                          .data$Age, .data$Rate)
  
  return(v_r_mort_year_age)
}

#' Load excess mortality data
#'
#' \code{load_mort_data} is used to load age-specific excess mortality due to
#'  heavy drinking from .csv file into vector.
#' @param file String with the location and name of the file with excess mortality 
#' data. If \code{NULL}, \code{v_r_mort_year_sex_age} will be used as default
#' @return 
#' A vector with mortality by age in 2010
#' A vector with hazard rate ratio with heavy drinking by age
#' @export
load_exmort_data <- function(file=NULL){
    # Load excess mortality data from file
    if(!is.null(file)) {
      df_hrr_mort_age <- read.csv(file = file)}
    else{
      # df_r_mort_by_age <- all_cause_mortality
      df_hrr_mort_age <- hrr_alcohol_age
    }
    # Vector with hazard rate ratios
    df_hrr_mort_age  <- dplyr::select(df_hrr_mort_age, 
                                        .data$Age, .data$Ratio)
    
    return(df_hrr_mort_age)
  }
#' Load all parameters
#'
#' \code{load_all_params} loads all parameters for the decision model from multiple sources and creates a list.
#'
#' @param file.init String with the location and name of the file with initial set of parameters
#' @param file.mort String with the location and name of the file with mortality data
#' @return 
#' A list of all parameters used for the decision model.
#' @export
load_all_params <- function(l_params_init = NULL,
                            file_init = NULL,
                            file_mort = NULL){ # User defined
  # Load initial set of initial parameters from .csv file -------------------
  if(is.null(l_params_init)) {
    l_params_init <- load_params_init()
  } else{
    l_params_init <- l_params_init
  }
  
  # All-cause age-specific mortality from .csv file -------------------------
  v_r_mort_year_age <- load_mort_data(file_mort)
  v_hr_mort_age <- load_exmort_data(file_exmort)
  
  l_params_all <- with(as.list(l_params_init), {
    
    # General setup -----------------------------------------------------------
    v_names_str <- c("No Treat", "Treat")# CEA strategies
    n_str       <- length(v_names_str) # Number of strategies
    v_age_names <- paste(rep(n_age_init:(n_age_max-1), each = n_cycles_year), 
                         1:n_cycles_year, 
                         sep = ".")
    # Vector with the health states of the model
    v_names_states <- c("ND", "CD",
                        "FD", "Dead") 
    n_states <- length(v_names_states) # number of health states
    
    # Within-cycle correction (WCC) using half-cycle rule
    v_wcc <- darthtools::gen_wcc(n_cycles = n_cycles, 
                                 # method = "Simpson1/3")
                                 method = "half-cycle")
    
    # Filter for selected sexes, ages and years
    v_r_mort_year_age <- v_r_mort_year_age %>%
      dplyr::filter(Age >= (n_age_init-1) & Age < n_age_max) %>%
      dplyr::select(Rate) %>%
      as.matrix()
    
    # Compute the mortality rates 
    v_r_mort_year_age_cycle <- c(rep(v_r_mort_year_age[1],4),
                                 rep(v_r_mort_year_age[2:length(v_r_mort_year_age)],each=5)
                                ) 
    # Excess mortality with drinking
    v_hr_mort_age <- v_hr_mort_age %>%
      dplyr::filter(Age >= n_age_init & Age < n_age_max) %>%
      dplyr::select(Ratio)%>%
      as.matrix()
    
    # Create list with all parameters -----------------------------------------
    l_params_all <- list(
      v_names_str = v_names_str,
      n_str       = n_str,
      n_age_init  = n_age_init, 
      n_cycles    = n_cycles, 
      v_age_names = v_age_names,
      v_names_states = v_names_states,
      n_states = n_states,
      n_cycles_year = n_cycles_year,
      v_r_mort_year_age = v_r_mort_year_age,
      v_r_mort_year_age_cycle = v_r_mort_year_age_cycle,
      v_hr_mort_age = v_hr_mort_age,
      v_wcc = v_wcc
    )
    return(l_params_all)
  }
  )
  
  l_params_all <- c(l_params_all, 
                    l_params_init) # Add initial set of parameters
  
  return(l_params_all)
}

#' Update parameters
#'
#' \code{update_param_list} is used to update list of all parameters with new 
#' values for specific parameters.
#'
#' @param l_params_all List with all parameters of decision model
#' @param params_updated Parameters for which values need to be updated
#' @return 
#' A list with all parameters updated.
#' @export
update_param_list <- function(l_params_all, params_updated){
  
  if (typeof(params_updated)!="list"){
    params_updated <- split(unname(params_updated), names(params_updated)) #converts the named vector to a list
  }
  l_params_all <- modifyList(l_params_all, params_updated) #update the values
  return(l_params_all)
}