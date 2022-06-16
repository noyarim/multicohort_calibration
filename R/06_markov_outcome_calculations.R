#' Calculate secondary outcomes from markov trace
#'
#' \code{cal_prev} calculates age-specific prevalence of current drinkers
#' 
#' @param mc_trace Markov trace outcome from a markov model
#' @return prevalence of drinkers by age
#' @export
cal_prev <- function(mc_trace){
  # Calculate prevalence of current drinkers
  prev_CD <- mc_trace[,2] / (1-mc_trace[,ncol(mc_trace)])
  prev_ND <- mc_trace[,1] / (1-mc_trace[,ncol(mc_trace)])
  prev_FD <- mc_trace[,3] / (1-mc_trace[,ncol(mc_trace)])
  # Group into 4-year age buckets
  age_groups <- c("16-19","20-23","24-27","28-31","32-35",
                  "36-40","44-47","48-51","52-55","56-59",
                  "60-63","64-67","68-71","72-")
  n_agegroups = length(age_groups)
  prev_CD_bygroup <- data.frame(Age=age_groups,
                            Prev=rep(0,n_agegroups))
  prev_ND_bygroup <- data.frame(Age=age_groups,
                                Prev=rep(0,n_agegroups))
  prev_FD_bygroup <- data.frame(Age=age_groups,
                                Prev=rep(0,n_agegroups))
  # Calculate the average prevalence by age groups
  for (i in 1:n_agegroups){
    if (i < n_agegroups){
      prev_CD_bygroup$Prev[i] <- mean(prev_CD[(4*(i-1)+1):(4*i)])
      prev_ND_bygroup$Prev[i] <- mean(prev_ND[(4*(i-1)+1):(4*i)])
      prev_FD_bygroup$Prev[i] <- mean(prev_FD[(4*(i-1)+1):(4*i)])
    }else{
      # age group older than 72
      prev_CD_bygroup$Prev[i] <- mean(prev_CD[(4*(i-1)+1):nrow(mc_trace)])
      prev_ND_bygroup$Prev[i] <- mean(prev_ND[(4*(i-1)+1):nrow(mc_trace)])
      prev_FD_bygroup$Prev[i] <- mean(prev_FD[(4*(i-1)+1):nrow(mc_trace)])
    }
  }
  
  model.res <- list(Prev_CD = prev_CD_bygroup,
                    Prev_ND = prev_ND_bygroup,
                    Prev_FD = prev_FD_bygroup)
  return(model.res)
  
}

#'
#' \code{cal_dwell} calculates dwelling time in the current drinker state
#' 
#' @param mc_trace Markov trace outcome from a markov model
#' @return dwelling time in the current drinker state
#' @export
cal_dwell <- function(mc_trace){
  
}