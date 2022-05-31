library(lhs)
library(IMIS)
library(matrixStats)

####################################################################
######  Load model as a function  ######
####################################################################
source("markov model here")


####################################################################
######  Specify calibration parameters  ######
####################################################################
# Specify seed (for reproducible sequence of random numbers)
set.seed(072218)

# number of random samples
n.resamp <- 1000

# names and number of input parameters to be calibrated
v.param.names <- c( #r_ND_CD: transition rate from ND to CD
                   "r_ND_CD_16", "r_ND_CD_20",, "r_ND_CD_24","r_ND_CD_28",
                   "r_ND_CD_32","r_ND_CD_36","r_ND_CD_40","r_ND_CD_44",
                   "r_ND_CD_48","r_ND_CD_52","r_ND_CD_56","r_ND_CD_60",
                   "r_ND_CD_64","r_ND_CD_68","r_ND_CD_72", 
                   # r_CD_FD: transition rate from CD to FD
                   "r_CD_FD_16", "r_CD_FD_20",, "r_CD_FD_24","r_CD_FD_28",
                   "r_CD_FD_32","r_CD_FD_36","r_CD_FD_40","r_CD_FD_44",
                   "r_CD_FD_48","r_CD_FD_52","r_CD_FD_56","r_CD_FD_60",
                   "r_CD_FD_64","r_CD_FD_68","r_CD_FD_72")
n.param <- length(v.param.names)

# range on input search space
lb <- c(r = 0.0) # lower bound
ub <- c(r = 1.0) # upper bound

# number of calibration targets
v.target.names <- c("Surv", "Prev", "PropSick")
n.target <- length(v.target.names)

####################################################################
######  Calibrate!  ######
####################################################################
# record start time of calibration
t.init <- Sys.time()

###  Write function to sample from prior ###
sample_prior <- function(n.samp){
  m.lhs.unit   <- randomLHS(n = n.samp, k = n.param)
  m.param.samp <- matrix(nrow = n.samp, ncol = n.param)
  colnames(m.param.samp) <- v.param.names
  for (i in 1:n.param){
    m.param.samp[, i] <- qunif(m.lhs.unit[,i],
                               min = lb,
                               max = ub)
    # ALTERNATIVE prior using beta (or other) distributions
    # m.param.samp[, i] <- qbeta(m.lhs.unit[,i],
    #                            min = 1,
    #                            max = 1)
  }
  return(m.param.samp)
}

# view resulting parameter set samples
pairs.panels(sample_prior(1000))

###  Write functions to evaluate log-prior and prior ###
calc_log_prior <- function(v.params){
  if(is.null(dim(v.params))) { # If vector, change to matrix
    v.params <- t(v.params) 
  }
  n.samp <- nrow(v.params)
  colnames(v.params) <- v.param.names
  lprior <- rep(0, n.samp)
  for (i in 1:n.param){
    lprior <- lprior + dunif(v.params[, i],
                             min = lb,
                             max = ub, 
                             log = T)
    # ALTERNATIVE prior using beta distributions
    # lprior <- lprior + dbeta(v.params[, i],
    #                          shape1 = 1,
    #                          shape2 = 1, 
    #                          log = T)
  }
  return(lprior)
}
# test calc_log_prior()
calc_log_prior(v.params = v.params.test)
calc_log_prior(v.params = sample_prior(10))

# function that calculates the (non-log) prior
calc_prior <- function(v.params) { 
  exp(calc_log_prior(v.params)) 
}
# test calc_prior()
calc_prior(v.params = v.params.test)
calc_prior(v.params = sample_prior(10))


###  Write log-likelihood and likelihood functions to pass to IMIS algorithm  ###
calc_log_lik <- function(v.params){
  # par_vector: a vector (or matrix) of model parameters 
  if(is.null(dim(v.params))) { # If vector, change to matrix
    v.params <- t(v.params) 
  }
  n.samp <- nrow(v.params)
  v.llik <- matrix(0, nrow = n.samp, ncol = n.target) 
  llik.overall <- numeric(n.samp)
  for(j in 1:n.samp) { # j=1
    jj <- tryCatch( { 
      ###   Include other parameters  ###
      v.params.full <- cbind(v.params, xxxxx)
      ###   Run model for parametr set "v.params" ###
      model.res <- run_sick_sicker_markov(v.params.full[j, ])
      
      ###  Calculate log-likelihood of model outputs to targets  ###
      # TARGET 1: Survival ("Surv")
      # log likelihood  
      v.llik[j, 1] <- sum(dnorm(x = SickSicker.targets$Surv$value,
                                mean = model.res$Surv,
                                sd = SickSicker.targets$Surv$se,
                                log = T))
      # TARGET 2: Prevalence ("Prev")
      # log likelihood
      v.llik[j, 2] <- sum(dnorm(x = SickSicker.targets$Prev$value,
                                mean = model.res$Prev,
                                sd = SickSicker.targets$Prev$se,
                                log = T))
      # TARGET 3: Proportion Sick+Sicker who are Sick
      # log likelihood
      v.llik[j, 3] <- sum(dnorm(x = SickSicker.targets$PropSick$value,
                                mean = model.res$PropSick,
                                sd = SickSicker.targets$PropSick$se,
                                log = T))
      # OVERALL 
      llik.overall[j] <- sum(v.llik[j, ])
    }, error = function(e) NA) 
    if(is.na(jj)) { llik.overall <- -Inf }
  } # End loop over sampled parameter sets
  # return LLIK
  return(llik.overall)
}
calc_log_lik(v.params = v.params.test)
calc_log_lik(v.params = sample_prior(10))

# write function to calculate the (non-log) likelihood
calc_likelihood <- function(v.params){ 
  exp(calc_log_lik(v.params)) 
}
calc_likelihood(v.params = v.params.test)
calc_likelihood(v.params = sample_prior(10))


###  Write function to evaluate log-posterior ###
calc_log_post <- function(v.params) { 
  lpost <- calc_log_prior(v.params) + calc_log_lik(v.params)
  return(lpost) 
}
calc_log_post(v.params = v.params.test)
calc_log_post(v.params = sample_prior(10))


###  Bayesian calibration using IMIS  ###
# define three functions needed by IMIS: prior(x), likelihood(x), sample.prior(n)
prior <- calc_prior
likelihood <- calc_likelihood
sample.prior <- sample_prior

# run IMIS
fit.imis <- IMIS(B = 1000, # the incremental sample size at each iteration of IMIS
                 B.re = n.resamp, # the desired posterior sample size
                 number_k = 10, # the maximum number of iterations in IMIS
                 D = 0)

# obtain draws from posterior
m.calib.post <- fit.imis$resample

# Calculate log-likelihood and posterior probability of each sample
m.calib.post <- cbind(m.calib.post, 
                      "Log_likelihood" = calc_log_lik(m.calib.post[,v.param.names]),
                      "Posterior_prob" = exp(calc_log_post(m.calib.post[,v.param.names])))

# normalize posterior probability
m.calib.post[,"Posterior_prob"] <- m.calib.post[,"Posterior_prob"]/sum(m.calib.post[,"Posterior_prob"])

# Calculate computation time
comp.time <- Sys.time() - t.init


