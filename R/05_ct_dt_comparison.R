library(dplyr)
library(msm)
library(reshape)
library(ggplot2)

# test whether loading parameters works
source("R/01_model_inputs_functions.R")
file_mort = "data/Annual_basemort_rate_2010.csv"
file_exmort = "data/heavydrinking_mort_rateratio.csv"
l_params_all<-load_all_params(file_mort=file_mort)

# run alcohol model with continuous-time
source("R/02_decision_model_functions.R")
source("R/03_alcohol_model_functions.R")
ct_trace <- run_drinkmc(l_params_all)


# run alcohol model with discrete-time
source("R/02_decision_model_functions_DT.R")
dt_trace <- run_drinkmc(l_params_all)

# compare outcomes between continuous and discrete-time models
source("R/06_markov_outcome_calculations.R")
with(as.list(l_params_all),{
  ## 1. Compare Markov trace
  ct_trace_dt <- data.frame(Age=seq(n_age_init,n_age_max-1,by=1), ct_trace)
  colnames(ct_trace_dt) <- c("Age",v_names_states)
  ct_trace_dt_t <- melt(ct_trace_dt, id.vars="Age")
  ct_trace_dt_t$type = 'continuous'

  dt_trace_dt <- data.frame(Age=seq(n_age_init,n_age_max-1,by=1), dt_trace)
  colnames(dt_trace_dt) <- c("Age",v_names_states)
  dt_trace_dt_t <- melt(dt_trace_dt, id.vars="Age")
  dt_trace_dt_t$type = 'discrete'
  
  trace_all <- rbind(ct_trace_dt_t, dt_trace_dt_t)
  
  ggplot(data=trace_all)+
    facet_wrap(~variable)+
    geom_line(aes(x=Age,y=value*100, color=type))+
    ylab("% of population")
  ggsave("figures/mc_trace_comparison.pdf")
  
  ## 2. Prevalence of current drinkers
  # Calculate prevalence
  ct_prev <- cal_prev(ct_trace)
  dt_prev <- cal_prev(dt_trace)
  # Add labels of model type
  ct_prev$type = 'continuous'
  dt_prev$type = 'discrete'
  
  prev_all <- rbind(ct_prev, dt_prev)
  
  ggplot(data = prev_all, aes(x=Age,y=Prev,fill=type))+
    geom_bar(position="dodge", stat="identity")+
    ylab("Prevalence of drinkers")+
    scale_y_continuous(limits=c(0,max(prev_all$Prev)+0.01), 
                                breaks = seq(0, (max(prev_all$Prev)+0.01), by = 0.02), 
                                expand = c(0,0))+
    theme_bw()
  ggsave("figures/prev_comparison.pdf")
  
  ## 3. Input parameters
  
  })

dt_trace_dt <- data.frame(dt_trace)
colnames(dt_trace_dt) <- l_params_all$v_names_states

ggplot(data = ct_trace)+
  