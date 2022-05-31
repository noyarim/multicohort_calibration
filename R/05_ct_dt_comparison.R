library(dplyr)

# test whether loading parameters works
source("R/01_model_inputs_functions.R")
file_mort = "data/Annual_basemort_rate_2010.csv"
load_all_params(file_mort=file_mort)

#