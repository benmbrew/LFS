### This script will analyze results from training and testing on cases 

##########
# initiate library
##########
library(tidyverse)

##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
results_folder <- paste0(project_folder, '/Scripts/predict_age/Results')
reg_folder <- paste0(results_folder, '/reg_results_05')

##########
# set fixed variables
##########
# no_transform
method = 'funnorm'
type = 'no_transform'

#
# list 4 = folds
# list 2 = mod, mod_resid
# list 6 = alpha, lambda_value, importance, cases_cor, age_cor, temp.non_zero_coeff
###########
# load in reg results 
###########
fold_results <- list()
reg_results <- list()
temp.dims <- list()
iter_length <- 30
folds <- 5
# i = 1
# j = 1

# loop through seed number and read in data
# the second element of the list can be removed (its a duplicate of the 4th element in each fold)

# get dims 
for (i in 1:iter_length) {
  # read in raw results
  temp.result <- readRDS(paste0(reg_folder, '/train_test_', i, '_', type, '_', method,'.rda' ))
  
  temp.dims[[i]] <- temp.result[[2]]
  
  temp.result_folds <- temp.result[[1]]

  for (j in 1:folds){
    # get fold - returns list of 4
    temp.fold <- temp.result_folds[[j]]

    
    # get results - list of 2 - normal and resid
    temp.fold_norm <- temp.fold[[1]]
    temp.fold_resid <- temp.fold[[2]]
    
    #alpha, lambda_value, importance, cases_cor, age_cor, controls_cor, valid_cor, temp.non_zero_coeff
    
    # remove 3rd element 
    temp.fold_norm[[3]] <- NULL
    temp.fold_resid[[3]] <- NULL
    

    # list of 5 - alpha, lambda_value, cases_cor, age_cor, temp.non_zero_coeff
    temp.fold_norm_dat <- as.data.frame(t(do.call(rbind, temp.fold_norm)))
    temp.fold_resid_dat <- as.data.frame(t(do.call(rbind, temp.fold_resid)))
    
    # add in resid and norm
    temp.fold_norm_dat$V8 <- 'norm'
    temp.fold_resid_dat$V8 <- 'resid'
    
    

    # get column names 
    colnames(temp.fold_norm_dat) <- c('alpha', 
                                      'lambda', 
                                      'onset_correlation', 
                                      'age_correlation', 
                                      'controls_cor',
                                      'valid_cor',
                                      'vars', 
                                      'type')
    
    colnames(temp.fold_resid_dat) <- c('alpha', 
                                       'lambda', 
                                       'onset_correlation', 
                                       'age_correlation', 
                                       'controls_cor',
                                       'valid_cor',
                                       'vars', 
                                       'type')

    # add in indicator for norm and resid
    temp.fold_norm_dat$seed_num <- i
    temp.fold_resid_dat$seed_num <- i
    

    # combine 
    temp.result_folds_dat <- rbind(temp.fold_norm_dat,
                               temp.fold_resid_dat)
    
    # store in list
    # fold_results[[j]] <- temp.result_folds
    fold_results[[j]] <- temp.result_folds_dat
    
    
  }

  # collapse fold_results
  temp.collpased <- do.call(rbind, fold_results)
  reg_results[[i]] <- temp.collpased
  print(i)
  
}

# get temp.dims
dim_dat <- unlist(temp.dims)

# collpase list into data frame 
result_table <- as.data.frame(do.call(rbind, reg_results))

# result_table make dim  column
result_table$mean_dim <- dim_dat

# unlist alpha 
result_table$alpha <- unlist(result_table$alpha)
##########
# analyze results 
##########

# first group alpha and get mean correlation 
alpha_result <- result_table %>%
  group_by(type, alpha) %>%
  summarise(mean_cor = mean(onset_correlation),
            mean_cor_age = mean(age_correlation),
            mean_cor_controls = mean(controls_cor),
            mean_cor_valid = mean(valid_cor),
            sd_cor = sd(onset_correlation),
            sd_cor_age = sd(age_correlation),
            sd_cor_controls = sd(controls_cor),
            sd_cor_valid = sd(valid_cor),
            mean_vars_import = mean(vars),
            mean_dim = mean(mean_dim))

