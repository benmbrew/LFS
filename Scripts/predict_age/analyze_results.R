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
reg_folder <- paste0(results_folder, '/reg_results_new')

# list 4 = folds
# list 2 = mod, mod_resid
# alpha, lambda_value, importance, cases_cor, age_cor, temp.non_zero_coeff
###########
# load in reg results 
###########
fold_results <- list()
reg_results <- list()
iter_length <- 50
folds <- 4
i = 1
j = 1



# loop through seed number and read in data
for (i in 1:iter_length) {
  # read in raw results
  temp.result <- readRDS(paste0(reg_folder, '/train_test_', i, '.rda' ))
  
  for (j in 1:folds){
    # get fold - returns list of 4
    temp.fold <- temp.result[[j]]
    
    # get results - list of 2 - normal and resid
    temp.fold_norm <- temp.fold[[1]]
    temp.fold_resid <- temp.fold[[2]]
    
    # remove 3rd element 
    temp.fold_norm[[3]] <- NULL
    temp.fold_resid[[3]] <- NULL
    
    # list of 5 - alpha, lambda_value, cases_cor, age_cor, temp.non_zero_coeff
    temp.fold_norm <- as.data.frame(t(do.call(rbind, temp.fold_norm)))
    temp.fold_resid <- as.data.frame(t(do.call(rbind, temp.fold_resid)))
    
    # get column names 
    colnames(temp.fold_norm) <- c('alpha', 'lambda', 'onset_correlation', 'age_correlation', 'vars')
    colnames(temp.fold_resid) <- c('alpha', 'lambda', 'onset_correlation', 'age_correlation', 'vars')
    
    # add in indicator for norm and resid
    temp.fold_norm$type <- 'normal'
    temp.fold_resid$type <- 'resid'
    
    # combine 
    temp.result_folds <- rbind(temp.fold_norm,
                               temp.fold_resid)
    
    # store in list
    fold_results[[j]] <- temp.result_folds
    
  }

  # collapse fold_results
  temp.collpased <- do.call(rbind, fold_results)
  reg_results[[i]] <- temp.collpased
  print(i)
}


# collpase list into data frame 
result_table <- as.data.frame(do.call(rbind, reg_results))

# unlist alpha 
result_table$alpha <- unlist(result_table$alpha)
##########
# analyze results 
##########

# first group alpha and get mean correlation 
alpha_result <- result_table %>%
  group_by(alpha, type) %>%
  summarise(mean_cor = mean(onset_correlation),
            mean_cor_age = mean(age_correlation))

