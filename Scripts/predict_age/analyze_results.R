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
data_folder <- paste0(project_folder, '/Data')
results_folder <- paste0(project_folder, '/Scripts/predict_age/Results')
reg_folder <- paste0(results_folder, '/reg_results_05')
model_data <- paste0(data_folder, '/model_data')


##########
# function to loop through folder and get all results for each method, type, control size combination
##########

get_results <- 
  function(method, 
           type, 
           control_size, 
           iter_length = 10, 
           folds = 5) {
    
    fold_results <- list()
    reg_results <- list()
    temp.dims <- list()
    
    # get dims 
    for (i in 1:iter_length) {
      # read in raw results
      temp.result <- readRDS(paste0(reg_folder, '/train_test_', i, '_', type, '_', method, control_size ,'.rda' ))
      
      temp.dims[[i]] <- temp.result[[2]]
      
      temp.result_folds <- temp.result[[1]]
      
      for (j in 1:folds){
        # get fold - returns list of 4
        temp.fold <- temp.result_folds[[j]]
        
        
        # get results - list of 2 - normal and resid
        temp.fold_norm <- temp.fold

        #alpha, lambda_value, importance, cases_cor, age_cor, controls_cor, valid_cor, temp.non_zero_coeff
        
        # remove 3rd element 
        temp.fold_norm[[3]] <- NULL
        # temp.fold_resid[[3]] <- NULL
        
        
        # list of 5 - alpha, lambda_value, cases_cor, age_cor, temp.non_zero_coeff
        temp.fold_norm_dat <- as.data.frame(t(do.call(rbind, temp.fold_norm)))
        # temp.fold_resid_dat <- as.data.frame(t(do.call(rbind, temp.fold_resid)))
        
        # add in resid and norm
        temp.fold_norm_dat$V8 <- 'norm'
        # temp.fold_resid_dat$V8 <- 'resid'
        
        
        
        # get column names 
        colnames(temp.fold_norm_dat) <- c('alpha', 
                                          'lambda', 
                                          'onset_correlation', 
                                          'age_correlation', 
                                          'controls_cor',
                                          'valid_cor',
                                          'vars', 
                                          'type')
        
        # colnames(temp.fold_resid_dat) <- c('alpha', 
        #                                    'lambda', 
        #                                    'onset_correlation', 
        #                                    'age_correlation', 
        #                                    'controls_cor',
        #                                    'valid_cor',
        #                                    'vars', 
        #                                    'type')
        
        # add in indicator for norm and resid
        temp.fold_norm_dat$seed_num <- i
        # temp.fold_resid_dat$seed_num <- i
        
        
        # combine 
        temp.result_folds_dat <- temp.fold_norm_dat
        
        # store in list
        # fold_results[[j]] <- temp.result_folds
        fold_results[[j]] <- temp.result_folds_dat
        
        
      }
      
      # collapse fold_results
      temp.collpased <- do.call(rbind, fold_results)
      reg_results[[i]] <- temp.collpased
      
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
                mean_vars_import = mean(vars),
                mean_dim = mean(mean_dim))
    
    # add identifier
    if(control_size == '') {
      controls_size <- 'new'
    }
    alpha_result$key <- paste0(method, '_', type, control_size)
    
    return(alpha_result)
    
    
}

#########
# load results 
#########

# loop through method types and control sizes and aggregate results
methods <- c('raw', 'quan', 'funnorm')
types <- c('original', 'transform', 'no_transform')
sizes <- c('', '_old', '_full')

temp_i <- 
  temp_j <- 
  temp_k <- list()

for (i in 1:length(methods)) {
  
  for (j in 1:length(types)) {
    
    for (k in 1:length(sizes)) {
      
      temp_k[[k]] <- get_results(method = methods[i],
                                 type = types[j],
                                 control_size = sizes[k],
                                 iter_length = 10, 
                                 folds = 5)
      
    }
    
    temp_j[[j]] <- do.call(rbind, temp_k)
    
  }
  
  temp_i[[i]] <- do.call(rbind, temp_j)
  
  print(paste0(methods[i]))
}



# get full results
full_results <- as.data.frame(do.call(rbind, temp_i))

# remove residual models
full_results <- full_results[full_results$type != 'resid',]
full_results$norm <- NULL

# subset to mean_cor, mean_cor_controls, and mean_valid, key
results <- full_results[, c('alpha', 
                            'mean_cor', 
                            'mean_cor_controls',
                            'mean_cor_valid', 
                            'key')]

# group by key and get max mean_cor with corresponding alpha
alpha_max <- full_results %>%
  group_by(key) %>% top_n(1, mean_cor)

saveRDS(alpha_max, paste0(model_data, '/alpha_max.rda'))

# group by key 
result_key <-
  results %>% 
  group_by(key) %>%
  summarise(mean_cor = mean(mean_cor, na.rm = T),
            mean_cor_con = mean(mean_cor_controls, na.rm = T),
            mean_cor_val = mean(mean_cor_valid, na.rm = T))
