
library(plotly)
library(scatterplot3d)
library(car)
library(rgl)
library(ROCR)
library(caret)
library(pROC)
library(dplyr)
library(grid)
library(broom)
library(tidyr)
library(scales)
library(gridExtra)
library(data.table)
library(glmnet)

source('helper_functions.R')

# source functions script
# source('all_functions.R')

# create plotting functions


# create fixed objects to model and pipeline inputs and saving  
methyl_type = 'beta'
combat = TRUE
gender = TRUE
tech = FALSE 
how_many_seeds = 100
how_many_folds = 5


# condition on fixed objects to get saving identifiers

if(methyl_type == 'beta'){
  which_methyl <- 'beta'
} else {
  which_methyl <- 'm'
}

if(combat){
  which_combat <- 'use_combat'
  beta_thresh = 0.05
} else {
  which_combat <- 'no_combat'
  beta_thresh = 0.5
}

if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}

num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds

# read in cases_450
cases_dat <- readRDS(paste0('validation_age_predictions/', 'cases_450',which_methyl, '_', which_combat, '_',
                          num_seeds, '_', num_folds, '_', is_gen, '.rda'))
# read in validation set 
valid_dat <- readRDS(paste0('validation_age_predictions/', 'valid_',which_methyl, '_', which_combat, '_',
                            num_seeds, '_', num_folds, '_', is_gen, '.rda'))

# save validation set 
con_dat <- readRDS(paste0('validation_age_predictions/', 'controls_',which_methyl, '_', which_combat, '_',
                          num_seeds, '_', num_folds, '_', is_gen, '.rda'))

# save associated lfs bumps
lfs_bump_probes <- readRDS(paste0('validation_age_predictions/', 'lfs_bumps_', which_methyl, '_', 
                                  which_combat, '.rda'))
# cases <- cases_dat
# controls <- controls_dat
# valid <- valid_dat
# age_cutoff <- 72
# gender = TRUE
# tech = FALSE
# train_lambda = TRUE
# alpha_value <- 0.8
# lambda_value <- 0.18
# bh_features <- lfs_bump_probes
test_model <- function(cases, 
                       controls, 
                       valid, 
                       test_controls,
                       train_lambda,
                       gender,
                       tech,
                       age_cutoff,
                       alpha_value, 
                       lambda_value,
                       bh_features) {
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(cases))
  
  if(gender) {
    intersected_feats <- c('Female', 'Male', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('batch_1', 'batch_2', intersected_feats)
  }
  
  
  
  # get outcomes
  cases_y <- as.factor(ifelse(cases$age_diagnosis < age_cutoff, 'positive', 'negative'))
  valid_y <- as.factor(ifelse(valid$age_diagnosis < age_cutoff, 'positive', 'negative')) 
  controls_y <-  as.factor(ifelse(controls$age_sample_collection < age_cutoff, 'positive', 'negative'))
  controls_clin <- controls[, !grepl('^cg', colnames(controls))]
  valid_clin <- valid[, !grepl('^cg', names(valid))]
  
  # get model data
  cases <- cases[, intersected_feats]
  controls <- controls[, intersected_feats]
  valid <- valid[, intersected_feats]
  
  # store fixed values
  best_alpha <- alpha_value
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  nfolds = 5
  
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(cases)
                                     , y =  cases_y
                                     , alpha = alpha_value
                                     , type.measure = type_measure
                                     , family = type_family
                                     , standardize=FALSE
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    
    
    # get outcome variables and clin variables
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    trained_lambda_index_set = elastic_net.cv_model$lambda[elastic_net.cv_model$lambda < 0.26 & elastic_net.cv_model$lambda > 0.24]
    trained_lambda_value <- trained_lambda_index_set[1]
    trained_lambda_index <- which(elastic_net.cv_model$lambda == trained_lambda_value) 
    
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    temp.loop_count = temp.loop_count + 1
    
    # set seed for next loop iteration
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    if (temp.loop_count > 10) {
      print("diverged")
      temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
      break
    }
  }# while loop ends 
  # print(temp.non_zero_coeff)  
  
  model  = glmnet(x = as.matrix(cases)
                  , y =  cases_y
                  ,alpha = best_alpha
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  
  # Predictions on controls data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_con <- predict(model, 
                                       data.matrix(controls),
                                       type = 'response')
  
  if(train_lambda){
    # get predictions with corresponding lambda.
    test.predictions_con <- temp_test.predictions_con[, trained_lambda_index]
    
  } else {
    # get predictions with corresponding lambda.
    test.predictions_con <- temp_test.predictions_con[, temp.min_lambda_index]
    
  }
  
  # combine predictions and real labels 
  temp_dat_con <- as.data.frame(cbind(controls_age_pred = test.predictions_con, controls_age_label = controls_y, controls_clin))
  temp_dat_con$alpha <- best_alpha
  
  
  # Predictions on validation data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_valid <- predict(model, 
                                         data.matrix(valid),
                                         type = 'response')
  
  if(train_lambda){
    # get predictions with corresponding lambda.
    test.predictions_valid <- temp_test.predictions_valid[, trained_lambda_index]
    
  } else {
    # get predictions with corresponding lambda.
    test.predictions_valid <- temp_test.predictions_valid[, temp.min_lambda_index]
    
  }
  
  # combine predictions and real labels 
  temp_dat_valid <- as.data.frame(cbind(valid_age_pred = test.predictions_valid, vallid_age_label = valid_y, valid_clin))
  temp_dat_valid$alpha <- best_alpha
  
  if(test_controls){
    return(temp_dat_con)
  } else{
    return(temp_dat_valid)
  }
  
  
}

# creat list to store results for alpha
result_list <- list()

alpha_values <- (1:10/10)

for(i in 1:length(alpha_values)){
  alpha_num <- alpha_values[i]
  
  message('working on alpha = ', alpha_num)
  result_list[[i]] <- test_model(cases = cases_dat,
                                 controls = con_dat,
                                 valid = valid_dat,
                                 test_controls = FALSE,
                                 age_cutoff = 72,
                                 gender = gender,
                                 tech = tech,
                                 train_lambda = FALSE,
                                 alpha_value = alpha_num,
                                 lambda_value = 0.15,
                                 bh_features = lfs_bump_probes)
}


temp <- do.call('rbind', result_list)


##########
# validation set
##########

# read in cases_450
saveRDS(temp, paste0('validation_age_predictions/', 'valid_test_untrained',alpha_num,'_' ,which_methyl, '_', which_combat, '_',
                            num_seeds, '_', num_folds, '_', is_gen, '.rda'))

