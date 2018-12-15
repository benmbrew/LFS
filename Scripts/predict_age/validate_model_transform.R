
# load libraries
library(plotly)
library(scatterplot3d)
library(preprocessCore)
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

# register other cpus
registerDoParallel(2)

# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  
methyl_type = 'm'
gender = TRUE
tech = FALSE 
how_many_seeds = 20
how_many_folds = 5
# use_offset = FALSE
remove_age  = TRUE
beta_thresh = 0.05
age_cutoff = 72


# condition on fixed objects to get saving identifiers

if(methyl_type == 'beta'){
  which_methyl <- 'beta'
} else {
  which_methyl <- 'm'
  beta_thresh= 0.5
}

if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}

# if(use_offset){
#   is_offset <- 'use_offset'
# } else {
#   is_offset <- 'no_offset'
# }

num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds

# read trainig  set 
cases_450 <- readRDS(paste0('transform_data/', 'cases_450_',which_methyl, '.rda'))
# read validation set 
con_transform <- readRDS(paste0('transform_data/', 'con_transform_',which_methyl, '.rda'))

# read validation set 
valid_transform <- readRDS(paste0('transform_data/', 'valid_transform_',which_methyl,'.rda'))
# read wt controls 
con_wt <- readRDS(paste0('transform_data/', 'con_wt_',which_methyl, '.rda'))
# read mut controls 
con_mut <- readRDS(paste0('transform_data/', 'con_mut_',which_methyl,'.rda'))

# load model_params
model_params <- readRDS(paste0('transform_data_test/', 'model_params_',which_methyl,'_',
                               num_seeds, '_',k_folds, '_', is_gen ,'.rda'))

# load optimal_cutoff
optimal_cutoff <- readRDS(paste0('transform_data_test/', 'optimal_cutoff_',which_methyl,'_',
                                 num_seeds, '_',k_folds, '_', is_gen ,'.rda'))

# read associated lfs bumps
lfs_bump_probes <- readRDS(paste0('transform_data/', 'lfs_bumps_', which_methyl, '_', '.rda'))
############################
HERE: make sure to load in model params and optimal_cutoff and implement them correctl in the model
# get intersecting features for con_wt

test_model <- function(cases, 
                       controls, 
                       valid, 
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
    intersected_feats <- c('450k', '850k', intersected_feats)
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

  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  nfolds = 5
  
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  elastic_net.cv_model = cv.glmnet(x = as.matrix(cases)
                                   , y =  cases_y
                                   , alpha = alpha_value
                                   , type.measure = type_measure
                                   , family = type_family
                                   , standardize=FALSE
                                   , nlambda = 100
                                   , nfolds = nfolds
                                   , parallel = TRUE)
  
    
  
    # get outcome variables and clin variables
    lambda_s <- elastic_net.cv_model$lambda.min
    lambda_s_train <- lambda_value
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 

  # print(temp.non_zero_coeff)  
  
  model  = glmnet(x = as.matrix(cases)
                  , y =  cases_y
                  ,alpha = best_alpha
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  
  # Predictions on controls data
  
  
  
  if(cv_lambda){
    # This returns 100 prediction with 1-100 lambdas
    test.predictions_con <- predict.glmnet(model, 
                                                data.matrix(controls),
                                                type = 'response',
                                                s = lambda_s_train)
    
  } else {
    # This returns 100 prediction with 1-100 lambdas
    temp_test.predictions_con <- predict(model, 
                                           data.matrix(controls),
                                           type = 'response')
    # get predictions with corresponding lambda.
    test.predictions_con <- temp_test.predictions_con[, temp.min_lambda_index]
    
  }
  
  # combine predictions and real labels 
  temp_dat_con <- as.data.frame(cbind(controls_age_pred = test.predictions_con, controls_age_label = controls_y, controls_clin))
  temp_dat_con$alpha <- best_alpha
  temp_dat_con$non_zero <- temp.non_zero_coeff
  
  if(cv_lambda){
    # This returns 100 prediction with 1-100 lambdas
    test.predictions_valid <- predict.glmnet(model, 
                                           data.matrix(valid),
                                           type = 'response',
                                           s = lambda_s_train)
    
  } else {
    # This returns 100 prediction with 1-100 lambdas
    temp_test.predictions_valid <- predict(model, 
                                         data.matrix(valid),
                                         type = 'response')
    # get predictions with corresponding lambda.
    test.predictions_valid <- temp_test.predictions_valid[, temp.min_lambda_index]
    
  }
  
  # combine predictions and real labels 
  temp_dat_valid <- as.data.frame(cbind(valid_age_pred = test.predictions_valid, valid_age_label = valid_y, valid_clin))
  temp_dat_valid$alpha <- best_alpha
  temp_dat_valid$non_zero <- temp.non_zero_coeff
  
  
  
  return(list(temp_dat_con, temp_dat_valid))
  
  
}

# creat list to store results for alpha
result_list <- list()
con_list <- list()
valid_list <- list()

alpha_values <- (1:10/10)

for(i in 1:length(alpha_values)){
  alpha_num <- alpha_values[i]
  
  message('working on alpha = ', alpha_num)
  result_list[[i]] <- test_model(cases = cases_450,
                                 controls = con_transform,
                                 valid = valid_transform,
                                 age_cutoff = 72,
                                 gender = gender,
                                 tech = tech,
                                 train_lambda = FALSE,
                                 alpha_value = alpha_num,
                                 lambda_value = s_value,
                                 bh_features = lfs_bump_probes)
  
  con_list[[i]] <- result_list[[i]][[1]]
  valid_list[[i]] <- result_list[[i]][[2]]
  
}


temp_con <- do.call('rbind', con_list)
temp_valid <- do.call('rbind', valid_list)


##########
# validation set
##########

# read in cases_450
saveRDS(temp_con, paste0('transform_data_test/', 'con_test_transform',which_methyl, '_', is_gen, '.rda'))

saveRDS(temp_valid, paste0('transform_data_test/', 'valid_test_transform',which_methyl, '_', is_gen, '.rda'))
