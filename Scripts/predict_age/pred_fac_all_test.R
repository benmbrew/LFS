##########
# initialize libraries
##########
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(ModelMetrics)
library(doParallel)
library(nnet)
library(dplyr)
library(bumphunter)
library(sqldf)
library(e1071)
library(reshape2)

registerDoParallel(1)

##########
# initialize folders
##########

home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
feat_data <- paste0(data_folder, '/feat_data')
results_data <- paste0(data_folder, '/results_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'funnorm'
k = 5
combat = T

##########
# load data
##########


if(combat) {
  
  mod_data <-  readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_combat.rda')))
  
  
} else {
  mod_data <- readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data.rda')))
  
}


##########
# get column names
##########
intersect_names <- colnames(betaFull)[9:ncol(betaFull)]

##########
# read in all features from feat data
##########
setwd(feat_data)
file_list_names = list.files()
file_list_names = file_list_names[1:5]

# store all raw rda feature lists in feat_list
feat_list = lapply(file_list_names, function(x) readRDS(x))

# get first 6 only 
model_names <- c('enet', 'rf', 'lasso')



l = j = m = 1
#
beta_full = betaFull
model_method = model_names[m]
mod_feats = feat_list[[l]]
feat_name = file_list_names[[l]]
k = k
class_age = 72
max_age = 850
data_type = 'standard'

#

trainTest <- function(beta_full,
                      model_method,
                      mod_feats,
                      feat_name,
                      data_type,
                      max_age,
                      remove_cases,
                      class_age,
                      k) {
  
  
  if(data_type == 'combat') {
    
    # load in data
    beta_cases <- beta_full[beta_full$type == 'cases',]
    beta_controls <- beta_full[beta_full$type == 'controls',]
    beta_valid <- beta_full[beta_full$type == 'valid',]
    
    # get mod data 
    beta_cases <- getModData(beta_cases)
    
    # get rid of cancer samples in controls 
    beta_controls <- beta_controls[grepl('Unaffected', beta_controls$cancer_diagnosis_diagnoses),]
    
    #subset valid
    beta_valid <- beta_valid[!beta_valid$ids %in% beta_cases$ids,]
    
  } else if(data_type == 'standard')  {
    
    beta_cases <- beta_full[beta_full$type == 'cases_450k',]
    beta_controls <- beta_full[beta_full$type == 'controls_850k',]
    beta_valid <- beta_full[beta_full$type == 'valid_850k',]
    
  } else if (data_type == 'full') {
    
    beta_cases <- beta_full[grepl('cases|valid', beta_full$type),]
    beta_controls <- beta_full[grepl('controls', beta_full$type),]
    beta_valid <- beta_full[beta_full$type == 'valid_850k',]
    
  } else if (data_type == 'old_controls') {
    
    beta_cases <- beta_full[beta_full$type == 'cases_450k',]
    beta_controls <- beta_full[beta_full$type == 'controls_850k',]
    beta_valid <- beta_full[beta_full$type == 'valid_850k',]
    
  }
  
  # remove samples that dont have an age of sample collection
  beta_cases <- beta_cases[complete.cases(beta_cases),]
  beta_valid <- beta_valid[complete.cases(beta_valid),]
  
  beta_controls <- beta_controls[!is.na(beta_controls$age_sample_collection),]
  beta_cases <- beta_cases[!is.na(beta_cases$age_sample_collection),]
  
  
  # if(remove_cases) {
  #   length_index <- 1:nrow(beta_cases)
  #   remove_index <- length_index[beta_cases$age_sample_collection >= 12 & beta_cases$age_sample_collection <= 24]
  #   remove_index <- sample(remove_index, 10, replace = F)
  #   
  #   beta_cases <- beta_cases[-remove_index,]
  #   
  # }
  
  beta_cases <- beta_cases[beta_cases$age_sample_collection <=max_age,]
  
  
  # list to store results
  temp_results <- list()
  
  # sample random 10 from each cluster 
  column_names <- colnames(beta_full)[9:ncol(beta_full)]
  sample_cols <- column_names[!column_names %in% mod_feats]
  rand_feats <- sample(sample_cols, length(mod_feats), replace = T)

  if(model_method == 'enet') {
    
    mod_result <- runEnetRandFac(training_dat = beta_cases, 
                                 controls_dat = beta_controls,
                                 valid_dat = beta_valid,
                                 test_dat = beta_cases, 
                                 age_cutoff = class_age,
                                 bh_features = mod_feats,
                                 rand_feats = rand_feats,
                                 gender = F)
  } 
  
  if(model_method == 'rf') {
    # # get residuals
    # cases_resid <- getResidual(data = cases, 
    #                            bh_features = bh_feat_sig)
    
    mod_result <- runRfRandFac(training_dat = beta_cases, 
                               controls_dat = beta_controls,
                               valid_dat = beta_valid,
                               test_dat = beta_cases, 
                               age_cutoff = class_age,
                               bh_features = mod_feats,
                               rand_feats = rand_feats,
                               pred_cutoff = .5,
                               gender = F)
    
    
  } 
  
  if(model_method == 'svm') {
    # # get residuals
    # cases_resid <- getResidual(data = cases, 
    #                            bh_features = bh_feat_sig)
    
    
    mod_result <- runSvmRandFac(training_dat = beta_cases[train_index,], 
                                controls_dat = beta_controls,
                                valid_dat = beta_valid,
                                test_dat = beta_cases[test_index,], 
                                age_cutoff = age_cutoff,
                                bh_features = mod_feats,
                                rand_feats = rand_feats,
                                gender = F)
    
  } 
  
  if(model_method == 'lasso') {
    mod_result <- runLassoL1RandFac(training_dat = beta_cases, 
                                    controls_dat = beta_controls,
                                    valid_dat = beta_valid,
                                    test_dat = beta_cases, 
                                    age_cutoff = class_age,
                                    bh_features = mod_feats,
                                    rand_feats = rand_feats,
                                    gender = F)
    
  } 
  
  temp_results <- mod_result
  
  return(temp_results)
  
}

temp_results_class <- list()
temp_results_mat <- list()
temp_results_class_2 <- list()
temp_results_mat_2 <- list()
full_results_class <- list()
full_results_mat <- list()
# set fixed variables 

for(m in 1:length(model_names)) {
  
  for(l in 1:length(feat_list)) {
    

      # return list of 2
      mod_results <- trainTest(beta_full = betaFull,
                               model_method = model_names[m],
                               mod_feats = feat_list[[l]],
                               feat_name = file_list_names[[l]],
                               data_type = 'combat',
                               max_age = 1000,
                               remove_cases = F,
                               class_age = 72,
                               k)
      
      
      # this is 5 observations (folds for each of the 7 age variables)
      # cases onset
      temp_results_class[[l]] <- get_class_results_test(mod_results, dims_of_dat = length(feat_list[[l]]), 
                                                        mod_name = model_names[m], feat_name = file_list_names[[l]])[[1]]
      temp_results_mat[[l]] <- get_class_results_test(mod_results, dims_of_dat = length(feat_list[[l]]), 
                                                      mod_name = model_names[m], feat_name = file_list_names[[l]])[[2]]
      
      # five folds, then first element is confusion matrix, and secod element is clss scores
      
      print(paste0('done with ', l, ' featues'))
    }
    temp_results_class_2[[m]] <- do.call(rbind, temp_results_class)
    temp_results_mat_2[[m]] <- temp_results_mat
    print(paste0('done with ', m, ' models'))
    
    
}

full_results <- do.call(rbind, temp_results_class_2)
full_results_mat <- temp_results_mat_2

# options(scipen=999)
# save.image('~/Desktop/result_class_temp_217.RData')



saveRDS(full_results, paste0(results_data, '/', method,'_','full_results_class_test_850_450_48.rda'))


colnames(full_results) <- tolower(colnames(full_results))
colnames(full_results) <- gsub(' ', replacement = '_', colnames(full_results))

# group by age outcome type, number of features in model, and the model name and get 
# the mean "sensitivity (true positive, recall)", "specificity (true negative)", "precision", "balanced_accuracy
# FNR = 1 - TPRs
# FPR = 1 - TNR
# TRP = TP/(TP + FN)
# Precision = TP/(TP + FP)

# Precision is 
# (TP)/(TP+FP)
# which tells us what proportion of patients we diagnosed as having cancer actually had cancer. 
# In other words, proportion of TP in the set of positive cancer diagnoses. This is given by the rightmost 
# column in the confusion matrix.

# Recall is
# (TP)/(TP+FN)
# which tells us what proportion of patients that actually had cancer were diagnosed by us as having cancer. 
# In other words, proportion of TP in the set of true cancer states. This is given by the bottom row in 
# the confusion matrix.

# Recall (TPR, sensitvity) is the probability that a (randomly selected) relevant document is retrieved in a search.
# Precision is the probability that a (randomly selected) retrieved document is relevant.

# In this representation, it is clearer that recall gives us information about a classifiers 
# performance with respect to false negatives (how many did we miss), while precision gives us 
# information about its performance with respect to false positives.
full_results <- as.data.frame(full_results)
temp <- 
  full_results %>%
  group_by(age_type, feature_num,feat_name, model_method) %>%
  summarise(mean_acc = mean(balanced_accuracy, na.rm = T),
            mean_prec = mean(precision, na.rm = T),
            mean_tpr = mean(sensitivity, na.rm = T),
            mean_tnr = mean(specificity, na.rm = T))

# remove svm and make fnr and fpr 
temp <- temp[!grepl('svm', temp$model_method),]
temp$mean_fnr <- 1 - temp$mean_tpr
temp$mean_fpr <- 1 - temp$mean_tnr

temp <- as.data.frame(temp)

temp <- temp[order(temp$feature_num, decreasing = F),]

saveRDS(temp, paste0(results_data, '/', method,'_','full_results_class_collapsed.rda'))



