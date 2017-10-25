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
method = 'noob'
k = 5
combat = T

##########
# load data
##########


if(combat) {
  
  betaFull <-  readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_combat_m.rda')))
  
  
} else {
  
  betaFull <- readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_m.rda')))
  
}

##########
# get column names
##########
intersect_names <- colnames(betaFull)[9:ncol(betaFull)]



##########
# read in all features from feat data
##########

lfs_feats_m <- readRDS(paste0(feat_data, paste0('/', method, '_', 'lfs_m.rda')))
no_cancer_feats_m <- readRDS(paste0(feat_data, paste0('/', method, '_', 'no_cancer_m.rda')))


setwd(feat_data)
file_list_names = list.files()


# store all raw rda feature lists in feat_list
feat_list <- lapply(file_list_names, function(x) readRDS(x))

# order feat list
feat_list  <- feat_list[order(sapply(feat_list, length), decreasing=F)]

# select first 10
feat_list <- feat_list[10]
# 
# # model_names <- c('enet', 'rf', 'lasso')
# # seeds <- c(1, 2, 3)

model_names <- c('enet')
file_list_names <- 'first'



l = j = m = 2
#
beta_full = betaFull
model_method = model_names[m]
mod_feats = feat_list[[l]]
feat_name = file_list_names[[l]]
k = k
class_age = 72
max_age = 850
data_type = 'standard'
temp_feat <- feat_list[[2]]
# temp <- betaFull[betaFull$batch == 'controls',]
# temp <- temp[, c('ids','gender','family_name','age_sample_collection',temp_feat)]
# #
# final <- cbind(test.predictions_controls, temp)
# final <- final[, 1:5]
# final$real <- ifelse(final$age_sample_collection < 72, 'yes', 'no')
# final$missed <- ifelse(final$test.predictions_controls == final$real, 'good', 'bad')
# final <- final[final$missed == 'bad',]
# test.predictions_controls
# temp$ids



trainTest <- function(beta_full,
                      model_method,
                      mod_feats,
                      feat_name,
                      data_type,
                      max_age,
                      remove_cases,
                      class_age,
                      k) {
  
  
  beta_full <- beta_full[beta_full$age_sample_collection < max_age,]
  if(data_type == 'standard') {
    
    beta_cases <- beta_full[beta_full$batch == 'cases',]
    beta_controls <- beta_full[beta_full$batch == 'controls',]
    beta_valid <- beta_full[beta_full$batch == 'valid',]
    
  } else {
    
    beta_cases <- beta_full[grepl('cases|valid', beta_full$batch),]
    beta_controls <- beta_full[grepl('controls', beta_full$batch),]
    beta_valid <- beta_full[beta_full$batch == 'valid',]
    
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
    
    mod_result <- runEnetRand(training_dat = beta_cases, 
                                 controls_dat = beta_controls,
                                 valid_dat = beta_valid,
                                 test_dat = beta_cases, 
                                 age_cutoff = class_age,
                                 bh_features = mod_feats,
                                 gender = T)
  } 
  
  if(model_method == 'rf') {
    # # get residuals
    # cases_resid <- getResidual(data = cases, 
    #                            bh_features = bh_feat_sig)
    
    mod_result <- runRfRand(training_dat = beta_cases, 
                               controls_dat = beta_controls,
                               valid_dat = beta_valid,
                               test_dat = beta_cases, 
                               age_cutoff = class_age,
                               bh_features = mod_feats,
                               rand_feats = rand_feats,
                               pred_cutoff = .5,
                               gender = F)
    
    
    
  } 
  
  
  
  temp_results <- mod_result
  
  return(temp_results)
  
}


results <- trainTest(beta_full = betaFull, 
                     model_method = 'enet', 
                     mod_feats = feat_list[[1]], 
                     feat_name = file_list_names[[1]], 
                     data_type = 'standard', 
                     max_age = 1000, 
                     remove_cases = F, 
                     class_age = 48) 


real_age <- results[[2]]
pred_age <- results[[1]]

#here
pred_age <- pred_age +80

temp_final <- as.data.frame(cbind(real_age, pred_age))

plot(temp_final$real_age, temp_final$pred_age, xlim = c(0, 800), ylim = c(0, 800))

ggplot(temp_final, aes(real_age, pred_age)) + 
  geom_point(size = 3, colour = 'darkblue', alpha = 0.7) + 
  xlim(c(0,800)) + ylim(c(0,800)) + 
  xlab('Real Age') + ylab('Predicted Age') +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + theme(axis.text.y = element_text(size = 12),
                     axis.text.x = element_text(size = 12)) 


cor(real_age, pred_age)
abline(0 ,1)

temp_final$ages <- ifelse(temp_final$pred_age > temp_final$real_age, 'good', 'bad')


# 1 test.predictions, 
# 2 test_y, 
# 3 patient_age, 
# 4 test.predictions_valid 
# 5 valid_y 
# 6 patient_age_valid
# 7 test.predictions_controls 
# 8 patient_age_controls
# 
