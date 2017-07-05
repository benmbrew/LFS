#########
# This script will train model on full data and test on controls and valid 
##########
# initialize libraries
##########
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(doParallel)
library(nnet)
library(dplyr)
library(bumphunter)
library(sqldf)

registerDoParallel(1)
##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'quan'
k = 4
seed_num = 1


##########
# load data
##########
betaCases <- readRDS(paste0(model_data, '/betaCases', method,'.rda'))
betaControls <- readRDS(paste0(model_data, '/betaControls', method,'.rda'))
betaControlsOld <- readRDS(paste0(model_data, '/betaControlsOld', method,'.rda'))
betaValid <- readRDS(paste0(model_data, '/betaValid', method,'.rda'))
# 
# # # TEMP
# betaCases <- betaCases[!grepl('9721365183', betaCases$sentrix_id),]

# load cg_locations
cg_locations <- read.csv(paste0(model_data, 
                                '/cg_locations.csv'))

cg_locations$X <- NULL
##########
# get intersecting colnames and prepare data for modeling
##########

intersect_names <- Reduce(intersect, list(colnames(betaCases)[8:ncol(betaCases)], 
                                          colnames(betaControls)[8:ncol(betaControls)], 
                                          colnames(betaValid)[8:ncol(betaValid)]))
# organize each data set accordling

# cases
betaCases <- betaCases[, c('age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           intersect_names)]
# controls
betaControls <- betaControls[, c('age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender', 
                                 intersect_names)]

# controls
betaControlsOld <- betaControlsOld[, c('age_diagnosis', 
                                       'age_sample_collection', 
                                       'cancer_diagnosis_diagnoses', 
                                       'gender', 
                                       intersect_names)]
#validation
betaValid <- betaValid[, c('age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           intersect_names)]

betaControlsFull <- rbind(betaControls, betaControlsOld)

# remove na in sample collection and duplicate ids
betaControlsFull <- betaControlsFull[!is.na(betaControlsFull$age_sample_collection),]

rm(betaControlsOld)

##########
# test for distribution sameness 
########## 
testKS(x = betaCases$age_sample_collection, y = betaControls$age_sample_collection)
testKS(x = betaCases$age_sample_collection, y = betaControlsFull$age_sample_collection)


##########
# run bumphunter
##########
bh_feat <- bumpHunterSurv(dat_cases = betaCases, dat_controls = betaControls)
# bh_feat_full <- bumpHunterSurv(dat_cases = betaCases, dat_controls = betaControlsFull)


##########
# get features
##########
bh_feat_all <- getProbe(bh_feat)
bh_feat_tot <- getRun(bh_feat_all[[1]], run_num = .15)
bh_feat_sig <- getRun(bh_feat_all[[2]], run_num = .15)


# saveRDS(bh_feat_all, paste0(model_data, '/bh_feat_all.rda'))
# saveRDS(bh_feat_all, paste0(model_data, '/bh_feat_all_quan.rda'))
bh_feat_all <- readRDS(paste0(model_data, '/bh_feat_all_quan.rda'))

# #raw
# save.image('/home/benbrew/Desktop/temp_full_test.RData')
# load('/home/benbrew/Desktop/temp_full_test.RData')
# #quan
# save.image('/home/benbrew/Desktop/temp_full_test_quan.RData')
# load('/home/benbrew/Desktop/temp_full_test_quan.RData')


##########
# test model
# ##########
# cases_dat <- betaCases
# controls_dat <- betaControls
# controls_dat_full <- betaControlsFull
# valid_dat <- betaValid
# bh_features <- bh_feat_sig
# alpha = 0.9
testModel <- function(cases_dat,
                      controls_dat,
                      controls_dat_full,
                      valid_dat,
                      bh_features,
                      alpha)
{
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(cases_dat))
  
  intersected_feats <- append('gender', intersected_feats)
  cases_dat$gender <- as.numeric(as.factor(cases_dat$gender))
  controls_dat$gender <- as.numeric(as.factor(controls_dat$gender))
  controls_dat_full$gender <- as.numeric(as.factor(controls_dat_full$gender))
  valid_dat$gender <- as.numeric(as.factor(valid_dat$gender))
  
  # get y
  train_y <- as.numeric(cases_dat$age_diagnosis)
  test_y_controls <- as.numeric(controls_dat$age_sample_collection)
  test_y_controls_full <- as.numeric(controls_dat_full$age_sample_collection)
  test_y_valid <- as.numeric(valid_dat$age_diagnosis)
  
  # get bumphunter features
  cases_dat <- cases_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  controls_dat_full <- controls_dat_full[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
  N_CV_REPEATS = 2
  nfolds = 3
  
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  elastic_net.cv_model = cv.glmnet(x = as.matrix(cases_dat)
                                   , y =  train_y
                                   , alpha = alpha
                                   , type.measure = 'deviance'
                                   , family = 'gaussian'
                                   , standardize=FALSE
                                   , nlambda = 100
                                   , nfolds = nfolds
                                   , parallel = TRUE
  )
    
    
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    print(temp.non_zero_coeff)  
  
  
    model  = glmnet(x = as.matrix(cases_dat)
                  , y =  train_y
                  ,alpha = alpha
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = 'gaussian')
  
  # get controls
  temp_test.predictions_controls <- predict(model,
                                            data.matrix(controls_dat),
                                            type = 'response')



  test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
  
  # get controls full
  temp_test.predictions_controls_full <- predict(model,
                                                data.matrix(controls_dat_full),
                                                type = 'response')
  
  
  
  test.predictions_controls_full <- temp_test.predictions_controls_full[, temp.min_lambda_index]
  


  # get validation
  temp_test.predictions_valid <- predict(model,
                                         data.matrix(valid_dat),
                                         type = 'response')



  test.predictions_valid  <- temp_test.predictions_valid[, temp.min_lambda_index]

  importance <- coef(model)
  
  lambda_value <- elastic_net.cv_model$lambda.min
  
  # # for each iteration, this should always be the same.
  controls_cor <- cor(test_y_controls, test.predictions_controls)
  
  controls_cor_full <- cor(test_y_controls_full, test.predictions_controls_full)
  
  valid_cor  <- cor(test_y_valid, test.predictions_valid)

  return(list(controls_cor, controls_cor_full, valid_cor))
  
  
}


##########
# calssification
##########

testModelFac <- function(cases_dat,
                         controls_dat,
                         controls_dat_full,
                         valid_dat,
                         bh_features,
                         cutoff,
                         alpha)
{
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(cases_dat))
  
  intersected_feats <- append('gender', intersected_feats)
  cases_dat$gender <- as.numeric(as.factor(cases_dat$gender))
  controls_dat$gender <- as.numeric(as.factor(controls_dat$gender))
  controls_dat_full$gender <- as.numeric(as.factor(controls_dat_full$gender))
  valid_dat$gender <- as.numeric(as.factor(valid_dat$gender))
  
  # # get y
  train_y <- factor(ifelse(cases_dat$age_diagnosis <= cutoff, 'yes', 'no'), 
                    levels = c('yes', 'no'))
  test_y_controls <- factor(ifelse(controls_dat$age_sample_collection <= cutoff, 'yes', 'no'), 
                            levels = c('yes', 'no'))
  test_y_controls_full <- factor(ifelse(controls_dat_full$age_sample_collection <= cutoff, 'yes', 'no'), 
                                 levels = c('yes', 'no'))
  test_y_valid <- factor(ifelse(valid_dat$age_sample_collection <= cutoff, 'yes', 'no'), 
                         levels = c('yes', 'no'))
  
  # get bumphunter features
  cases_dat <- cases_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  controls_dat_full <- controls_dat_full[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
  N_CV_REPEATS = 2
  nfolds = 3
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  elastic_net.cv_model = cv.glmnet(x = as.matrix(cases_dat)
                                   , y =  train_y
                                   , alpha = alpha
                                   , type.measure = 'auc'
                                   , family = type_family
                                   , standardize=FALSE
                                   , nlambda = 100
                                   , nfolds = nfolds
                                   , parallel = TRUE
  )
  
  
  temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
  
  # # number of non zero coefficients at that lambda    
  temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
  print(temp.non_zero_coeff)  
  
  
  model  = glmnet(x = as.matrix(cases_dat)
                  , y =  train_y
                  ,alpha = alpha
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  # get controls
  temp_test.predictions_controls <- predict(model,
                                            data.matrix(controls_dat),
                                            type = 'class')
  
  
  
  test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
  
  test.predictions_controls <- factor(test.predictions_controls, levels = c('yes', 'no'))
  
  
  test_stats_controls <- confusionMatrix(test_y_controls, test.predictions_controls)
  
  
  # get controls full
  temp_test.predictions_controls_full <- predict(model,
                                                 data.matrix(controls_dat_full),
                                                 type = 'class')
  
  
  
  test.predictions_controls_full <- temp_test.predictions_controls_full[, temp.min_lambda_index]
  
  test.predictions_controls_full <- factor(test.predictions_controls_full, levels = c('yes', 'no'))
  
  
  test_stats_controls_full <- confusionMatrix(test_y_controls_full, test.predictions_controls_full)
  
  
  # get validation
  temp_test.predictions_valid <- predict(model,
                                         data.matrix(valid_dat),
                                         type = 'class')
  
  
  
  test.predictions_valid  <- temp_test.predictions_valid[, temp.min_lambda_index]
  
  test.predictions_valid <- factor(test.predictions_valid, levels = c('yes', 'no'))
  
  
  test_stats_valid <- confusionMatrix(test_y_valid, test.predictions_valid)
  
  
  importance <- coef(model)
  
  lambda_value <- elastic_net.cv_model$lambda.min
  
  
  
  return(list(test_stats_controls, test_stats_controls_full, test_stats_valid))
  
  
}



testModel(betaCases,
          betaControls,
          betaControls,
          betaValid,
          bh_feat_tot,
          alpha = 0.9)



testModelFac(betaCases,
             betaControls,
             betaControlsFull,
             betaValid,
             bh_feat_sig,
             cutoff = 60,
             alpha = 0.9)
