####### Script will run models

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

registerDoParallel(1)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')
results_folder <- paste0(project_folder, '/Results')
quan_folder <- paste0(results_folder, '/quan')
rand_folder <- paste0(results_folder, '/rand')
scripts_folder <- paste0(project_folder, '/Scripts')


##########
# load data batch data
##########
# read cases
quan_controls <- readRDS(paste0(model_data, '/quan_controls.rda'))

quan_controls_gen <- readRDS(paste0(model_data, '/quan_controls_gen.rda'))


# load features
load(paste0(model_data, '/bh_feat.RData'))

##########
# source model_functions to get functions to run models 
##########
source(paste0(scripts_folder, '/predict_age/functions.R'))

# ##########
# # set parameters 
# ##########
# data_thresholds <- c(48, 60, 72)
# p53 <- c('Mut', 'WT')


##########
# function to run models - subset, get residuals, get categorical, predict with regression and fac. 
##########

runModels <- function(data,
                      model,
                      random = F,
                      bump_hunter = F,
                      bump_hunter_data,
                      num_feat = NULL,
                      gender,
                      residual,
                      seed_num) 
{
  
  # get differenct variations of data
  if (random) {
    
    data <- subsetDat(data, 
                      random = T, 
                      num_feat,
                      seed_num)
  } else {
    
    data <- subsetDat(data, 
                      random = F, 
                      gender,
                      num_feat = NULL)
  }
  
  
  #########
  # get regression data and run regression
  #########
  # get data
  if (bump_hunter) {
    
    data <- bhSubset(data, bh_data = bump_hunter_data, gender)
    
  }
  
  if (residual) {
    data <- getResidual(data, gender)
  }
  
  
  
  if (model == 'rf') {
    
    data_result <- rfPredictReg(data, cutoff = .7, iterations = 10, control = T)
    # data_resid_result <- rfPredictReg(data_resid, cutoff = .7, iterations = 5)
  }
  
  if (model == 'enet') {
    
    data_result <- enetPredReg(data, N_CV_REPEATS = 2, nfolds = 5,cutoff = .7,gender = T ,iterations = 10, control =T)
    # data_resid_result <- enetPredReg(data_resid, N_CV_REPEATS = 2, nfolds = 5,cutoff = .7, iterations = 5)
    
  }
  
  # if (model == 'lasso') {
  #   
  #   data_result <- regPred(data, alpha = 1, nfolds = 5, cutoff = .7, iterations = 5)
  #   # data_resid_result <- regPred(data_resid, alpha = 1, nfolds = 5, cutoff = .7, iterations = 5)
  #   
  # }
  # 
  # #########
  # # get classification data and run classifier
  # #########
  # data_fac <- makeFac(data, threshold = 48)
  # data_resid_fac <- makeFac(data_resid, threshold = 48)
  # 
  # data_fac_result  <- rfPredictFac(data_fac, cutoff = .7, iterations = 5)
  # data_resid_fac_result <- rfPredictFac(data_resid_fac, cutoff = .7, iterations = 5)
  
  return (data_result) 
  # data_fac_result, 
  # data_resid_fac_result))
  
  
}

########################################################################

# 
# # full data sam
# even_full_sam_rf <- runModels(quan_cases_sam, 
#                               model = 'rf',
#                               rand = F,
#                               bump_hunter = F,
#                               bump_hunter_data)
# 
# even_full_sam_table_rf <- extractResults(even_full_sam_rf, data_name = 'even sam_full rf', regularize = F)
# 
# even_full_sam_enet <- runModels(quan_cases_sam, 
#                                 model = 'enet',
#                                 rand = F,
#                                 bump_hunter = F,
#                                 bump_hunter_data)
# 
# even_full_sam_table_enet <- extractResults(even_full_sam_enet, data_name = 'even sam_full enet', regularize = F)
# 
# #here
# # get table 
# even_full_sam <- rbind(even_full_sam_table_rf,
#                        even_full_sam_table_enet)
# 
# 
# rm(even_full_sam_table_rf,
#    even_full_sam_table_enet)
# 
# 
# 
# saveRDS(even_full_sam, paste0(model_data, '/even_full_2.rds'))

########################################################################3
##########
# run models
##########

#### 10
sam_fwer_10_rf <- runModels(quan_controls, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_fwer_10,
                         gender= F,
                         residual = F)

sam_fwer_10_table_rf <- extractResults(sam_fwer_10_rf, 
                                    data_name = 'sam_10_fwer_rf', 
                                    regularize = F,
                                    bh_data = quan_even_fwer_10)


#### 10
sam_fwer_10_enet <- runModels(quan_controls, 
                         model = 'enet',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_fwer_10,
                         gender= F,
                         residual = F)

sam_fwer_10_table_enet <- extractResults(sam_fwer_10_enet, 
                                    data_name = 'sam_10_fwer_enet', 
                                    regularize = T,
                                    bh_data = quan_even_fwer_10)

#############################################################################################################################3
###################################################################################################################################

