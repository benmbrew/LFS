####### This script will predict age of onset and save results
# this is 7th step in pipeline

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
raw_folder <- paste0(results_folder, '/raw')
swan_folder <- paste0(results_folder, '/swan')
quan_folder <- paste0(results_folder, '/quan')
funnorm_folder <- paste0(results_folder, '/funnorm')
rand_folder <- paste0(results_folder, '/rand')
scripts_folder <- paste0(project_folder, '/Scripts')


##########
# source model_functions to get functions to run models 
##########
source(paste0(scripts_folder, '/predict_age/model_functions_short.R'))

##########
# set parameters 
##########
data_thresholds <- 48

##########
# Read in cases and bumphunter features
##########
# cases - raw, swan, quan, funnorm, and clin
# load(paste0(model_data, '/model_data_cases.RData'))

# load cases
load(paste0(model_data, '/model_data.RData'))

# bh_features
load(paste0(model_data, '/bh_feat.RData'))

##########
# Read in controls (if needed)
##########
# load(paste0(model_data, '/model_data_controls.RData'))
##########
# function to run models - subset, get residuals, get categorical, predict with regression and fac. 
##########
runModels <- function(data,
                      model,
                      random = F,
                      bump_hunter = F,
                      bump_hunter_data,
                      num_feat = NULL,
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
                      num_feat = NULL)
  }
  
  #########
  # get regression data and run regression
  #########
  # get data
  if (bump_hunter) {
    
    data <- bhSubset(data, bh_data = bump_hunter_data)
    
  }
  
  # data_resid <- getResidual(data)
  
  if (model == 'rf') {
    
    data_result <- rfPredictReg(data, cutoff = .7, iterations = 5)
    # data_resid_result <- rfPredictReg(data_resid, cutoff = .7, iterations = 5)
  }
  
  if (model == 'enet') {
    
    data_result <- enetPredReg(data, N_CV_REPEATS = 2, nfolds = 5,cutoff = .7, iterations = 5)
    # data_resid_result <- enetPredReg(data_resid, N_CV_REPEATS = 2, nfolds = 5,cutoff = .7, iterations = 5)
    
  }
  
  if (model == 'lasso') {
    
    data_result <- regPred(data, alpha = 1, nfolds = 5, cutoff = .7, iterations = 5)
    # data_resid_result <- regPred(data_resid, alpha = 1, nfolds = 5, cutoff = .7, iterations = 5)
    
  }
  
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

##########
# get stats of data used in model
##########
# get data info from 
data_stats <- dataStats(raw_cases)

saveRDS(data_stats, 
        file = paste0(results_folder, 
                      '/model_data_stats.rda'))

# remove object
rm(data_stats)
