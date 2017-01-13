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
source(paste0(scripts_folder, '/predict_age/model_functions.R'))

##########
# set parameters 
##########
data_thresholds <- c(48, 60, 72, 84)
p53 <- c('Mut', 'WT')

##########
# Read in cases and bumphunter features
##########
# cases - raw, swan, quan, funnorm, and clin
load(paste0(model_data, '/model_data_cases.RData'))
# bh_features
load(paste0(model_data, '/bh_features.RData'))

##########
# Read in controls (if needed)
##########
# load(paste0(model_data, '/model_data_controls.RData'))

##########
# function to run models - subset, get residuals, get categorical, predict with regression and fac. 
##########
runModels <- function(data,
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
  
  if (bump_hunter) {
    
    data <- bhSubset(data, bh_data = bump_hunter_data)
    
  }
  
  data_resid <- getResidual(data)
  
  data_fac <- list()
  data_resid_fac <- list()
  
  for (thresh in 1:length(data_thresholds)) {
    
    data_fac[[thresh]] <- makeFac(data, threshold = data_thresholds[thresh])
    data_resid_fac[[thresh]] <- makeFac(data_resid, threshold = data_thresholds[thresh])
    
  }
  
  # run regressions
  data_result <- list()
  data_resid_result <- list()
  
  for (dat in 1:length(data)) {
    
    sub_dat <- data[[dat]]
    sub_data_resid <- data_resid[[dat]]
    data_result[[dat]] <- rfPredictReg(sub_dat, cutoff = .7, iterations = 10)
    data_resid_result[[dat]] <- rfPredictReg(sub_data_resid, cutoff = .7, iterations = 10)
  }
  
  # run classification
  data_fac_result <- list()
  data_resid_fac_result <- list()
  
  for (dat in 1:length(data_fac)) {
    
    sub_dat_fac <- data_fac[[dat]]
    sub_dat_resid_fac <- data_resid_fac[[dat]]
    
    temp.data_fac_result <- list()
    temp.data_resid_fac_result <- list()
    
    for (sub_dat in 1:length(sub_dat_fac)) {
      
      temp.sub_dat_fac <- sub_dat_fac[[sub_dat]]
      temp.sub_dat_resid_fac <- sub_dat_resid_fac[[sub_dat]]
      
      temp.data_fac_result[[sub_dat]] <- rfPredictFac(temp.sub_dat_fac, cutoff = .7, iterations = 10)
      temp.data_resid_fac_result[[sub_dat]] <- rfPredictFac(temp.sub_dat_resid_fac, cutoff = .7, iterations = 10)
      
    }
    
    data_fac_result[[dat]] <- temp.data_fac_result
    data_resid_fac_result[[dat]] <- temp.data_resid_fac_result
    
  }
  
  return (list(data_result, 
               data_resid_result, 
               data_fac_result, 
               data_resid_fac_result))
  
}

##########
# get stats of data used in model
##########
# get data info from 
data_stats <- dataStats(beta_raw)

saveRDS(data_stats, 
        file = paste0(results_folder, 
                      '/model_data_stats.rda'))

# remove object
rm(data_stats)
