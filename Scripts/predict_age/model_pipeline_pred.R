#!/hpf/tools/centos6/R/3.2.3/bin/Rscript

# Script to evaluate the how each imputation method affects the
# performance of the clustering methods

argv <- as.numeric(commandArgs(T))

##########
# This script will get cleaned data from data saved in clean_data.R
# this will be the modeling pipeline script where we select features on 
# our training set and fit model. On the test set we test our model 
# and compare our predictions to values. 

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
home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
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
method = 'raw'
k = 4
seed_num <- 1
seed_num <- argv[1]


##########
# load data
##########
betaCases <- readRDS(paste0(model_data, '/betaCases', method,'.rda'))
betaControls <- readRDS(paste0(model_data, '/betaControls', method,'.rda'))
betaControlsWT <- readRDS(paste0(model_data, '/betaControlsWT', method,'.rda'))
betaControlsOld <- readRDS(paste0(model_data, '/betaControlsOld', method,'.rda'))


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
                                          colnames(betaControls)[8:ncol(betaControls)]))
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

# controls wild type
betaControlsWT <- betaControlsWT[, c('age_diagnosis', 
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

# #TEMP
# betaControlsFull <- rbind(betaControls, betaControlsOld)
# 
# # remove na in sample collection and duplicate ids
# betaControlsFull <- betaControlsFull[!is.na(betaControlsFull$age_sample_collection),]


###########################################################################
# Next part of the pipline selects regions of the genome that are most differentially methylated 
# between 2 groups

# get a column for each dataset indicating the fold
betaCases <- getFolds(betaCases, seed_number = seed_num, k_num = k)
betaControls <- getFolds(betaControls, seed_number = seed_num, k_num = k)
betaControlsWT <- getFolds(betaControlsWT, seed_number = seed_num, k_num = k)

# betaCases <- betaCases[, c(1:5000, ncol(betaCases))]
# betaControls <- betaControls[, c(1:5000, ncol(betaControls))]
# betaControlsWT <- betaControlsWT[, c(1:5000, ncol(betaControlsWT))]


# high pvalue, no evidence they are different
# print(testKS(cases$age_sample_collection[train_index], controls$age_sample_collection)$p.value)
print(testKS(betaControlsWT$age_sample_collection, betaControls$age_sample_collection)$p.value)

# subset betaControlsWT by removing a sample with age of sample collection between 300 and 400
# iteratively, each time running until p.value of testKS is above 0.05
p_value =  0
loop_counts=  0

while (p_value < 0.05) {
  
  if (loop_counts == 0) {
    
    remove_index <- which(betaControlsWT$age_sample_collection >= 300 & 
                            betaControlsWT$age_sample_collection <= 400)
    remove_index <- remove_index[1]
    
    
    betaControlsWTSub <- betaControlsWT[-remove_index,]
  } else {
    
    remove_index <- which(betaControlsWTSub$age_sample_collection >= 300 & 
                            betaControlsWTSub$age_sample_collection <= 400)
    
    remove_index <- remove_index[1]
    
    
    betaControlsWTSub <- betaControlsWTSub[-remove_index,]
  }
 
  p_value <- testKS(betaControlsWTSub$age_sample_collection, betaControls$age_sample_collection)$p.value
  
  loop_counts <- loop_counts + 1
}
# # adjust for age
# hist(betaControlsWTSub$age_sample_collection)
# hist(betaControls$age_sample_collection)
bh_feat <- bumpHunterPred(dat_controls_wt = betaControlsWTSub, dat_controls_mut = betaControls)

# get probes with regions
bh_feat_3 <- getProbe(bh_feat)

# saveRDS(bh_feat_3, paste0(model_data, '/bh_feat_pred.rda'))

# get all data sets from bh_feat_3
# bh_feat_all <- getRun(bh_feat_3[[1]], run_num = .20)
bh_feat_sig <- getRun(bh_feat_3[[2]], run_num = .20)
# bh_feat_fwer <- getRun(bh_feat_3[[3]], run_num = seed_num)

trainTest <- function(cases, 
                      bh_feat_sig,
                      k) 
{
  
  # remove samples that dont have an age of sample collection
  cases <- cases[complete.cases(cases),]
  
  # get residuals
  cases_resid <- getResidual(data = cases, 
                             bh_features = bh_feat_sig)
  
  # list to store results
  bh_feat <- list()
  model_results <- list()
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, cases$folds)
    test_index <- !train_index
    
    # run regression model
    mod_result <- runEnet(training_dat = cases[train_index,], 
                          test_dat = cases[test_index,], 
                          bh_features = bh_feat_sig,
                          gender = T)
    
    
    # run classification model
    mod_result_resid <- runEnet(training_dat = cases_resid[train_index,], 
                                test_dat = cases_resid[test_index,], 
                                bh_features = bh_feat_sig,
                                gender = T)
    
    
    
    model_results[[i]] <- getResults(mod_result, mod_result_resid)
    
    
  }
  
  return(model_results)
  
}

mod_results <- trainTest(cases = betaCases,
                         bh_feat_sig = bh_feat_sig,
                         k = 4)

seed_num
# change pred to nothing if doing surv
saveRDS(mod_results, 
        paste0('/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/LFS/Scripts/predict_age/Results/reg_results/train_test', '_', 'pred' , '_' ,seed_num, '.rda'))

