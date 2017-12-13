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
library(ROCR)
library(ModelMetrics)
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
age = 48
seed_num <- 1
bh_method <- 'pred'
seed_num <- argv[1]



##########
# load data
##########
betaCases <- readRDS(paste0(model_data, '/betaCases', method,'.rda'))
betaControls <- readRDS(paste0(model_data, '/betaControls', method,'.rda'))
betaControlsWT <- readRDS(paste0(model_data, '/betaControlsWT', method,'.rda'))
betaControlsOld <- readRDS(paste0(model_data, '/betaControlsOld', method,'.rda'))
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

betaCases <- betaCases[, c(1:3000, ncol(betaCases))]
betaControls <- betaControls[, c(1:3000, ncol(betaControls))]
betaControlsWT <- betaControlsWT[, c(1:3000, ncol(betaControlsWT))]

# read in bh data
bh_dat <- readRDS(paste0(model_data, '/bh_feat_pred.rda'))

# bh_sig
bh_feat_sig <- getRun(bh_dat[[2]], run_num = .20)


trainTestFac <- function(cases, 
                         age_cutoff,
                         k) 
{
  
  # list to store results
  bh_feat <- list()
  test_stats <- list()
  alpha_score <- list()
  lambda_num <- list()
  import <- list()
  models <- list()
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, cases$folds)
    test_index <- !train_index
    
    mod_result <- runEnetFac(training_dat = cases[train_index,],
                             test_dat = cases[test_index,], 
                             bh_features = bh_feat_sig,
                             gender = T,
                             cutoff = age_cutoff)
    
    
    alpha_score[[i]] <- mod_result[[1]]
    lambda_num[[i]] <- mod_result[[2]]
    import[[i]] <- mod_result[[3]]
    test_stats[[i]] <- mod_result[[4]]
    models[[i]] <- mod_result[[5]]
    
    
  }
  
  return(list(alpha_score, lambda_num, import, test_stats, models))
  
}

mod_results <- trainTestFac(cases = betaCases,
                            age_cutoff = age,
                            k = 4)

saveRDS(mod_results, paste0('/hpf/largeprojects/agoldenb/ben/Projects/LFS/Scripts/
                            predict_age/Results/class_results/train_test', '_' , bh_method, '_',
                            seed_num, '_' ,age,'.rds'))

