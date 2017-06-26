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
method = 'raw'
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

#TEMP
betaControlsFull <- rbind(betaControls, betaControlsOld)

# remove na in sample collection and duplicate ids
betaControlsFull <- betaControlsFull[!is.na(betaControlsFull$age_sample_collection),]


###########################################################################
# Next part of the pipline selects regions of the genome that are most differentially methylated 
# between 2 groups

# get a column for each dataset indicating the fold
betaCases <- getFolds(betaCases)
betaControls <- getFolds(betaControls)
betaValid <- getFolds(betaValid)


trainTest <- function(betaCases, 
                      betaControls, 
                      betaValid, 
                      k) 
{
  
  # list to store results
  bh_feat <- list()
  
  # now write forloop to 
  for (i in 1:k) {

    # get x 
    train_index <- !grepl(i, betaCases$folds)
    test_index <- !train_index
  
    # high pvalue, no evidence they are different
    print(testKS(betaCases$age_sample_collection[train_index], betaControls$age_sample_collection)$p.value)

    # use bumphunter surveillance function to get first set of regions
    bh_feat[[i]] <- bumpHunterSurv(dat_cases = betaCases[train_index,], dat_controls = betaControls)
   
    # get probes with regions
    bh_feat_3 <- getProbe(bh_feat[[i]])
    
    # get all data sets from bh_feat_3
    bh_feat_all <- getRun(bh_feat_3[[2]], run_num = .10)
    # bh_feat_sig <- getRun(bh_feat_3[[2]], run_num = seed_num)
    # bh_feat_fwer <- getRun(bh_feat_3[[3]], run_num = seed_num)
    
    # save.image('/home/benbrew/Desktop/temp_pipeline_model.RData')
    # load('/home/benbrew/Desktop/temp_pipeline_model.RData')
    
   
    
    
  }
  
}
  
  


