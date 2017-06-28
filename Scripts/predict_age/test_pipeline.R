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
bh_feat_tot <- getRun(bh_feat_all[[1]], run_num = .20)
bh_feat_sig <- getRun(bh_feat_all[[2]], run_num = .20)

# save.image('/home/benbrew/Desktop/temp_full_test.RData')
# load('/home/benbrew/Desktop/temp_full_test.RData')

##########
# train on full data and test on controls and valid 
##########








