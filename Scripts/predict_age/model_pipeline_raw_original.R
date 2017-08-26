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
home_folder <- '/hpf/largeprojects/agoldenb/ben/Projects'
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
k = 5
seed_num <- argv[1]
type = 'original'

##########
# load data
##########
if (type == 'original') {
  
  betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new.rda')))
  betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new.rda'))) #34 449936
  betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new.rda'))) #35 449783
  
} else if (type == 'transform') {
  
  betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new.rda')))
  betaControls <- readRDS(paste0(model_data, paste0('/controls_transform','_' , method, '.rda'))) # 34 435813
  betaValid <- readRDS(paste0(model_data, paste0('/valid_transform','_' , method, '.rda'))) # 45 400842
  
} else {
  
  betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new.rda')))
  betaControls <- readRDS(paste0(model_data, paste0('/controls_no_transform','_' , method, '.rda'))) # 34 435813
  betaValid <- readRDS(paste0(model_data, paste0('/valid_no_transform','_' , method, '.rda')))# 45 400842
  
}

colnames(betaCases)[1] <- 'ids'
colnames(betaControls)[1] <- 'ids'
colnames(betaValid)[1] <- 'ids'

###########
# get model data
###########
betaCases <- getModData(betaCases)

# get rid of cancer samples in controls 
betaControls <- betaControls[grepl('Unaffected', betaControls$cancer_diagnosis_diagnoses),]


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
betaCases <- betaCases[, c('ids',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           intersect_names)]
# controls
betaControls <- betaControls[, c('ids',
                                 'age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender', 
                                 intersect_names)]

#validation
betaValid <- betaValid[, c('ids',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           intersect_names)]

betaValid <- betaValid[!betaValid$ids %in% betaCases$ids,]


###########################################################################
# Next part of the pipline selects regions of the genome that are most differentially methylated 
# between 2 groups

# get a column for each dataset indicating the fold
betaCases <- getFolds(betaCases, seed_number = seed_num, k = k)
betaControls <- getFolds(betaControls, seed_number = seed_num, k = k)

# get gender 
# get gender dummy variable
betaCases <- cbind(as.data.frame(class.ind(betaCases$gender)), betaCases)
betaControls <- cbind(as.data.frame(class.ind(betaControls$gender)), betaControls)
betaValid <- cbind(as.data.frame(class.ind(betaValid$gender)), betaValid)


# betaCases <- betaCases[, c(1:3000, ncol(betaCases))]
# betaControls <- betaControls[, c(1:3000, ncol(betaControls))]
# betaValid <- betaValid[, c(1:3000, ncol(betaValid))]

trainTest <- function(cases, 
                      controls,
                      valid,
                      k) 
{
  
  # remove samples that dont have an age of sample collection
  cases <- cases[complete.cases(cases),]
  
  # list to store results
  bh_feat <- list()
  model_results <- list()
  bh_dim <- list()
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, cases$folds)
    test_index <- !train_index
    
    # high pvalue, no evidence they are different
    # print(testKS(cases$age_sample_collection[train_index], controls$age_sample_collection)$p.value)
    # print(testKS(controls_wt$age_sample_collection[train_index], controls$age_sample_collection)$p.value)
    
    # use bumphunter surveillance function to get first set of regions
    bh_feat[[i]] <- bumpHunterSurv(dat_cases = cases[train_index,], dat_controls = controls)
    
    # get probes with regions
    bh_feat_3 <- getProbe(bh_feat[[i]])
    
    # get all data sets from bh_feat_3
    bh_feat_all <- getRun(bh_feat_3[[1]], run_num = .05)
    # bh_feat_sig <- getRun(bh_feat_3[[2]], run_num = .10)
    bh_dim[[i]] <- length(bh_feat_all)
    # bh_feat_fwer <- getRun(bh_feat_3[[3]], run_num = seed_num)
    
    # get residuals
    # cases_resid <- getResidual(data = cases,
    #                            bh_features = bh_feat_all)
    
    mod_result <- runEnet(training_dat = cases[train_index,], 
                          test_dat = cases[test_index,], 
                          controls_dat = controls,
                          valid_dat = valid,
                          bh_features = bh_feat_all,
                          gender = T)
    
    
    # mod_result_resid <- runEnet(training_dat = cases_resid[train_index,], 
    #                             test_dat = cases_resid[test_index,], 
    #                             controls_dat = controls,
    #                             valid_dat = valid,
    #                             bh_features = bh_feat_all,
    #                             gender = T)
    # 
    
    
    model_results[[i]] <- getResults(mod_result)
    
    
  }
  
  return(list(model_results, bh_dim))
  
}

mod_results <- trainTest(cases = betaCases,
                         controls = betaControls,
                         valid = betaValid,
                         k = k)

# change pred to nothing if doing surv
saveRDS(mod_results, paste0('/hpf/largeprojects/agoldenb/ben/Projects/LFS/Scripts/predict_age/Results/reg_results/train_test', '_' , seed_num,'_', type, '_', method, '.rda'))

