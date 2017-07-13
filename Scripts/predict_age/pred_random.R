#!/hpf/tools/centos6/R/3.2.3/bin/Rscript
argv <- as.numeric(commandArgs(T))

##########
# This script will predict with random features, increasing.

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
method = 'funnorm'
k = 5
seed_num <- 1
seed_num <- argv[1]


##########
# load data
##########
betaCases <- readRDS(paste0(model_data, '/betaCases', method,'.rda'))
betaControls <- readRDS(paste0(model_data, '/betaControls', method,'.rda'))


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

###########################################################################
# Next part of the pipline selects regions of the genome that are most differentially methylated 
# between 2 groups

# get a column for each dataset indicating the fold
betaCases <- getFolds(betaCases, 
                      seed_number = seed_num, 
                      k_num = k)


# betaCases <- betaCases[, c(1:5000, ncol(betaCases))]
# betaControls <- betaControls[, c(1:5000, ncol(betaControls))]
# betaControlsWT <- betaControlsWT[, c(1:5000, ncol(betaControlsWT))]
# cases <- betaCases
# num_feats <- 100
trainTestRand <- function(cases, 
                          num_feats,
                          k) 
{
  
  # remove samples that dont have an age of sample collection
  cases <- cases[complete.cases(cases),]
  
  # get random sample of column names
  features <- colnames(cases)[5:ncol(cases)]
  
  # sample num_feats
  set.seed(seed_num)
  sub_feats <- sample(features, num_feats)
  
  # list to store results
  model_results <- list()
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, cases$folds)
    test_index <- !train_index
    
    # run regression model
    mod_result <- runEnet(training_dat = cases[train_index,], 
                          test_dat = cases[test_index,], 
                          bh_features = sub_feats,
                          gender = T)
  
    model_results[[i]] <- getResultsCancer(mod_result)
    
    
  }
  
  return(model_results)
  
}

mod_results <- trainTestRand(cases = betaCases,
                             num_feats = 800,
                             k = k)

# change pred to nothing if doing surv
saveRDS(mod_results, 
        paste0('/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/LFS/Scripts/predict_age/Results/reg_results/train_test', '_', 'pred' , '_' ,seed_num, '.rda'))


