## this script will read in batch data from get_cases, get_controls, or get_valid and explore potential batches and outliers

##########
# initialize libraries
##########
library(tidyverse)
library(sva)
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(ModelMetrics)
library(doParallel)
library(nnet)
library(dplyr)
library(bumphunter)
library(sqldf)
library(e1071)
library(reshape2)

registerDoParallel(1)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
clin_data <- paste0(data_folder, '/clin_data')

# get method 
method = 'noob'

##########
# load data
##########
  
# read in data
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'beta_cases.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'beta_controls.rda')))
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'beta_valid.rda')))

  
##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))



##########
# if you use transform data, homogenize here
##########
intersect_names <- Reduce(intersect, list(colnames(betaCases)[9:ncol(betaCases)], 
                                          colnames(betaControls)[9:ncol(betaControls)], 
                                          colnames(betaValid)[9:ncol(betaValid)]))



# add column in
betaCases$batch <- 'cases'
betaControls$batch <- 'controls'
betaValid$batch <- 'valid'


# cases
betaCases <- betaCases[, c('ids',
                           'p53_germline',
                     'age_diagnosis', 
                     'age_sample_collection', 
                     'cancer_diagnosis_diagnoses', 
                     'gender',
                     'batch',
                     'family_name',
                     intersect_names)]


# controls
betaControls <- betaControls[, c('ids',
                                 'p53_germline',
                     'age_diagnosis', 
                     'age_sample_collection', 
                     'cancer_diagnosis_diagnoses', 
                     'gender', 
                     'batch',
                     'family_name',
                     intersect_names)]

#validation
betaValid <- betaValid[, c('ids',
                           'p53_germline',
                     'age_diagnosis', 
                     'age_sample_collection', 
                     'cancer_diagnosis_diagnoses', 
                     'gender',
                     'batch',
                     'family_name',
                     intersect_names)]


# combine for model data
full_data <- rbind(betaCases,
                   betaControls,
                   betaValid)

full_data_combat <- run_combat(full_data)


##########
# get mod data
##########

sep_combat_dat <- function(combat_data){
  

  beta_cases <- combat_data[combat_data$batch == 'cases',]
  beta_controls <- combat_data[combat_data$batch == 'controls',]
  beta_valid <- combat_data[combat_data$batch == 'valid',]
  
  
  beta_cases <- getModData(beta_cases)
  
  beta_controls <- subset(beta_controls, p53_germline == 'Mut' & 
                            cancer_diagnosis_diagnoses == 'Unaffected')
  
  beta_valid <- beta_valid[!beta_valid$ids %in% beta_cases$ids,]
  
  combined_data <- rbind(beta_cases,
                         beta_controls,
                         beta_valid)
  
  
  return(combined_data)
  
}

mod_data_combat <- sep_combat_dat(full_data_combat, analysis = 'normal' )
mod_data <- sep_combat_dat(full_data, analysis = 'normal')

saveRDS(mod_data_combat, paste0(model_data, paste0('/', method, '_', 'mod_data_combat.rda')))
saveRDS(mod_data, paste0(model_data, paste0('/', method, '_', 'mod_data.rda')))






