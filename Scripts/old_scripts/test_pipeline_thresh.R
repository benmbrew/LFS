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
# function for model
##########

run_model <- 
  function(method, type, thresh, set) {
    
    ##########
    # load data
    ##########
    
    if (method != 'raw') {
      
      # not raw, so load in raw columns to subset quan or funnorm
      raw_cols <- colnames(readRDS(paste0(model_data, paste0('/controls_no_transform','_', set,'_','raw','_' , thresh, '.rda'))))
      betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new.rda')))[, raw_cols]
      
      # if (type == 'transform') {
      #   
      #   betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new.rda')))[, raw_cols]
      #   betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new.rda')))[, raw_cols]
      
      
      betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new.rda')))[, raw_cols]
      betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new.rda')))[, raw_cols]
      
    } else {
      
      betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new.rda')))
      # 
      # if (type == 'transform') {
      #   
      #   betaControls <- readRDS(paste0(model_data, paste0('/controls_transform','_' , method,'_' , thresh, '.rda'))) # 34 435813
      #   betaValid <- readRDS(paste0(model_data, paste0('/valid_transform','_' , method,'_' , thresh, '.rda'))) # 45 400842
      #   
      
      betaControls <- readRDS(paste0(model_data, paste0('/controls_no_transform','_' ,set,'_', method,'_' , thresh, '.rda'))) # 34 435813
      betaValid <- readRDS(paste0(model_data, paste0('/valid_no_transform','_' ,set,'_' , method,'_' , thresh, '.rda')))# 45 400842
      
      
    }
    
    ###########
    # make id into ids
    ###########
    colnames(betaCases)[1] <- 'ids'
    colnames(betaControls)[1] <- 'ids'
    colnames(betaValid)[1] <- 'ids'
    
    ###########
    # get extra controls from 450k 
    ###########
    betaControlsOld <- getControls(betaCases, mut = T)
    
    ###########
    # get controls wt from betaCases
    ###########
    betaControlsWT <- subset(betaCases, p53_germline == 'WT' & cancer_diagnosis_diagnoses == 'Unaffected')
    
    
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
    
    # controlswt
    betaControlsWT <- betaControlsWT[, c('ids',
                                         'age_diagnosis', 
                                         'age_sample_collection', 
                                         'cancer_diagnosis_diagnoses', 
                                         'gender', 
                                         intersect_names)]
    
    # controls
    betaControlsOld <- betaControlsOld[, c('ids', 
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
    
    betaControls$ids <- as.character(betaControls$ids) 
    betaControls$cancer_diagnosis_diagnoses <- as.character(betaControls$cancer_diagnosis_diagnoses) 
    betaControls$gender <- as.character(betaControls$gender) 
    betaControls$age_diagnosis <- as.numeric(as.character(betaControls$age_diagnosis))
    betaControls$age_sample_collection <- as.numeric(as.character(betaControls$age_sample_collection))
    
    
    # get controls full
    betaControlsFull <- rbind(betaControls,
                              betaControlsOld)
    
    # remove duplicates from betaControlsFull
    length(which(duplicated(betaControlsFull$ids)))
    betaControlsFull <- betaControlsFull[!duplicated(betaControlsFull$ids),]
    
    # betaControlsWT
    length(which(duplicated(betaControlsWT$ids)))
    
    ##########
    # remove samples with no age data
    ##########
    betaControls <- betaControls[!is.na(betaControls$age_sample_collection),]
    betaControlsOld <- betaControlsOld[!is.na(betaControlsOld$age_sample_collection),]
    betaControlsFull <- betaControlsFull[!is.na(betaControlsFull$age_sample_collection),]
    
    # get gender dummy variable
    betaCases <- cbind(as.data.frame(class.ind(betaCases$gender)), betaCases)
    betaControls <- cbind(as.data.frame(class.ind(betaControls$gender)), betaControls)
    betaControlsOld <- cbind(as.data.frame(class.ind(betaControlsOld$gender)), betaControlsOld)
    betaControlsFull <- cbind(as.data.frame(class.ind(betaControlsFull$gender)), betaControlsFull)
    betaValid <- cbind(as.data.frame(class.ind(betaValid$gender)), betaValid)
    
    #subset valid
    betaValid <- betaValid[!betaValid$ids %in% betaCases$ids,]
    
    # betaCases <- betaCases[, c(1:500, ncol(betaCases))]
    # betaControls <- betaControls[, c(1:500, ncol(betaControls))]
    # betaValid<- betaValid[, c(1:500, ncol(betaValid))]
    
    alphas <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    
    result_list <- list()
    
    for(i in 1:length(alphas)) {
      
      alpha_num <- alphas[i]
      
      temp_results  <- testModelThresh(cases_dat = betaCases,
                                      controls_dat = betaControls,
                                      controls_dat_old = betaControlsOld,
                                      controls_dat_full  = betaControlsFull,
                                      valid_dat = betaValid,
                                      features = intersect_names,
                                      alpha = alpha_num)
      
      temp_results <- as.data.frame(t(unlist(temp_results)))
      temp_results$alpha <- alpha_num
      colnames(temp_results) <- c('controls', 'controls_old', 'controls_full', 'valid', 'alpha')
      
      result_list[[i]] <- temp_results
      
      print(i)
      
    }
    
    result_dat <- do.call(rbind, result_list)
    
    # assign identifiers
    result_dat$thresh <- thresh
    result_dat$set <- set
    result_dat$type <- 'no_transform'
    
    return(result_dat)
    
}



#########
# load results 
#########

# loop through method types and control sizes and aggregate results
method <- c('raw')
sets <- c('int', 'union')
thresh_holds <- c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09)


temp_i <- 
  temp_j <-
  temp_l <-list()

for (i in 1:length(method)) {
  
  for (j in 1:length(sets)) {
    
    for(l in 1:length(thresh_holds)) {
      
      temp_l[[l]] <- run_model(method = method[i], 
                               type = 'no_transform', 
                               set = sets[j],
                               thresh = thresh_holds[l])
      
      print(paste('finished', method[i], sets[j], thresh_holds[l] ,collapse = ' '))
    }
    
    temp_j[[j]] <- do.call(rbind, temp_l)
    
  }
  
  temp_i[[i]] <- do.call(rbind, temp_j)
  
}

# get full results
full_results <- do.call(rbind, temp_i)

# group by alpha

saveRDS(full_results, paste0(project_folder, paste0('/Scripts/predict_age/Results/reg_results_thresh/full_results_pipeline.rda')))
