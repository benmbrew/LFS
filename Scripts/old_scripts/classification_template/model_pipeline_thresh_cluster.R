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
# function for model
##########

run_model <- 
  function(method, k, seed_num, type, thresh, set) {
    
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
    
    
    colnames(betaCases)[1] <- 'ids'
    colnames(betaControls)[1] <- 'ids'
    colnames(betaValid)[1] <- 'ids'
    
    ###########
    # get model data
    ###########
    betaCases <- getModData(betaCases)
    
    # get rid of cancer samples in controls 
    betaControls <- betaControls[grepl('Unaffected', betaControls$cancer_diagnosis_diagnoses),]
    
    # remove valid samples that are already present in betaCases
    betaValid <- betaValid[!betaValid$ids %in% betaCases$ids,]
    
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
    
    # get gender dummy variable
    betaCases <- cbind(as.data.frame(class.ind(betaCases$gender)), betaCases)
    betaControls <- cbind(as.data.frame(class.ind(betaControls$gender)), betaControls)
    betaValid <- cbind(as.data.frame(class.ind(betaValid$gender)), betaValid)
    # # 
    # betaCases <- betaCases[, c(1:500, ncol(betaCases))]
    # betaControls <- betaControls[, c(1:500, ncol(betaControls))]
    # betaValid<- betaValid[, c(1:500, ncol(betaValid))]
    
    ##########
    # function for training and testing
    ##########
    
    trainTest <- function(cases, 
                          controls, 
                          valid,
                          features,
                          k) {
      
      # remove samples that dont have an age of sample collection
      cases <- cases[complete.cases(cases),]
      
      dimensions <- ncol(cases)
      
      # list to store results
      model_results <- list()
      
      # now write forloop to 
      for (i in 1:k) {
        
        # get x 
        train_index <- !grepl(i, cases$folds)
        test_index <- !train_index
        
        
        mod_result <- runEnet(training_dat = cases[train_index,], 
                              test_dat = cases[test_index,], 
                              controls_dat = controls,
                              valid_dat = valid,
                              bh_features = features,
                              gender = T)
        
        
        model_results[[i]] <- getResults(mod_result)
        
        
      }
      
      return(list(model_results, dimensions))
      
    }
    
    mod_results <- trainTest(cases = betaCases,
                             controls = betaControls,
                             valid = betaValid,
                             features = intersect_names,
                             k = k)
    
    return(mod_results)
    
  }

# method = 'raw'
# type = 'transform'
# iter_length = 2
# k = 5
# thresh = 0.05

get_results <- 
  function(method, 
           type, 
           set,
           seed_num,
           k,
           thresh) {
    
    fold_results <- list()
    reg_results <- list()

    # read in raw results
    temp_results <- run_model(method = method, 
                              k = k, 
                              seed_num = seed_num, 
                              type = type, 
                              thresh = thresh,
                              set = set)
    
    temp.result_folds <- temp_results[[1]]
    
    temp.dims <- temp_results[[2]]
      
      
      for (j in 1:k){
        # get fold - returns list of 4
        temp.fold <- temp.result_folds[[j]]
        
        
        # get results - list of 2 - normal and resid
        temp.fold_norm <- temp.fold
        # temp.fold_resid <- temp.fold[[2]]
        
        #alpha, lambda_value, importance, cases_cor, age_cor, controls_cor, valid_cor, temp.non_zero_coeff
        
        # remove 3rd element 
        temp.fold_norm[[3]] <- NULL
        # temp.fold_resid[[3]] <- NULL
        
        
        # list of 5 - alpha, lambda_value, cases_cor, age_cor, temp.non_zero_coeff
        temp.fold_norm_dat <- as.data.frame(t(do.call(rbind, temp.fold_norm)))
        # temp.fold_resid_dat <- as.data.frame(t(do.call(rbind, temp.fold_resid)))
        
        # add in resid and norm
        temp.fold_norm_dat$V8 <- 'norm'
        # temp.fold_resid_dat$V8 <- 'resid'
        
        # get column names 
        colnames(temp.fold_norm_dat) <- c('alpha', 
                                          'lambda', 
                                          'onset_correlation', 
                                          'age_correlation', 
                                          'controls_cor',
                                          'valid_cor',
                                          'vars', 
                                          'type')
        
        # colnames(temp.fold_resid_dat) <- c('alpha', 
        #                                    'lambda', 
        #                                    'onset_correlation', 
        #                                    'age_correlation', 
        #                                    'controls_cor',
        #                                    'valid_cor',
        #                                    'vars', 
        #                                    'type')
        
        # add in indicator for norm and resid
        temp.fold_norm_dat$seed_num <- i
        temp.fold_norm_dat$method <- paste0(method, '_', 'no_transform','_' ,thresh, '_', set)
        # temp.fold_resid_dat$seed_num <- i
        
        
        # combine 
        temp.result_folds_dat <- temp.fold_norm_dat
        
        # store in list
        # fold_results[[j]] <- temp.result_folds
        fold_results[[j]] <- temp.result_folds_dat
        
        
      }
      
    # collapse fold_results
    temp.collpased <- as.data.frame(do.call(rbind, fold_results))

    
    # get temp.dims
    dim_dat <- unlist(temp.dims)[1]
    
    # result_table make dim  column
    temp.collpased$mean_dim <- dim_dat
    
    # unlist alpha 
    temp.collpased$alpha <- unlist(temp.collpased$alpha)
    ##########
    # analyze results 
    ##########
    
    # # first group alpha and get mean correlation 
    # alpha_result <- result_table %>%
    #   group_by(method, alpha) %>%
    #   summarise(mean_cor = mean(onset_correlation),
    #             mean_cor_age = mean(age_correlation),
    #             mean_cor_controls = mean(controls_cor),
    #             mean_cor_valid = mean(valid_cor),
    #             mean_vars_import = mean(vars),
    #             mean_dim = mean(mean_dim))
    # 
    return(temp.collpased)
    
  }

#########
# load results 
#########


# seed number
seed_num <- argv[1]

# loop through method types and control sizes and aggregate results
method <- c('raw', 'quan', 'funnorm')
sets <- c('int', 'union')
thresh_holds <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.5)

temp_i <- 
  temp_j <-
  temp_l <-list()

for (i in 1:length(method)) {
  
  for (j in 1:length(sets)) {
    
    for(l in 1:length(thresh_holds)) {
      
      temp_l[[l]] <- get_results(method = method[i], 
                                 type = 'no_transform', 
                                 set = sets[j],
                                 seed_num = seed_num, 
                                 k = 5,
                                 thresh = thresh_holds[l])
      
      print(paste('finished', method[i], sets[j], thresh_holds[l] ,collapse = ' '))
    }
    
    temp_j[[j]] <- do.call(rbind, temp_l)
    
  }
  
  temp_i[[i]] <- do.call(rbind, temp_j)
  
}

# get full results
full_results <- do.call(rbind, temp_i)

# # group by alpha 
# alpha_results <- full_results %>%
#   group_by(alpha) %>%
#   summarise(cor = mean(mean_cor),
#             cor_age = mean(mean_cor_age),
#             cor_con = mean(mean_cor_controls),
#             cor_valid = mean(mean_cor_valid))


saveRDS(full_results, paste0(project_folder, paste0('/Scripts/predict_age/Results/reg_results_thresh/full_results_thresh', '_',seed_num,'.rda')))
