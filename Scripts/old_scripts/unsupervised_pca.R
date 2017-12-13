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
k = 5
seed_num <- 1
set = 'int'
thresh = 0.08
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


# homogenize ids column name
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


# get folds
betaCases <- getFolds(betaCases, seed_number = seed_num, k = k)

###########################################################################
# Next part of the pipline selects regions of the genome that are most differentially methylated 
# between 2 groups


# get gender 
# get gender dummy variable
# betaCases <- cbind(as.data.frame(class.ind(betaCases$gender)), betaCases)
# betaControls <- cbind(as.data.frame(class.ind(betaControls$gender)), betaControls)
# betaValid <- cbind(as.data.frame(class.ind(betaValid$gender)), betaValid)


# betaCases <- betaCases[, c(1:300, ncol(betaCases))]
# betaControls <- betaControls[, c(1:300, ncol(betaControls))]
# betaValid <- betaValid[, c(1:300, ncol(betaValid))]

# 
# cases <- betaCases
# controls <- betaControls
# valid <- betaValid
# probe_start = 6

trainTest <- function(cases, 
                      controls,
                      valid,
                      probe_start,
                      k) {
  
  # remove samples that dont have an age of sample collection
  cases <- cases[complete.cases(cases),]
  
  # get clin data
  clin_dat <- cases[, 1:probe_start -1]
  # creat features
  feature_names <- colnames(cases)[probe_start:(ncol(cases) - 1)]
  stopifnot(all(grepl('cg', feature_names)))
  
  # list to store results
  bh_feat <- list()
  model_results <- list()
  bh_dim <- list()
  cor_pred <- list()
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, cases$folds)
    test_index <- !train_index
    
    # get train x and train_y
    train_x <- data.matrix(cases[train_index, feature_names])
    train_y <- as.vector(cases$age_diagnosis[train_index])
    
    # get test x and test y
    test_x <- data.matrix(cases[!train_index, feature_names])
    test_y <- as.vector(cases$age_diagnosis[!train_index])
    
    # gender train and test 
    
    gen_train <- ifelse(grepl("M", cases$gender[train_index]), 1, 0)
    gen_test <- ifelse(grepl("M", cases$gender[!train_index]), 1, 0)
    
    # check sizes
    stopifnot(dim(train_x)[1] == length(train_y))
    stopifnot(dim(test_x)[1] == length(test_y))
    
    # get pcas of train_x
    pca_train_x <- prcomp(train_x)
    pca_test_x <- prcomp(test_x)
    
    # get first PCs
    get_pc <- function(pc_obj) {
      temp_col <- list()
      for(j in 1:5) {
        temp_col[[j]] <- pc_obj$x[, j]
      }
      
      pc_data <- do.call(cbind, temp_col)
      
      return(pc_data)
    }
    
    pc_train_x <- get_pc(pca_train_x)
    pc_test_x <- get_pc(pca_test_x)
    
    # fit model on train y with pcax
    mod_data <- as.data.frame(cbind(y = train_y ,pc_train_x,gen = gen_train))
    mod_obj <- lm(y~., mod_data)
    
    # predict
    test_data <- as.data.frame(cbind(y= test_y, pc_test_x, gen = gen_test))
    pred_y <- predict(mod_obj, test_data[-1])
    
    # correlate
    cor_pred[[i]] <- cor(pred_y, test_data$y)
    
    print(i)
  }
  
  return(cor_pred)
  
}

mod_results <- trainTest(cases = betaCases,
                         controls = betaControls,
                         valid = betaValid,
                         probe_start = 6,
                         k = k)

mean(unlist(mod_results))
# change pred to nothing if doing surv
saveRDS(mod_results, paste0('/hpf/largeprojects/agoldenb/ben/Projects/LFS/Scripts/predict_age/Results/reg_results/train_test', '_' , seed_num,'_', type, '_', method, '.rda'))

