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
library(DAAG)
library(superpc)
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
probe_start = 6
k = 3
seed_num <- 5

##########
# load data
##########
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new.rda'))) #34 449936
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new.rda'))) #35 449783

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
# # controls
# betaControls <- betaControls[, c('ids',
#                                  'age_diagnosis', 
#                                  'age_sample_collection', 
#                                  'cancer_diagnosis_diagnoses', 
#                                  'gender', 
#                                  intersect_names)]
# 
# #validation
# betaValid <- betaValid[, c('ids',
#                            'age_diagnosis', 
#                            'age_sample_collection', 
#                            'cancer_diagnosis_diagnoses', 
#                            'gender', 
#                            intersect_names)]
# 
# betaValid <- betaValid[!betaValid$ids %in% betaCases$ids,]
# # 
# set.seed(464)
#
# #
# x<-matrix(rnorm(1000*100),ncol=100)
# v1<- svd(x[1:80,])$v[,1]
# 
# y<-2+5*v1+ .05*rnorm(100)
# #
# xtest<-x
# ytest<-2+5*v1+ .05*rnorm(100)
# 
# featurenames <- paste("feature",as.character(1:1000),sep="")
# 
# data<-list(x=x,y=y, featurenames=featurenames)
# data.test<-list(x=xtest,y=ytest,  featurenames= featurenames)#
# 

#  get complete cases 
betaCases <- betaCases[complete.cases(betaCases),]

# order data by id
betaCases <- betaCases[order(betaCases$ids), ]


##########
# get methylation matrix and clinical data
##########

# get model data list
mod_data<- get_model_dat(betaCases, probe_start = 6,
                         seed_num = seed_num, k = 5)

# get beta matrix, clinical variables, fold vector, and feature names for the supervised PCA
beta_methyl <- mod_data[[1]] # a transposed (p x n) methylation matrix
clin_data <- mod_data[[2]] # c
fold_vec <- mod_data[[3]]
feats <- mod_data[[4]]



##########
# function 
##########
train_test_pca <- 
  function(cases,
           clin, 
           fold_numbers,
           feature_names,
           cv_type) {
  
  
  # # creat features
  # feature_names <- colnames(cases)[probe_start:(ncol(cases) - 1)]
  # stopifnot(all(grepl('cg', feature_names)))
  
  # list to store results
  bh_feat <- list()
  model_results <- list()
  bh_dim <- list()
  test_results <- list()
  cv_object <- list()
  opt_feats <- list()
  opt_thresh <- list()
  
  
  # now write forloop to 
  for (i in 1:max(fold_numbers)) {
    
    # get x 
    train_index <- !grepl(i, fold_numbers)
    test_index <- !train_index
    
    # get train x and train_y
    train_x <- beta_methyl[,train_index]
    train_y <- as.numeric(clin$age_diagnosis[train_index])
    
    # get test x and test y
    test_x <- beta_methyl[,!train_index]
    test_y <- as.numeric(clin$age_diagnosis[!train_index])
    
    # check sizes
    stopifnot(dim(train_x)[2] == length(train_y))
    stopifnot(dim(test_x)[2] == length(test_y))
    
    # check name length
    stopifnot(dim(train_x)[1] == length(feature_names))
    
    # put into data lists for model 
    train_dat <- list(x = train_x, y = train_y, feature_names = feature_names)
    test_dat <- list(x = test_x, y = test_y,  feature_names = feature_names)
    
    # This step just computes the  scores for each feature and creates a training object
    train_object <- superpc.train(train_dat, type="regression")
    
    # cross-validate the model to get optimal threshold to be used in next function
    cv_object[[i]] <-superpc.cv(train_object, train_dat, n.components = 5, n.fold = 5)
    superpc.plotcv(cv_object[[i]])
    
    # round cv_object scores 
    lrst <- as.data.frame(round(cv_object[[i]]$scor, 2))
    
    # get number of features used at each threshold 
    used_feats <- as.data.frame(cv_object[[i]]$nonzero)
    names(used_feats) <- 'used_feats'
    
    # get the max from each column (pc with highest lrst)
    lrst_scores <- apply(lrst, 2, function(x) max(x, na.rm = T))
    
    # combine lrst_scores with threshold data
    thresh_dat <- as.data.frame(cbind(lrst_score = lrst_scores, thresh = cv_object[[i]]$thresholds, used_features = used_feats))
    
    # get index 
    thresh_index <- which(thresh_dat$lrst_score == max(thresh_dat$lrst_score))
    
    # get threshold with highest lrst score
    opt_thresh[[i]] <- thresh_dat$thresh[thresh_index]
    
    # get optimal feats
    opt_feats[[i]] <- thresh_dat$used_feats[thresh_index]
    
    pred_pca <- superpc.predict(train_object, 
                                train_dat, 
                                test_dat, 
                                threshold = opt_thresh[[i]], 
                                n.components=3, 
                                prediction.type="continuous")$v.pred
    
    
    # will return correlation on test data with pred_pca predictors 
    test_results[[i]] <- superpc_fit_lm_cv(y = test_dat$y, 
                                           score = pred_pca,
                                           cv = cv_type)[[1]]
    
    # pring status 
    print(paste0('finished with fold', ' ', i))
    
  }
  
  
  return(test_results)
  
}

mod_results <- train_test_pca(cases = beta_methyl, 
                              clin = clin_data, 
                              fold_numbers = fold_vec, 
                              feature_names = feats,
                              cv_type = 'auto')



# change pred to nothing if doing surv
saveRDS(mod_results, paste0('/hpf/largeprojects/agoldenb/ben/Projects/LFS/Scripts/predict_age/Results/reg_results/train_test', '_' , seed_num,'_', type, '_', method, '.rda'))

