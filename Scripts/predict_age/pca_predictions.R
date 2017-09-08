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
k = 5
seed_num <- 1

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
                         seed_num = seed_num, k = k)

# get beta matrix, clinical variables, fold vector, and feature names for the supervised PCA
beta_methyl <- mod_data[[1]] # a transposed (p x n) methylation matrix
clin_data <- mod_data[[2]] # c
fold_vec <- mod_data[[3]]
feats <- mod_data[[4]]


# beta_methyl <- beta_methyl[1:1000, ]
cases <- beta_methyl
clin <- clin_data
fold_numbers <- fold_vec
feature_names <- feats[1:1000]
i = 1

trainTest <- function(cases, 
                      feature_names,
                      clin,
                      controls,
                      valid,
                      probe_start,
                      k) {
  
  
  # # creat features
  # feature_names <- colnames(cases)[probe_start:(ncol(cases) - 1)]
  # stopifnot(all(grepl('cg', feature_names)))
  
  # list to store results
  bh_feat <- list()
  model_results <- list()
  bh_dim <- list()
  
  # now write forloop to 
  for (i in 1:k) {
    
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
    
    # put into data lists for model 
    train_dat <- list(x = train_x, y = train_y, feature_names = feature_names)
    test_dat <- list(x = test_x, y = test_y,  feature_names = feature_names)
    
    # This step just computes the  scores for each feature and creates a training object
    train.obj<- superpc.train(train_dat, type="regression")
    
    # cross-validate the model
    cv.obj <-superpc.cv(train.obj, train_dat)
    superpc.plotcv(cv.obj)
    
    # # here we have the luxury of  test data, so we can compute the  likelihood ratio statistic
    lrtest.obj <- superpc.lrtest.curv(train.obj, train_dat, test_dat)
    # 
    superpc.plot.lrtest(lrtest.obj)
    # 
    
    fit.cts <- superpc.predict(train.obj, 
                               train_dat, 
                               test_dat, 
                               threshold =1.5, 
                               n.components=3, 
                               prediction.type="continuous")
    
    
    superpc.fit.to.outcome(train.obj, test_dat, fit.cts$v.pred)
    
    
  }
  
  return(list(model_results, bh_dim))
  
}

mod_results <- trainTest(cases = betaCases,
                         controls = betaControls,
                         valid = betaValid,
                         k = k)

# change pred to nothing if doing surv
saveRDS(mod_results, paste0('/hpf/largeprojects/agoldenb/ben/Projects/LFS/Scripts/predict_age/Results/reg_results/train_test', '_' , seed_num,'_', type, '_', method, '.rda'))

