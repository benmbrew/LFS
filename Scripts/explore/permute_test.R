####### This script will look more in depth at models - class balance, permute outcome

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

registerDoParallel(1)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')
scripts_folder <- paste0(project_folder, '/Scripts')

##########
# load data 
##########
# cases - raw, swan, quan, funnorm, and clin
load(paste0(model_data, '/model_data_cases.RData'))

# bh_features
load(paste0(model_data, '/bh_features.RData'))

##########
# function taht gets data you want and removes other
##########
removeDat <- function(keep)
{
  data_set <- c('raw', 'quan', 'swan', 'funnorm')
  
  remove <- data_set[data_set != keep]
  
  # remove unwated
  rm(list=ls(pattern=remove[[1]], envir = .GlobalEnv), envir = .GlobalEnv)
  rm(list=ls(pattern=remove[[2]], envir = .GlobalEnv), envir = .GlobalEnv)
  rm(list=ls(pattern=remove[[3]], envir = .GlobalEnv), envir = .GlobalEnv)
  
}

removeDat(keep = 'funnorm')

##########
# get distribtuion of age for each seed
##########
# subset data by not na in age of diagnosis and mut
dat <- beta_funnorm[!is.na(beta_funnorm$age_diagnosis),]
dat <- dat[dat$p53_germline == 'Mut',]

# make histograms for seed 1-10
for (seed in 1:10) {
  set.seed(seed)
  train_index <- sample(nrow(dat), nrow(dat) *.7, replace = F)
  
  hist(dat[train_index, 'age_diagnosis'], main = seed)
  hist(dat[-train_index, 'age_diagnosis'], main = seed)
  print(seed)
}

##########
# get residuals
##########
# funnorm_bal_p53
# bh_data <- beta_funnorm_bal_p53_features
# data <- dat
# Function that takes model model_data and obstains residual features based on regressin each gene/probe on age of sample collection
getResidual <- function(data, bh_data) 
{
  # subet data by bh
  bh_features <- bh_data$probe
  
  # get features
  feature_names <- colnames(data)[6:ncol(data)]
  
  # intersect the two
  model_features <- intersect(bh_features, feature_names)
  
  # rearrange
  data <- data[, c('age_diagnosis','age_sample_collection', model_features)]
  
  # complete cases 
  data <- data[complete.cases(data),]
  
  resid <- list()
    
    for (i in 3:ncol(data)) {
      
      temp_response <- data[, i]
      temp_var <- data$age_sample_collection
      
      resid[[i]] <- lm(temp_response ~ temp_var)$residuals
      
      print(i)
      
    }
  
  resid_data <- as.data.frame(do.call('cbind', resid))
  data <- cbind(data$age_diagnosis, data$age_sample_collection, resid_data)
  colnames(data) <- c('age_diagnosis', 'age_sample_collection', model_features)
  
  return(data)

  
  return(data_list)
}

dat_resid <- getResidual(dat, beta_funnorm_p53_union_features)

##########
# function that predicts regression, gets class counts, and permutes results 
##########

rfPredictReg <- function(model_data,
                         cutoff,
                         iterations) 
{
  
  model <- list()
  best_features <- list()
  importance <- list()
  test.predictions <- list()
  test.ground_truth <- list()
  test_acc <- list()
  test_stats  <- list()

  
  dims <- dim(model_data)
  selected_features <- names(model_data[, 3:ncol(model_data)])
  
  
  for (i in 1:iterations) {
    
    set.seed(i)
    
    train_index <- sample(nrow(model_data), nrow(model_data) *cutoff, replace = F)
    
    NFOLDS <- 2
    fitControl <- trainControl( 
      method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
      number = min(10, NFOLDS),      
      repeats = 1,
      allowParallel = TRUE
    )
    
    # get y
    y <- model_data$age_diagnosis[train_index]
    
    # permute y
    y <- sample(y, length(y))
    
    # mtry: Number of variables randomly sampled as candidates at each split.
    # ntree: Number of trees to grow.
    mtry <- sqrt(ncol(model_data[train_index, selected_features]))
    tunegrid <- expand.grid(.mtry=mtry)
    
    model[[i]] <- train(x = model_data[train_index, selected_features]
                        , y = y
                        , method = "rf"
                        , trControl = fitControl
                        , tuneGrid = tunegrid
                        , importance = T
                        , verbose = FALSE)
    
    temp <- varImp(model[[i]])[[1]]
    importance[[i]] <- cbind(rownames(temp), temp$Overall)
    
    
    test.predictions[[i]] <- predict(model[[i]] 
                                     , newdata = model_data[-train_index, selected_features])
    
    
    test.ground_truth[[i]] <- model_data$age_diagnosis[-train_index]
    
    print(i)
    
  }
  
  return(list(test.predictions, test.ground_truth, test_acc, 
              test_stats, model, importance, dims))
  
}

result_list <- rfPredictReg(model_data = dat_resid,
                            cutoff = 0.7,
                            iterations = 5)

plot(unlist(result_list[[1]]), unlist(result_list[[2]]), 
     xlab = 'Predictions', ylab = 'Real Age', main = 'Permuted outcome')

