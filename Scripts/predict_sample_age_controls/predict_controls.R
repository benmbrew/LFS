#######################################
# This script will predict age of sample collection from 
# control data - that is samples with mut status, but no cancer diagnosis. 
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(doParallel)
library(nnet)

registerDoParallel(1)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
scripts_folder <- paste0(project_folder, '/Scripts')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
idat_data <- paste0(methyl_data, '/raw_files')
clin_data <- paste0(data_folder, '/clin_data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')
model_data <- paste0(data_folder, '/model_data')


##########
# load controls dat - ex:beta_raw_transformation, beta_raw_controls_overlap
##########
load(paste0(model_data, '/controls.RData'))

##########
# load bumphunter features
##########
load(paste0(model_data, '/bh_features_idat.RData'))

##########
# function to get bumphunter features
##########
getFeatures <- function(data, bh_features) {
  
  bh_features <- bh_features$probe
  model_features <- colnames(data)
  bh_intersect <- intersect(bh_features, model_features)
  data <- data[, c('age_sample_collection', bh_intersect)]
  return(data)
}


##########
# function to get random features
##########
getRand <- function(data) {
  features <- colnames(data)[5:ncol(data)]
  rand_features <- sample(features, 900, replace = F)
  data <-  data[, c('age_sample_collection', rand_features)]
  return(data)
}

##########
# make numeric 
##########

makeNum <- function(model_data) {
  
  model_data[, 1:ncol(model_data)] <- apply(model_data[, 1:ncol(model_data)], 2, function(x) as.numeric(as.character(x)))
  
  return(model_data)
}

#########################################################################################
# Prepare data - for each beta, get features on the transformed data and on original data

##########
# beta raw - transformation
##########
beta_raw_tran_global_unbal_bh <- getFeatures(beta_raw_transformation, beta_raw_global_unbal_features)
beta_raw_tran_cancer_unbal_bh <- getFeatures(beta_raw_transformation, beta_raw_cancer_unbal_features)
beta_raw_tran_global_bal_bh <- getFeatures(beta_raw_transformation, beta_raw_global_bal_features)
beta_raw_tran_cancer_bal_bh <- getFeatures(beta_raw_transformation, beta_raw_cancer_bal_features)
beta_raw_tran_rand <- getRand(beta_raw_transformation)

beta_raw_tran_global_unbal_bh <- makeNum(beta_raw_tran_global_unbal_bh)
beta_raw_tran_cancer_unbal_bh <- makeNum(beta_raw_tran_cancer_unbal_bh)
beta_raw_tran_global_bal_bh <- makeNum(beta_raw_tran_global_bal_bh)
beta_raw_tran_cancer_bal_bh <- makeNum(beta_raw_tran_cancer_bal_bh)
beta_raw_tran_rand <- makeNum(beta_raw_tran_rand)


##########
# beta raw - original
##########
beta_raw_orig_global_unbal_bh <- getFeatures(beta_raw_controls_overlap, beta_raw_global_unbal_features)
beta_raw_orig_cancer_unbal_bh <- getFeatures(beta_raw_controls_overlap, beta_raw_cancer_unbal_features)
beta_raw_orig_global_bal_bh <- getFeatures(beta_raw_controls_overlap, beta_raw_global_bal_features)
beta_raw_orig_cancer_bal_bh <- getFeatures(beta_raw_controls_overlap, beta_raw_cancer_bal_features)
beta_raw_orig_rand <- getRand(beta_raw_controls_overlap)

beta_raw_orig_global_unbal_bh <- makeNum(beta_raw_orig_global_unbal_bh)
beta_raw_orig_cancer_unbal_bh <- makeNum(beta_raw_orig_cancer_unbal_bh)
beta_raw_orig_global_bal_bh <- makeNum(beta_raw_orig_global_bal_bh)
beta_raw_orig_cancer_bal_bh <- makeNum(beta_raw_orig_cancer_bal_bh)
beta_raw_orig_rand <- makeNum(beta_raw_orig_rand)

rm(beta_raw_transformation, beta_raw_controls_overlap)

##############################################################################################

##########
# beta swan - transformation
##########
beta_swan_tran_global_unbal_bh <- getFeatures(beta_swan_transformation, beta_swan_global_unbal_features)
beta_swan_tran_cancer_unbal_bh <- getFeatures(beta_swan_transformation, beta_swan_cancer_unbal_features)
beta_swan_tran_global_bal_bh <- getFeatures(beta_swan_transformation, beta_swan_global_bal_features)
beta_swan_tran_cancer_bal_bh <- getFeatures(beta_swan_transformation, beta_swan_cancer_bal_features)
beta_swan_tran_rand <- getRand(beta_swan_transformation)

beta_swan_tran_global_unbal_bh <- makeNum(beta_swan_tran_global_unbal_bh)
beta_swan_tran_cancer_unbal_bh <- makeNum(beta_swan_tran_cancer_unbal_bh)
beta_swan_tran_global_bal_bh <- makeNum(beta_swan_tran_global_bal_bh)
beta_swan_tran_cancer_bal_bh <- makeNum(beta_swan_tran_cancer_bal_bh)
beta_swan_tran_rand <- makeNum(beta_swan_tran_rand)


##########
# beta swan - original
##########
beta_swan_orig_global_unbal_bh <- getFeatures(beta_swan_controls_overlap, beta_swan_global_unbal_features)
beta_swan_orig_cancer_unbal_bh <- getFeatures(beta_swan_controls_overlap, beta_swan_cancer_unbal_features)
beta_swan_orig_global_bal_bh <- getFeatures(beta_swan_controls_overlap, beta_swan_global_bal_features)
beta_swan_orig_cancer_bal_bh <- getFeatures(beta_swan_controls_overlap, beta_swan_cancer_bal_features)
beta_swan_orig_rand <- getRand(beta_swan_controls_overlap)

beta_swan_orig_global_unbal_bh <- makeNum(beta_swan_orig_global_unbal_bh)
beta_swan_orig_cancer_unbal_bh <- makeNum(beta_swan_orig_cancer_unbal_bh)
beta_swan_orig_global_bal_bh <- makeNum(beta_swan_orig_global_bal_bh)
beta_swan_orig_cancer_bal_bh <- makeNum(beta_swan_orig_cancer_bal_bh)
beta_swan_orig_rand <- makeNum(beta_swan_orig_rand)

rm(beta_swan_transformation, beta_swan_controls_overlap)

##############################################################################################

##########
# beta quan - transformation
##########
beta_quan_tran_global_unbal_bh <- getFeatures(beta_quan_transformation, beta_quan_global_unbal_features)
beta_quan_tran_cancer_unbal_bh <- getFeatures(beta_quan_transformation, beta_quan_cancer_unbal_features)
beta_quan_tran_global_bal_bh <- getFeatures(beta_quan_transformation, beta_quan_global_bal_features)
beta_quan_tran_cancer_bal_bh <- getFeatures(beta_quan_transformation, beta_quan_cancer_bal_features)
beta_quan_tran_rand <- getRand(beta_quan_transformation)

beta_quan_tran_global_unbal_bh <- makeNum(beta_quan_tran_global_unbal_bh)
beta_quan_tran_cancer_unbal_bh <- makeNum(beta_quan_tran_cancer_unbal_bh)
beta_quan_tran_global_bal_bh <- makeNum(beta_quan_tran_global_bal_bh)
beta_quan_tran_cancer_bal_bh <- makeNum(beta_quan_tran_cancer_bal_bh)
beta_quan_tran_rand <- makeNum(beta_quan_tran_rand)


##########
# beta quan - original
##########
beta_quan_orig_global_unbal_bh <- getFeatures(beta_quan_controls_overlap, beta_quan_global_unbal_features)
beta_quan_orig_cancer_unbal_bh <- getFeatures(beta_quan_controls_overlap, beta_quan_cancer_unbal_features)
beta_quan_orig_global_bal_bh <- getFeatures(beta_quan_controls_overlap, beta_quan_global_bal_features)
beta_quan_orig_cancer_bal_bh <- getFeatures(beta_quan_controls_overlap, beta_quan_cancer_bal_features)
beta_quan_orig_rand <- getRand(beta_quan_controls_overlap)

beta_quan_orig_global_unbal_bh <- makeNum(beta_quan_orig_global_unbal_bh)
beta_quan_orig_cancer_unbal_bh <- makeNum(beta_quan_orig_cancer_unbal_bh)
beta_quan_orig_global_bal_bh <- makeNum(beta_quan_orig_global_bal_bh)
beta_quan_orig_cancer_bal_bh <- makeNum(beta_quan_orig_cancer_bal_bh)
beta_quan_orig_rand <- makeNum(beta_quan_orig_rand)

rm(beta_quan_transformation, beta_quan_controls_overlap)

##############################################################################################

##########
# beta funnorm - transformation
##########
beta_funnorm_tran_global_unbal_bh <- getFeatures(beta_funnorm_transformation, beta_funnorm_global_unbal_features)
beta_funnorm_tran_cancer_unbal_bh <- getFeatures(beta_funnorm_transformation, beta_funnorm_cancer_unbal_features)
beta_funnorm_tran_global_bal_bh <- getFeatures(beta_funnorm_transformation, beta_funnorm_global_bal_features)
beta_funnorm_tran_cancer_bal_bh <- getFeatures(beta_funnorm_transformation, beta_funnorm_cancer_bal_features)
beta_funnorm_tran_rand <- getRand(beta_funnorm_transformation)

beta_funnorm_tran_global_unbal_bh <- makeNum(beta_funnorm_tran_global_unbal_bh)
beta_funnorm_tran_cancer_unbal_bh <- makeNum(beta_funnorm_tran_cancer_unbal_bh)
beta_funnorm_tran_global_bal_bh <- makeNum(beta_funnorm_tran_global_bal_bh)
beta_funnorm_tran_cancer_bal_bh <- makeNum(beta_funnorm_tran_cancer_bal_bh)
beta_funnorm_tran_rand <- makeNum(beta_funnorm_tran_rand)


##########
# beta funnorm - original
##########
beta_funnorm_orig_global_unbal_bh <- getFeatures(beta_funnorm_controls_overlap, beta_funnorm_global_unbal_features)
beta_funnorm_orig_cancer_unbal_bh <- getFeatures(beta_funnorm_controls_overlap, beta_funnorm_cancer_unbal_features)
beta_funnorm_orig_global_bal_bh <- getFeatures(beta_funnorm_controls_overlap, beta_funnorm_global_bal_features)
beta_funnorm_orig_cancer_bal_bh <- getFeatures(beta_funnorm_controls_overlap, beta_funnorm_cancer_bal_features)
beta_funnorm_orig_rand <- getRand(beta_funnorm_controls_overlap)

beta_funnorm_orig_global_unbal_bh <- makeNum(beta_funnorm_orig_global_unbal_bh)
beta_funnorm_orig_cancer_unbal_bh <- makeNum(beta_funnorm_orig_cancer_unbal_bh)
beta_funnorm_orig_global_bal_bh <- makeNum(beta_funnorm_orig_global_bal_bh)
beta_funnorm_orig_cancer_bal_bh <- makeNum(beta_funnorm_orig_cancer_bal_bh)
beta_funnorm_orig_rand <- makeNum(beta_funnorm_orig_rand)

rm(beta_funnorm_transformation, beta_funnorm_controls_overlap)

##############################################################################################

##########
# function that runs random forest as a classification
##########
rfPredictFac <- function(model_data,
                         cutoff,
                         iterations) {
  
  model <- list()
  best_features <- list()
  importance <- list()
  train.predictions <- list()
  test.predictions <- list()
  train.ground_truth <- list()
  test.ground_truth <- list()
  test_acc <- list()
  test_stats  <- list()
  
  
  dims <- dim(model_data)
  
  selected_features <- names(model_data[, 2:ncol(model_data)])
  
  model_data$age_sample_fac <- ifelse(model_data$age_sample_collection > 140, 1,2)
  model_data$age_sample_collection <- NULL
  
  for (i in 1:iterations) {
    
    set.seed(i)
    
    train_index <- sample(nrow(model_data), nrow(model_data) *cutoff, replace = F)
    
    # determines how you train the model.
    NFOLDS <- 2
    fitControl <- trainControl( 
      method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
      number = min(10, NFOLDS),
      classProbs = TRUE,
      repeats = 1,
      allowParallel = TRUE,
      summaryFunction = twoClassSummary
      
    )
    
    y = make.names(as.factor(model_data$age_sample_fac[train_index]))
    
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
                                     , newdata = model_data[-train_index, selected_features]
                                     , type = "prob")
    
    train.predictions[[i]] <- predict(model[[i]] 
                                      , newdata = model_data[train_index, selected_features]
                                      ,type = "prob")
    
    
    
    
    
    train.ground_truth[[i]] <- as.factor(make.names(model_data$age_sample_fac[train_index]))
    test.ground_truth[[i]] <- as.factor(make.names(model_data$age_sample_fac[-train_index]))
    
    
    # For age of diagnosis
    # Accuracy
    test_acc[[i]] <- sum(levels(test.ground_truth[[i]])[apply(test.predictions[[i]], 1, which.is.max)] == test.ground_truth[[i]]) / dim(test.predictions[[i]])[1]
    # Confustion Matrix
    test_stats[[i]] <- confusionMatrix(levels(test.ground_truth[[i]])[apply(test.predictions[[i]], 1, which.is.max)], test.ground_truth[[i]])
    # print(rf.test_stats
    
    
    print(i)
    
  }
  
  return(list(train.predictions, test.predictions, train.ground_truth, test.ground_truth, test_acc, test_stats, model, importance, dims))
  
}

###############################
# Apply function to beta_funnorm 
###############################
# function that runs random forest as a regression

rfPredictReg <- function(model_data,
                         cutoff,
                         features,
                         iterations) {
  
  model <- list()
  best_features <- list()
  importance <- list()
  train.mse <- list()
  test.mse <- list()
  train.predictions <- list()
  test.predictions <- list()
  train.ground_truth <- list()
  test.ground_truth <- list()
  test_acc <- list()
  test_stats  <- list()
  
  
  dims <- dim(model_data)
  selected_features <- names(model_data[, 2:ncol(model_data)])
  
  
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
    
    y <- model_data$age_sample_collection[train_index]
    
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
    
    train.predictions[[i]] <- predict(model[[i]] 
                                      , newdata = model_data[train_index, selected_features])
    
    
    
    
    
    train.ground_truth[[i]] <- model_data$age_sample_collection[train_index]
    test.ground_truth[[i]] <- model_data$age_sample_collection[-train_index]
    
    train.mse[[i]] <- rmse(unlist(train.predictions[[i]]), unlist(train.ground_truth[[i]]))
    test.mse[[i]] <- rmse(unlist(test.predictions[[i]]), unlist(test.ground_truth[[i]]))
    
    
    print(i)
    
  }
  
  return(list(train.mse, test.mse,  train.predictions, test.predictions, train.ground_truth, test.ground_truth,
              test_acc, test_stats,  model, importance, dims))
  
}



# Function that takes results list from regression and plots predictions against ground truth

plotModel <- function(result_list,
                      main,
                      xlim,
                      ylim) {
  
  # plot predictions against ground truth
  plot(unlist(result_list[[4]]), unlist(result_list[[6]]), 
       xlim = xlim,
       ylim = ylim,
       bty = 'n',
       col = adjustcolor('blue', alpha.f = 0.6),
       pch = 16,
       xlab = 'Predictions',
       ylab = 'Real Age of Diagnosis',
       main = main)
  abline(0,1, lty = 3)
  corr <- round(cor(unlist(result_list[[4]]), unlist(result_list[[6]])), 2)
  legend("topleft", legend = paste0('correlation = ', corr), cex = 1, bty = 'n')
  legend("bottomright", legend = paste0('# dims = ', result_list[[11]]), cex = 1, bty = 'n')
  
  
  
}

conMatrix <- function(results) {
  
  # test acc for age of diagnosis
  acc_age <- mean(unlist(results[[5]]))
  
  # confustion matrix age of diagnosis 10
  iterations <- 10
  temp <- list()
  for (i in 1:10){
    temp[[i]] <- results[[6]][[i]]$table
  }
  mat <- unlist(temp)
  new_mat <- matrix(, 2, 2)
  
  mat_index <- seq(1, length(mat), 4)
  
  new_mat[1,1] <- sum(mat[mat_index])/iterations
  new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
  new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
  new_mat[2,2] <- sum(mat[mat_index + 3])/iterations
  
  return(list(new_mat, acc_age))
  
}




##################################################################################################
# apply functions to all preprocessing methods and transformed and orig
########################################

##########
# BETA RAW - transformed
##########

##########
# global unbal
##########
beta_raw_tran_global_unbal_model_fac <- rfPredictFac(beta_raw_tran_global_unbal_bh, 
                                                cutoff = .7, 
                                                iterations = 10)

conMatrix(beta_raw_tran_global_unbal_model_fac)

beta_raw_tran_global_unbal_model_reg <- rfPredictReg(beta_raw_tran_global_unbal_bh,
                                                cutoff = .7,
                                                iterations = 10)

plotModel(beta_raw_tran_global_unbal_model_reg, main = 'raw prepocessing, tran, unbal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# global bal
##########

beta_raw_tran_global_bal_model_fac <- rfPredictFac(beta_raw_tran_global_bal_bh, 
                                              cutoff = .7, 
                                              iterations = 10)

conMatrix(beta_raw_tran_global_bal_model_fac)

beta_raw_tran_global_bal_model_reg <- rfPredictReg(beta_raw_tran_global_bal_bh,
                                              cutoff = .7,
                                              iterations = 10)

plotModel(beta_raw_tran_global_bal_model_reg, main = 'raw prepocessing, tran, bal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# cancer unbal
##########
beta_raw_tran_cancer_unbal_model_fac <- rfPredictFac(beta_raw_tran_cancer_unbal_bh, 
                                                cutoff = .7, 
                                                iterations = 10)

conMatrix(beta_raw_tran_cancer_unbal_model_fac)

beta_raw_tran_cancer_unbal_model_reg <- rfPredictReg(beta_raw_tran_cancer_unbal_bh,
                                                cutoff = .7,
                                                iterations = 10)

plotModel(beta_raw_tran_cancer_unbal_model_reg, main = 'raw prepocessing, tran, unbal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))


##########
# cancer bal
##########
beta_raw_tran_cancer_bal_model_fac <- rfPredictFac(beta_raw_tran_cancer_bal_bh, 
                                              cutoff = .7, 
                                              iterations = 10)

conMatrix(beta_raw_tran_cancer_bal_model_fac)

beta_raw_tran_cancer_bal_model_reg <- rfPredictReg(beta_raw_tran_cancer_bal_bh,
                                              cutoff = .7,
                                              iterations = 10)

plotModel(beta_raw_tran_cancer_bal_model_reg, main = 'raw prepocessing, tran, bal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########  
# random
##########
beta_raw_tran_rand_fac <- rfPredictFac(beta_raw_tran_rand, 
                                  cutoff = .7, 
                                  iterations = 10)

conMatrix(beta_raw_tran_rand_fac)

beta_raw_tran_rand_reg <- rfPredictReg(beta_raw_tran_rand,
                                  cutoff = .7,
                                  iterations = 10)

plotModel(beta_raw_tran_rand_reg, main = 'raw random, tran',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

######################################################################################################

##########
# BETA swan - transformed
##########

##########
# global unbal
##########
beta_swan_tran_global_unbal_model_fac <- rfPredictFac(beta_swan_tran_global_unbal_bh, 
                                                     cutoff = .7, 
                                                     iterations = 10)

conMatrix(beta_swan_tran_global_unbal_model_fac)

beta_swan_tran_global_unbal_model_reg <- rfPredictReg(beta_swan_tran_global_unbal_bh,
                                                     cutoff = .7,
                                                     iterations = 10)

plotModel(beta_swan_tran_global_unbal_model_reg, main = 'swan prepocessing, tran, unbal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# global bal
##########

beta_swan_tran_global_bal_model_fac <- rfPredictFac(beta_swan_tran_global_bal_bh, 
                                                   cutoff = .7, 
                                                   iterations = 10)

conMatrix(beta_swan_tran_global_bal_model_fac)

beta_swan_tran_global_bal_model_reg <- rfPredictReg(beta_swan_tran_global_bal_bh,
                                                   cutoff = .7,
                                                   iterations = 10)

plotModel(beta_swan_tran_global_bal_model_reg, main = 'swan prepocessing, tran, bal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# cancer unbal
##########
beta_swan_tran_cancer_unbal_model_fac <- rfPredictFac(beta_swan_tran_cancer_unbal_bh, 
                                                     cutoff = .7, 
                                                     iterations = 10)

conMatrix(beta_swan_tran_cancer_unbal_model_fac)

beta_swan_tran_cancer_unbal_model_reg <- rfPredictReg(beta_swan_tran_cancer_unbal_bh,
                                                     cutoff = .7,
                                                     iterations = 10)

plotModel(beta_swan_tran_cancer_unbal_model_reg, main = 'swan prepocessing, tran, unbal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))


##########
# cancer bal
##########
beta_swan_tran_cancer_bal_model_fac <- rfPredictFac(beta_swan_tran_cancer_bal_bh, 
                                                   cutoff = .7, 
                                                   iterations = 10)

conMatrix(beta_swan_tran_cancer_bal_model_fac)

beta_swan_tran_cancer_bal_model_reg <- rfPredictReg(beta_swan_tran_cancer_bal_bh,
                                                   cutoff = .7,
                                                   iterations = 10)

plotModel(beta_swan_tran_cancer_bal_model_reg, main = 'swan prepocessing, tran, bal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# random
##########
beta_swan_tran_rand_fac <- rfPredictFac(beta_swan_tran_rand, 
                                       cutoff = .7, 
                                       iterations = 10)

conMatrix(beta_swan_tran_rand_fac)

beta_swan_tran_rand_reg <- rfPredictReg(beta_swan_tran_rand,
                                       cutoff = .7,
                                       iterations = 10)

plotModel(beta_swan_tran_rand_reg, main = 'swan random, tran',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

######################################################################################################

##########
# BETA quan - transformed
##########

##########
# global unbal
##########
beta_quan_tran_global_unbal_model_fac <- rfPredictFac(beta_quan_tran_global_unbal_bh, 
                                                     cutoff = .7, 
                                                     iterations = 10)

conMatrix(beta_quan_tran_global_unbal_model_fac)

beta_quan_tran_global_unbal_model_reg <- rfPredictReg(beta_quan_tran_global_unbal_bh,
                                                     cutoff = .7,
                                                     iterations = 10)

plotModel(beta_quan_tran_global_unbal_model_reg, main = 'quan prepocessing, tran, unbal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# global bal
##########

beta_quan_tran_global_bal_model_fac <- rfPredictFac(beta_quan_tran_global_bal_bh, 
                                                   cutoff = .7, 
                                                   iterations = 10)

conMatrix(beta_quan_tran_global_bal_model_fac)

beta_quan_tran_global_bal_model_reg <- rfPredictReg(beta_quan_tran_global_bal_bh,
                                                   cutoff = .7,
                                                   iterations = 10)

plotModel(beta_quan_tran_global_bal_model_reg, main = 'quan prepocessing, tran, bal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# cancer unbal
##########
beta_quan_tran_cancer_unbal_model_fac <- rfPredictFac(beta_quan_tran_cancer_unbal_bh, 
                                                     cutoff = .7, 
                                                     iterations = 10)

conMatrix(beta_quan_tran_cancer_unbal_model_fac)

beta_quan_tran_cancer_unbal_model_reg <- rfPredictReg(beta_quan_tran_cancer_unbal_bh,
                                                     cutoff = .7,
                                                     iterations = 10)

plotModel(beta_quan_tran_cancer_unbal_model_reg, main = 'quan prepocessing, tran, unbal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))


##########
# cancer bal
##########
beta_quan_tran_cancer_bal_model_fac <- rfPredictFac(beta_quan_tran_cancer_bal_bh, 
                                                   cutoff = .7, 
                                                   iterations = 10)

conMatrix(beta_quan_tran_cancer_bal_model_fac)

beta_quan_tran_cancer_bal_model_reg <- rfPredictReg(beta_quan_tran_cancer_bal_bh,
                                                   cutoff = .7,
                                                   iterations = 10)

plotModel(beta_quan_tran_cancer_bal_model_reg, main = 'quan prepocessing, tran, bal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# random
##########
beta_quan_tran_rand_fac <- rfPredictFac(beta_quan_tran_rand, 
                                       cutoff = .7, 
                                       iterations = 10)

conMatrix(beta_quan_tran_rand_fac)

beta_quan_tran_rand_reg <- rfPredictReg(beta_quan_tran_rand,
                                       cutoff = .7,
                                       iterations = 10)

plotModel(beta_quan_tran_rand_reg, main = 'quan random, tran',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

######################################################################################################

##########
# BETA funnorm - transformed
##########

##########
# global unbal
##########
beta_funnorm_tran_global_unbal_model_fac <- rfPredictFac(beta_funnorm_tran_global_unbal_bh, 
                                                     cutoff = .7, 
                                                     iterations = 10)

conMatrix(beta_funnorm_tran_global_unbal_model_fac)

beta_funnorm_tran_global_unbal_model_reg <- rfPredictReg(beta_funnorm_tran_global_unbal_bh,
                                                     cutoff = .7,
                                                     iterations = 10)

plotModel(beta_funnorm_tran_global_unbal_model_reg, main = 'funnorm prepocessing, tran, unbal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# global bal
##########

beta_funnorm_tran_global_bal_model_fac <- rfPredictFac(beta_funnorm_tran_global_bal_bh, 
                                                   cutoff = .7, 
                                                   iterations = 10)

conMatrix(beta_funnorm_tran_global_bal_model_fac)

beta_funnorm_tran_global_bal_model_reg <- rfPredictReg(beta_funnorm_tran_global_bal_bh,
                                                   cutoff = .7,
                                                   iterations = 10)

plotModel(beta_funnorm_tran_global_bal_model_reg, main = 'funnorm prepocessing, tran, bal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# cancer unbal
##########
beta_funnorm_tran_cancer_unbal_model_fac <- rfPredictFac(beta_funnorm_tran_cancer_unbal_bh, 
                                                     cutoff = .7, 
                                                     iterations = 10)

conMatrix(beta_funnorm_tran_cancer_unbal_model_fac)

beta_funnorm_tran_cancer_unbal_model_reg <- rfPredictReg(beta_funnorm_tran_cancer_unbal_bh,
                                                     cutoff = .7,
                                                     iterations = 10)

plotModel(beta_funnorm_tran_cancer_unbal_model_reg, main = 'funnorm prepocessing, tran, unbal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))


##########
# cancer bal
##########
beta_funnorm_tran_cancer_bal_model_fac <- rfPredictFac(beta_funnorm_tran_cancer_bal_bh, 
                                                   cutoff = .7, 
                                                   iterations = 10)

conMatrix(beta_funnorm_tran_cancer_bal_model_fac)

beta_funnorm_tran_cancer_bal_model_reg <- rfPredictReg(beta_funnorm_tran_cancer_bal_bh,
                                                   cutoff = .7,
                                                   iterations = 10)

plotModel(beta_funnorm_tran_cancer_bal_model_reg, main = 'funnorm prepocessing, tran, bal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# random
##########
beta_funnorm_tran_rand_fac <- rfPredictFac(beta_funnorm_tran_rand, 
                                       cutoff = .7, 
                                       iterations = 10)

conMatrix(beta_funnorm_tran_rand_fac)

beta_funnorm_tran_rand_reg <- rfPredictReg(beta_funnorm_tran_rand,
                                       cutoff = .7,
                                       iterations = 10)

plotModel(beta_funnorm_tran_rand_reg, main = 'funnorm random, tran',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

###############################################################################################

##########
# BETA RAW - original
##########

##########
# global unbal
##########
beta_raw_orig_global_unbal_model_fac <- rfPredictFac(beta_raw_orig_global_unbal_bh, 
                                                     cutoff = .7, 
                                                     iterations = 10)

conMatrix(beta_raw_orig_global_unbal_model_fac)

beta_raw_orig_global_unbal_model_reg <- rfPredictReg(beta_raw_orig_global_unbal_bh,
                                                     cutoff = .7,
                                                     iterations = 10)

plotModel(beta_raw_orig_global_unbal_model_reg, main = 'raw prepocessing, orig, unbal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# global bal
##########

beta_raw_orig_global_bal_model_fac <- rfPredictFac(beta_raw_orig_global_bal_bh, 
                                                   cutoff = .7, 
                                                   iterations = 10)

conMatrix(beta_raw_orig_global_bal_model_fac)

beta_raw_orig_global_bal_model_reg <- rfPredictReg(beta_raw_orig_global_bal_bh,
                                                   cutoff = .7,
                                                   iterations = 10)

plotModel(beta_raw_orig_global_bal_model_reg, main = 'raw prepocessing, orig, bal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# cancer unbal
##########
beta_raw_orig_cancer_unbal_model_fac <- rfPredictFac(beta_raw_orig_cancer_unbal_bh, 
                                                     cutoff = .7, 
                                                     iterations = 10)

conMatrix(beta_raw_orig_cancer_unbal_model_fac)

beta_raw_orig_cancer_unbal_model_reg <- rfPredictReg(beta_raw_orig_cancer_unbal_bh,
                                                     cutoff = .7,
                                                     iterations = 10)

plotModel(beta_raw_orig_cancer_unbal_model_reg, main = 'raw prepocessing, orig, unbal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))


##########
# cancer bal
##########
beta_raw_orig_cancer_bal_model_fac <- rfPredictFac(beta_raw_orig_cancer_bal_bh, 
                                                   cutoff = .7, 
                                                   iterations = 10)

conMatrix(beta_raw_orig_cancer_bal_model_fac)

beta_raw_orig_cancer_bal_model_reg <- rfPredictReg(beta_raw_orig_cancer_bal_bh,
                                                   cutoff = .7,
                                                   iterations = 10)

plotModel(beta_raw_orig_cancer_bal_model_reg, main = 'raw prepocessing, orig, bal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########  
# random
##########
beta_raw_orig_rand_fac <- rfPredictFac(beta_raw_orig_rand, 
                                       cutoff = .7, 
                                       iterations = 10)

conMatrix(beta_raw_orig_rand_fac)

beta_raw_orig_rand_reg <- rfPredictReg(beta_raw_orig_rand,
                                       cutoff = .7,
                                       iterations = 10)

plotModel(beta_raw_orig_rand_reg, main = 'raw random, orig',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

######################################################################################################

##########
# BETA swan - original
##########

##########
# global unbal
##########
beta_swan_orig_global_unbal_model_fac <- rfPredictFac(beta_swan_orig_global_unbal_bh, 
                                                      cutoff = .7, 
                                                      iterations = 10)

conMatrix(beta_swan_orig_global_unbal_model_fac)

beta_swan_orig_global_unbal_model_reg <- rfPredictReg(beta_swan_orig_global_unbal_bh,
                                                      cutoff = .7,
                                                      iterations = 10)

plotModel(beta_swan_orig_global_unbal_model_reg, main = 'swan prepocessing, orig, unbal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# global bal
##########

beta_swan_orig_global_bal_model_fac <- rfPredictFac(beta_swan_orig_global_bal_bh, 
                                                    cutoff = .7, 
                                                    iterations = 10)

conMatrix(beta_swan_orig_global_bal_model_fac)

beta_swan_orig_global_bal_model_reg <- rfPredictReg(beta_swan_orig_global_bal_bh,
                                                    cutoff = .7,
                                                    iterations = 10)

plotModel(beta_swan_orig_global_bal_model_reg, main = 'swan prepocessing, orig, bal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# cancer unbal
##########
beta_swan_orig_cancer_unbal_model_fac <- rfPredictFac(beta_swan_orig_cancer_unbal_bh, 
                                                      cutoff = .7, 
                                                      iterations = 10)

conMatrix(beta_swan_orig_cancer_unbal_model_fac)

beta_swan_orig_cancer_unbal_model_reg <- rfPredictReg(beta_swan_orig_cancer_unbal_bh,
                                                      cutoff = .7,
                                                      iterations = 10)

plotModel(beta_swan_orig_cancer_unbal_model_reg, main = 'swan prepocessing, orig, unbal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))


##########
# cancer bal
##########
beta_swan_orig_cancer_bal_model_fac <- rfPredictFac(beta_swan_orig_cancer_bal_bh, 
                                                    cutoff = .7, 
                                                    iterations = 10)

conMatrix(beta_swan_orig_cancer_bal_model_fac)

beta_swan_orig_cancer_bal_model_reg <- rfPredictReg(beta_swan_orig_cancer_bal_bh,
                                                    cutoff = .7,
                                                    iterations = 10)

plotModel(beta_swan_orig_cancer_bal_model_reg, main = 'swan prepocessing, orig, bal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# random
##########
beta_swan_orig_rand_fac <- rfPredictFac(beta_swan_orig_rand, 
                                        cutoff = .7, 
                                        iterations = 10)

conMatrix(beta_swan_orig_rand_fac)

beta_swan_orig_rand_reg <- rfPredictReg(beta_swan_orig_rand,
                                        cutoff = .7,
                                        iterations = 10)

plotModel(beta_swan_orig_rand_reg, main = 'swan random, orig',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

######################################################################################################

##########
# BETA quan - origsformed
##########

##########
# global unbal
##########
beta_quan_orig_global_unbal_model_fac <- rfPredictFac(beta_quan_orig_global_unbal_bh, 
                                                      cutoff = .7, 
                                                      iterations = 10)

conMatrix(beta_quan_orig_global_unbal_model_fac)

beta_quan_orig_global_unbal_model_reg <- rfPredictReg(beta_quan_orig_global_unbal_bh,
                                                      cutoff = .7,
                                                      iterations = 10)

plotModel(beta_quan_orig_global_unbal_model_reg, main = 'quan prepocessing, orig, unbal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# global bal
##########

beta_quan_orig_global_bal_model_fac <- rfPredictFac(beta_quan_orig_global_bal_bh, 
                                                    cutoff = .7, 
                                                    iterations = 10)

conMatrix(beta_quan_orig_global_bal_model_fac)

beta_quan_orig_global_bal_model_reg <- rfPredictReg(beta_quan_orig_global_bal_bh,
                                                    cutoff = .7,
                                                    iterations = 10)

plotModel(beta_quan_orig_global_bal_model_reg, main = 'quan prepocessing, orig, bal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# cancer unbal
##########
beta_quan_orig_cancer_unbal_model_fac <- rfPredictFac(beta_quan_orig_cancer_unbal_bh, 
                                                      cutoff = .7, 
                                                      iterations = 10)

conMatrix(beta_quan_orig_cancer_unbal_model_fac)

beta_quan_orig_cancer_unbal_model_reg <- rfPredictReg(beta_quan_orig_cancer_unbal_bh,
                                                      cutoff = .7,
                                                      iterations = 10)

plotModel(beta_quan_orig_cancer_unbal_model_reg, main = 'quan prepocessing, orig, unbal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))


##########
# cancer bal
##########
beta_quan_orig_cancer_bal_model_fac <- rfPredictFac(beta_quan_orig_cancer_bal_bh, 
                                                    cutoff = .7, 
                                                    iterations = 10)

conMatrix(beta_quan_orig_cancer_bal_model_fac)

beta_quan_orig_cancer_bal_model_reg <- rfPredictReg(beta_quan_orig_cancer_bal_bh,
                                                    cutoff = .7,
                                                    iterations = 10)

plotModel(beta_quan_orig_cancer_bal_model_reg, main = 'quan prepocessing, orig, bal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# random
##########
beta_quan_orig_rand_fac <- rfPredictFac(beta_quan_orig_rand, 
                                        cutoff = .7, 
                                        iterations = 10)

conMatrix(beta_quan_orig_rand_fac)

beta_quan_orig_rand_reg <- rfPredictReg(beta_quan_orig_rand,
                                        cutoff = .7,
                                        iterations = 10)

plotModel(beta_quan_orig_rand_reg, main = 'quan random, orig',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

######################################################################################################

##########
# BETA funnorm - origsformed
##########

##########
# global unbal
##########
beta_funnorm_orig_global_unbal_model_fac <- rfPredictFac(beta_funnorm_orig_global_unbal_bh, 
                                                         cutoff = .7, 
                                                         iterations = 10)

conMatrix(beta_funnorm_orig_global_unbal_model_fac)

beta_funnorm_orig_global_unbal_model_reg <- rfPredictReg(beta_funnorm_orig_global_unbal_bh,
                                                         cutoff = .7,
                                                         iterations = 10)

plotModel(beta_funnorm_orig_global_unbal_model_reg, main = 'funnorm prepocessing, orig, unbal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# global bal
##########

beta_funnorm_orig_global_bal_model_fac <- rfPredictFac(beta_funnorm_orig_global_bal_bh, 
                                                       cutoff = .7, 
                                                       iterations = 10)

conMatrix(beta_funnorm_orig_global_bal_model_fac)

beta_funnorm_orig_global_bal_model_reg <- rfPredictReg(beta_funnorm_orig_global_bal_bh,
                                                       cutoff = .7,
                                                       iterations = 10)

plotModel(beta_funnorm_orig_global_bal_model_reg, main = 'funnorm prepocessing, orig, bal, global',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# cancer unbal
##########
beta_funnorm_orig_cancer_unbal_model_fac <- rfPredictFac(beta_funnorm_orig_cancer_unbal_bh, 
                                                         cutoff = .7, 
                                                         iterations = 10)

conMatrix(beta_funnorm_orig_cancer_unbal_model_fac)

beta_funnorm_orig_cancer_unbal_model_reg <- rfPredictReg(beta_funnorm_orig_cancer_unbal_bh,
                                                         cutoff = .7,
                                                         iterations = 10)

plotModel(beta_funnorm_orig_cancer_unbal_model_reg, main = 'funnorm prepocessing, orig, unbal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))


##########
# cancer bal
##########
beta_funnorm_orig_cancer_bal_model_fac <- rfPredictFac(beta_funnorm_orig_cancer_bal_bh, 
                                                       cutoff = .7, 
                                                       iterations = 10)

conMatrix(beta_funnorm_orig_cancer_bal_model_fac)

beta_funnorm_orig_cancer_bal_model_reg <- rfPredictReg(beta_funnorm_orig_cancer_bal_bh,
                                                       cutoff = .7,
                                                       iterations = 10)

plotModel(beta_funnorm_orig_cancer_bal_model_reg, main = 'funnorm prepocessing, orig, bal, cancer',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

##########
# random
##########
beta_funnorm_orig_rand_fac <- rfPredictFac(beta_funnorm_orig_rand, 
                                           cutoff = .7, 
                                           iterations = 10)

conMatrix(beta_funnorm_orig_rand_fac)

beta_funnorm_orig_rand_reg <- rfPredictReg(beta_funnorm_orig_rand,
                                           cutoff = .7,
                                           iterations = 10)

plotModel(beta_funnorm_orig_rand_reg, main = 'funnorm random, orig',
          xlim = c(0, 1000),
          ylim = c(0, 1000))

######################################################################################################
# ggplot for beta_funnorm cancer bal, transformed and raw
# 4 and 6 are predictions and real age 
pred_tran <- unlist(beta_funnorm_tran_cancer_bal_model_reg[[4]])
real_tran <- unlist(beta_funnorm_tran_cancer_bal_model_reg[[6]])
plotDataTran <- as.data.frame(cbind(predictions = pred_tran, real = real_tran))

pred_orig <- unlist(beta_funnorm_orig_cancer_bal_model_reg[[4]])
real_orig  <- unlist(beta_funnorm_orig_cancer_bal_model_reg[[6]])
plotDataOrig <- as.data.frame(cbind(predictions = pred_orig, real = real_orig))

ggplot(data = plotDataOrig, aes(pred_tran, real_tran)) + 
  geom_point(colour = 'black', alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.3) +
  xlim(0, 800) +
  ylim(0, 800) +
  xlab('Model predictions') +
  ylab('Real age of sample collection') +
  ggtitle('Age of sample collection') + theme(panel.background=element_rect(fill="#F0F0F0"), 
                                  plot.background=element_rect(fill="#F0F0F0"), 
                                  # panel.border=element_rect(colour="#F0F0F0"),
                                  panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
                                  legend.position="right",  plot.title=element_text(face="bold",hjust=-.08,vjust=2,colour="#3C3C3C",size=20),
                                  axis.text.x=element_text(size=11,colour="#535353",face="bold"),
                                  axis.text.y=element_text(size=11,colour="#535353",face="bold"),
                                  axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
                                  axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5))

###############################
# MAYBE INSTEAD OF PLOTTING EACH ONE, CREAT A TABLE SIMILAR TO RESULTS TABLE FROM IDAT MODELS
