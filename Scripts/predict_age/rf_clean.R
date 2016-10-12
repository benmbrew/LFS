#############################################################################################################
# This script will be a clean version of the random forest function predicting age of diagnosis (log and not log) from 
# methylation, using full_data, and the two correlation data sets. It will also be able to run 
# the predictions with the residual from regressing each gene on age of sample collection for the three data 
# sets

##################################################################################################
# this script will read in different versions of subsetted full data and run random forest. 
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
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')
model_data <- paste0(data_folder, '/model_data')


###########################################
# Read in data- gene_knn, gene_lsa, 
# probe_knn, probe_lsa, and bh_features
###########################################

load(paste0(model_data, '/model_data.RData'))
load(paste0(model_data, '/bh_features.RData'))
rm(cg_locations)

# Data types: 
# 1) gene_knn, gene_lsa, probe_knn, probe_lsa all methylation data that needs
# to be joined with clinical data
# 2) probe_knn_cancer_features, probe_lsa_cancer_features,
#    probe_knn_global_features, probe_lsa_global_features
# simpy have the probes chosen as significant in the bumphunter analysis. 

###################
# function that combines model_data with 
# clinical data
###################
joinData <- function(model_data) {
  
  features <- colnames(model_data)[2:(length(colnames(model_data)))]
  model_data <- inner_join(model_data, clin, by = 'id')
  model_data <- model_data[!is.na(model_data$p53_germline),]
  model_data <- model_data[!duplicated(model_data$id),]
  model_data <- model_data[!duplicated(model_data$tm_donor_),]
  model_data <- model_data[, c('p53_germline', 'age_diagnosis', 'age_sample_collection', 'cancer_diagnosis_diagnoses', features)]
  return(model_data)
}

# take p53 germline column and relevel factors to get rid of NA level
relevelFactor <- function(model_data) {
  model_data$p53_germline <- factor(model_data$p53_germline, levels = c('Mut', 'WT'))
  return(model_data)
}

###################
# Function that takes model model_data and obstains residual features based
# on regressin each gene/probe on age of sample collection
###################
getResidual <- function(model_data) {
  
  feature_names <- colnames(model_data)[2:ncol(model_data)]
  
  resid <- list()
  
  for (i in 2:ncol(model_data)){
    
    resid[[i]] <- lm(model_data[, i] ~ model_data$age_sample_collection, data = model_data)$residuals
    
    print(i)
    
  }
  
  resid_data <- as.data.frame(do.call('cbind', resid))
  model_data <- cbind(model_data$age_diagnosis, model_data$age_sample_collection, resid_data)
  colnames(model_data) <- c('age_diagnosis', 'age_sample_collection', feature_names)
  
}


#########################################
# function that takes data and arguments for regression, cutoff, and selected features. 
#########################################
# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       fac,
                       reg,
                       log,
                       resid,
                       cutoff,
                       random,
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
  train.sample_collection <- list()
  test.sample_collection <- list()
  test_acc <- list()
  test_stats  <- list()
  test_acc_samp <- list()
  test_stats_samp <- list()
  
  # get features, and subset by complete age of diagnosis 
  selected_features <- colnames(data)[3:(ncol(data) -2)]
  
  if (random){
    
    selected_features <- sample(selected_features, 100, replace = F)
    
  }
  
  if (fac) {
    
    # remove first two columns
    data <- data[, c('age_diagnosis_fac', 'age_sample_fac', selected_features)]

  }
  
  # set log transformation
  if (reg) {
    
    data <- data[, c('age_diagnosis', 'age_sample_collection', selected_features)]
    
    if (log & !resid) {
      
      data <- log(data)
      data <- as.data.frame(data)
      
    }
  
    if (log & resid) {
    
      data[, 3:ncol(data)] <-  data[, 3:ncol(data)] + abs(min(data[, 3:ncol(data)])) + .1
      data <- log(data)
    
    }
  
  }

  obs <- nrow(data)
  
  for (i in 1:iterations) {
    
    set.seed(i)
    
    train_index <- sample(nrow(data), nrow(data) *cutoff, replace = F)
    
    if (fac) {
      
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
      
      y = make.names(as.factor(data$age_diagnosis_fac[train_index]))
      
      
    } else {
        
      NFOLDS <- 2
        fitControl <- trainControl( 
          method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
          number = min(10, NFOLDS),      
          repeats = 1,
          allowParallel = TRUE
          )
      
        y <- data$age_diagnosis[train_index]
        
    }
    
    # mtry: Number of variables randomly sampled as candidates at each split.
    # ntree: Number of trees to grow.
  
    
      mtry <- sqrt(ncol(data[train_index, selected_features]))
      tunegrid <- expand.grid(.mtry=mtry)
      
      model[[i]] <- train(x = data[train_index, selected_features]
                          , y = y
                          , method = "rf"
                          , trControl = fitControl
                          , tuneGrid = tunegrid
                          , importance = T
                          , verbose = FALSE)
      
      temp <- varImp(model[[i]])[[1]]
      importance[[i]] <- cbind(rownames(temp), temp$Overall)
      
      if (fac) {
        test.predictions[[i]] <- predict(model[[i]] 
                                         , newdata = data[-train_index, selected_features]
                                         , type = "prob")
        
        train.predictions[[i]] <- predict(model[[i]] 
                                          , newdata = data[train_index, selected_features]
                                          ,type = "prob")
        
      } else {
        test.predictions[[i]] <- predict(model[[i]] 
                                         , newdata = data[-train_index, selected_features])
        
        train.predictions[[i]] <- predict(model[[i]] 
                                          , newdata = data[train_index, selected_features])
        
      }
      
      
      if (fac) {
        
        train.ground_truth[[i]] <- as.factor(make.names(data$age_diagnosis_fac[train_index]))
        test.ground_truth[[i]] <- as.factor(make.names(data$age_diagnosis_fac[-train_index]))
        train.sample_collection[[i]] = as.factor(make.names(data$age_sample_fac[train_index]))
        train.sample_collection[[i]] <- factor(train.sample_collection[[i]], levels = c("X1", "X2"))
        # temp <- as.factor(make.names(data$age_sample_fac[-train_index]))
        test.sample_collection[[i]] = as.factor(make.names(data$age_sample_fac[-train_index]))
        test.sample_collection[[i]] <- factor(test.sample_collection[[i]], levels = c("X1", "X2"))
        
      
        # For age of diagnosis
        # Accuracy
        test_acc[[i]] <- sum(levels(test.ground_truth[[i]])[apply(test.predictions[[i]], 1, which.is.max)] == test.ground_truth[[i]]) / dim(test.predictions[[i]])[1]
        # Confustion Matrix
        test_stats[[i]] <- confusionMatrix(levels(test.ground_truth[[i]])[apply(test.predictions[[i]], 1, which.is.max)], test.ground_truth[[i]])
        # print(rf.test_stats)
        missing_ind <- !is.na(unlist(test.sample_collection[[i]]))
        test.sample_collection[[i]] <- unlist(test.sample_collection[[i]])[missing_ind]
        test.predictions[[i]] <- test.predictions[[i]][missing_ind,]
        # for age of collection
        # subset test.predictions by missing index in age.sample_collection, and remove the NAs in age.sample collection
        
        # Accuracy
        test_acc_samp[[i]] <- sum(levels(test.sample_collection[[i]])[apply(test.predictions[[i]], 1, which.is.max)] == test.sample_collection[[i]], na.rm = T) / dim(test.predictions[[i]])[1]
        # Confustion Matrix
        # subset to remove NAs
        test_stats_samp[[i]] <- confusionMatrix(levels(test.sample_collection[[i]])[apply(test.predictions[[i]], 1, which.is.max)], test.sample_collection[[i]])
        # print(rf.test_stats)
      } else {
        
        train.ground_truth[[i]] <- data$age_diagnosis[train_index]
        test.ground_truth[[i]] <- data$age_diagnosis[-train_index]
        train.sample_collection[[i]] = data$age_sample_collection[train_index]
        test.sample_collection[[i]] = data$age_sample_collection[-train_index]
        train.mse[[i]] <- rmse(unlist(train.predictions[[i]]), unlist(train.ground_truth[[i]]))
        test.mse[[i]] <- rmse(unlist(test.predictions[[i]]), unlist(test.ground_truth[[i]]))
        
      }
      
    
    print(i)
    
  }
  
  return(list(train.mse, test.mse,  train.predictions, test.predictions, train.ground_truth, test.ground_truth, train.sample_collection,
              test.sample_collection, test_acc, test_stats, test_acc_samp, test_stats_samp, model, importance, obs))
  
}

# 1 - train.mse
# 2 - test.mse
# 3 - train.predictions
# 4 - test.predictions
# 5 - train.ground_truth
# 6 - test.ground_truth
# 7 - train.sample_collection
# 8 - test.sample_collection
# 9 - test_acc
# 10 - test_stats
# 11 - test_acc_samp
# 12 - test_stats_samp
# 13 - model
# 14 - importance
# 15 - obs

# plot function 
plotModel <- function(result_list,
                      main1,
                      main2,
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
       main = main1)
  abline(0,1, lty = 3)
  corr <- round(cor(unlist(result_list[[4]]), unlist(result_list[[6]])), 2)
  legend("topleft", legend = paste0('correlation = ', corr), cex = 1, bty = 'n')
  legend("bottomright", legend = paste0('# obs = ', result_list[[15]]), cex = 1, bty = 'n')
  
  
  # plot predictions against ground truth
  plot(unlist(result_list[[4]]), unlist(result_list[[8]]), 
       xlim = xlim,
       ylim = ylim,
       bty = 'n',
       col = adjustcolor('blue', alpha.f = 0.6),
       pch = 16,
       xlab = 'Predictions',
       ylab = 'Real Age of Sample Collection',
       main = main2)
  abline(0,1, lty = 3)
  corr <- round(cor(unlist(result_list[[4]]), unlist(result_list[[8]]), use = "complete.obs"), 2)
  legend("topleft", legend = paste0('correlation = ', corr), cex = 1, bty = 'n')
   legend("bottomright", legend = paste0('# obs = ', result_list[[15]]), cex = 1, bty = 'n')
  
}

#######################################################################################
# Methylation
#######################################################################################

###################
# Regression - all data
###################

####
# regression, not log, not residual
####
methyl_reg <- predictAll(data = full_data_probe_mut,
                         random = F,
                         reg = T,
                         fac = F,
                         log = F,
                         cutoff = .7,
                         resid = F,
                         iterations = 10)

####
# regression, log, not resid
####
methyl_reg_log <- predictAll(data = full_data_probe,
                             random = F,
                             reg = T,
                             fac = F,
                             log = T,
                             cutoff = .7,
                             resid = F,
                             iterations = 10)
####
# regression, not log, resid
####
methyl_reg_resid <- predictAll(data = resid_full_probe,
                               random = F,
                               reg = T,
                               fac = F,
                               log = F,
                               cutoff = .7,
                               resid = T,
                               iterations = 10)

####
# regression, log, resid
####
methyl_reg_log_resid <- predictAll(data = resid_full_probe,
                                   random = F,
                                   reg = T,
                                   fac = F,
                                   log = T,
                                   cutoff = .7,
                                   resid = T,
                                   iterations = 10)


plotModel(methyl_reg,
          'Predicting Age of Onset',
          'Predicting Age of Sample Collection',
          xlim = c(0,1000),
          ylim = c(0,1000))

plotModel(methyl_reg_log,
          'Predicting Age of Onset (Log)',
          'Predicting Age of Sample Collection (Log)',
          xlim = c(0,10),
          ylim = c(0,10))

plotModel(methyl_reg_resid,
          'Predicting Age of Onset Residuals',
          'Predicting Age of Sample Collection Residuals',
          xlim = c(0,1000),
          ylim = c(0,1000))

plotModel(methyl_reg_log_resid,
          'Predicting Age of Onset Residuals (Log)',
          'Predicting Age of Sample Collection Residuals (Log)',
          xlim = c(0,15),
          ylim = c(0,15))


###################
# classification
###################
load('/home/benbrew/Desktop/confusion.matrix.RData')

# age of diagnosis, classification, not log
methyl_fac <- predictAll(data = full_data_probe,
                         random = F,
                         reg = F,
                         fac = T,
                         log = F,
                         cutoff = .7,
                         resid = F,
                         iterations = 10)

# test acc for age of diagnosis
mean(unlist(methyl_fac[[9]]))

# test acc for age of sample collection
mean(unlist(methyl_fac[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- methyl_fac[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations


# confustion matrix age of sample
temp <- list()
for (i in 1:10){
  temp[[i]] <- methyl_fac[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations

