#############################################################################################################
# This script will be a clean version of the random forest function predicting age of diagnosis (log and not log) from 
# methylation and trained on observaitons that have very little difference between age of diagnosis and age of sample
# collection and test on observations that a large difference between age of diagnosis and age of sample collection.

##################################################################################################
# Load libraries
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
test <- paste0(project_folder, '/Scripts/Analyze')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

# Read in 3 different data sets 
full_data <- read.csv(paste0(data_folder, '/full_data.csv'), stringsAsFactors = F)
full_data_cor <- read.csv(paste0(data_folder, '/full_data_cor.csv'), stringsAsFactors = F)
full_data_rf <- read.csv(paste0(data_folder, '/full_data_rf.csv'), stringsAsFactors = F)

full_data$X <- NULL
full_data_cor$X <- NULL
full_data_rf$X <- NULL


# make categroica variable from age of methylaion and age of sample collection
full_data_rf$age_diagnosis_fac <- as.integer(ifelse(full_data_rf$age_diagnosis <= 48, 1, 2))

full_data_rf$age_sample_fac <- as.integer(ifelse(full_data_rf$age_sample_collection <= 48, 1, 2))

# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       fac,
                       subset, 
                       selected_features,
                       cutoff,
                       log,
                       iterations) {
  
  model <- list()
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
  
  
  # set log transformation
  if(log) {
    
    data[,c(6,8, 27:(ncol(data) - 2))]  <- log(data[,c(6,8,27:(ncol(data) -2))])
  }
  
  
  
  genes <- colnames(data)[27:(ncol(data) - 2)]
  data <- data[, c(subset, genes)]
  data <- data[complete.cases(data),]
  
  
  
  
  for ( i in 3:ncol(data)){
    
    if(typeof(data[,i]) == 'character' || typeof(data[,i]) == 'integer') {
      data[,i] <- as.factor(data[,i])
    }
  }
  
  
  
  obs <- nrow(data)
  
  
  for (i in 1:iterations){
    
    set.seed(i)
    train_index <- abs(data$age_sample_collection - data$age_diagnosis) < cutoff
    
    if(fac){
      rf_y = make.names(as.factor(data$age_diagnosis_fac[train_index]))
    }else {
      rf_y = data$age_diagnosis[train_index]
    }
    
    if(fac) {
      
      if (length(levels(rf_y)) == 2) {
        summaryFunc <- twoClassSummary
      } else {
        summaryFunc <- multiClassSummary
      }
      # determines how you train the model.
      NFOLDS <- 2
      fitControl <- trainControl( 
        method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
        number = min(10, NFOLDS),
        classProbs = TRUE,
        repeats = 1,
        allowParallel = TRUE,
        summaryFunction = summaryFunc
      )
    } else {
      NFOLDS <- 2
      fitControl <- trainControl( 
        method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
        number = min(10, NFOLDS),      
        repeats = 1,
        allowParallel = TRUE
      )
    }
    
    # mtry: Number of variables randomly sampled as candidates at each split.
    # ntree: Number of trees to grow.
    
    mtry <- sqrt(ncol(data))
    tunegrid <- expand.grid(.mtry=mtry)
    
    model[[i]] <- train(x = data[train_index, c(selected_features, genes)]
                        , y = rf_y
                        , method = "rf"
                        , trControl = fitControl
                        , tuneGrid = tunegrid
                        , importance = T
                        , verbose = FALSE)
    
    temp <- varImp(model[[i]])[[1]]
    importance[[i]] <- cbind(rownames(temp), temp$Overall)
    
    if(fac){
      test.predictions[[i]] <- predict(model[[i]] 
                                       , newdata = data[-train_index, c(selected_features, genes)]
                                       , type = "prob")
      
      train.predictions[[i]] <- predict(model[[i]] 
                                        , newdata = data[train_index, c(selected_features, genes)]
                                        ,type = "prob")
      
    } else {
      test.predictions[[i]] <- predict(model[[i]] 
                                       , newdata = data[-train_index, c(selected_features, genes)])
      
      train.predictions[[i]] <- predict(model[[i]] 
                                        , newdata = data[train_index, c(selected_features, genes)])
      
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

# age of diagnosis, regression, not log
methyl_reg <- predictAll(data = full_data_rf,
                         fac = F,
                         log = F,
                         subset = c('age_diagnosis', 'age_sample_collection'),
                         selected_features = NULL,
                         cutoff = 10,
                         iterations = 10)

# plot predictions against ground truth
plot(unlist(methyl_reg[[4]]), unlist(methyl_reg[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'Age in Months')
abline(0,1)
r_squared <- round(summary(lm(unlist(methyl_reg[[4]]) ~ unlist(methyl_reg[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', methyl_reg[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(methyl_reg[[4]]), unlist(methyl_reg[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'Age in Months')
abline(0,1)
r_squared <- round(summary(lm(unlist(methyl_reg[[4]]) ~ unlist(methyl_reg[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', methyl_reg[[15]]), cex = 0.7)



# age of diagnosis, regression with log transform
methyl_reg_log <- predictAll(data = full_data_rf,
                             fac = F,
                             log = T,
                             subset = c('age_diagnosis', 'age_sample_collection'),
                             selected_features = NULL,
                             cutoff = .7,
                             iterations = 10)

# plot predictions against ground truth
plot(unlist(methyl_reg_log[[4]]), unlist(methyl_reg_log[[6]]), 
     xlim = c(0, 8),
     ylim = c(0, 8),
     xlab = 'Log Predictions',
     ylab = 'Log Real Age of Diagnosis',
     main = 'Log Months')
abline(0,1)
r_squared <- round(summary(lm(unlist(methyl_reg_log[[4]]) ~ unlist(methyl_reg_log[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', methyl_reg_log[[15]]), cex = 0.7)


# plot predictions against age of sample collection
plot(unlist(methyl_reg_log[[4]]), unlist(methyl_reg_log[[8]]), 
     xlim = c(0, 8),
     ylim = c(0, 8),
     xlab = 'Log Predictions',
     ylab = 'Log Real Age of Sample Collection',
     main = 'Log Months')
abline(0,1)
r_squared <- round(summary(lm(unlist(methyl_reg_log[[4]]) ~ unlist(methyl_reg_log[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', methyl_reg_log[[15]]), cex = 0.7)


