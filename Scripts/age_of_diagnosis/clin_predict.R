###########################################################################
# this script will predict clinical variables using methylation data.

library(dplyr)
library(stringr)
library(impute)
library(mlbench)
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(dplyr)
library(Metrics)
library(doParallel)
library(ROCR)
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

# class.ind function
class.ind <- function(cl)
{
  n <- length(cl)
  cl <- as.factor(cl)
  x <- matrix(0, n, length(levels(cl)) )
  x[(1:n) + n*(unclass(cl)-1)] <- 1
  dimnames(x) <- list(names(cl), levels(cl))
  x
}


# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       variable, 
                       selected_features,
                       iterations) {
  
  model <- list()
  predictions <- list()
  mse <- list()
  importance <- list()
  test.ground_truth <- list()

  genes <- colnames(data)[28:ncol(data)]
  
  data <- data[, c(variable, genes)]
  data[, variable] <- as.factor(data[, variable])
  
  # Try the model with all different selection of features based on number of missinginess. 
  data <- data[complete.cases(data),]
  
  obs <- nrow(data)
  
  for (i in 1:iterations){
    
    set.seed(i)
    train_index <- sample(nrow(data), nrow(data) *.7, replace = F)
    
    # 4) Random Forest 
    rf_y = make.names(as.factor(data[, variable])[train_index])

    if (length(levels(rf_y)) == 2) {
      summaryFunc <- twoClassSummary
    } else {
      summaryFunc <- multiClassSummary
    }
    
    # determines how you train the model.
    fitControl <- trainControl( 
      method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
      number = 5, 
      classProbs = TRUE,     
      repeats = 1,
      allowParallel = TRUE,
      summaryFunction = summaryFunc)
    
    rf.model <- train(x = data[train_index, c(selected_features, genes)]
                      , y = rf_y
                      , method = "rf"
                      , trControl = fitControl                   
                      , verbose = FALSE
                      , metric = "ROC")
    
    rf.predictions <- predict(rf.model 
                              , newdata = data[-train_index, c(selected_features, genes)]
                              , type = "prob")
    
    test.ground_truth <- data[, variable][-train_index]
    test.ground_truth_1inN <- class.ind(data[, variable][-train_index])
    #print(table(levels(test.ground_truth)[apply(rf.predictions, 1, which.is.max)], test.ground_truth))
    # AUC
    #create ROCR prediction object
    temp.predict <- prediction(rf.predictions, test.ground_truth_1inN)  
    rf.test_auc <- unlist(slot(performance(temp.predict, "auc"), "y.values"))  
    print(paste("RF test AUC:", rf.test_auc))  
    # Accuracy
    rf.test_acc[[i]] <- sum(levels(test.ground_truth)[apply(rf.predictions, 1, which.is.max)] == test.ground_truth) / dim(rf.predictions)[1]
    print(paste("RF test acc:", rf.test_acc))  
    # Compute Confusion Matrix and Statistics
    #confusionMatrix(pred, truth)
    rf.test_stats[[i]] <- confusionMatrix(levels(test.ground_truth)[apply(rf.predictions, 1, which.is.max)], test.ground_truth)
    print(rf.test_stats)
    
    print(i)
    
  }
  
  return(list(rf.test_stats, predictions, model, rf.test_acc, test.ground_truth, obs))
  
}

# predict mdm2.nG
rf_mdm2.nG <- predictAll(data = full_data_rf,
                        variable = "mdm2.nG",
                        selected_features = NULL, iterations = 50)









#
#
# plot(unlist(rf_methyl[[2]]), unlist(rf_methyl[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = 'Just methylation')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_methyl[[2]]) ~ unlist(rf_methyl[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_methyl[[6]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))
#
#
#
# plot(unlist(rf_methyl[[2]]), unlist(rf_methyl[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = 'Just methylation')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_methyl[[6]]))


