##################################################################################################
# this script will subset data based on biggest difference between age of diagnosis and age of sample collection and predict
# on thresholds

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

# #split the data by how close together age of diagnosis and age of sample collection are
# diff <- full_data_rf$age_sample_collection - full_data_rf$age_diagnosis
# histogram(diff)
# median(diff, na.rm = T)
# mean(diff, na.rm = T)
# diff <- diff[!is.na(diff)]
# 
# length(diff)
# length(diff[diff > 20])

# # with log
# data <- full_data_rf
# data[,c(6,8, 27:ncol(data))]  <- log(data[,c(6,8, 27:ncol(data))])
# diff <- data$age_sample_collection - data$age_diagnosis
# histogram(diff)
# median(diff, na.rm = T)
# mean(diff, na.rm = T)
# diff <- diff[!is.na(diff)]
# 
# length(diff)
# length(diff[diff > 20])

# use 20 months as cutoff

# Create binary variables for age of diagnosis and age of sample collection on both data sets
full_data$age_diagnosis_fac <- ifelse(full_data$age_diagnosis <= 48, 1, 2)

full_data$age_sample_fac<- ifelse(full_data$age_sample_collection <= 48, 1, 2)

full_data_rf$age_diagnosis_fac <- ifelse(full_data_rf$age_diagnosis <= 48, 1, 2)

full_data_rf$age_sample_fac <- ifelse(full_data_rf$age_sample_collection <= 48, 1, 2)

# Create binomial variables for age of diagnosis and age of sample collection
full_data$age_diagnosis_multi <- ifelse(full_data$age_diagnosis <= 30, 1, 
                                        ifelse(full_data$age_diagnosis > 30 & full_data$age_diagnosis <= 50, 2,
                                               ifelse(full_data$age_diagnosis > 50 & full_data$age_diagnosis <= 300, 3, 4)))

full_data$age_sample_multi <- ifelse(full_data$age_sample_collection <= 30, 1, 
                                     ifelse(full_data$age_sample_collection > 30 & full_data$age_sample_collection <= 50, 2,
                                            ifelse(full_data$age_sample_collection > 50 & full_data$age_sample_collection <= 300, 3, 4)))

full_data_rf$age_diagnosis_multi <- ifelse(full_data_rf$age_diagnosis <= 30, 1, 
                                           ifelse(full_data_rf$age_diagnosis > 30 & full_data_rf$age_diagnosis <= 50, 2,
                                                  ifelse(full_data_rf$age_diagnosis > 50 & full_data_rf$age_diagnosis <= 300, 3, 4)))

full_data_rf$age_sample_multi <- ifelse(full_data_rf$age_sample_collection <= 30, 1, 
                                        ifelse(full_data_rf$age_sample_collection > 30 & full_data_rf$age_sample_collection <= 50, 2,
                                               ifelse(full_data_rf$age_sample_collection > 50 & full_data_rf$age_sample_collection <= 300, 3, 4)))

summary(as.factor(full_data_rf$age_diagnosis_multi))
summary(as.factor(full_data_rf$age_sample_multi))

# Create multinomal variables for age of diagnosis and age of sample collection
full_data$age_diagnosis_year <- ifelse(full_data$age_diagnosis <= 24, 1, 
                                       ifelse(full_data$age_diagnosis > 24 & full_data$age_diagnosis <= 48, 2,
                                              ifelse(full_data$age_diagnosis > 48 & full_data$age_diagnosis <= 180, 3,
                                                     ifelse(full_data$age_diagnosis > 180 & full_data$age_diagnosis <= 360, 4, 5))))

full_data$age_sample_year <- ifelse(full_data$age_diagnosis <= 24, 1, 
                                    ifelse(full_data$age_diagnosis > 24 & full_data$age_diagnosis <= 48, 2,
                                           ifelse(full_data$age_diagnosis > 48 & full_data$age_diagnosis <= 180, 3,
                                                  ifelse(full_data$age_diagnosis > 180 & full_data$age_diagnosis <= 360, 4, 5))))

full_data_rf$age_diagnosis_year <- ifelse(full_data_rf$age_diagnosis <= 24, 1, 
                                          ifelse(full_data_rf$age_diagnosis > 24 & full_data_rf$age_diagnosis <= 48, 2,
                                                 ifelse(full_data_rf$age_diagnosis > 48 & full_data_rf$age_diagnosis <= 180, 3,
                                                        ifelse(full_data_rf$age_diagnosis > 180 & full_data_rf$age_diagnosis <= 360, 4, 5))))

full_data_rf$age_sample_year <- ifelse(full_data_rf$age_diagnosis <= 24, 1, 
                                       ifelse(full_data_rf$age_diagnosis > 24 & full_data_rf$age_diagnosis <= 48, 2,
                                              ifelse(full_data_rf$age_diagnosis > 48 & full_data_rf$age_diagnosis <= 180, 3,
                                                     ifelse(full_data_rf$age_diagnosis > 180 & full_data_rf$age_diagnosis <= 360, 4, 5))))

summary(as.factor(full_data_rf$age_diagnosis_year))
summary(as.factor(full_data_rf$age_sample_year))

# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       subset, 
                       selected_features,
                       binary,
                       multi,
                       log,
                       cutoff,
                       iterations) {
  
  if(log) {
    
    data[,27:ncol(data)]  <- log(data[, 27:ncol(data)])
  }
  
  genes <- colnames(data)[27:(ncol(data)-6)]
  
  data <- data[, c(subset, genes)]
  
  # Try the model with all different selection of features based on number of missinginess. 
  data <- data[complete.cases(data),]
  
  # set train index to difference less than 20
  
  train_index <- data[,2]  == data[,1]
  
  obs <- nrow(data)

  if (binary){
    level_count <- 2
  } else if (multi) {
    level_count <- 4
  } else {
    level_count <- 5
  }
  
  rf_y = as.factor(make.names(data[,1][train_index]))
  
  # 4) Random Forest 
  
  if (length(levels(rf_y)) == 2) {
    summaryFunc <- twoClassSummary
  } else {
    summaryFunc <- multiClassSummary
  }
  
  # determines how you train the model.
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = 2, 
    classProbs = TRUE,     
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = summaryFunc)
  
  mtry <- sqrt(ncol(data))
  tunegrid <- expand.grid(.mtry=mtry)
  
  model <- train(x = data[train_index, c(selected_features, genes)]
                      , y = rf_y
                      , method = "rf"
                      , trControl = fitControl                   
                      , verbose = FALSE
                      , metric = "logLoss")
  
  predictions <- predict(model 
                              , newdata = data[!train_index, c(selected_features, genes)]
                              , type = "prob")
  
  test.ground_truth <- as.factor(data[, 1][!train_index])
  test.ground_truth_samp <- as.factor(data[,2][!train_index])
  
  # test.ground_truth_1inN <- as.factor(class.ind(data[, 2][!train_index]))
  #print(table(levels(test.ground_truth)[apply(rf.predictions, 1, which.is.max)], test.ground_truth))
  # AUC
  # #create ROCR prediction object
  # temp.predict <- prediction(unlist(predictions[[i]]), test.ground_truth_1inN)  
  # 
  # if(auc){
  #   rf.test_auc <- unlist(slot(performance(temp.predict, "auc"), "y.values"))  
  #   print(paste("RF test AUC:", rf.test_auc))  
  # }
  
  # Accuracy
  rf.test_acc <- sum(levels(test.ground_truth)[apply(predictions, 1, which.is.max)] == test.ground_truth) / dim(predictions)[1]
  # print(paste("RF test acc:", rf.test_acc))  
  # Compute Confusion Matrix and Statistics
  #confusionMatrix(pred, truth)
  rf.test_stats <- confusionMatrix(levels(test.ground_truth)[apply(predictions, 1, which.is.max)], test.ground_truth)
  rf.test_stats_samp <- confusionMatrix(levels(test.ground_truth)[apply(predictions, 1, which.is.max)], test.ground_truth)
  
  # print(rf.test_stats)
  

  return(list(predictions, test.ground_truth, rf.test_acc, model, rf.test_stats, obs))
  
}

# just methylation
rf_methyl <- predictAll(data = full_data_rf,
                        subset <- c("age_diagnosis_fac", "age_sample_fac"),
                        selected_features = NULL,
                        binary = T,
                        log = F,
                        cutoff = 10)

# just methylation
rf_methyl_multi <- predictAll(data = full_data_rf,
                        subset <- c("age_diagnosis_multi", "age_sample_multi"),
                        selected_features = NULL,
                        binary = F,
                        multi = T,
                        log = F,
                        cutoff = 10)

# just methylation
rf_methyl_year <- predictAll(data = full_data_rf,
                              subset <- c("age_diagnosis", "age_sample_collection"),
                              selected_features = NULL,
                              binary = F,
                              multi = F,
                              log = F,
                              cutoff = 10)



###############################################################################################
# just methylation log
rf_methyl_log <- predictAll(data = full_data_rf,
                            subset <- c("age_diagnosis", "age_sample_collection"),
                            selected_features = NULL, 
                            binary = T,
                            log = T,
                            cutoff = 0.15)


