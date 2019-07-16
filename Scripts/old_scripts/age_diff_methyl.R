##################################################################################################
# this script will subset data based on biggest difference between age of diagnosis and age of sample collection.
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

#split the data by how close together age of diagnosis and age of sample collection are
diff <- full_data_rf$age_sample_collection - full_data_rf$age_diagnosis
histogram(diff)
median(diff, na.rm = T)
mean(diff, na.rm = T)
diff <- diff[!is.na(diff)]

length(diff)
length(diff[diff > 20])

# with log
data <- full_data_rf
data[,c(6,8, 27:ncol(data))]  <- log(data[,c(6,8, 27:ncol(data))])
diff <- data$age_sample_collection - data$age_diagnosis
histogram(diff)
median(diff, na.rm = T)
mean(diff, na.rm = T)
diff <- diff[!is.na(diff)]

length(diff)
length(diff[diff > 20])

# use 20 months as cutoff

# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       subset,
                       selected_features,
                       log,
                       cutoff) {
  
  
  # set log transformation
  if(log) {

    data[,c(6,8, 27:ncol(data))]  <- log(data[,c(6,8,27:ncol(data))])
  }
  genes <- colnames(data)[27:ncol(data)]
  
  data <- data[, c(subset, genes)]
  
  # Try the model with all different selection of features based on number of missinginess. 
  data <- data[complete.cases(data),]
  
  # set train index to difference less than 20
  
  train_index <- abs(data$age_sample_collection - data$age_diagnosis) < cutoff
  
  obs <- nrow(data)

  # determines how you train the model.
  NFOLDS <- 2
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = min(10, NFOLDS),      
    repeats = 1,
    allowParallel = TRUE
    #summaryFunction = summaryFunc
  )
  
  # mtry: Number of variables randomly sampled as candidates at each split.
  # ntree: Number of trees to grow.

    mtry <- sqrt(ncol(data))
    tunegrid <- expand.grid(.mtry=mtry)
    
    rf_y = data$age_diagnosis[train_index]
    
    
    model <- train(x = data[train_index, c(selected_features, genes)]
                        , y = rf_y
                        , method = "rf"
                        , trControl = fitControl
                        , tuneGrid = tunegrid
                        , importance = T
                        , verbose = FALSE)
    
    temp <- varImp(model)[[1]]
    importance<- cbind(rownames(temp), temp$Overall)
    
    
    test.predictions <- predict(model
                                , newdata = data[!train_index, c(selected_features, genes)])
    
    train.predictions <- predict(model 
                                , newdata = data[train_index, c(selected_features, genes)])
    
    train.ground_truth <- data$age_diagnosis[train_index]
    test.ground_truth <- data$age_diagnosis[!train_index]
    train.sample_collection = data$age_sample_collection[train_index]
    test.sample_collection = data$age_sample_collection[!train_index]
    train.mse <- rmse(unlist(train.predictions), unlist(train.ground_truth))
    test.mse <- rmse(unlist(test.predictions), unlist(test.ground_truth))
  
  return(list(train.mse, test.mse,  train.predictions, test.predictions, train.ground_truth, 
              test.ground_truth, train.sample_collection, test.sample_collection, model, importance, obs))
  
}


# just methylation
rf_methyl <- predictAll(data = full_data_rf,
                        subset <- c("age_diagnosis", "age_sample_collection"),
                        selected_features = NULL, 
                        log = F,
                        cutoff = 10)


###############################
# plot test age diagnosis
plot(unlist(rf_methyl[[4]]), unlist(rf_methyl[[6]]), 
     xlab = 'Test Predictions', 
     ylab = 'Test Age of Diagnosis',
     xlim = c(0,800),
     ylim = c(0,800),
     main = 'Methylation')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_methyl[[4]]) ~ unlist(rf_methyl[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_methyl[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))


# plot test age sample collection
plot(unlist(rf_methyl[[4]]), unlist(rf_methyl[[8]]), 
     xlab = 'Test Predictions', 
     ylab = 'Test Age Sample Collection',
     xlim = c(0,800),
     ylim = c(0,800),
     main = 'Methylation')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_methyl[[4]]) ~ unlist(rf_methyl[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_methyl[[11]]))
# legend("topright", legend = paste0('r_squared = ', r_squared))

###############################################################################################
# just methylation log
rf_methyl_log <- predictAll(data = full_data_rf,
                        subset <- c("age_diagnosis", "age_sample_collection"),
                        selected_features = NULL, 
                        log = T,
                        cutoff = 0.15)


###############################
# plot test age diagnosis
plot(unlist(rf_methyl_log[[4]]), unlist(rf_methyl_log[[6]]), 
     xlab = 'Test Predictions', 
     ylab = 'Test Age of Diagnosis',
     xlim = c(0,8),
     ylim = c(0,8),
     main = 'Methylation Log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_methyl_log[[4]]) ~ unlist(rf_methyl_log[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_methyl_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


# plot test age sample collection
plot(unlist(rf_methyl_log[[4]]), unlist(rf_methyl_log[[8]]), 
     xlab = 'Test Predictions', 
     ylab = 'Test Age Sample Collection',
     xlim = c(0,8),
     ylim = c(0,8),
     main = 'Methylation Log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_methyl_log[[4]]) ~ unlist(rf_methyl_log[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_methyl_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))

