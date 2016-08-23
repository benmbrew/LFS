#############################################################################################################
# This script will be a clean version of the random forest function predicting age of diagnosis (log and not log) from 
# methylation. The script will also be able to add in clinical variables, and just do clinical if needed. It will also have
# an option for using residuals as predictors
# 1) clin and methylation
# 2) class and regresson
# 3) log and not log
# 4) residual and not residual
# 5) add in WT for clin

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
test <- paste0(project_folder, '/Scripts/Analyze')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

# Read in 3 different data sets 
full_data <- read.csv(paste0(data_folder, '/full_data.csv'), stringsAsFactors = F)


full_data$X <- NULL

# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       svm,
                       subset, 
                       selected_features,
                       cutoff,
                       log,
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
  
  
  # set log transformation
  if(log) {
    
    data[,c(6,8, 30:(ncol(data) - 2))]  <- log(data[,c(6,8,30:(ncol(data) -2))])
  }

  genes <- colnames(data)[30:(ncol(data) - 2)]
  data <- data[, c(subset, genes)]
  data <- data[!(is.na(data$age_diagnosis)),]
  

  
  obs <- nrow(data)
  
  
  for (i in 1:iterations){
    
    set.seed(i)
    train_index <- sample(nrow(data), nrow(data) *cutoff, replace = F)
    
    y = data$age_diagnosis[train_index]
    
    
  
    control <- trainControl(method="repeatedcv", number=2, repeats=1)
    
    # train the model
    best_model <- train(x = data[train_index, c(selected_features, genes)],
                   y = y,
                   importance = TRUE,
                   trControl = control)
    
    # estimate variable best
    best_features[[i]] <- varImp(best_model)
    
    # get vector of best 
    best <- best_features[[i]]$importance
    
    # make rownames a column
    best$gene <- rownames(best)
    rownames(best) <- NULL
    
    # sort best vector
    best <- best[order(best$Overall, decreasing = T),]
    
    # subset data by top features 
    final <- best[best$Overall > 35,]
    final_genes <- final$gene
    
    # mtry: Number of variables randomly sampled as candidates at each split.
    # ntree: Number of trees to grow.
    
    mtry <- sqrt(ncol(data[train_index, c(selected_features, final_genes)]))
    tunegrid <- expand.grid(.mtry=mtry)
    
  
    
    if(svm) {
      
      ctrl <- trainControl(
        
        method="repeatedcv", # 10fold cross validation
        repeats=5 # do 5 repititions of cv
                           # summaryFunction=twoClassSummary # Use AUC to pick the best model
        # classProbs=TRUE)
      )
      
      #Train and Tune the SVM
      model[[i]] <- train(x=data[train_index, c(selected_features, final_genes)], 
                          y = y,
                          method = "svmRadial",
                          tuneLength = 9, # 9 values of the cost function
                          preProc = c("center","scale"),
                          trControl=ctrl) 
    } else {
      
      NFOLDS <- 2
      fitControl <- trainControl( 
        
        method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
        number = min(10, NFOLDS),      
        repeats = 1,
        allowParallel = TRUE
        
      )
      
      model[[i]] <- train(x = data[train_index, c(selected_features, final_genes)]
                          , y = y
                          , method = "rf"
                          , trControl = fitControl
                          , tuneGrid = tunegrid
                          , importance = T
                          , verbose = FALSE)
      
      temp <- varImp(model[[i]])[[1]]
      importance[[i]] <- cbind(rownames(temp), temp$Overall)
    }
    
    
   
    
    
    test.predictions[[i]] <- predict(model[[i]] 
                                     , newdata = data[-train_index, c(selected_features, final_genes)])
    
    train.predictions[[i]] <- predict(model[[i]] 
                                      , newdata = data[train_index, c(selected_features, final_genes)])
    
  
  
  
      
    train.ground_truth[[i]] <- data$age_diagnosis[train_index]
    test.ground_truth[[i]] <- data$age_diagnosis[-train_index]
    train.sample_collection[[i]] = data$age_sample_collection[train_index]
    test.sample_collection[[i]] = data$age_sample_collection[-train_index]
    train.mse[[i]] <- rmse(unlist(train.predictions[[i]]), unlist(train.ground_truth[[i]]))
    test.mse[[i]] <- rmse(unlist(test.predictions[[i]]), unlist(test.ground_truth[[i]]))
      
    
    
    print(i)
    
  }
  
  return(list(train.mse, test.mse,  train.predictions, test.predictions, train.ground_truth, test.ground_truth, train.sample_collection,
              test.sample_collection, test_acc, test_stats, test_acc_samp, test_stats_samp, model, importance, obs, best_features))
  
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
# 16 - best_model

#######################################################################################
# Methylation
#######################################################################################

###################
# Regression
###################

# age of diagnosis, regression, not log
methyl_reg <- predictAll(data = full_data,
                         svm = F,
                         log = F,
                         subset = c('age_diagnosis', 'age_sample_collection'),
                         selected_features = NULL,
                         cutoff = .7,
                         iterations = 15)

# plot predictions against ground truth
plot(unlist(methyl_reg[[4]]), unlist(methyl_reg[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'Age of Diagnosis (Months)')
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
     main = 'Age of Sample Collection (Months)')
abline(0,1)
r_squared <- round(summary(lm(unlist(methyl_reg[[4]]) ~ unlist(methyl_reg[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', methyl_reg[[15]]), cex = 0.7)



# age of diagnosis, regression with log transform
methyl_reg_log <- predictAll(data = full_data,
                             svm = F,
                             log = T,
                             subset = c('age_diagnosis', 'age_sample_collection'),
                             selected_features = NULL,
                             cutoff = .7,
                             iterations = 15)

# plot predictions against ground truth
plot(unlist(methyl_reg_log[[4]]), unlist(methyl_reg_log[[6]]), 
     xlim = c(0, 8),
     ylim = c(0, 8),
     xlab = 'Log Predictions',
     ylab = 'Log Real Age of Diagnosis',
     main = 'Log Age of Diagnosis (Months)')
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
     main = 'Log Age of Sample Collection (Months)')
abline(0,1)
r_squared <- round(summary(lm(unlist(methyl_reg_log[[4]]) ~ unlist(methyl_reg_log[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', methyl_reg_log[[15]]), cex = 0.7)

# save.image(file = '/home/benbrew/Desktop/rf_clean_feature.RData')
