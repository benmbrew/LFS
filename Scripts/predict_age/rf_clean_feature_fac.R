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

# Load in clinical data
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# Read in 3 different data sets 
full_data <- read.csv(paste0(data_folder, '/full_data_rf.csv'), stringsAsFactors = F)


full_data$X <- NULL


# make categroical variable from age of methylaion and age of sample collection
full_data$age_diagnosis_fac <- as.integer(ifelse(full_data$age_diagnosis <= 48, 1, 2))

full_data$age_sample_fac <- as.integer(ifelse(full_data$age_sample_collection <= 48, 1, 2))


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
    
    y = make.names(as.factor(data$age_diagnosis_fac[train_index]))
    
    control <- trainControl(method="repeatedcv", 
                            number=2, 
                            repeats=1,
                            classProbs=TRUE,
                            summaryFunction=twoClassSummary # Use AUC to pick the best model
    )
    
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
    best <- best[order(best$X1, decreasing = T),]
    
    # subset data by top features 
    final <- best[best$X1 > 35,]
    final_genes <- final$gene
    
    # mtry: Number of variables randomly sampled as candidates at each split.
    # ntree: Number of trees to grow.
    
    mtry <- sqrt(ncol(data[train_index, c(selected_features, final_genes)]))
    tunegrid <- expand.grid(.mtry=mtry)
    
    
    
    if(svm) {
      
      ctrl <- trainControl(
        
        method="repeatedcv", # 10fold cross validation
        repeats=5, # do 5 repititions of cv
        summaryFunction=twoClassSummary, # Use AUC to pick the best model
        classProbs=TRUE
      )
      
      #Train and Tune the SVM
      model[[i]] <- train(x=data[train_index, c(selected_features, final_genes)], 
                          y = y,
                          method = "svmRadial",
                          tuneLength = 9, # 9 values of the cost function
                          metric="ROC",
                          trControl=ctrl) 
    } else {
      
      NFOLDS <- 2
      fitControl <- trainControl( 
        
        method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
        number = min(10, NFOLDS),      
        repeats = 1,
        classProbs=TRUE,
        summaryFunction=twoClassSummary, # Use AUC to pick the best model
        allowParallel = TRUE
        
      )
      model[[i]] <- train(x = data[train_index, c(selected_features, final_genes)]
                          , y = y
                          , method = "rf"
                          , trControl = fitControl
                          , tuneGrid = tunegrid
                          , importance = T
                          , metric="ROC"
                          , verbose = FALSE)
      
      temp <- varImp(model[[i]])[[1]]
      importance[[i]] <- cbind(rownames(temp), temp$Overall)
      
    }
     

      test.predictions[[i]] <- predict(model[[i]] 
                                       , newdata = data[-train_index, c(selected_features, final_genes)]
                                       , type = "prob")
      
      train.predictions[[i]] <- predict(model[[i]] 
                                        , newdata = data[train_index, c(selected_features, final_genes)]
                                        ,type = "prob")
      
  
    

      train.ground_truth[[i]] <- as.factor(make.names(data$age_diagnosis_fac[train_index]))
      test.ground_truth[[i]] <- as.factor(make.names(data$age_diagnosis_fac[-train_index]))
      train.sample_collection[[i]] = as.factor(make.names(data$age_sample_fac[train_index]))
      train.sample_collection[[i]] <- factor(train.sample_collection[[i]], levels = c("X1", "X2"))
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
# 16 - best_features

#######################################################################################
# Methylation
#######################################################################################

###################
# classification
###################

# age of diagnosis, regression, not log
methyl_fac <- predictAll(data = full_data,
                         svm = T,
                         log = F,
                         subset = c('age_diagnosis_fac', 'age_sample_fac'),
                         selected_features = NULL,
                         cutoff = .7,
                         iterations = 15)

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


# confustion matrix age of diagnosis 
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

# age of diagnosis, classification, log
methyl_fac_log <- predictAll(data = full_data_rf,
                             fac = T,
                             clin_only =  F,
                             clin_methyl = F,
                             log = T,
                             subset = c('age_diagnosis_fac', 'age_sample_fac'),
                             selected_features = NULL,
                             cutoff = .7,
                             resid = F,
                             iterations = 10)

# test acc for age of diagnosis
mean(unlist(methyl_fac_log[[9]]))

# test acc for age of sample collection
mean(unlist(methyl_fac_log[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- methyl_fac_log[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations


# confustion matrix age of diagnosis 
temp <- list()
for (i in 1:10){
  temp[[i]] <- methyl_fac_log[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations

