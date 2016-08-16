#############################################################################################################
# This script will be a clean version of the ridge regression function predicting age of diagnosis (log and not log) from 
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
full_data <- read.csv(paste0(data_folder, '/full_data.csv'), stringsAsFactors = F)
full_data_cor <- read.csv(paste0(data_folder, '/full_data_cor.csv'), stringsAsFactors = F)
full_data_rf <- read.csv(paste0(data_folder, '/full_data_rf.csv'), stringsAsFactors = F)
# Read in 2 smaller cor and rf
full_data_cor_small <- read.csv(paste0(data_folder, '/full_data_cor_small.csv'), stringsAsFactors = F)
full_data_rf_small <- read.csv(paste0(data_folder, '/full_data_rf_small.csv'), stringsAsFactors = F)


full_data$X <- NULL
full_data_cor$X <- NULL
full_data_rf$X <- NULL
full_data_cor$X <- NULL
full_data_rf_small$X <- NULL

# Read in residuals from regressing gene on age of sample collection
resid_rf <- read.csv(paste0(data_folder, '/resid_rf.csv'), stringsAsFactors = F)
resid_cor <- read.csv(paste0(data_folder, '/resid_cor.csv'), stringsAsFactors = F)
resid_full <- read.csv(paste0(data_folder, '/resid_full.csv'), stringsAsFactors = F)

# make categorifcal variable for clinical data
clin$age_diagnosis_fac <- as.integer(ifelse(clin$age_diagnosis <= 48, 1, 2))

clin$age_sample_fac <- as.integer(ifelse(clin$age_sample_collection <= 48, 1, 2))


# make categroical variable from age of methylaion and age of sample collection
full_data_rf$age_diagnosis_fac <- as.integer(ifelse(full_data_rf$age_diagnosis <= 48, 1, 2))

full_data_rf$age_sample_fac <- as.integer(ifelse(full_data_rf$age_sample_collection <= 48, 1, 2))

# make categroical variable from age of methylaion and age of sample collection
full_data_cor$age_diagnosis_fac <- as.integer(ifelse(full_data_cor$age_diagnosis <= 48, 1, 2))

full_data_cor$age_sample_fac <- as.integer(ifelse(full_data_cor$age_sample_collection <= 48, 1, 2))

# make categroical variable from age of methylaion and age of sample collection
full_data_rf_small$age_diagnosis_fac <- as.integer(ifelse(full_data_rf_small$age_diagnosis <= 48, 1, 2))

full_data_rf_small$age_sample_fac <- as.integer(ifelse(full_data_rf_small$age_sample_collection <= 48, 1, 2))

# make categroical variable from age of methylaion and age of sample collection
full_data_cor_small$age_diagnosis_fac <- as.integer(ifelse(full_data_rf$age_diagnosis <= 48, 1, 2))

full_data_cor_small$age_sample_fac <- as.integer(ifelse(full_data_rf$age_sample_collection <= 48, 1, 2))

# make categroical variable from age of methylaion and age of sample collection in residual data
resid_rf$age_diagnosis_fac <- as.integer(ifelse(resid_rf$age_diagnosis <= 48, 1, 2))

resid_rf$age_sample_fac <- as.integer(ifelse(resid_rf$age_sample_collection <= 48, 1, 2))

resid_cor$age_diagnosis_fac <- as.integer(ifelse(resid_cor$age_diagnosis <= 48, 1, 2))

resid_cor$age_sample_fac <- as.integer(ifelse(resid_cor$age_sample_collection <= 48, 1, 2))

resid_full$age_diagnosis_fac <- as.integer(ifelse(resid_full$age_diagnosis <= 48, 1, 2))

resid_full$age_sample_fac <- as.integer(ifelse(resid_full$age_sample_collection <= 48, 1, 2))



# remove variable 
resid_rf$X <- NULL
resid_cor$X <- NULL
resid_full$X <- NULL



# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       clin_only,
                       clin_methyl,
                       fac,
                       subset, 
                       selected_features,
                       cutoff,
                       log,
                       resid,
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
  if(log & !resid) {
    
    data[,c(6,8, 30:(ncol(data) - 2))]  <- log(data[,c(6,8,30:(ncol(data) -2))])
  }
  
  if(log & resid) {
    
    data[,c(1:(ncol(data) - 2))]  <- log(data[,c(1:(ncol(data) -2))])
  }
  
  
  if(clin_only) {
    genes <- NULL
    data <- data[, c(subset, genes)]
    data <- data[complete.cases(data),]
  } else if (clin_methyl) {
    genes <- colnames(data)[30:(ncol(data) - 2)]
    data <- data[, c(subset, genes)]
    data <- data[complete.cases(data),]
  } else {
    genes <- colnames(data)[30:(ncol(data) - 2)]
    data <- data[, c(subset, genes)]
    data <- data[!(is.na(data$age_diagnosis)),]
  }
  
  
  
  
  for ( i in 3:ncol(data)){
    
    if(typeof(data[,i]) == 'character' || typeof(data[,i]) == 'integer') {
      data[,i] <- as.factor(data[,i])
    }
  }
  
  
  
  obs <- nrow(data)
  
  
  for (i in 1:iterations){
    
    set.seed(i)
    train_index <- sample(nrow(data), nrow(data) *cutoff, replace = F)
    
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

#######################################################################################
# Methylation
#######################################################################################

###################
# Regression
###################

# age of diagnosis, regression, not log
methyl_reg <- predictAll(data = full_data_rf,
                         fac = F,
                         clin_only =  F,
                         clin_methyl = F,
                         log = F,
                         subset = c('age_diagnosis', 'age_sample_collection'),
                         selected_features = NULL,
                         cutoff = .7,
                         resid = F,
                         iterations = 10)

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
methyl_reg_log <- predictAll(data = full_data_rf,
                             fac = F,
                             clin_only =  F,
                             clin_methyl = F,
                             log = T,
                             subset = c('age_diagnosis', 'age_sample_collection'),
                             selected_features = NULL,
                             cutoff = .7,
                             resid = F,
                             iterations = 10)

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

###################
# classification
###################

# age of diagnosis, classification, not log
methyl_fac <- predictAll(data = full_data_rf,
                         fac = T,
                         clin_only =  F,
                         clin_methyl = F,
                         log = F,
                         subset = c('age_diagnosis_fac', 'age_sample_fac'),
                         selected_features = NULL,
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



#####################################################################
# Regression resid
#####################################################################

# age of diagnosis, regression, not log
methyl_reg_resid <- predictAll(data = resid_rf,
                               fac = F,
                               clin_only =  F,
                               clin_methyl = F,
                               log = F,
                               subset = c('age_diagnosis', 'age_sample_collection'),
                               selected_features = NULL,
                               cutoff = .7,
                               resid = T,
                               iterations = 10)

# plot predictions against ground truth
plot(unlist(methyl_reg_resid[[4]]), unlist(methyl_reg_resid[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions with residuals',
     ylab = 'Real Age of Diagnosis',
     main = 'Age of Diagnosis (Months)')
abline(0,1)
r_squared <- round(summary(lm(unlist(methyl_reg_resid[[4]]) ~ unlist(methyl_reg_resid[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', methyl_reg_resid[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(methyl_reg_resid[[4]]), unlist(methyl_reg_resid[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions with residuals',
     ylab = 'Real Age of Sample Collection',
     main = 'Age of Sample Collection (Months)')
abline(0,1)
r_squared <- round(summary(lm(unlist(methyl_reg_resid[[4]]) ~ unlist(methyl_reg_resid[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', methyl_reg_resid[[15]]), cex = 0.7)


###################
# classification
###################

# age of diagnosis, classification, not log
methyl_fac_resid <- predictAll(data = resid_rf,
                               fac = T,
                               clin_only =  F,
                               clin_methyl = F,
                               log = F,
                               subset = c('age_diagnosis_fac', 'age_sample_fac'),
                               selected_features = NULL,
                               cutoff = .7,
                               resid = T,
                               iterations = 10)

# test acc for age of diagnosis
mean(unlist(methyl_fac_resid[[9]]))

# test acc for age of sample collection
mean(unlist(methyl_fac_resid[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- methyl_fac_resid[[10]][[i]]$table
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
  temp[[i]] <- methyl_fac_resid[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations


############################################################################
# combine methyl and clin
############################################################################

# best clinical variables (accuracy) - gdna.exon.intron(fac) 65, mdm2 60, p53 62

###################
# Regression
###################

# age of diagnosis, regression, not log gdna.exon.intron
methyl_gdna.exon.intron <- predictAll(data = full_data_rf,
                                      fac = F,
                                      clin_only =  F,
                                      clin_methyl = T,
                                      log = F,
                                      subset = c('age_diagnosis', 'age_sample_collection', 'gdna.exon.intron.'),
                                      selected_features = c('gdna.exon.intron'),
                                      cutoff = .7,
                                      resid = F,
                                      iterations = 10)

# plot predictions against ground truth
plot(unlist(methyl_gdna.exon.intron[[4]]), unlist(methyl_gdna.exon.intron[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'Methylation and gdna.exon.intron')
abline(0,1)
r_squared <- round(summary(lm(unlist(methyl_gdna.exon.intron[[4]]) ~ unlist(methyl_gdna.exon.intron[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', methyl_gdna.exon.intron[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(methyl_gdna.exon.intron[[4]]), unlist(methyl_gdna.exon.intron[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'Methylation and gdna.exon.intron')
abline(0,1)
r_squared <- round(summary(lm(unlist(methyl_gdna.exon.intron[[4]]) ~ unlist(methyl_gdna.exon.intron[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', methyl_gdna.exon.intron[[15]]), cex = 0.7)


# age of diagnosis, regression, not log
methyl_gdna.exon.intron_log <- predictAll(data = full_data_rf,
                                          fac = F,
                                          clin_only =  F,
                                          clin_methyl = T,
                                          log = T,
                                          subset = c('age_diagnosis', 'age_sample_collection', 'gdna.exon.intron'),
                                          selected_features = c('gdna.exon.intron'),
                                          cutoff = .7,
                                          resid = F,
                                          iterations = 10)

# plot predictions against ground truth
plot(unlist(methyl_gdna.exon.intron_log[[4]]), unlist(methyl_gdna.exon.intron_log[[6]]), 
     xlim = c(0, 8),
     ylim = c(0, 8),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'Methylation and gdna.exon.intron Log')
abline(0,1)
r_squared <- round(summary(lm(unlist(methyl_gdna.exon.intron_log[[4]]) ~ unlist(methyl_gdna.exon.intron_log[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', methyl_gdna.exon.intron_log[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(methyl_gdna.exon.intron_log[[4]]), unlist(methyl_gdna.exon.intron_log[[8]]), 
     xlim = c(0, 8),
     ylim = c(0, 8),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'Methylation and gdna.exon.intron Log')
abline(0,1)
r_squared <- round(summary(lm(unlist(methyl_gdna.exon.intron_log[[4]]) ~ unlist(methyl_gdna.exon.intron_log[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', methyl_gdna.exon.intron_log[[15]]), cex = 0.7)



###################
# classification
###################

# age of diagnosis, classification, not log
methyl_fac <- predictAll(data = full_data_rf,
                         fac = T,
                         clin_only =  F,
                         clin_methyl = F,
                         log = F,
                         subset = c('age_diagnosis_fac', 'age_sample_fac'),
                         selected_features = NULL,
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




#######################################################################################
# Clincial
#######################################################################################

##########################
# Regression
##########################
# age of diagnosis, classification, not log - using gender and p53_germline
mod_p53 <- predictAll(data = clin,
                      fac = F,
                      clin_only =  T,
                      clin_methyl = F,
                      log = F,
                      subset = c('age_diagnosis', 'age_sample_collection', 'gender', 'p53_germline'),
                      selected_features = c('p53_germline', 'gender'),
                      cutoff = .7,
                      resid = F,
                      iterations = 10)

# plot predictions against ground truth
plot(unlist(mod_p53[[4]]), unlist(mod_p53[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'p53_germline')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_p53[[4]]) ~ unlist(mod_p53[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_p53[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(mod_p53[[4]]), unlist(mod_p53[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'p53_germline')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_p53[[4]]) ~ unlist(mod_p53[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_p53[[15]]), cex = 0.7)


# age of diagnosis, classification, not log - using gender and gdna.exon.intron
mod_gdna.exon.intron <- predictAll(data = clin,
                                   fac = F,
                                   clin_only =  T,
                                   clin_methyl = F,
                                   log = F,
                                   subset = c('age_diagnosis', 'age_sample_collection', 'gender', 'gdna.exon.intron'),
                                   selected_features = c('gdna.exon.intron', 'gender'),
                                   cutoff = .7,
                                   resid = F,
                                   iterations = 10)

# plot predictions against ground truth
plot(unlist(mod_gdna.exon.intron[[4]]), unlist(mod_gdna.exon.intron[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'gdna.exon.intron')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_gdna.exon.intron[[4]]) ~ unlist(mod_gdna.exon.intron[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_gdna.exon[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(mod_gdna.exon.intron[[4]]), unlist(mod_gdna.exon.intron[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'gdna.exon.intron')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_gdna.exon.intron[[4]]) ~ unlist(mod_gdna.exon.intron[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_gdna.exon.intron[[15]]), cex = 0.7)


# age of diagnosis, classification, not log - using gender and gdna.exon.intron.fac
mod_gdna.exon.intron.fac <- predictAll(data = clin,
                                       fac = F,
                                       clin_only =  T,
                                       clin_methyl = F,
                                       log = F,
                                       subset = c('age_diagnosis', 'age_sample_collection', 'gender', 'gdna.exon.intron.fac'),
                                       selected_features = c('gdna.exon.intron.fac', 'gender'),
                                       cutoff = .7,
                                       resid = F,
                                       iterations = 10)

# plot predictions against ground truth
plot(unlist(mod_gdna.exon.intron.fac[[4]]), unlist(mod_gdna.exon.intron.fac[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'gdna.exon.intron.fac')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_gdna.exon.intron.fac[[4]]) ~ unlist(mod_gdna.exon.intron.fac[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_gdna.exon.intron.fac[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(mod_gdna.exon.intron.fac[[4]]), unlist(mod_gdna.exon.intron.fac[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'gdna.exon.intron.fac')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_gdna.exon.intron.fac[[4]]) ~ unlist(mod_gdna.exon.intron.fac[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_gdna.exon.intron.fac[[15]]), cex = 0.7)




# age of diagnosis, classification, not log - using gender and gdna.base.change
mod_gdna.base.change <- predictAll(data = clin,
                                   fac = F,
                                   clin_only =  T,
                                   clin_methyl = F,
                                   log = F,
                                   subset = c('age_diagnosis', 'age_sample_collection', 'gender', 'gdna.base.change'),
                                   selected_features = c('gdna.base.change', 'gender'),
                                   cutoff = .7,
                                   resid = F,
                                   iterations = 10)

# plot predictions against ground truth
plot(unlist(mod_gdna.base.change[[4]]), unlist(mod_gdna.base.change[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'gdna.base.change')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_gdna.base.change[[4]]) ~ unlist(mod_gdna.base.change[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_gdna.base.change[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(mod_gdna.base.change[[4]]), unlist(mod_gdna.base.change[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'gdna.base.change')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_gdna.base.change[[4]]) ~ unlist(mod_gdna.base.change[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_gdna.base.change[[15]]), cex = 0.7)



# age of diagnosis, classification, not log - using gender and gdna.base.change.fac
mod_gdna.base.change.fac <- predictAll(data = clin,
                                       fac = F,
                                       clin_only =  T,
                                       clin_methyl = F,
                                       log = F,
                                       subset = c('age_diagnosis', 'age_sample_collection', 'gender', 'gdna.base.change.fac'),
                                       selected_features = c('gdna.base.change.fac', 'gender'),
                                       cutoff = .7,
                                       resid = F,
                                       iterations = 10)

# plot predictions against ground truth
plot(unlist(mod_gdna.base.change.fac[[4]]), unlist(mod_gdna.base.change.fac[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'gdna.base.change.fac')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_gdna.base.change.fac[[4]]) ~ unlist(mod_gdna.base.change.fac[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_gdna.base.change.fac[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(mod_gdna.base.change.fac[[4]]), unlist(mod_gdna.base.change.fac[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'gdna.base.change.fac')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_gdna.base.change.fac[[4]]) ~ unlist(mod_gdna.base.change.fac[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_gdna.base.change.fac[[15]]), cex = 0.7)


#HERE
# age of diagnosis, classification, not log - using gender and gdna.base.change.fac
mod_gdna.codon<- predictAll(data = clin,
                            fac = F,
                            clin_only =  T,
                            clin_methyl = F,
                            log = F,
                            subset = c('age_diagnosis', 'age_sample_collection', 'gender', 'gdna.codon'),
                            selected_features = c('gdna.codon', 'gender'),
                            cutoff = .7,
                            resid = F,
                            iterations = 10)

# plot predictions against ground truth
plot(unlist(mod_gdna.codon[[4]]), unlist(mod_gdna.codon[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'gdna.codon')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_gdna.codon[[4]]) ~ unlist(mod_gdna.codon[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_gdna.codon[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(mod_gdna.codon[[4]]), unlist(mod_gdna.codon[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'gdna.codon')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_gdna.codon[[4]]) ~ unlist(mod_gdna.codon[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_gdna.codon[[15]]), cex = 0.7)




# age of diagnosis, classification, not log - using gender and protein.codon.change
mod_protein.codon.change <- predictAll(data = clin,
                                       fac = F,
                                       clin_only =  T,
                                       clin_methyl = F,
                                       log = F,
                                       subset = c('age_diagnosis', 'age_sample_collection', 'gender', 'protein.codon.change'),
                                       selected_features = c('protein.codon.change', 'gender'),
                                       cutoff = .7,
                                       resid = F,
                                       iterations = 10)

# plot predictions against ground truth
plot(unlist(mod_protein.codon.change[[4]]), unlist(mod_protein.codon.change[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'protein.codon.change')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_protein.codon.change[[4]]) ~ unlist(mod_protein.codon.change[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_protein.codon.change[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(mod_protein.codon.change[[4]]), unlist(mod_protein.codon.change[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'protein.codon.change')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_protein.codon.change[[4]]) ~ unlist(mod_protein.codon.change[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_protein.codon.change[[15]]), cex = 0.7)



# age of diagnosis, classification, not log - using gender and protein.codon.change.fac
mod_protein.codon.change.fac <- predictAll(data = clin,
                                           fac = F,
                                           clin_only =  T,
                                           clin_methyl = F,
                                           log = F,
                                           subset = c('age_diagnosis', 'age_sample_collection', 'gender', 'protein.codon.change.fac'),
                                           selected_features = c('protein.codon.change.fac', 'gender'),
                                           cutoff = .7,
                                           resid = F,
                                           iterations = 10)

# plot predictions against ground truth
plot(unlist(mod_protein.codon.change.fac[[4]]), unlist(mod_protein.codon.change.fac[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'protein.codon.change.fac')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_protein.codon.change.fac[[4]]) ~ unlist(mod_protein.codon.change.fac[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_protein.codon.change.fac[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(mod_protein.codon.change.fac[[4]]), unlist(mod_protein.codon.change.fac[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'protein.codon.change.fac')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_protein.codon.change.fac[[4]]) ~ unlist(mod_protein.codon.change.fac[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_protein.codon.change.fac[[15]]), cex = 0.7)





# age of diagnosis, classification, not log - using gender and protein.codon.change
mod_protein.codon.num <- predictAll(data = clin,
                                    fac = F,
                                    clin_only =  T,
                                    clin_methyl = F,
                                    log = F,
                                    subset = c('age_diagnosis', 'age_sample_collection', 'gender', 'protein.codon.num'),
                                    selected_features = c('protein.codon.num', 'gender'),
                                    cutoff = .7,
                                    resid = F,
                                    iterations = 10)

# plot predictions against ground truth
plot(unlist(mod_protein.codon.num[[4]]), unlist(mod_protein.codon.num[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'protein.codon.num')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_protein.codon.num[[4]]) ~ unlist(mod_protein.codon.num[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_protein.codon.num[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(mod_protein.codon.num[[4]]), unlist(mod_protein.codon.num[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'protein.codon.num')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_protein.codon.num[[4]]) ~ unlist(mod_protein.codon.num[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_protein.codon.num[[15]]), cex = 0.7)




# age of diagnosis, classification, not log - using gender and protein.codon.change
mod_splice.delins.snv <- predictAll(data = clin,
                                    fac = F,
                                    clin_only =  T,
                                    clin_methyl = F,
                                    log = F,
                                    subset = c('age_diagnosis', 'age_sample_collection', 'gender', 'splice.delins.snv'),
                                    selected_features = c('splice.delins.snv', 'gender'),
                                    cutoff = .7,
                                    resid = F,
                                    iterations = 10)

# plot predictions against ground truth
plot(unlist(mod_splice.delins.snv[[4]]), unlist(mod_splice.delins.snv[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'splice.delins.snv')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_splice.delins.snv[[4]]) ~ unlist(mod_splice.delins.snv[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_splice.delins.snv[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(mod_splice.delins.snv[[4]]), unlist(mod_splice.delins.snv[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'splice.delins.snv')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_splice.delins.snv[[4]]) ~ unlist(mod_splice.delins.snv[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_splice.delins.snv[[15]]), cex = 0.7)




# age of diagnosis, classification, not log - using gender and protein.codon.change
mod_codon72.npro <- predictAll(data = clin,
                               fac = F,
                               clin_only =  T,
                               clin_methyl = F,
                               log = F,
                               subset = c('age_diagnosis', 'age_sample_collection', 'gender', 'codon72.npro'),
                               selected_features = c('codon72.npro', 'gender'),
                               cutoff = .7,
                               resid = F,
                               iterations = 10)

# plot predictions against ground truth
plot(unlist(mod_codon72.npro[[4]]), unlist(mod_codon72.npro[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'codon72.npro')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_codon72.npro[[4]]) ~ unlist(mod_codon72.npro[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_codon72.npro[[15]]), cex = 0.7)


# plot predictions against ground truth
plot(unlist(mod_codon72.npro[[4]]), unlist(mod_codon72.npro[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'codon72.npro')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_codon72.npro[[4]]) ~ unlist(mod_codon72.npro[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_codon72.npro[[15]]), cex = 0.7)



# age of diagnosis, classification, not log - using gender and protein.codon.change
mod_mdm2.nG <- predictAll(data = clin,
                          fac = F,
                          clin_only =  T,
                          clin_methyl = F,
                          log = F,
                          subset = c('age_diagnosis', 'age_sample_collection', 'gender', 'mdm2.nG'),
                          selected_features = c('mdm2.nG', 'gender'),
                          cutoff = .7,
                          resid = F,
                          iterations = 10)

# plot predictions against ground truth
plot(unlist(mod_mdm2.nG[[4]]), unlist(mod_mdm2.nG[[6]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Diagnosis',
     main = 'mdm2.nG')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_mdm2.nG[[4]]) ~ unlist(mod_mdm2.nG[[6]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_mdm2.nG[[15]]), cex = 0.7)

# plot predictions against ground truth
plot(unlist(mod_mdm2.nG[[4]]), unlist(mod_mdm2.nG[[8]]), 
     xlim = c(0, 1000),
     ylim = c(0, 1000),
     xlab = 'Predictions',
     ylab = 'Real Age of Sample Collection',
     main = 'mdm2.nG')
abline(0,1)
r_squared <- round(summary(lm(unlist(mod_mdm2.nG[[4]]) ~ unlist(mod_mdm2.nG[[8]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared), cex = 0.7)
legend("bottomright", legend = paste0('# obs = ', mod_mdm2.nG[[15]]), cex = 0.7)



##########################
# Classification
##########################


# age of diagnosis, classification, not log - using gender and p53_germline
mod_p53 <- predictAll(data = clin,
                      fac = T,
                      clin_only =  T,
                      clin_methyl = F,
                      log = F,
                      subset = c('age_diagnosis_fac', 'age_sample_fac', 'gender', 'p53_germline'),
                      selected_features = c('gender', 'p53_germline'),
                      cutoff = .7,
                      resid = F,
                      iterations = 10)

# test acc for age of diagnosis
mean(unlist(mod_p53[[9]]))

# test acc for age of sample collection
mean(unlist(mod_p53[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_p53[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat)
# confustion matrix age of sample collection
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_p53[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat_sample)

# age of diagnosis, classification, not log - using gender and gdna.exon.intron
mod_gdna.exon.intron <- predictAll(data = clin,
                                   fac = T,
                                   clin_only =  T,
                                   clin_methyl = F,
                                   log = F,
                                   subset = c('age_diagnosis_fac', 'age_sample_fac', 'gender', 'gdna.exon.intron'),
                                   selected_features = c('gdna.exon.intron', 'gender'),
                                   cutoff = .7,
                                   resid = F,
                                   iterations = 10)

# test acc for age of diagnosis
mean(unlist(mod_gdna.exon.intron[[9]]))

# test acc for age of sample collection
mean(unlist(mod_gdna.exon.intron[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_gdna.exon.intron[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat)


# confustion matrix age of sample collection
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_gdna.exon.intron[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations


print(new_mat_sample)


# age of diagnosis, classification, not log - using gender and gdna.exon.intron.fac
mod_gdna.exon.intron.fac <- predictAll(data = clin,
                                       fac = T,
                                       clin_only =  T,
                                       clin_methyl = F,
                                       log = F,
                                       subset = c('age_diagnosis_fac', 'age_sample_fac', 'gender', 'gdna.exon.intron.fac'),
                                       selected_features = c('gdna.exon.intron.fac', 'gender'),
                                       cutoff = .7,
                                       resid = F,
                                       iterations = 10)

# test acc for age of diagnosis
mean(unlist(mod_gdna.exon.intron.fac[[9]]))

# test acc for age of sample collection
mean(unlist(mod_gdna.exon.intron.fac[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_gdna.exon.intron.fac[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat)

# confustion matrix age of sample collection
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_gdna.exon.intron.fac[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat_sample)


# age of diagnosis, classification, not log - using gender and gdna.base.change
mod_gdna.base.change <- predictAll(data = clin,
                                   fac = T,
                                   clin_only =  T,
                                   clin_methyl = F,
                                   log = F,
                                   subset = c('age_diagnosis_fac', 'age_sample_fac', 'gender', 'gdna.base.change'),
                                   selected_features = c('gdna.base.change', 'gender'),
                                   cutoff = .7,
                                   resid = F,
                                   iterations = 10)


# test acc for age of diagnosis
mean(unlist(mod_gdna.base.change[[9]]))

# test acc for age of sample collection
mean(unlist(mod_gdna.base.change[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_gdna.base.change[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat)

# confustion matrix age of sample collection
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_gdna.base.change[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat_sample)

# age of diagnosis, classification, not log - using gender and gdna.base.change.fac
mod_gdna.base.change.fac <- predictAll(data = clin,
                                       fac = T,
                                       clin_only =  T,
                                       clin_methyl = F,
                                       log = F,
                                       subset = c('age_diagnosis_fac', 'age_sample_fac', 'gender', 'gdna.base.change.fac'),
                                       selected_features = c('gdna.base.change.fac', 'gender'),
                                       cutoff = .7,
                                       resid = F,
                                       iterations = 10)


# test acc for age of diagnosis
mean(unlist(mod_gdna.base.change.fac[[9]]))

# test acc for age of sample collection
mean(unlist(mod_gdna.base.change.fac[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_gdna.base.change.fac[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat)


# confustion matrix age of sample collection
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_gdna.base.change.fac[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat_sample)


# age of diagnosis, classification, not log - using gender and gdna.base.change.fac
mod_gdna.codon<- predictAll(data = clin,
                            fac = T,
                            clin_only =  T,
                            clin_methyl = F,
                            log = F,
                            subset = c('age_diagnosis_fac', 'age_sample_fac', 'gender', 'gdna.codon'),
                            selected_features = c('gdna.codon', 'gender'),
                            cutoff = .7,
                            resid = F,
                            iterations = 10)


# test acc for age of diagnosis
mean(unlist(mod_gdna.codon[[9]]))

# test acc for age of sample collection
mean(unlist(mod_gdna.codon[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_gdna.codon[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat)

# confustion matrix age of sample collection
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_gdna.codon[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations


print(new_mat_sample)


# age of diagnosis, classification, not log - using gender and protein.codon.change
mod_protein.codon.change <- predictAll(data = clin,
                                       fac = T,
                                       clin_only =  T,
                                       clin_methyl = F,
                                       log = F,
                                       subset = c('age_diagnosis_fac', 'age_sample_fac', 'gender', 'protein.codon.change'),
                                       selected_features = c('protein.codon.change', 'gender'),
                                       cutoff = .7,
                                       resid = F,
                                       iterations = 10)


# test acc for age of diagnosis
mean(unlist(mod_protein.codon.change[[9]]))

# test acc for age of sample collection
mean(unlist(mod_protein.codon.change[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_protein.codon.change[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat)

# confustion matrix age of sample collection
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_protein.codon.change[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat_sample)

# age of diagnosis, classification, not log - using gender and protein.codon.change.fac
mod_protein.codon.change.fac <- predictAll(data = clin,
                                           fac = T,
                                           clin_only =  T,
                                           clin_methyl = F,
                                           log = F,
                                           subset = c('age_diagnosis_fac', 'age_sample_fac', 'gender', 'protein.codon.change.fac'),
                                           selected_features = c('protein.codon.change.fac', 'gender'),
                                           cutoff = .7,
                                           resid = F,
                                           iterations = 10)


# test acc for age of diagnosis
mean(unlist(mod_protein.codon.change.fac[[9]]))

# test acc for age of sample collection
mean(unlist(mod_protein.codon.change.fac[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_protein.codon.change.fac[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat)

# confustion matrix age of sample collection
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_protein.codon.change.fac[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat_sample)

# age of diagnosis, classification, not log - using gender and protein.codon.change
mod_protein.codon.num <- predictAll(data = clin,
                                    fac = T,
                                    clin_only =  T,
                                    clin_methyl = F,
                                    log = F,
                                    subset = c('age_diagnosis_fac', 'age_sample_fac', 'gender', 'protein.codon.num'),
                                    selected_features = c('protein.codon.num', 'gender'),
                                    cutoff = .7,
                                    resid = F,
                                    iterations = 10)

# test acc for age of diagnosis
mean(unlist(mod_protein.codon.num[[9]]))

# test acc for age of sample collection
mean(unlist(mod_protein.codon.num[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_protein.codon.num[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations
print(new_mat)


# confustion matrix age of sample collection
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_protein.codon.num[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat_sample)

# age of diagnosis, classification, not log - using gender and protein.codon.change
mod_splice.delins.snv <- predictAll(data = clin,
                                    fac = T,
                                    clin_only =  T,
                                    clin_methyl = F,
                                    log = F,
                                    subset = c('age_diagnosis_fac', 'age_sample_fac', 'gender', 'splice.delins.snv'),
                                    selected_features = c('splice.delins.snv', 'gender'),
                                    cutoff = .7,
                                    resid = F,
                                    iterations = 10)



# test acc for age of diagnosis
mean(unlist(mod_splice.delins.snv[[9]]))

# test acc for age of sample collection
mean(unlist(mod_splice.delins.snv[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_splice.delins.snv[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat)

# confustion matrix age of sample collection
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_splice.delins.snv[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat_sample)

# age of diagnosis, classification, not log - using gender and protein.codon.change
mod_codon72.npro <- predictAll(data = clin,
                               fac = T,
                               clin_only =  T,
                               clin_methyl = F,
                               log = F,
                               subset = c('age_diagnosis_fac', 'age_sample_fac', 'gender', 'codon72.npro'),
                               selected_features = c('codon72.npro', 'gender'),
                               cutoff = .7,
                               resid = F,
                               iterations = 10)



# test acc for age of diagnosis
mean(unlist(mod_codon72.npro[[9]]))

# test acc for age of sample collection
mean(unlist(mod_codon72.npro[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_codon72.npro[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat)

# confustion matrix age of sample collection
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_codon72.npro[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations
print(new_mat_sample)


# age of diagnosis, classification, not log - using gender and protein.codon.change
mod_mdm2.nG <- predictAll(data = clin,
                          fac = T,
                          clin_only =  T,
                          clin_methyl = F,
                          log = F,
                          subset = c('age_diagnosis_fac', 'age_sample_fac', 'gender', 'mdm2.nG'),
                          selected_features = c('mdm2.nG', 'gender'),
                          cutoff = .7,
                          resid = F,
                          iterations = 10)


# test acc for age of diagnosis
mean(unlist(mod_mdm2.nG[[9]]))

# test acc for age of sample collection
mean(unlist(mod_mdm2.nG[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_mdm2.nG[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations

print(new_mat)

# confustion matrix age of sample collection
temp <- list()
for (i in 1:10){
  temp[[i]] <- mod_mdm2.nG[[12]][[i]]$table
}
mat <- unlist(temp)
new_mat_sample <- matrix(, 2, 2)

new_mat_sample[1,1] <- sum(mat[mat_index])/iterations
new_mat_sample[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat_sample[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat_sample[2,2] <- sum(mat[mat_index + 3])/iterations


print(new_mat_sample)
