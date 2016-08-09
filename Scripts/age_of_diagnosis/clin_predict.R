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
# Load in clinical data
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)


full_data$X <- NULL
full_data_cor$X <- NULL
full_data_rf$X <- NULL


# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       subset, 
                       selected_features,
                       use_genes,
                       iterations) {
  
  model <- list()
  predictions <- list()
  rf.test_stats <- list()
  rf.test_acc <- list()
  predict_test <- list()
  test.ground_truth <- list()
  rf_rest_results <- list()


  if (use_genes) {
    genes <- colnames(data)[27:ncol(data)]
    data <- data[, c(subset, genes)]
    data[, subset] <- as.factor(data[,subset])
  }else{
    genes <- NULL
    data <- data[, subset]
    # convert characters to factors 
    for ( i in 1:ncol(data)){
      
      if(typeof(data[,i]) == 'character' || typeof(data[,i]) == 'integer') {
        data[,i] <- as.factor(data[,i])
        print(i)
      } 
    }
  }
  
  # Try the model with all different selection of features based on number of missinginess. 
  test_index <- is.na(data$codon72.npro)
  test_data <- data
  
  data <- data[complete.cases(data),]
  
  obs <- nrow(data)
  
  for (i in 1:iterations){
    
    set.seed(i)
    train_index <- sample(nrow(data), nrow(data) *.7, replace = F)
    
    # 4) Random Forest 
    variable <- setdiff(subset, selected_features)
    rf_y = make.names(as.factor(data[, variable])[train_index])

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
    
    model[[i]] <- train(x = data[train_index, c(selected_features, genes)]
                      , y = rf_y
                      , method = "rf"
                      , trControl = fitControl                   
                      , verbose = FALSE
                      , metric = "logLoss")
    
    predictions[[i]] <- predict(model[[i]] 
                              , newdata = data[-train_index, c(selected_features, genes)]
                              , type = "prob")
    
    predict_test[[i]] <- predict(model[[i]] 
                                , newdata = test_data[test_index, c(selected_features, genes)]
                                , type = "prob")
    
    test.ground_truth[[i]] <- data[, variable][-train_index]
    test.ground_truth_1inN <- class.ind(data[, variable][-train_index])
    #print(table(levels(test.ground_truth)[apply(rf.predictions, 1, which.is.max)], test.ground_truth))
    # AUC
    #create ROCR prediction object
    temp.predict <- prediction(predictions[[i]], test.ground_truth_1inN)  
    # rf.test_auc <- unlist(slot(performance(temp.predict, "auc"), "y.values"))  
    # print(paste("RF test AUC:", rf.test_auc))  
    # Accuracy
    rf.test_acc[[i]] <- sum(levels(test.ground_truth[[i]])[apply(predictions[[i]], 1, which.is.max)] == test.ground_truth[[i]]) / dim(predictions[[i]])[1]
    rf_rest_results[[i]] <- levels(test.ground_truth[[i]])[apply(predict_test[[i]], 1, which.is.max)]
    # print(paste("RF test acc:", rf.test_acc))  
    # Compute Confusion Matrix and Statistics
    #confusionMatrix(pred, truth)
    # rf.test_stats[[i]] <- confusionMatrix(levels(test.ground_truth[[i]])[apply(predictions[[i]], 1, which.is.max)], test.ground_truth[[i]])
    # print(rf.test_stats)
    
    print(i)
    
  }
  
  return(list(rf.test_stats, predictions, model, rf.test_acc, test.ground_truth, obs))
  
}


##################################################################################################
# MDM2.NG
##########################
# Use Methylation
# predict mdm2.nG with subset 
rf_mdm2.nG_reg <- predictAll(data = full_data_rf,
                             use_genes = TRUE,
                             subset = "mdm2.nG",
                             selected_features = NULL, iterations = 5)

# get avg accuracy 
mean(unlist(rf_mdm2.nG_reg[[4]]))


##########################
# predict mdm2.nG with correlation data 
rf_mdm2.nG_cor <- predictAll(data = full_data_cor,
                            use_genes = TRUE,
                            subset = "mdm2.nG",
                            selected_features = NULL, iterations = 5)

# get avg accuracy 
mean(unlist(rf_mdm2.nG_cor[[4]]))


##########################
# predict mdm2.nG with full data 
rf_mdm2.nG_full <- predictAll(data = full_data_cor,
                              use_genes = TRUE,
                              subset = "mdm2.nG",
                              selected_features = NULL, iterations = 5)

# get avg accuracy 
mean(unlist(rf_mdm2.nG_full[[4]]))

######################################################################################################
# codon72.npro
##########################
# Use Methylation
# predict codon72.npro with subset 
rf_codon72.npro_reg <- predictAll(data = full_data_rf,
                             use_genes = TRUE,
                             subset = "codon72.npro",
                             selected_features = NULL, iterations = 5)

# get avg accuracy 
mean(unlist(rf_codon72.npro_reg[[4]]))


##########################
# predict codon72.npro with correlation data 
rf_codon72.npro_cor <- predictAll(data = full_data_cor,
                             use_genes = TRUE,
                             subset = "codon72.npro",
                             selected_features = NULL, iterations = 5)

# get avg accuracy 
mean(unlist(rf_codon72.npro_cor[[4]]))


##########################
# predict codon72.npro with full data 
rf_codon72.npro_full <- predictAll(data = full_data_cor,
                              use_genes = TRUE,
                              subset = "codon72.npro",
                              selected_features = NULL, iterations = 5)

# get avg accuracy 
mean(unlist(rf_codon72.npro_full[[4]]))


#######################################################################################
# Use Clinical
###########################
# predict codon72.npro with clinical data
rf_codon72.npro_clin_full <- predictAll(data = clin,
                                   subset = c('codon72.npro', 'gdna.exon.intron', 'gdna.base.change', 'gdna.codon',
                                              'protein.codon.change', 'protein.codon.num', 'splice.delins.snv',
                                              'mdm2.nG'),
                                   selected_features = c('gdna.exon.intron', 'gdna.base.change', 'gdna.codon',
                                                    'protein.codon.change', 'protein.codon.num', 'splice.delins.snv',
                                                    'mdm2.nG'),
                                   use_genes = FALSE,
                                   iterations = 5)

# get avg accuracy 
mean(unlist(rf_codon72.npro_clin_full[[4]]))



######################################################################################################
# splice.delins.snv
##########################
# Use Methylation
# predict splice.delins.snv with subset 
rf_splice.delins.snv_reg <- predictAll(data = full_data_rf,
                                  use_genes = TRUE,
                                  subset = "splice.delins.snv",
                                  selected_features = NULL, iterations = 5)

# get avg accuracy 
mean(unlist(rf_splice.delins.snv_reg[[4]]))


##########################
# predict splice.delins.snv with correlation data 
rf_splice.delins.snv_cor <- predictAll(data = full_data_cor,
                                  use_genes = TRUE,
                                  subset = "splice.delins.snv",
                                  selected_features = NULL, iterations = 5)

# get avg accuracy 
mean(unlist(rf_splice.delins.snv_cor[[4]]))


##########################
# predict splice.delins.snv with full data 
rf_splice.delins.snv_full <- predictAll(data = full_data_cor,
                                   use_genes = TRUE,
                                   subset = "splice.delins.snv",
                                   selected_features = NULL, iterations = 20)

# get avg accuracy 
mean(unlist(rf_splice.delins.snv_full[[4]]))


#######################################################################################
# Use Clinical
###########################
# predict splice.delins.snv with clinical data
rf_splice.delins.snv_clin_full <- predictAll(data = clin,
                                        subset = c('splice.delins.snv', 'gdna.exon.intron', 'gdna.base.change', 'gdna.codon',
                                                   'protein.codon.change', 'protein.codon.num', 'splice.delins.snv',
                                                   'mdm2.nG'),
                                        selected_features = c('gdna.exon.intron', 'gdna.base.change', 'gdna.codon',
                                                              'protein.codon.change', 'protein.codon.num', 'splice.delins.snv',
                                                              'splice.delins.snv'),
                                        use_genes = FALSE,
                                        iterations = 20)

# get avg accuracy 
mean(unlist(rf_splice.delins.snv_clin_full[[4]]))



######################################################################################################3
# predict onto 


