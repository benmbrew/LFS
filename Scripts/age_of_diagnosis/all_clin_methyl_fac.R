####################################################################################################3
# This script will predict a multinomal age model ibrary(dplyr)
library(doParallel)
library(randomForest)
library(caret)


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
# remove extra columns
full_data$X <- NULL
full_data_cor$X <- NULL
full_data_rf$X <- NULL

# Create binary variables for age of diagnosis and age of sample collection on both data sets
full_data$age_diagnosis_fac <- ifelse(full_data$age_diagnosis <= 50, 1, 2)

full_data$age_diagnosis_fac_samp <- ifelse(full_data$age_sample_collection <= 50, 1, 2)

full_data_rf$age_diagnosis_fac <- ifelse(full_data_rf$age_diagnosis <= 50, 1, 2)

full_data_rf$age_diagnosis_fac_samp <- ifelse(full_data_rf$age_sample_collection <= 50, 1, 2)

# Create multinomal variables for age of diagnosis and age of sample collection
full_data$age_diagnosis_multi <- ifelse(full_data$age_diagnosis <= 30, 1, 
                                        ifelse(full_data$age_diagnosis > 30 & full_data$age_diagnosis <= 50, 2,
                                               ifelse(full_data$age_diagnosis > 50 & full_data$age_diagnosis <= 300, 3, 4)))

full_data$age_diagnosis_multi_samp <- ifelse(full_data$age_sample_collection <= 30, 1, 
                                        ifelse(full_data$age_sample_collection > 30 & full_data$age_sample_collection <= 50, 2,
                                               ifelse(full_data$age_sample_collection > 50 & full_data$age_sample_collection <= 300, 3, 4)))


full_data_rf$age_diagnosis_multi <- ifelse(full_data_rf$age_diagnosis <= 30, 1, 
                                        ifelse(full_data_rf$age_diagnosis > 30 & full_data_rf$age_diagnosis <= 50, 2,
                                               ifelse(full_data_rf$age_diagnosis > 50 & full_data_rf$age_diagnosis <= 300, 3, 4)))

full_data_rf$age_diagnosis_multi_samp <- ifelse(full_data_rf$age_sample_collection <= 30, 1, 
                                        ifelse(full_data_rf$age_sample_collection > 30 & full_data_rf$age_sample_collection <= 50, 2,
                                               ifelse(full_data_rf$age_sample_collection > 50 & full_data_rf$age_sample_collection <= 300, 3, 4)))

summary(as.factor(full_data$age_diagnosis_multi))
summary(as.factor(full_data_rf$age_diagnosis_multi))


# Create more granular multinomial variables for age of diagnosis and age of sample collection
full_data$age_diagnosis_multi_more <- ifelse(full_data$age_diagnosis <= 18, 1, 
                                        ifelse(full_data$age_diagnosis > 18 & full_data$age_diagnosis <= 36, 2,
                                               ifelse(full_data$age_diagnosis > 36 & full_data$age_diagnosis <= 288, 3,
                                                      ifelse(full_data$age_diagnosis > 288 & full_data$age_diagnosis <= 500, 4, 5))))

full_data$age_diagnosis_multi_more_samp <- ifelse(full_data$age_sample_collection <= 18, 1, 
                                                  ifelse(full_data$age_sample_collection > 18 & full_data$age_sample_collection <= 36, 2,
                                                         ifelse(full_data$age_sample_collection > 36 & full_data$age_sample_collection <= 288, 3,
                                                                ifelse(full_data$age_sample_collection > 288 & full_data$age_sample_collection <= 500, 4, 5))))

full_data_rf$age_diagnosis_multi_more <- ifelse(full_data_rf$age_diagnosis <= 18, 1, 
                                             ifelse(full_data_rf$age_diagnosis > 18 & full_data_rf$age_diagnosis <= 36, 2,
                                                    ifelse(full_data_rf$age_diagnosis > 36 & full_data_rf$age_diagnosis <= 288, 3,
                                                           ifelse(full_data_rf$age_diagnosis > 288 & full_data_rf$age_diagnosis <= 500, 4, 5))))

full_data_rf$age_diagnosis_multi_more_samp <- ifelse(full_data_rf$age_sample_collection <= 18, 1, 
                                                  ifelse(full_data_rf$age_sample_collection > 18 & full_data_rf$age_sample_collection <= 36, 2,
                                                         ifelse(full_data_rf$age_sample_collection > 36 & full_data_rf$age_sample_collection <= 288, 3,
                                                                ifelse(full_data_rf$age_sample_collection > 288 & full_data_rf$age_sample_collection <= 500, 4, 5))))



summary(as.factor(full_data$age_diagnosis_multi_more))
summary(as.factor(full_data_rf$age_diagnosis_multi_more))

# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       subset, 
                       selected_features,
                       binary,
                       iterations) {
  
  model <- list()
  predictions <- list()
  test.ground_truth <- list()
  rf.test_acc <- list()
  rf.test_auc <- list()

  genes <- colnames(data)[27:ncol(data)]
  
  data <- data[, c(subset, genes)]
  
  # Try the model with all different selection of features based on number of missinginess. 
  data <- data[complete.cases(data),]
  
  obs <- nrow(data)
  
  for (i in 1:iterations){
    
    set.seed(i)
    temp_levels <- 1
    loop_count <- i
    
    variable <- setdiff(subset, selected_features)
    
    if (binary){
      level_count <- 2
    } else if (multi) {
      level_count <- 4
    } else {
      level_count <- 5
    }
    
    while (length(temp_levels) < level_count) {
      set.seed(loop_count)
      train_index <- sample(nrow(data), nrow(data) *.7, replace = F)
      temp_levels <- levels(as.factor(data[, variable][-train_index]))
      loop_count <- i + 1
    }
    
    rf_y = make.names(as.factor(data[, variable])[train_index])

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
    
    model[[i]] <- train(x = data[train_index, c(selected_features, genes)]
                        , y = rf_y
                        , method = "rf"
                        , trControl = fitControl                   
                        , verbose = FALSE
                        , metric = "logLoss")
    
    predictions[[i]] <- predict(model[[i]] 
                                , newdata = data[-train_index, c(selected_features, genes)]
                                , type = "prob")
    
    test.ground_truth[[i]] <- as.factor(data[, variable][-train_index])
    test.ground_truth_1inN <- as.factor(class.ind(data[, variable][-train_index]))
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
    rf.test_acc[[i]] <- sum(levels(test.ground_truth[[i]])[apply(predictions[[i]], 1, which.is.max)] == test.ground_truth[[i]]) / dim(predictions[[i]])[1]
    # print(paste("RF test acc:", rf.test_acc))  
    # Compute Confusion Matrix and Statistics
    #confusionMatrix(pred, truth)
    # rf.test_stats[[i]] <- confusionMatrix(levels(test.ground_truth[[i]])[apply(predictions[[i]], 1, which.is.max)], test.ground_truth[[i]])
    # print(rf.test_stats)
    
    print(i)
    
  }
  
  return(list(predictions, test.ground_truth, rf.test_acc, model, rf.test_auc, obs))
  
}


##################################################################################################
# MDM2.NG
##########################
# Use Methylation
rf_methyl_fac <- predictAll(data = full_data_rf,
                        subset = c('age_diagnosis_fac'),
                        selected_features = NULL,
                        binary = TRUE,
                        iterations = 20)

rf_methyl_fac_samp <- predictAll(data = full_data_rf,
                            subset = c('age_diagnosis_fac_samp'),
                            selected_features = NULL,
                            binary = TRUE,
                            iterations = 20)

rf_methyl_multi <- predictAll(data = full_data_rf,
                            subset = c('age_diagnosis_multi'),
                            selected_features = NULL,
                            binary = FALSE,
                            iterations = 20)

rf_methyl_multi_samp <- predictAll(data = full_data_rf,
                              subset = c('age_diagnosis_multi_samp'),
                              selected_features = NULL,
                              binary = FALSE,
                              iterations = 20)

###########################
# unlist and combine 

all <- rbind (
  append('binary_50', mean(unlist(rf_methyl_fac[[3]]))),
  append('binary_50_sample', mean(unlist(rf_methyl_fac_samp[[3]]))) ,
  append('multi_30_50_300', mean(unlist(rf_methyl_multi[[3]]))),
  append('multi_50_50_300_sample', mean(unlist(rf_methyl_multi[[3]])))

)


# rf_methyl_multi_more <- predictAll(data = full_data_rf,
#                               subset = c('age_diagnosis_multi_more'),
#                               selected_features = NULL,
#                               binary = FALSE,
#                               multi = FALSE,
#                               iterations = 20)
# 
# rf_methyl_multi_samp <- predictAll(data = full_data_rf,
#                                    subset = c('age_diagnosis_multi_more_samp'),
#                                    selected_features = NULL,
#                                    binary = FALSE,
#                                    iterations = 20)

# # unlist accuracy results and put into data frame 
# 
# rf_mdm2.nG_clin_reg <- predictAll(data = full_data_rf,
#                                   use_genes = FALSE,
#                                   subset = c('age_diagnosis_fac', 'codon72.npro', 'gdna.exon.intron', 'gdna.base.change', 'gdna.codon',
#                                              'protein.codon.change', 'protein.codon.num', 'splice.delins.snv',
#                                              'mdm2.nG'),
#                                   selected_features = c('codon72.npro', 'gdna.exon.intron', 'gdna.base.change', 'gdna.codon',
#                                                         'protein.codon.change', 'protein.codon.num', 'splice.delins.snv',
#                                                         'mdm2.nG'),
#                                   iterations = 5)
