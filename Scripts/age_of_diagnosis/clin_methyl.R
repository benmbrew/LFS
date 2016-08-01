### 
# this script will run the models

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/regression_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


# Load Libraries and Dependencies ---------------------------------
library("igraph")
library("foreach")
library("ROCR")
library("doParallel")
library("caret")
library("glmnet")
library("randomForest")
library("kernlab")
library("pROC")
library(dplyr)

# Load in clinical data
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# subset to those variables 
# subset <- c("family_name", "relationship", "age_diagnosis", "p53_germline","gdna", 
#             "protein", "codon72", "mdm2", "gender", 'methyl_indicator')

subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
            "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
            "mdm2.nG")

clin <- clin[, subset]

# Try the model with all different selection of features based on number of missinginess. 
clin <- clin[complete.cases(clin),]

# convert characters to factors 
# convert characters to factors 
for ( i in 1:ncol(clin)){
  
  if(!grepl('age', colnames(clin[i]))) {
    clin[,i] <- as.factor(clin[,i])
    
  } 
}
##################################################################
# Random Forest - this is training on no methylation and test on methylation

train_index <- clin$methyl_indicator == 'No'
test_index <- clin$methyl_indicator == 'Yes'

rf_y = clin$age_diagnosis[train_index]

# determines how you train the model.
NFOLDS <- 5
fitControl <- trainControl( 
  method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
  number = min(10, NFOLDS),      
  repeats = 1,
  allowParallel = TRUE
  #summaryFunction = summaryFunc
)

selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                      "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                      "mdm2.nG")

rf.model <- train(x = clin[train_index, selected_features]
                  , y = rf_y
                  , method = "rf"
                  , trControl = fitControl                   
                  , verbose = FALSE)


rf.predictions <- predict(rf.model 
                          , newdata = clin[test_index, selected_features])

test.ground_truth <- clin$age_diagnosis[test_index]
test.rf_mse_no_methyl <- rmse(rf.predictions, test.ground_truth)

# 34.8 months - lower than all other rmse 

# ##################################################################
# # Random Forest - this is training on no methylation - split into 70% and predict on the 30% 
# # then it predicts on yes methylation
# # Set training set to no methylation
# train <- clin[clin$methyl_indicator == 'No',]
# train_index <- sample(nrow(train), nrow(train) *.7, replace = F)
# 
# for (i in clin$methyl_indicator) {
#   
#   test <- clin[clin$methyl_indicator == i,]  
#   
#   rf_y = train$age_diagnosis[train_index]
#   
#   # determines how you train the model.
#   NFOLDS <- 5
#   fitControl <- trainControl( 
#     method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
#     number = min(10, NFOLDS),      
#     repeats = 1,
#     allowParallel = TRUE
#     #summaryFunction = summaryFunc
#   )
#   
#   selected_features <- c("p53_germline","gdna", 
#                          "protein", "codon72", "mdm2", "gender")
#   
#   rf.model <- train(x = train[train_index, selected_features]
#                     , y = rf_y
#                     , method = "rf"
#                     , trControl = fitControl                   
#                     , verbose = FALSE)
#   
#   if (i == 'No') {
#     
#     rf.predictions <- predict(rf.model 
#                               , newdata = train[-train_index, selected_features])
#     
#     test.ground_truth <- train$age_diagnosis[-train_index]
#     test.rf_mse_no_methyl <- mean((rf.predictions - test.ground_truth))^2
#     
#   } else {
#     
#     rf.predictions <- predict(rf.model 
#                               , newdata = test[1:nrow(test), selected_features])
#     
#     test.ground_truth <- test$age_diagnosis
#     test.rf_mse_methyl <- mean((rf.predictions - test.ground_truth))^2
#   }
#   
# }
# 
