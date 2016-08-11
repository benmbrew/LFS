# this script will read in different versions of subsetted full data and run random forest. 
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
full_data_cor <- read.csv(paste0(data_folder, '/full_data_cor.csv'), stringsAsFactors = F)
full_data_rf <- read.csv(paste0(data_folder, '/full_data_rf.csv'), stringsAsFactors = F)

full_data_cor$X <- NULL
full_data_rf$X <- NULL

data <- full_data_cor

# get genes 
genes <- colnames(data)[27:ncol(data)]

## get matrix of sample collection and methylation
data <- data[, c(6, 8, 27:ncol(data))]
data <- data[complete.cases(data),]

# predict age of sample collection with each gene and store residuals
resid <- list()

for (i in 3:ncol(data)){
  
  resid[[i]] <- lm(data[, i] ~ data$age_sample_collection, data = data)$residuals
  print(i)
  
}

methyl_resid <- do.call('cbind', resid)
colnames(methyl_resid) <- genes


write.csv(methyl_resid, paste0(data_folder, '/methyl_resid_cor.csv'))
