#############################################################################################################
# This script will be predict age of sample collection as a function of each genes methylation score and save the residuals
# to be used as predictors for age of diagnosis. 

##################################################################################################
# Load libraries
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

full_data_probe <- read.csv(paste0(data_folder, '/full_data_probe.csv'), stringsAsFactors = F)
full_data_probe$X <- NULL

load('/home/benbrew/Desktop/probe_data.RData')

getResidual <- function(data) {
  
  # get genes 
  genes <- colnames(data)[4:ncol(data)]

  # get matrix of sample collection and methylation
  data <- data[complete.cases(data),]
  
  # data <-  apply(data, 2, function(x) round(x, 3))
  # data <- as.data.frame(data)

  # predict each gene's methylation as a function of age of sample collection and save residuals - this will get the 
  # part of age of sample collection that is not predictive of gene methylation. your "instrument" is left hand side. 
  # some of the variance in the gene is predicted by age of sample collection, we want the variance that is independent 
  # of age of sample collection. if we did the age of sample collection, then we would have the variance of age of sample
  # collection that is indepent of gene, but we want to keep the gene component. So below will give us the variance of the 
  # gene that is not predicted (independent) of the age of sample collection and use that to predict age of diagnosis. 
  resid <- list()
  
  for (i in 4:ncol(data)){
    
    temp <- data[, i]
    temp1 <- data$age_sample_collection
    
    resid[[i]] <- lm(temp ~ temp1)$residuals

    print(i)
    
  }
  
  # general association between age and probe - some threshold 
  # treating age as a batch effect (combat)
  
  methyl_resid <- do.call('cbind', resid)
  methyl_resid <- cbind(data$age_diagnosis, data$age_sample_collection, methyl_resid)
  colnames(methyl_resid) <- c('age_diagnosis', 'age_sample_collection',genes)
  
  return(methyl_resid)
  
}

resid_full <- getResidual(full_data)
resid_full_probe <- getResidual(full_data_probe)


write.csv(resid_full, paste0(data_folder, '/resid_full.csv'))
write.csv(resid_full_probe, paste0(data_folder, '/resid_full_probe.csv'))






