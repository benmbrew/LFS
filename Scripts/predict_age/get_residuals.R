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
full_data_cor <- read.csv(paste0(data_folder, '/full_data_cor.csv'), stringsAsFactors = F)
full_data_rf <- read.csv(paste0(data_folder, '/full_data_rf.csv'), stringsAsFactors = F)
# load dmr methylatio data
methyl_gene_dmr <- read.csv(paste0(data_folder, '/methyl_gene_dmr.csv'))
methyl_dmr <- read.csv(paste0(data_folder, '/methyl_dmr.csv'))

full_data$X <- NULL
full_data_cor$X <- NULL
full_data_rf$X <- NULL
methyl_gene_dmr$X <- NULL
methyl_dmr$X <- NULL



getResidual <- function(data) {
  
  # get genes 
  genes <- colnames(data)[30:ncol(data)]

  # get matrix of sample collection and methylation
  data <- data[, c(6, 8, 30:ncol(data))]
  data <- data[complete.cases(data),]
  
  # predict each gene's methylation as a function of age of sample collection and save residuals - this will get the 
  # part of age of sample collection that is not predictive of gene methylation. your "instrument" is left hand side. 
  # some of the variance in the gene is predicted by age of sample collection, we want the variance that is independent 
  # of age of sample collection. if we did the age of sample collection, then we would have the variance of age of sample
  # collection that is indepent of gene, but we want to keep the gene component. So below will give us the variance of the 
  # gene that is not predicted (independent) of the age of sample collection and use that to predict age of diagnosis. 
  resid <- list()
  
  for (i in 3:ncol(data)){
    
    resid[[i]] <- lm(data[, i] ~ data$age_sample_collection, data = data)$residuals

    print(i)
    
  }
  
  methyl_resid <- do.call('cbind', resid)
  methyl_resid <- cbind(data$age_diagnosis, data$age_sample_collection, methyl_resid)
  colnames(methyl_resid) <- c('age_diagnosis', 'age_sample_collection',genes)
  
  return(methyl_resid)
  
}

resid_rf <- getResidual(full_data_rf)
resid_cor <- getResidual(full_data_cor)
resid_full <- getResidual(full_data)
resid_dmr <- getResidual(methyl_dmr)
resid_gene_dmr <- getResidual(methyl_gene_dmr)

write.csv(resid_rf, paste0(data_folder, '/resid_rf.csv'))
write.csv(resid_cor, paste0(data_folder, '/resid_cor.csv'))
write.csv(resid_full, paste0(data_folder, '/resid_full.csv'))
write.csv(resid_dmr, paste0(data_folder, '/resid_dmr.csv'))
write.csv(resid_gene_dmr, paste0(data_folder, '/resid_gene_dmr.csv'))




