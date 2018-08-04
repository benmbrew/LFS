
#!/hpf/tools/centos6/R/3.2.3/bin/Rscript

# Script to evaluate the how each imputation method affects the
# performance of the clustering methods

# argv <- as.numeric(commandArgs(T))

##########
# This script will predict cancer with cases and controls, both 450k and 850k

# predict cancer = methylation (p53 mutants)

##########
# initialize libraries
##########
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(doParallel)
library(nnet)
library(dplyr)
library(bumphunter)
library(sqldf)
library(e1071)

registerDoParallel(1)
##########
# initialize folders
##########
home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'noob'
k = 5
combat = F

##########
# load data
##########

if (combat) {
  
  betaFull <-  readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_combat.rda')))
  
  
} else {
  
  betaFull <- readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data.rda')))
  
}



# summary of data types 
summary(as.factor(betaFull$batch))

# remove duplicates 
betaFull <- betaFull[!duplicated(betaFull$ids),]

# cancer, no cancer
summary(as.factor(betaFull$cancer_diagnosis_diagnoses))

# get folds 
betaFull <- getFolds(betaFull, seed_number = 1, k = 5)
summary(as.factor(betaFull$folds))

beta_dat <- betaFull
# lets have this model either do 450 only, 850 only, or both, no combinations
trainTest <- function(beta_dat,
                      methyl_tech,
                      k) 
{
  
  
  # if(methyl_tech == '450k'){
  #   
  #   beta_dat <- dat[grepl('450k', dat$type),]
  #   
  # } else if (methyl_tech =='850k') {
  #   
  #   beta_dat <- dat[grepl('850k', dat$type),]
  #   
  # } else if(methyl_tech == '450_850') {
  #   
  #   beta_dat <- subset(dat, type == '450k' |  (type == '850k' & cancer_diagnosis_diagnoses == 'Unaffected'))
  #   
  # } else if (methyl_tech =='850_450') {
  #   
  #   beta_dat <- subset(dat, type == '850k' |  (type == '450k' & cancer_diagnosis_diagnoses == 'Unaffected'))
  #   
  # } else {
  #   
  #   beta_dat <- dat
  #   
  # }
  # 
  
  # list to store results
  bh_feat <- list()
  model_results <- list()
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, beta_dat$folds)
    test_index <- !train_index
    
    mod_feats <- colnames(beta_dat)[9:(ncol(beta_dat) -1)]
    
    mod_feats <- sample(mod_feats, 10000, replace = T)
    
    mod_result <- predCancer(training_dat = beta_dat[train_index,], 
                             test_dat = beta_dat[test_index,], 
                             bh_features = mod_feats)
    
    
    model_results[[i]] <- getResultsCancer(mod_result)
    
    
  }
  
  return(model_results)
  
}
set.seed(4)
mod_results <- trainTest(beta_dat = betaFull, methyl_tech = 'all', k = 5)


getClassResutls <- function(temp.result) {
  
  # creat matrix
  mat <- matrix(NA, nrow = 2, ncol = 2)

  for (j in 1:4) {
    
    if (j > 1) {
      
      mat <- mat + temp.result[[j]][[4]]$table[1:2,1:2]
      acc_temp <- acc_temp + temp.result[[j]][[4]]$overall[1]
      tpr_temp <- tpr_temp + temp.result[[j]][[4]]$byClass[1]
      tnr_temp <- tnr_temp + temp.result[[j]][[4]]$byClass[2]
      alpha_temp <- alpha_temp + temp.result[[j]][[1]]
      lambda_temp <- lambda_temp + temp.result[[j]][[1]]
      
      
    } else{
      
      mat <- temp.result[[j]][[4]]$table[1:2,1:2]
      acc_temp <- as.numeric(temp.result[[j]][[4]]$overall[1])
      tpr_temp <- as.numeric(temp.result[[j]][[4]]$byClass[1])
      tnr_temp <- as.numeric(temp.result[[j]][[4]]$byClass[2])
      alpha_temp <- as.numeric(temp.result[[j]][[1]])
      lambda_temp <- as.numeric(temp.result[[j]][[1]])
      
    }
    
  }
  
  # get average from j
  confusion_table <- mat/4
  acc_table <- acc_temp/4
  tnr_table <- tnr_temp/4
  tpr_table <- tpr_temp/4
  alpha_table <- alpha_temp/4
  lambda_table <- lambda_temp/4
  
  # store at i 
  table_results <- confusion_table
  
  table_results <- as.data.frame(table_results)
  
  table_scores <- as.data.frame(cbind(alpha_table, lambda_table, acc_table, tnr_table, tpr_table))
  
  # remove NAs - bad predictions dont yield a value for precision, etc
  table_scores <- table_scores[complete.cases(table_scores),]
  
  # fpr, fnr
  table_scores$fpr <- 1 - table_scores$tnr
  table_scores$fnr <- 1 - table_scores$tpr
  

  return(list(table_results, table_scores))
}

# get results for age 48
result_48 <- getClassResutls(mod_results)
conMat_48 <- result_48[[1]]
resultTab_48 <- result_48[[2]]
