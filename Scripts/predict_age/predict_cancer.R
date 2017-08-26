#######
# this script will use p53 cases and controls (full) to predict cancer or not to get an 
# idea as to what extent methylation has a strong cancer signature.

#!/hpf/tools/centos6/R/3.2.3/bin/Rscript

# Script to evaluate the how each imputation method affects the
# performance of the clustering methods

# argv <- as.numeric(commandArgs(T))

##########
# This script will get cleaned data from data saved in clean_data.R
# this will be the modeling pipeline script where we select features on 
# our training set and fit model. On the test set we test our model 
# and compare our predictions to values. 

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
method = 'raw'
k = 4
# seed_num <- 1
# seed_num <- argv[1]


##########
# load data
##########
betaCases <- readRDS(paste0(model_data, '/betaCases', method,'.rda'))
betaControls <- readRDS(paste0(model_data, '/betaControls', method,'.rda'))
betaControlsOld <- readRDS(paste0(model_data, '/betaControlsOld', method,'.rda'))

# load features
bh_feat_all <- readRDS(paste0(model_data, '/bh_feat_all.rda'))

# 
# # # TEMP
# betaCases <- betaCases[!grepl('9721365183', betaCases$sentrix_id),]

# load cg_locations
cg_locations <- read.csv(paste0(model_data, 
                                '/cg_locations.csv'))

cg_locations$X <- NULL
##########
# get intersecting colnames and prepare data for modeling
##########

intersect_names <- Reduce(intersect, list(colnames(betaCases)[8:ncol(betaCases)], 
                                          colnames(betaControls)[8:ncol(betaControls)]))
# organize each data set accordling

# cases
betaCases <- betaCases[, c('age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           intersect_names)]
# controls
betaControls <- betaControls[, c('age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender', 
                                 intersect_names)]

# controls
betaControlsOld <- betaControlsOld[, c('age_diagnosis', 
                                       'age_sample_collection', 
                                       'cancer_diagnosis_diagnoses', 
                                       'gender', 
                                       intersect_names)]

# #TEMP
betaFull <- rbind(betaCases, betaControls, betaControlsOld)

rm(betaCases, betaControls, betaControlsOld)

# create outcome variable cancer or not 
betaFull$cancer_diagnosis_diagnoses <- ifelse(grepl('Unaffected', 
                                                    betaFull$cancer_diagnosis_diagnoses), 'no', 'yes')


# get gender dummy variable
betaCases <- cbind(as.data.frame(class.ind(betaCases$gender)), betaCases)
betaControls <- cbind(as.data.frame(class.ind(betaControls$gender)), betaControls)


###########################################################################
# Next part of the pipline selects regions of the genome that are most differentially methylated 
# between 2 groups

# get a column for each dataset indicating the fold
seed_num <- 1
  
betaFull <- getFolds(betaFull, seed_number = seed_num, k_num = k)

bh_feat_sig <- getRun(bh_feat_all[[2]], run_num = .20)


trainTest <- function(cases, 
                      k) 
{
  
  # list to store results
  bh_feat <- list()
  model_results <- list()
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, cases$folds)
    test_index <- !train_index
    
    mod_result <- predCancer(training_dat = cases[train_index,], 
                             test_dat = cases[test_index,], 
                             bh_features = bh_feat_sig,
                             gender = T)
    
    
    model_results[[i]] <- getResultsCancer(mod_result)
    
    
  }
  
  return(model_results)
  
}
set.seed(4)
mod_results <- trainTest(cases = betaFull,
                         k = 4)


getClassResutls <- function(temp.result) {
  
  # creat matrix
  mat <- matrix(NA, nrow = 2, ncol = 2)

  for (j in 1:4) {
    
    if(j > 1) {
      
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
