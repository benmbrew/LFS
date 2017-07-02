### This script will analyze results from training and testing on cases 
# specificity - TNR
# sensitivity - TPR
# FPR (1 -TNR)
# FNR (1- TPR)
##########
# initiate library
##########
library(tidyverse)

##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
results_folder <- paste0(project_folder, '/Scripts/predict_age/Results')
class_folder <- paste0(results_folder, '/class_results')


###########
# load in reg results 
###########
iter_length <- 50


getClassResutls <- function(age) {
  
  # creat matrix
  mat <- matrix(NA, nrow = 2, ncol = 2)
  table_results <- list()
  alpha_scores <- list()
  lambda_num <- list()
  acc <- list()
  tnr <- list()
  tpr <- list()
  
  for (i in 1:iter_length) {
    # read in raw results
    temp.result <- readRDS(paste0(class_folder, '/train_test_', i, '_', age,  '.rda' ))
    
    # remove importance and model
    temp.result[[3]] <- temp.result[[5]] <- NULL
    
    for (j in 1:4) {

      if(j > 1) {
        
        mat <- mat + temp.result[[3]][[j]]$table[1:2,1:2]
        acc_temp <- acc_temp + temp.result[[3]][[j]]$overall[1]
        tpr_temp <- tpr_temp + temp.result[[3]][[j]]$byClass[1]
        tnr_temp <- tnr_temp + temp.result[[3]][[j]]$byClass[2]
      } else{
        
        mat <- temp.result[[3]][[j]]$table[1:2,1:2]
        acc_temp <- as.numeric(temp.result[[3]][[j]]$overall[1])
        tpr_temp <- as.numeric(temp.result[[3]][[j]]$byClass[1])
        tnr_temp <- as.numeric(temp.result[[3]][[j]]$byClass[2])

      }
      
    }
    
    # get average from j
    confusion_table <- mat/4
    acc_table <- acc_temp/4
    tnr_table <- tnr_temp/4
    tpr_table <- tpr_temp/4
    
    # store at i 
    table_results[[i]] <- confusion_table
    alpha_scores[[i]] <- unlist(temp.result[[1]])
    lambda_num[[i]] <- unlist(temp.result[[2]])
    acc[[i]] <- acc_table
    tnr[[i]] <- tnr_table
    tpr[[i]] <- tpr_table
    
    
  }
  
  table_results <- Reduce("+",table_results)/50
  table_results <- as.data.frame(table_results)
  
  # get scores by alpha parameter
  alpha_scores <- unlist(alpha_scores)
  lambda_num <- unlist(lambda_num)
  acc <- unlist(acc)
  tnr <- unlist(tnr)
  tpr <- unlist(tpr)
  
  table_scores <- as.data.frame(cbind(alpha_scores, lambda_num, acc, tnr, tpr))
  
  # remove NAs - bad predictions dont yield a value for precision, etc
  table_scores <- table_scores[complete.cases(table_scores),]
  
  # fpr, fnr
  table_scores$fpr <- 1 - table_scores$tnr
  table_scores$fnr <- 1 - table_scores$tpr
  
  # get mean
  # group by alpha get means
  mean_result <- table_scores %>%
    group_by(alpha_scores) %>%
    summarise(mean_acc = mean(acc),
              mean_tpr = mean(tpr),
              mean_tnr = mean(tnr),
              mean_fpr = mean(fpr),
              mean_fnr = mean(fnr))
  
  # get indicator for age 
  mean_result$age <- age
  
  
  return(list(table_results, mean_result))
}

# get results for age 48
result_48 <- getClassResutls(48)
conMat_48 <- result_48[[1]]
resultTab_48 <- result_48[[2]]

# # get results for age 60
# result_60 <- getClassResutls(60)
# conMat_60 <- result_60[[1]]
# resultTab_60 <- result_60[[2]]

# get results for age 72
result_72 <- getClassResutls(72)
conMat_72 <- result_72[[1]]
resultTab_72 <- result_72[[2]]

# combine all ages 
results <- rbind(resultTab_48,
                 resultTab_72)

