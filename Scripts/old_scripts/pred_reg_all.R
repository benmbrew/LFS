##########
# initialize libraries
##########
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(ModelMetrics)
library(doParallel)
library(nnet)
library(dplyr)
library(bumphunter)
library(sqldf)
library(e1071)
library(reshape2)

registerDoParallel(1)

##########
# initialize folders
##########

home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
feat_data <- paste0(data_folder, '/feat_data')
results_data <- paste0(data_folder, '/results_data')

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
  
  betaFull <-  readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_combat_m.rda')))
  
  
} else {
  betaFull <- readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_m.rda')))
  
}


##########
# get column names
##########
intersect_names <- colnames(betaFull)[9:ncol(betaFull)]

##########
# read in all features from feat data
##########

lfs_feats_m <- readRDS(paste0(feat_data, paste0('/', method, '_', 'lfs_m.rda')))
no_cancer_feats_m <- readRDS(paste0(feat_data, paste0('/', method, '_', 'no_cancer_m.rda')))


setwd(feat_data)
file_list_names = list.files()


# store all raw rda feature lists in feat_list
feat_list <- lapply(file_list_names, function(x) readRDS(x))

# order feat list
feat_list  <- feat_list[order(sapply(feat_list, length), decreasing=F)]

# select first 10
feat_list <- feat_list[10]

# model_names <- c('enet', 'rf', 'lasso')
# seeds <- c(1, 2, 3)

model_names <- c('enet')
seeds <- c(1,2,3)
feat_list <- no_cancer_feats_m
file_list_names <- 'lfs_m'

betaFull <- betaFull[, 1:10000]

# 
l = j = m = 1
#
beta_full = betaFull
model_method = model_names[m]
mod_feats = feat_list
feat_name = file_list_names[[l]]
seed_num = seeds[j]
k = k
max_age = 1000
data_type = 'standard'



trainTest <- function(beta_full,
                      model_method,
                      mod_feats,
                      seed_num,
                      feat_name,
                      data_type,
                      max_age,
                      remove_cases,
                      k) {
  
  
  if(data_type == 'standard') {
    
    beta_cases <- beta_full[beta_full$batch == 'cases',]
    beta_controls <- beta_full[beta_full$batch == 'controls',]
    beta_valid <- beta_full[beta_full$batch == 'valid',]
    
  } else {
    
    beta_cases <- beta_cases[grepl('cases|valid', beta_cases$batch)]
    beta_controls <- beta_controls[grepl('controls', beta_cases$batch)]
    beta_valid <- beta_full[beta_full$batch == 'valid',]
    
  }
  
  # remove samples that dont have an age of sample collection
  beta_cases <- beta_cases[complete.cases(beta_cases),]
  beta_valid <- beta_valid[complete.cases(beta_valid),]
  
  beta_controls <- beta_controls[!is.na(beta_controls$age_sample_collection),]
  beta_cases <- beta_cases[!is.na(beta_cases$age_sample_collection),]
  
  # subet by max age
  beta_cases <- beta_cases[beta_cases$age_sample_collection <= max_age,]
  
  if(remove_cases) {
    length_index <- 1:nrow(beta_cases)
    remove_index <- length_index[beta_cases$age_sample_collection >= 12 & beta_cases$age_sample_collection <= 24]
    remove_index <- sample(remove_index, 8, replace = F)
    
    beta_cases <- beta_cases[-remove_index,]
    
  }
  
  set.seed(seed_num)
  # get a column for each dataset indicating the fold
  beta_cases <- getFolds(beta_cases, seed_number = seed_num, k = k)
  
  # list to store results
  temp_results <- list()
  
  # sample random 10 from each cluster 
  column_names <- colnames(beta_full)[9:ncol(beta_full)]
  sample_cols <- column_names[!column_names %in% mod_feats]
  rand_feats <- sample(sample_cols, length(mod_feats), replace = T)
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, beta_cases$folds)
    test_index <- !train_index
    
    if(model_method == 'enet') {
      
      mod_result <- runEnetRand(training_dat = beta_cases[train_index,], 
                                controls_dat = beta_controls,
                                valid_dat = beta_valid,
                                test_dat = beta_cases[test_index,], 
                                age_cutoff = class_age,
                                bh_features = mod_feats,
                                gender = T)
      
    } 
    
    if(model_method == 'rf') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      mod_result <- runRfRandFac(training_dat = beta_cases[train_index,], 
                                 controls_dat = beta_controls,
                                 valid_dat = beta_valid,
                                 test_dat = beta_cases[test_index,], 
                                 age_cutoff = class_age,
                                 bh_features = mod_feats,
                                 rand_feats = rand_feats,
                                 pred_cutoff = .5,
                                 gender = gender)
      
      
    } 
    # get results 
    temp_results[[i]] <- mod_result
    
  }
  
  return(temp_results)
  
}

seeds <- c(1,2, 3, 4, 5)
con_pred_age_temp <- list()
con_real_age_temp <- list()
con_pred_age <- list()
con_real_age <- list()
case_pred_age_temp <- list()
case_real_age_temp <- list()
case_pred_age <- list()
case_real_age <- list()
val_pred_age_temp <- list()
val_real_age_temp <- list()
val_pred_age <- list()
val_real_age <- list()



# set fixed variables 
    
for (j in 1:length(seeds)) {
  
  # return list of 2
  mod_results  <- trainTest(beta_full = betaFull,
                           model_method = 'enet',
                           mod_feats = feat_list,
                           seed_num = seeds[j],
                           data_type = 'standard',
                           max_age = 10000,
                           remove_cases = F,
                           k = 5)
  
  for (i in 1:k) {
    
    con_pred_age_temp[[i]]  <- mod_results[[i]][7]
    con_real_age_temp[[i]] <- mod_results[[i]][8]
    case_pred_age_temp[[i]]  <- mod_results[[i]][1]
    case_real_age_temp[[i]] <- mod_results[[i]][2]
    val_pred_age_temp[[i]]  <- mod_results[[i]][4]
    val_real_age_temp[[i]] <- mod_results[[i]][5]
    
  }
  
  con_pred_age[[j]] <- unlist(con_pred_age_temp)
  con_real_age[[j]] <- unlist(con_real_age_temp)
  case_pred_age[[j]] <- unlist(case_pred_age_temp)
  case_real_age[[j]] <- unlist(case_real_age_temp)
  val_pred_age[[j]] <- unlist(val_pred_age_temp)
  val_real_age[[j]] <- unlist(val_real_age_temp)
  


  print(paste0('done with ', j, ' seed'))
}

# comibine 
con_pred_age_final <- do.call(rbind, con_pred_age)
con_real_age_final <- do.call(rbind, con_real_age)
case_pred_age_final <- do.call(rbind, case_pred_age)
case_real_age_final <- do.call(rbind, case_real_age)
val_pred_age_final <- do.call(rbind, val_pred_age)
val_real_age_final <- do.call(rbind, val_real_age)

# get mean across seeds
con_pred_age_final_1 <- apply(con_pred_age_final, 2, function(x) mean(x))
con_real_age_final_1 <- apply(con_real_age_final, 2, function(x) mean(x))
case_pred_age_final_1 <- apply(case_pred_age_final, 2, function(x) mean(x))
# case_real_age_final_1 <- apply(case_real_age_final, 2, function(x) mean(x))
val_pred_age_final_1 <- apply(val_pred_age_final, 2, function(x) mean(x))
val_real_age_final_1 <- apply(val_real_age_final, 2, function(x) mean(x))

# get mean across folds - sets of 30 by 5
con_pred_age_final_1  <- as.data.frame(cbind(con_pred_age_final_1, c(seq(1, 30, 1), 
                                          seq(1, 30,1), 
                                          seq(1, 30,1), 
                                          seq(1, 30,1), 
                                          seq(1, 30,1))))

# get mean across folds - sets of 30 by 5
con_real_age_final_1  <- as.data.frame(cbind(con_real_age_final_1, c(seq(1, 30, 1), 
                                                             seq(1, 30,1), 
                                                             seq(1, 30,1), 
                                                             seq(1, 30,1), 
                                                             seq(1, 30,1))))
# # get mean across folds - sets of 30 by 5
# case_pred_age_final_1  <- as.data.frame(cbind(case_pred_age_final_1, c(seq(1, 30, 1), 
#                                                                      seq(1, 30,1), 
#                                                                      seq(1, 30,1), 
#                                                                      seq(1, 30,1), 
#                                                                      seq(1, 30,1))))
# 
# # get mean across folds - sets of 30 by 5
# case_real_age_final_1  <- as.data.frame(cbind(case_real_age_final_1, c(seq(1, 30, 1), 
#                                                                      seq(1, 30,1), 
#                                                                      seq(1, 30,1), 
#                                                                      seq(1, 30,1), 
#                                                                      seq(1, 30,1))))

# get mean across folds - sets of 30 by 5
val_pred_age_final_1  <- as.data.frame(cbind(val_pred_age_final_1, c(seq(1, 37, 1), 
                                                                     seq(1, 37,1), 
                                                                     seq(1, 37,1), 
                                                                     seq(1, 37,1), 
                                                                     seq(1, 37,1))))

# get mean across folds - sets of 30 by 5
val_real_age_final_1  <- as.data.frame(cbind(val_real_age_final_1, c(seq(1, 37, 1), 
                                                                     seq(1, 37,1), 
                                                                     seq(1, 37,1), 
                                                                     seq(1, 37,1), 
                                                                     seq(1, 37,1))))


# make factors
con_pred_age_final_1$V2 <- as.factor(con_pred_age_final_1$V2)
con_real_age_final_1$V2 <- as.factor(con_real_age_final_1$V2)

val_pred_age_final_1$V2 <- as.factor(val_pred_age_final_1$V2)
val_real_age_final_1$V2 <- as.factor(val_real_age_final_1$V2)


# group by fold and get mean
con_pred_final <- con_pred_age_final_1 %>%
  group_by(V2) %>%
  summarise(mean_x = mean(con_pred_age_final_1))

# group by fold and get mean
con_real_final <- con_real_age_final_1 %>%
  group_by(V2) %>%
  summarise(mean_x = mean(con_real_age_final_1))

# # group by fold and get mean
# case_pred_final <- case_pred_age_final_1 %>%
#   group_by(V2) %>%
#   summarise(mean_x = mean(case_pred_age_final_1))
# 
# # group by fold and get mean
# case_real_final <- case_real_age_final_1 %>%
#   group_by(V2) %>%
#   summarise(mean_x = mean(case_real_age_final_1))

# group by fold and get mean
val_pred_final <- val_pred_age_final_1 %>%
  group_by(V2) %>%
  summarise(mean_x = mean(val_pred_age_final_1))

# group by fold and get mean
val_real_final <- val_real_age_final_1 %>%
  group_by(V2) %>%
  summarise(mean_x = mean(val_real_age_final_1))


# plot
plot(con_real_final$mean_x, con_pred_final$mean_x, xlim = c(0, 800), ylim = c(0, 800))
abline(0, 1)
cor(con_real_final$mean_x, con_pred_final$mean_x)

# plot
plot(case_pred_age_final_1, case_real_age_final[1,], xlim = c(0, 800), ylim = c(0, 800))
abline(0, 1)
cor(case_pred_age_final_1, case_real_age_final_1)

# plot
plot(val_real_final$mean_x, val_pred_final$mean_x, xlim = c(0, 800), ylim = c(0, 800))
abline(0, 1)
cor(val_real_final$mean_x, val_pred_final$mean_x)


# 1 test.predictions
# 2 test_y
# 3 patient_age
# 4 test.predictions_valid
# 5 valid_y
# 6 patient_age_valid
# 7 test.predictions_controls
# 8 patient_age_controls


