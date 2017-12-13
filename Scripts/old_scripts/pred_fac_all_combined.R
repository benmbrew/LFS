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
method = 'funnorm'
k = 5
combat = F

##########
# load data
##########

if(combat) {
  
  betaFullCon <- readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_con_combat.rda')))
  betaFullVal <- readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_val_combat.rda')))
  
  
} else {
  
  betaFullCon <- readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_con.rda')))
  betaFullVal <- readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_val.rda')))
  
}


##########
# get column names
##########
intersect_names_con <- colnames(betaFullCon)[9:ncol(betaFullCon)]
intersect_names_val <- colnames(betaFullVal)[9:ncol(betaFullVal)]


##########
# read in all features from feat data
##########

lfs_feats_m <- readRDS(paste0(feat_data, paste0('/', 'noob', '_', 'lfs_m.rda')))
no_cancer_feats_m <- readRDS(paste0(feat_data, paste0('/', 'noob', '_', 'no_cancer_m.rda')))


# setwd(feat_data)
# file_list_names = list.files()
# 
# 
# # store all raw rda feature lists in feat_list
# feat_list <- lapply(file_list_names, function(x) readRDS(x))
# 
# # order feat list
# feat_list  <- feat_list[order(sapply(feat_list, length), decreasing=F)]
# 
# # select first 10 
# feat_list <- feat_list[5]
# 
# # model_names <- c('enet', 'rf', 'lasso')
# # seeds <- c(1, 2, 3)

model_names <- c('enet')
seeds <- c(1,2)
feat_list <- c(no_cancer_feats_m, lfs_feats_m)
file_list_names <- c('no_cancer_m', 'lfs')



l = j = m = 1
#
beta_full = betaFullCon
model_method = model_names[m]
mod_feats = feat_list[[2]]
feat_name = file_list_names[[l]]
seed_num = seeds[j]
k = k
class_age = 72
min_age = 5
max_age = 850
data_type = 'controls'



trainTest <- function(beta_full,
                      model_method,
                      mod_feats,
                      feat_name,
                      seed_num,
                      data_type,
                      max_age,
                      remove_cases,
                      class_age,
                      k) {
  
  
  if(data_type == 'controls') {
    
    beta_cases <- beta_full[beta_full$batch == 'cases',]
    beta_controls <- beta_full[beta_full$batch == 'controls',]
    beta_valid <- beta_full[beta_full$batch == 'controls',]
    
  } else {
    
    beta_cases <- beta_full[beta_full$batch == 'cases',]
    beta_controls <- beta_full[beta_full$batch == 'valid',]
    beta_valid <- beta_full[beta_full$batch == 'valid',]
    
  }
  
  # remove samples that dont have an age of sample collection
  beta_cases <- beta_cases[complete.cases(beta_cases),]
  # beta_valid <- beta_valid[complete.cases(beta_valid),]
  
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
      
      mod_result <- runEnetRandFac(training_dat = beta_cases[train_index,], 
                                   controls_dat = beta_controls,
                                   valid_dat = beta_valid,
                                   test_dat = beta_cases[test_index,], 
                                   age_cutoff = class_age,
                                   bh_features = mod_feats,
                                   rand_feats = rand_feats,
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
                                 gender = T)
      
      
    } 
    
    if(model_method == 'svm') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      
      mod_result <- runSvmRandFac(training_dat = beta_cases[train_index,], 
                                  controls_dat = beta_controls,
                                  valid_dat = beta_valid,
                                  test_dat = beta_cases[test_index,], 
                                  age_cutoff = age_cutoff,
                                  bh_features = mod_feats,
                                  rand_feats = rand_feats,
                                  gender = F)
      
    } 
    
    if(model_method == 'lasso') {
      mod_result <- runLassoL1RandFac(training_dat = beta_cases[train_index,], 
                                      controls_dat = beta_controls,
                                      valid_dat = beta_valid,
                                      test_dat = beta_cases[test_index,], 
                                      age_cutoff = class_age,
                                      bh_features = mod_feats,
                                      rand_feats = rand_feats,
                                      gender = F)
      
    } 
    
    # get results 
    temp_results[[i]] <- mod_result
    
  }
  
  return(temp_results)
  
}

temp_results_class <- list()
temp_results_mat <- list()
temp_results_class_2 <- list()
temp_results_mat_2 <- list()
full_results_class <- list()
full_results_mat <- list()
# set fixed variables 

for(m in 1:length(model_names)) {
  
  for(l in 1:length(feat_list)) {
    
    for (j in 1:length(seeds)) {
      
      # return list of 2
      mod_results <- trainTest(beta_full = betaFullVal[, 1:50000],
                               model_method = model_names[m],
                               mod_feats = feat_list[[l]],
                               feat_name = file_list_names[[l]],
                               seed_num = seeds[j],
                               data_type = 'valid',
                               max_age = 1000,
                               remove_cases = F,
                               class_age = 72,
                               k)
      
      
      # this is 5 observations (folds for each of the 7 age variables)
      # cases onset
      temp_results_class[[j]] <- get_class_results(mod_results, dims_of_dat = length(feat_list[[l]]), 
                                                   mod_name = model_names[m], feat_name = file_list_names[[l]],
                                                   seed_number = seeds[j])[[1]]
      temp_results_mat[[j]] <- get_class_results(mod_results, dims_of_dat = length(feat_list[[l]]), 
                                                 mod_name = model_names[m], feat_name = file_list_names[[l]],
                                                 seed_number = seeds[j])[[2]]
      
      # five folds, then first element is confusion matrix, and secod element is clss scores
      
      print(paste0('done with ', j, ' seed'))
    }
    temp_results_class_2[[l]] <- do.call(rbind, temp_results_class)
    temp_results_mat_2[[l]] <- temp_results_mat
    print(paste0('done with ', l, ' features'))
    
    
  }
  full_results_class[[m]] <- do.call(rbind, temp_results_class_2)
  full_results_mat[[m]] <- temp_results_mat_2
  print(paste0('done with ', m, ' model'))
  
}

# oct 19th 11 am.

# options(scipen=999)

# unscaled 
#########

# # no cluster

full_results <- do.call(rbind, full_results_class)


saveRDS(full_results, paste0(results_data, '/', method,'_','full_results_class_combined_m.rda'))


colnames(full_results) <- tolower(colnames(full_results))
colnames(full_results) <- gsub(' ', replacement = '_', colnames(full_results))

# group by age outcome type, number of features in model, and the model name and get 
# the mean "sensitivity (true positive, recall)", "specificity (true negative)", "precision", "balanced_accuracy
# FNR = 1 - TPRs
# FPR = 1 - TNR
# TRP = TP/(TP + FN)
# Precision = TP/(TP + FP)

# Precision is 
# (TP)/(TP+FP)
# which tells us what proportion of patients we diagnosed as having cancer actually had cancer. 
# In other words, proportion of TP in the set of positive cancer diagnoses. This is given by the rightmost 
# column in the confusion matrix.

# Recall is
# (TP)/(TP+FN)
# which tells us what proportion of patients that actually had cancer were diagnosed by us as having cancer. 
# In other words, proportion of TP in the set of true cancer states. This is given by the bottom row in 
# the confusion matrix.

# Recall (TPR, sensitvity) is the probability that a (randomly selected) relevant document is retrieved in a search.
# Precision is the probability that a (randomly selected) retrieved document is relevant.

# In this representation, it is clearer that recall gives us information about a classifiers 
# performance with respect to false negatives (how many did we miss), while precision gives us 
# information about its performance with respect to false positives.

temp <- 
  full_results %>%
  group_by(age_type, feature_num, feat_name, model_method) %>%
  summarise(mean_acc = mean(balanced_accuracy, na.rm = T),
            mean_prec = mean(precision, na.rm = T),
            mean_tpr = mean(sensitivity, na.rm = T),
            mean_tnr = mean(specificity, na.rm = T))

# remove svm and make fnr and fpr 
temp <- temp[!grepl('svm', temp$model_method),]
temp$mean_fnr <- 1 - temp$mean_tpr
temp$mean_fpr <- 1 - temp$mean_tnr

temp <- as.data.frame(temp)

temp <- temp[order(temp$feature_num, decreasing = F),]

saveRDS(temp, paste0(results_data, '/', method,'_','full_results_class_collapsed_under_18_48.rda'))
temp <- temp[order(temp$mean_acc, decreasing = T),]


# function for plotting class results 
plot_class_results <- function(results_data, age_type, data_type, plot_column) {
  
  max_feat <- max(results_data$feature_num)
  plot_temp <- as.data.frame(results_data[grepl(age_type, results_data$age_type),])
  plot_temp$age_type <- NULL
  plot_temp$data_type <- NULL
  
  plot_melt <- melt(plot_temp, id.vars = c('feature_num', 'model_method'))
  plot_dat <- plot_melt[grepl(plot_column, plot_melt$variable),]
  plot_dat$value <- as.numeric(plot_dat$value)
  
  ggplot(plot_dat, aes(feature_num, value, group = model_method, colour = model_method)) +
    geom_point(size = 4) + 
    geom_line(size = 2, alpha = 0.4) + 
    xlab('number of features') + 
    ylab(unique(plot_dat$variable)) +
    xlim(c(0, max_feat)) + ylim(c(0,1)) +
    scale_fill_manual(name = 'Model',
                      breaks = c('enet', 'rf'),
                      labels = c('Enet', 'Random forest'),
                      values = c('blue', 'red')) +
    theme_bw() + theme(axis.text.y = element_text(size = 12),
                       axis.text.x = element_text(size = 12)) + ggtitle(paste0(age_type, '_', data_type, '_', plot_column))
}


summary(as.factor(temp$age_type))
summary(as.factor(temp$data_type))

plot_pdf <- function(data_type){
  
  temp <- temp[grepl(data_type, temp$data_type),]
  
  pdf(paste0('~/Desktop/', data_type ,'.pdf'))
  
  for (i in unique(temp$age_type)) {
    
    print(plot_class_results(results_data = temp , age_type = i, data_type = data_type, plot_column = 'mean_fpr'))
    print(plot_class_results(results_data = temp , age_type = i, data_type = data_type, plot_column = 'mean_fnr'))
    print(plot_class_results(results_data = temp , age_type = i, data_type = data_type, plot_column = 'mean_acc'))
    
  }
  
  dev.off()
  
}




plot_pdf('no_clust_no_scale')
plot_pdf('no_clust_yes_scale')
plot_pdf('yes_clust_no_scale')
plot_pdf('yes_clust_yes_scale')


