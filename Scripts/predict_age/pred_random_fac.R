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

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'raw'
k = 5

##########
# load data
##########
# read each data set in
betaFull <- readRDS(paste0(model_data, paste0('/', method, '_', 'full_mod.rda')))

##########
# read in cluster labels
##########
kmeans_lab <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_kmeans_labs.rda')))


# set fixed variables 
model_names <- c('enet', 'rf')
feature_length <- c(20, 50, 100, 300, 500, 1000, 5000, 10000)
seeds <- c(1, 2, 3, 4, 5)


# i = j = m = l = 1
# #
# beta_full = betaFull
# model_method = model_names[m]
# cluster_feats = kmeans_lab
# mod_feats = intersect_names
# feature_set = feature_length[l]
# seed_num = seeds[j]
# k = k
# pred_cutoff = 0.5
# age_cutoff = 72


trainTest <- function(beta_full,
                      pred_cutoff,
                      age_cutoff,
                      model_method,
                      mod_feats,
                      feature_set,
                      seed_num, 
                      k) 
{
  
  
  
  set.seed(seed_num)
  # get a column for each dataset indicating the fold
  beta_cases <- getFolds(beta_cases, seed_number = seed_num, k = k)
  
  # remove samples that dont have an age of sample collection
  beta_cases <- beta_cases[complete.cases(beta_cases),]
  beta_valid <- beta_valid[complete.cases(beta_valid),]
  
  beta_controls <- beta_controls[!is.na(beta_controls$age_sample_collection),]
  beta_controls_old <- beta_controls_old[!is.na(beta_controls_old$age_sample_collection),]
  beta_controls_full <- beta_controls_full[!is.na(beta_controls_full$age_sample_collection),]
  
  # list to store results
  temp_results <- list()
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, beta_cases$folds)
    test_index <- !train_index
    
    # sample random 10 from each cluster 
    mod_feats <- (cluster_feats %>% group_by(count) %>% sample_n(feature_set/10))$probe

    }
    
    if(model_method == 'enet') {
      
      mod_result <- runEnetRandFac(training_dat = beta_cases[train_index,], 
                                   controls_dat = beta_controls,
                                   controls_dat_old = beta_controls_old,
                                   controls_dat_full = beta_controls_full,
                                   valid_dat = beta_valid,
                                   test_dat = beta_cases[test_index,], 
                                   age_cutoff = age_cutoff,
                                   pred_cutoff = pred_cutoff,
                                   bh_features = mod_feats,
                                   gender = F)
    } 
    
    if(model_method == 'rf') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      mod_result <- runRfRandFac(training_dat = beta_cases[train_index,], 
                                 controls_dat = beta_controls,
                                 controls_dat_old = beta_controls_old,
                                 controls_dat_full = beta_controls_full,
                                 valid_dat = beta_valid,
                                 test_dat = beta_cases[test_index,], 
                                 age_cutoff = age_cutoff,
                                 pred_cutoff = pred_cutoff,
                                 bh_features = mod_feats,
                                 gender = F)
      
      
    } 
    
    if(model_method == 'svm') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      
      mod_result <- runSvmRandFac(training_dat = beta_cases[train_index,], 
                               controls_dat = beta_controls,
                               controls_dat_old = beta_controls_old,
                               controls_dat_full = beta_controls_full,
                               valid_dat = beta_valid,
                               test_dat = beta_cases[test_index,], 
                               age_cutoff = age_cutoff,
                               pred_cutoff = pred_cutoff,
                               bh_features = mod_feats,
                               gender = F)
      
    } 
    
    if(model_method == 'lasso') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      
      mod_result <- runLassoRandFac(training_dat = beta_cases[train_index,], 
                                    controls_dat = beta_controls,
                                    controls_dat_old = beta_controls_old,
                                    controls_dat_full = beta_controls_full,
                                    valid_dat = beta_valid,
                                    test_dat = beta_cases[test_index,], 
                                    age_cutoff = age_cutoff,
                                    pred_cutoff = pred_cutoff,
                                    bh_features = mod_feats,
                                    gender = F)
      
    } 
    
    
    if(model_method == 'ridge') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      
      mod_result <- runRidgeRandFac(training_dat = beta_cases[train_index,],
                                 controls_dat = beta_controls,
                                 controls_dat_old = beta_controls_old,
                                 controls_dat_full = beta_controls_full,
                                 valid_dat = beta_valid,
                                 test_dat = beta_cases[test_index,], 
                                 age_cutoff = age_cutoff,
                                 pred_cutoff = pred_cutoff,
                                 bh_features = mod_feats,
                                 gender = F)
      
    } 
    
    
    
    # get results 
    temp_results[[i]] <- mod_result
    
    
    # test_stats, test_stats_age,
    # test_stats_controls, test_stats_controls_old, test_stats_controls_full, 
    # test_stats_valid, test_stats_age_vali
   
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
  
  for(l in 1:length(feature_length)) {
    
    for (j in 1:length(seeds)) {
      
      # return list of 2
      mod_results <- trainTest(beta_full = betaFull,
                               pred_cutoff = 0.5,
                               age_cutoff = 72,
                               model_method = model_names[m],
                               cluster_feats = kmeans_lab_scaled,
                               mod_feats = intersect_names,
                               feature_set = feature_length[l],
                               seed_num = seeds[j],
                               k = k)
      
      
      # this is 5 observations (folds for each of the 7 age variables)
      # cases onset
      temp_results_class[[j]] <- get_class_results(mod_results)[[1]]
      temp_results_mat[[j]] <- get_class_results(mod_results)[[2]]

      # five folds, then first element is confusion matrix, and secod element is clss scores

      print(paste0('done with ', j, ' seed'))
    }
    temp_results_class_2[[l]] <- do.call(rbind, temp_results_class)
    temp_results_mat_2[[l]] <- temp_results_mat
    print(paste0('done with ', l, ' random feature set'))
    
  }
  
  full_results_class[[m]] <- do.call(rbind, temp_results_class_2)
  full_results_mat[[m]] <- temp_results_mat_2
  
  print(paste0('done with ', model_names[m], ' models'))
  
}

options(scipen=999)


#########
# unscaled 
#########

# # no cluster
# full_results <- do.call(rbind, full_results_class)
# saveRDS(full_results, paste0(model_data, '/full_results_class.rda'))
full_results <- readRDS(paste0(model_data, '/full_results_class.rda'))
#DONE

# # cluster
# full_results <- do.call(rbind, full_results_class)
# saveRDS(full_results, paste0(model_data, '/full_results_class_clust.rda'))
full_results_c <- readRDS(paste0(model_data, '/full_results_class_clust.rda'))
#DONE

#########
# scaled
#########
# # 
# # no cluster
# full_results <- do.call(rbind, full_results_class)
# saveRDS(full_results, paste0(model_data, '/full_results_class_scaled.rda'))
full_results_scaled <- readRDS(paste0(model_data, '/full_results_class_scaled.rda'))
# 
# # cluster
# full_results <- do.call(rbind, full_results_class)
# saveRDS(full_results, paste0(model_data, '/full_results_class_clust_scaled.rda'))
full_results_c_scaled <- readRDS(paste0(model_data, '/full_results_class_clust_scaled.rda'))

# combine to make full results, but first add indicator
full_results$data_type <- 'no_clust_no_scale'
full_results_c$data_type <- 'yes_clust_no_scale'
full_results_scaled$data_type <- 'no_clust_yes_scale'
full_results_c_scaled$data_type <- 'yes_clust_yes_scale'

full_results <- rbind(full_results,
                      full_results_c,
                      full_results_scaled,
                      full_results_c_scaled)


# remove unneccesarry obejcs 
rm(list = ls(pattern = "*beta|kmeans")) 

# this will return full results will be a list the size of 
# model_length =  3, 
# features size = 10, 
# seed = 5. 
# folds = 5
# result = 2 class table or matrix

colnames(full_results) <- tolower(colnames(full_results))
colnames(full_results) <- gsub(' ', replacement = '_', colnames(full_results))

# group by age outcome type, number of features in model, and the model name and get 
# the mean "sensitivity (true positive, recall)", "specificity (true negative)", "precision", "balanced_accuracy
# FNR = 1 - TPR
# FPR = 1 - TNR
# TRP = TP/(TP + FN)
# Precision = TP/(TP + FP)

# Precision is 
# (TP)/(TP+FP)(TP)/(TP+FP)
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
  group_by(data_type,age_type, feature_num, model_method) %>%
  summarise(mean_acc = mean(balanced_accuracy, na.rm = T),
            mean_prec = mean(precision, na.rm = T),
            mean_tpr = mean(sensitivity, na.rm = T),
            mean_tnr = mean(specificity, na.rm = T))

# remove svm and make fnr and fpr 
temp <- temp[!grepl('svm', temp$model_method),]
temp$mean_fnr <- 1 - temp$mean_tpr
temp$mean_fpr <- 1 - temp$mean_tnr

temp <- as.data.frame(temp)

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


  