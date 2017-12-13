
##########
# load libraries
##########
library(tidyverse)

##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
idat_data_case <- paste0(methyl_data, '/raw_files')
idat_data_con <- paste0(methyl_data, '/controls')
idat_data_val <- paste0(data_folder, '/methyl_data/validation/idat_files')
model_data <- paste0(data_folder, '/model_data')
map_data <- paste0(data_folder, '/methyl_data/validation')


##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

# result names 
# raw_5_100_50000_0.5_FALSE_FALSE_full_results,
# raw_5_100_50000_0.5_TRUE_controls_full_results,
# raw_5_100_50000_0.5_TRUE_valid_full_results,
# funnorm_5_100_50000_0.5_FALSE_FALSE_full_results,
# noob_5_100_50000_0.5_FALSE_FALSE_full_results,
# noob_5_100_50000_0.5_TRUE_valid_full_results,

# read in results
# list of 2 - results and clin
# 

# raw, funnorm, noob, swan
# fixed variables
method = 'noob'
combined = T
combined_data_type = 'controls'

# raw no combined
full_results <- readRDS(paste0('~/Desktop/',method, '_5_1000_50000_0.5_' , combined,'_', combined_data_type, '_no_gender_full_results.rda'))[[1]]
full_results_clin <- readRDS(paste0('~/Desktop/',method, '_5_1000_50000_0.5_' , combined,'_', combined_data_type, '_no_gender_full_results.rda'))[[2]]

# create real value column
full_results_clin$ground_truth[full_results_clin$batch == 'test_data' & full_results_clin$age_diagnosis < 72] <- 'a'
full_results_clin$ground_truth[full_results_clin$batch == 'test_data' & full_results_clin$age_diagnosis > 72] <- 'b'

full_results_clin$ground_truth[full_results_clin$batch == 'controls_data' & full_results_clin$age_sample_collection < 72] <- 'a'
full_results_clin$ground_truth[full_results_clin$batch == 'controls_data' & full_results_clin$age_sample_collection > 72] <- 'b'

full_results_clin$ground_truth[full_results_clin$batch == 'valid_data' & full_results_clin$age_diagnosis < 72] <- 'a'
full_results_clin$ground_truth[full_results_clin$batch == 'valid_data' & full_results_clin$age_diagnosis > 72] <- 'b'

# keep only necessary columns
full_results_clin$F <- full_results_clin$M <- full_results_clin$p53_germline <-
  full_results_clin$ids <- full_results_clin$sentrix_id <- NULL

# create column that indicates of good or bad pred
full_results_clin$result <- ifelse(full_results_clin$prediction == full_results_clin$ground_truth, 'good', 'bad')

# subset by batch and look at pred vs real with age
sub_batch <- subset(full_results_clin, batch == 'controls_data')

# make indicator for fold 

sub_batch$folds[1:30] <- '1'
sub_batch$folds[31:60] <- '2'
sub_batch$folds[61:90] <- '3'
sub_batch$folds[91:120] <- '4'
sub_batch$folds[121:150] <- '5'

# keep only relevant columns 
sub_batch$cancer_diagnosis_diagnoses <- sub_batch$age_diagnosis <- NULL 

# group by fold and ground truth and get sum of goods and bads
good_bad <- sub_batch %>%
  group_by(folds, ground_truth, prediction) %>%
  summarise(counts = n())



# edit results and clin results
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
  group_by(age_type, model_method) %>%
  summarise(mean_feat = mean(feature_num, na.rm = T),
            mean_acc = mean(balanced_accuracy, na.rm = T),
            mean_prec = mean(precision, na.rm = T),
            mean_tpr = mean(sensitivity, na.rm = T),
            mean_tnr = mean(specificity, na.rm = T))

# remove svm and make fnr and fpr 
temp$mean_fnr <- 1 - temp$mean_tpr
temp$mean_fpr <- 1 - temp$mean_tnr

temp <- as.data.frame(temp)

