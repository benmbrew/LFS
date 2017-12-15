library(tidyverse)
library(reshape2)

method = 'noob'
age_cutoff = 72
gender = T
tech = T
k_folds = 4
beta_thresh = 0.01
control_type = 'full'


# load result data for 450_850 analysis
temp_results <- readRDS(paste0('../../Data/results_data/',method, '_', 
                               age_cutoff, '_', gender, '_', tech,'_' , 
                               control_type, '_', beta_thresh, '.rda'))



# get results from list 
temp_cases <- temp_results[[1]]
temp_controls <- temp_results[[2]]

# cases are ok, but controls need to be avg over 4 folds 
# create indicator for fold (1-47, by 4)
temp_controls$fold <- 
temp_controls <- temp_controls %>%
  group_by()

# get cases, controls and valid 
cases <- full_results[[1]]
controls <- full_results[[2]]
valid <- full_results[[3]]

# valid and controls are multiplied by 5 
controls <- controls[1:30,]
valid <- valid[1:36,]

# remove unneeded objects
rm(training_clin, training_dat, temp.cv_error_matrix, temp.cv_error_mean, test_dat, test_clin)
rm(beta_controls, beta_controls_mod, beta_controls_mod_old, beta_test, beta_test_cases, beta_train)
rm(beta_valid, beta_train_cases, bh_feats, bh_features, beta_valid_mod)
rm(rg_cases, rg_controls, rg_valid, rgCasesM, rgCasesT, rgCases, rgControls, rgValid)

#########
# cases controls and valid 
#########
cases$score <- ifelse(round(cases$test_pred,.1) == cases$test_label, 'good', 'bad')
caret::confusionMatrix(round(cases$test_pred, 0.1), cases$test_label)

# get pred object from predictions and real values 
temp_pred_cases <- prediction(cases$test_pred, cases$test_label)
temp_roc <- performance(temp_pred_cases, measure = 'tpr', x.measure = 'fpr')
temp_pr <- performance(temp_pred_cases, measure = 'prec', x.measure = 'rec')
temp_lc <- performance(temp_pred_cases, measure="lift")




# get score 
controls$score <- ifelse(controls$test_pred == controls$test_label, 'good', 'bad')

# get young samples 
young <- controls[controls$age_sample_collection <216,]

# get bad samples 
bad <- young[young$score == 'bad',]


# get confusion matrix
caret::confusionMatrix(controls$test_pred, controls$test_label)

# get pred object from predictions and real values 

temp_pred_controls <- prediction(controls$test_pred, controls$test_label)
temp_roc <- performance(temp_pred_controls, measure = 'tpr', x.measure = 'fpr')
temp_pr <- performance(temp_pred_controls, measure = 'prec', x.measure = 'rec')
temp_lc <- performance(temp_pred_controls, measure="lift")

plot(temp_roc)
