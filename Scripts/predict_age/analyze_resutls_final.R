library(tidyverse)
library(reshape2)
library(ROCR)

method = 'funnorm'
age_cutoff = 72
gender = F
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

# get person identfier 
temp_controls$p_id <- rep.int(seq(1, 47, 1), 4)

# group by fold and get mean 
temp_controls_pred <- temp_controls %>%
  group_by(p_id) %>%
  summarise(mean_pred = mean(controls_age_pred, na.rm =T)) %>%
  cbind(temp_controls[1:47,]) 

# remove original prediction 
temp_controls_pred$controls_age_pred <- NULL
rm(temp_controls,temp_results)


##########
# examin cases with prediction objects (ROC, TRP, etc)
##########

temp_pred_cases <- prediction(temp_cases$test_pred, temp_cases$test_label)
temp_roc <- performance(temp_pred_cases, measure = 'tpr', x.measure = 'fpr')
temp_pr <- performance(temp_pred_cases, measure = 'prec', x.measure = 'rec')
temp_lc <- performance(temp_pred_cases, measure="lift")

plot(temp_roc)
plot(temp_pr)
plot(temp_lc)
caret::confusionMatrix(round(temp_cases$test_pred, 0.1), temp_cases$test_label)

##########
# examine controls 
##########

# assigan to each row an indicator if the prediction was good or bad 
temp_controls_pred$pred_label <- ifelse(temp_controls_pred$mean_pred > 0.5, 1, 0)
temp_controls_pred$pred_is <- ifelse(temp_controls_pred$pred_label == temp_controls_pred$controls_age_label, 
                                     'good',
                                     'bad')

table(temp_controls_pred$pred_is)

