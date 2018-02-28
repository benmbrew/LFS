library(tidyverse)
library(reshape2)
library(ROCR)

# 72_noob_body_false_true_false_true_true_false_false_false

#1stExon   3'UTR   5'UTR    Body TSS1500  TSS200
age_cutoff = 72
method = 'noob'
cg_gene_regions <- "Body"
survival = F
remove_age_lit = T
remove_age_cgs = T
gender = T
tech = T
base_change = F
exon_intron = F
control_for_family = F
k_folds = 5
beta_thresh = 0.05

# 72_noob_body_false_true_false_true_true_false_false_false

temp_results <- readRDS(paste0('../../Data/results_data/',age_cutoff,'_',method,'_', cg_gene_regions,'_',
                               survival ,'_', remove_age_lit,'_', remove_age_cgs ,'_',gender, '_', tech, '_', 
                               base_change,'_',exon_intron, '_', control_for_family,'.rda'))


# get results from list 
temp_cases <- temp_results[[1]]

temp_cases <- temp_cases[ , c('test_pred', 'test_label', 
                                 'age_diagnosis' ,  'age_sample_collection')]
temp_cases$test_pred <- ifelse(temp_cases$test_pred > .5, 1, 0)

temp_cases$pred_is <- ifelse(temp_cases$test_pred == temp_cases$test_label, 
                                'good',
                                'bad')


temp_controls <- temp_results[[2]]

temp_controls <- temp_controls[ , c('controls_age_pred', 'controls_age_label', 
                                    'age_sample_collection')]
temp_controls$controls_pred_label <- ifelse(temp_controls$controls_age_pred > .5, 1, 0)

temp_controls <- temp_controls[ 1:44,]

temp_controls$pred_is <- ifelse(temp_controls$controls_pred_label == temp_controls$controls_age_label, 
                                     'good',
                                     'bad')

temp_controls <- temp_controls[temp_controls$age_sample_collection < 73,]
# get person identfier 
temp_controls$p_id <- rep.int(seq(1, 44, 1), 5)

# group by fold and get mean 
temp_controls_pred <- temp_controls %>%
  group_by(p_id) %>%
  summarise(mean_pred = mean(controls_age_pred, na.rm =T)) %>%
  cbind(temp_controls[1:44,])

# remove original prediction 
temp_controls_pred$controls_age_pred <- NULL
rm(temp_controls,temp_results)


##########
# examin cases with prediction objects (ROC, TRP, etc)
##########

temp_pred_cases <- prediction(temp_cases$test_pred, temp_cases$test_label)
temp_roc <- performance(temp_pred_cases, measure = 'tpr', x.measure = 'tnr')
cutoffs <- data.frame(cut=temp_roc@alpha.values[[1]], tpr=temp_roc@x.values[[1]], 
                      tnr=temp_roc@y.values[[1]])

plot(cutoffs$cut, cutoffs$tpr,type="l",col="red")
par(new=TRUE)
plot(cutoffs$cut, cutoffs$tnr,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")

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

