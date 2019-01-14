# source functions script
source('helper_functions.R')

# load libraries
library(plotly)
library(scatterplot3d)
library(tidyverse)
library(grid)
library(broom)
library(scales)
library(gridExtra)
library(data.table)
library(doParallel)

# register other cpus
registerDoParallel(2)

# set fixed variables
size = 'full'
model_type = 'rf'
gender = FALSE
method = 'noob'
methyl_type = 'beta'
beta_thresh = 0.01
optimal_cutoff = 0.5

# create objects to indicate method and model details when saving
age_cutoff = 72
trained_lambda = FALSE
tech = FALSE



if(trained_lambda){
  is_lambda <- 'lambda_test'
} else {
  is_lambda <- 'lambda_train'
}


if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}

# read in cases_450
if(model_type == 'enet'){
  temp_con <- readRDS(paste0('transform_data_test/', 'con_test_transform_',method,'_',size,'_',is_gen, '_', model_type,'.rda'))
  
  temp_valid <- readRDS(paste0('transform_data_test/', 'valid_test_transform_',method,'_',size,'_',is_gen, '_', model_type,'.rda'))
  
} else {
  
  temp_con <- readRDS(paste0('transform_data_test/', 'con_test_transform',method,'_',size,'_',is_gen, '_', model_type,'.rda'))
  
  temp_valid <- readRDS(paste0('transform_data_test/', 'valid_test_transform',method,'_',size,'_',is_gen, '_', model_type,'.rda'))
  
  temp_importance <- readRDS(paste0('transform_data_test/', 'importance_transform',method,'_',size,'_',is_gen, '_', model_type,'.rda'))
  
  
}

# CONTROLS
if(model_type == 'enet'){
  temp_group <- temp_con %>% group_by(tm_donor, alpha) %>% 
    filter(p53_germline == 'MUT') %>%
    mutate(mean_pred = mean(controls_age_pred, na.rm = TRUE))
   
} else {
  temp_group <- temp_con %>% group_by(tm_donor) %>% 
    filter(p53_germline == 'MUT') %>%
    mutate(mean_pred = mean(positive, na.rm = TRUE))
  
}

write.csv(temp_group, paste0('~/Desktop/lfs_plots_jan_2019/test_transform/','con_test_', method,'_',size,'_',is_gen, '_', is_lambda, '_', model_type,'.csv'))

if(model_type == 'enet'){
  # remove duplicates
  temp <- get_acc_val(temp_valid, thresh = 0.5)
  temp <- temp[order(temp$acc, decreasing = TRUE),]
  write.csv(temp_group, paste0('~/Desktop/lfs_plots_jan_2019/test_transform/','val_test_', method,'_',size,'_',is_gen, '_', is_lambda,'_', model_type,'.csv'))
  
} else {
  write.csv(temp_valid, paste0('~/Desktop/lfs_plots_jan_2019/test_transform/','val_test_', method,'_',size,'_',is_gen, '_', is_lambda,'_', model_type,'.csv'))
  
}

# temp_dedup <- temp[!duplicated(temp$tm_donor),]
# max_alpha = temp_dedup$alpha[which(temp_dedup$acc == max(temp_dedup$acc))]

# 
# temp <- temp[temp$alpha == 0.7,]
# temp <- temp[!duplicated(temp$tm_donor),]
# write.csv(temp, '~/Desktop/strategy_1_log/valid_table_normal.csv')
# 
# # get dataset of predictions and labels for both small and large data
# pred_short <- prediction(temp$valid_age_pred, temp$valid_age_label)
# # pred_long <- prediction($preds, final_dat$real)
# 
# # get performace objects
# perf_s <- performance(prediction.obj = pred_short, measure = 'tpr', x.measure = 'fpr')
# # perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')
# 
# # plot mean preds
# plot(perf_s)
# abline(a = 0, b =1)
# 
# # plot all preds
# plot(perf_l)
# abline(a = 0, b =1)
# 
# 
# # plot optimal cutoff (TPR vs TNR)
# temp_roc <- performance(pred_short, measure = 'tpr', x.measure = 'tnr')
# cutoffs <- data.frame(cut=temp_roc@alpha.values[[1]], tpr=temp_roc@x.values[[1]],
#                       tnr=temp_roc@y.values[[1]])
# 
# # combine plots 
# plot(cutoffs$cut, cutoffs$tpr,type="l",col="red")
# par(new=TRUE)
# plot(cutoffs$cut, cutoffs$tnr,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
# 
# # get precision recall curves 
# temp_pr <- performance(pred_short, measure = 'prec', x.measure = 'rec')
# plot(temp_pr)
# 
# # get lift curve
# temp_lc <- performance(pred_short, measure="lift")
# plot(temp_lc)
# 
# 
# ###############
# # get optimal cutoff and confusion matrix info 
# ###############
# 
# # optimal cutoff
# cost.perf = performance(pred_short, "cost", cost.fp = 1, cost.fn = 1)
# optimal_cutoff <- pred_short@cutoffs[[1]][which.min(cost.perf@y.values[[1]])]
# 
# # get confusion matrix function for plotting 
# ConfusionMatrixInfo(data = temp, 
#                     predict = 'valid_age_pred', 
#                     actual = 'valid_age_label', 
#                     cutoff = .6,
#                     get_plot = TRUE)
# 
# 
# 
# dev.off()
# dat_sample$real_label <- as.character(dat_sample$real_label)
# # get cost function by specifying cost before
# cost_fp <- 1
# cost_fn <- 1
# roc_info <- ROCInfo(data = dat_sample, 
#                     predict = 'mean_preds', 
#                     actual = 'real_label', 
#                     cost.fp = cost_fp,
#                     cost.fn = cost_fn)
# 
# grid.draw(roc_info$plot)
# 
# ###############
# # get acc optimal cutoff plot
# ###############
# dev.off()
# acc.perf = performance(pred_short, measure = "acc")
# plot(acc.perf)
# abline(v = optimal_cutoff)
# 
# ###############
# # Other plots
# ###############
# 
# # distribution of preds
# ggplot(dat_sample, aes(mean_preds)) + 
#   geom_histogram(aes(y=..density..), bins = 15, colour="black", fill="white")+
#   geom_density(adjust = 1, fill= 'black', alpha = 0.4) +
#   labs(title = 'Distribution of risk scores',
#        x = 'Mean risk',
#        y = 'Density')
# 
# # plot age as a function of prediction
# ggplot(dat_sample, 
#        aes(age, mean_preds)) + 
#   geom_point() +
#   labs(title = 'Age in months vs risk score',
#        x = 'Age (months)',
#        y = 'Mean risk score') +
#   geom_smooth(method = 'loess')
# 
# 
# t_dat <- dat[dat$alpha == 0.5,]
# # get confusion matrix function for plotting 
# t <-ConfusionMatrixInfo(data = t_dat, 
#                         predict = 'valid_age_pred', 
#                         actual = 'valid_age_label', 
#                         cutoff = .6,
#                         get_plot = FALSE)
# 
# val_dat <- as.data.frame(cbind(t_dat, t))
# val_dat$valid_age_pred <- val_dat$vallid_age_label <- NULL
# 
# 
# temp <-read.csv('~/Desktop/strategy_1_log/con_table_log.csv')
# plot(temp$age_sample_collection, temp$controls_age_pred)
# a
# 
