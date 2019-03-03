
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
library(pROC)
# register other cpus
registerDoParallel(2)
# Build a ROC object and compute the AUC

# set fixed variables
size = 'used_bh'
model_type = 'rf'
# standardize = FALSE
gender = FALSE
method = 'quan'
combat = 'combat_1'
which_methyl = 'beta'
beta_thresh = 0.05
optimal_cutoff = 0.5
remove_leading_pcs = 'first'
tech = FALSE


# create objects to indicate method and model details when saving
age_cutoff = 72
trained_lambda = FALSE
tech = FALSE

source('helper_functions.R')

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

if(model_type == 'enet'){
  # read in cases_450
  temp_con <- readRDS(paste0('age_combat_data_test/', 'con_test_pc_2',method,'_',size,'_',is_gen,'_',combat,'_', model_type,'.rda'))
  
  temp_valid <- readRDS(paste0('age_combat_data_test/', 'valid_test_pc_2',method,'_',size,'_',is_gen, '_',combat, '_', model_type,'.rda'))
  
  
  temp_valid <- temp_valid[, c('tm_donor','p53_germline','valid_age_label',  'age_diagnosis', 'age_sample_collection', 'valid_age_pred','non_zero', 'alpha', 'tech')]
  temp_con <- temp_con[, c('tm_donor','p53_germline','controls_age_label',  'age_sample_collection', 'controls_age_pred','non_zero', 'alpha', 'tech')]
  
} else {
  temp_con <- readRDS(paste0('age_combat_data_test/', 'con_test_pc_2', method,'_',size,'_',is_gen,'_',combat,'_', model_type,'.rda'))
  
  temp_valid <- readRDS(paste0('age_combat_data_test/', 'valid_test_pc_2',method,'_',size,'_',is_gen,'_',combat,'_', model_type,'.rda'))
  
  temp_importance <- readRDS(paste0('age_combat_data_test/', 'importance_pc_2',method,'_',size,'_',is_gen,'_',combat,'_', model_type,'.rda'))
  temp_valid <- temp_valid[, c('tm_donor','p53_germline','positive','real' ,'age_diagnosis', 'age_sample_collection', 'tech')]
  temp_con <- temp_con[, c('tm_donor','p53_germline','positive', 'real', 'age_sample_collection', 'tech')]
  
}

# CONTROLS
if(model_type == 'enet'){
  
  
  temp_group <- temp_con %>% group_by(tm_donor, alpha) %>% 
    filter(p53_germline == 'MUT') %>%
    mutate(mean_pred = mean(controls_age_pred, na.rm = TRUE))
  
  temp_450 <- temp_group[temp_group$tech == '450k',]
  temp_850 <- temp_group[temp_group$tech == '850k',]
  
  
  temp_alpha <- temp_con %>% group_by(alpha) %>% 
    filter(p53_germline == 'MUT') %>%
    summarise(mean_pred_under = mean(controls_age_pred[age_sample_collection < 72], na.rm = TRUE),
              mean_pred_over = mean(controls_age_pred[age_sample_collection > 72], na.rm = TRUE),
              mean_pred_450 = mean(controls_age_pred[tech == '450k'], na.rm = TRUE),
              mean_pred_850 = mean(controls_age_pred[tech == '850k'], na.rm = TRUE))
  
  
} else {
  temp_group <- temp_con %>% group_by(tm_donor) %>% 
    filter(p53_germline == 'MUT') %>%
    mutate(mean_pred = mean(positive, na.rm = TRUE))
  
  temp_450 <- temp_group[temp_group$tech == '450k',]
  temp_850 <- temp_group[temp_group$tech == '850k',]
  temp_450$real_label <- as.factor(ifelse(temp_450$age_sample_collection > 72, 'negative', 'positive'))
  temp_850$real_label <- as.factor(ifelse(temp_850$age_sample_collection > 72, 'negative', 'positive'))
  temp_450$real_label <- factor(temp_450$real_label, levels = c('positive', 'negative'))
  temp_850$real_label <- factor(temp_850$real_label, levels = c('positive', 'negative'))
  
}

# write.csv(temp_group, paste0('~/Desktop/lfs_plots_jan_2019/test_pc/','con_test_', method,'_',size,'_',is_gen, '_', is_lambda, '_',combat,'_', model_type,'.csv'))

if(model_type == 'enet'){
  # remove duplicates
  temp <- get_acc_val(temp_valid, thresh = 0.5)
  temp <- temp[order(temp$acc, decreasing = TRUE),]
  # write.csv(temp_group, paste0('~/Desktop/lfs_plots_jan_2019/test_pc/','val_test_', method,'_',size,'_',is_gen, '_', is_lambda, '_',combat,'_', model_type,'.csv'))
  
} else {
  # write.csv(temp_valid, paste0('~/Desktop/lfs_plots_jan_2019/test_pc/','val_test_', method,'_',size,'_',is_gen, '_', is_lambda, '_',combat,'_', model_type,'.csv'))
  # size = 'used_bh'
  # model_type = 'rf'
  # gender = FALSE
  # method = 'swan'
  # combat = 'normal'
  # optimal cutoff
  pred_val <- prediction(temp_valid$positive, temp_valid$real)
  
  cost.perf = performance(pred_val, "cost", cost.fp = 1, cost.fn = 1)
  optimal_cutoff <- pred_val@cutoffs[[1]][which.min(cost.perf@y.values[[1]])]
  
  # get confusion matrix function for plotting
  ConfusionMatrixInfo(data = temp_valid,
                      other_title = 'validation set - x axis positive if age of onset is less than 72 months',
                      predict = 'positive',
                      actual = 'real',
                      cutoff = .55,
                      get_plot = TRUE)
  
  ConfusionMatrixInfo(data = temp_850,
                      other_title = 'null set 850 - x axis positive if age of collection is less than 72 months',
                      predict = 'mean_pred',
                      actual = 'real_label',
                      cutoff = .55,
                      get_plot = TRUE)
  
  
  ConfusionMatrixInfo(data = temp_450,
                      other_title = 'null set 450 - x axis positive if age of collection is less than 72 months',
                      predict = 'mean_pred',
                      actual = 'real_label',
                      cutoff = .55,
                      get_plot = TRUE)
  
  temp_all <- rbind(temp_450, temp_850)
  ConfusionMatrixInfo(data = temp_all,
                      other_title = 'null set 450 - x axis positive if age of collection is less than 72 months',
                      predict = 'mean_pred',
                      actual = 'real_label',
                      cutoff = .55,
                      get_plot = TRUE)
  
}



# plots
alpha = 0.5

# subset data by 0.2
temp <- temp[temp$alpha == alpha,]
temp_850 <- temp_850[temp_850$alpha == alpha,]

# get dataset of predictions and labels for both small and large data
pred_val <- prediction(temp$valid_age_pred, temp$valid_age_label)
pred_con <- prediction(temp_850$controls_age_pred, temp_850$controls_age_label)

# get performace objects
perf_s <- performance(prediction.obj = pred_val, measure = 'tpr', x.measure = 'fpr')
# perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')

# plot mean preds
dev.off()
plot(perf_s)
abline(a = 0, b =1)
roc(temp$valid_age_label, temp$valid_age_pred)
legend(x = 0.9, y = 0.2, legend = .86)



# plot optimal cutoff (TPR vs TNR)
temp_roc <- performance(pred_val, measure = 'tpr', x.measure = 'tnr')
cutoffs <- data.frame(cut=temp_roc@alpha.values[[1]], tpr=temp_roc@x.values[[1]],
                      tnr=temp_roc@y.values[[1]])

# combine plots
plot(cutoffs$cut, cutoffs$tpr,type="l",col="red")
par(new=TRUE)
plot(cutoffs$cut, cutoffs$tnr,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")

# # get precision recall curves
# temp_pr <- performance(pred_short, measure = 'prec', x.measure = 'rec')
# plot(temp_pr)
# 
# # get lift curve
# temp_lc <- performance(pred_short, measure="lift")
# plot(temp_lc)


###############
# get optimal cutoff and confusion matrix info

# optimal cutoff
cost.perf = performance(pred_val, "cost", cost.fp = 1, cost.fn = 1)
optimal_cutoff <- pred_val@cutoffs[[1]][which.min(cost.perf@y.values[[1]])]

# get confusion matrix function for plotting
ConfusionMatrixInfo(data = temp,
                    other_title = 'validation set',
                    predict = 'valid_age_pred',
                    actual = 'valid_age_label',
                    cutoff = .53,
                    get_plot = TRUE)

ConfusionMatrixInfo(data = temp_850,
                    other_title = 'cross_validation_null set',
                    predict = 'mean_pred',
                    actual = 'controls_age_label',
                    cutoff = optimal_cutoff,
                    get_plot = TRUE)

# get confusion matrix function for plotting
ConfusionMatrixInfo(data = temp,
                    other_title = 'cross validation validation set',
                    predict = 'valid_age_pred',
                    actual = 'valid_age_label',
                    cutoff = .53,
                    get_plot = TRUE)

ConfusionMatrixInfo(data = temp_850,
                    other_title = 'cross_validation_null set',
                    predict = 'mean_pred',
                    actual = 'controls_age_label',
                    cutoff = .53,
                    get_plot = TRUE)



dev.off()
dat_sample$real_label <- as.character(dat_sample$real_label)
# get cost function by specifying cost before
cost_fp <- 1
cost_fn <- 1
roc_info <- ROCInfo(data = temp,
                    predict = 'valid_age_pred',
                    actual = 'valid_age_label',
                    other_title = 'validation set',
                    cost.fp = cost_fp,
                    cost.fn = cost_fn)

grid.draw(roc_info$plot)

###############
# get acc optimal cutoff plot
###############
dev.off()
acc.perf = performance(pred_short, measure = "acc")
plot(acc.perf)
abline(v = optimal_cutoff)

###############
# Other plots
###############

# distribution of preds
ggplot(dat_sample, aes(mean_preds)) +
  geom_histogram(aes(y=..density..), bins = 15, colour="black", fill="white")+
  geom_density(adjust = 1, fill= 'black', alpha = 0.4) +
  labs(title = 'Distribution of risk scores',
       x = 'Mean risk',
       y = 'Density')

# plot age as a function of prediction
ggplot(temp,
       aes(age_sample_collection, valid_age_pred)) +
  geom_point() +
  labs(title = 'Age in months vs risk score',
       x = 'Age (months)',
       y = 'Mean risk score') +
  geom_smooth(method = 'loess')


t_dat <- dat[dat$alpha == 0.5,]
# get confusion matrix function for plotting
t <-ConfusionMatrixInfo(data = t_dat,
                        predict = 'valid_age_pred',
                        actual = 'valid_age_label',
                        cutoff = .6,
                        get_plot = FALSE)

val_dat <- as.data.frame(cbind(t_dat, t))
val_dat$valid_age_pred <- val_dat$vallid_age_label <- NULL


