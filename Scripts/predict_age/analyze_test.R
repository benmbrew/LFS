# This scrip will analyze cross validation error and find optimal cut off for classification
# https://ethen8181.github.io/machine-learning/unbalanced/unbalanced.html
# https://hopstat.wordpress.com/2014/12/19/a-small-introduction-to-the-rocr-package/
# https://pat-s.github.io/mlr/articles/tutorial/devel/roc_analysis.html
library(plotly)
library(scatterplot3d)
library(car)
library(rgl)
library(ROCR)
library(caret)
library(pROC)
library(dplyr)
library(grid)
library(broom)
library(tidyr)
library(scales)
library(gridExtra)
library(data.table)

source('helper_functions.R')

# source functions script
# source('all_functions.R')

# create plotting functions


# create fixed objects to model and pipeline inputs and saving  
methyl_type = 'beta'
combat = TRUE
gender = FALSE
tech = FALSE 
how_many_seeds = 100
how_many_folds = 5


# condition on fixed objects to get saving identifiers

if(methyl_type == 'beta'){
  which_methyl <- 'beta'
} else {
  which_methyl <- 'm'
}

if(combat){
  which_combat <- 'use_combat'
  beta_thresh = 0.05
} else {
  which_combat <- 'no_combat'
  beta_thresh = 0.5
}

if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}

num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds
alpha_num = 1

# read in cases_450
final_dat <- readRDS(paste0('validation_age_predictions/', 'valid_test_untrained',alpha_num,'_' ,which_methyl, '_', which_combat, '_',
                     num_seeds, '_', num_folds, '_', is_gen, '.rda'))

# subset data
dat <- final_dat[, c('tm_donor','cancer_diagnosis_diagnoses','valid_age_pred', 'vallid_age_label', 'age_sample_collection',
                         'alpha')]


dat <- dat[dat$alpha == 1,]
# get confusion matrix function for plotting 
t <-ConfusionMatrixInfo(data = dat, 
                    predict = 'valid_age_pred', 
                    actual = 'vallid_age_label', 
                    cutoff = .5,
                    get_plot = FALSE)


# create a function for controls that gets the young predictions for each alpha level
temp_dat <- dat
thresh = .5
i = 1

age = 72
get_young_labels <- function(temp_dat, thresh, age){
  all_alphas <- (1:10)/10
  result_list <- list()
  for(i in 1:length(all_alphas)){
    this_alpha <- all_alphas[i]
    sub_dat <- temp_dat[temp_dat$alpha == this_alpha,]
    sub_dat <- sub_dat[sub_dat$age_sample_collection < age,]
    sub_dat$pred_label <-as.factor(ifelse(sub_dat$controls_age_pred > thresh, 'positive', 'negative'))
    sub_dat$pred_label <- factor(sub_dat$pred_label, levels = c('positive', 'negative'))
    sub_dat$controls_age_label <- factor(sub_dat$controls_age_label, levels = c('positive', 'negative'))
    sub_dat <- sub_dat[, c('tm_donor','alpha','age_sample_collection','valid_age_pred', 'vallid_age_label')]
    result_list[[i]] <- sub_dat
    print(i)
  }
  temp <- do.call('rbind', result_list)
  return(temp)
}


temp <- get_young_labels(dat, thresh = 0.5, age = 72)

# remove duplicates
temp <- temp[!duplicated(temp$tm_donor, fromLast = TRUE),]

# creat a function that gets accuracy for each alpha level
get_acc_val <- function(temp_dat, thresh){
  all_alphas <- (1:10)/10
  result_list <- list()
  for(i in 1:length(all_alphas)){
    this_alpha <- all_alphas[i]
    sub_dat <- temp_dat[temp_dat$alpha == this_alpha,]
    sub_dat$pred_label <-as.factor(ifelse(sub_dat$valid_age_pred > thresh, 'positive', 'negative'))
    sub_dat$pred_label <- factor(sub_dat$pred_label, levels = c('positive', 'negative'))
    sub_dat$vallid_age_label <- factor(sub_dat$vallid_age_label, levels = c('positive', 'negative'))
    sub_dat$acc <- caret::confusionMatrix(sub_dat$pred_label, sub_dat$vallid_age_label)$overall[1]
    result_list[[i]] <- sub_dat
    print(i)
  }
  temp <- do.call('rbind', result_list)
  return(temp)
}
temp <- get_acc_val(dat, thresh = 0.55)

# alpha - 1.0 is best
temp <- temp[temp$alpha == 1,]

# get dataset of predictions and labels for both small and large data
pred_short <- prediction(sub$mean_preds, dat_sample$real_label)
pred_long <- prediction(sub_dat$preds, final_dat$real)

# get performace objects
perf_s <- performance(prediction.obj = pred_short, measure = 'tpr', x.measure = 'fpr')
perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')

# plot mean preds
plot(perf_s)
abline(a = 0, b =1)

# plot all preds
plot(perf_l)
abline(a = 0, b =1)


# plot optimal cutoff (TPR vs TNR)
temp_roc <- performance(pred_short, measure = 'tpr', x.measure = 'tnr')
cutoffs <- data.frame(cut=temp_roc@alpha.values[[1]], tpr=temp_roc@x.values[[1]],
                      tnr=temp_roc@y.values[[1]])

# combine plots 
plot(cutoffs$cut, cutoffs$tpr,type="l",col="red")
par(new=TRUE)
plot(cutoffs$cut, cutoffs$tnr,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")

# get precision recall curves 
temp_pr <- performance(pred_short, measure = 'prec', x.measure = 'rec')
plot(temp_pr)

# get lift curve
temp_lc <- performance(pred_short, measure="lift")
plot(temp_lc)


###############
# get optimal cutoff and confusion matrix info 
###############

# optimal cutoff
cost.perf = performance(pred_short, "cost", cost.fp = 1, cost.fn = 1)
optimal_cutoff <- pred_short@cutoffs[[1]][which.min(cost.perf@y.values[[1]])]

# get confusion matrix function for plotting 
ConfusionMatrixInfo(data = dat_sample, 
                    predict = 'mean_preds', 
                    actual = 'real_label', 
                    cutoff = optimal_cutoff,
                    get_plot = TRUE)



dev.off()
dat_sample$real_label <- as.character(dat_sample$real_label)
# get cost function by specifying cost before
cost_fp <- 1
cost_fn <- 1
roc_info <- ROCInfo(data = dat_sample, 
                    predict = 'mean_preds', 
                    actual = 'real_label', 
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
ggplot(dat_sample, 
       aes(age, mean_preds)) + 
  geom_point() +
  labs(title = 'Age in months vs risk score',
       x = 'Age (months)',
       y = 'Mean risk score') +
  geom_smooth(method = 'loess')


# beta and gender .45
# beta and no gender .45
# m and gender .49
# m and no gender .49


