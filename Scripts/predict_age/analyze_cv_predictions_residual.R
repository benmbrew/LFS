source('helper_functions.R')
#
# create fixed objects to model and pipeline inputs and saving  
data_type = 'beta'
combat = TRUE

if(combat){
  used_combat <- 'used_combat'
} else {
  used_combat <- 'no_combat'
}
if(data_type == 'beta'){
  beta_thresh = 0.0001
} else {
  beta_thresh = 0.001
}

#create objects to indicate method and model details when saving
gender = FALSE
tech = FALSE
how_many_seeds = 20
how_many_folds = 5
age_cutoff = 72
model_type = 'enet'

if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}

num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds

# save data
final_dat <- readRDS(paste0('residual_data_cv/results/', data_type, '_',
                          num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',used_combat,'.rda'))


# subset data
sub_dat <- final_dat[, c('tm_donor','cancer_diagnosis_diagnoses','preds', 'pred_class', 'real', 'accuracy', 'lambda', 
                         'alpha', 'tot_probes', 'fold', 'seed', 'non_zero')]

# plot alpha against accuracy
plot_acc(sub_dat, acc_column = 'accuracy', column = 'alpha', bar = FALSE)

# plot lambda against accuracy
plot_acc(sub_dat, acc_column = 'accuracy', column = 'lambda', bar = FALSE)

# plot cancer_diagnosis against accuracy
sub_cancer <- sub_dat %>% group_by(cancer_diagnosis_diagnoses) %>% summarise(mean_acc = mean(accuracy))
plot_acc(sub_cancer, acc_column = 'mean_acc', column = 'cancer_diagnosis_diagnoses', bar = TRUE)

# plot tot_probes against accuracy
plot_acc(sub_dat, acc_column = 'accuracy', column = 'tot_probes', bar = FALSE)

# plot non zero against accuracy
plot_acc(sub_dat, acc_column = 'accuracy', column = 'non_zero', bar = FALSE)

# plot non zero against accuracy
plot_acc(sub_dat,acc_column = 'accuracy', column = 'fold', bar = FALSE)

plot_3d_model(sub_dat, type = 'h')


###############
# get avg over seeds for model analysis
###############

# get the average over seeds
final_dat$seed <- as.character(final_dat$seed)
dat <- final_dat %>% group_by(seed) %>% summarise(mean_acc = mean(accuracy, na.rm = TRUE),
                                                  mean_lambda = mean(lambda, na.rm = TRUE),
                                                  mean_alpha = mean(alpha, na.rm = TRUE),
                                                  sd_acc = sd(accuracy, na.rm = TRUE),
                                                  sd_lambda = sd(lambda, na.rm = TRUE),
                                                  sd_alpha= sd(alpha, na.rm = TRUE),
                                                  n = n()) %>% 
  mutate(se_acc = sd_acc/sqrt(n),
         acc_lower_ci = mean_acc - qt(1 - (0.05 / 2), n - 1)* se_acc,
         acc_upper_ci = mean_acc + qt(1 - (0.05 / 2), n - 1)* se_acc,
         se_lambda = sd_lambda/sqrt(n),
         lambda_lower_ci = mean_lambda - qt(1 - (0.05 / 2), n - 1)* se_lambda,
         lambda_upper_ci = mean_lambda + qt(1 - (0.05 / 2), n - 1)* se_lambda,
         se_alpha = sd_alpha/sqrt(n),
         alpha_lower_ci = mean_alpha - qt(1 - (0.05 / 2), n - 1)* se_alpha,
         alpha_upper_ci = mean_alpha + qt(1 - (0.05 / 2), n - 1)* se_alpha)

# remove sd_ variables 
dat$sd_acc <- dat$sd_alpha <- dat$se_lambda <- NULL


# plot alpha against accuracy
plot_acc(dat, acc_column = 'mean_acc', column = 'mean_alpha', bar = FALSE)

# plot lambda against accuracy
plot_acc(dat, acc_column = 'mean_acc', column = 'mean_lambda', bar = FALSE)

# 3d scatter plot with mean values 
plot_3d_model_means(dat)

# save mean lambda  and mean alpha for highest accuracy
mean_lambda <- mean(dat$mean_lambda[dat$mean_acc == max(dat$mean_acc)])
mean_alpha <- mean(dat$mean_alpha[dat$mean_acc == max(dat$mean_acc)])
model_params <- c(mean_lambda = mean_lambda, mean_alpha = mean_alpha)
saveRDS(model_params, paste0('residual_data_cv/results/model_params_', data_type, '_',
                             num_seeds, '_', k_folds, '_', is_gen, '_', model_type,'.rda'))
####

###############
# get mean predictions for each sample
###############

# first group by tm donor and get avg preds with errorbars
dat_sample <- final_dat %>% group_by(tm_donor) %>% summarise(sum_positive = sum(real == 'positive', na.rm = TRUE),
                                                             sum_negative = sum(real == 'negative', na.rm = TRUE),
                                                             mean_preds = mean(preds, na.rm = TRUE),
                                                             sd_preds= sd(preds, na.rm = TRUE),
                                                             age = mean(age_sample_collection, na.rm = TRUE),
                                                             n = n()) %>% 
  mutate(se_preds = sd_preds/sqrt(n),
         lower_ci = mean_preds - qt(1 - (0.05 / 2), n - 1)* se_preds,
         upper_ci = mean_preds + qt(1 - (0.05 / 2), n - 1)* se_preds)

# get official label
dat_sample$real_label <- as.factor(ifelse(dat_sample$sum_positive == 20, 'positive', 'negative'))
dat_sample$sum_negative <- dat_sample$sum_positive <- NULL
dat_sample$real_label <- factor(dat_sample$real_label, levels = c('positive', 'negative'))

###############
# get performance scores and plots
###############

# get dataset of predictions and labels for both small and large data
pred_short <- prediction(dat_sample$mean_preds, dat_sample$real_label)
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
dev.off()
plot(cutoffs$cut, cutoffs$tpr,type="l",col="red", cex = 2,  xlab="Cutoff",ylab="Accuracy")
par(new=TRUE)
plot(cutoffs$cut, cutoffs$tnr,type="l", lub = 2,cex = 2,col="blue",xaxt="n",yaxt="n",xlab="",ylab="")

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
saveRDS(optimal_cutoff, paste0('residual_data_cv/results/optimal_cutoff_', data_type, '_',
                               num_seeds, '_', k_folds, '_', is_gen, '_', model_type,'.rda'))
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



