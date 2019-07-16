#source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  
data_type = 'beta'
combat = FALSE

if(combat){
  used_combat <- 'used_combat'
} else {
  used_combat <- 'no_combat'
}

# create objects to indicate method and model details when saving
gender = FALSE
tech = FALSE
how_many_seeds = 20
how_many_folds = 5
age_cutoff = 72
model_type = 'enet'
trained_lambda = FALSE

if(trained_lambda == TRUE){
  is_lambda <- 'no_cv_lambda'
} else {
  is_lambda <- 'used_cv_lambda'
}


if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}

num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds


if(model_type == 'rf'){
  temp_con <- readRDS(paste0('residual_data_test/', 'con_test_resid',data_type, '_', is_gen, '_',used_combat,'_', model_type,'.rda'))
  temp_valid <- readRDS(paste0('residual_data_test/', 'valid_test_resid',data_type, '_', is_gen,'_',used_combat, '_', model_type,'.rda'))
  temp_importance <- readRDS( paste0('residual_data_test/', 'importance_resid',data_type, '_', is_gen,'_',used_combat, '_', model_type,'.rda'))
  temp_con <- temp_con[, c('tm_donor', 'positive', 'tech','p53_germline','age_sample_collection', 'tot_probes', 'accuracy')]
  temp_valid <- temp_valid[, c('tm_donor', 'positive', 'tech','p53_germline','age_sample_collection', 'tot_probes', 'accuracy')]
  
} else {
  temp_con <- readRDS(paste0('residual_data_test/', 'con_test_resid',data_type, '_', is_gen, '_',is_lambda,'_',used_combat,'_', model_type,'.rda'))
  temp_valid <- readRDS(paste0('residual_data_test/', 'valid_test_resid',data_type, '_', is_gen,'_',is_lambda,'_',used_combat, '_', model_type,'.rda'))
  temp_con <- temp_con[, c('tm_donor', 'p53_germline', 'tech','controls_age_pred', 'controls_age_label', 'age_sample_collection', 'non_zero', 'alpha')]
  temp_valid <- temp_valid[, c('tm_donor', 'p53_germline','tech','valid_age_pred', 'valid_age_label', 'age_sample_collection', 'non_zero', 'alpha')]
}

# subset data
temp_group <- temp_con %>% group_by(tm_donor, alpha) %>% 
  filter(p53_germline == 'MUT') %>%
  mutate(mean_pred = mean(controls_age_pred, na.rm = TRUE))
temp_1 = temp_group[temp_group$alpha == 0.1,]
temp_2 = temp_group[temp_group$alpha == 0.2,]
temp_3 = temp_group[temp_group$alpha == 0.3,]
temp_9 = temp_group[temp_group$alpha == 0.9,]

temp_3 <- temp_3[!duplicated(temp_3$tm_donor),]
write.csv(temp_1, '~/Desktop/strategy_3_residual/con_table.csv')

# remove duplicates
temp <- get_acc_val(temp_valid, thresh = 0.5)
temp <- temp[order(temp$acc, decreasing = TRUE),]
temp_dedup <- temp[!duplicated(temp$tm_donor),]
max_alpha = temp_dedup$alpha[which(temp_dedup$acc == max(temp_dedup$acc))]

temp <- temp[temp$alpha == 0.9,]
temp <- temp[!duplicated(temp$tm_donor),]
write.csv(temp, '~/Desktop/strategy_4_pc/valid_table.csv')

# get dataset of predictions and labels for both small and large data
pred_short <- prediction(temp$valid_age_pred, temp$valid_age_label)
# pred_long <- prediction($preds, final_dat$real)

# get performace objects
perf_s <- performance(prediction.obj = pred_short, measure = 'tpr', x.measure = 'fpr')
# perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')

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
ConfusionMatrixInfo(data = temp, 
                    predict = 'valid_age_pred', 
                    actual = 'valid_age_label', 
                    cutoff = optimal_cutoff,
                    get_plot = TRUE)



dev.off()
dat_sample$real_label <- as.character(dat_sample$real_label)
# get cost function by specifying cost before
cost_fp <- 1
cost_fn <- 1
roc_info <- ROCInfo(data = temp, 
                    predict = 'valid_age_pred', 
                    actual = 'valid_age_label', 
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


t_dat <- dat[dat$alpha == 0.5,]
# get confusion matrix function for plotting 
t <-ConfusionMatrixInfo(data = t_dat, 
                        predict = 'valid_age_pred', 
                        actual = 'valid_age_label', 
                        cutoff = .6,
                        get_plot = FALSE)

val_dat <- as.data.frame(cbind(t_dat, t))
val_dat$valid_age_pred <- val_dat$vallid_age_label <- NULL



