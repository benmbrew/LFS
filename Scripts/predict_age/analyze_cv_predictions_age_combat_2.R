
source('helper_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


# set fixed variables
size = 'full'
model_type = 'rf'
standardize = FALSE
gender = TRUE
method = 'noob'
combat = 'combat_1'
which_methyl = 'beta'
beta_thresh = 0.05
optimal_cutoff = 0.5
remove_leading_pcs = 'first'
tech = FALSE
# standardized =FALSE

# create objects to indicate method and model details when saving
age_cutoff = 72
trained_lambda = FALSE

if(standardize){
  standardize_data <- 'standardized'
} else {
  standardize_data <- 'not_standardized'
}

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
how_many_seeds =50
how_many_folds = 5


num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds


if(model_type == 'rf'){
  
  # # read data
  final_dat <- readRDS(paste0('age_combat_data_cv/results/',combat,'_' , method, '_', size, '_',
                              num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))
  # # read data
  final_importance <- readRDS(paste0('age_combat_data_cv/results/importance_',combat,'_' , method, '_',size, '_',
                                     num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff, '.rda'))
  temp_imp <- final_importance %>% group_by(probe) %>% summarise(score = mean(score),
                                                                 counts = n())
  names(final_dat)[2] <- 'preds'
  
  
} else {
  
  # # read data
  final_dat <- readRDS(paste0('age_combat_data_cv/results/',combat,'_' , method, '_', size, '_',
                              num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',standardize_data,'.rda'))
}

###############
# get avg over seeds for model analysis
###############

# creat real label
final_dat$real <- as.factor(ifelse(final_dat$age_diagnosis <= 72, 'positive', 'negative'))
final_dat$real <- factor(final_dat$real, levels = c('positive', 'negative'))

roc_curve = roc(real ~ preds , data = final_dat)

library(DescTools)
# coords(roc=roc_curve, x = "local maximas", ret='threshold')
thresh <- as.data.frame(t(coords(roc=roc_curve, x = "all")))
thresh <- thresh[order(thresh$threshold, decreasing = FALSE),]
optimal_thresh <- thresh$threshold[which(thresh$sensitivity == Closest(thresh$sensitivity, thresh$specificity))][1]
saveRDS(optimal_thresh, paste0('pc_data_cv/results/optimal_cutoff_', combat,'_' , method, '_', size, '_',
                               num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',standardize_data,'.rda'))
basic_thresh = 0.5
# create pred class wit optimalthresholl
final_dat$pred_class <- as.factor(ifelse(final_dat$preds > basic_thresh, 'positive', 'negative'))
final_dat$pred_class <- factor(final_dat$pred_class, c('positive', 'negative'))
final_dat$real <- factor(final_dat$real, c('positive', 'negative'))
stopifnot(levels(final_dat$pred_class) == levels(final_dat$real))

if(model_type == 'rf'){
  temp <- final_dat %>% 
    group_by(tm_donor) %>% 
    summarise(mean_preds = mean(preds),
              age_onset = mean(age_diagnosis))
  
  temp$pred_class <- as.factor(ifelse(temp$mean_preds > .5, 'positive', 'negative'))
  temp$real <- as.factor(ifelse(temp$age_onset <= 72, 'positive', 'negative'))
  
  temp$real <- factor(temp$real, c('positive', 'negative'))
  temp$pred_class<- factor(temp$pred_class, c('positive', 'negative'))
  temp$acc <- caret::confusionMatrix(table(temp$pred_class, temp$real))$overall[[1]]
}

###############
# get mean predictions for each sample
###############
# best_lambda_round <- dat$lambda_value_round[dat$mean_acc == max(dat$mean_acc)]
# 
# temp_alpha <- final_dat[final_dat$alpha == best_alpha,]
# temp_alpha <- temp_alpha[temp_alpha$lambda_value_round == best_lambda_round,]

# create an indicator for each group 
final_dat$lambda_round <- round(final_dat$lambda_value, 2)
final_dat$alpha_round <- round(final_dat$alpha, 2)

best_model <- final_dat %>% 
  group_by(alpha_round, seed) %>% 
  summarise(mean_acc = caret::confusionMatrix(table(pred_class, real))$overall[1],
            mean_lambda = round(mean(lambda_value), 2),
            tot_probes = mean(tot_probes),
            non_zero = mean(non_zero),
            counts = n())

best_model <- best_model %>% 
  group_by(alpha_round) %>% 
  summarise(mean_acc = mean(mean_acc),
            mean_lambda = mean(mean_lambda),
            tot_probes = mean(tot_probes),
            non_zero = mean(non_zero),
            counts = n())

best_lambda <- best_model$mean_lambda[best_model$mean_acc == max(best_model$mean_acc)]
best_alpha <- best_model$alpha_round[best_model$mean_acc == max(best_model$mean_acc)]
model_params <- c(mean_lambda = best_lambda, mean_alpha = best_alpha)
saveRDS(model_params, paste0('pc_data_cv/results/model_params_',combat,'_' , method, '_', size, '_',
                             num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',standardize_data,'.rda'))

dat_sample <- final_dat[final_dat$alpha == best_alpha,]
dat_sample <- dat_sample[dat_sample$lambda_round == round(best_lambda,2 ),]
dat_sample$acc <- caret::confusionMatrix(table(dat_sample$pred_class, dat_sample$real))$overall[1]


# # t group by tm donor and get avg preds with errorbars
# dat <- final_dat %>% group_by(tm_donor) %>% summarise(positives = sum(real == 'positive'),
#                                                              negatives = sum(real == 'negative'),
#                                                              mean_preds = mean(preds, na.rm = TRUE),
#                                                              positives = sum(pred_class == 'positive'),
#                                                              mean_alpha = mean(alpha),
#                                                              mean_lambda = mean(lambda),
#                                                              mean_tot_probes = mean(tot_probes),
#                                                              mean_non_zero = mean(non_zero),
#                                                              age_onset = mean(age_diagnosis, na.rm = TRUE),
#                                                              age = mean(age_sample_collection, na.rm = TRUE),
#                                                              n = n())
# 
# dat$real_label <- as.factor(ifelse(dat$negatives == 0, 'positive', 'negative'))
# dat$pred_class <- as.factor(ifelse(dat$mean_preds > optimal_thresh, 'positive', 'negative'))
# dat$pred_class <- factor(dat$pred_class, c('positive', 'negative'))
# dat$real_label <- factor(dat$real_label, c('positive', 'negative'))
# stopifnot(levels(dat$pred_class) == levels(dat$real_label))
# 
# dat$acc <- caret::confusionMatrix(table(dat$pred_class, dat$real_label))$overall[1]
# 
# # get official label
# dat$real_label <- as.factor(ifelse(dat$age_onset <= 72, 'positive', 'negative'))
# dat$real_label <- factor(dat$real_label, levels = c('positive', 'negative'))
# 
# 
# 
# best_lambda <- dat$mean_lambda[dat$acc == max(dat$acc)]
# best_alpha <- dat$alpha[dat$mean_acc == max(dat$mean_acc)]
# model_params <- c(mean_lambda = best_lambda, mean_alpha = best_alpha)
# saveRDS(model_params, paste0('pc_data_cv/results/model_params_',combat,'_' , method, '_', size, '_',
#                              num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',standardize_data,'.rda'))

##############
# get performance scores and plots
###############

# get dataset of predictions and labels for both small and large data
pred_short <- prediction(dat_sample$preds, dat_sample$real)
# pred_long <- prediction(sub_dat$preds, final_dat$real)

# get performace objects
perf_s <- performance(prediction.obj = pred_short, measure = 'tpr', x.measure = 'fpr')
# perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')

# plot mean preds
plot(perf_s)
abline(a = 0, b =1)
temp <- roc(dat_sample$real, dat_sample$preds)[[9]]
legend(x = 0.8, y = 0.2, legend = temp[[1]])


# # plot all preds
# plot(perf_l)
# abline(a = 0, b =1)


# plot optimal cutoff (TPR vs TNR)
temp_roc <- performance(pred_short, measure = 'tpr', x.measure = 'tnr')
cutoffs <- data.frame(cut=temp_roc@alpha.values[[1]], tpr=temp_roc@x.values[[1]],
                      tnr=temp_roc@y.values[[1]])

# combine plots
plot(cutoffs$cut, cutoffs$tpr,type="l",col="red", cex = 2,  xlab="Cutoff",ylab="Accuracy")
par(new=TRUE)
plot(cutoffs$cut, cutoffs$tnr,type="l", lub = 2,cex = 2,col="blue",xaxt="n",yaxt="n",xlab="",ylab="")



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
                    cutoff = .46,
                    other_title = 'cv',
                    get_plot = TRUE)



dat_sample$real_label <- as.character(dat_sample$real_label)
# get cost function by specifying cost before
cost_fp <- 1
cost_fn <- 1
roc_info <- ROCInfo(data = dat_sample, 
                    predict = 'mean_preds', 
                    actual = 'real_label', 
                    other_title = 'cv',
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


# save mean lambda  and mean alpha for highest accuracy


# beta and gender .45
# beta and no gender .45
# m and gender .49
# m and no gender .49


