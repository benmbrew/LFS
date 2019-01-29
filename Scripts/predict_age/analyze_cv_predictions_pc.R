
source('helper_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


# next do swan with standardize and not
# set fixed variables
size = 'used_bh'
model_type = 'rf'
gender = FALSE
method = 'noob'
combat = 'normal'
standardize = TRUE
which_methyl = 'beta'
beta_thresh = 0.05
optimal_cutoff = 0.5
tech = FALSE

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
how_many_seeds = 20
how_many_folds = 5


num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds


if(model_type == 'rf'){
  
  final_dat <- readRDS(paste0('pc_data_cv/new_results/', combat,'_' , method, '_',size, '_',
                            num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))
  final_dat_con <- readRDS(paste0('pc_data_cv/new_results/', 'con_', combat,'_' , method, '_',size, '_',
                                num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))
  
  final_dat_val <- readRDS(paste0('pc_data_cv/new_results/', 'valid_', combat,'_' , method, '_',size, '_',
                                num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))
  
  # # read data
  final_importance <- readRDS(paste0('pc_data_cv/new_results/', 'importance_', combat,'_' ,method, '_',size, '_',
                                   num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))
  
  library(DescTools)
  # get optimal threshold 
  roc_curve = roc(real ~ positive , data = final_dat)
  
  # coords(roc=roc_curve, x = "local maximas", ret='threshold')
  thresh <- as.data.frame(t(coords(roc=roc_curve, x = "all")))
  thresh <- thresh[order(thresh$threshold, decreasing = FALSE),]
  optimal_thresh <- thresh$threshold[which(thresh$sensitivity == Closest(thresh$sensitivity, thresh$specificity))][1]
  final_dat$optimal_thresh <- optimal_thresh
  
  final_dat <- final_dat %>% group_by(tm_donor) %>% summarise(preds = mean(positive),
                                                              age_diagnosis = mean(age_diagnosis),
                                                              age_sample_collection = mean(age_sample_collection),
                                                              optimal_thresh = mean(optimal_thresh),
                                                              counts = n())
  final_dat$real <- as.factor(ifelse(final_dat$age_diagnosis <= 72, 'positive', 'negative'))
  final_dat$real <- factor(final_dat$real, levels = c('positive', 'negative'))
  final_dat$pred_class <- as.factor(ifelse(final_dat$preds > .5, 'positive', 'negative'))
  
  final_dat$pred_class<- factor(final_dat$pred_class, c('positive', 'negative'))
  final_dat$acc <- caret::confusionMatrix(table(final_dat$pred_class, final_dat$real))$overall[[1]]
  
  final_dat$pred_class_opt <- as.factor(ifelse(final_dat$preds > unique(final_dat$optimal_thresh), 'positive', 'negative'))
  
  final_dat$pred_class_opt<- factor(final_dat$pred_class_opt, c('positive', 'negative'))
  final_dat$acc_opt <- caret::confusionMatrix(table(final_dat$pred_class_opt, final_dat$real))$overall[[1]]
  
  
  saveRDS(unique(final_dat$optimal_thresh), paste0('pc_data_cv/optimal_cutoff_', combat,'_' , method, '_',size, '_',
                                                   num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh,'.rda'))
  acc_opt = round(unique(final_dat$acc_opt), 2)
  acc_normal = round(unique(final_dat$acc), 2)
  # get confusion matrix function for plotting 
  ConfusionMatrixInfo(data = final_dat, 
                      predict = 'preds', 
                      actual = 'real', 
                      cutoff = optimal_thresh,
                      other_title = paste0('optimal_thresh_cv', '_',acc_opt,'_',combat,'_' , method, '_', size, '_',
                                           num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh),
                      get_plot = TRUE)
  
  # get confusion matrix function for plotting 
  ConfusionMatrixInfo(data = final_dat, 
                      predict = 'preds', 
                      actual = 'real', 
                      cutoff = 0.5,
                      other_title = paste0('normal_thresh','_', acc_normal,'_',combat,'_' , method, '_', size, '_',
                                           num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh),
                      get_plot = TRUE)
  
  # get dataset of predictions and labels for both small and large data
  pred_short <- prediction(final_dat$preds, final_dat$real)
  # pred_long <- prediction(sub_dat$preds, final_dat$real)
  
  # get performace objects
  perf_s <- performance(prediction.obj = pred_short, measure = 'tpr', x.measure = 'fpr')
  # perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')
  
  # plot mean preds
  plot(perf_s)
  abline(a = 0, b =1)
  temp <- round(roc(final_dat$real, final_dat$preds)[[9]], 2)
  legend(x = 0.8, y = 0.2, legend = paste0('AUC = ', temp[[1]]))
  
  # plot age as a function of prediction
  ggplot(final_dat, 
         aes(age_sample_collection, preds)) + 
    geom_point() +
    labs(title = 'Age in months vs risk score',
         x = 'Age (months)',
         y = 'Mean risk score') +
    geom_smooth(method = 'loess')
  
  
} else {
  
  # # read data
  final_dat <- readRDS(paste0('pc_data_cv/',combat,'_' , method, '_', size, '_',
                              num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',standardize_data,'.rda'))
  
  # remove duplicates
  final_dat <- get_acc_val(final_dat, thresh = 0.5)
  final_dat <- final_dat[order(final_dat$acc, decreasing = TRUE),]
  
  optimal_alpha = final_dat$alpha[final_dat$acc == max(final_dat$acc)]
  optimal_alpha <- unique(optimal_alpha)
  
  optimal_lambda = as.data.frame(cbind(final_dat$lambda[final_dat$acc == max(final_dat$acc)],
                                       final_dat$non_zero[final_dat$acc == max(final_dat$acc)]))
  optimal_lambda <- optimal_lambda$V1[optimal_lambda$V2 == max(optimal_lambda$V2)]
  optimal_lambda <- unique(optimal_lambda)
  saveRDS(unique(optimal_alpha),paste0('pc_data_cv/optimal_alpha',combat,'_' , method, '_', size, '_',
                                        num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
  
  sub_dat = final_dat[final_dat$alpha == optimal_alpha,]
  # get optmal lamda 
  sub_dat <- final_dat %>% 
    group_by(tm_donor) %>% 
    summarise(mean_preds = mean(preds),
              mean_lambda = mean(lambda),
              mean_acc = mean(acc),
              age_onset = mean(age_diagnosis),
              counts = n())
  
  sub_dat$real <- as.factor(ifelse(sub_dat$age_onset <= 72, 'positive', 'negative'))
  sub_dat$real <- factor(sub_dat$real, levels = c('positive', 'negative'))
  
  # coords(roc=roc_curve, x = "local maximas", ret='threshold')
  roc_curve = roc(real ~ mean_preds , data = sub_dat)
  thresh <- as.data.frame(t(coords(roc=roc_curve, x = "all")))
  thresh <- thresh[order(thresh$threshold, decreasing = FALSE),]
  optimal_thresh <- thresh$threshold[which(thresh$sensitivity == Closest(thresh$sensitivity, thresh$specificity))][1]
  sub_dat$optimal_thresh <- optimal_thresh
  sub_dat$pred_class <- ifelse(sub_dat$mean_preds > optimal_thresh, 'positive', 'negative')
  sub_dat$pred_class<- factor(sub_dat$pred_class, c('positive', 'negative'))
  sub_dat$acc <- caret::confusionMatrix(table(sub_dat$pred_class, sub_dat$real))$overall[[1]]
  
  saveRDS(unique(final_dat$optimal_thresh), paste0('pc_data_cv/optimal_cutoff_', combat,'_' , method, '_', size, '_',
                                                   num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
  
  
  saveRDS(unique(optimal_lambda),paste0('pc_data_cv/optimal_lambda',combat,'_' , method, '_', size, '_',
                                     num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
  
  # group by tm donor 
  final_dat <- final_dat %>% 
    group_by(tm_donor) %>% 
    summarise(preds = mean(preds),
              mean_acc = mean(accuracy),
              mean_lambda = mean(lambda_value),
              age_diagnosis = mean(age_diagnosis),
              age_sample_collection = mean(age_sample_collection),
              counts = n())
  
  # get mean lambda 

  final_dat$real <- as.factor(ifelse(final_dat$age_diagnosis <= 72, 'positive', 'negative'))
  final_dat$real <- factor(final_dat$real, levels = c('positive', 'negative'))
  final_dat$pred_class <- as.factor(ifelse(final_dat$preds > .5, 'positive', 'negative'))
  roc_curve = roc(real ~ preds , data = final_dat)
  
  # coords(roc=roc_curve, x = "local maximas", ret='threshold')
  thresh <- as.data.frame(t(coords(roc=roc_curve, x = "all")))
  thresh <- thresh[order(thresh$threshold, decreasing = FALSE),]
  mean_thresh <- thresh$threshold[which(thresh$sensitivity == Closest(thresh$sensitivity, thresh$specificity))][1]
  final_dat$mean_thresh <- mean_thresh
  
  final_dat$pred_class<- factor(final_dat$pred_class, c('positive', 'negative'))
  final_dat$acc <- caret::confusionMatrix(table(final_dat$pred_class, final_dat$real))$overall[[1]]
  
  
  
  final_dat$pred_class_opt <- as.factor(ifelse(final_dat$preds > unique(final_dat$mean_thresh), 'positive', 'negative'))
  
  final_dat$pred_class_opt<- factor(final_dat$pred_class_opt, c('positive', 'negative'))
  final_dat$acc_opt <- caret::confusionMatrix(table(final_dat$pred_class_opt, final_dat$real))$overall[[1]]
  
  mean_lambda = round(mean(final_dat$mean_lambda), 4)
  saveRDS(unique(final_dat$mean_thresh), paste0('pc_data_cv/mean_cutoff_', combat,'_' , method, '_', size, '_',
                                                   num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
  
  
  saveRDS(unique(mean_lambda),paste0('pc_data_cv/mean_lambda',combat,'_' , method, '_', size, '_',
                                     num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
  
  
  
  acc_opt = round(unique(final_dat$acc_opt), 2)
  acc_normal = round(unique(final_dat$acc), 2)
  # get confusion matrix function for plotting 
  ConfusionMatrixInfo(data = final_dat, 
                      predict = 'preds', 
                      actual = 'real', 
                      cutoff = optimal_thresh,
                      other_title = paste0('optimal_thresh_cv', '_',acc_opt,'_',combat,'_' , method, '_', size, '_',
                                           num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh, '_', standardize_data),
                      get_plot = TRUE)
  
  # get confusion matrix function for plotting 
  ConfusionMatrixInfo(data = final_dat, 
                      predict = 'preds', 
                      actual = 'real', 
                      cutoff = 0.5,
                      other_title = paste0('normal_thresh','_',acc_normal,'_',combat,'_' , method, '_', size, '_',
                                           num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh, '_', standardize_data),
                      get_plot = TRUE)
  
  # get dataset of predictions and labels for both small and large data
  pred_short <- prediction(final_dat$preds, final_dat$real)
  # pred_long <- prediction(sub_dat$preds, final_dat$real)
  
  # get performace objects
  perf_s <- performance(prediction.obj = pred_short, measure = 'tpr', x.measure = 'fpr')
  # perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')
  
  # plot mean preds
  plot(perf_s)
  abline(a = 0, b =1)
  temp <- round(roc(final_dat$real, final_dat$preds)[[9]], 2)
  legend(x = 0.8, y = 0.2, legend = paste0('AUC = ', temp[[1]]))
  
  # plot age as a function of prediction
  ggplot(final_dat, 
         aes(age_sample_collection, preds)) + 
    geom_point() +
    labs(title = 'Age in months vs risk score',
         x = 'Age (months)',
         y = 'Mean risk score') +
    geom_smooth(method = 'loess')
  
  
}

###############
# get avg over seeds for model analysis
###############





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
saveRDS(model_params, paste0('pc_data_cv/model_params_',combat,'_' , method, '_', size, '_',
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



# save mean lambda  and mean alpha for highest accuracy


# beta and gender .45
# beta and no gender .45
# m and gender .49
# m and no gender .49


