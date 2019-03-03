
source('helper_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


# set fixed variables
size = 'used_bh'
model_type = 'rf'
gender = FALSE
method = 'quan'
combat = 'normal'
which_methyl = 'beta'
standardize = FALSE
beta_thresh = 0.05
alpha_val = 0.1
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
how_many_seeds =10
how_many_folds = 5


num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds


if(model_type == 'rf'){
  
  # # read data
  final_dat <- readRDS(paste0('final_age_combat_cv/new_results/', combat,'_' , method, '_',size, '_',
                              num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh,'.rda'))
  
  # # read data
  final_dat_con <- readRDS(paste0('final_age_combat_cv/new_results/con_450_', combat,'_' , method, '_',size, '_',
                              num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh,'.rda'))
  # # read data
  # # read data
  final_importance <- readRDS(paste0('final_age_combat_cv/new_results/importance_', combat,'_' , method, '_',size, '_',
                            num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh,'.rda'))
  
  library(DescTools)
  # get optimal threshold 
  roc_curve = roc(real ~ positive , data = final_dat)
  
  # coords(roc=roc_curve, x = "local maximas", ret='threshold')
  thresh <- as.data.frame(t(coords(roc=roc_curve, x = "all")))
  thresh <- thresh[order(thresh$threshold, decreasing = FALSE),]
  optimal_thresh <- thresh$threshold[which(thresh$sensitivity == Closest(thresh$sensitivity, thresh$specificity))][1]
  final_dat$optimal_thresh <- optimal_thresh
  
  final_dat <- final_dat %>% group_by(tm_donor) %>% summarise(preds = mean(positive),
                                                              mean_acc = mean(accuracy),
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
  
  
  final_dat_con_mut <- final_dat_con %>% 
    filter(p53_germline == 'MUT') %>%
    group_by(tm_donor) %>% 
    summarise(preds = mean(positive),
              age_sample_collection = mean(age_sample_collection),
              optimal_thresh = mean(optimal_thresh),
              counts = n())
  final_dat_con_mut$real <- as.factor(ifelse(final_dat_con_mut$age_sample_collection <= 72, 'positive', 'negative'))
  final_dat_con_mut$real <- factor(final_dat_con_mut$real, levels = c('positive', 'negative'))
  final_dat_con_mut$pred_class <- as.factor(ifelse(final_dat_con_mut$preds > .5, 'positive', 'negative'))
  final_dat_con_mut$pred_class<- factor(final_dat_con_mut$pred_class, c('positive', 'negative'))  
  
  ConfusionMatrixInfo(data = final_dat_con_mut, 
                      predict = 'preds', 
                      actual = 'real', 
                      cutoff = optimal_thresh,
                      other_title = paste0('optimal_thresh_null_450_mut', '_',acc_opt,'_',combat,'_' , method, '_', size, '_',
                                           num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh),
                      get_plot = TRUE)
  
  # get confusion matrix function for plotting 
  ConfusionMatrixInfo(data = final_dat_con_mut, 
                      predict = 'preds', 
                      actual = 'real', 
                      cutoff = 0.5,
                      other_title = paste0('normal_thresh_null_450_mut','_', acc_normal,'_',combat,'_' , method, '_', size, '_',
                                           num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh),
                      get_plot = TRUE)
  
  
  
  
  final_dat_con_wt <- final_dat_con %>% 
    filter(p53_germline == 'WT') %>%
    group_by(tm_donor) %>% 
    summarise(preds = mean(positive),
              age_sample_collection = mean(age_sample_collection),
              optimal_thresh = mean(optimal_thresh),
              counts = n())
  final_dat_con_wt$real <- as.factor(ifelse(final_dat_con_wt$age_sample_collection <= 72, 'positive', 'negative'))
  final_dat_con_wt$real <- factor(final_dat_con_wt$real, levels = c('positive', 'negative'))
  final_dat_con_wt$pred_class <- as.factor(ifelse(final_dat_con_wt$preds > .5, 'positive', 'negative'))
  final_dat_con_wt$pred_class<- factor(final_dat_con_wt$pred_class, c('positive', 'negative'))  
  
  ConfusionMatrixInfo(data = final_dat_con_wt, 
                      predict = 'preds', 
                      actual = 'real', 
                      cutoff = optimal_thresh,
                      other_title = paste0('optimal_thresh_null_450_wt', '_',acc_opt,'_',combat,'_' , method, '_', size, '_',
                                           num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh),
                      get_plot = TRUE)
  
  # get confusion matrix function for plotting 
  ConfusionMatrixInfo(data = final_dat_con_wt, 
                      predict = 'preds', 
                      actual = 'real', 
                      cutoff = 0.5,
                      other_title = paste0('normal_thresh_null_450_wt','_', acc_normal,'_',combat,'_' , method, '_', size, '_',
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
  saveRDS(unique(final_dat$optimal_thresh), paste0('final_age_combat_cv/newer_results/optimal_cutoff_', combat,'_' , method, '_',size, '_',
                                                   num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh,'.rda'))
  
} 





if(model_type == 'enet') {
  
  # # read data
  final_dat <- readRDS(paste0('final_age_combat_cv/new_results/',combat,'_' ,method, '_', size, '_',
                              num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_', alpha_val,'_',beta_thresh,'.rda'))
  
  # get optimal lambda 
  optimal_lambda <- as.data.frame(cbind(final_dat$lambda[final_dat$accuracy == max(final_dat$accuracy)], 
                                        final_dat$non_zero[final_dat$accuracy == max(final_dat$accuracy)]))
  optimal_lambda <- optimal_lambda$V1[optimal_lambda$V2 == max(optimal_lambda$V2)]
  saveRDS(unique(optimal_lambda),paste0('final_age_combat_cv/optimal_lambda',combat,'_' , method, '_', size, '_',
                                        num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
  
  
  # group by tm donor 
  final_dat <- final_dat %>% 
    group_by(tm_donor) %>% 
    summarise(preds = mean(preds),
              mean_acc = mean(accuracy),
              mean_lambda = mean(lambda_value),
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
  
  mean_lambda = round(mean(final_dat$mean_lambda), 4)
  saveRDS(unique(final_dat$optimal_thresh), paste0('final_age_combat_cv/optimal_cutoff_', combat,'_' , method, '_', size, '_',
                                                   num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
  
  
  saveRDS(unique(mean_lambda),paste0('final_age_combat_cv/mean_lambda',combat,'_' , method, '_', size, '_',
                                     num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
  
  
  # get confusion matrix function for plotting 
  ConfusionMatrixInfo(data = final_dat, 
                      predict = 'preds', 
                      actual = 'real', 
                      cutoff = optimal_thresh,
                      other_title = paste0('optimal_thresh_cv', combat,'_' , method, '_', size, '_',
                                           num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh),
                      get_plot = TRUE)
  
  # get confusion matrix function for plotting 
  ConfusionMatrixInfo(data = final_dat, 
                      predict = 'preds', 
                      actual = 'real', 
                      cutoff = 0.5,
                      other_title = paste0('normal_thresh',combat,'_' , method, '_', size, '_',
                                           num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh),
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


