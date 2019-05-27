source('all_functions.R')


# next do swan with standardize and not
# set fixed variables
# set fixed variables
size = 'full'
model_type = 'rf'
null_450 = TRUE
null_450_all = FALSE
use_p53 = FALSE
gender = TRUE
use_cancer = TRUE
method = 'noob'
include_under_6 = FALSE
combat = 'normal'
alpha_val = 0.1
train_alpha = FALSE
train_lambda = FALSE
train_cutoff = TRUE
mean_lambda = FALSE
which_methyl = 'beta'
beta_thresh = 0.1
age_cutoff =  72
tech = FALSE
standardize = FALSE
boot_num = 50


if(null_450 & !null_450_all){
  use_null_450 <- 'used_null_450_mut'
} else  if(!null_450 & null_450_all){
  use_null_450 <- 'used_null_450_all'
} else {
  use_null_450 <- 'no_null_450'
  
}

if(include_under_6){
  used_under_6 <- 'used_under_6'
} else {
  used_under_6 <- 'no_under_6'
}
if(use_cancer){
  removed_cancer <- 'removed_cancer'
} else {
  removed_cancer <- 'no_removed_cancer'
}

if(use_p53){
  used_p53 <- 'used_p53'
} else {
  used_p53 <- 'no_p53'
}

if(use_cancer){
  removed_cancer <- 'removed_cancer'
} else {
  removed_cancer <- 'no_removed_cancer'
}
if(standardize){
  standardize_data <- 'standardized'
} else {
  standardize_data <- 'not_standardized'
}

if(train_lambda){
  is_lambda <- 'lambda_test'
} else {
  is_lambda <- 'lambda_train'
}

if(train_alpha){
  is_alpha <- 'alpha_test'
} else {
  is_alpha <- 'alpha_train'
}


if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}
how_many_seeds = 50
how_many_folds = 5


num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds

# create range for random sampling for cross validation
seed_range <- c(1:how_many_seeds)

# create list to store model
all_test_results <- list()
importance_results <- list()
rf_pred_results <- list()
rf_important_results <- list()


if(size =='used_bh'){
  cases_450 <-  readRDS(paste0('../../Data/', method,'/cases_450_small_cv', combat,'.rda'))
  cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_small_cv', combat,'.rda'))
  con_mut <- readRDS(paste0('../../Data/', method,'/con_mut_small_cv', combat,'.rda'))
  con_850 <- readRDS( paste0('../../Data/', method,'/con_850_small_cv', combat,'.rda'))
  con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_small_cv', combat,'.rda'))
  
} else {
  cases_450 <-  readRDS(paste0('../../Data/', method,'/cases_450_cv', combat,'.rda'))
  cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_cv', combat,'.rda'))
  con_mut <- readRDS(paste0('../../Data/', method,'/con_mut_cv', combat,'.rda'))
  con_850 <- readRDS( paste0('../../Data/', method,'/con_850_cv', combat,'.rda'))
  con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_cv', combat,'.rda'))
  
  
  # # randomly sample from all cgs
  # clin_names <- names(cases_450)[1:12]
  # r_cgs <- sample(names(cases_450)[13:ncol(cases_450)], 5000)
  # cases_450 <- cases_450[c(clin_names, r_cgs)]
  # con_wt <- con_wt[c(clin_names, r_cgs)]
  # con_mut <- con_mut[c(clin_names, r_cgs)]
  # con_850 <- con_850[c(clin_names, r_cgs)]
  # cases_850 <- cases_850[c(clin_names, r_cgs)]
  # 
}

optimal_thesh = 0.5
lambda_val = 0.05


if(model_type == 'enet'){
  if(train_lambda){
    
    lambda_val <- 0.05
    
  } else if (mean_lambda) {
    lambda_val <-  readRDS(paste0('pc_data_cv/mean_lambda',combat,'_' , method, '_', size, '_',
                                  num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
    
  } else {
    lambda_val <-  readRDS(paste0('pc_data_cv/optimal_lambda',combat,'_' , method, '_', size, '_',
                                  num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
  }
} 

if(model_type == 'enet'){
  
  if(train_cutoff){
    optimal_thresh <- readRDS(paste0('pc_data_cv/optimal_cutoff_', combat,'_' , method, '_', size, '_',
                                     num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
    
    # readRDS(paste0('pc_data_cv/mean_cutoff_', combat,'_' , method, '_', size, '_',
    #                num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
    # 
    
  } else {
    optimal_thresh = 0.5
  }
  
  if(train_alpha){
    optimal_alpha <- readRDS(paste0('pc_data_cv/optimal_alpha',combat,'_' , method, '_', size, '_',
                                    num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
    
  } else {
    optimal_alpha = alpha_val
  }
  
} else {
  
  
  optimal_thresh = 0.5

}


##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
##########
g_ranges <- readRDS('../../Data/g_ranges.rda')

# get probes from rownames
g_ranges$probe <- rownames(g_ranges)

# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]

names(g_ranges)[1] <- 'chr'

con_wt$tech <- '450k'
if(use_cancer){
  # remove cancer signature
  bh_feats <- bump_hunter(dat_1 = cases_450, 
                          dat_2 = con_mut , 
                          bump = 'cancer', 
                          boot_num = boot_num, 
                          beta_thresh = beta_thresh,
                          methyl_type = methyl_type,
                          g_ranges = g_ranges)
  
  # combine call controls 
  con_all<- rbind(con_850, 
                  con_mut,
                  con_wt)
  rm(con_mut)
  rm(con_wt)
  rm(con_850)
  
  # get intersect_names
  intersect_names <- names(cases_450)[grepl('^cg', names(cases_450))]
  
  # get feature list
  colnames(bh_feats)[1] <- 'chr'
  remove_features <- inner_join(bh_feats, g_ranges)$probe
  
  # take remove features out of colnames 
  bh_features <- intersect_names[!intersect_names %in% remove_features]
  
  # subset all data by bh_features
  cases_450 <- remove_cancer_feats(cases_450, bh_feats = bh_features)
  con_all <- remove_cancer_feats(con_all, bh_feats = bh_features)
  cases_850 <- remove_cancer_feats(cases_850, bh_feats = bh_features)
} else {
  
  # combine call controls 
  con_all<- rbind(con_850, 
                  con_mut,
                  con_wt)
  rm(con_mut)
  rm(con_wt)
  rm(con_850)
  
  bh_features <- names(cases_450)[grepl('^cg', names(cases_450))]
  
}
# load cases
cases_450 <- as.data.frame(cbind(WT = 0, as.data.frame(class.ind(cases_450$p53_germline)), 
                   cases_450))


# load cases
cases_850 <- as.data.frame(cbind(WT = 0, as.data.frame(class.ind(cases_850$p53_germline)), 
                                 cases_850))



# load cases
con_all <- cbind(as.data.frame(class.ind(con_all$p53_germline)), 
                                 con_all)



optimal_thresh = 0.5
lambda_val = 0.5
# get s_num and alpha_value
if(model_type == 'enet'){
  # s_num = mean_lambda
  s_num = lambda_val
  # s_num <- round(model_params[[1]], 3)
  # alpha_num <- round(model_params[[2]], 2)
  # 
  # creat list to store results for alpha
  
  alpha_num <- alpha_val
  
  
  message('working on alpha = ', alpha_num)
  result_list <- test_model_enet(cases = cases_450,
                                 controls = con_all,
                                 valid = cases_850,
                                 null_450 = use_null_450,
                                 use_6 = include_under_6,
                                 use_p53 = use_p53,
                                 age_cutoff = age_cutoff,
                                 gender = gender,
                                 tech = tech,
                                 test_lambda = train_lambda,
                                 alpha_value = alpha_num,
                                 lambda_value = s_num,
                                 control_age = FALSE,
                                 bh_features = bh_features)
  
  con_dat <- result_list[[1]]
  valid_dat <- result_list[[2]]
  mod_dat <- result_list[[3]]
  
  
  
  
  # read in cases_450
  saveRDS(con_dat, paste0('pc_data_test/', 'con_test_', method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_', is_lambda,'_',standardize_data,'_',alpha_num,'_',s_num,'_', is_alpha,'_',optimal_thresh,'_',use_null_450,'_',removed_cancer,'_',use_p53,'_', used_under_6,'.rda'))
  
  saveRDS(valid_dat, paste0('pc_data_test/', 'valid_test_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_', is_lambda,'_',standardize_data,'_',alpha_num,'_',s_num,'_', is_alpha,'_',optimal_thresh,'_',use_null_450,'_',removed_cancer,'_',use_p53,'_', used_under_6,'.rda'))
  saveRDS(mod_dat, paste0('pc_data_test/', 'model_test_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_', is_lambda,'_',standardize_data,'_',alpha_num,'_',s_num,'_',is_alpha,'_' ,optimal_thresh,'_',use_null_450,'_',removed_cancer,'_',use_p53,'_', used_under_6,'.rda'))
  
  
  
} else {
  
  
  result_list <- test_model_rf(cases = cases_450,
                               controls = con_all,
                               valid = cases_850,
                               null_450 = use_null_450,
                               use_6 = include_under_6,
                               use_p53 = use_p53,
                               age_cutoff = age_cutoff,
                               gender = gender,
                               tech = tech,
                               control_age = FALSE,
                               optimal_cutoff = optimal_thresh,
                               bh_features = bh_features)
  
  temp_valid <- result_list[[1]]
  temp_con <- result_list[[2]]
  temp_importance  <- result_list[[3]]
  temp_model <- result_list[[4]]
  
  
  saveRDS(temp_con, paste0('pc_data_test/', 'con_test_',method,'_',size,'_',is_gen,'_',combat,'_', model_type,'_',optimal_thresh,'_',use_null_450,'_',removed_cancer,'_',removed_cancer,'_',use_p53,'_', used_under_6,'.rda'))
  
  saveRDS(temp_valid, paste0('pc_data_test/', 'valid_test_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_',optimal_thresh,'_',use_null_450,'_',removed_cancer,'_', used_under_6,'.rda'))
  
  saveRDS(temp_importance, paste0('pc_data_test/', 'importance_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_',optimal_thresh,'_',use_null_450,'_',removed_cancer,'_',use_p53,'_', used_under_6,'.rda'))
  saveRDS(temp_model, paste0('pc_data_test/', 'model_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_',optimal_thresh,'_',use_null_450,'_',removed_cancer, '_',use_p53,'_', used_under_6,'.rda'))
  
}


if (model_type == 'rf'){
  
  # optimal_thresh = 0.5
  temp_valid <- temp_valid[, c('MUT','WT','tm_donor','positive','real' ,'age_diagnosis', 'age_sample_collection', 'tech')]
  temp_con <- temp_con[, c('MUT','WT','tm_donor','positive', 'real', 'age_sample_collection', 'tech')]
  temp_con$pred_class <- as.factor(ifelse(temp_con$positive > optimal_thresh, 'positive', 'negative'))
  temp_con$pred_class <- factor(temp_con$pred_class, levels = c('positive', 'negative'))
  
  # temp_group <- temp_con %>% group_by(tm_donor) %>% 
  #   filter(p53_germline == 'MUT') %>%
  #   mutate(mean_pred = mean(positive, na.rm = TRUE))
  # temp_850 <- temp_con[temp_con$tech == '850k',]
  # any(temp_850$p53_germline == 'WT')
  # if(any(temp_850$p53_germline == 'WT')){
  #   
  #   temp_450_wt <- temp_850[temp_850$p53_germline == 'WT',]
  #   temp_850 <- temp_850[temp_850$p53_germline == 'MUT',]
  #   temp_450 <- temp_con[temp_con$tech == '450k',]
  #   
  #   
  # } else {
  #   temp_450 <- temp_con[temp_con$tech == '450k',]
  #   
  #   temp_450_wt <- temp_450[temp_450$p53_germline == 'WT',]
  #   temp_450 <- temp_450[temp_450$p53_germline == 'MUT',]
  #   
  # }
  # temp_450$real_label <- as.factor(ifelse(temp_450$age_sample_collection > 72, 'negative', 'positive'))
  # temp_450_wt$real_label <- as.factor(ifelse(temp_450_wt$age_sample_collection > 72, 'negative', 'positive'))
  # temp_850$real_label <- as.factor(ifelse(temp_850$age_sample_collection > 72, 'negative', 'positive'))
  # temp_450$real_label <- factor(temp_450$real_label, levels = c('positive', 'negative'))
  # temp_850$real_label <- factor(temp_850$real_label, levels = c('positive', 'negative'))
  # temp_450_wt$real_label <- factor(temp_450_wt$real_label, levels = c('positive', 'negative'))
  # 
  # temp_valid$pred_class <- as.factor(ifelse(temp_valid$positive > optimal_thresh, 'positive', 'negative'))
  # temp_valid$pred_class <- factor(temp_valid$pred_class, levels = c('positive', 'negative'))
  # temp_valid$real <- factor(temp_valid$real, levels = c('positive', 'negative'))
  
  
  temp_valid$pred_class <- as.factor(ifelse(temp_valid$positive > optimal_thresh, 'positive', 'negative'))
  temp_valid$pred_class <- factor(temp_valid$pred_class, levels = c('positive', 'negative'))
  temp_valid$acc <- caret::confusionMatrix(table(temp_valid$pred_class, temp_valid$real))$overall[[1]]
  acc <-  round(caret::confusionMatrix(table(temp_valid$pred_class, temp_valid$real))$overall[[1]], 2)
  # get confusion matrix function for plotting
  library(data.table)
  ConfusionMatrixInfo(data = temp_valid,
                      other_title = '',
                      predict = 'positive',
                      actual = 'real',
                      cutoff = optimal_thresh,
                      get_plot = TRUE,
                      data_type = 'valid',
                      points = FALSE)
  
  ConfusionMatrixInfo(data = temp_con,
                      other_title = paste0('null all','_',combat,'_' , method, '_', size, '_',
                                           num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh, '_',optimal_thresh),
                      predict = 'positive',
                      actual = 'real',
                      cutoff = optimal_thresh,
                      get_plot = TRUE,
                      data_type = 'null',
                      points = FALSE)
  
  temp_con$acc <- caret::confusionMatrix(table(temp_con$pred_class, temp_con$real))$overall[[1]]
  acc <-  round(caret::confusionMatrix(table(temp_con$pred_class, temp_con$real))$overall[[1]], 2)
  
  t.test(x = temp_con$positive[temp_con$real == 'positive'], 
         y = temp_con$positive[temp_con$real == 'negative'])
  
  if(use_null_450 == 'used_null_450_all'){
    temp_con_wt <- temp_con[temp_con$WT == 1,]
    temp_con_mut <- temp_con[temp_con$MUT == 1,]
    
    ConfusionMatrixInfo(data = temp_con_wt,
                        other_title = '',
                        predict = 'positive',
                        actual = 'real',
                        cutoff = optimal_thresh,
                        get_plot = TRUE,
                        data_type = 'null',
                        points = FALSE)
    
    ConfusionMatrixInfo(data = temp_con_mut,
                        other_title = '',
                        predict = 'positive',
                        actual = 'real',
                        cutoff = optimal_thresh,
                        get_plot = TRUE,
                        data_type = 'null',
                        points = FALSE)
  }
  
  t1 <- temp_con_mut[temp_con_mut$age_sample_collection < 72,,]
  t2 <- temp_con_mut[temp_con_mut$positive  > 0.5,]
  
  # ConfusionMatrixInfo(data = temp_450_wt,
  #                     other_title = paste0('null 450 wt','_',combat,'_' , method, '_', size, '_',
  #                                          num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh, '_',optimal_thresh),
  #                     predict = 'positive',
  #                     actual = 'real_label',
  #                     cutoff = optimal_thresh,
  #                     get_plot = TRUE)
  # 
  # 
  # 
  # ConfusionMatrixInfo(data = temp_850,
  #                     other_title = paste0('null 850', '_',combat,'_' , method, '_', size, '_',
  #                                          num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh, '_',optimal_thresh),
  #                     predict = 'positive',
  #                     actual = 'real_label',
  #                     cutoff = optimal_thresh,
  #                     get_plot = TRUE)
  # 
  # get dataset of predictions and labels for both small and large data
  pred_val <- prediction(temp_valid$positive, temp_valid$real)
  # pred_con <- prediction(temp_valid_850$controls_age_pred, temp_valid_850$controls_age_label)
  
  # get performace objects
  perf_s <- performance(prediction.obj = pred_val, measure = 'tpr', x.measure = 'fpr')
  # perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')
  
  # plot mean preds
  dev.off()
  plot(perf_s)
  abline(a = 0, b =1)
  temp <- round(roc(temp_valid$real, temp_valid$positive)[[9]][1], 2)
  legend(x = 0.8, y = 0.2, legend = paste0('AUC = ', temp))
  
  
  # plot age as a function of prediction
  ggplot(temp_valid,
         aes(age_sample_collection, positive)) +
    geom_point() +
    labs(title = 'Age in months vs risk score',
         x = 'Age (months)',
         y = 'Mean risk score') +
    geom_smooth(method = 'loess')
  
  
  
}



# CONTROLS
if(model_type == 'enet'){
  
  
  temp_valid <- valid_dat[, c('tm_donor','valid_age_label',  'age_diagnosis', 'age_sample_collection', 'valid_age_pred','non_zero', 'alpha', 'tech')]
  temp_con <- con_dat[, c('tm_donor','controls_age_label',  'age_sample_collection', 'controls_age_pred','non_zero', 'alpha', 'tech')]
  temp_con$pred_class <- as.factor(ifelse(temp_con$controls_age_pred > optimal_thresh, 'positive', 'negative'))
  temp_con$pred_class <- factor(temp_con$pred_class, levels = c('positive', 'negative'))
  
 
  
  temp_valid$real <- as.factor(ifelse(temp_valid$age_diagnosis <= age_cutoff, 'positive', 'negative'))
  temp_valid$real <- factor(temp_valid$real, levels = c('positive', 'negative'))
  temp_valid$pred_class <- as.factor(ifelse(temp_valid$valid_age_pred > optimal_thresh, 'positive', 'negative'))
  
  temp_valid$pred_class<- factor(temp_valid$pred_class, c('positive', 'negative'))
  temp_valid$acc <- caret::confusionMatrix(table(temp_valid$pred_class, temp_valid$real))$overall[[1]]
  acc <- round(unique(temp_valid$acc), 3)
  library(data.table)
  # get confusion matrix function for plotting 
  ConfusionMatrixInfo(data = temp_valid, 
                      predict = 'valid_age_pred', 
                      actual = 'real', 
                      cutoff = optimal_thresh,
                      other_title = paste0('valid','_',acc ,'_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_', is_lambda,'_',standardize_data,'_',alpha_num,'_',round(s_num, 4),'_', is_alpha),
                      get_plot = TRUE,
                      data_type = 'valid')
  
  # get confusion matrix function for plotting 
  ConfusionMatrixInfo(data = temp_450, 
                      predict = 'mean_pred', 
                      actual = 'controls_age_label', 
                      cutoff = optimal_thresh,
                      other_title = paste0('null 450','_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_', is_lambda,'_',standardize_data,'_',alpha_num,'_',round(s_num, 4),'_', is_alpha),
                      get_plot = TRUE,
                      data_type = 'null')
  
  # get confusion matrix function for plotting 
  ConfusionMatrixInfo(data = temp_850, 
                      predict = 'mean_pred', 
                      actual = 'controls_age_label', 
                      cutoff = optimal_thresh,
                      other_title = paste0('null 850','_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_', is_lambda,'_',standardize_data,'_',alpha_num,'_',round(s_num, 4),'_', is_alpha),
                      get_plot = TRUE,
                      data_type = 'null')
  
  # get dataset of predictions and labels for both small and large data
  pred_short <- prediction(temp_valid$valid_age_pred, temp_valid$real)
  # pred_long <- prediction(sub_dat$preds, temp_valid$real)
  
  # get performace objects
  perf_s <- performance(prediction.obj = pred_short, measure = 'tpr', x.measure = 'fpr')
  # perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')
  
  # plot mean preds
  plot(perf_s)
  abline(a = 0, b =1)
  temp <- round(roc(temp_valid$real, temp_valid$valid_age_pred)[[9]], 2)
  legend(x = 0.8, y = 0.2, legend = paste0('AUC = ', temp[[1]]))
  
  # plot age as a function of prediction
  ggplot(temp_valid, 
         aes(age_sample_collection, valid_age_pred)) + 
    geom_point() +
    labs(title = 'Age in months vs risk score',
         x = 'Age (months)',
         y = 'Mean risk score') +
    geom_smooth(method = 'loess')
  
  
} 



