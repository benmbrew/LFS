source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}

# strategy 
# all methods, both models, both sizes
# combat only normal, combat_1

size = 'used_bh'
model_type = 'rf'
null_450= TRUE
null_450_all = FALSE
use_p53 = FALSE
gender = TRUE
use_cancer = TRUE
method = 'noob'
include_under_6 = FALSE
combat = 'combat_1'
alpha_val = 0.1
train_alpha = FALSE
train_lambda = FALSE
train_cutoff = TRUE
mean_lambda = FALSE
which_methyl = 'beta'
beta_thresh = 0.05
age_cutoff = 72
tech = FALSE
standardize = FALSE


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
# create list to store model
all_test_results <- list()
importance_results <- list()
rf_pred_results <- list()
rf_pred_results_con <- list()
rf_pred_results_valid <- list()
rf_important_results <- list()

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

if(use_null_450 == 'used_null_450_all'){
  null_450_all <- rbind(con_mut,
                        con_wt)
  null_450_all$age <- null_450_all$age_sample_collection
  if(!include_under_6){
    null_450_all <- null_450_all[null_450_all$age > 72, ]
    
  }
  cases_450$age <- cases_450$age_diagnosis
  cases_450 <- rbind(cases_450,
                 null_450_all)

  cases_450 <- cases_450[!duplicated(cases_450$tm_donor),]
  
  
  
} else  if(use_null_450 == 'used_null_450_mut'){
  null_450_all <- con_mut
  null_450_all$age <- null_450_all$age_sample_collection
  if(!include_under_6){
    null_450_all <- null_450_all[null_450_all$age > 72, ]
    
  }
  cases_450$age <- cases_450$age_diagnosis
  cases_450 <- rbind(cases_450,
                     null_450_all)
  
  cases_450 <- cases_450[!duplicated(cases_450$tm_donor),]
  
  
}
  
  



# create objects to indicate method and model details when saving

# save image
# save.image('~/Desktop/temp_valid.RData')
cases_full <- cases_450
controls_full <- con_850
valid_full <- cases_850
for(random_seed in 1:length(seed_range)) {
  # prepare data sets for modelling
  
  run_model <- function(cases_full,
                        controls_full,
                        valid_full,
                        model_type = model_type,
                        age_cutoff = age_cutoff,
                        k_folds = k_folds,
                        tech = tech,
                        gender = gender,
                        beta_thresh = beta_thresh,
                        methyl_type = methyl_type,
                        g_ranges = g_ranges,
                        standardize = standardize) {
    
    
    
    set.seed(random_seed)
    
    # get vector of random folds
    fold_vec <- sample(1:k_folds, nrow(cases_full), replace = T)
    
    # test_data_results
    
    importance_rf <- list()
    temp_results <- list()
    temp_results_con <- list()
    
    
    # combine 
    for(i in 1:k_folds) {
      
      # get train and test index
      train_index <- !grepl(i, fold_vec)
      test_index <- !train_index
      
      # subset to get training and test data
      train_cases <- cases_full[train_index, ]
      test_cases <- cases_full[test_index, ]
      
      # get controls 450k mut
      temp_con <- controls_full[controls_full$tech == '850k' & controls_full$p53_germline == 'MUT',]
      temp_con <- temp_con[!duplicated(temp_con$tm_donor),]
      temp_con$age <- temp_con$age_sample_collection
      
      # use cases training and controls to get bumphunter features
      bh_feats <- bump_hunter(dat_1 = train_cases, 
                              dat_2 = temp_con , 
                              bump = 'cancer', 
                              boot_num = 5, 
                              beta_thresh = beta_thresh,
                              methyl_type = methyl_type,
                              g_ranges = g_ranges)
      rm(temp_con)
      
      
      # get intersect_names
      intersect_names <- names(train_cases)[grepl('^cg', names(train_cases))]
      
      # get feature list
      colnames(bh_feats)[1] <- 'chr'
      remove_features <- inner_join(bh_feats, g_ranges)$probe
      
      # take remove features out of colnames 
      bh_features <- intersect_names[!intersect_names %in% remove_features]
      
      if(model_type == 'enet'){
        # function to predict with all test, controls, controls old, and valid
        
        run_enet_test(training_dat = train_cases,
                      test_dat = test_cases,
                      controls_dat = controls_full,
                      valid_dat = cases_850,
                      age_cutoff = age_cutoff,
                      alpha_num = alpha_val,
                      gender = gender,
                      tech = tech,
                      bh_features = bh_features,
                      standardize = standardize)
        
        temp_results[[i]] <- mod_results[[i]][[1]]
        temp_results_con[[i]] <- mod_results[[i]][[2]]
      } else {
        # mod_result  <- run_rf_all_test(training_dat = train_cases,
        #                                test_dat = test_cases,
        #                                controls_dat = controls_full,
        #                                valid_dat = cases_850,
        #                                age_cutoff = age_cutoff,
        #                                gender = gender,
        #                                tech = tech,
        #                                bh_features = bh_features)
        
        mod_results  <- run_rf_test(training_dat = train_cases,
                                         test_dat = test_cases,
                                         controls_dat = controls_full,
                                         valid_dat = cases_850,
                                         age_cutoff = age_cutoff,
                                         gender = gender,
                                         tech = tech,
                                         bh_features = bh_features)
        
        importance_rf[[i]] <- mod_results[[2]]
        temp_results[[i]] <- mod_results[[1]]
        temp_results_con[[i]] <- mod_results[[3]]
        
      
      }
      
      
      
      print(i)
    }
    
    test_result <- do.call('rbind',temp_results)
    con_result <- do.call('rbind',temp_results_con)
    importance_result <- do.call('rbind', importance_rf)
    test_result$seed <- random_seed
    con_result$seed <- random_seed

   
    
    
    # combine list of case and control result data frames and return all result objects (two dataframes and 4 lists)
    
   return(list(test_result, con_result, importance_result))
  }
  
  if(model_type == 'rf'){
    
    
    # run model with 5 k fold cross validation
    temp_test_results <- run_model(cases_full = cases_450,
                                   controls_full = con_850,
                                   valid_full = cases_850,
                                   model_type = model_type,
                                   age_cutoff = age_cutoff,
                                   k_folds = k_folds,
                                   tech = tech,
                                   gender = gender,
                                   beta_thresh = beta_thresh,
                                   g_ranges = g_ranges,
                                   standardize = standardize)
    rf_pred_results[[random_seed]] <- temp_test_results[[1]]
    rf_pred_results_con[[random_seed]] <- temp_test_results[[2]]
    rf_important_results[[random_seed]] <- temp_test_results[[3]]
    
    
    message('finished working on random seed = ', random_seed)
  } else {
    
    # run model with 5 k fold cross validation
    all_test_results[[random_seed]] <- run_model(cases_full = cases_450,
                                                 controls_full = con_850,
                                                 valid_full = cases_850,
                                                 model_type = model_type,
                                                 age_cutoff = age_cutoff,
                                                 k_folds = k_folds,
                                                 tech = tech,
                                                 gender = gender,
                                                 beta_thresh = beta_thresh,
                                                 g_ranges = g_ranges,
                                                 standardize = standardize)
    
    message('finished working on random seed = ', random_seed)
  }
  
  
}


if(model_type == 'rf'){
  final_dat <- do.call(rbind, rf_pred_results)
  final_dat_con <- do.call(rbind, rf_pred_results_con)
  # final_dat_val <- do.call(rbind, rf_pred_results_valid)
  
  final_importance <- do.call(rbind, rf_important_results )
  # # save data
  saveRDS(final_dat, paste0('pc_data_cv/new_results/', combat,'_' , method, '_',size, '_',
                            num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))
  saveRDS(final_dat_con, paste0('pc_data_cv/new_results/', 'con_', combat,'_' , method, '_',size, '_',
                            num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))
  
  # saveRDS(final_dat_val, paste0('pc_data_cv/new_results/', 'valid_', combat,'_' , method, '_',size, '_',
                                # num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))
  # 
  # # save data
  saveRDS(final_importance, paste0('pc_data_cv/new_results/', 'importance_', combat,'_' ,method, '_',size, '_',
                                   num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))
  
  
} else {
  final_dat <- do.call(rbind, all_test_results)
  
  # # save data
  saveRDS(final_dat, paste0('pc_data_cv/new_results/',combat,'_' ,method, '_', size, '_',
                            num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_', standardize_data,'.rda'))
  
}


temp <- final_dat %>% 
  group_by(tm_donor) %>% 
  summarise(positive = mean(positive),
            age = mean(age))

temp$real <- as.factor(ifelse(temp$age < 72, 'positive', 'negative'))
temp$real <- factor(temp$real, levels = c('positive', 'negative'))
temp$pred_class <- as.factor(ifelse(temp$positive > 0.5, 'positive', 'negative'))
temp$pred_class <- factor(temp$pred_class, levels = c('positive', 'negative'))
temp$acc <- caret::confusionMatrix(table(temp$pred_class, temp$real))$overall[[1]]


ConfusionMatrixInfo(data = temp, 
                    predict = 'positive', 
                    actual = 'real', 
                    cutoff = 0.5, 
                    get_plot = TRUE, 
                    other_title = '', 
                    data_type = 'valid')

# get dataset of predictions and labels for both small and large data
pred_val <- prediction(temp$positive, temp$real)
# pred_con <- prediction(temp_valid_850$controls_age_pred, temp_valid_850$controls_age_label)

# get performace objects
perf_s <- performance(prediction.obj = pred_val, measure = 'tpr', x.measure = 'fpr')
# perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')

# plot mean preds
dev.off()
plot(perf_s)
abline(a = 0, b =1)
temp_roc <- round(roc(temp$real, temp$positive)[[9]][1], 2)
legend(x = 0.8, y = 0.2, legend = paste0('AUC = ', temp_roc))


# plot age as a function of prediction
ggplot(temp,
       aes(age, positive)) +
  geom_point() +
  labs(title = 'Age in months vs risk score',
       x = 'Age (months)',
       y = 'Mean risk score') +
  geom_smooth(method = 'loess')
