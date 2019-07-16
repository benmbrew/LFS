source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}

# strategy 
# all methods, both models, both sizes
# combat only normal, combat_1

# next do swan with standardize and not
# set fixed variables
size = 'used_bh'
model_type = 'rf'
control_age = FALSE
gender = FALSE
method = 'swan'
combat = 'normal'
which_methyl = 'beta'
beta_thresh = 0.05
alpha_val= 0.1
tech = FALSE
standardize = FALSE

# create objects to indicate method and model details when saving
age_cutoff = 72
trained_lambda = FALSE


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
how_many_seeds = 10
how_many_folds = 5


num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds

# create range for random sampling for cross validation
seed_range <- c(1:how_many_seeds)

# create list to store model
all_test_results <- list()
importance_results <- list()
pred_results <- list()
pred_results_con <- list()

important_results <- list()


if(size =='used_bh'){
  cases_450 <-  readRDS(paste0('../../Data/', method,'/cases_450_small_cv_age_combat_', combat,'.rda'))
  cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_small_cv_age_combat_', combat,'.rda'))
  con_mut <- readRDS(paste0('../../Data/', method,'/con_mut_small_cv_age_combat_', combat,'.rda'))
  con_850 <- readRDS( paste0('../../Data/', method,'/con_850_small_cv_age_combat_', combat,'.rda'))
  con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_small_cv_age_combat_', combat,'.rda'))
  
} else {
  cases_450 <-  readRDS(paste0('../../Data/', method,'/cases_450_cv_age_combat_', combat,'.rda'))
  cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_cv_age_combat_', combat,'.rda'))
  con_mut <- readRDS(paste0('../../Data/', method,'/con_mut_cv_age_combat_', combat,'.rda'))
  con_850 <- readRDS( paste0('../../Data/', method,'/con_850_cv_age_combat_', combat,'.rda'))
  con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_cv_age_combat_', combat,'.rda'))
  
  # if(model_type == 'enet'){
  #   num_probes = 10000
  # } else {
  #   num_probes = 5000
  # }
  # 
  # # randomly sample from all cgs
  # clin_names <- names(cases_450)[1:12]
  # r_cgs <- sample(names(cases_450)[13:ncol(cases_450)], num_probes)
  # cases_450 <- cases_450[c(clin_names, r_cgs)]
  # con_wt <- con_wt[c(clin_names, r_cgs)]
  # con_mut <- con_mut[c(clin_names, r_cgs)]
  # con_850 <- con_850[c(clin_names, r_cgs)]
  # cases_850 <- cases_850[c(clin_names, r_cgs)]
  # 
}

# combine controls
controls_full <- rbind(con_850,
                       con_mut,
                       con_wt)

rm(con_850, con_mut, con_wt)
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

# create objects to indicate method and model details when saving

# save image
# save.image('~/Desktop/temp_valid.RData')
cases_full = cases_450
controls_full = controls_full
valid_full = cases_850
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
                        alpha_val = alpha_val) {
    
    
    
    set.seed(random_seed)
    
    # get vector of random folds
    fold_vec <- sample(1:k_folds, nrow(cases_full), replace = T)
    
    # get nrow for tot probes in rf
    tot_probes <- ncol(cases_full)
    # test_data_results
    test_data_results <- list()
    mod_results <- list()
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
        mod_results[[i]]  <- run_enet_test(training_dat = train_cases,
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
        
        
        mod_results[[i]]  <- run_rf_test(training_dat = train_cases,
                                       test_dat = test_cases,
                                       controls_dat = controls_full,
                                       valid_dat = cases_850,
                                       age_cutoff = age_cutoff,
                                       gender = gender,
                                       tech = tech,
                                       bh_features = bh_features)
        
        importance_rf[[i]] <- mod_results[[i]][[2]]
        temp_results[[i]] <- mod_results[[i]][[1]]
        temp_results_con[[i]] <- mod_results[[i]][[3]]
        
      }
      
      
     
      
    }
    library(DescTools)
    
    
    if(model_type == 'rf'){
      test_results <- do.call(rbind, temp_results)
      test_results_con <- do.call(rbind, temp_results_con)
      
      final_imp <- do.call(rbind, importance_rf)
      final_imp$seed <- random_seed
      test_results$seed <- random_seed
      test_results_con$seed <- random_seed
      
      roc_curve = roc(real ~ positive , data = test_results)
      
      # coords(roc=roc_curve, x = "local maximas", ret='threshold')
      thresh <- as.data.frame(t(coords(roc=roc_curve, x = "all")))
      thresh <- thresh[order(thresh$threshold, decreasing = FALSE),]
      optimal_thresh <- thresh$threshold[which(thresh$sensitivity == Closest(thresh$sensitivity, thresh$specificity))][1]
      
      test_results$pred_class <- as.factor(ifelse(test_results$positive > 0.5, 'positive', 'negative'))
      
      # relevel both factore
      test_results$pred_class <- factor(test_results$pred_class, c('positive', 'negative'))
      test_results$real <- factor(test_results$real, c('positive', 'negative'))
      
      test_results$accuracy <- caret::confusionMatrix(table(test_results$pred_class, test_results$real))$overall[1]
      
      
      
      test_results$pred_class_opt <- as.factor(ifelse(test_results$positive > optimal_thresh, 'positive', 'negative'))
      
      # relevel both factore
      test_results$pred_class_opt <- factor(test_results$pred_class_opt, c('positive', 'negative'))
      test_results$real <- factor(test_results$real, c('positive', 'negative'))
      
      test_results$accuracy_opt <- caret::confusionMatrix(table(test_results$pred_class_opt, test_results$real))$overall[1]
      test_results$tot_probes <- tot_probes
    } else {
      # collaspe model_resutl
      test_results <- do.call(rbind, temp_results)
      test_results_con <- do.call(rbind, temp_results_con)
      
      test_results$seed <- random_seed
      test_results_con$seed <- random_seed
      
      
      
      roc_curve = roc(real ~ preds , data = test_results)
      
      # coords(roc=roc_curve, x = "local maximas", ret='threshold')
      thresh <- as.data.frame(t(coords(roc=roc_curve, x = "all")))
      thresh <- thresh[order(thresh$threshold, decreasing = FALSE),]
      optimal_thresh <- thresh$threshold[which(thresh$sensitivity == Closest(thresh$sensitivity, thresh$specificity))][1]
      
      test_results$pred_class <- as.factor(ifelse(test_results$preds > .5, 'positive', 'negative'))
      
      # relevel both factore
      test_results$pred_class <- factor(test_results$pred_class, c('positive', 'negative'))
      test_results$real <- factor(test_results$real, c('positive', 'negative'))
      test_results$accuracy <- caret::confusionMatrix(table(test_results$pred_class, test_results$real))$overall[1]
      
      test_results$pred_class_opt <- as.factor(ifelse(test_results$preds > optimal_thresh, 'positive', 'negative'))
      
      # relevel both factore
      test_results$pred_class_opt <- factor(test_results$pred_class_opt, c('positive', 'negative'))
      test_results$accuracy_opt <- caret::confusionMatrix(table(test_results$pred_class_opt, test_results$real))$overall[1]
      
      
      test_results_con$pred_class <- as.factor(ifelse(test_results_con$preds > optimal_thresh, 'positive', 'negative'))
      
      # relevel both factore
      test_results_con$pred_class <- factor(test_results_con$pred_class, c('positive', 'negative'))

      test_results_con$optimal_thresh = optimal_thresh
      
    }
    
    if(model_type == 'rf'){
      return(list(test_results, final_imp, test_results_con))
    } else{
      return(list(test_results, test_results_con))
    }
    
  }
  
  if(model_type == 'rf'){
    
    
    # run model with 5 k fold cross validation
    temp_test_results <- run_model(cases_full = cases_450,
                                   controls_full = controls_full,
                                   valid_full = cases_850,
                                   model_type = model_type,
                                   age_cutoff = age_cutoff,
                                   k_folds = k_folds,
                                   tech = tech,
                                   gender = gender,
                                   beta_thresh = beta_thresh,
                                   g_ranges = g_ranges,
                                   alpha_val = alpha_val)
    pred_results[[random_seed]] <- temp_test_results[[1]]
    important_results[[random_seed]] <- temp_test_results[[2]]
    pred_results_con[[random_seed]] <- temp_test_results[[3]]
    
    
    
    message('finished working on random seed = ', random_seed)
  } else {
    
    # run model with 5 k fold cross validation
    temp_test_results  <- run_model(cases_full = cases_450,
                                    controls_full = controls_full,
                                    valid_full = cases_850,
                                    model_type = model_type,
                                    age_cutoff = age_cutoff,
                                    k_folds = k_folds,
                                    tech = tech,
                                    gender = gender,
                                    beta_thresh = beta_thresh,
                                    g_ranges = g_ranges,
                                    alpha_val = alpha_val)
    
    pred_results[[random_seed]] <- temp_test_results[[1]]
    pred_results_con[[random_seed]] <- temp_test_results[[2]]
    
    message('finished working on random seed = ', random_seed)
  }
  
  
}


if(model_type == 'rf'){
  final_dat <- do.call(rbind, pred_results)
  final_dat_con <- do.call(rbind, pred_results_con)
  
  final_importance <- do.call(rbind, important_results)
  # # save data
  saveRDS(final_dat, paste0('final_age_combat_cv/new_results/', combat,'_' , method, '_',size, '_',
                            num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh,'.rda'))
  saveRDS(final_dat_con, paste0('final_age_combat_cv/new_results/con_450_', combat,'_' , method, '_',size, '_',
                            num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh,'.rda'))
  # # save data
  saveRDS(final_importance, paste0('final_age_combat_cv/new_results/', 'importance_', combat,'_' ,method, '_',size, '_',
                                   num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh,'.rda'))
  
  
} else {
  final_dat <- do.call(rbind, pred_results)
  final_dat_con <- do.call(rbind, pred_results_con)
  
  
  # # save data
  saveRDS(final_dat, paste0('final_age_combat_cv/new_results/',combat,'_' ,method, '_', size, '_',
                            num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_', alpha_val,'_',beta_thresh,'.rda'))
  saveRDS(final_dat_con, paste0('final_age_combat_cv/new_results/','con_450_', combat,'_' ,method, '_', size, '_',
                            num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_', alpha_val,'_',beta_thresh,'.rda'))
  
}



