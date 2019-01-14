# prepare pc data
source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


# set fixed variables
data_type = 'beta'
combat = TRUE
remove_leading_pcs = 'first'


# set saving objects
if(combat){
  used_combat <- 'used_combat'
} else {
  used_combat <- 'no_combat'
}

# condition on fixed objects to get saving identifiers
if(data_type == 'beta'){
  which_methyl <- 'beta'
  beta_thresh = 0.03
  
} else {
  which_methyl <- 'm'
  beta_thresh= 0.3
}

# read in data
cases_450 <- readRDS(paste0('pc_data_cv/cases_450', remove_leading_pcs,'_', data_type,'_', 
                            used_combat, '.rda'))
cases_850 <- readRDS(paste0('pc_data_cv/cases_850',  remove_leading_pcs, '_',data_type,  '_', 
                            used_combat, '.rda'))
con_mut <- readRDS(paste0('pc_data_cv/con_mut',remove_leading_pcs,'_',data_type,  '_', 
                          used_combat, '.rda'))
con_850 <- readRDS(paste0('pc_data_cv/con_850', remove_leading_pcs, '_', data_type, '_', 
                          used_combat, '.rda'))
con_wt <- readRDS(paste0('pc_data_cv/con_wt',  remove_leading_pcs, '_',data_type,  '_', 
                         used_combat, '.rda'))
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

# create range for random sampling for cross validation
seed_range <- c(1:how_many_seeds)

# create list to store model
all_test_results <- list()
importance_results <- list()
rf_pred_results <- list()
rf_important_results <- list()
# save image
# save.image('~/Desktop/temp_valid.RData')

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
                        g_ranges = g_ranges) {
    
    
    
    set.seed(random_seed)
    
    # get vector of random folds
    fold_vec <- sample(1:k_folds, nrow(cases_full), replace = T)
    
    # test_data_results
    test_data_results <- list()
    
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
        mod_result  <- run_enet_all_test(training_dat = train_cases,
                                         test_dat = test_cases,
                                         controls_dat = con_transform,
                                         valid_dat = valid_transform,
                                         age_cutoff = age_cutoff,
                                         gender = gender,
                                         tech = tech,
                                         bh_features = bh_features)
      } else {
        mod_result  <- run_rf_all_test(training_dat = train_cases,
                                       test_dat = test_cases,
                                       controls_dat = con_transform,
                                       valid_dat = valid_transform,
                                       age_cutoff = age_cutoff,
                                       gender = gender,
                                       tech = tech,
                                       bh_features = bh_features)
        
        importance_rf <- mod_result[[2]]
        mod_result <- mod_result[[1]]
      }
      
      
      
      if(model_type == 'rf'){
        importance_results[[i]] <- importance_rf
      }
      
      mod_result$seed <- random_seed
      mod_result$fold <- i
      test_data_results[[i]] <- mod_result
      print(i)
    }
    
    
    # combine list of case and control result data frames and return all result objects (two dataframes and 4 lists)
    cv_testing_results <- do.call(rbind, test_data_results)
    
    if(model_type == 'rf'){
      cv_importance_results <- do.call(rbind, importance_results)
      return(list(cv_testing_results, cv_importance_results))
    } else {
      return(cv_testing_results)
      
    }
  }
  
  if(model_type == 'rf'){
    
    
    # run model with 5 k fold cross validation
    temp_test_results <- run_model(cases_full = cases_450,
                                   controls_full = con_850,
                                   valid_full = valid_850,
                                   model_type = model_type,
                                   age_cutoff = age_cutoff,
                                   k_folds = k_folds,
                                   tech = tech,
                                   gender = gender,
                                   beta_thresh = beta_thresh,
                                   g_ranges = g_ranges)
    rf_pred_results[[random_seed]] <- temp_test_results[[1]]
    rf_important_results[[random_seed]] <- temp_test_results[[2]]
    
    
    message('finished working on random seed = ', random_seed)
  } else {
    
    # run model with 5 k fold cross validation
    all_test_results[[random_seed]] <- run_model(cases_full = cases_450,
                                                 controls_full = con_850,
                                                 valid_full = valid_850,
                                                 model_type = model_type,
                                                 age_cutoff = age_cutoff,
                                                 k_folds = k_folds,
                                                 tech = tech,
                                                 gender = gender,
                                                 beta_thresh = beta_thresh,
                                                 g_ranges = g_ranges)
    
    message('finished working on random seed = ', random_seed)
  }
  
  
}


if(model_type == 'rf'){
  final_dat <- do.call(rbind, rf_pred_results)
  final_importance <- do.call(rbind, rf_important_results )
  # # save data
  saveRDS(final_dat, paste0('pc_data_cv/results/', data_type, '_', remove_leading_pcs,'_' ,
                            num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))
  # # save data
  saveRDS(final_importance, paste0('pc_data_cv/results/', 'importance_',data_type, '_',remove_leading_pcs,'_' ,
                                   num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))
  
  
} else {
  final_dat <- do.call(rbind, all_test_results)
  
  # # save data
  saveRDS(final_dat, paste0('pc_data_cv/results', data_type, '_',remove_leading_pcs,'_' ,
                            num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))
  
}



