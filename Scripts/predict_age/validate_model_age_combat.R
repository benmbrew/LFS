source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}

# set fixed variables
size = 'used_bh'
model_type = 'rf'
standardize = FALSE
gender = TRUE
method = 'quan'
combat = 'combat_sen'
train_lambda = TRUE
train_cutoff = FALSE
mean_lambda = FALSE
which_methyl = 'beta'
beta_thresh = 0.01
alpha_val = 0.7
age_cutoff = 72
tech = FALSE


# standardized =FALSE

# create objects to indicate method and model details when saving

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

if(mean_lambda){
  is_mean_lambda <- 'used_mean_lambda'
} else {
  is_mean_lambda <- 'no_mean_lambda'
  
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
  
  # 
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



if(model_type == 'enet'){
  if(train_lambda){
   
   lambda_val <- 'lambda_on_test'
   
  } else if (mean_lambda) {
   lambda_val <-  readRDS(paste0('final_age_combat_cv/mean_lambda',combat,'_' , method, '_', size, '_',
                                       num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
  } else {
  lambda_val <-  readRDS(paste0('final_age_combat_cv/optimal_lambda',combat,'_' , method, '_', size, '_',
                                  num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
  }
} 
    
if(model_type == 'enet'){
  
  if(train_cutoff){
    optimal_thresh <- readRDS(paste0('final_age_combat_cv/optimal_cutoff_', combat,'_' , method, '_', size, '_',
                                     num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',alpha_val,'_',beta_thresh,'.rda'))
    
  } else {
    optimal_thresh = 0.5
  }
  
} else {
  
  if(train_cutoff){
    optimal_thresh <- readRDS(paste0('final_age_combat_cv/optimal_cutoff_', combat,'_' , method, '_', size, '_',
                                     num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'_',beta_thresh,'.rda'))
    
  } else {
    optimal_thresh = 0.5
  }
  
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


# remove cancer signature
bh_feats <- bump_hunter(dat_1 = cases_450, 
                        dat_2 = con_mut , 
                        bump = 'cancer', 
                        boot_num = 5, 
                        beta_thresh = beta_thresh,
                        methyl_type = methyl_type,
                        g_ranges = g_ranges)

# combine call controls 
con_850 <- rbind(con_850, 
                 con_mut,
                 con_wt)
rm(con_mut)
rm(con_wt)

# get intersect_names
intersect_names <- names(cases_450)[grepl('^cg', names(cases_450))]

# get feature list
colnames(bh_feats)[1] <- 'chr'
remove_features <- inner_join(bh_feats, g_ranges)$probe

# take remove features out of colnames 
bh_features <- intersect_names[!intersect_names %in% remove_features]

# subset all data by bh_features
cases_450 <- remove_cancer_feats(cases_450, bh_feats = bh_features)
con_850 <- remove_cancer_feats(con_850, bh_feats = bh_features)
cases_850 <- remove_cancer_feats(cases_850, bh_feats = bh_features)



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
                                 controls = con_850,
                                 valid = cases_850,
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
  saveRDS(con_dat, paste0('final_age_combat_test/', 'con_test_', method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_', is_lambda,'_',standardize_data,'_',alpha_num,'_',s_num,'_', optimal_thresh,'.rda'))
  
  saveRDS(valid_dat, paste0('final_age_combat_test/', 'valid_test_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_', is_lambda,'_',standardize_data,'_',alpha_num,'_',s_num,'_', optimal_thresh,'.rda'))
  saveRDS(mod_dat, paste0('final_age_combat_test/', 'model_test_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_', is_lambda,'_',standardize_data,'_',alpha_num,'_',s_num,'_', optimal_thresh,'.rda'))
  
 
  
} else {
  
  
  result_list <- test_model_rf(cases = cases_450,
                               controls = con_850,
                               valid = cases_850,
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
  
  
  saveRDS(temp_con, paste0('final_age_combat_test/', 'con_test_',method,'_',size,'_',is_gen,'_',combat,'_', model_type,'_',optimal_thresh,'.rda'))
  
  saveRDS(temp_valid, paste0('final_age_combat_test/', 'valid_test_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_',optimal_thresh,'.rda'))
  
  saveRDS(temp_importance, paste0('final_age_combat_test/', 'importance_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_',optimal_thresh,'.rda'))
  saveRDS(temp_model, paste0('final_age_combat_test/', 'model_',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'_',optimal_thresh,'.rda'))
  
}

