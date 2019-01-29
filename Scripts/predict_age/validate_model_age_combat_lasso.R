source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


# set fixed variables
size = 'full'
model_type = 'rf'
gender = FALSE
method = 'noob'
combat = 'combat_1'
which_methyl = 'beta'
beta_thresh = 0.01
optimal_cutoff = 0.5

# create objects to indicate method and model details when saving
age_cutoff = 72
trained_lambda = FALSE
tech = FALSE

if(control_age){
  age_control <- 'age_control'
} else {
  age_control <- 'no_age_control'
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
  
}


# # load lambda
# model_params <- readRDS(paste0('pc_data_cv/results/model_params_', data_type, '_', remove_leading_pcs, '_',
#                                num_seeds, '_', k_folds, '_', is_gen, '_', model_type, '.rda'))
# 
# mean_lambda <- model_params[1]
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
if(model_type == 'lasso'){
  # s_num = mean_lambda
  s_num = 0.0154
  # s_num <- round(model_params[[1]], 3)
  # alpha_num <- round(model_params[[2]], 2)
  # 
  # creat list to store results for alpha
  result_list <- list()
  con_list <- list()
  valid_list <- list()
  
  alpha_values <- 1
  
  for(i in 1:length(alpha_values)){
    alpha_num <- alpha_values[i]
    
    message('working on alpha = ', alpha_num)
    result_list[[i]] <- test_model_enet(cases = cases_450,
                                        controls = con_850,
                                        valid = cases_850,
                                        age_cutoff = age_cutoff,
                                        gender = gender,
                                        tech = tech,
                                        cv_lambda = trained_lambda,
                                        alpha_value = alpha_num,
                                        lambda_value = s_num,
                                        bh_features = bh_features)
    
    con_list[[i]] <- result_list[[i]][[1]]
    valid_list[[i]] <- result_list[[i]][[2]]
    
  }
  
  
  
  temp_con <- do.call('rbind', con_list)
  temp_valid <- do.call('rbind', valid_list)
  
  # read in cases_450
  saveRDS(temp_con, paste0('age_combat_data_test/', 'con_test', method,'_',size,'_',is_gen, '_',combat,'_', model_type,'.rda'))
  
  saveRDS(temp_valid, paste0('age_combat_data_test/', 'valid_test',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'.rda'))
  
  
} else {
  
  
  result_list <- test_model_rf(cases = cases_450,
                               controls = con_850,
                               valid = cases_850,
                               age_cutoff = age_cutoff,
                               gender = gender,
                               tech = tech,
                               bh_features = bh_features)
  
  temp_valid <- result_list[[1]]
  temp_con <- result_list[[2]]
  temp_importance  <- result_list[[3]]
  
  
  saveRDS(temp_con, paste0('age_combat_data_test/', 'con_test_pc',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'.rda'))
  
  saveRDS(temp_valid, paste0('age_combat_data_test/', 'valid_test_pc',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'.rda'))
  
  saveRDS(temp_importance, paste0('age_combat_data_test/', 'importance_pc',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'.rda'))
  
}

