#source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  
data_type = 'm'
combat = FALSE

if(combat){
  used_combat <- 'used_combat'
} else {
  used_combat <- 'no_combat'
}
if(data_type == 'beta'){
  beta_thresh = 0.0001
} else {
  beta_thresh = 0.001
}


# read in all data
cases_450 <- readRDS(paste0('residual_data_cv/', 'cases_450_',data_type, '_',used_combat,'.rda'))
con_850 <- readRDS(paste0('residual_data_cv/', 'con_850_',data_type,'_',used_combat, '.rda'))
cases_850 <- readRDS(paste0('residual_data_cv/', 'cases_850_',data_type,'_',used_combat,'.rda'))
con_wt <-  readRDS(paste0('residual_data_cv/', 'con_wt_',data_type,'_',used_combat, '.rda'))
lfs_bump_probes <- readRDS(paste0('transform_data_cv/', 'lfs_bumps_', data_type,'_.rda'))

# subset by lfs_bump_probes
cases_clin <- names(cases_450)[!grepl('^cg', names(cases_450))]
cases_450 <- cases_450[c(cases_clin, lfs_bump_probes)]

cases_clin_450 <- names(cases_850)[!grepl('^cg', names(cases_850))]
cases_850 <- cases_850[c(cases_clin, lfs_bump_probes)]

con_clin_850 <- names(con_850)[!grepl('^cg', names(con_850))]
con_850 <- cases_450[c(cases_clin, lfs_bump_probes)]

con_wt <- con_wt[c(con_clin_850, lfs_bump_probes)]


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
model_type = 'rf'
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

optimal_cutoff = 0.5
bh_features <- lfs_bump_probes
# get s_num and alpha_value
if(model_type == 'enet'){
  s_num = 1
  # s_num <- round(model_params[[1]], 3)
  # alpha_num <- round(model_params[[2]], 2)
  # 
  # creat list to store results for alpha
  result_list <- list()
  con_list <- list()
  valid_list <- list()
  
  alpha_values <- (1:10/10)
  
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
                                        bh_features = lfs_bump_probes)
    
    con_list[[i]] <- result_list[[i]][[1]]
    valid_list[[i]] <- result_list[[i]][[2]]
    
  }
  
  
  
  temp_con <- do.call('rbind', con_list)
  temp_valid <- do.call('rbind', valid_list)
  
  # read in cases_450
  saveRDS(temp_con, paste0('residual_data_test/', 'con_test_resid',data_type, '_', is_gen, '_',is_lambda, '_',used_combat,'_', model_type,'.rda'))
  
  saveRDS(temp_valid, paste0('residual_data_test/', 'valid_test_resid',data_type, '_', is_gen,'_',is_lambda, '_',used_combat, '_', model_type,'.rda'))
  
  
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
  
  saveRDS(temp_con, paste0('residual_data_test/', 'con_test_resid',data_type, '_', is_gen, '_',used_combat,'_', model_type,'.rda'))
  
  saveRDS(temp_valid, paste0('residual_data_test/', 'valid_test_resid',data_type, '_', is_gen,'_',used_combat, '_', model_type,'.rda'))
  
  saveRDS(temp_importance, paste0('residual_data_test/', 'importance_resid',data_type, '_', is_gen,'_',used_combat, '_', model_type,'.rda'))
  
}





