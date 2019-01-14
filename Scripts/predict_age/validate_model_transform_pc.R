# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  
data_type = 'beta'

if(data_type == 'beta'){
  beta_thresh = 0.05
} else {
  beta_thresh = 0.5
}
remove_leading_pcs = 'first'

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

# load model params
# load lambda
model_params <- readRDS(paste0('transform_pc_data_cv/results/model_params_', data_type, '_', remove_leading_pcs,'_',
                               num_seeds, '_', k_folds, '_', is_gen, '_', model_type,'.rda'))

mean_lambda <- model_params[1]

# read in all data
cases_450 <- readRDS(paste0('transform_pc_data/', 'cases_450_pc_', remove_leading_pcs,'_',data_type,'.rda'))
con_850 <- readRDS( paste0('transform_pc_data/', 'con_850_pc_',remove_leading_pcs,'_',data_type, '.rda'))
cases_850 <- readRDS(paste0('transform_pc_data/', 'cases_850_pc_',remove_leading_pcs,'_',data_type,'.rda'))
con_wt <-  readRDS(paste0('transform_pc_data/', 'con_wt_pc_',remove_leading_pcs,'_',data_type, '.rda'))

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
rm(con_mut, con_wt)

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

# create objects to indicate method and model details when saving
gender = FALSE
tech = FALSE
how_many_seeds = 20
how_many_folds = 5
age_cutoff = 72
model_type = 'enet'
trained_lambda = TRUE

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

# get s_num and alpha_value
if(model_type == 'enet'){
  s_num = mean_lambda
  s_num = .128
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
                                        bh_features = bh_features)
    
    con_list[[i]] <- result_list[[i]][[1]]
    valid_list[[i]] <- result_list[[i]][[2]]
    
  }
  
  
  
  temp_con <- do.call('rbind', con_list)
  temp_valid <- do.call('rbind', valid_list)
  
  # read in cases_450
  saveRDS(temp_con, paste0('transform_pc_data_test/', 'con_test',data_type, '_', is_gen, '_', remove_leading_pcs, '_',is_lambda,'_', model_type,'.rda'))
  
  saveRDS(temp_valid, paste0('transform_pc_data_test/', 'valid_test',data_type, '_', is_gen,'_', remove_leading_pcs,'_',is_lambda, '_', model_type,'.rda'))
  
  
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
  
  saveRDS(temp_con, paste0('data_test/', 'con_test',data_type, '_', is_gen, '_',base_num,'_',used_combat,'_', model_type,'.rda'))
  
  saveRDS(temp_valid, paste0('data_test/', 'valid_test',data_type, '_', is_gen, '_',base_num,'_',used_combat,'_', model_type,'.rda'))
  
  saveRDS(temp_importance, paste0('data_test/', 'importance',data_type, '_', is_gen, '_',base_num,'_',used_combat,'_', model_type,'.rda'))
  
}

