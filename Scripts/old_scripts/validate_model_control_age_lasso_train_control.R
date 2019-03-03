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
combat = 'normal'
control_age = TRUE
which_methyl = 'beta'
beta_thresh = 0.01
optimal_cutoff = 0.5

# create objects to indicate method and model details when saving
age_cutoff = 72
trained_lambda = FALSE
tech = FALSE



if(trained_lambda){
  is_lambda <- 'lambda_test'
} else {
  is_lambda <- 'lambda_train'
}

if(control_age){
  age_control <- 'age_control'
} else {
  age_control <- 'no_age_control'
}


if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}




if(size =='used_bh'){
  cases_450 <- readRDS( paste0('../../Data/', method,'/cases_450_small_norm_', combat,'.rda'))
  cases_850<- readRDS(paste0('../../Data/', method,'/cases_850_small_norm_', combat,'.rda'))
  con_mut <- readRDS(paste0('../../Data/', method,'/con_mut_small_norm_', combat,'.rda'))
  con_850 <- readRDS(paste0('../../Data/', method,'/con_850_small_norm_', combat,'.rda'))
  con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_small_norm_', combat,'.rda'))
  
} else {
  cases_450 <- readRDS( paste0('../../Data/', method,'/cases_450_norm_', combat,'.rda'))
  cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_norm_', combat,'.rda'))
  con_mut <- readRDS(paste0('../../Data/', method,'/con_mut_norm_', combat,'.rda'))
  con_850 <- readRDS(paste0('../../Data/', method,'/con_850_norm_', combat,'.rda'))
  con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_norm_', combat,'.rda'))
  
  # add dummy tech variable for data sets with only one, replace family_name
  names(cases_450)[9] <- 'tech'
  names(con_850)[9] <- 'tech'
  names(cases_850)[9] <- 'tech'
  
  # fill them with Zero
  cases_450$tech <- '450k'
  con_850$tech <- '850k'
  cases_850$tech <- '850k'
  
  # do the same to con_mut and con_wt
  names(con_mut)[9] <- 'tech'
  names(con_wt)[9] <- 'tech'
  
  # fill new variable with right tech indication
  con_mut$tech <- '450k'
  con_wt$tech <- '450k'
  # 
  # # randomly sample from all cgs
  # clin_names <- names(cases_450)[1:11]
  # r_cgs <- sample(names(cases_450)[12:ncol(cases_450)], 3000)
  # cases_450 <- cases_450[c(clin_names, r_cgs)]
  # cases_850 <- cases_850[c(clin_names, r_cgs)]
  # con_850 <- con_850[c(clin_names, r_cgs)]
  # con_mut <- con_mut[c(clin_names, r_cgs)]
  # con_wt <- con_wt[c(clin_names, r_cgs)]
  # 
  
  
}


# # load lambda
# model_params <- readRDS(paste0('data_cv/results/model_params_', data_type, '_', base_num, '_',
#                              num_seeds, '_', k_folds, '_', is_gen, '_', model_type, '.rda'))
# 
# mean_lambda <- model_params[1]
####
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
con_all <- rbind(con_850, 
                 con_mut,
                 con_wt)
rm(con_mut, con_wt, con_850)

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

# get age factor to control for age.
get_age_fac <- function(temp_dat){
  temp_dat$age_fac <- ifelse(temp_dat$age_sample_collection >= 140, 'age_1', 'age_2')
  return(temp_dat)
}

cases_450 <- get_age_fac(cases_450)
cases_850 <- get_age_fac(cases_850)
con_all <- get_age_fac(con_all)


# age
cases_450 <- cbind(as.data.frame(class.ind(cases_450$age_fac)), 
                   cases_450)

# rempove old agevariable 
cases_450$age_fac <- NULL

# age
cases_850 <- cbind(as.data.frame(class.ind(cases_850$age_fac)), 
                   cases_850)

# rempove old agevariable 
cases_850$age_fac <- NULL

# age
con_all <- cbind(as.data.frame(class.ind(con_all$age_fac)), 
                 con_all)

# rempove old agevariable 
con_all$age_fac <- NULL

# get s_num and alpha_value
if(model_type == 'lasso'){
  # s_num = mean_lambda
  s_num = .0681
  # creat list to store results for alpha
  result_list <- list()
  con_list <- list()
  valid_list <- list()
  
  alpha_values <- 1
  
  for(i in 1:length(alpha_values)){
    alpha_num <- alpha_values[i]
    
    message('working on alpha = ', alpha_num)
    result_list[[i]] <- test_model_enet(cases = cases_450,
                                        controls = con_all,
                                        valid = cases_850,
                                        age_cutoff = age_cutoff,
                                        gender = gender,
                                        tech = tech,
                                        control_age = control_age,
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
  saveRDS(temp_con, paste0('data_test_train_control/', 'con_test_age',method,'_',size,'_',is_gen, '_', is_lambda, '_',combat,'_', model_type,'.rda'))
  
  saveRDS(temp_valid, paste0('data_test_train_control/', 'valid_test_age',method,'_',size,'_',is_gen,'_', is_lambda, '_',combat, '_', model_type,'.rda'))
  
  
} else {
  
  
  result_list <- test_model_rf(cases = cases_450,
                               controls = con_all,
                               valid = cases_850,
                               age_cutoff = age_cutoff,
                               gender = gender,
                               tech = tech,
                               control_age = control_age,
                               bh_features = bh_features)
  
  temp_valid <- result_list[[1]]
  temp_con <- result_list[[2]]
  temp_importance  <- result_list[[3]]
  
  saveRDS(temp_con, paste0('data_test_train_control/', 'con_test_control_age', method,'_',size,'_',is_gen, '_',combat,'_', model_type,'.rda'))
  
  saveRDS(temp_valid, paste0('data_test_train_control/', 'valid_test_control_age',method,'_',size,'_',is_gen, '_',combat,'_', model_type,'.rda'))
  
  saveRDS(temp_importance, paste0('data_test_train_control/', 'importance_control_age', method,'_',size,'_',is_gen,'_',combat,'_', model_type,'.rda'))
  
}

