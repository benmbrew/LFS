
# old_m_processed_outlier_fun, 
# new_m_processed_outlier, new_m_processed_outliner_more, 
# new_m_processed_outlier_more_funnorm
# validation or combined data 
image_ids <- read_csv('~/Documents/LFS_shiny/data/image_ids.csv')
data_used <- 'new'

# get cg regions
cg_gene_regions = 'Body'

# set preprocessing method
method <- 'noob'

# set type of data, beta or m
methyl_type <- 'beta'

# set data directory
data_dir <- '../../Data/'

# 
# # new_m_final_outlier
# # new_m_final_outlier_more
# # new_m_final_outlier_more_bad
# # new_m_final_outlier_funnorm_more
# temp_cancer <- full_data[!grepl('Unaffected', full_data$cancer_diagnosis_diagnoses),]
# temp <- as.data.frame(cbind(onset = temp_cancer$age_diagnosis, sample = temp_cancer$age_sample_collection))
# 
# 
# 
# 
# 
# temp_controls <- full_data[grepl('Unaffected', full_data$cancer_diagnosis_diagnoses),]
# length(which(duplicated(temp_contorls$ids)))
# 
# temp <- cbind(temp_cancer$age_diagnosis, temp_cancer$age_sample_collection)
# length(which(temp_cancer$age_diagnosis >72))
# get data
load(paste0(data_dir,paste0(data_used,'_',methyl_type, '_final_beta_first_last', '.RData')))

full_data <- full_data_first
full_data <- full_data[!duplicated(full_data$ids),]


length(which(duplicated(full_data$tm_donor_)))






# source all_functions.R to load libraries and my functions
source('all_functions.R')


# read in wt data
data_wt <- readRDS(paste0(data_dir,paste0('new','_',methyl_type, '_wild_type', '.rda')))

# full_data <- full_data[full_data$tm_donor_ != '3955',]
# full_data_fun_fromfirst <- full_data_fun_fromfirst[full_data_fun_fromfirst$tm_donor_ != '3955',]
# full_data_fun_fromlast <- full_data_fun_fromlast[full_data_fun_fromlast$tm_donor_ != '3955',]
full_data_first <- full_data_first[full_data_first$tm_donor_ != '3955',]
full_data_first_combat <- full_data_first_combat[full_data_first_combat$tm_donor_ != '3955',]

full_data_last <- full_data_last[full_data_last$tm_donor_ != '3955',]
full_data_last_combat <- full_data_last_combat[full_data_last_combat$tm_donor_ != '3955',]

# temp <- get_diff_dups(full_data)
# full_data_first <- temp[[1]]
# full_data_last <- temp[[2]]

# temp <- get_diff_dups(full_data_combat)
# full_data_first_combat <- temp[[1]]
# full_data_last_combat <- temp[[2]]

temp <- get_diff_dups(full_data_first)
full_data_first <- temp[[1]]
full_data_last <- temp[[2]]
temp_combat <- get_diff_dups(full_data_first_combat)
full_data_first_combat <- temp_combat[[1]]
full_data_last_combat <- temp_combat[[2]]
# get controls - HERE
# get controls that are deduped from both sides

# # get age numeric category
# full_data_first <- get_age_cat(full_data_first)
# full_data_first_combat <- get_age_cat(full_data_first_combat)
# 
# full_data_last <- get_age_cat(full_data_last)
# full_data_last_combat <- get_age_cat(full_data_last_combat)
# 
# data_wt <- get_age_cat(data_wt)


# get age dummy cat
full_data_first <- get_age_cat_dummy(full_data_first)
full_data_first_combat <- get_age_cat_dummy(full_data_first_combat)

full_data_last <- get_age_cat_dummy(full_data_last)
full_data_last_combat <- get_age_cat_dummy(full_data_last_combat)

data_wt <- get_age_cat_dummy(data_wt)
# 
# temp_cancer <- full_data_first[!grepl('Unaffected', full_data_first$cancer_diagnosis_diagnoses),]
# length(which(temp_cancer$age_diagnosis<= 72))
# # 
# data_full <- full_data_first
# wt_data <- data_wt
# beta_thresh = 0.05


full_data <- full_data_first[!grepl('Unaffected', full_data_first$cancer_diagnosis_diagnoses),]
##########
# run model
##########
# https://rdrr.io/rforge/glmmixedlasso/man/glmmlasso.html

image_ids <- as.character(image_ids)

image_ids[image_ids %in% full_data_first$ids]
full_data_first$ids

run_model <- function(data_full,
                      enet,
                      wt_data,
                      bump_type,
                      gender,
                      tech,
                      fam_num,
                      fam_ratio,
                      age_dum,
                      remove_age_cgs_lit,
                      remove_age_cgs_lm,
                      beta_thresh = beta_thresh) {
  
  
  intersect_names <- colnames(data_full)
  clin_data_names <- names(data_full)[!grepl('^cg', colnames(data_full))]
  # read in probes associated with age
  if (remove_age_cgs_lit) {
    age_cgs <- readRDS('../../Data/age_probes.rda')
    intersect_names <- intersect_names[!intersect_names %in% age_cgs]
    data_full <- data_full[, colnames(data_full) %in% intersect_names]
  }
  
  
  # remove age
  if (remove_age_cgs_lm) {
    age_cgs_lm <- readRDS('../../Data/age_cgs_lm.rda')
    intersect_names <- intersect_names[!intersect_names %in% age_cgs_lm]
    data_full <- data_full[, colnames(data_full) %in% intersect_names]
  }
  
  
  # get two data sets from data_full - cases and controls (containing both technologies)
  cases_full <- data_full[data_full$a == 1,]
  controls_full <- data_full[data_full$a == 0,]
  
  remaining_features <- colnames(cases_full)[23:ncol(cases_full)]
  
  # cases_full <- data_full[data_full$cancer_diagnosis_diagnoses != 'Unaffected',]
  # controls_full <- data_full[data_full$cancer_diagnosis_diagnoses == 'Unaffected',]
  
  # if(bump_type == 'both') {
  #   # use cases training and controls to get bumphunter features
  #   bh_feats <- bump_hunter(dat_1 = cases_full, 
  #                           dat_2 = controls_full, 
  #                           wild_type = wt_data,
  #                           bump = 'lfs', 
  #                           boot_num = 5, 
  #                           thresh = beta_thresh,
  #                           g_ranges = g_ranges)
  #   
  #   # get feature list
  #   colnames(bh_feats)[1] <- 'chr'
  #   these_features <- inner_join(bh_feats, g_ranges)$probe
  #   
  #   # subset data by these_features
  #   cases_full <- cases_full[, c(clin_data_names, these_features)]
  #   controls_temp <- controls_full[, c(clin_data_names, these_features)]
  #   
  #   # use cases training and controls to get bumphunter features
  #   bh_feats <- bump_hunter(dat_1 = cases_full, 
  #                           dat_2 = controls_temp,
  #                           wild_type = NULL,
  #                           bump = 'cancer', 
  #                           boot_num = 5, 
  #                           thresh = 0.05,
  #                           g_ranges = g_ranges)
  #   
  #   
  #   # get feature list
  #   colnames(bh_feats)[1] <- 'chr'
  #   remove_features <- inner_join(bh_feats, g_ranges)$probe
  #   
  #   # get features
  #   temp_features <- colnames(cases_full)[23:ncol(cases_full)]
  #   # take remove features out of colnames 
  #   remaining_features <- temp_features[!temp_features %in% remove_features]
  #   
  # } else if(bump_type == 'lfs_only') {
  #   
  #   # read in lfs data 
  #   lfs_features <- readRDS('../../Data/lfs_features_05.rda')
  #   
  #   remaining_features <- lfs_features
  #   
  # } else {
  #   # use cases training and controls to get bumphunter features
  #   bh_feats <- bump_hunter(dat_1 = cases_full, 
  #                           dat_2 = controls_full,
  #                           wild_type = NULL,
  #                           bump = 'cancer', 
  #                           boot_num = 5, 
  #                           thresh = 0.1,
  #                           g_ranges = g_ranges)
  #   
  #   
  #   # get feature list
  #   colnames(bh_feats)[1] <- 'chr'
  #   remove_features <- inner_join(bh_feats, g_ranges)$probe
  #   
  #   # get features
  #   temp_features <- colnames(cases_full)[23:ncol(cases_full)]
  #   # take remove features out of colnames 
  #   remaining_features <- temp_features[!temp_features %in% remove_features]
  # }
  # 
  if(enet){
    num_after_bh <- ncol(cases_full)
    
    # function to predict with all test, controls, controls old, and valid
    mod_result  <- run_enet_test(cases_dat = cases_full,
                                 controls_dat = controls_full,
                                 age_cutoff = 72,
                                 age_dum = age_dum,
                                 gender = gender,
                                 tech = tech,
                                 fam_num = fam_num,
                                 fam_ratio = fam_ratio,
                                 bh_features = remaining_features)
    
    
    temp_controls <- mod_result[[2]]
    temp_model <- mod_result[[1]]
    temp_non_zero_min <- mod_result[[3]]
    temp_non_zero_1se <- mod_result[[4]]
    temp_alpha <- mod_result[[5]]
    
    return(list(temp_controls, temp_model, temp_non_zero_min, temp_non_zero_1se, temp_alpha, num_after_bh))
    
    
  } else {
    num_after_bh <- ncol(cases_full)
    
    # function to predict with all test, controls, controls old, and valid
    mod_result  <- run_rf_test(cases_dat = cases_full,
                               controls_dat = controls_full,
                               age_dum = age_dum,
                               age_cutoff = 72,
                               gender = gender,
                               tech = tech,
                               fam_num = fam_num,
                               fam_ratio = fam_ratio,
                               bh_features = remaining_features)
    
    
    temp_controls <- mod_result[[2]]
    temp_model <- mod_result[[1]]
    temp_importance<- mod_result[[3]]
    temp_data_size <- mod_result[[4]]
    
    return(list(temp_controls, temp_model, temp_importance, temp_data_size))
    
    
  }
  
  
}

##########
# fixed variables
##########


# run full pipeline
full_results <- run_model(full_data,
                          enet = T,
                          wt_data = data_wt,
                          bump_type = 'both',
                          gender = T,
                          tech = F,
                          fam_num = F,
                          fam_ratio = F,
                          age_dum = T,
                          remove_age_cgs_lit = T,
                          remove_age_cgs_lm = T,
                          beta_thresh = 0.05)

# (1) temp_controls, (2) temp_model, (3) temp_non_zero_coef, (4) temp_non_zero_1se, (5) num_after_bh, (6) temp_alpha
# (1) temp_controls, (2) temp_model, (3) temp_importance, (4) temp_data_size

tots <- full_results[[6]]
temp_1se <- full_results[[4]]
temp_min <- full_results[[3]]

temp_controls <- full_results[[1]]
temp_controls <- temp_controls[ , c('controls_age_pred', 'age_sample_collection',
                                    'age_diagnosis')]
# sort by age
temp_controls$pred_label <- ifelse(temp_controls$controls_age_pred >= 0.5, 'yes', 'no')
temp_controls$real_label <- ifelse(temp_controls$age_diagnosis <= 72, 'yes', 'no')

temp_over <- temp_controls[temp_controls$age_sample_collection > 72,]
temp_under <- temp_controls[temp_controls$age_sample_collection <= 72,]

caret::confusionMatrix(temp_under$pred_label, temp_under$real_label)

temp_controls <- temp_controls[order(temp_controls$age_diagnosis, decreasing = TRUE),]

# saveRDS(temp_controls,'~/Desktop/decent_controls.rda')
# funnorm first (0.1): good
# funnorm first (0.5): good 
# outlier from last (0.1) everything marked : 

temp_controls <- read_csv('~/Desktop/temp_controls.csv')

# sort by risk score
temp_controls <- temp_controls[order(temp_controls$controls_age_pred, decreasing = TRUE),]
