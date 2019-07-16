
# old_m_processed_outlier_fun, 
# new_m_processed_outlier, new_m_processed_outliner_more, 
# new_m_processed_outlier_more_funnorm
# validation or combined data 
data_used <- 'old'

# get cg regions
cg_gene_regions = 'Body'

# set preprocessing method
method <- 'noob'

# set type of data, beta or m
methyl_type <- 'm'

# set data directory
data_dir <- '../../Data/'

# get data
load(paste0(data_dir,paste0(data_used,'_',methyl_type, '_final_outlier', '.RData')))

full_data <- rbind(data_cases_full, data_controls_full)
# read in wt data
data_wt <- readRDS(paste0(data_dir,paste0('new','_','beta', '_wild_type', '.rda')))

data_cases_full <- data_cases_full[data_cases_full$tm_donor_ != '3955',]

data_full <- full_data

# # read in wt data
# data_wt <- readRDS(paste0(data_dir,paste0('new','_',methyl_type, '_wild_type', '.rda')))
# 
# full_data <- full_data[full_data$tm_donor_ != '3955',]
# full_data_combat <- full_data[full_data_combat$tm_donor_ != '3955',]

# source all_functions.R to load libraries and my functions
source('all_functions.R')
##########
# run model
##########
# https://rdrr.io/rforge/glmmixedlasso/man/glmmlasso.html

data_full <- full_data
run_model <- function(data_full,
                      enet,
                      wt_data,
                      bump_type,
                      gender = gender,
                      tech = tech,
                      fam_num = fam_num,
                      fam_ratio = fam_ratio,
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
  cases  <- data_full[data_full$a == 1,]
  controls  <- data_full[data_full$a == 0 & full_data$cancer_diagnosis_diagnoses == 'Unaffected',]
  valid  <- data_full[data_full$a == 0 & full_data$cancer_diagnosis_diagnoses != 'Unaffected',]
  
  

  
  if(bump_type == 'both') {
    # use cases training and controls to get bumphunter features
    bh_feats <- bump_hunter(dat_1 = cases, 
                            dat_2 = controls, 
                            wild_type = wt_data,
                            bump = 'lfs', 
                            boot_num = 5, 
                            thresh = beta_thresh,
                            g_ranges = g_ranges)
    
    # get feature list
    colnames(bh_feats)[1] <- 'chr'
    these_features <- inner_join(bh_feats, g_ranges)$probe
    
    # subset data by these_features
    cases <- cases[, c(clin_data_names, these_features)]
    valid  <- valid[, c(clin_data_names, these_features)]
    controls_temp <- controls[, c(clin_data_names, these_features)]
    
    # use cases training and controls to get bumphunter features
    bh_feats <- bump_hunter(dat_1 =cases, 
                            dat_2 = controls_temp,
                            wild_type = NULL,
                            bump = 'cancer', 
                            boot_num = 5, 
                            thresh = 0.1,
                            g_ranges = g_ranges)
    
    
    # get feature list
    colnames(bh_feats)[1] <- 'chr'
    remove_features <- inner_join(bh_feats, g_ranges)$probe
    
    # get features
    temp_features <- colnames(cases)[grepl('^cg', colnames(cases))]
    # take remove features out of colnames 
    remaining_features <- temp_features[!temp_features %in% remove_features]
    
  } else {
    
    # use cases training and controls to get bumphunter features
    bh_feats <- bump_hunter(dat_1 = cases, 
                            dat_2 = controls,
                            wild_type = NULL,
                            bump = 'cancer', 
                            boot_num = 5, 
                            thresh = beta_thresh,
                            g_ranges = g_ranges)
    
    
    # get feature list
    colnames(bh_feats)[1] <- 'chr'
    remove_features <- inner_join(bh_feats, g_ranges)$probe
    
    # get features
    temp_features <- colnames(cases)[grepl('^cg', colnames(cases))]
    # take remove features out of colnames 
    remaining_features <- temp_features[!temp_features %in% remove_features]
    
  }
  if(enet){
    
    # function to predict with all test, controls, controls old, and valid
    mod_result  <- run_enet_all_test(cases_dat = cases,
                                     controls_dat = controls,
                                     valid_dat = valid,
                                     age_cutoff = 72,
                                     gender = gender,
                                     tech = tech,
                                     fam_num = fam_num,
                                     fam_ratio = fam_ratio,
                                     bh_features = remaining_features)
    
    
    temp_model <- mod_result[[1]]
    temp_controls <- mod_result[[2]]
    temp_valid <- mod_result[[3]]
    temp_alpha <- mod_result[[4]]
    
    
    return(list(temp_controls, temp_valid, temp_model, temp_alpha))
    
    
  } else {
    # function to predict with all test, controls, controls old, and valid
    mod_result  <- run_rf_all_test(cases_dat = cases,
                                   controls_dat = controls,
                                   valid_dat = valid,
                                   test_dat = test_cases,
                                   age_cutoff = 72,
                                   gender = gender,
                                   tech = tech,
                                   fam_num = fam_num,
                                   fam_ratio = fam_ratio,
                                   bh_features = remaining_features)
    
    temp_model <- mod_result[[1]]
    temp_controls <- mod_result[[2]]
    temp_valid <- mod_result[[3]]
    temp_importance <- mod_result[[4]]
    
    return(list(temp_controls, temp_valid, temp_importance, temp_importance))
    
    
  }
  
 
}

##########
# fixed variables
##########


# run full pipeline
full_results <- run_model(full_data_combat,
                          enet = T,
                          wt_data = data_wt,
                          bump_type = 'cancer',
                          gender = T,
                          tech = F,
                          fam_num  = F,
                          fam_ratio = F,
                          remove_age_cgs_lit = T,
                          remove_age_cgs_lm = F,
                          beta_thresh = 0.1)


temp_controls <- full_results[[1]]

temp_controls <- temp_controls[ , c('controls_age_pred', 'controls_age_label', 
                                    'age_sample_collection')]

# sort by age
temp_controls <- temp_controls[order(temp_controls$age_sample_collection, decreasing = TRUE),]

# sort by risk score
temp_controls <- temp_controls[order(temp_controls$controls_age_pred, decreasing = TRUE),]


temp_valid <- full_results[[2]]

temp_valid <- temp_valid[ , c('valid_age_pred', 'valide_age_label', 
                                    'age_diagnosis')]

temp_valid$valid_pred_label <- ifelse(temp_valid$valid_age_pred > .5, 1, 0)
temp_valid$valid_onset_label <-  ifelse(temp_valid$age_diagnosis <= 72, 1, 0)

caret::confusionMatrix(temp_valid$valid_pred_label, temp_valid$valid_onset_label)


# sort by age
temp_valid <- temp_valid[order(temp_valid$age_diagnosis, decreasing = TRUE),]

# sort by risk score
temp_valid <- temp_valid[order(temp_valid$valid_age_pred, decreasing = TRUE),]
