# old_m_processed_outlier_fun, 
# new_m_processed_outlier, new_m_processed_outliner_more, 
# new_m_processed_outlier_more_funnorm

# this script will get combat data on 2 tech types
method <- 'noob'

# set type of data, beta or m
methyl_type <- 'beta'

# get the type of data
data_used <- 'new'

# source all_functions.R to load libraries and my functions
source('all_functions.R')

# set data directory
data_dir <- '../../Data/'

# get data
if(paste0(data_used,'_',methyl_type, '_processed_beta_first_last_gen', '.RData') %in% dir(data_dir)) {
  load(paste0(data_dir,paste0(data_used,'_',methyl_type, '_processed_beta_first_last', '.RData')))
} else {
  source(paste0('get_data_', data_used , '.R'))
}

if(data_used == 'new') {
  # combine data
  full_data_last <- rbind(data_cases_full_last,
                     data_controls_full)
  
  full_data_first <- rbind(data_cases_full_first,
                          data_controls_full)
  
  rm(data_cases_full_first, data_cases_full_last, data_controls_full)
} else {
  # remove 'a' and 'b' from columns
  
  # combine data
  full_data <- rbind(data_cases,
                     data_controls_mod,
                     data_valid_mod)
  rm(data_cases, data_controls_mod, data_valid_mod)
}

# run combat on technology
full_data_first_combat <- run_combat(full_data_first)
full_data_last_combat <- run_combat(full_data_last)


# remove tech variable from full_data_combat
full_data_first_combat$tech <- NULL
full_data_last_combat$tech <- NULL


save.image(paste0(data_dir,paste0(data_used,'_',methyl_type, '_final_beta_first_last_gen', '.RData')))


summary(as.factor(full_data_first$a))
temp_cancer <- full_data_first[!grepl('Unaffected', full_data_first$cancer_diagnosis_diagnoses),]
temp_controls <- full_data_first[grepl('Unaffected', full_data_first$cancer_diagnosis_diagnoses),]

temp_controls_first <- temp_controls[!duplicated(temp_controls$tm_donor_, fromLast = F),]
temp_controls_last <- temp_controls[!duplicated(temp_controls$tm_donor_, fromLast = T),]


