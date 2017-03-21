##### This script will run all variations of rf model on batch_sen
#HERE
##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
mod_data_folder <- paste0(data_folder, '/model_data')
results_folder <- paste0(project_folder, '/Results')
quan_folder <- paste0(results_folder, '/quan')
rand_folder <- paste0(results_folder, '/rand')
scripts_folder <- paste0(project_folder, '/Scripts')


# soucce final_model.R
source(paste0(scripts_folder, '/predict_age/final_model.R'))

# read in cases and controls

# read cases
quan_cases_sam <- readRDS(paste0(mod_data_folder, '/quan_cases_sam.rda'))
# read controls
quan_controls <- readRDS(paste0(mod_data_folder, '/quan_controls.rda'))

## batch sam

# load features
load(paste0(mod_data_folder, '/bh_feat.RData'))


# rf, no gender, no resid
quan_rf_sam <- runModels(quan_cases_sam,
                         quan_controls,
                         resid = F,
                         mod = 'rf',
                         gender = F,
                         num_it = 10,
                         DELTA_BETA_THRESH = .10)

# extract resutls
quan_rf_sam_table <- extractResults(quan_rf_sam, 
                                    'quan_rf_sam',
                                    regularize = F)

# rf, gender, no resid
quan_rf_gen_sam <- runModels(quan_cases_sam,
                             quan_controls,
                             resid = F,
                             mod = 'rf',
                             gender = T,
                             num_it = 10,
                             DELTA_BETA_THRESH = .20)

# extract resutls
quan_rf_gen_sam_table <- extractResults(quan_rf_gen_sam, 
                                        'quan_rf_gen_same',
                                        regularize = F)

# 
# rf, no gender, resid
quan_rf_resid_sam <- runModels(temp1,
                                   temp2,
                                   resid = T,
                                   mod = 'rf',
                                   gender = F,
                                   num_it = 5,
                                   DELTA_BETA_THRESH = 0.01)

# extract resutls
quan_rf_resid_sam_table <- extractResults(quan_rf_resid_sam,
                                              'quan_rf_resid_sam',
                                              regularize = F)
# 
# # rf, no gender, no rresid
# quan_rf_gen_resid_sam <- runModels(quan_cases_sam,
#                                        quan_controls,
#                                        resid = T,
#                                        mod = 'rf',
#                                        gender = T,
#                                        num_it = 10,
#                                        DELTA_BETA_THRESH = .20)
# 
# # extract resutls
# quan_rf_gen_resid_sam_table <- extractResults(quan_rf_gen_resid_sam,
#                                                   'quan_rf_gen_resid_sam',
#                                                   regularize = F)

##########
# combine rf table 
##########

rf_table <- rbind(quan_rf_sam_table,
                  quan_rf_gen_sam_table)
                  # quan_rf_resid_sam_table,
                  # quan_rf_gen_resid_sam_table)

saveRDS(rf_table, paste0(quan_folder, '/rf_sam.rda'))

