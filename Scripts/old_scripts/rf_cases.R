##### This script will run all variations of rf model on cases 

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
quan_cases <- readRDS(paste0(mod_data_folder, '/quan_cases.rda'))
# read controls
quan_controls <- readRDS(paste0(mod_data_folder, '/quan_controls.rda'))


## no batch

# rf, no gender, no resid
quan_rf <- runModels(quan_cases,
                     quan_controls,
                     resid = F,
                     mod = "rf",
                     gender = F,
                     num_it = 3,
                     DELTA_BETA_THRESH = .20)

# extract resutls
quan_rf_table <- extractResults(quan_rf, 
                                'quan_rf',
                                regularize = F)


# rf, gender, no resid
quan_rf_gen <- runModels(quan_cases,
                         quan_controls,
                         resid = F,
                         mod = 'rf',
                         gender = T,
                         num_it = 10,
                         DELTA_BETA_THRESH = .20)

# extract resutls
quan_rf_gen_table <- extractResults(quan_rf_gen, 
                                    'quan_rf_gen',
                                    regularize = F)

# rf, no gender, resid
quan_rf_resid <- runModels(quan_cases,
                           quan_controls,
                           resid = T,
                           mod = 'rf',
                           gender = F,
                           num_it = 10,
                           DELTA_BETA_THRESH = .20)

# extract resutls
quan_rf_resid_table <- extractResults(quan_rf_resid,
                                      'quan_resid',
                                      regularize = F)

# # rf, no gender, no rresid
# quan_rf_gen_resid <- runModels(quan_cases,
#                                quan_controls,
#                                resid = T,
#                                mod = 'rf',
#                                gender = T,
#                                num_it = 10,
#                                DELTA_BETA_THRESH = .20)
# 
# # extract resutls
# quan_rf_gen_resid_table <- extractResults(quan_rf_gen_resid, 
#                                           'quan_rf_gen_resid',
#                                           regularize = F)


##########
# combine rf table 
##########

rf_table <- rbind(quan_rf_table,
                  quan_rf_gen_table)
                  # quan_rf_resid_table,
                  # quan_rf_gen_resid_table)


saveRDS(rf_table, paste0(quan_folder, '/rf_cases.rda'))



