##### This script will run all variations of enet model on batch_sen

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
quan_cases_sen <- readRDS(paste0(mod_data_folder, '/quan_cases_sen.rda'))
# read controls
quan_controls <- readRDS(paste0(mod_data_folder, '/quan_controls.rda'))


## batch sen

# enet, no gender, no resid
quan_enet_sen <- runModels(quan_cases_sen,
                         quan_controls,
                         resid = F,
                         mod = 'enet',
                         gender = F,
                         num_it = 10,
                         DELTA_BETA_THRESH = .20)

# extract resutls
quan_enet_sen_table <- extractResults(quan_enet_sen, 
                                    'quan_enet_sen',
                                    regularize = T)

# enet, gender, no resid
quan_enet_gen_sen <- runModels(quan_cases_sen,
                             quan_controls,
                             resid = F,
                             mod = 'enet',
                             gender = T,
                             num_it = 10,
                             DELTA_BETA_THRESH = .20)

# extract resutls
quan_enet_gen_sen_table <- extractResults(quan_enet_gen_sen, 
                                        'quan_enet_gen_sen',
                                        regularize = T)

# # enet, no gender, resid
# quan_enet_resid_sen <- runModels(quan_cases_sen,
#                                quan_controls,
#                                resid = T,
#                                mod = 'enet',
#                                gender = F,
#                                num_it = 10,
#                                DELTA_BETA_THRESH = .20)
# 
# # extract resutls
# quan_enet_resid_sen_table <- extractResults(quan_enet_resid_sen, 
#                                           'quan_enet_resid_sen',
#                                           regularize = T)
# 
# # enet, no gender, no rresid
# quan_enet_gen_resid_sen <- runModels(quan_cases_sen,
#                                    quan_controls,
#                                    resid = T,
#                                    mod = 'enet',
#                                    gender = T,
#                                    num_it = 10,
#                                    DELTA_BETA_THRESH = .20)
# 
# # extract resutls
# quan_enet_gen_resid_sen_table <- extractResults(quan_enet_gen_resid_sen, 
#                                               'quan_enet_gen_resid_sen',
#                                               regularize = T)
# 

##########
# combine enet table 
##########

enet_table <- rbind(quan_enet_sen_table,
                  quan_enet_gen_sen_table)
                  # quan_enet_resid_sen_table,
                  # quan_enet_gen_resid_sen_table)

saveRDS(enet_table, paste0(quan_folder, '/enet_sen.rda'))

