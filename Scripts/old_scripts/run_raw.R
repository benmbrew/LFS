####### This script will predict age of onset with raw preprocessed data
# this is part of 7th step in pipleline


##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
scripts_folder <- paste0(project_folder, '/Scripts')

##########
# source model_functions.R and run_models.R
##########
source(paste0(scripts_folder, '/predict_age/model_functions_short.R'))
source(paste0(scripts_folder, '/predict_age/run_models_short.R'))

##########
# remove raw, swan, and funnorm objects
##########
rm(list = ls(pattern = "quan_*"))
rm(list = ls(pattern = "swan_*"))
rm(list = ls(pattern = "funnorm_*"))

###################################################################################################################################
## beta_raw

##########
# modesl with full data
##########

#########
# not batch

# Random forest
raw_result_full_rf <- runModels(raw_cases,
                                 model = 'rf',
                                 bump_hunter = F)

raw_table_full_rf <- extractResults(raw_result_full_rf, 
                                     data_name = 'raw_full_rf',
                                     regularize = F)

# Elastic Net
raw_result_full_enet <- runModels(raw_cases,
                                   model = 'enet',
                                   bump_hunter = F)

raw_table_full_enet <- extractResults(raw_result_full_enet, 
                                       data_name = 'raw_full_enet',
                                       regularize = T)

# Lasso
raw_result_full_lasso <- runModels(raw_cases,
                                    model = 'lasso',
                                    bump_hunter = F)

raw_table_full_lasso <- extractResults(raw_result_full_lasso, 
                                        data_name = 'raw_full_lasso',
                                        regularize = T)


#########
# batch

# Random forest
raw_result_batch_full_rf <- runModels(raw_cases_batch,
                                       model = 'rf',
                                       bump_hunter = F)

raw_table_batch_full_rf <- extractResults(raw_result_batch_full_rf, 
                                           data_name = 'raw_full_batch_rf',
                                           regularize = F)

# Elastic Net
raw_result_batch_full_enet <- runModels(raw_cases_batch,
                                         model = 'enet',
                                         bump_hunter = F)

raw_table_batch_full_enet <- extractResults(raw_result_batch_full_enet, 
                                             data_name = 'raw_full_bath_enet',
                                             regularize = T)

# Lasso
raw_result_batch_full_lasso <- runModels(raw_cases_batch,
                                          model = 'lasso',
                                          bump_hunter = F)

raw_table_batch_full_lasso <- extractResults(raw_result_batch_full_lasso, 
                                              data_name = 'raw_full_batch_lasso',
                                              regularize = T)
# get full table
raw_full <- rbind(raw_table_full_rf,
                   raw_table_full_enet,
                   raw_table_full_lasso,
                   raw_table_batch_full_rf,
                   raw_table_batch_full_enet,
                   raw_table_batch_full_lasso)

#save table 
saveRDS(raw_full, 
        file = paste0(raw_folder, '/raw_table_full.rda'))

# remove
rm(raw_table_full_rf,
   raw_table_full_enet,
   raw_table_full_lasso,
   raw_table_batch_full_rf,
   raw_table_batch_full_enet,
   raw_table_batch_full_lasso,
   raw_result_full_rf,
   raw_result_full_enet,
   raw_result_full_lasso,
   raw_result_batch_full_rf,
   raw_result_batch_full_enet,
   raw_result_batch_full_lasso)



########################################################################################################################
# RANDOM FOREST WITH BH
###########################################################################################################################

#################################################################################
# No batch

##################
# Even counts

##########
#normal

#### raw 10 
raw_result_unbatch_even_rf_10 <- runModels(raw_cases,
                                            model = 'rf',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_even_10)

raw_table_unbatch_even_rf_10 <- extractResults(raw_result_unbatch_even_rf_10, 
                                                data_name = 'raw_unbatch_even_rf_10',
                                                regularize = F)


#### raw 20 
raw_result_unbatch_even_rf_20 <- runModels(raw_cases,
                                            model = 'rf',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_even_20)

raw_table_unbatch_even_rf_20 <- extractResults(raw_result_unbatch_even_rf_20, 
                                                data_name = 'raw_unbatch_even_rf_20',
                                                regularize = F)


#### raw 30 
raw_result_unbatch_even_rf_30 <- runModels(raw_cases,
                                            model = 'rf',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_even_30)

raw_table_unbatch_even_rf_30 <- extractResults(raw_result_unbatch_even_rf_30, 
                                                data_name = 'raw_unbatch_even_rf_30',
                                                regularize = F)


#### raw 40 
raw_result_unbatch_even_rf_40 <- runModels(raw_cases,
                                            model = 'rf',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_even_40)

raw_table_unbatch_even_rf_40 <- extractResults(raw_result_unbatch_even_rf_40, 
                                                data_name = 'raw_unbatch_even_rf_40',
                                                regularize = F)


#### raw 10 
raw_result_unbatch_even_rf_50 <- runModels(raw_cases,
                                            model = 'rf',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_even_50)

raw_table_unbatch_even_rf_50 <- extractResults(raw_result_unbatch_even_rf_50, 
                                                data_name = 'raw_unbatch_even_rf_50',
                                                regularize = F)

###########################################################################################################################
# No batch

##################
# Even counts

########
# fwer

raw_result_unbatch_even_rf_fwer_10 <- runModels(raw_cases,
                                                 model = 'rf',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_even_fwer_10)

raw_table_unbatch_even_rf_fwer_10 <- extractResults(raw_result_unbatch_even_rf_fwer_10 , 
                                                     data_name = 'raw_unbatch_rf_even_fwer_10',
                                                     regularize = F)



raw_result_unbatch_even_rf_fwer_20 <- runModels(raw_cases,
                                                 model = 'rf',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_even_fwer_20)

raw_table_unbatch_even_rf_fwer_20 <- extractResults(raw_result_unbatch_even_rf_fwer_20 , 
                                                     data_name = 'raw_unbatch_rf_even_fwer_20',
                                                     regularize = F)


raw_result_unbatch_even_rf_fwer_30 <- runModels(raw_cases,
                                                 model = 'rf',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_even_fwer_30)

raw_table_unbatch_even_rf_fwer_30 <- extractResults(raw_result_unbatch_even_rf_fwer_30 , 
                                                     data_name = 'raw_unbatch_rf_even_fwer_30',
                                                     regularize = F)


raw_result_unbatch_even_rf_fwer_40 <- runModels(raw_cases,
                                                 model = 'rf',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_even_fwer_40)

raw_table_unbatch_even_rf_fwer_40 <- extractResults(raw_result_unbatch_even_rf_fwer_40 , 
                                                     data_name = 'raw_unbatch_rf_even_fwer_40',
                                                     regularize = F)


raw_result_unbatch_even_rf_fwer_50 <- runModels(raw_cases,
                                                 model = 'rf',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_even_fwer_50)

raw_table_unbatch_even_rf_fwer_50 <- extractResults(raw_result_unbatch_even_rf_fwer_50 , 
                                                     data_name = 'raw_unbatch_even_rf_fwer_50',
                                                     regularize = F)


###########################################################################################################################
# No batch

##################
# Even counts

##########
#sig

raw_result_unbatch_even_rf_sig_10 <- runModels(raw_cases,
                                                model = 'rf',
                                                bump_hunter = T, 
                                                bump_hunter_data = raw_even_sig_10)

raw_table_unbatch_even_rf_sig_10 <- extractResults(raw_result_unbatch_even_rf_sig_10 , 
                                                    data_name = 'raw_unbatch_even_rf_sig_10',
                                                    regularize = F)

###########################################################################################################################
# No batch

##################
# Uneven counts

##########
#normal

#### raw 10 
raw_result_unbatch_uneven_rf_10 <- runModels(raw_cases,
                                              model = 'rf',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_uneven_10)

raw_table_unbatch_uneven_rf_10 <- extractResults(raw_result_unbatch_uneven_rf_10, 
                                                  data_name = 'raw_unbatch_uneven_rf_10',
                                                  regularize = F)


#### raw 20 
raw_result_unbatch_uneven_rf_20 <- runModels(raw_cases,
                                              model = 'rf',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_uneven_20)

raw_table_unbatch_uneven_rf_20 <- extractResults(raw_result_unbatch_uneven_rf_20, 
                                                  data_name = 'raw_unbatch_uneven_rf_20',
                                                  regularize = F)


#### raw 30 
raw_result_unbatch_uneven_rf_30 <- runModels(raw_cases,
                                              model = 'rf',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_uneven_30)

raw_table_unbatch_uneven_rf_30 <- extractResults(raw_result_unbatch_uneven_rf_30, 
                                                  data_name = 'raw_unbatch_uneven_rf_30',
                                                  regularize = F)


#### raw 40 
raw_result_unbatch_uneven_rf_40 <- runModels(raw_cases,
                                              model = 'rf',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_uneven_40)

raw_table_unbatch_uneven_rf_40 <- extractResults(raw_result_unbatch_uneven_rf_40, 
                                                  data_name = 'raw_unbatch_uneven_rf_40',
                                                  regularize = F)


#### raw 10 
raw_result_unbatch_uneven_rf_50 <- runModels(raw_cases,
                                              model = 'rf',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_uneven_50)

raw_table_unbatch_uneven_rf_50 <- extractResults(raw_result_unbatch_uneven_rf_50, 
                                                  data_name = 'raw_unbatch_uneven_rf_50',
                                                  regularize = F)

###########################################################################################################################
# No batch

##################
# uneven counts

########
# fwer

raw_result_unbatch_uneven_rf_fwer_10 <- runModels(raw_cases,
                                                   model = 'rf',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_uneven_fwer_10)

raw_table_unbatch_uneven_rf_fwer_10 <- extractResults(raw_result_unbatch_uneven_rf_fwer_10 , 
                                                       data_name = 'raw_unbatch_rf_uneven_fwer_10',
                                                       regularize = F)



raw_result_unbatch_uneven_rf_fwer_20 <- runModels(raw_cases,
                                                   model = 'rf',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_uneven_fwer_20)

raw_table_unbatch_uneven_rf_fwer_20 <- extractResults(raw_result_unbatch_uneven_rf_fwer_20 , 
                                                       data_name = 'raw_unbatch_rf_uneven_fwer_20',
                                                       regularize = F)


raw_result_unbatch_uneven_rf_fwer_30 <- runModels(raw_cases,
                                                   model = 'rf',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_uneven_fwer_30)

raw_table_unbatch_uneven_rf_fwer_30 <- extractResults(raw_result_unbatch_uneven_rf_fwer_30 , 
                                                       data_name = 'raw_unbatch_rf_uneven_fwer_30',
                                                       regularize = F)


raw_result_unbatch_uneven_rf_fwer_40 <- runModels(raw_cases,
                                                   model = 'rf',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_uneven_fwer_40)

raw_table_unbatch_uneven_rf_fwer_40 <- extractResults(raw_result_unbatch_uneven_rf_fwer_40 , 
                                                       data_name = 'raw_unbatch_rf_uneven_fwer_40',
                                                       regularize = F)


raw_result_unbatch_uneven_rf_fwer_50 <- runModels(raw_cases,
                                                   model = 'rf',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_uneven_fwer_50)

raw_table_unbatch_uneven_rf_fwer_50 <- extractResults(raw_result_unbatch_uneven_rf_fwer_50 , 
                                                       data_name = 'raw_unbatch_uneven_rf_fwer_50',
                                                       regularize = F)


###########################################################################################################################
# No batch

##################
# uneven counts

##########
#sig

raw_result_unbatch_uneven_rf_sig_10 <- runModels(raw_cases,
                                                  model = 'rf',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = raw_uneven_sig_10)

raw_table_unbatch_uneven_rf_sig_10 <- extractResults(raw_result_unbatch_uneven_rf_sig_10 , 
                                                      data_name = 'raw_unbatch_rf_uneven_sig_10',
                                                      regularize = F)


# get full table
raw_unbatch_rf <- rbind(raw_table_unbatch_even_rf_10,
                      raw_table_unbatch_even_rf_20,
                      raw_table_unbatch_even_rf_30,
                      raw_table_unbatch_even_rf_40,
                      raw_table_unbatch_even_rf_50,
                      raw_table_unbatch_even_rf_fwer_10,
                      raw_table_unbatch_even_rf_fwer_20,
                      raw_table_unbatch_even_rf_fwer_30,
                      raw_table_unbatch_even_rf_fwer_40,
                      raw_table_unbatch_even_rf_fwer_50,
                      raw_table_unbatch_even_rf_sig_10,
                      raw_table_unbatch_uneven_rf_10,
                      raw_table_unbatch_uneven_rf_20,
                      raw_table_unbatch_uneven_rf_30,
                      raw_table_unbatch_uneven_rf_40,
                      raw_table_unbatch_uneven_rf_50,
                      raw_table_unbatch_uneven_rf_fwer_10,
                      raw_table_unbatch_uneven_rf_fwer_20,
                      raw_table_unbatch_uneven_rf_fwer_30,
                      raw_table_unbatch_uneven_rf_fwer_40,
                      raw_table_unbatch_uneven_rf_fwer_50,
                      raw_table_unbatch_uneven_rf_sig_10)

#save table 
saveRDS(raw_unbatch_rf, 
        file = paste0(raw_folder, '/raw_table_unbatch_rf.rda'))

# remove
rm(raw_table_unbatch_even_rf_10,
   raw_table_unbatch_even_rf_20,
   raw_table_unbatch_even_rf_30,
   raw_table_unbatch_even_rf_40,
   raw_table_unbatch_even_rf_50,
   raw_table_unbatch_even_rf_fwer_10,
   raw_table_unbatch_even_rf_fwer_20,
   raw_table_unbatch_even_rf_fwer_30,
   raw_table_unbatch_even_rf_fwer_40,
   raw_table_unbatch_even_rf_fwer_50,
   raw_table_unbatch_even_rf_sig_10,
   raw_table_unbatch_uneven_rf_10,
   raw_table_unbatch_uneven_rf_20,
   raw_table_unbatch_uneven_rf_30,
   raw_table_unbatch_uneven_rf_40,
   raw_table_unbatch_uneven_rf_50,
   raw_table_unbatch_uneven_rf_fwer_10,
   raw_table_unbatch_uneven_rf_fwer_20,
   raw_table_unbatch_uneven_rf_fwer_30,
   raw_table_unbatch_uneven_rf_fwer_40,
   raw_table_unbatch_uneven_rf_fwer_50,
   raw_table_unbatch_uneven_rf_sig_10)



########################################################################################################################
# RANDOM FOREST WITH BH
###########################################################################################################################

#################################################################################
# batch

##################
# Even counts

##########
#normal

#### raw 10 
raw_result_batch_even_rf_10 <- runModels(raw_cases,
                                          model = 'rf',
                                          bump_hunter = T, 
                                          bump_hunter_data = raw_even_10)

raw_table_batch_even_rf_10 <- extractResults(raw_result_batch_even_rf_10, 
                                              data_name = 'raw_batch_even_rf_10',
                                              regularize = F)


#### raw 20 
raw_result_batch_even_rf_20 <- runModels(raw_cases,
                                          model = 'rf',
                                          bump_hunter = T, 
                                          bump_hunter_data = raw_even_20)

raw_table_batch_even_rf_20 <- extractResults(raw_result_batch_even_rf_20, 
                                              data_name = 'raw_batch_even_rf_20',
                                              regularize = F)


#### raw 30 
raw_result_batch_even_rf_30 <- runModels(raw_cases,
                                          model = 'rf',
                                          bump_hunter = T, 
                                          bump_hunter_data = raw_even_30)

raw_table_batch_even_rf_30 <- extractResults(raw_result_batch_even_rf_30, 
                                              data_name = 'raw_batch_even_rf_30',
                                              regularize = F)


#### raw 40 
raw_result_batch_even_rf_40 <- runModels(raw_cases,
                                          model = 'rf',
                                          bump_hunter = T, 
                                          bump_hunter_data = raw_even_40)

raw_table_batch_even_rf_40 <- extractResults(raw_result_batch_even_rf_40, 
                                              data_name = 'raw_batch_even_rf_40',
                                              regularize = F)


#### raw 10 
raw_result_batch_even_rf_50 <- runModels(raw_cases,
                                          model = 'rf',
                                          bump_hunter = T, 
                                          bump_hunter_data = raw_even_50)

raw_table_batch_even_rf_50 <- extractResults(raw_result_batch_even_rf_50, 
                                              data_name = 'raw_batch_even_rf_50',
                                              regularize = F)

###########################################################################################################################
# batch

##################
# Even counts

########
# fwer

raw_result_batch_even_rf_fwer_10 <- runModels(raw_cases,
                                               model = 'rf',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_even_fwer_10)

raw_table_batch_even_rf_fwer_10 <- extractResults(raw_result_batch_even_rf_fwer_10 , 
                                                   data_name = 'raw_batch_rf_even_fwer_10',
                                                   regularize = F)



raw_result_batch_even_rf_fwer_20 <- runModels(raw_cases,
                                               model = 'rf',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_even_fwer_20)

raw_table_batch_even_rf_fwer_20 <- extractResults(raw_result_batch_even_rf_fwer_20 , 
                                                   data_name = 'raw_batch_rf_even_fwer_20',
                                                   regularize = F)


raw_result_batch_even_rf_fwer_30 <- runModels(raw_cases,
                                               model = 'rf',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_even_fwer_30)

raw_table_batch_even_rf_fwer_30 <- extractResults(raw_result_batch_even_rf_fwer_30 , 
                                                   data_name = 'raw_batch_rf_even_fwer_30',
                                                   regularize = F)


raw_result_batch_even_rf_fwer_40 <- runModels(raw_cases,
                                               model = 'rf',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_even_fwer_40)

raw_table_batch_even_rf_fwer_40 <- extractResults(raw_result_batch_even_rf_fwer_40 , 
                                                   data_name = 'raw_batch_rf_even_fwer_40',
                                                   regularize = F)


raw_result_batch_even_rf_fwer_50 <- runModels(raw_cases,
                                               model = 'rf',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_even_fwer_50)

raw_table_batch_even_rf_fwer_50 <- extractResults(raw_result_batch_even_rf_fwer_50 , 
                                                   data_name = 'raw_batch_even_rf_fwer_50',
                                                   regularize = F)


###########################################################################################################################
# batch

##################
# Even counts

##########
#sig

raw_result_batch_even_rf_sig_10 <- runModels(raw_cases,
                                              model = 'rf',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_even_sig_10)

raw_table_batch_even_rf_sig_10 <- extractResults(raw_result_batch_even_rf_sig_10 , 
                                                  data_name = 'raw_batch_even_rf_sig_10',
                                                  regularize = F)

###########################################################################################################################
# batch

##################
# Uneven counts

##########
#normal

#### raw 10 
raw_result_batch_uneven_rf_10 <- runModels(raw_cases,
                                            model = 'rf',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_uneven_10)

raw_table_batch_uneven_rf_10 <- extractResults(raw_result_batch_uneven_rf_10, 
                                                data_name = 'raw_batch_uneven_rf_10',
                                                regularize = F)


#### raw 20 
raw_result_batch_uneven_rf_20 <- runModels(raw_cases,
                                            model = 'rf',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_uneven_20)

raw_table_batch_uneven_rf_20 <- extractResults(raw_result_batch_uneven_rf_20, 
                                                data_name = 'raw_batch_uneven_rf_20',
                                                regularize = F)


#### raw 30 
raw_result_batch_uneven_rf_30 <- runModels(raw_cases,
                                            model = 'rf',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_uneven_30)

raw_table_batch_uneven_rf_30 <- extractResults(raw_result_batch_uneven_rf_30, 
                                                data_name = 'raw_batch_uneven_rf_30',
                                                regularize = F)


#### raw 40 
raw_result_batch_uneven_rf_40 <- runModels(raw_cases,
                                            model = 'rf',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_uneven_40)

raw_table_batch_uneven_rf_40 <- extractResults(raw_result_batch_uneven_rf_40, 
                                                data_name = 'raw_batch_uneven_rf_40',
                                                regularize = F)


#### raw 10 
raw_result_batch_uneven_rf_50 <- runModels(raw_cases,
                                            model = 'rf',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_uneven_50)

raw_table_batch_uneven_rf_50 <- extractResults(raw_result_batch_uneven_rf_50, 
                                                data_name = 'raw_batch_uneven_rf_50',
                                                regularize = F)

###########################################################################################################################
# batch

##################
# uneven counts

########
# fwer

raw_result_batch_uneven_rf_fwer_10 <- runModels(raw_cases,
                                                 model = 'rf',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_uneven_fwer_10)

raw_table_batch_uneven_rf_fwer_10 <- extractResults(raw_result_batch_uneven_rf_fwer_10 , 
                                                     data_name = 'raw_batch_rf_uneven_fwer_10',
                                                     regularize = F)



raw_result_batch_uneven_rf_fwer_20 <- runModels(raw_cases,
                                                 model = 'rf',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_uneven_fwer_20)

raw_table_batch_uneven_rf_fwer_20 <- extractResults(raw_result_batch_uneven_rf_fwer_20 , 
                                                     data_name = 'raw_batch_rf_uneven_fwer_20',
                                                     regularize = F)


raw_result_batch_uneven_rf_fwer_30 <- runModels(raw_cases,
                                                 model = 'rf',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_uneven_fwer_30)

raw_table_batch_uneven_rf_fwer_30 <- extractResults(raw_result_batch_uneven_rf_fwer_30 , 
                                                     data_name = 'raw_batch_rf_uneven_fwer_30',
                                                     regularize = F)


raw_result_batch_uneven_rf_fwer_40 <- runModels(raw_cases,
                                                 model = 'rf',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_uneven_fwer_40)

raw_table_batch_uneven_rf_fwer_40 <- extractResults(raw_result_batch_uneven_rf_fwer_40 , 
                                                     data_name = 'raw_batch_rf_uneven_fwer_40',
                                                     regularize = F)


raw_result_batch_uneven_rf_fwer_50 <- runModels(raw_cases,
                                                 model = 'rf',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_uneven_fwer_50)

raw_table_batch_uneven_rf_fwer_50 <- extractResults(raw_result_batch_uneven_rf_fwer_50 , 
                                                     data_name = 'raw_batch_uneven_rf_fwer_50',
                                                     regularize = F)


###########################################################################################################################
# batch

##################
# uneven counts

##########
#sig

raw_result_batch_uneven_rf_sig_10 <- runModels(raw_cases,
                                                model = 'rf',
                                                bump_hunter = T, 
                                                bump_hunter_data = raw_uneven_sig_10)

raw_table_batch_uneven_rf_sig_10 <- extractResults(raw_result_batch_uneven_rf_sig_10 , 
                                                    data_name = 'raw_batch_rf_uneven_sig_10',
                                                    regularize = F)


# get full table
raw_batch_rf <- rbind(raw_table_batch_even_rf_10,
                    raw_table_batch_even_rf_20,
                    raw_table_batch_even_rf_30,
                    raw_table_batch_even_rf_40,
                    raw_table_batch_even_rf_50,
                    raw_table_batch_even_rf_fwer_10,
                    raw_table_batch_even_rf_fwer_20,
                    raw_table_batch_even_rf_fwer_30,
                    raw_table_batch_even_rf_fwer_40,
                    raw_table_batch_even_rf_fwer_50,
                    raw_table_batch_even_rf_sig_10,
                    raw_table_batch_uneven_rf_10,
                    raw_table_batch_uneven_rf_20,
                    raw_table_batch_uneven_rf_30,
                    raw_table_batch_uneven_rf_40,
                    raw_table_batch_uneven_rf_50,
                    raw_table_batch_uneven_rf_fwer_10,
                    raw_table_batch_uneven_rf_fwer_20,
                    raw_table_batch_uneven_rf_fwer_30,
                    raw_table_batch_uneven_rf_fwer_40,
                    raw_table_batch_uneven_rf_fwer_50,
                    raw_table_batch_uneven_rf_sig_10)

#save table 
saveRDS(raw_batch_rf, 
        file = paste0(raw_folder, '/raw_table_batch_rf.rda'))

# remove
rm(raw_table_batch_even_rf_10,
   raw_table_batch_even_rf_20,
   raw_table_batch_even_rf_30,
   raw_table_batch_even_rf_40,
   raw_table_batch_even_rf_50,
   raw_table_batch_even_rf_fwer_10,
   raw_table_batch_even_rf_fwer_20,
   raw_table_batch_even_rf_fwer_30,
   raw_table_batch_even_rf_fwer_40,
   raw_table_batch_even_rf_fwer_50,
   raw_table_batch_even_rf_sig_10,
   raw_table_batch_uneven_rf_10,
   raw_table_batch_uneven_rf_20,
   raw_table_batch_uneven_rf_30,
   raw_table_batch_uneven_rf_40,
   raw_table_batch_uneven_rf_50,
   raw_table_batch_uneven_rf_fwer_10,
   raw_table_batch_uneven_rf_fwer_20,
   raw_table_batch_uneven_rf_fwer_30,
   raw_table_batch_uneven_rf_fwer_40,
   raw_table_batch_uneven_rf_fwer_50,
   raw_table_batch_uneven_rf_sig_10)

#######################################################################################################################
#######################################################################################################################

########################################################################################################################
# ENET WITH BH
###########################################################################################################################

#################################################################################
# No batch

##################
# Even counts

##########
#normal

#### raw 10 
raw_result_unbatch_even_enet_10 <- runModels(raw_cases,
                                              model = 'enet',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_even_10)

raw_table_unbatch_even_enet_10 <- extractResults(raw_result_unbatch_even_enet_10, 
                                                  data_name = 'raw_unbatch_even_enet_10',
                                                  regularize = F)


#### raw 20 
raw_result_unbatch_even_enet_20 <- runModels(raw_cases,
                                              model = 'enet',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_even_20)

raw_table_unbatch_even_enet_20 <- extractResults(raw_result_unbatch_even_enet_20, 
                                                  data_name = 'raw_unbatch_even_enet_20',
                                                  regularize = F)


#### raw 30 
raw_result_unbatch_even_enet_30 <- runModels(raw_cases,
                                              model = 'enet',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_even_30)

raw_table_unbatch_even_enet_30 <- extractResults(raw_result_unbatch_even_enet_30, 
                                                  data_name = 'raw_unbatch_even_enet_30',
                                                  regularize = F)


#### raw 40 
raw_result_unbatch_even_enet_40 <- runModels(raw_cases,
                                              model = 'enet',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_even_40)

raw_table_unbatch_even_enet_40 <- extractResults(raw_result_unbatch_even_enet_40, 
                                                  data_name = 'raw_unbatch_even_enet_40',
                                                  regularize = F)


#### raw 10 
raw_result_unbatch_even_enet_50 <- runModels(raw_cases,
                                              model = 'enet',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_even_50)

raw_table_unbatch_even_enet_50 <- extractResults(raw_result_unbatch_even_enet_50, 
                                                  data_name = 'raw_unbatch_even_enet_50',
                                                  regularize = F)

###########################################################################################################################
# No batch

##################
# Even counts

########
# fwer

raw_result_unbatch_even_enet_fwer_10 <- runModels(raw_cases,
                                                   model = 'enet',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_even_fwer_10)

raw_table_unbatch_even_enet_fwer_10 <- extractResults(raw_result_unbatch_even_enet_fwer_10 , 
                                                       data_name = 'raw_unbatch_enet_even_fwer_10',
                                                       regularize = F)



raw_result_unbatch_even_enet_fwer_20 <- runModels(raw_cases,
                                                   model = 'enet',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_even_fwer_20)

raw_table_unbatch_even_enet_fwer_20 <- extractResults(raw_result_unbatch_even_enet_fwer_20 , 
                                                       data_name = 'raw_unbatch_enet_even_fwer_20',
                                                       regularize = F)


raw_result_unbatch_even_enet_fwer_30 <- runModels(raw_cases,
                                                   model = 'enet',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_even_fwer_30)

raw_table_unbatch_even_enet_fwer_30 <- extractResults(raw_result_unbatch_even_enet_fwer_30 , 
                                                       data_name = 'raw_unbatch_enet_even_fwer_30',
                                                       regularize = F)


raw_result_unbatch_even_enet_fwer_40 <- runModels(raw_cases,
                                                   model = 'enet',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_even_fwer_40)

raw_table_unbatch_even_enet_fwer_40 <- extractResults(raw_result_unbatch_even_enet_fwer_40 , 
                                                       data_name = 'raw_unbatch_enet_even_fwer_40',
                                                       regularize = F)


raw_result_unbatch_even_enet_fwer_50 <- runModels(raw_cases,
                                                   model = 'enet',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_even_fwer_50)

raw_table_unbatch_even_enet_fwer_50 <- extractResults(raw_result_unbatch_even_enet_fwer_50 , 
                                                       data_name = 'raw_unbatch_even_enet_fwer_50',
                                                       regularize = F)


###########################################################################################################################
# No batch

##################
# Even counts

##########
#sig

raw_result_unbatch_even_enet_sig_10 <- runModels(raw_cases,
                                                  model = 'enet',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = raw_even_sig_10)

raw_table_unbatch_even_enet_sig_10 <- extractResults(raw_result_unbatch_even_enet_sig_10 , 
                                                      data_name = 'raw_unbatch_even_enet_sig_10',
                                                      regularize = F)

###########################################################################################################################
# No batch

##################
# Uneven counts

##########
#normal

#### raw 10 
raw_result_unbatch_uneven_enet_10 <- runModels(raw_cases,
                                                model = 'enet',
                                                bump_hunter = T, 
                                                bump_hunter_data = raw_uneven_10)

raw_table_unbatch_uneven_enet_10 <- extractResults(raw_result_unbatch_uneven_enet_10, 
                                                    data_name = 'raw_unbatch_uneven_enet_10',
                                                    regularize = F)


#### raw 20 
raw_result_unbatch_uneven_enet_20 <- runModels(raw_cases,
                                                model = 'enet',
                                                bump_hunter = T, 
                                                bump_hunter_data = raw_uneven_20)

raw_table_unbatch_uneven_enet_20 <- extractResults(raw_result_unbatch_uneven_enet_20, 
                                                    data_name = 'raw_unbatch_uneven_enet_20',
                                                    regularize = F)


#### raw 30 
raw_result_unbatch_uneven_enet_30 <- runModels(raw_cases,
                                                model = 'enet',
                                                bump_hunter = T, 
                                                bump_hunter_data = raw_uneven_30)

raw_table_unbatch_uneven_enet_30 <- extractResults(raw_result_unbatch_uneven_enet_30, 
                                                    data_name = 'raw_unbatch_uneven_enet_30',
                                                    regularize = F)


#### raw 40 
raw_result_unbatch_uneven_enet_40 <- runModels(raw_cases,
                                                model = 'enet',
                                                bump_hunter = T, 
                                                bump_hunter_data = raw_uneven_40)

raw_table_unbatch_uneven_enet_40 <- extractResults(raw_result_unbatch_uneven_enet_40, 
                                                    data_name = 'raw_unbatch_uneven_enet_40',
                                                    regularize = F)


#### raw 10 
raw_result_unbatch_uneven_enet_50 <- runModels(raw_cases,
                                                model = 'enet',
                                                bump_hunter = T, 
                                                bump_hunter_data = raw_uneven_50)

raw_table_unbatch_uneven_enet_50 <- extractResults(raw_result_unbatch_uneven_enet_50, 
                                                    data_name = 'raw_unbatch_uneven_enet_50',
                                                    regularize = F)

###########################################################################################################################
# No batch

##################
# uneven counts

########
# fwer

raw_result_unbatch_uneven_enet_fwer_10 <- runModels(raw_cases,
                                                     model = 'enet',
                                                     bump_hunter = T, 
                                                     bump_hunter_data = raw_uneven_fwer_10)

raw_table_unbatch_uneven_enet_fwer_10 <- extractResults(raw_result_unbatch_uneven_enet_fwer_10 , 
                                                         data_name = 'raw_unbatch_enet_uneven_fwer_10',
                                                         regularize = F)



raw_result_unbatch_uneven_enet_fwer_20 <- runModels(raw_cases,
                                                     model = 'enet',
                                                     bump_hunter = T, 
                                                     bump_hunter_data = raw_uneven_fwer_20)

raw_table_unbatch_uneven_enet_fwer_20 <- extractResults(raw_result_unbatch_uneven_enet_fwer_20 , 
                                                         data_name = 'raw_unbatch_enet_uneven_fwer_20',
                                                         regularize = F)


raw_result_unbatch_uneven_enet_fwer_30 <- runModels(raw_cases,
                                                     model = 'enet',
                                                     bump_hunter = T, 
                                                     bump_hunter_data = raw_uneven_fwer_30)

raw_table_unbatch_uneven_enet_fwer_30 <- extractResults(raw_result_unbatch_uneven_enet_fwer_30 , 
                                                         data_name = 'raw_unbatch_enet_uneven_fwer_30',
                                                         regularize = F)


raw_result_unbatch_uneven_enet_fwer_40 <- runModels(raw_cases,
                                                     model = 'enet',
                                                     bump_hunter = T, 
                                                     bump_hunter_data = raw_uneven_fwer_40)

raw_table_unbatch_uneven_enet_fwer_40 <- extractResults(raw_result_unbatch_uneven_enet_fwer_40 , 
                                                         data_name = 'raw_unbatch_enet_uneven_fwer_40',
                                                         regularize = F)


raw_result_unbatch_uneven_enet_fwer_50 <- runModels(raw_cases,
                                                     model = 'enet',
                                                     bump_hunter = T, 
                                                     bump_hunter_data = raw_uneven_fwer_50)

raw_table_unbatch_uneven_enet_fwer_50 <- extractResults(raw_result_unbatch_uneven_enet_fwer_50 , 
                                                         data_name = 'raw_unbatch_uneven_enet_fwer_50',
                                                         regularize = F)


###########################################################################################################################
# No batch

##################
# uneven counts

##########
#sig

raw_result_unbatch_uneven_enet_sig_10 <- runModels(raw_cases,
                                                    model = 'enet',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = raw_uneven_sig_10)

raw_table_unbatch_uneven_enet_sig_10 <- extractResults(raw_result_unbatch_uneven_enet_sig_10 , 
                                                        data_name = 'raw_unbatch_enet_uneven_sig_10',
                                                        regularize = F)


# get full table
raw_unbatch_enet <- rbind(raw_table_unbatch_even_enet_10,
                      raw_table_unbatch_even_enet_20,
                      raw_table_unbatch_even_enet_30,
                      raw_table_unbatch_even_enet_40,
                      raw_table_unbatch_even_enet_50,
                      raw_table_unbatch_even_enet_fwer_10,
                      raw_table_unbatch_even_enet_fwer_20,
                      raw_table_unbatch_even_enet_fwer_30,
                      raw_table_unbatch_even_enet_fwer_40,
                      raw_table_unbatch_even_enet_fwer_50,
                      raw_table_unbatch_even_enet_sig_10,
                      raw_table_unbatch_uneven_enet_10,
                      raw_table_unbatch_uneven_enet_20,
                      raw_table_unbatch_uneven_enet_30,
                      raw_table_unbatch_uneven_enet_40,
                      raw_table_unbatch_uneven_enet_50,
                      raw_table_unbatch_uneven_enet_fwer_10,
                      raw_table_unbatch_uneven_enet_fwer_20,
                      raw_table_unbatch_uneven_enet_fwer_30,
                      raw_table_unbatch_uneven_enet_fwer_40,
                      raw_table_unbatch_uneven_enet_fwer_50,
                      raw_table_unbatch_uneven_enet_sig_10)

#save table 
saveRDS(raw_unbatch_enet, 
        file = paste0(raw_folder, '/raw_table_unbatch_enet.rda'))

# remove
rm(raw_table_unbatch_even_enet_10,
   raw_table_unbatch_even_enet_20,
   raw_table_unbatch_even_enet_30,
   raw_table_unbatch_even_enet_40,
   raw_table_unbatch_even_enet_50,
   raw_table_unbatch_even_enet_fwer_10,
   raw_table_unbatch_even_enet_fwer_20,
   raw_table_unbatch_even_enet_fwer_30,
   raw_table_unbatch_even_enet_fwer_40,
   raw_table_unbatch_even_enet_fwer_50,
   raw_table_unbatch_even_enet_sig_10,
   raw_table_unbatch_uneven_enet_10,
   raw_table_unbatch_uneven_enet_20,
   raw_table_unbatch_uneven_enet_30,
   raw_table_unbatch_uneven_enet_40,
   raw_table_unbatch_uneven_enet_50,
   raw_table_unbatch_uneven_enet_fwer_10,
   raw_table_unbatch_uneven_enet_fwer_20,
   raw_table_unbatch_uneven_enet_fwer_30,
   raw_table_unbatch_uneven_enet_fwer_40,
   raw_table_unbatch_uneven_enet_fwer_50,
   raw_table_unbatch_uneven_enet_sig_10)



########################################################################################################################
# RANDOM FOREST WITH BH
###########################################################################################################################

#################################################################################
# batch

##################
# Even counts

##########
#normal

#### raw 10 
raw_result_batch_even_enet_10 <- runModels(raw_cases,
                                            model = 'enet',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_even_10)

raw_table_batch_even_enet_10 <- extractResults(raw_result_batch_even_enet_10, 
                                                data_name = 'raw_batch_even_enet_10',
                                                regularize = F)


#### raw 20 
raw_result_batch_even_enet_20 <- runModels(raw_cases,
                                            model = 'enet',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_even_20)

raw_table_batch_even_enet_20 <- extractResults(raw_result_batch_even_enet_20, 
                                                data_name = 'raw_batch_even_enet_20',
                                                regularize = F)


#### raw 30 
raw_result_batch_even_enet_30 <- runModels(raw_cases,
                                            model = 'enet',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_even_30)

raw_table_batch_even_enet_30 <- extractResults(raw_result_batch_even_enet_30, 
                                                data_name = 'raw_batch_even_enet_30',
                                                regularize = F)


#### raw 40 
raw_result_batch_even_enet_40 <- runModels(raw_cases,
                                            model = 'enet',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_even_40)

raw_table_batch_even_enet_40 <- extractResults(raw_result_batch_even_enet_40, 
                                                data_name = 'raw_batch_even_enet_40',
                                                regularize = F)


#### raw 10 
raw_result_batch_even_enet_50 <- runModels(raw_cases,
                                            model = 'enet',
                                            bump_hunter = T, 
                                            bump_hunter_data = raw_even_50)

raw_table_batch_even_enet_50 <- extractResults(raw_result_batch_even_enet_50, 
                                                data_name = 'raw_batch_even_enet_50',
                                                regularize = F)

###########################################################################################################################
# batch

##################
# Even counts

########
# fwer

raw_result_batch_even_enet_fwer_10 <- runModels(raw_cases,
                                                 model = 'enet',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_even_fwer_10)

raw_table_batch_even_enet_fwer_10 <- extractResults(raw_result_batch_even_enet_fwer_10 , 
                                                     data_name = 'raw_batch_enet_even_fwer_10',
                                                     regularize = F)



raw_result_batch_even_enet_fwer_20 <- runModels(raw_cases,
                                                 model = 'enet',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_even_fwer_20)

raw_table_batch_even_enet_fwer_20 <- extractResults(raw_result_batch_even_enet_fwer_20 , 
                                                     data_name = 'raw_batch_enet_even_fwer_20',
                                                     regularize = F)


raw_result_batch_even_enet_fwer_30 <- runModels(raw_cases,
                                                 model = 'enet',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_even_fwer_30)

raw_table_batch_even_enet_fwer_30 <- extractResults(raw_result_batch_even_enet_fwer_30 , 
                                                     data_name = 'raw_batch_enet_even_fwer_30',
                                                     regularize = F)


raw_result_batch_even_enet_fwer_40 <- runModels(raw_cases,
                                                 model = 'enet',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_even_fwer_40)

raw_table_batch_even_enet_fwer_40 <- extractResults(raw_result_batch_even_enet_fwer_40 , 
                                                     data_name = 'raw_batch_enet_even_fwer_40',
                                                     regularize = F)


raw_result_batch_even_enet_fwer_50 <- runModels(raw_cases,
                                                 model = 'enet',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_even_fwer_50)

raw_table_batch_even_enet_fwer_50 <- extractResults(raw_result_batch_even_enet_fwer_50 , 
                                                     data_name = 'raw_batch_even_enet_fwer_50',
                                                     regularize = F)


###########################################################################################################################
# batch

##################
# Even counts

##########
#sig

raw_result_batch_even_enet_sig_10 <- runModels(raw_cases,
                                                model = 'enet',
                                                bump_hunter = T, 
                                                bump_hunter_data = raw_even_sig_10)

raw_table_batch_even_enet_sig_10 <- extractResults(raw_result_batch_even_enet_sig_10 , 
                                                    data_name = 'raw_batch_even_enet_sig_10',
                                                    regularize = F)

###########################################################################################################################
# batch

##################
# Uneven counts

##########
#normal

#### raw 10 
raw_result_batch_uneven_enet_10 <- runModels(raw_cases,
                                              model = 'enet',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_uneven_10)

raw_table_batch_uneven_enet_10 <- extractResults(raw_result_batch_uneven_enet_10, 
                                                  data_name = 'raw_batch_uneven_enet_10',
                                                  regularize = F)


#### raw 20 
raw_result_batch_uneven_enet_20 <- runModels(raw_cases,
                                              model = 'enet',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_uneven_20)

raw_table_batch_uneven_enet_20 <- extractResults(raw_result_batch_uneven_enet_20, 
                                                  data_name = 'raw_batch_uneven_enet_20',
                                                  regularize = F)


#### raw 30 
raw_result_batch_uneven_enet_30 <- runModels(raw_cases,
                                              model = 'enet',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_uneven_30)

raw_table_batch_uneven_enet_30 <- extractResults(raw_result_batch_uneven_enet_30, 
                                                  data_name = 'raw_batch_uneven_enet_30',
                                                  regularize = F)


#### raw 40 
raw_result_batch_uneven_enet_40 <- runModels(raw_cases,
                                              model = 'enet',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_uneven_40)

raw_table_batch_uneven_enet_40 <- extractResults(raw_result_batch_uneven_enet_40, 
                                                  data_name = 'raw_batch_uneven_enet_40',
                                                  regularize = F)


#### raw 10 
raw_result_batch_uneven_enet_50 <- runModels(raw_cases,
                                              model = 'enet',
                                              bump_hunter = T, 
                                              bump_hunter_data = raw_uneven_50)

raw_table_batch_uneven_enet_50 <- extractResults(raw_result_batch_uneven_enet_50, 
                                                  data_name = 'raw_batch_uneven_enet_50',
                                                  regularize = F)

###########################################################################################################################
# batch

##################
# uneven counts

########
# fwer

raw_result_batch_uneven_enet_fwer_10 <- runModels(raw_cases,
                                                   model = 'enet',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_uneven_fwer_10)

raw_table_batch_uneven_enet_fwer_10 <- extractResults(raw_result_batch_uneven_enet_fwer_10 , 
                                                       data_name = 'raw_batch_enet_uneven_fwer_10',
                                                       regularize = F)



raw_result_batch_uneven_enet_fwer_20 <- runModels(raw_cases,
                                                   model = 'enet',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_uneven_fwer_20)

raw_table_batch_uneven_enet_fwer_20 <- extractResults(raw_result_batch_uneven_enet_fwer_20 , 
                                                       data_name = 'raw_batch_enet_uneven_fwer_20',
                                                       regularize = F)


raw_result_batch_uneven_enet_fwer_30 <- runModels(raw_cases,
                                                   model = 'enet',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_uneven_fwer_30)

raw_table_batch_uneven_enet_fwer_30 <- extractResults(raw_result_batch_uneven_enet_fwer_30 , 
                                                       data_name = 'raw_batch_enet_uneven_fwer_30',
                                                       regularize = F)


raw_result_batch_uneven_enet_fwer_40 <- runModels(raw_cases,
                                                   model = 'enet',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_uneven_fwer_40)

raw_table_batch_uneven_enet_fwer_40 <- extractResults(raw_result_batch_uneven_enet_fwer_40 , 
                                                       data_name = 'raw_batch_enet_uneven_fwer_40',
                                                       regularize = F)


raw_result_batch_uneven_enet_fwer_50 <- runModels(raw_cases,
                                                   model = 'enet',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_uneven_fwer_50)

raw_table_batch_uneven_enet_fwer_50 <- extractResults(raw_result_batch_uneven_enet_fwer_50 , 
                                                       data_name = 'raw_batch_uneven_enet_fwer_50',
                                                       regularize = F)


###########################################################################################################################
# batch

##################
# uneven counts

##########
#sig

raw_result_batch_uneven_enet_sig_10 <- runModels(raw_cases,
                                                  model = 'enet',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = raw_uneven_sig_10)

raw_table_batch_uneven_enet_sig_10 <- extractResults(raw_result_batch_uneven_enet_sig_10 , 
                                                      data_name = 'raw_batch_enet_uneven_sig_10',
                                                      regularize = F)


# get full table
raw_batch_enet <- rbind(raw_table_batch_even_enet_10,
                    raw_table_batch_even_enet_20,
                    raw_table_batch_even_enet_30,
                    raw_table_batch_even_enet_40,
                    raw_table_batch_even_enet_50,
                    raw_table_batch_even_enet_fwer_10,
                    raw_table_batch_even_enet_fwer_20,
                    raw_table_batch_even_enet_fwer_30,
                    raw_table_batch_even_enet_fwer_40,
                    raw_table_batch_even_enet_fwer_50,
                    raw_table_batch_even_enet_sig_10,
                    raw_table_batch_uneven_enet_10,
                    raw_table_batch_uneven_enet_20,
                    raw_table_batch_uneven_enet_30,
                    raw_table_batch_uneven_enet_40,
                    raw_table_batch_uneven_enet_50,
                    raw_table_batch_uneven_enet_fwer_10,
                    raw_table_batch_uneven_enet_fwer_20,
                    raw_table_batch_uneven_enet_fwer_30,
                    raw_table_batch_uneven_enet_fwer_40,
                    raw_table_batch_uneven_enet_fwer_50,
                    raw_table_batch_uneven_enet_sig_10)

#save table 
saveRDS(raw_batch_enet, 
        file = paste0(raw_folder, '/raw_table_batch_enet.rda'))

# remove
rm(raw_table_batch_even_enet_10,
   raw_table_batch_even_enet_20,
   raw_table_batch_even_enet_30,
   raw_table_batch_even_enet_40,
   raw_table_batch_even_enet_50,
   raw_table_batch_even_enet_fwer_10,
   raw_table_batch_even_enet_fwer_20,
   raw_table_batch_even_enet_fwer_30,
   raw_table_batch_even_enet_fwer_40,
   raw_table_batch_even_enet_fwer_50,
   raw_table_batch_even_enet_sig_10,
   raw_table_batch_uneven_enet_10,
   raw_table_batch_uneven_enet_20,
   raw_table_batch_uneven_enet_30,
   raw_table_batch_uneven_enet_40,
   raw_table_batch_uneven_enet_50,
   raw_table_batch_uneven_enet_fwer_10,
   raw_table_batch_uneven_enet_fwer_20,
   raw_table_batch_uneven_enet_fwer_30,
   raw_table_batch_uneven_enet_fwer_40,
   raw_table_batch_uneven_enet_fwer_50,
   raw_table_batch_uneven_enet_sig_10)


#######################################################################################################################
#######################################################################################################################

########################################################################################################################
# lasso WITH BH
###########################################################################################################################

#################################################################################
# No batch

##################
# Even counts

##########
#normal

#### raw 10 
raw_result_unbatch_even_lasso_10 <- runModels(raw_cases,
                                               model = 'lasso',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_even_10)

raw_table_unbatch_even_lasso_10 <- extractResults(raw_result_unbatch_even_lasso_10, 
                                                   data_name = 'raw_unbatch_even_lasso_10',
                                                   regularize = F)


#### raw 20 
raw_result_unbatch_even_lasso_20 <- runModels(raw_cases,
                                               model = 'lasso',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_even_20)

raw_table_unbatch_even_lasso_20 <- extractResults(raw_result_unbatch_even_lasso_20, 
                                                   data_name = 'raw_unbatch_even_lasso_20',
                                                   regularize = F)


#### raw 30 
raw_result_unbatch_even_lasso_30 <- runModels(raw_cases,
                                               model = 'lasso',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_even_30)

raw_table_unbatch_even_lasso_30 <- extractResults(raw_result_unbatch_even_lasso_30, 
                                                   data_name = 'raw_unbatch_even_lasso_30',
                                                   regularize = F)


#### raw 40 
raw_result_unbatch_even_lasso_40 <- runModels(raw_cases,
                                               model = 'lasso',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_even_40)

raw_table_unbatch_even_lasso_40 <- extractResults(raw_result_unbatch_even_lasso_40, 
                                                   data_name = 'raw_unbatch_even_lasso_40',
                                                   regularize = F)


#### raw 10 
raw_result_unbatch_even_lasso_50 <- runModels(raw_cases,
                                               model = 'lasso',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_even_50)

raw_table_unbatch_even_lasso_50 <- extractResults(raw_result_unbatch_even_lasso_50, 
                                                   data_name = 'raw_unbatch_even_lasso_50',
                                                   regularize = F)

###########################################################################################################################
# No batch

##################
# Even counts

########
# fwer

raw_result_unbatch_even_lasso_fwer_10 <- runModels(raw_cases,
                                                    model = 'lasso',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = raw_even_fwer_10)

raw_table_unbatch_even_lasso_fwer_10 <- extractResults(raw_result_unbatch_even_lasso_fwer_10 , 
                                                        data_name = 'raw_unbatch_lasso_even_fwer_10',
                                                        regularize = F)



raw_result_unbatch_even_lasso_fwer_20 <- runModels(raw_cases,
                                                    model = 'lasso',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = raw_even_fwer_20)

raw_table_unbatch_even_lasso_fwer_20 <- extractResults(raw_result_unbatch_even_lasso_fwer_20 , 
                                                        data_name = 'raw_unbatch_lasso_even_fwer_20',
                                                        regularize = F)


raw_result_unbatch_even_lasso_fwer_30 <- runModels(raw_cases,
                                                    model = 'lasso',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = raw_even_fwer_30)

raw_table_unbatch_even_lasso_fwer_30 <- extractResults(raw_result_unbatch_even_lasso_fwer_30 , 
                                                        data_name = 'raw_unbatch_lasso_even_fwer_30',
                                                        regularize = F)


raw_result_unbatch_even_lasso_fwer_40 <- runModels(raw_cases,
                                                    model = 'lasso',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = raw_even_fwer_40)

raw_table_unbatch_even_lasso_fwer_40 <- extractResults(raw_result_unbatch_even_lasso_fwer_40 , 
                                                        data_name = 'raw_unbatch_lasso_even_fwer_40',
                                                        regularize = F)


raw_result_unbatch_even_lasso_fwer_50 <- runModels(raw_cases,
                                                    model = 'lasso',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = raw_even_fwer_50)

raw_table_unbatch_even_lasso_fwer_50 <- extractResults(raw_result_unbatch_even_lasso_fwer_50 , 
                                                        data_name = 'raw_unbatch_even_lasso_fwer_50',
                                                        regularize = F)


###########################################################################################################################
# No batch

##################
# Even counts

##########
#sig

raw_result_unbatch_even_lasso_sig_10 <- runModels(raw_cases,
                                                   model = 'lasso',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_even_sig_10)

raw_table_unbatch_even_lasso_sig_10 <- extractResults(raw_result_unbatch_even_lasso_sig_10 , 
                                                       data_name = 'raw_unbatch_even_lasso_sig_10',
                                                       regularize = F)

###########################################################################################################################
# No batch

##################
# Uneven counts

##########
#normal

#### raw 10 
raw_result_unbatch_uneven_lasso_10 <- runModels(raw_cases,
                                                 model = 'lasso',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_uneven_10)

raw_table_unbatch_uneven_lasso_10 <- extractResults(raw_result_unbatch_uneven_lasso_10, 
                                                     data_name = 'raw_unbatch_uneven_lasso_10',
                                                     regularize = F)


#### raw 20 
raw_result_unbatch_uneven_lasso_20 <- runModels(raw_cases,
                                                 model = 'lasso',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_uneven_20)

raw_table_unbatch_uneven_lasso_20 <- extractResults(raw_result_unbatch_uneven_lasso_20, 
                                                     data_name = 'raw_unbatch_uneven_lasso_20',
                                                     regularize = F)


#### raw 30 
raw_result_unbatch_uneven_lasso_30 <- runModels(raw_cases,
                                                 model = 'lasso',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_uneven_30)

raw_table_unbatch_uneven_lasso_30 <- extractResults(raw_result_unbatch_uneven_lasso_30, 
                                                     data_name = 'raw_unbatch_uneven_lasso_30',
                                                     regularize = F)


#### raw 40 
raw_result_unbatch_uneven_lasso_40 <- runModels(raw_cases,
                                                 model = 'lasso',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_uneven_40)

raw_table_unbatch_uneven_lasso_40 <- extractResults(raw_result_unbatch_uneven_lasso_40, 
                                                     data_name = 'raw_unbatch_uneven_lasso_40',
                                                     regularize = F)


#### raw 10 
raw_result_unbatch_uneven_lasso_50 <- runModels(raw_cases,
                                                 model = 'lasso',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_uneven_50)

raw_table_unbatch_uneven_lasso_50 <- extractResults(raw_result_unbatch_uneven_lasso_50, 
                                                     data_name = 'raw_unbatch_uneven_lasso_50',
                                                     regularize = F)

###########################################################################################################################
# No batch

##################
# uneven counts

########
# fwer

raw_result_unbatch_uneven_lasso_fwer_10 <- runModels(raw_cases,
                                                      model = 'lasso',
                                                      bump_hunter = T, 
                                                      bump_hunter_data = raw_uneven_fwer_10)

raw_table_unbatch_uneven_lasso_fwer_10 <- extractResults(raw_result_unbatch_uneven_lasso_fwer_10 , 
                                                          data_name = 'raw_unbatch_lasso_uneven_fwer_10',
                                                          regularize = F)



raw_result_unbatch_uneven_lasso_fwer_20 <- runModels(raw_cases,
                                                      model = 'lasso',
                                                      bump_hunter = T, 
                                                      bump_hunter_data = raw_uneven_fwer_20)

raw_table_unbatch_uneven_lasso_fwer_20 <- extractResults(raw_result_unbatch_uneven_lasso_fwer_20 , 
                                                          data_name = 'raw_unbatch_lasso_uneven_fwer_20',
                                                          regularize = F)


raw_result_unbatch_uneven_lasso_fwer_30 <- runModels(raw_cases,
                                                      model = 'lasso',
                                                      bump_hunter = T, 
                                                      bump_hunter_data = raw_uneven_fwer_30)

raw_table_unbatch_uneven_lasso_fwer_30 <- extractResults(raw_result_unbatch_uneven_lasso_fwer_30 , 
                                                          data_name = 'raw_unbatch_lasso_uneven_fwer_30',
                                                          regularize = F)


raw_result_unbatch_uneven_lasso_fwer_40 <- runModels(raw_cases,
                                                      model = 'lasso',
                                                      bump_hunter = T, 
                                                      bump_hunter_data = raw_uneven_fwer_40)

raw_table_unbatch_uneven_lasso_fwer_40 <- extractResults(raw_result_unbatch_uneven_lasso_fwer_40 , 
                                                          data_name = 'raw_unbatch_lasso_uneven_fwer_40',
                                                          regularize = F)


raw_result_unbatch_uneven_lasso_fwer_50 <- runModels(raw_cases,
                                                      model = 'lasso',
                                                      bump_hunter = T, 
                                                      bump_hunter_data = raw_uneven_fwer_50)

raw_table_unbatch_uneven_lasso_fwer_50 <- extractResults(raw_result_unbatch_uneven_lasso_fwer_50 , 
                                                          data_name = 'raw_unbatch_uneven_lasso_fwer_50',
                                                          regularize = F)


###########################################################################################################################
# No batch

##################
# uneven counts

##########
#sig

raw_result_unbatch_uneven_lasso_sig_10 <- runModels(raw_cases,
                                                     model = 'lasso',
                                                     bump_hunter = T, 
                                                     bump_hunter_data = raw_uneven_sig_10)

raw_table_unbatch_uneven_lasso_sig_10 <- extractResults(raw_result_unbatch_uneven_lasso_sig_10 , 
                                                         data_name = 'raw_unbatch_lasso_uneven_sig_10',
                                                         regularize = F)


# get full table
raw_unbatch_lasso <- rbind(raw_table_unbatch_even_lasso_10,
                      raw_table_unbatch_even_lasso_20,
                      raw_table_unbatch_even_lasso_30,
                      raw_table_unbatch_even_lasso_40,
                      raw_table_unbatch_even_lasso_50,
                      raw_table_unbatch_even_lasso_fwer_10,
                      raw_table_unbatch_even_lasso_fwer_20,
                      raw_table_unbatch_even_lasso_fwer_30,
                      raw_table_unbatch_even_lasso_fwer_40,
                      raw_table_unbatch_even_lasso_fwer_50,
                      raw_table_unbatch_even_lasso_sig_10,
                      raw_table_unbatch_uneven_lasso_10,
                      raw_table_unbatch_uneven_lasso_20,
                      raw_table_unbatch_uneven_lasso_30,
                      raw_table_unbatch_uneven_lasso_40,
                      raw_table_unbatch_uneven_lasso_50,
                      raw_table_unbatch_uneven_lasso_fwer_10,
                      raw_table_unbatch_uneven_lasso_fwer_20,
                      raw_table_unbatch_uneven_lasso_fwer_30,
                      raw_table_unbatch_uneven_lasso_fwer_40,
                      raw_table_unbatch_uneven_lasso_fwer_50,
                      raw_table_unbatch_uneven_lasso_sig_10)

#save table 
saveRDS(raw_unbatch_lasso, 
        file = paste0(raw_folder, '/raw_table_unbatch_lasso.rda'))

# remove
rm(raw_table_unbatch_even_lasso_10,
   raw_table_unbatch_even_lasso_20,
   raw_table_unbatch_even_lasso_30,
   raw_table_unbatch_even_lasso_40,
   raw_table_unbatch_even_lasso_50,
   raw_table_unbatch_even_lasso_fwer_10,
   raw_table_unbatch_even_lasso_fwer_20,
   raw_table_unbatch_even_lasso_fwer_30,
   raw_table_unbatch_even_lasso_fwer_40,
   raw_table_unbatch_even_lasso_fwer_50,
   raw_table_unbatch_even_lasso_sig_10,
   raw_table_unbatch_uneven_lasso_10,
   raw_table_unbatch_uneven_lasso_20,
   raw_table_unbatch_uneven_lasso_30,
   raw_table_unbatch_uneven_lasso_40,
   raw_table_unbatch_uneven_lasso_50,
   raw_table_unbatch_uneven_lasso_fwer_10,
   raw_table_unbatch_uneven_lasso_fwer_20,
   raw_table_unbatch_uneven_lasso_fwer_30,
   raw_table_unbatch_uneven_lasso_fwer_40,
   raw_table_unbatch_uneven_lasso_fwer_50,
   raw_table_unbatch_uneven_lasso_sig_10)



########################################################################################################################
# RANDOM FOREST WITH BH
###########################################################################################################################

#################################################################################
# batch

##################
# Even counts

##########
#normal

#### raw 10 
raw_result_batch_even_lasso_10 <- runModels(raw_cases,
                                             model = 'lasso',
                                             bump_hunter = T, 
                                             bump_hunter_data = raw_even_10)

raw_table_batch_even_lasso_10 <- extractResults(raw_result_batch_even_lasso_10, 
                                                 data_name = 'raw_batch_even_lasso_10',
                                                 regularize = F)


#### raw 20 
raw_result_batch_even_lasso_20 <- runModels(raw_cases,
                                             model = 'lasso',
                                             bump_hunter = T, 
                                             bump_hunter_data = raw_even_20)

raw_table_batch_even_lasso_20 <- extractResults(raw_result_batch_even_lasso_20, 
                                                 data_name = 'raw_batch_even_lasso_20',
                                                 regularize = F)


#### raw 30 
raw_result_batch_even_lasso_30 <- runModels(raw_cases,
                                             model = 'lasso',
                                             bump_hunter = T, 
                                             bump_hunter_data = raw_even_30)

raw_table_batch_even_lasso_30 <- extractResults(raw_result_batch_even_lasso_30, 
                                                 data_name = 'raw_batch_even_lasso_30',
                                                 regularize = F)


#### raw 40 
raw_result_batch_even_lasso_40 <- runModels(raw_cases,
                                             model = 'lasso',
                                             bump_hunter = T, 
                                             bump_hunter_data = raw_even_40)

raw_table_batch_even_lasso_40 <- extractResults(raw_result_batch_even_lasso_40, 
                                                 data_name = 'raw_batch_even_lasso_40',
                                                 regularize = F)


#### raw 10 
raw_result_batch_even_lasso_50 <- runModels(raw_cases,
                                             model = 'lasso',
                                             bump_hunter = T, 
                                             bump_hunter_data = raw_even_50)

raw_table_batch_even_lasso_50 <- extractResults(raw_result_batch_even_lasso_50, 
                                                 data_name = 'raw_batch_even_lasso_50',
                                                 regularize = F)

###########################################################################################################################
# batch

##################
# Even counts

########
# fwer

raw_result_batch_even_lasso_fwer_10 <- runModels(raw_cases,
                                                  model = 'lasso',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = raw_even_fwer_10)

raw_table_batch_even_lasso_fwer_10 <- extractResults(raw_result_batch_even_lasso_fwer_10 , 
                                                      data_name = 'raw_batch_lasso_even_fwer_10',
                                                      regularize = F)



raw_result_batch_even_lasso_fwer_20 <- runModels(raw_cases,
                                                  model = 'lasso',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = raw_even_fwer_20)

raw_table_batch_even_lasso_fwer_20 <- extractResults(raw_result_batch_even_lasso_fwer_20 , 
                                                      data_name = 'raw_batch_lasso_even_fwer_20',
                                                      regularize = F)


raw_result_batch_even_lasso_fwer_30 <- runModels(raw_cases,
                                                  model = 'lasso',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = raw_even_fwer_30)

raw_table_batch_even_lasso_fwer_30 <- extractResults(raw_result_batch_even_lasso_fwer_30 , 
                                                      data_name = 'raw_batch_lasso_even_fwer_30',
                                                      regularize = F)


raw_result_batch_even_lasso_fwer_40 <- runModels(raw_cases,
                                                  model = 'lasso',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = raw_even_fwer_40)

raw_table_batch_even_lasso_fwer_40 <- extractResults(raw_result_batch_even_lasso_fwer_40 , 
                                                      data_name = 'raw_batch_lasso_even_fwer_40',
                                                      regularize = F)


raw_result_batch_even_lasso_fwer_50 <- runModels(raw_cases,
                                                  model = 'lasso',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = raw_even_fwer_50)

raw_table_batch_even_lasso_fwer_50 <- extractResults(raw_result_batch_even_lasso_fwer_50 , 
                                                      data_name = 'raw_batch_even_lasso_fwer_50',
                                                      regularize = F)


###########################################################################################################################
# batch

##################
# Even counts

##########
#sig

raw_result_batch_even_lasso_sig_10 <- runModels(raw_cases,
                                                 model = 'lasso',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = raw_even_sig_10)

raw_table_batch_even_lasso_sig_10 <- extractResults(raw_result_batch_even_lasso_sig_10 , 
                                                     data_name = 'raw_batch_even_lasso_sig_10',
                                                     regularize = F)

###########################################################################################################################
# batch

##################
# Uneven counts

##########
#normal

#### raw 10 
raw_result_batch_uneven_lasso_10 <- runModels(raw_cases,
                                               model = 'lasso',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_uneven_10)

raw_table_batch_uneven_lasso_10 <- extractResults(raw_result_batch_uneven_lasso_10, 
                                                   data_name = 'raw_batch_uneven_lasso_10',
                                                   regularize = F)


#### raw 20 
raw_result_batch_uneven_lasso_20 <- runModels(raw_cases,
                                               model = 'lasso',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_uneven_20)

raw_table_batch_uneven_lasso_20 <- extractResults(raw_result_batch_uneven_lasso_20, 
                                                   data_name = 'raw_batch_uneven_lasso_20',
                                                   regularize = F)


#### raw 30 
raw_result_batch_uneven_lasso_30 <- runModels(raw_cases,
                                               model = 'lasso',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_uneven_30)

raw_table_batch_uneven_lasso_30 <- extractResults(raw_result_batch_uneven_lasso_30, 
                                                   data_name = 'raw_batch_uneven_lasso_30',
                                                   regularize = F)


#### raw 40 
raw_result_batch_uneven_lasso_40 <- runModels(raw_cases,
                                               model = 'lasso',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_uneven_40)

raw_table_batch_uneven_lasso_40 <- extractResults(raw_result_batch_uneven_lasso_40, 
                                                   data_name = 'raw_batch_uneven_lasso_40',
                                                   regularize = F)


#### raw 10 
raw_result_batch_uneven_lasso_50 <- runModels(raw_cases,
                                               model = 'lasso',
                                               bump_hunter = T, 
                                               bump_hunter_data = raw_uneven_50)

raw_table_batch_uneven_lasso_50 <- extractResults(raw_result_batch_uneven_lasso_50, 
                                                   data_name = 'raw_batch_uneven_lasso_50',
                                                   regularize = F)

###########################################################################################################################
# batch

##################
# uneven counts

########
# fwer

raw_result_batch_uneven_lasso_fwer_10 <- runModels(raw_cases,
                                                    model = 'lasso',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = raw_uneven_fwer_10)

raw_table_batch_uneven_lasso_fwer_10 <- extractResults(raw_result_batch_uneven_lasso_fwer_10 , 
                                                        data_name = 'raw_batch_lasso_uneven_fwer_10',
                                                        regularize = F)



raw_result_batch_uneven_lasso_fwer_20 <- runModels(raw_cases,
                                                    model = 'lasso',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = raw_uneven_fwer_20)

raw_table_batch_uneven_lasso_fwer_20 <- extractResults(raw_result_batch_uneven_lasso_fwer_20 , 
                                                        data_name = 'raw_batch_lasso_uneven_fwer_20',
                                                        regularize = F)


raw_result_batch_uneven_lasso_fwer_30 <- runModels(raw_cases,
                                                    model = 'lasso',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = raw_uneven_fwer_30)

raw_table_batch_uneven_lasso_fwer_30 <- extractResults(raw_result_batch_uneven_lasso_fwer_30 , 
                                                        data_name = 'raw_batch_lasso_uneven_fwer_30',
                                                        regularize = F)


raw_result_batch_uneven_lasso_fwer_40 <- runModels(raw_cases,
                                                    model = 'lasso',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = raw_uneven_fwer_40)

raw_table_batch_uneven_lasso_fwer_40 <- extractResults(raw_result_batch_uneven_lasso_fwer_40 , 
                                                        data_name = 'raw_batch_lasso_uneven_fwer_40',
                                                        regularize = F)


raw_result_batch_uneven_lasso_fwer_50 <- runModels(raw_cases,
                                                    model = 'lasso',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = raw_uneven_fwer_50)

raw_table_batch_uneven_lasso_fwer_50 <- extractResults(raw_result_batch_uneven_lasso_fwer_50 , 
                                                        data_name = 'raw_batch_uneven_lasso_fwer_50',
                                                        regularize = F)


###########################################################################################################################
# batch

##################
# uneven counts

##########
#sig

raw_result_batch_uneven_lasso_sig_10 <- runModels(raw_cases,
                                                   model = 'lasso',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = raw_uneven_sig_10)

raw_table_batch_uneven_lasso_sig_10 <- extractResults(raw_result_batch_uneven_lasso_sig_10 , 
                                                       data_name = 'raw_batch_lasso_uneven_sig_10',
                                                       regularize = F)


# get full table
raw_batch_lasso <- rbind(raw_table_batch_even_lasso_10,
                    raw_table_batch_even_lasso_20,
                    raw_table_batch_even_lasso_30,
                    raw_table_batch_even_lasso_40,
                    raw_table_batch_even_lasso_50,
                    raw_table_batch_even_lasso_fwer_10,
                    raw_table_batch_even_lasso_fwer_20,
                    raw_table_batch_even_lasso_fwer_30,
                    raw_table_batch_even_lasso_fwer_40,
                    raw_table_batch_even_lasso_fwer_50,
                    raw_table_batch_even_lasso_sig_10,
                    raw_table_batch_uneven_lasso_10,
                    raw_table_batch_uneven_lasso_20,
                    raw_table_batch_uneven_lasso_30,
                    raw_table_batch_uneven_lasso_40,
                    raw_table_batch_uneven_lasso_50,
                    raw_table_batch_uneven_lasso_fwer_10,
                    raw_table_batch_uneven_lasso_fwer_20,
                    raw_table_batch_uneven_lasso_fwer_30,
                    raw_table_batch_uneven_lasso_fwer_40,
                    raw_table_batch_uneven_lasso_fwer_50,
                    raw_table_batch_uneven_lasso_sig_10)

#save table 
saveRDS(raw_batch_lasso, 
        file = paste0(raw_folder, '/raw_table_batch_lasso.rda'))

# remove
rm(raw_table_batch_even_lasso_10,
   raw_table_batch_even_lasso_20,
   raw_table_batch_even_lasso_30,
   raw_table_batch_even_lasso_40,
   raw_table_batch_even_lasso_50,
   raw_table_batch_even_lasso_fwer_10,
   raw_table_batch_even_lasso_fwer_20,
   raw_table_batch_even_lasso_fwer_30,
   raw_table_batch_even_lasso_fwer_40,
   raw_table_batch_even_lasso_fwer_50,
   raw_table_batch_even_lasso_sig_10,
   raw_table_batch_uneven_lasso_10,
   raw_table_batch_uneven_lasso_20,
   raw_table_batch_uneven_lasso_30,
   raw_table_batch_uneven_lasso_40,
   raw_table_batch_uneven_lasso_50,
   raw_table_batch_uneven_lasso_fwer_10,
   raw_table_batch_uneven_lasso_fwer_20,
   raw_table_batch_uneven_lasso_fwer_30,
   raw_table_batch_uneven_lasso_fwer_40,
   raw_table_batch_uneven_lasso_fwer_50,
   raw_table_batch_uneven_lasso_sig_10)


