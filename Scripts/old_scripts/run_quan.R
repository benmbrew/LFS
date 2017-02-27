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
rm(list = ls(pattern = "raw_*"))
rm(list = ls(pattern = "swan_*"))
rm(list = ls(pattern = "funnorm_*"))

###################################################################################################################################
## beta_quan

##########
# modesl with full data
##########

#########
# not batch

# Random forest
quan_result_full_rf <- runModels(quan_cases,
                                model = 'rf',
                                bump_hunter = F)

quan_table_full_rf <- extractResults(quan_result_full_rf, 
                                    data_name = 'quan_full_rf',
                                    regularize = F)

# Elastic Net
quan_result_full_enet <- runModels(quan_cases,
                                  model = 'enet',
                                  bump_hunter = F)

quan_table_full_enet <- extractResults(quan_result_full_enet, 
                                      data_name = 'quan_full_enet',
                                      regularize = T)

# Lasso
quan_result_full_lasso <- runModels(quan_cases,
                                   model = 'lasso',
                                   bump_hunter = F)

quan_table_full_lasso <- extractResults(quan_result_full_lasso, 
                                       data_name = 'quan_full_lasso',
                                       regularize = T)


#########
# batch

# Random forest
quan_result_batch_full_rf <- runModels(quan_cases_batch,
                                      model = 'rf',
                                      bump_hunter = F)

quan_table_batch_full_rf <- extractResults(quan_result_batch_full_rf, 
                                          data_name = 'quan_full_batch_rf',
                                          regularize = F)

# Elastic Net
quan_result_batch_full_enet <- runModels(quan_cases_batch,
                                        model = 'enet',
                                        bump_hunter = F)

quan_table_batch_full_enet <- extractResults(quan_result_batch_full_enet, 
                                            data_name = 'quan_full_bath_enet',
                                            regularize = T)

# Lasso
quan_result_batch_full_lasso <- runModels(quan_cases_batch,
                                         model = 'lasso',
                                         bump_hunter = F)

quan_table_batch_full_lasso <- extractResults(quan_result_batch_full_lasso, 
                                             data_name = 'quan_full_batch_lasso',
                                             regularize = T)
# get full table
quan_full <- rbind(quan_table_full_rf,
                  quan_table_full_enet,
                  quan_table_full_lasso,
                  quan_table_batch_full_rf,
                  quan_table_batch_full_enet,
                  quan_table_batch_full_lasso)

#save table 
saveRDS(quan_full, 
        file = paste0(quan_folder, '/quan_table_full.rda'))

# remove
rm(quan_table_full_rf,
   quan_table_full_enet,
   quan_table_full_lasso,
   quan_table_batch_full_rf,
   quan_table_batch_full_enet,
   quan_table_batch_full_lasso,
   quan_result_full_rf,
   quan_result_full_enet,
   quan_result_full_lasso,
   quan_result_batch_full_rf,
   quan_result_batch_full_enet,
   quan_result_batch_full_lasso)



########################################################################################################################
# RANDOM FOREST WITH BH
###########################################################################################################################

#################################################################################
# No batch

##################
# Even counts

##########
#normal

#### quan 10 
quan_result_unbatch_even_rf_10 <- runModels(quan_cases,
                                           model = 'rf',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_even_10)

quan_table_unbatch_even_rf_10 <- extractResults(quan_result_unbatch_even_rf_10, 
                                               data_name = 'quan_unbatch_even_rf_10',
                                               regularize = F)


#### quan 20 
quan_result_unbatch_even_rf_20 <- runModels(quan_cases,
                                           model = 'rf',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_even_20)

quan_table_unbatch_even_rf_20 <- extractResults(quan_result_unbatch_even_rf_20, 
                                               data_name = 'quan_unbatch_even_rf_20',
                                               regularize = F)


#### quan 30 
quan_result_unbatch_even_rf_30 <- runModels(quan_cases,
                                           model = 'rf',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_even_30)

quan_table_unbatch_even_rf_30 <- extractResults(quan_result_unbatch_even_rf_30, 
                                               data_name = 'quan_unbatch_even_rf_30',
                                               regularize = F)


#### quan 40 
quan_result_unbatch_even_rf_40 <- runModels(quan_cases,
                                           model = 'rf',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_even_40)

quan_table_unbatch_even_rf_40 <- extractResults(quan_result_unbatch_even_rf_40, 
                                               data_name = 'quan_unbatch_even_rf_40',
                                               regularize = F)


#### quan 10 
quan_result_unbatch_even_rf_50 <- runModels(quan_cases,
                                           model = 'rf',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_even_50)

quan_table_unbatch_even_rf_50 <- extractResults(quan_result_unbatch_even_rf_50, 
                                               data_name = 'quan_unbatch_even_rf_50',
                                               regularize = F)

###########################################################################################################################
# No batch

##################
# Even counts

########
# fwer

quan_result_unbatch_even_rf_fwer_10 <- runModels(quan_cases,
                                                model = 'rf',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_even_fwer_10)

quan_table_unbatch_even_rf_fwer_10 <- extractResults(quan_result_unbatch_even_rf_fwer_10 , 
                                                    data_name = 'quan_unbatch_rf_even_fwer_10',
                                                    regularize = F)



quan_result_unbatch_even_rf_fwer_20 <- runModels(quan_cases,
                                                model = 'rf',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_even_fwer_20)

quan_table_unbatch_even_rf_fwer_20 <- extractResults(quan_result_unbatch_even_rf_fwer_20 , 
                                                    data_name = 'quan_unbatch_rf_even_fwer_20',
                                                    regularize = F)


quan_result_unbatch_even_rf_fwer_30 <- runModels(quan_cases,
                                                model = 'rf',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_even_fwer_30)

quan_table_unbatch_even_rf_fwer_30 <- extractResults(quan_result_unbatch_even_rf_fwer_30 , 
                                                    data_name = 'quan_unbatch_rf_even_fwer_30',
                                                    regularize = F)


quan_result_unbatch_even_rf_fwer_40 <- runModels(quan_cases,
                                                model = 'rf',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_even_fwer_40)

quan_table_unbatch_even_rf_fwer_40 <- extractResults(quan_result_unbatch_even_rf_fwer_40 , 
                                                    data_name = 'quan_unbatch_rf_even_fwer_40',
                                                    regularize = F)


quan_result_unbatch_even_rf_fwer_50 <- runModels(quan_cases,
                                                model = 'rf',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_even_fwer_50)

quan_table_unbatch_even_rf_fwer_50 <- extractResults(quan_result_unbatch_even_rf_fwer_50 , 
                                                    data_name = 'quan_unbatch_even_rf_fwer_50',
                                                    regularize = F)


###########################################################################################################################
# No batch

##################
# Even counts

##########
#sig

quan_result_unbatch_even_rf_sig_10 <- runModels(quan_cases,
                                               model = 'rf',
                                               bump_hunter = T, 
                                               bump_hunter_data = quan_even_sig_10)

quan_table_unbatch_even_rf_sig_10 <- extractResults(quan_result_unbatch_even_rf_sig_10 , 
                                                   data_name = 'quan_unbatch_even_rf_sig_10',
                                                   regularize = F)

###########################################################################################################################
# No batch

##################
# Uneven counts

##########
#normal

#### quan 10 
quan_result_unbatch_uneven_rf_10 <- runModels(quan_cases,
                                             model = 'rf',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_uneven_10)

quan_table_unbatch_uneven_rf_10 <- extractResults(quan_result_unbatch_uneven_rf_10, 
                                                 data_name = 'quan_unbatch_uneven_rf_10',
                                                 regularize = F)


#### quan 20 
quan_result_unbatch_uneven_rf_20 <- runModels(quan_cases,
                                             model = 'rf',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_uneven_20)

quan_table_unbatch_uneven_rf_20 <- extractResults(quan_result_unbatch_uneven_rf_20, 
                                                 data_name = 'quan_unbatch_uneven_rf_20',
                                                 regularize = F)


#### quan 30 
quan_result_unbatch_uneven_rf_30 <- runModels(quan_cases,
                                             model = 'rf',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_uneven_30)

quan_table_unbatch_uneven_rf_30 <- extractResults(quan_result_unbatch_uneven_rf_30, 
                                                 data_name = 'quan_unbatch_uneven_rf_30',
                                                 regularize = F)


#### quan 40 
quan_result_unbatch_uneven_rf_40 <- runModels(quan_cases,
                                             model = 'rf',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_uneven_40)

quan_table_unbatch_uneven_rf_40 <- extractResults(quan_result_unbatch_uneven_rf_40, 
                                                 data_name = 'quan_unbatch_uneven_rf_40',
                                                 regularize = F)


#### quan 10 
quan_result_unbatch_uneven_rf_50 <- runModels(quan_cases,
                                             model = 'rf',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_uneven_50)

quan_table_unbatch_uneven_rf_50 <- extractResults(quan_result_unbatch_uneven_rf_50, 
                                                 data_name = 'quan_unbatch_uneven_rf_50',
                                                 regularize = F)

###########################################################################################################################
# No batch

##################
# uneven counts

########
# fwer

quan_result_unbatch_uneven_rf_fwer_10 <- runModels(quan_cases,
                                                  model = 'rf',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_uneven_fwer_10)

quan_table_unbatch_uneven_rf_fwer_10 <- extractResults(quan_result_unbatch_uneven_rf_fwer_10 , 
                                                      data_name = 'quan_unbatch_rf_uneven_fwer_10',
                                                      regularize = F)



quan_result_unbatch_uneven_rf_fwer_20 <- runModels(quan_cases,
                                                  model = 'rf',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_uneven_fwer_20)

quan_table_unbatch_uneven_rf_fwer_20 <- extractResults(quan_result_unbatch_uneven_rf_fwer_20 , 
                                                      data_name = 'quan_unbatch_rf_uneven_fwer_20',
                                                      regularize = F)


quan_result_unbatch_uneven_rf_fwer_30 <- runModels(quan_cases,
                                                  model = 'rf',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_uneven_fwer_30)

quan_table_unbatch_uneven_rf_fwer_30 <- extractResults(quan_result_unbatch_uneven_rf_fwer_30 , 
                                                      data_name = 'quan_unbatch_rf_uneven_fwer_30',
                                                      regularize = F)


quan_result_unbatch_uneven_rf_fwer_40 <- runModels(quan_cases,
                                                  model = 'rf',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_uneven_fwer_40)

quan_table_unbatch_uneven_rf_fwer_40 <- extractResults(quan_result_unbatch_uneven_rf_fwer_40 , 
                                                      data_name = 'quan_unbatch_rf_uneven_fwer_40',
                                                      regularize = F)


quan_result_unbatch_uneven_rf_fwer_50 <- runModels(quan_cases,
                                                  model = 'rf',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_uneven_fwer_50)

quan_table_unbatch_uneven_rf_fwer_50 <- extractResults(quan_result_unbatch_uneven_rf_fwer_50 , 
                                                      data_name = 'quan_unbatch_uneven_rf_fwer_50',
                                                      regularize = F)


###########################################################################################################################
# No batch

##################
# uneven counts

##########
#sig

quan_result_unbatch_uneven_rf_sig_10 <- runModels(quan_cases,
                                                 model = 'rf',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = quan_uneven_sig_10)

quan_table_unbatch_uneven_rf_sig_10 <- extractResults(quan_result_unbatch_uneven_rf_sig_10 , 
                                                     data_name = 'quan_unbatch_rf_uneven_sig_10',
                                                     regularize = F)


# get full table
quan_unbatch_rf <- rbind(quan_table_unbatch_even_rf_10,
                        quan_table_unbatch_even_rf_20,
                        quan_table_unbatch_even_rf_30,
                        quan_table_unbatch_even_rf_40,
                        quan_table_unbatch_even_rf_50,
                        quan_table_unbatch_even_rf_fwer_10,
                        quan_table_unbatch_even_rf_fwer_20,
                        quan_table_unbatch_even_rf_fwer_30,
                        quan_table_unbatch_even_rf_fwer_40,
                        quan_table_unbatch_even_rf_fwer_50,
                        quan_table_unbatch_even_rf_sig_10,
                        quan_table_unbatch_uneven_rf_10,
                        quan_table_unbatch_uneven_rf_20,
                        quan_table_unbatch_uneven_rf_30,
                        quan_table_unbatch_uneven_rf_40,
                        quan_table_unbatch_uneven_rf_50,
                        quan_table_unbatch_uneven_rf_fwer_10,
                        quan_table_unbatch_uneven_rf_fwer_20,
                        quan_table_unbatch_uneven_rf_fwer_30,
                        quan_table_unbatch_uneven_rf_fwer_40,
                        quan_table_unbatch_uneven_rf_fwer_50,
                        quan_table_unbatch_uneven_rf_sig_10)

#save table 
saveRDS(quan_unbatch_rf, 
        file = paste0(quan_folder, '/quan_table_unbatch_rf.rda'))

# remove
rm(quan_table_unbatch_even_rf_10,
   quan_table_unbatch_even_rf_20,
   quan_table_unbatch_even_rf_30,
   quan_table_unbatch_even_rf_40,
   quan_table_unbatch_even_rf_50,
   quan_table_unbatch_even_rf_fwer_10,
   quan_table_unbatch_even_rf_fwer_20,
   quan_table_unbatch_even_rf_fwer_30,
   quan_table_unbatch_even_rf_fwer_40,
   quan_table_unbatch_even_rf_fwer_50,
   quan_table_unbatch_even_rf_sig_10,
   quan_table_unbatch_uneven_rf_10,
   quan_table_unbatch_uneven_rf_20,
   quan_table_unbatch_uneven_rf_30,
   quan_table_unbatch_uneven_rf_40,
   quan_table_unbatch_uneven_rf_50,
   quan_table_unbatch_uneven_rf_fwer_10,
   quan_table_unbatch_uneven_rf_fwer_20,
   quan_table_unbatch_uneven_rf_fwer_30,
   quan_table_unbatch_uneven_rf_fwer_40,
   quan_table_unbatch_uneven_rf_fwer_50,
   quan_table_unbatch_uneven_rf_sig_10)



########################################################################################################################
# RANDOM FOREST WITH BH
###########################################################################################################################

#################################################################################
# batch

##################
# Even counts

##########
#normal

#### quan 10 
quan_result_batch_even_rf_10 <- runModels(quan_cases,
                                         model = 'rf',
                                         bump_hunter = T, 
                                         bump_hunter_data = quan_even_10)

quan_table_batch_even_rf_10 <- extractResults(quan_result_batch_even_rf_10, 
                                             data_name = 'quan_batch_even_rf_10',
                                             regularize = F)


#### quan 20 
quan_result_batch_even_rf_20 <- runModels(quan_cases,
                                         model = 'rf',
                                         bump_hunter = T, 
                                         bump_hunter_data = quan_even_20)

quan_table_batch_even_rf_20 <- extractResults(quan_result_batch_even_rf_20, 
                                             data_name = 'quan_batch_even_rf_20',
                                             regularize = F)


#### quan 30 
quan_result_batch_even_rf_30 <- runModels(quan_cases,
                                         model = 'rf',
                                         bump_hunter = T, 
                                         bump_hunter_data = quan_even_30)

quan_table_batch_even_rf_30 <- extractResults(quan_result_batch_even_rf_30, 
                                             data_name = 'quan_batch_even_rf_30',
                                             regularize = F)


#### quan 40 
quan_result_batch_even_rf_40 <- runModels(quan_cases,
                                         model = 'rf',
                                         bump_hunter = T, 
                                         bump_hunter_data = quan_even_40)

quan_table_batch_even_rf_40 <- extractResults(quan_result_batch_even_rf_40, 
                                             data_name = 'quan_batch_even_rf_40',
                                             regularize = F)


#### quan 10 
quan_result_batch_even_rf_50 <- runModels(quan_cases,
                                         model = 'rf',
                                         bump_hunter = T, 
                                         bump_hunter_data = quan_even_50)

quan_table_batch_even_rf_50 <- extractResults(quan_result_batch_even_rf_50, 
                                             data_name = 'quan_batch_even_rf_50',
                                             regularize = F)

###########################################################################################################################
# batch

##################
# Even counts

########
# fwer

quan_result_batch_even_rf_fwer_10 <- runModels(quan_cases,
                                              model = 'rf',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_even_fwer_10)

quan_table_batch_even_rf_fwer_10 <- extractResults(quan_result_batch_even_rf_fwer_10 , 
                                                  data_name = 'quan_batch_rf_even_fwer_10',
                                                  regularize = F)



quan_result_batch_even_rf_fwer_20 <- runModels(quan_cases,
                                              model = 'rf',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_even_fwer_20)

quan_table_batch_even_rf_fwer_20 <- extractResults(quan_result_batch_even_rf_fwer_20 , 
                                                  data_name = 'quan_batch_rf_even_fwer_20',
                                                  regularize = F)


quan_result_batch_even_rf_fwer_30 <- runModels(quan_cases,
                                              model = 'rf',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_even_fwer_30)

quan_table_batch_even_rf_fwer_30 <- extractResults(quan_result_batch_even_rf_fwer_30 , 
                                                  data_name = 'quan_batch_rf_even_fwer_30',
                                                  regularize = F)


quan_result_batch_even_rf_fwer_40 <- runModels(quan_cases,
                                              model = 'rf',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_even_fwer_40)

quan_table_batch_even_rf_fwer_40 <- extractResults(quan_result_batch_even_rf_fwer_40 , 
                                                  data_name = 'quan_batch_rf_even_fwer_40',
                                                  regularize = F)


quan_result_batch_even_rf_fwer_50 <- runModels(quan_cases,
                                              model = 'rf',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_even_fwer_50)

quan_table_batch_even_rf_fwer_50 <- extractResults(quan_result_batch_even_rf_fwer_50 , 
                                                  data_name = 'quan_batch_even_rf_fwer_50',
                                                  regularize = F)


###########################################################################################################################
# batch

##################
# Even counts

##########
#sig

quan_result_batch_even_rf_sig_10 <- runModels(quan_cases,
                                             model = 'rf',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_even_sig_10)

quan_table_batch_even_rf_sig_10 <- extractResults(quan_result_batch_even_rf_sig_10 , 
                                                 data_name = 'quan_batch_even_rf_sig_10',
                                                 regularize = F)

###########################################################################################################################
# batch

##################
# Uneven counts

##########
#normal

#### quan 10 
quan_result_batch_uneven_rf_10 <- runModels(quan_cases,
                                           model = 'rf',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_uneven_10)

quan_table_batch_uneven_rf_10 <- extractResults(quan_result_batch_uneven_rf_10, 
                                               data_name = 'quan_batch_uneven_rf_10',
                                               regularize = F)


#### quan 20 
quan_result_batch_uneven_rf_20 <- runModels(quan_cases,
                                           model = 'rf',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_uneven_20)

quan_table_batch_uneven_rf_20 <- extractResults(quan_result_batch_uneven_rf_20, 
                                               data_name = 'quan_batch_uneven_rf_20',
                                               regularize = F)


#### quan 30 
quan_result_batch_uneven_rf_30 <- runModels(quan_cases,
                                           model = 'rf',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_uneven_30)

quan_table_batch_uneven_rf_30 <- extractResults(quan_result_batch_uneven_rf_30, 
                                               data_name = 'quan_batch_uneven_rf_30',
                                               regularize = F)


#### quan 40 
quan_result_batch_uneven_rf_40 <- runModels(quan_cases,
                                           model = 'rf',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_uneven_40)

quan_table_batch_uneven_rf_40 <- extractResults(quan_result_batch_uneven_rf_40, 
                                               data_name = 'quan_batch_uneven_rf_40',
                                               regularize = F)


#### quan 10 
quan_result_batch_uneven_rf_50 <- runModels(quan_cases,
                                           model = 'rf',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_uneven_50)

quan_table_batch_uneven_rf_50 <- extractResults(quan_result_batch_uneven_rf_50, 
                                               data_name = 'quan_batch_uneven_rf_50',
                                               regularize = F)

###########################################################################################################################
# batch

##################
# uneven counts

########
# fwer

quan_result_batch_uneven_rf_fwer_10 <- runModels(quan_cases,
                                                model = 'rf',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_uneven_fwer_10)

quan_table_batch_uneven_rf_fwer_10 <- extractResults(quan_result_batch_uneven_rf_fwer_10 , 
                                                    data_name = 'quan_batch_rf_uneven_fwer_10',
                                                    regularize = F)



quan_result_batch_uneven_rf_fwer_20 <- runModels(quan_cases,
                                                model = 'rf',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_uneven_fwer_20)

quan_table_batch_uneven_rf_fwer_20 <- extractResults(quan_result_batch_uneven_rf_fwer_20 , 
                                                    data_name = 'quan_batch_rf_uneven_fwer_20',
                                                    regularize = F)


quan_result_batch_uneven_rf_fwer_30 <- runModels(quan_cases,
                                                model = 'rf',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_uneven_fwer_30)

quan_table_batch_uneven_rf_fwer_30 <- extractResults(quan_result_batch_uneven_rf_fwer_30 , 
                                                    data_name = 'quan_batch_rf_uneven_fwer_30',
                                                    regularize = F)


quan_result_batch_uneven_rf_fwer_40 <- runModels(quan_cases,
                                                model = 'rf',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_uneven_fwer_40)

quan_table_batch_uneven_rf_fwer_40 <- extractResults(quan_result_batch_uneven_rf_fwer_40 , 
                                                    data_name = 'quan_batch_rf_uneven_fwer_40',
                                                    regularize = F)


quan_result_batch_uneven_rf_fwer_50 <- runModels(quan_cases,
                                                model = 'rf',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_uneven_fwer_50)

quan_table_batch_uneven_rf_fwer_50 <- extractResults(quan_result_batch_uneven_rf_fwer_50 , 
                                                    data_name = 'quan_batch_uneven_rf_fwer_50',
                                                    regularize = F)


###########################################################################################################################
# batch

##################
# uneven counts

##########
#sig

quan_result_batch_uneven_rf_sig_10 <- runModels(quan_cases,
                                               model = 'rf',
                                               bump_hunter = T, 
                                               bump_hunter_data = quan_uneven_sig_10)

quan_table_batch_uneven_rf_sig_10 <- extractResults(quan_result_batch_uneven_rf_sig_10 , 
                                                   data_name = 'quan_batch_rf_uneven_sig_10',
                                                   regularize = F)


# get full table
quan_batch_rf <- rbind(quan_table_batch_even_rf_10,
                      quan_table_batch_even_rf_20,
                      quan_table_batch_even_rf_30,
                      quan_table_batch_even_rf_40,
                      quan_table_batch_even_rf_50,
                      quan_table_batch_even_rf_fwer_10,
                      quan_table_batch_even_rf_fwer_20,
                      quan_table_batch_even_rf_fwer_30,
                      quan_table_batch_even_rf_fwer_40,
                      quan_table_batch_even_rf_fwer_50,
                      quan_table_batch_even_rf_sig_10,
                      quan_table_batch_uneven_rf_10,
                      quan_table_batch_uneven_rf_20,
                      quan_table_batch_uneven_rf_30,
                      quan_table_batch_uneven_rf_40,
                      quan_table_batch_uneven_rf_50,
                      quan_table_batch_uneven_rf_fwer_10,
                      quan_table_batch_uneven_rf_fwer_20,
                      quan_table_batch_uneven_rf_fwer_30,
                      quan_table_batch_uneven_rf_fwer_40,
                      quan_table_batch_uneven_rf_fwer_50,
                      quan_table_batch_uneven_rf_sig_10)

#save table 
saveRDS(quan_batch_rf, 
        file = paste0(quan_folder, '/quan_table_batch_rf.rda'))

# remove
rm(quan_table_batch_even_rf_10,
   quan_table_batch_even_rf_20,
   quan_table_batch_even_rf_30,
   quan_table_batch_even_rf_40,
   quan_table_batch_even_rf_50,
   quan_table_batch_even_rf_fwer_10,
   quan_table_batch_even_rf_fwer_20,
   quan_table_batch_even_rf_fwer_30,
   quan_table_batch_even_rf_fwer_40,
   quan_table_batch_even_rf_fwer_50,
   quan_table_batch_even_rf_sig_10,
   quan_table_batch_uneven_rf_10,
   quan_table_batch_uneven_rf_20,
   quan_table_batch_uneven_rf_30,
   quan_table_batch_uneven_rf_40,
   quan_table_batch_uneven_rf_50,
   quan_table_batch_uneven_rf_fwer_10,
   quan_table_batch_uneven_rf_fwer_20,
   quan_table_batch_uneven_rf_fwer_30,
   quan_table_batch_uneven_rf_fwer_40,
   quan_table_batch_uneven_rf_fwer_50,
   quan_table_batch_uneven_rf_sig_10)

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

#### quan 10 
quan_result_unbatch_even_enet_10 <- runModels(quan_cases,
                                             model = 'enet',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_even_10)

quan_table_unbatch_even_enet_10 <- extractResults(quan_result_unbatch_even_enet_10, 
                                                 data_name = 'quan_unbatch_even_enet_10',
                                                 regularize = F)


#### quan 20 
quan_result_unbatch_even_enet_20 <- runModels(quan_cases,
                                             model = 'enet',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_even_20)

quan_table_unbatch_even_enet_20 <- extractResults(quan_result_unbatch_even_enet_20, 
                                                 data_name = 'quan_unbatch_even_enet_20',
                                                 regularize = F)


#### quan 30 
quan_result_unbatch_even_enet_30 <- runModels(quan_cases,
                                             model = 'enet',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_even_30)

quan_table_unbatch_even_enet_30 <- extractResults(quan_result_unbatch_even_enet_30, 
                                                 data_name = 'quan_unbatch_even_enet_30',
                                                 regularize = F)


#### quan 40 
quan_result_unbatch_even_enet_40 <- runModels(quan_cases,
                                             model = 'enet',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_even_40)

quan_table_unbatch_even_enet_40 <- extractResults(quan_result_unbatch_even_enet_40, 
                                                 data_name = 'quan_unbatch_even_enet_40',
                                                 regularize = F)


#### quan 10 
quan_result_unbatch_even_enet_50 <- runModels(quan_cases,
                                             model = 'enet',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_even_50)

quan_table_unbatch_even_enet_50 <- extractResults(quan_result_unbatch_even_enet_50, 
                                                 data_name = 'quan_unbatch_even_enet_50',
                                                 regularize = F)

###########################################################################################################################
# No batch

##################
# Even counts

########
# fwer

quan_result_unbatch_even_enet_fwer_10 <- runModels(quan_cases,
                                                  model = 'enet',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_even_fwer_10)

quan_table_unbatch_even_enet_fwer_10 <- extractResults(quan_result_unbatch_even_enet_fwer_10 , 
                                                      data_name = 'quan_unbatch_enet_even_fwer_10',
                                                      regularize = F)



quan_result_unbatch_even_enet_fwer_20 <- runModels(quan_cases,
                                                  model = 'enet',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_even_fwer_20)

quan_table_unbatch_even_enet_fwer_20 <- extractResults(quan_result_unbatch_even_enet_fwer_20 , 
                                                      data_name = 'quan_unbatch_enet_even_fwer_20',
                                                      regularize = F)


quan_result_unbatch_even_enet_fwer_30 <- runModels(quan_cases,
                                                  model = 'enet',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_even_fwer_30)

quan_table_unbatch_even_enet_fwer_30 <- extractResults(quan_result_unbatch_even_enet_fwer_30 , 
                                                      data_name = 'quan_unbatch_enet_even_fwer_30',
                                                      regularize = F)


quan_result_unbatch_even_enet_fwer_40 <- runModels(quan_cases,
                                                  model = 'enet',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_even_fwer_40)

quan_table_unbatch_even_enet_fwer_40 <- extractResults(quan_result_unbatch_even_enet_fwer_40 , 
                                                      data_name = 'quan_unbatch_enet_even_fwer_40',
                                                      regularize = F)


quan_result_unbatch_even_enet_fwer_50 <- runModels(quan_cases,
                                                  model = 'enet',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_even_fwer_50)

quan_table_unbatch_even_enet_fwer_50 <- extractResults(quan_result_unbatch_even_enet_fwer_50 , 
                                                      data_name = 'quan_unbatch_even_enet_fwer_50',
                                                      regularize = F)


###########################################################################################################################
# No batch

##################
# Even counts

##########
#sig

quan_result_unbatch_even_enet_sig_10 <- runModels(quan_cases,
                                                 model = 'enet',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = quan_even_sig_10)

quan_table_unbatch_even_enet_sig_10 <- extractResults(quan_result_unbatch_even_enet_sig_10 , 
                                                     data_name = 'quan_unbatch_even_enet_sig_10',
                                                     regularize = F)

###########################################################################################################################
# No batch

##################
# Uneven counts

##########
#normal

#### quan 10 
quan_result_unbatch_uneven_enet_10 <- runModels(quan_cases,
                                               model = 'enet',
                                               bump_hunter = T, 
                                               bump_hunter_data = quan_uneven_10)

quan_table_unbatch_uneven_enet_10 <- extractResults(quan_result_unbatch_uneven_enet_10, 
                                                   data_name = 'quan_unbatch_uneven_enet_10',
                                                   regularize = F)


#### quan 20 
quan_result_unbatch_uneven_enet_20 <- runModels(quan_cases,
                                               model = 'enet',
                                               bump_hunter = T, 
                                               bump_hunter_data = quan_uneven_20)

quan_table_unbatch_uneven_enet_20 <- extractResults(quan_result_unbatch_uneven_enet_20, 
                                                   data_name = 'quan_unbatch_uneven_enet_20',
                                                   regularize = F)


#### quan 30 
quan_result_unbatch_uneven_enet_30 <- runModels(quan_cases,
                                               model = 'enet',
                                               bump_hunter = T, 
                                               bump_hunter_data = quan_uneven_30)

quan_table_unbatch_uneven_enet_30 <- extractResults(quan_result_unbatch_uneven_enet_30, 
                                                   data_name = 'quan_unbatch_uneven_enet_30',
                                                   regularize = F)


#### quan 40 
quan_result_unbatch_uneven_enet_40 <- runModels(quan_cases,
                                               model = 'enet',
                                               bump_hunter = T, 
                                               bump_hunter_data = quan_uneven_40)

quan_table_unbatch_uneven_enet_40 <- extractResults(quan_result_unbatch_uneven_enet_40, 
                                                   data_name = 'quan_unbatch_uneven_enet_40',
                                                   regularize = F)


#### quan 10 
quan_result_unbatch_uneven_enet_50 <- runModels(quan_cases,
                                               model = 'enet',
                                               bump_hunter = T, 
                                               bump_hunter_data = quan_uneven_50)

quan_table_unbatch_uneven_enet_50 <- extractResults(quan_result_unbatch_uneven_enet_50, 
                                                   data_name = 'quan_unbatch_uneven_enet_50',
                                                   regularize = F)

###########################################################################################################################
# No batch

##################
# uneven counts

########
# fwer

quan_result_unbatch_uneven_enet_fwer_10 <- runModels(quan_cases,
                                                    model = 'enet',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = quan_uneven_fwer_10)

quan_table_unbatch_uneven_enet_fwer_10 <- extractResults(quan_result_unbatch_uneven_enet_fwer_10 , 
                                                        data_name = 'quan_unbatch_enet_uneven_fwer_10',
                                                        regularize = F)



quan_result_unbatch_uneven_enet_fwer_20 <- runModels(quan_cases,
                                                    model = 'enet',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = quan_uneven_fwer_20)

quan_table_unbatch_uneven_enet_fwer_20 <- extractResults(quan_result_unbatch_uneven_enet_fwer_20 , 
                                                        data_name = 'quan_unbatch_enet_uneven_fwer_20',
                                                        regularize = F)


quan_result_unbatch_uneven_enet_fwer_30 <- runModels(quan_cases,
                                                    model = 'enet',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = quan_uneven_fwer_30)

quan_table_unbatch_uneven_enet_fwer_30 <- extractResults(quan_result_unbatch_uneven_enet_fwer_30 , 
                                                        data_name = 'quan_unbatch_enet_uneven_fwer_30',
                                                        regularize = F)


quan_result_unbatch_uneven_enet_fwer_40 <- runModels(quan_cases,
                                                    model = 'enet',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = quan_uneven_fwer_40)

quan_table_unbatch_uneven_enet_fwer_40 <- extractResults(quan_result_unbatch_uneven_enet_fwer_40 , 
                                                        data_name = 'quan_unbatch_enet_uneven_fwer_40',
                                                        regularize = F)


quan_result_unbatch_uneven_enet_fwer_50 <- runModels(quan_cases,
                                                    model = 'enet',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = quan_uneven_fwer_50)

quan_table_unbatch_uneven_enet_fwer_50 <- extractResults(quan_result_unbatch_uneven_enet_fwer_50 , 
                                                        data_name = 'quan_unbatch_uneven_enet_fwer_50',
                                                        regularize = F)


###########################################################################################################################
# No batch

##################
# uneven counts

##########
#sig

quan_result_unbatch_uneven_enet_sig_10 <- runModels(quan_cases,
                                                   model = 'enet',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = quan_uneven_sig_10)

quan_table_unbatch_uneven_enet_sig_10 <- extractResults(quan_result_unbatch_uneven_enet_sig_10 , 
                                                       data_name = 'quan_unbatch_enet_uneven_sig_10',
                                                       regularize = F)


# get full table
quan_unbatch_enet <- rbind(quan_table_unbatch_even_enet_10,
                          quan_table_unbatch_even_enet_20,
                          quan_table_unbatch_even_enet_30,
                          quan_table_unbatch_even_enet_40,
                          quan_table_unbatch_even_enet_50,
                          quan_table_unbatch_even_enet_fwer_10,
                          quan_table_unbatch_even_enet_fwer_20,
                          quan_table_unbatch_even_enet_fwer_30,
                          quan_table_unbatch_even_enet_fwer_40,
                          quan_table_unbatch_even_enet_fwer_50,
                          quan_table_unbatch_even_enet_sig_10,
                          quan_table_unbatch_uneven_enet_10,
                          quan_table_unbatch_uneven_enet_20,
                          quan_table_unbatch_uneven_enet_30,
                          quan_table_unbatch_uneven_enet_40,
                          quan_table_unbatch_uneven_enet_50,
                          quan_table_unbatch_uneven_enet_fwer_10,
                          quan_table_unbatch_uneven_enet_fwer_20,
                          quan_table_unbatch_uneven_enet_fwer_30,
                          quan_table_unbatch_uneven_enet_fwer_40,
                          quan_table_unbatch_uneven_enet_fwer_50,
                          quan_table_unbatch_uneven_enet_sig_10)

#save table 
saveRDS(quan_unbatch_enet, 
        file = paste0(quan_folder, '/quan_table_unbatch_enet.rda'))

# remove
rm(quan_table_unbatch_even_enet_10,
   quan_table_unbatch_even_enet_20,
   quan_table_unbatch_even_enet_30,
   quan_table_unbatch_even_enet_40,
   quan_table_unbatch_even_enet_50,
   quan_table_unbatch_even_enet_fwer_10,
   quan_table_unbatch_even_enet_fwer_20,
   quan_table_unbatch_even_enet_fwer_30,
   quan_table_unbatch_even_enet_fwer_40,
   quan_table_unbatch_even_enet_fwer_50,
   quan_table_unbatch_even_enet_sig_10,
   quan_table_unbatch_uneven_enet_10,
   quan_table_unbatch_uneven_enet_20,
   quan_table_unbatch_uneven_enet_30,
   quan_table_unbatch_uneven_enet_40,
   quan_table_unbatch_uneven_enet_50,
   quan_table_unbatch_uneven_enet_fwer_10,
   quan_table_unbatch_uneven_enet_fwer_20,
   quan_table_unbatch_uneven_enet_fwer_30,
   quan_table_unbatch_uneven_enet_fwer_40,
   quan_table_unbatch_uneven_enet_fwer_50,
   quan_table_unbatch_uneven_enet_sig_10)



########################################################################################################################
# RANDOM FOREST WITH BH
###########################################################################################################################

#################################################################################
# batch

##################
# Even counts

##########
#normal

#### quan 10 
quan_result_batch_even_enet_10 <- runModels(quan_cases,
                                           model = 'enet',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_even_10)

quan_table_batch_even_enet_10 <- extractResults(quan_result_batch_even_enet_10, 
                                               data_name = 'quan_batch_even_enet_10',
                                               regularize = F)


#### quan 20 
quan_result_batch_even_enet_20 <- runModels(quan_cases,
                                           model = 'enet',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_even_20)

quan_table_batch_even_enet_20 <- extractResults(quan_result_batch_even_enet_20, 
                                               data_name = 'quan_batch_even_enet_20',
                                               regularize = F)


#### quan 30 
quan_result_batch_even_enet_30 <- runModels(quan_cases,
                                           model = 'enet',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_even_30)

quan_table_batch_even_enet_30 <- extractResults(quan_result_batch_even_enet_30, 
                                               data_name = 'quan_batch_even_enet_30',
                                               regularize = F)


#### quan 40 
quan_result_batch_even_enet_40 <- runModels(quan_cases,
                                           model = 'enet',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_even_40)

quan_table_batch_even_enet_40 <- extractResults(quan_result_batch_even_enet_40, 
                                               data_name = 'quan_batch_even_enet_40',
                                               regularize = F)


#### quan 10 
quan_result_batch_even_enet_50 <- runModels(quan_cases,
                                           model = 'enet',
                                           bump_hunter = T, 
                                           bump_hunter_data = quan_even_50)

quan_table_batch_even_enet_50 <- extractResults(quan_result_batch_even_enet_50, 
                                               data_name = 'quan_batch_even_enet_50',
                                               regularize = F)

###########################################################################################################################
# batch

##################
# Even counts

########
# fwer

quan_result_batch_even_enet_fwer_10 <- runModels(quan_cases,
                                                model = 'enet',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_even_fwer_10)

quan_table_batch_even_enet_fwer_10 <- extractResults(quan_result_batch_even_enet_fwer_10 , 
                                                    data_name = 'quan_batch_enet_even_fwer_10',
                                                    regularize = F)



quan_result_batch_even_enet_fwer_20 <- runModels(quan_cases,
                                                model = 'enet',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_even_fwer_20)

quan_table_batch_even_enet_fwer_20 <- extractResults(quan_result_batch_even_enet_fwer_20 , 
                                                    data_name = 'quan_batch_enet_even_fwer_20',
                                                    regularize = F)


quan_result_batch_even_enet_fwer_30 <- runModels(quan_cases,
                                                model = 'enet',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_even_fwer_30)

quan_table_batch_even_enet_fwer_30 <- extractResults(quan_result_batch_even_enet_fwer_30 , 
                                                    data_name = 'quan_batch_enet_even_fwer_30',
                                                    regularize = F)


quan_result_batch_even_enet_fwer_40 <- runModels(quan_cases,
                                                model = 'enet',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_even_fwer_40)

quan_table_batch_even_enet_fwer_40 <- extractResults(quan_result_batch_even_enet_fwer_40 , 
                                                    data_name = 'quan_batch_enet_even_fwer_40',
                                                    regularize = F)


quan_result_batch_even_enet_fwer_50 <- runModels(quan_cases,
                                                model = 'enet',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_even_fwer_50)

quan_table_batch_even_enet_fwer_50 <- extractResults(quan_result_batch_even_enet_fwer_50 , 
                                                    data_name = 'quan_batch_even_enet_fwer_50',
                                                    regularize = F)


###########################################################################################################################
# batch

##################
# Even counts

##########
#sig

quan_result_batch_even_enet_sig_10 <- runModels(quan_cases,
                                               model = 'enet',
                                               bump_hunter = T, 
                                               bump_hunter_data = quan_even_sig_10)

quan_table_batch_even_enet_sig_10 <- extractResults(quan_result_batch_even_enet_sig_10 , 
                                                   data_name = 'quan_batch_even_enet_sig_10',
                                                   regularize = F)

###########################################################################################################################
# batch

##################
# Uneven counts

##########
#normal

#### quan 10 
quan_result_batch_uneven_enet_10 <- runModels(quan_cases,
                                             model = 'enet',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_uneven_10)

quan_table_batch_uneven_enet_10 <- extractResults(quan_result_batch_uneven_enet_10, 
                                                 data_name = 'quan_batch_uneven_enet_10',
                                                 regularize = F)


#### quan 20 
quan_result_batch_uneven_enet_20 <- runModels(quan_cases,
                                             model = 'enet',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_uneven_20)

quan_table_batch_uneven_enet_20 <- extractResults(quan_result_batch_uneven_enet_20, 
                                                 data_name = 'quan_batch_uneven_enet_20',
                                                 regularize = F)


#### quan 30 
quan_result_batch_uneven_enet_30 <- runModels(quan_cases,
                                             model = 'enet',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_uneven_30)

quan_table_batch_uneven_enet_30 <- extractResults(quan_result_batch_uneven_enet_30, 
                                                 data_name = 'quan_batch_uneven_enet_30',
                                                 regularize = F)


#### quan 40 
quan_result_batch_uneven_enet_40 <- runModels(quan_cases,
                                             model = 'enet',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_uneven_40)

quan_table_batch_uneven_enet_40 <- extractResults(quan_result_batch_uneven_enet_40, 
                                                 data_name = 'quan_batch_uneven_enet_40',
                                                 regularize = F)


#### quan 10 
quan_result_batch_uneven_enet_50 <- runModels(quan_cases,
                                             model = 'enet',
                                             bump_hunter = T, 
                                             bump_hunter_data = quan_uneven_50)

quan_table_batch_uneven_enet_50 <- extractResults(quan_result_batch_uneven_enet_50, 
                                                 data_name = 'quan_batch_uneven_enet_50',
                                                 regularize = F)

###########################################################################################################################
# batch

##################
# uneven counts

########
# fwer

quan_result_batch_uneven_enet_fwer_10 <- runModels(quan_cases,
                                                  model = 'enet',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_uneven_fwer_10)

quan_table_batch_uneven_enet_fwer_10 <- extractResults(quan_result_batch_uneven_enet_fwer_10 , 
                                                      data_name = 'quan_batch_enet_uneven_fwer_10',
                                                      regularize = F)



quan_result_batch_uneven_enet_fwer_20 <- runModels(quan_cases,
                                                  model = 'enet',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_uneven_fwer_20)

quan_table_batch_uneven_enet_fwer_20 <- extractResults(quan_result_batch_uneven_enet_fwer_20 , 
                                                      data_name = 'quan_batch_enet_uneven_fwer_20',
                                                      regularize = F)


quan_result_batch_uneven_enet_fwer_30 <- runModels(quan_cases,
                                                  model = 'enet',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_uneven_fwer_30)

quan_table_batch_uneven_enet_fwer_30 <- extractResults(quan_result_batch_uneven_enet_fwer_30 , 
                                                      data_name = 'quan_batch_enet_uneven_fwer_30',
                                                      regularize = F)


quan_result_batch_uneven_enet_fwer_40 <- runModels(quan_cases,
                                                  model = 'enet',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_uneven_fwer_40)

quan_table_batch_uneven_enet_fwer_40 <- extractResults(quan_result_batch_uneven_enet_fwer_40 , 
                                                      data_name = 'quan_batch_enet_uneven_fwer_40',
                                                      regularize = F)


quan_result_batch_uneven_enet_fwer_50 <- runModels(quan_cases,
                                                  model = 'enet',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_uneven_fwer_50)

quan_table_batch_uneven_enet_fwer_50 <- extractResults(quan_result_batch_uneven_enet_fwer_50 , 
                                                      data_name = 'quan_batch_uneven_enet_fwer_50',
                                                      regularize = F)


###########################################################################################################################
# batch

##################
# uneven counts

##########
#sig

quan_result_batch_uneven_enet_sig_10 <- runModels(quan_cases,
                                                 model = 'enet',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = quan_uneven_sig_10)

quan_table_batch_uneven_enet_sig_10 <- extractResults(quan_result_batch_uneven_enet_sig_10 , 
                                                     data_name = 'quan_batch_enet_uneven_sig_10',
                                                     regularize = F)


# get full table
quan_batch_enet <- rbind(quan_table_batch_even_enet_10,
                        quan_table_batch_even_enet_20,
                        quan_table_batch_even_enet_30,
                        quan_table_batch_even_enet_40,
                        quan_table_batch_even_enet_50,
                        quan_table_batch_even_enet_fwer_10,
                        quan_table_batch_even_enet_fwer_20,
                        quan_table_batch_even_enet_fwer_30,
                        quan_table_batch_even_enet_fwer_40,
                        quan_table_batch_even_enet_fwer_50,
                        quan_table_batch_even_enet_sig_10,
                        quan_table_batch_uneven_enet_10,
                        quan_table_batch_uneven_enet_20,
                        quan_table_batch_uneven_enet_30,
                        quan_table_batch_uneven_enet_40,
                        quan_table_batch_uneven_enet_50,
                        quan_table_batch_uneven_enet_fwer_10,
                        quan_table_batch_uneven_enet_fwer_20,
                        quan_table_batch_uneven_enet_fwer_30,
                        quan_table_batch_uneven_enet_fwer_40,
                        quan_table_batch_uneven_enet_fwer_50,
                        quan_table_batch_uneven_enet_sig_10)

#save table 
saveRDS(quan_batch_enet, 
        file = paste0(quan_folder, '/quan_table_batch_enet.rda'))

# remove
rm(quan_table_batch_even_enet_10,
   quan_table_batch_even_enet_20,
   quan_table_batch_even_enet_30,
   quan_table_batch_even_enet_40,
   quan_table_batch_even_enet_50,
   quan_table_batch_even_enet_fwer_10,
   quan_table_batch_even_enet_fwer_20,
   quan_table_batch_even_enet_fwer_30,
   quan_table_batch_even_enet_fwer_40,
   quan_table_batch_even_enet_fwer_50,
   quan_table_batch_even_enet_sig_10,
   quan_table_batch_uneven_enet_10,
   quan_table_batch_uneven_enet_20,
   quan_table_batch_uneven_enet_30,
   quan_table_batch_uneven_enet_40,
   quan_table_batch_uneven_enet_50,
   quan_table_batch_uneven_enet_fwer_10,
   quan_table_batch_uneven_enet_fwer_20,
   quan_table_batch_uneven_enet_fwer_30,
   quan_table_batch_uneven_enet_fwer_40,
   quan_table_batch_uneven_enet_fwer_50,
   quan_table_batch_uneven_enet_sig_10)


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

#### quan 10 
quan_result_unbatch_even_lasso_10 <- runModels(quan_cases,
                                              model = 'lasso',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_even_10)

quan_table_unbatch_even_lasso_10 <- extractResults(quan_result_unbatch_even_lasso_10, 
                                                  data_name = 'quan_unbatch_even_lasso_10',
                                                  regularize = F)


#### quan 20 
quan_result_unbatch_even_lasso_20 <- runModels(quan_cases,
                                              model = 'lasso',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_even_20)

quan_table_unbatch_even_lasso_20 <- extractResults(quan_result_unbatch_even_lasso_20, 
                                                  data_name = 'quan_unbatch_even_lasso_20',
                                                  regularize = F)


#### quan 30 
quan_result_unbatch_even_lasso_30 <- runModels(quan_cases,
                                              model = 'lasso',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_even_30)

quan_table_unbatch_even_lasso_30 <- extractResults(quan_result_unbatch_even_lasso_30, 
                                                  data_name = 'quan_unbatch_even_lasso_30',
                                                  regularize = F)


#### quan 40 
quan_result_unbatch_even_lasso_40 <- runModels(quan_cases,
                                              model = 'lasso',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_even_40)

quan_table_unbatch_even_lasso_40 <- extractResults(quan_result_unbatch_even_lasso_40, 
                                                  data_name = 'quan_unbatch_even_lasso_40',
                                                  regularize = F)


#### quan 10 
quan_result_unbatch_even_lasso_50 <- runModels(quan_cases,
                                              model = 'lasso',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_even_50)

quan_table_unbatch_even_lasso_50 <- extractResults(quan_result_unbatch_even_lasso_50, 
                                                  data_name = 'quan_unbatch_even_lasso_50',
                                                  regularize = F)

###########################################################################################################################
# No batch

##################
# Even counts

########
# fwer

quan_result_unbatch_even_lasso_fwer_10 <- runModels(quan_cases,
                                                   model = 'lasso',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = quan_even_fwer_10)

quan_table_unbatch_even_lasso_fwer_10 <- extractResults(quan_result_unbatch_even_lasso_fwer_10 , 
                                                       data_name = 'quan_unbatch_lasso_even_fwer_10',
                                                       regularize = F)



quan_result_unbatch_even_lasso_fwer_20 <- runModels(quan_cases,
                                                   model = 'lasso',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = quan_even_fwer_20)

quan_table_unbatch_even_lasso_fwer_20 <- extractResults(quan_result_unbatch_even_lasso_fwer_20 , 
                                                       data_name = 'quan_unbatch_lasso_even_fwer_20',
                                                       regularize = F)


quan_result_unbatch_even_lasso_fwer_30 <- runModels(quan_cases,
                                                   model = 'lasso',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = quan_even_fwer_30)

quan_table_unbatch_even_lasso_fwer_30 <- extractResults(quan_result_unbatch_even_lasso_fwer_30 , 
                                                       data_name = 'quan_unbatch_lasso_even_fwer_30',
                                                       regularize = F)


quan_result_unbatch_even_lasso_fwer_40 <- runModels(quan_cases,
                                                   model = 'lasso',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = quan_even_fwer_40)

quan_table_unbatch_even_lasso_fwer_40 <- extractResults(quan_result_unbatch_even_lasso_fwer_40 , 
                                                       data_name = 'quan_unbatch_lasso_even_fwer_40',
                                                       regularize = F)


quan_result_unbatch_even_lasso_fwer_50 <- runModels(quan_cases,
                                                   model = 'lasso',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = quan_even_fwer_50)

quan_table_unbatch_even_lasso_fwer_50 <- extractResults(quan_result_unbatch_even_lasso_fwer_50 , 
                                                       data_name = 'quan_unbatch_even_lasso_fwer_50',
                                                       regularize = F)


###########################################################################################################################
# No batch

##################
# Even counts

##########
#sig

quan_result_unbatch_even_lasso_sig_10 <- runModels(quan_cases,
                                                  model = 'lasso',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_even_sig_10)

quan_table_unbatch_even_lasso_sig_10 <- extractResults(quan_result_unbatch_even_lasso_sig_10 , 
                                                      data_name = 'quan_unbatch_even_lasso_sig_10',
                                                      regularize = F)

###########################################################################################################################
# No batch

##################
# Uneven counts

##########
#normal

#### quan 10 
quan_result_unbatch_uneven_lasso_10 <- runModels(quan_cases,
                                                model = 'lasso',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_uneven_10)

quan_table_unbatch_uneven_lasso_10 <- extractResults(quan_result_unbatch_uneven_lasso_10, 
                                                    data_name = 'quan_unbatch_uneven_lasso_10',
                                                    regularize = F)


#### quan 20 
quan_result_unbatch_uneven_lasso_20 <- runModels(quan_cases,
                                                model = 'lasso',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_uneven_20)

quan_table_unbatch_uneven_lasso_20 <- extractResults(quan_result_unbatch_uneven_lasso_20, 
                                                    data_name = 'quan_unbatch_uneven_lasso_20',
                                                    regularize = F)


#### quan 30 
quan_result_unbatch_uneven_lasso_30 <- runModels(quan_cases,
                                                model = 'lasso',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_uneven_30)

quan_table_unbatch_uneven_lasso_30 <- extractResults(quan_result_unbatch_uneven_lasso_30, 
                                                    data_name = 'quan_unbatch_uneven_lasso_30',
                                                    regularize = F)


#### quan 40 
quan_result_unbatch_uneven_lasso_40 <- runModels(quan_cases,
                                                model = 'lasso',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_uneven_40)

quan_table_unbatch_uneven_lasso_40 <- extractResults(quan_result_unbatch_uneven_lasso_40, 
                                                    data_name = 'quan_unbatch_uneven_lasso_40',
                                                    regularize = F)


#### quan 10 
quan_result_unbatch_uneven_lasso_50 <- runModels(quan_cases,
                                                model = 'lasso',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_uneven_50)

quan_table_unbatch_uneven_lasso_50 <- extractResults(quan_result_unbatch_uneven_lasso_50, 
                                                    data_name = 'quan_unbatch_uneven_lasso_50',
                                                    regularize = F)

###########################################################################################################################
# No batch

##################
# uneven counts

########
# fwer

quan_result_unbatch_uneven_lasso_fwer_10 <- runModels(quan_cases,
                                                     model = 'lasso',
                                                     bump_hunter = T, 
                                                     bump_hunter_data = quan_uneven_fwer_10)

quan_table_unbatch_uneven_lasso_fwer_10 <- extractResults(quan_result_unbatch_uneven_lasso_fwer_10 , 
                                                         data_name = 'quan_unbatch_lasso_uneven_fwer_10',
                                                         regularize = F)



quan_result_unbatch_uneven_lasso_fwer_20 <- runModels(quan_cases,
                                                     model = 'lasso',
                                                     bump_hunter = T, 
                                                     bump_hunter_data = quan_uneven_fwer_20)

quan_table_unbatch_uneven_lasso_fwer_20 <- extractResults(quan_result_unbatch_uneven_lasso_fwer_20 , 
                                                         data_name = 'quan_unbatch_lasso_uneven_fwer_20',
                                                         regularize = F)


quan_result_unbatch_uneven_lasso_fwer_30 <- runModels(quan_cases,
                                                     model = 'lasso',
                                                     bump_hunter = T, 
                                                     bump_hunter_data = quan_uneven_fwer_30)

quan_table_unbatch_uneven_lasso_fwer_30 <- extractResults(quan_result_unbatch_uneven_lasso_fwer_30 , 
                                                         data_name = 'quan_unbatch_lasso_uneven_fwer_30',
                                                         regularize = F)


quan_result_unbatch_uneven_lasso_fwer_40 <- runModels(quan_cases,
                                                     model = 'lasso',
                                                     bump_hunter = T, 
                                                     bump_hunter_data = quan_uneven_fwer_40)

quan_table_unbatch_uneven_lasso_fwer_40 <- extractResults(quan_result_unbatch_uneven_lasso_fwer_40 , 
                                                         data_name = 'quan_unbatch_lasso_uneven_fwer_40',
                                                         regularize = F)


quan_result_unbatch_uneven_lasso_fwer_50 <- runModels(quan_cases,
                                                     model = 'lasso',
                                                     bump_hunter = T, 
                                                     bump_hunter_data = quan_uneven_fwer_50)

quan_table_unbatch_uneven_lasso_fwer_50 <- extractResults(quan_result_unbatch_uneven_lasso_fwer_50 , 
                                                         data_name = 'quan_unbatch_uneven_lasso_fwer_50',
                                                         regularize = F)


###########################################################################################################################
# No batch

##################
# uneven counts

##########
#sig

quan_result_unbatch_uneven_lasso_sig_10 <- runModels(quan_cases,
                                                    model = 'lasso',
                                                    bump_hunter = T, 
                                                    bump_hunter_data = quan_uneven_sig_10)

quan_table_unbatch_uneven_lasso_sig_10 <- extractResults(quan_result_unbatch_uneven_lasso_sig_10 , 
                                                        data_name = 'quan_unbatch_lasso_uneven_sig_10',
                                                        regularize = F)


# get full table
quan_unbatch_lasso <- rbind(quan_table_unbatch_even_lasso_10,
                           quan_table_unbatch_even_lasso_20,
                           quan_table_unbatch_even_lasso_30,
                           quan_table_unbatch_even_lasso_40,
                           quan_table_unbatch_even_lasso_50,
                           quan_table_unbatch_even_lasso_fwer_10,
                           quan_table_unbatch_even_lasso_fwer_20,
                           quan_table_unbatch_even_lasso_fwer_30,
                           quan_table_unbatch_even_lasso_fwer_40,
                           quan_table_unbatch_even_lasso_fwer_50,
                           quan_table_unbatch_even_lasso_sig_10,
                           quan_table_unbatch_uneven_lasso_10,
                           quan_table_unbatch_uneven_lasso_20,
                           quan_table_unbatch_uneven_lasso_30,
                           quan_table_unbatch_uneven_lasso_40,
                           quan_table_unbatch_uneven_lasso_50,
                           quan_table_unbatch_uneven_lasso_fwer_10,
                           quan_table_unbatch_uneven_lasso_fwer_20,
                           quan_table_unbatch_uneven_lasso_fwer_30,
                           quan_table_unbatch_uneven_lasso_fwer_40,
                           quan_table_unbatch_uneven_lasso_fwer_50,
                           quan_table_unbatch_uneven_lasso_sig_10)

#save table 
saveRDS(quan_unbatch_lasso, 
        file = paste0(quan_folder, '/quan_table_unbatch_lasso.rda'))

# remove
rm(quan_table_unbatch_even_lasso_10,
   quan_table_unbatch_even_lasso_20,
   quan_table_unbatch_even_lasso_30,
   quan_table_unbatch_even_lasso_40,
   quan_table_unbatch_even_lasso_50,
   quan_table_unbatch_even_lasso_fwer_10,
   quan_table_unbatch_even_lasso_fwer_20,
   quan_table_unbatch_even_lasso_fwer_30,
   quan_table_unbatch_even_lasso_fwer_40,
   quan_table_unbatch_even_lasso_fwer_50,
   quan_table_unbatch_even_lasso_sig_10,
   quan_table_unbatch_uneven_lasso_10,
   quan_table_unbatch_uneven_lasso_20,
   quan_table_unbatch_uneven_lasso_30,
   quan_table_unbatch_uneven_lasso_40,
   quan_table_unbatch_uneven_lasso_50,
   quan_table_unbatch_uneven_lasso_fwer_10,
   quan_table_unbatch_uneven_lasso_fwer_20,
   quan_table_unbatch_uneven_lasso_fwer_30,
   quan_table_unbatch_uneven_lasso_fwer_40,
   quan_table_unbatch_uneven_lasso_fwer_50,
   quan_table_unbatch_uneven_lasso_sig_10)



########################################################################################################################
# RANDOM FOREST WITH BH
###########################################################################################################################

#################################################################################
# batch

##################
# Even counts

##########
#normal

#### quan 10 
quan_result_batch_even_lasso_10 <- runModels(quan_cases,
                                            model = 'lasso',
                                            bump_hunter = T, 
                                            bump_hunter_data = quan_even_10)

quan_table_batch_even_lasso_10 <- extractResults(quan_result_batch_even_lasso_10, 
                                                data_name = 'quan_batch_even_lasso_10',
                                                regularize = F)


#### quan 20 
quan_result_batch_even_lasso_20 <- runModels(quan_cases,
                                            model = 'lasso',
                                            bump_hunter = T, 
                                            bump_hunter_data = quan_even_20)

quan_table_batch_even_lasso_20 <- extractResults(quan_result_batch_even_lasso_20, 
                                                data_name = 'quan_batch_even_lasso_20',
                                                regularize = F)


#### quan 30 
quan_result_batch_even_lasso_30 <- runModels(quan_cases,
                                            model = 'lasso',
                                            bump_hunter = T, 
                                            bump_hunter_data = quan_even_30)

quan_table_batch_even_lasso_30 <- extractResults(quan_result_batch_even_lasso_30, 
                                                data_name = 'quan_batch_even_lasso_30',
                                                regularize = F)


#### quan 40 
quan_result_batch_even_lasso_40 <- runModels(quan_cases,
                                            model = 'lasso',
                                            bump_hunter = T, 
                                            bump_hunter_data = quan_even_40)

quan_table_batch_even_lasso_40 <- extractResults(quan_result_batch_even_lasso_40, 
                                                data_name = 'quan_batch_even_lasso_40',
                                                regularize = F)


#### quan 10 
quan_result_batch_even_lasso_50 <- runModels(quan_cases,
                                            model = 'lasso',
                                            bump_hunter = T, 
                                            bump_hunter_data = quan_even_50)

quan_table_batch_even_lasso_50 <- extractResults(quan_result_batch_even_lasso_50, 
                                                data_name = 'quan_batch_even_lasso_50',
                                                regularize = F)

###########################################################################################################################
# batch

##################
# Even counts

########
# fwer

quan_result_batch_even_lasso_fwer_10 <- runModels(quan_cases,
                                                 model = 'lasso',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = quan_even_fwer_10)

quan_table_batch_even_lasso_fwer_10 <- extractResults(quan_result_batch_even_lasso_fwer_10 , 
                                                     data_name = 'quan_batch_lasso_even_fwer_10',
                                                     regularize = F)



quan_result_batch_even_lasso_fwer_20 <- runModels(quan_cases,
                                                 model = 'lasso',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = quan_even_fwer_20)

quan_table_batch_even_lasso_fwer_20 <- extractResults(quan_result_batch_even_lasso_fwer_20 , 
                                                     data_name = 'quan_batch_lasso_even_fwer_20',
                                                     regularize = F)


quan_result_batch_even_lasso_fwer_30 <- runModels(quan_cases,
                                                 model = 'lasso',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = quan_even_fwer_30)

quan_table_batch_even_lasso_fwer_30 <- extractResults(quan_result_batch_even_lasso_fwer_30 , 
                                                     data_name = 'quan_batch_lasso_even_fwer_30',
                                                     regularize = F)


quan_result_batch_even_lasso_fwer_40 <- runModels(quan_cases,
                                                 model = 'lasso',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = quan_even_fwer_40)

quan_table_batch_even_lasso_fwer_40 <- extractResults(quan_result_batch_even_lasso_fwer_40 , 
                                                     data_name = 'quan_batch_lasso_even_fwer_40',
                                                     regularize = F)


quan_result_batch_even_lasso_fwer_50 <- runModels(quan_cases,
                                                 model = 'lasso',
                                                 bump_hunter = T, 
                                                 bump_hunter_data = quan_even_fwer_50)

quan_table_batch_even_lasso_fwer_50 <- extractResults(quan_result_batch_even_lasso_fwer_50 , 
                                                     data_name = 'quan_batch_even_lasso_fwer_50',
                                                     regularize = F)


###########################################################################################################################
# batch

##################
# Even counts

##########
#sig

quan_result_batch_even_lasso_sig_10 <- runModels(quan_cases,
                                                model = 'lasso',
                                                bump_hunter = T, 
                                                bump_hunter_data = quan_even_sig_10)

quan_table_batch_even_lasso_sig_10 <- extractResults(quan_result_batch_even_lasso_sig_10 , 
                                                    data_name = 'quan_batch_even_lasso_sig_10',
                                                    regularize = F)

###########################################################################################################################
# batch

##################
# Uneven counts

##########
#normal

#### quan 10 
quan_result_batch_uneven_lasso_10 <- runModels(quan_cases,
                                              model = 'lasso',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_uneven_10)

quan_table_batch_uneven_lasso_10 <- extractResults(quan_result_batch_uneven_lasso_10, 
                                                  data_name = 'quan_batch_uneven_lasso_10',
                                                  regularize = F)


#### quan 20 
quan_result_batch_uneven_lasso_20 <- runModels(quan_cases,
                                              model = 'lasso',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_uneven_20)

quan_table_batch_uneven_lasso_20 <- extractResults(quan_result_batch_uneven_lasso_20, 
                                                  data_name = 'quan_batch_uneven_lasso_20',
                                                  regularize = F)


#### quan 30 
quan_result_batch_uneven_lasso_30 <- runModels(quan_cases,
                                              model = 'lasso',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_uneven_30)

quan_table_batch_uneven_lasso_30 <- extractResults(quan_result_batch_uneven_lasso_30, 
                                                  data_name = 'quan_batch_uneven_lasso_30',
                                                  regularize = F)


#### quan 40 
quan_result_batch_uneven_lasso_40 <- runModels(quan_cases,
                                              model = 'lasso',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_uneven_40)

quan_table_batch_uneven_lasso_40 <- extractResults(quan_result_batch_uneven_lasso_40, 
                                                  data_name = 'quan_batch_uneven_lasso_40',
                                                  regularize = F)


#### quan 10 
quan_result_batch_uneven_lasso_50 <- runModels(quan_cases,
                                              model = 'lasso',
                                              bump_hunter = T, 
                                              bump_hunter_data = quan_uneven_50)

quan_table_batch_uneven_lasso_50 <- extractResults(quan_result_batch_uneven_lasso_50, 
                                                  data_name = 'quan_batch_uneven_lasso_50',
                                                  regularize = F)

###########################################################################################################################
# batch

##################
# uneven counts

########
# fwer

quan_result_batch_uneven_lasso_fwer_10 <- runModels(quan_cases,
                                                   model = 'lasso',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = quan_uneven_fwer_10)

quan_table_batch_uneven_lasso_fwer_10 <- extractResults(quan_result_batch_uneven_lasso_fwer_10 , 
                                                       data_name = 'quan_batch_lasso_uneven_fwer_10',
                                                       regularize = F)



quan_result_batch_uneven_lasso_fwer_20 <- runModels(quan_cases,
                                                   model = 'lasso',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = quan_uneven_fwer_20)

quan_table_batch_uneven_lasso_fwer_20 <- extractResults(quan_result_batch_uneven_lasso_fwer_20 , 
                                                       data_name = 'quan_batch_lasso_uneven_fwer_20',
                                                       regularize = F)


quan_result_batch_uneven_lasso_fwer_30 <- runModels(quan_cases,
                                                   model = 'lasso',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = quan_uneven_fwer_30)

quan_table_batch_uneven_lasso_fwer_30 <- extractResults(quan_result_batch_uneven_lasso_fwer_30 , 
                                                       data_name = 'quan_batch_lasso_uneven_fwer_30',
                                                       regularize = F)


quan_result_batch_uneven_lasso_fwer_40 <- runModels(quan_cases,
                                                   model = 'lasso',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = quan_uneven_fwer_40)

quan_table_batch_uneven_lasso_fwer_40 <- extractResults(quan_result_batch_uneven_lasso_fwer_40 , 
                                                       data_name = 'quan_batch_lasso_uneven_fwer_40',
                                                       regularize = F)


quan_result_batch_uneven_lasso_fwer_50 <- runModels(quan_cases,
                                                   model = 'lasso',
                                                   bump_hunter = T, 
                                                   bump_hunter_data = quan_uneven_fwer_50)

quan_table_batch_uneven_lasso_fwer_50 <- extractResults(quan_result_batch_uneven_lasso_fwer_50 , 
                                                       data_name = 'quan_batch_uneven_lasso_fwer_50',
                                                       regularize = F)


###########################################################################################################################
# batch

##################
# uneven counts

##########
#sig

quan_result_batch_uneven_lasso_sig_10 <- runModels(quan_cases,
                                                  model = 'lasso',
                                                  bump_hunter = T, 
                                                  bump_hunter_data = quan_uneven_sig_10)

quan_table_batch_uneven_lasso_sig_10 <- extractResults(quan_result_batch_uneven_lasso_sig_10 , 
                                                      data_name = 'quan_batch_lasso_uneven_sig_10',
                                                      regularize = F)


# get full table
quan_batch_lasso <- rbind(quan_table_batch_even_lasso_10,
                         quan_table_batch_even_lasso_20,
                         quan_table_batch_even_lasso_30,
                         quan_table_batch_even_lasso_40,
                         quan_table_batch_even_lasso_50,
                         quan_table_batch_even_lasso_fwer_10,
                         quan_table_batch_even_lasso_fwer_20,
                         quan_table_batch_even_lasso_fwer_30,
                         quan_table_batch_even_lasso_fwer_40,
                         quan_table_batch_even_lasso_fwer_50,
                         quan_table_batch_even_lasso_sig_10,
                         quan_table_batch_uneven_lasso_10,
                         quan_table_batch_uneven_lasso_20,
                         quan_table_batch_uneven_lasso_30,
                         quan_table_batch_uneven_lasso_40,
                         quan_table_batch_uneven_lasso_50,
                         quan_table_batch_uneven_lasso_fwer_10,
                         quan_table_batch_uneven_lasso_fwer_20,
                         quan_table_batch_uneven_lasso_fwer_30,
                         quan_table_batch_uneven_lasso_fwer_40,
                         quan_table_batch_uneven_lasso_fwer_50,
                         quan_table_batch_uneven_lasso_sig_10)

#save table 
saveRDS(quan_batch_lasso, 
        file = paste0(quan_folder, '/quan_table_batch_lasso.rda'))

# remove
rm(quan_table_batch_even_lasso_10,
   quan_table_batch_even_lasso_20,
   quan_table_batch_even_lasso_30,
   quan_table_batch_even_lasso_40,
   quan_table_batch_even_lasso_50,
   quan_table_batch_even_lasso_fwer_10,
   quan_table_batch_even_lasso_fwer_20,
   quan_table_batch_even_lasso_fwer_30,
   quan_table_batch_even_lasso_fwer_40,
   quan_table_batch_even_lasso_fwer_50,
   quan_table_batch_even_lasso_sig_10,
   quan_table_batch_uneven_lasso_10,
   quan_table_batch_uneven_lasso_20,
   quan_table_batch_uneven_lasso_30,
   quan_table_batch_uneven_lasso_40,
   quan_table_batch_uneven_lasso_50,
   quan_table_batch_uneven_lasso_fwer_10,
   quan_table_batch_uneven_lasso_fwer_20,
   quan_table_batch_uneven_lasso_fwer_30,
   quan_table_batch_uneven_lasso_fwer_40,
   quan_table_batch_uneven_lasso_fwer_50,
   quan_table_batch_uneven_lasso_sig_10)


