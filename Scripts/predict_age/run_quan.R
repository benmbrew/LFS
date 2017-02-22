####### This script will predict age of onset with quan preprocessed data
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
# remove quan, swan, and funnorm objects
##########
rm(list = ls(pattern = "quan_*"))
rm(list = ls(pattern = "swan_*"))
rm(list = ls(pattern = "funnorm_*"))

###################################################################################################################################
## beta_quan

##########
# threshold 0.10, 0.15, 0.20, 0.30
##########

#NORMAL
# quan 10 
quan_result_10 <- runModels(quan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = quan_bh_10)

quan_table_10 <- extractResults(quan_result_10, 
                               data_name = 'quan_10')


# quan 15
quan_result_15 <- runModels(quan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = quan_bh_15)

quan_table_15 <- extractResults(quan_result_15, 
                               data_name = 'quan_15')

# quan 20
quan_result_20 <- runModels(quan_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = quan_bh_20)

quan_table_20 <- extractResults(quan_result_20, 
                                data_name = 'quan_20')


#UNBAL
# quan 10 
quan_result_unbal_10 <- runModels(quan_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = quan_unbal_bh_10)

quan_table_unbal_10 <- extractResults(quan_result_unbal_10, 
                                data_name = 'quan_unbal_10')


# quan 15
quan_result_unbal_15 <- runModels(quan_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = quan_unbal_bh_15)

quan_table_unbal_15 <- extractResults(quan_result_unbal_15, 
                                data_name = 'quan_unbal_15')

# quan 20
quan_result_unbal_20 <- runModels(quan_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = quan_unbal_bh_20)

quan_table_unbal_20 <- extractResults(quan_result_unbal_20, 
                                data_name = 'quan_unbal_20')

# quan intersection
quan_result_int <- runModels(quan_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = quan_int_feat)

quan_table_int <- extractResults(quan_result_int, 
                                data_name = 'quan_int')


##########
# threshold 0.07 - 0.15 sig
##########

#NORMAL
# quan 10 
quan_sig_result_10 <- runModels(quan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = quan_bh_sig_10)

quan_sig_table_10 <- extractResults(quan_sig_result_10, 
                                   data_name = 'quan_sig_10')


# quan 15
quan_sig_result_15 <- runModels(quan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = quan_bh_sig_15)

quan_sig_table_15 <- extractResults(quan_sig_result_15, 
                                   data_name = 'quan_sig_15')

# quan 20
quan_sig_result_20 <- runModels(quan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = quan_bh_sig_20)

quan_sig_table_20 <- extractResults(quan_sig_result_20, 
                                   data_name = 'quan_sig_20')



# quan intersection
quan_sig_result_int <- runModels(quan_cases, 
                                bump_hunter = T, 
                                bump_hunter_data = quan_int_sig_feat)

quan_sig_table_int <- extractResults(quan_sig_result_int, 
                                    data_name = 'quan_sig_int')



#UNBAL
# quan 10 
quan_sig_result_unbal_10 <- runModels(quan_cases, 
                                     bump_hunter = T, 
                                     bump_hunter_data = quan_unbal_bh_sig_10)

quan_sig_table_unbal_10 <- extractResults(quan_sig_result_unbal_10, 
                                         data_name = 'quan_unbal_sig_10')


# quan 15
quan_sig_result_unbal_15 <- runModels(quan_cases, 
                                     bump_hunter = T, 
                                     bump_hunter_data = quan_unbal_bh_sig_15)

quan_sig_table_unbal_15 <- extractResults(quan_sig_result_unbal_15, 
                                         data_name = 'quan_unbal_sig_15')

# # quan 20
# quan_sig_result_unbal_20 <- runModels(quan_cases, 
#                                      bump_hunter = T, 
#                                      bump_hunter_data = quan_unbal_bh_sig_20)
# 
# quan_sig_table_unbal_20 <- extractResults(quan_sig_result_unbal_20, 
#                                          data_name = 'quan_unbal_sig_20')
# 

# # # quan intersection
# quan_sig_result_unbal_int <- runModels(quan_cases,
#                                       bump_hunter = T,
#                                       bump_hunter_data = quan_unbal_int_sig_feat)
# #
# quan_sig_table_unbal_int <- extractResults(quan_sig_result_unbal_int,
#                                           data_name = 'quan_unbal_sig_int')
# 



###########
# rbind tables and save RDA file
###########
quan_table <- rbind(quan_table_10, quan_table_15,
                   quan_table_20,
                   quan_sig_table_10, quan_sig_table_15,
                   quan_sig_table_20,
                   quan_table_int, quan_sig_table_int,
                   quan_table_unbal_10, quan_table_unbal_15,
                   quan_table_unbal_20, 
                   quan_sig_table_unbal_10, quan_sig_table_unbal_15
                   #quan_sig_table_unbal_20, 
                   #quan_table_unbal_int, 
                   #quan_sig_table_unbal_int)
)


# # remove data 
rm(quan_table_10, quan_table_15,
   quan_table_20,
   quan_sig_table_10, quan_sig_table_15,
   quan_sig_table_20,
   quan_table_int, quan_sig_table_int,
   quan_table_unbal_10, quan_table_unbal_15,
   quan_table_unbal_20, 
   quan_sig_table_unbal_10, quan_sig_table_unbal_15)
   #quan_sig_table_unbal_20, 
   #quan_table_unbal_int, 
   #quan_sig_table_unbal_int)

#save table 
saveRDS(quan_table, 
        file = paste0(quan_folder, '/quan_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(quan_result_10, 
        file = paste0(quan_folder, '/quan_result_10.rda'))

saveRDS(quan_result_15, 
        file = paste0(quan_folder, '/quan_result_15.rda'))

saveRDS(quan_result_20, 
        file = paste0(quan_folder, '/quan_result_20.rda'))

saveRDS(quan_result_unbal_10, 
        file = paste0(quan_folder, '/quan_result_unbal_10.rda'))

saveRDS(quan_result_unbal_15, 
        file = paste0(quan_folder, '/quan_result_unbal_15.rda'))

saveRDS(quan_result_unbal_20, 
        file = paste0(quan_folder, '/quan_result_unbal_20.rda'))
