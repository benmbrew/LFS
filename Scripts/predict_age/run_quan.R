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
# remove raw, swan, and funnorm objects
##########
rm(list = ls(pattern = "raw_*"))
rm(list = ls(pattern = "swan_*"))
rm(list = ls(pattern = "funnorm_*"))

###################################################################################################################################
## beta_quan

##########
# threshold 0.07 - 0.15 normal
##########

# quan 07 
quan_result_07 <- runModels(quan_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = quan_bh_07)

quan_table_07 <- extractResults(quan_result_07, 
                               data_name = 'quan_07')

# quan 08
quan_result_08 <- runModels(quan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = quan_bh_08)

quan_table_08 <- extractResults(quan_result_08, 
                               data_name = 'quan_08')

# quan 09 
quan_result_09 <- runModels(quan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = quan_bh_09)

quan_table_09 <- extractResults(quan_result_09, 
                               data_name = 'quan_09')

# quan 10 
quan_result_10 <- runModels(quan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = quan_bh_10)

quan_table_10 <- extractResults(quan_result_10, 
                               data_name = 'quan_10')

# quan 11
quan_result_11 <- runModels(quan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = quan_bh_11)

quan_table_11 <- extractResults(quan_result_11, 
                               data_name = 'quan_11')

# quan 12
quan_result_12 <- runModels(quan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = quan_bh_12)

quan_table_12 <- extractResults(quan_result_12, 
                               data_name = 'quan_12')

# quan 13
quan_result_13 <- runModels(quan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = quan_bh_13)

quan_table_13 <- extractResults(quan_result_13, 
                               data_name = 'quan_13')

# quan 14
quan_result_14 <- runModels(quan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = quan_bh_14)

quan_table_14 <- extractResults(quan_result_14, 
                               data_name = 'quan_14')

# quan 15
quan_result_15 <- runModels(quan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = quan_bh_15)

quan_table_15 <- extractResults(quan_result_15, 
                               data_name = 'quan_15')

# quan intersection
quan_result_int <- runModels(quan_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = quan_int_feat)

quan_table_int <- extractResults(quan_result_int, 
                                data_name = 'quan_int')

# quan union
quan_result_union <- runModels(quan_cases, 
                              bump_hunter = T, 
                              bump_hunter_data = quan_union_feat)

quan_table_union <- extractResults(quan_result_union, 
                                  data_name = 'quan_union')


##########
# threshold 0.07 - 0.15 sig
##########

# quan sig 07 
quan_result_sig_07 <- runModels(quan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = quan_bh_sig_07)

quan_table_sig_07 <- extractResults(quan_result_sig_07, 
                                   data_name = 'quan_sig_07')

# quan sig 08
quan_result_sig_08 <- runModels(quan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = quan_bh_sig_08)

quan_table_sig_08 <- extractResults(quan_result_sig_08, 
                                   data_name = 'quan_sig_08')

# quan sig 09 
quan_result_sig_09 <- runModels(quan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = quan_bh_sig_09)

quan_table_sig_09 <- extractResults(quan_result_sig_09, 
                                   data_name = 'quan_sig_09')

# quan sig 10 
quan_result_sig_10 <- runModels(quan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = quan_bh_sig_10)

quan_table_sig_10 <- extractResults(quan_result_sig_10, 
                                   data_name = 'quan_sig_10')

# quan sig 11 
quan_result_sig_11 <- runModels(quan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = quan_bh_sig_11)

quan_table_sig_11 <- extractResults(quan_result_sig_11, 
                                   data_name = 'quan_sig_11')

# quan sig 12 
quan_result_sig_12 <- runModels(quan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = quan_bh_sig_12)

quan_table_sig_12 <- extractResults(quan_result_sig_12, 
                                   data_name = 'quan_sig_12')

# quan sig 13 
quan_result_sig_13 <- runModels(quan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = quan_bh_sig_13)

quan_table_sig_13 <- extractResults(quan_result_sig_13, 
                                   data_name = 'quan_sig_13')

# quan sig 14 
quan_result_sig_14 <- runModels(quan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = quan_bh_sig_14)

quan_table_sig_14 <- extractResults(quan_result_sig_14, 
                                   data_name = 'quan_sig_14')

# quan sig 15 
quan_result_sig_15 <- runModels(quan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = quan_bh_sig_15)

quan_table_sig_15 <- extractResults(quan_result_sig_15, 
                                   data_name = 'quan_sig_15')


# quan intersection
quan_result_sig_int <- runModels(quan_cases, 
                                bump_hunter = T, 
                                bump_hunter_data = quan_int_sig_feat)

quan_table_sig_int <- extractResults(quan_result_sig_int, 
                                    data_name = 'quan_sig_int')

# quan union
quan_result_sig_union <- runModels(quan_cases, 
                                  bump_hunter = T, 
                                  bump_hunter_data = quan_union_sig_feat)

quan_table_sig_union <- extractResults(quan_result_sig_union, 
                                      data_name = 'quan_sig_union')


###########
# rbind tables and save RDA file
###########
quan_table <- rbind(quan_table_07, quan_table_sig_07,
                   quan_table_08, quan_table_sig_08,
                   quan_table_09, quan_table_sig_09,
                   quan_table_10, quan_table_sig_10,
                   quan_table_11, quan_table_sig_11,
                   quan_table_12, quan_table_sig_12,
                   quan_table_13, quan_table_sig_13,
                   quan_table_14, quan_table_sig_14,
                   quan_table_15, quan_table_sig_15,
                   quan_table_int, quan_table_sig_int,
                   quan_table_union, quan_table_union)

# remove data 
rm(quan_table_07, quan_table_sig_07,
   quan_table_08, quan_table_sig_08,
   quan_table_09, quan_table_sig_09,
   quan_table_10, quan_table_sig_10,
   quan_table_11, quan_table_sig_11,
   quan_table_12, quan_table_sig_12,
   quan_table_13, quan_table_sig_13,
   quan_table_14, quan_table_sig_14,
   quan_table_15, quan_table_sig_15,
   quan_table_int, quan_table_sig_int,
   quan_table_union, quan_table_sig_union)


#save table 
saveRDS(quan_table, 
        file = paste0(quan_folder, '/quan_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(quan_result_07, 
        file = paste0(quan_folder, '/quan_result_07.rda'))
saveRDS(quan_result_sig_07, 
        file = paste0(quan_folder, '/quan_result_sig_07.rda'))

saveRDS(quan_result_08, 
        file = paste0(quan_folder, '/quan_result_08.rda'))
saveRDS(quan_result_sig_08, 
        file = paste0(quan_folder, '/quan_result_sig_08.rda'))

saveRDS(quan_result_09, 
        file = paste0(quan_folder, '/quan_result_09.rda'))
saveRDS(quan_result_sig_09, 
        file = paste0(quan_folder, '/quan_result_sig_09.rda'))

saveRDS(quan_result_10, 
        file = paste0(quan_folder, '/quan_result_10.rda'))
saveRDS(quan_result_sig_10, 
        file = paste0(quan_folder, '/quan_result_sig_10.rda'))

saveRDS(quan_result_11, 
        file = paste0(quan_folder, '/quan_result_11.rda'))
saveRDS(quan_result_sig_11, 
        file = paste0(quan_folder, '/quan_result_sig_11.rda'))

saveRDS(quan_result_12, 
        file = paste0(quan_folder, '/quan_result_12.rda'))
saveRDS(quan_result_sig_12, 
        file = paste0(quan_folder, '/quan_result_sig_12.rda'))

saveRDS(quan_result_13, 
        file = paste0(quan_folder, '/quan_result_13.rda'))
saveRDS(quan_result_sig_13, 
        file = paste0(quan_folder, '/quan_result_sig_13.rda'))

saveRDS(quan_result_14, 
        file = paste0(quan_folder, '/quan_result_14.rda'))
saveRDS(quan_result_sig_14, 
        file = paste0(quan_folder, '/quan_result_sig_14.rda'))

saveRDS(quan_result_15, 
        file = paste0(quan_folder, '/quan_result_15.rda'))
saveRDS(quan_result_sig_15, 
        file = paste0(quan_folder, '/quan_result_sig_15.rda'))

saveRDS(quan_result_int, 
        file = paste0(quan_folder, '/quan_result_int.rda'))
saveRDS(quan_result_sig_int, 
        file = paste0(quan_folder, '/quan_result_sig_int.rda'))

saveRDS(quan_result_union, 
        file = paste0(quan_folder, '/quan_result_union.rda'))
saveRDS(quan_result_sig_union, 
        file = paste0(quan_folder, '/quan_result_sig_union.rda'))


