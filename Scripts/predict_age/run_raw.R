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
# remove quan, swan, and funnorm objects
##########
rm(list = ls(pattern = "quan_*"))
rm(list = ls(pattern = "swan_*"))
rm(list = ls(pattern = "funnorm_*"))


###################################################################################################################################
## beta_raw

##########
# threshold 0.07 - 0.15 normal
##########

# raw 07 
raw_result_07 <- runModels(raw_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = raw_bh_07)

raw_table_07 <- extractResults(raw_result_07, 
                               data_name = 'raw_07')

# raw 08
raw_result_08 <- runModels(raw_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = raw_bh_08)

raw_table_08 <- extractResults(raw_result_08, 
                               data_name = 'raw_08')

# raw 09 
raw_result_09 <- runModels(raw_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = raw_bh_09)

raw_table_09 <- extractResults(raw_result_09, 
                               data_name = 'raw_09')

# raw 10 
raw_result_10 <- runModels(raw_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = raw_bh_10)

raw_table_10 <- extractResults(raw_result_10, 
                               data_name = 'raw_10')

# raw 11
raw_result_11 <- runModels(raw_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = raw_bh_11)

raw_table_11 <- extractResults(raw_result_11, 
                               data_name = 'raw_11')

# raw 12
raw_result_12 <- runModels(raw_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = raw_bh_12)

raw_table_12 <- extractResults(raw_result_12, 
                               data_name = 'raw_12')

# raw 13
raw_result_13 <- runModels(raw_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = raw_bh_13)

raw_table_13 <- extractResults(raw_result_13, 
                               data_name = 'raw_13')

# raw 14
raw_result_14 <- runModels(raw_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = raw_bh_14)

raw_table_14 <- extractResults(raw_result_14, 
                               data_name = 'raw_14')

# raw 15
raw_result_15 <- runModels(raw_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = raw_bh_15)

raw_table_15 <- extractResults(raw_result_15, 
                               data_name = 'raw_15')

# raw intersection
raw_result_int <- runModels(raw_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = raw_int_feat)

raw_table_int <- extractResults(raw_result_int, 
                                data_name = 'raw_int')

# raw union
raw_result_union <- runModels(raw_cases, 
                              bump_hunter = T, 
                              bump_hunter_data = raw_union_feat)

raw_table_union <- extractResults(raw_result_union, 
                                  data_name = 'raw_union')


##########
# threshold 0.07 - 0.15 sig
##########

# raw sig 07 
raw_result_sig_07 <- runModels(raw_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = raw_bh_sig_07)

raw_table_sig_07 <- extractResults(raw_result_sig_07, 
                                   data_name = 'raw_sig_07')

# raw sig 08
raw_result_sig_08 <- runModels(raw_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = raw_bh_sig_08)

raw_table_sig_08 <- extractResults(raw_result_sig_08, 
                                   data_name = 'raw_sig_08')

# raw sig 09 
raw_result_sig_09 <- runModels(raw_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = raw_bh_sig_09)

raw_table_sig_09 <- extractResults(raw_result_sig_09, 
                                   data_name = 'raw_sig_09')

# raw sig 10 
raw_result_sig_10 <- runModels(raw_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = raw_bh_sig_10)

raw_table_sig_10 <- extractResults(raw_result_sig_10, 
                                   data_name = 'raw_sig_10')

# raw sig 11 
raw_result_sig_11 <- runModels(raw_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = raw_bh_sig_11)

raw_table_sig_11 <- extractResults(raw_result_sig_11, 
                                   data_name = 'raw_sig_11')

# raw sig 12 
raw_result_sig_12 <- runModels(raw_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = raw_bh_sig_12)

raw_table_sig_12 <- extractResults(raw_result_sig_12, 
                                   data_name = 'raw_sig_12')

# raw sig 13 
raw_result_sig_13 <- runModels(raw_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = raw_bh_sig_13)

raw_table_sig_13 <- extractResults(raw_result_sig_13, 
                                   data_name = 'raw_sig_13')

# raw sig 14 
raw_result_sig_14 <- runModels(raw_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = raw_bh_sig_14)

raw_table_sig_14 <- extractResults(raw_result_sig_14, 
                                   data_name = 'raw_sig_14')

# raw sig 15 
raw_result_sig_15 <- runModels(raw_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = raw_bh_sig_15)

raw_table_sig_15 <- extractResults(raw_result_sig_15, 
                                   data_name = 'raw_sig_15')


# raw intersection
raw_result_sig_int <- runModels(raw_cases, 
                                bump_hunter = T, 
                                bump_hunter_data = raw_int_sig_feat)

raw_table_sig_int <- extractResults(raw_result_sig_int, 
                                    data_name = 'raw_sig_int')

# raw union
raw_result_sig_union <- runModels(raw_cases, 
                                  bump_hunter = T, 
                                  bump_hunter_data = raw_union_sig_feat)

raw_table_sig_union <- extractResults(raw_result_sig_union, 
                                      data_name = 'raw_sig_union')


###########
# rbind tables and save RDA file
###########
raw_table <- rbind(raw_table_07, raw_table_sig_07,
                   raw_table_08, raw_table_sig_08,
                   raw_table_09, raw_table_sig_09,
                   raw_table_10, raw_table_sig_10,
                   raw_table_11, raw_table_sig_11,
                   raw_table_12, raw_table_sig_12,
                   raw_table_13, raw_table_sig_13,
                   raw_table_14, raw_table_sig_14,
                   raw_table_15, raw_table_sig_15,
                   raw_table_int, raw_table_sig_int,
                   raw_table_union, raw_table_sig_union)

# remove data 
rm(raw_table_07, raw_table_sig_07,
   raw_table_08, raw_table_sig_08,
   raw_table_09, raw_table_sig_09,
   raw_table_10, raw_table_sig_10,
   raw_table_11, raw_table_sig_11,
   raw_table_12, raw_table_sig_12,
   raw_table_13, raw_table_sig_13,
   raw_table_14, raw_table_sig_14,
   raw_table_15, raw_table_sig_15,
   raw_table_int, raw_table_sig_int,
   raw_table_union, raw_table_sig_union)


#save table 
saveRDS(raw_table, 
        file = paste0(raw_folder, '/raw_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(raw_result_07, 
        file = paste0(raw_folder, '/raw_result_07.rda'))
saveRDS(raw_result_sig_07, 
        file = paste0(raw_folder, '/raw_result_sig_07.rda'))

saveRDS(raw_result_08, 
        file = paste0(raw_folder, '/raw_result_08.rda'))
saveRDS(raw_result_sig_08, 
        file = paste0(raw_folder, '/raw_result_sig_08.rda'))

saveRDS(raw_result_09, 
        file = paste0(raw_folder, '/raw_result_09.rda'))
saveRDS(raw_result_sig_09, 
        file = paste0(raw_folder, '/raw_result_sig_09.rda'))

saveRDS(raw_result_10, 
        file = paste0(raw_folder, '/raw_result_10.rda'))
saveRDS(raw_result_sig_10, 
        file = paste0(raw_folder, '/raw_result_sig_10.rda'))

saveRDS(raw_result_11, 
        file = paste0(raw_folder, '/raw_result_11.rda'))
saveRDS(raw_result_sig_11, 
        file = paste0(raw_folder, '/raw_result_sig_11.rda'))

saveRDS(raw_result_12, 
        file = paste0(raw_folder, '/raw_result_12.rda'))
saveRDS(raw_result_sig_12, 
        file = paste0(raw_folder, '/raw_result_sig_12.rda'))

saveRDS(raw_result_13, 
        file = paste0(raw_folder, '/raw_result_13.rda'))
saveRDS(raw_result_sig_13, 
        file = paste0(raw_folder, '/raw_result_sig_13.rda'))

saveRDS(raw_result_14, 
        file = paste0(raw_folder, '/raw_result_14.rda'))
saveRDS(raw_result_sig_14, 
        file = paste0(raw_folder, '/raw_result_sig_14.rda'))

saveRDS(raw_result_15, 
        file = paste0(raw_folder, '/raw_result_15.rda'))
saveRDS(raw_result_sig_15, 
        file = paste0(raw_folder, '/raw_result_sig_15.rda'))

saveRDS(raw_result_int, 
        file = paste0(raw_folder, '/raw_result_inr.rda'))
saveRDS(raw_result_sig_int, 
        file = paste0(raw_folder, '/raw_result_sig_int.rda'))

saveRDS(raw_result_union, 
        file = paste0(raw_folder, '/raw_result_union.rda'))
saveRDS(raw_result_sig_union, 
        file = paste0(raw_folder, '/raw_result_sig_union.rda'))


