####### This script will predict age of onset with swan preprocessed data
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
# remove raw, quan, and funnorm objects
##########
rm(list = ls(pattern = "raw_*"))
rm(list = ls(pattern = "quan_*"))
rm(list = ls(pattern = "funnorm_*"))


###################################################################################################################################
## beta_swan

##########
# threshold 0.07 - 0.15 normal
##########

# swan 07 
swan_result_07 <- runModels(swan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = swan_bh_07)

swan_table_07 <- extractResults(swan_result_07, 
                               data_name = 'swan_07')

# swan 08
swan_result_08 <- runModels(swan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = swan_bh_08)

swan_table_08 <- extractResults(swan_result_08, 
                               data_name = 'swan_08')

# swan 09 
swan_result_09 <- runModels(swan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = swan_bh_09)

swan_table_09 <- extractResults(swan_result_09, 
                               data_name = 'swan_09')

# swan 10 
swan_result_10 <- runModels(swan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = swan_bh_10)

swan_table_10 <- extractResults(swan_result_10, 
                               data_name = 'swan_10')

# swan 11
swan_result_11 <- runModels(swan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = swan_bh_11)

swan_table_11 <- extractResults(swan_result_11, 
                               data_name = 'swan_11')

# swan 12
swan_result_12 <- runModels(swan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = swan_bh_12)

swan_table_12 <- extractResults(swan_result_12, 
                               data_name = 'swan_12')

# swan 13
swan_result_13 <- runModels(swan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = swan_bh_13)

swan_table_13 <- extractResults(swan_result_13, 
                               data_name = 'swan_13')

# swan 14
swan_result_14 <- runModels(swan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = swan_bh_14)

swan_table_14 <- extractResults(swan_result_14, 
                               data_name = 'swan_14')

# swan 15
swan_result_15 <- runModels(swan_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = swan_bh_15)

swan_table_15 <- extractResults(swan_result_15, 
                               data_name = 'swan_15')

# swan intersection
swan_result_int <- runModels(swan_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = swan_int_feat)

swan_table_int <- extractResults(swan_result_int, 
                                data_name = 'swan_int')

# swan union
swan_result_union <- runModels(swan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = swan_union_feat)

swan_table_union <- extractResults(swan_result_union, 
                                   data_name = 'swan_union')


##########
# threshold 0.07 - 0.15 sig
##########

# swan sig 07 
swan_result_sig_07 <- runModels(swan_cases, 
                                bump_hunter = T, 
                                bump_hunter_data = swan_bh_sig_07)

swan_table_sig_07 <- extractResults(swan_result_sig_07, 
                                    data_name = 'swan_sig_07')

# swan sig 08
swan_result_sig_08 <- runModels(swan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = swan_bh_sig_08)

swan_table_sig_08 <- extractResults(swan_result_sig_08, 
                                   data_name = 'swan_sig_08')

# swan sig 09 
swan_result_sig_09 <- runModels(swan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = swan_bh_sig_09)

swan_table_sig_09 <- extractResults(swan_result_sig_09, 
                                   data_name = 'swan_sig_09')

# swan sig 10 
swan_result_sig_10 <- runModels(swan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = swan_bh_sig_10)

swan_table_sig_10 <- extractResults(swan_result_sig_10, 
                                   data_name = 'swan_sig_10')

# swan sig 11 
swan_result_sig_11 <- runModels(swan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = swan_bh_sig_11)

swan_table_sig_11 <- extractResults(swan_result_sig_11, 
                                   data_name = 'swan_sig_11')

# swan sig 12 
swan_result_sig_12 <- runModels(swan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = swan_bh_sig_12)

swan_table_sig_12 <- extractResults(swan_result_sig_12, 
                                   data_name = 'swan_sig_12')

# swan sig 13 
swan_result_sig_13 <- runModels(swan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = swan_bh_sig_13)

swan_table_sig_13 <- extractResults(swan_result_sig_13, 
                                   data_name = 'swan_sig_13')

# swan sig 14 
swan_result_sig_14 <- runModels(swan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = swan_bh_sig_14)

swan_table_sig_14 <- extractResults(swan_result_sig_14, 
                                   data_name = 'swan_sig_14')

# swan sig 15 
swan_result_sig_15 <- runModels(swan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = swan_bh_sig_15)

swan_table_sig_15 <- extractResults(swan_result_sig_15, 
                                   data_name = 'swan_sig_15')


# swan intersection
swan_result_sig_int <- runModels(swan_cases, 
                                bump_hunter = T, 
                                bump_hunter_data = swan_int_sig_feat)

swan_table_sig_int <- extractResults(swan_result_sig_int, 
                                    data_name = 'swan_sig_int')

# swan union
swan_result_sig_union <- runModels(swan_cases, 
                                  bump_hunter = T, 
                                  bump_hunter_data = swan_union_sig_feat)

swan_table_sig_union <- extractResults(swan_result_sig_union, 
                                      data_name = 'swan_sig_union')



###########
# rbind tables and save RDA file
###########
swan_table <- rbind(swan_table_07, swan_table_sig_07,
                   swan_table_08, swan_table_sig_08,
                   swan_table_09, swan_table_sig_09,
                   swan_table_10, swan_table_sig_10,
                   swan_table_11, swan_table_sig_11,
                   swan_table_12, swan_table_sig_12,
                   swan_table_13, swan_table_sig_13,
                   swan_table_14, swan_table_sig_14,
                   swan_table_15, swan_table_sig_15,
                   swan_table_int, swan_table_sig_int,
                   swan_table_union, swan_table_sig_union)

# remove data 
rm(swan_table_07, swan_table_sig_07,
   swan_table_08, swan_table_sig_08,
   swan_table_09, swan_table_sig_09,
   swan_table_10, swan_table_sig_10,
   swan_table_11, swan_table_sig_11,
   swan_table_12, swan_table_sig_12,
   swan_table_13, swan_table_sig_13,
   swan_table_14, swan_table_sig_14,
   swan_table_15, swan_table_sig_15,
   swan_table_int, swan_table_sig_int,
   swan_table_union, swan_table_sig_union)


#save table 
saveRDS(swan_table, 
        file = paste0(swan_folder, '/swan_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(swan_result_07, 
        file = paste0(swan_folder, '/swan_result_07.rda'))
saveRDS(swan_result_sig_07, 
        file = paste0(swan_folder, '/swan_result_sig_07.rda'))

saveRDS(swan_result_08, 
        file = paste0(swan_folder, '/swan_result_08.rda'))
saveRDS(swan_result_sig_08, 
        file = paste0(swan_folder, '/swan_result_sig_08.rda'))

saveRDS(swan_result_09, 
        file = paste0(swan_folder, '/swan_result_09.rda'))
saveRDS(swan_result_sig_09, 
        file = paste0(swan_folder, '/swan_result_sig_09.rda'))

saveRDS(swan_result_10, 
        file = paste0(swan_folder, '/swan_result_10.rda'))
saveRDS(swan_result_sig_10, 
        file = paste0(swan_folder, '/swan_result_sig_10.rda'))

saveRDS(swan_result_11, 
        file = paste0(swan_folder, '/swan_result_11.rda'))
saveRDS(swan_result_sig_11, 
        file = paste0(swan_folder, '/swan_result_sig_11.rda'))

saveRDS(swan_result_12, 
        file = paste0(swan_folder, '/swan_result_12.rda'))
saveRDS(swan_result_sig_12, 
        file = paste0(swan_folder, '/swan_result_sig_12.rda'))

saveRDS(swan_result_13, 
        file = paste0(swan_folder, '/swan_result_13.rda'))
saveRDS(swan_result_sig_13, 
        file = paste0(swan_folder, '/swan_result_sig_13.rda'))

saveRDS(swan_result_14, 
        file = paste0(swan_folder, '/swan_result_14.rda'))
saveRDS(swan_result_sig_14, 
        file = paste0(swan_folder, '/swan_result_sig_14.rda'))

saveRDS(swan_result_15, 
        file = paste0(swan_folder, '/swan_result_15.rda'))
saveRDS(swan_result_sig_15, 
        file = paste0(swan_folder, '/swan_result_sig_15.rda'))

saveRDS(swan_result_int, 
        file = paste0(swan_folder, '/swan_result_inr.rda'))
saveRDS(swan_result_sig_int, 
        file = paste0(swan_folder, '/swan_result_sig_int.rda'))

saveRDS(swan_result_union, 
        file = paste0(swan_folder, '/swan_result_union.rda'))
saveRDS(swan_result_sig_union, 
        file = paste0(swan_folder, '/swan_result_sig_union.rda'))


