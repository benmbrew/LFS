####### This script will predict age of onset with funnorm preprocessed data
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
# remove raw, swan, and quan objects
##########
rm(list = ls(pattern = "raw_*"))
rm(list = ls(pattern = "swan_*"))
rm(list = ls(pattern = "quan_*"))


###################################################################################################################################
## beta_funnorm

##########
# threshold 0.07 - 0.15 normal
##########

# funnorm 07 
funnorm_result_07 <- runModels(funnorm_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = funnorm_bh_07)

funnorm_table_07 <- extractResults(funnorm_result_07, 
                               data_name = 'funnorm_07')

# funnorm 08
funnorm_result_08 <- runModels(funnorm_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = funnorm_bh_08)

funnorm_table_08 <- extractResults(funnorm_result_08, 
                               data_name = 'funnorm_08')

# funnorm 09 
funnorm_result_09 <- runModels(funnorm_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = funnorm_bh_09)

funnorm_table_09 <- extractResults(funnorm_result_09, 
                               data_name = 'funnorm_09')

# funnorm 10 
funnorm_result_10 <- runModels(funnorm_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = funnorm_bh_10)

funnorm_table_10 <- extractResults(funnorm_result_10, 
                               data_name = 'funnorm_10')

# funnorm 11
funnorm_result_11 <- runModels(funnorm_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = funnorm_bh_11)

funnorm_table_11 <- extractResults(funnorm_result_11, 
                               data_name = 'funnorm_11')

# funnorm 12
funnorm_result_12 <- runModels(funnorm_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = funnorm_bh_12)

funnorm_table_12 <- extractResults(funnorm_result_12, 
                               data_name = 'funnorm_12')

# funnorm 13
funnorm_result_13 <- runModels(funnorm_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = funnorm_bh_13)

funnorm_table_13 <- extractResults(funnorm_result_13, 
                               data_name = 'funnorm_13')

# funnorm 14
funnorm_result_14 <- runModels(funnorm_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = funnorm_bh_14)

funnorm_table_14 <- extractResults(funnorm_result_14, 
                               data_name = 'funnorm_14')

# funnorm 15
funnorm_result_15 <- runModels(funnorm_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = funnorm_bh_15)

funnorm_table_15 <- extractResults(funnorm_result_15, 
                               data_name = 'funnorm_15')

# funnorm intersection
funnorm_result_int <- runModels(funnorm_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = funnorm_int_feat)

funnorm_table_int <- extractResults(funnorm_result_int, 
                                data_name = 'funnorm_int')

# funnorm union
funnorm_result_union <- runModels(funnorm_cases, 
                              bump_hunter = T, 
                              bump_hunter_data = funnorm_union_feat)

funnorm_table_union <- extractResults(funnorm_result_union, 
                                  data_name = 'funnorm_union')


##########
# threshold 0.07 - 0.15 sig
##########

# funnorm sig 07 
funnorm_result_sig_07 <- runModels(funnorm_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = funnorm_bh_sig_07)

funnorm_table_sig_07 <- extractResults(funnorm_result_sig_07, 
                                   data_name = 'funnorm_sig_07')

# funnorm sig 08
funnorm_result_sig_08 <- runModels(funnorm_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = funnorm_bh_sig_08)

funnorm_table_sig_08 <- extractResults(funnorm_result_sig_08, 
                                   data_name = 'funnorm_sig_08')

# funnorm sig 09 
funnorm_result_sig_09 <- runModels(funnorm_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = funnorm_bh_sig_09)

funnorm_table_sig_09 <- extractResults(funnorm_result_sig_09, 
                                   data_name = 'funnorm_sig_09')

# funnorm sig 10 
funnorm_result_sig_10 <- runModels(funnorm_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = funnorm_bh_sig_10)

funnorm_table_sig_10 <- extractResults(funnorm_result_sig_10, 
                                   data_name = 'funnorm_sig_10')

# funnorm sig 11 
funnorm_result_sig_11 <- runModels(funnorm_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = funnorm_bh_sig_11)

funnorm_table_sig_11 <- extractResults(funnorm_result_sig_11, 
                                   data_name = 'funnorm_sig_11')

# funnorm sig 12 
funnorm_result_sig_12 <- runModels(funnorm_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = funnorm_bh_sig_12)

funnorm_table_sig_12 <- extractResults(funnorm_result_sig_12, 
                                   data_name = 'funnorm_sig_12')

# funnorm sig 13 
funnorm_result_sig_13 <- runModels(funnorm_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = funnorm_bh_sig_13)

funnorm_table_sig_13 <- extractResults(funnorm_result_sig_13, 
                                   data_name = 'funnorm_sig_13')

# funnorm sig 14 
funnorm_result_sig_14 <- runModels(funnorm_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = funnorm_bh_sig_14)

funnorm_table_sig_14 <- extractResults(funnorm_result_sig_14, 
                                   data_name = 'funnorm_sig_14')

# funnorm sig 15 
funnorm_result_sig_15 <- runModels(funnorm_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = funnorm_bh_sig_15)

funnorm_table_sig_15 <- extractResults(funnorm_result_sig_15, 
                                   data_name = 'funnorm_sig_15')


# funnorm intersection
funnorm_result_sig_int <- runModels(funnorm_cases, 
                                bump_hunter = T, 
                                bump_hunter_data = funnorm_int_sig_feat)

funnorm_table_sig_int <- extractResults(funnorm_result_sig_int, 
                                    data_name = 'funnorm_sig_int')

# funnorm union
funnorm_result_sig_union <- runModels(funnorm_cases, 
                                  bump_hunter = T, 
                                  bump_hunter_data = funnorm_union_sig_feat)

funnorm_table_sig_union <- extractResults(funnorm_result_sig_union, 
                                      data_name = 'funnorm_sig_union')


###########
# rbind tables and save RDA file
###########
funnorm_table <- rbind(funnorm_table_07, funnorm_table_sig_07,
                   funnorm_table_08, funnorm_table_sig_08,
                   funnorm_table_09, funnorm_table_sig_09,
                   funnorm_table_10, funnorm_table_sig_10,
                   funnorm_table_11, funnorm_table_sig_11,
                   funnorm_table_12, funnorm_table_sig_12,
                   funnorm_table_13, funnorm_table_sig_13,
                   funnorm_table_14, funnorm_table_sig_14,
                   funnorm_table_15, funnorm_table_sig_15,
                   funnorm_table_int, funnorm_table_sig_int,
                   funnorm_table_union, funnorm_table_sig_union)

# remove data 
rm(funnorm_table_07, funnorm_table_sig_07,
   funnorm_table_08, funnorm_table_sig_08,
   funnorm_table_09, funnorm_table_sig_09,
   funnorm_table_10, funnorm_table_sig_10,
   funnorm_table_11, funnorm_table_sig_11,
   funnorm_table_12, funnorm_table_sig_12,
   funnorm_table_13, funnorm_table_sig_13,
   funnorm_table_14, funnorm_table_sig_14,
   funnorm_table_15, funnorm_table_sig_15,
   funnorm_table_int, funnorm_table_sig_int,
   funnorm_table_union, funnorm_table_sig_union)


#save table 
saveRDS(funnorm_table, 
        file = paste0(funnorm_folder, '/funnorm_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(funnorm_result_07, 
        file = paste0(funnorm_folder, '/funnorm_result_07.rda'))
saveRDS(funnorm_result_sig_07, 
        file = paste0(funnorm_folder, '/funnorm_result_sig_07.rda'))

saveRDS(funnorm_result_08, 
        file = paste0(funnorm_folder, '/funnorm_result_08.rda'))
saveRDS(funnorm_result_sig_08, 
        file = paste0(funnorm_folder, '/funnorm_result_sig_08.rda'))

saveRDS(funnorm_result_09, 
        file = paste0(funnorm_folder, '/funnorm_result_09.rda'))
saveRDS(funnorm_result_sig_09, 
        file = paste0(funnorm_folder, '/funnorm_result_sig_09.rda'))

saveRDS(funnorm_result_10, 
        file = paste0(funnorm_folder, '/funnorm_result_10.rda'))
saveRDS(funnorm_result_sig_10, 
        file = paste0(funnorm_folder, '/funnorm_result_sig_10.rda'))

saveRDS(funnorm_result_11, 
        file = paste0(funnorm_folder, '/funnorm_result_11.rda'))
saveRDS(funnorm_result_sig_11, 
        file = paste0(funnorm_folder, '/funnorm_result_sig_11.rda'))

saveRDS(funnorm_result_12, 
        file = paste0(funnorm_folder, '/funnorm_result_12.rda'))
saveRDS(funnorm_result_sig_12, 
        file = paste0(funnorm_folder, '/funnorm_result_sig_12.rda'))

saveRDS(funnorm_result_13, 
        file = paste0(funnorm_folder, '/funnorm_result_13.rda'))
saveRDS(funnorm_result_sig_13, 
        file = paste0(funnorm_folder, '/funnorm_result_sig_13.rda'))

saveRDS(funnorm_result_14, 
        file = paste0(funnorm_folder, '/funnorm_result_14.rda'))
saveRDS(funnorm_result_sig_14, 
        file = paste0(funnorm_folder, '/funnorm_result_sig_14.rda'))

saveRDS(funnorm_result_15, 
        file = paste0(funnorm_folder, '/funnorm_result_15.rda'))
saveRDS(funnorm_result_sig_15, 
        file = paste0(funnorm_folder, '/funnorm_result_sig_15.rda'))

saveRDS(funnorm_result_int, 
        file = paste0(funnorm_folder, '/funnorm_result_int.rda'))
saveRDS(funnorm_result_sig_int, 
        file = paste0(funnorm_folder, '/funnorm_result_sig_int.rda'))

saveRDS(funnorm_result_union, 
        file = paste0(funnorm_folder, '/funnorm_result_union.rda'))
saveRDS(funnorm_result_sig_union, 
        file = paste0(funnorm_folder, '/funnorm_result_sig_union.rda'))


