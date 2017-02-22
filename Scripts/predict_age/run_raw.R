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
# threshold 0.10, 0.15, 0.20, 0.30
##########

#NORMAL
# raw 10 
raw_result_10 <- runModels(raw_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = raw_bh_10)

raw_table_10 <- extractResults(raw_result_10, 
                                data_name = 'raw_10')


# raw 15
raw_result_15 <- runModels(raw_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = raw_bh_15)

raw_table_15 <- extractResults(raw_result_15, 
                                data_name = 'raw_15')

# raw 20
raw_result_20 <- runModels(raw_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = raw_bh_20)

raw_table_20 <- extractResults(raw_result_20, 
                                data_name = 'raw_20')


#UNBAL
# raw 10 
raw_result_unbal_10 <- runModels(raw_cases, 
                                  bump_hunter = T, 
                                  bump_hunter_data = raw_unbal_bh_10)

raw_table_unbal_10 <- extractResults(raw_result_unbal_10, 
                                      data_name = 'raw_unbal_10')


# raw 15
raw_result_unbal_15 <- runModels(raw_cases, 
                                  bump_hunter = T, 
                                  bump_hunter_data = raw_unbal_bh_15)

raw_table_unbal_15 <- extractResults(raw_result_unbal_15, 
                                      data_name = 'raw_unbal_15')

# raw 20
raw_result_unbal_20 <- runModels(raw_cases, 
                                  bump_hunter = T, 
                                  bump_hunter_data = raw_unbal_bh_20)

raw_table_unbal_20 <- extractResults(raw_result_unbal_20, 
                                      data_name = 'raw_unbal_20')

# raw intersection
raw_result_int <- runModels(raw_cases, 
                             bump_hunter = T, 
                             bump_hunter_data = raw_int_feat)

raw_table_int <- extractResults(raw_result_int, 
                                 data_name = 'raw_int')

raw_result_unbal_int <- runModels(raw_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = raw_unbal_int_feat)

raw_table_unbal_int <- extractResults(raw_result_unbal_int, 
                                data_name = 'raw_unbal_int')



##########
# threshold 0.07 - 0.15 sig
##########

#NORMAL
# raw 10 
raw_sig_result_10 <- runModels(raw_cases, 
                                bump_hunter = T, 
                                bump_hunter_data = raw_bh_sig_10)

raw_sig_table_10 <- extractResults(raw_sig_result_10, 
                                    data_name = 'raw_sig_10')


# raw 15
raw_sig_result_15 <- runModels(raw_cases, 
                                bump_hunter = T, 
                                bump_hunter_data = raw_bh_sig_15)

raw_sig_table_15 <- extractResults(raw_sig_result_15, 
                                    data_name = 'raw_sig_15')

# raw 20
raw_sig_result_20 <- runModels(raw_cases, 
                                bump_hunter = T, 
                                bump_hunter_data = raw_bh_sig_20)

raw_sig_table_20 <- extractResults(raw_sig_result_20, 
                                    data_name = 'raw_sig_20')



# raw intersection
raw_sig_result_int <- runModels(raw_cases, 
                                 bump_hunter = T, 
                                 bump_hunter_data = raw_int_sig_feat)

raw_sig_table_int <- extractResults(raw_sig_result_int, 
                                     data_name = 'raw_sig_int')



#UNBAL
# raw 10 
raw_sig_result_unbal_10 <- runModels(raw_cases, 
                                      bump_hunter = T, 
                                      bump_hunter_data = raw_unbal_bh_sig_10)

raw_sig_table_unbal_10 <- extractResults(raw_sig_result_unbal_10, 
                                          data_name = 'raw_unbal_sig_10')


# raw 15
raw_sig_result_unbal_15 <- runModels(raw_cases, 
                                      bump_hunter = T, 
                                      bump_hunter_data = raw_unbal_bh_sig_15)

raw_sig_table_unbal_15 <- extractResults(raw_sig_result_unbal_15, 
                                          data_name = 'raw_unbal_sig_15')

# # raw 20
raw_sig_result_unbal_20 <- runModels(raw_cases,
                                     bump_hunter = T,
                                     bump_hunter_data = raw_unbal_bh_sig_20)

raw_sig_table_unbal_20 <- extractResults(raw_sig_result_unbal_20,
                                         data_name = 'raw_unbal_sig_20')


# # # raw intersection
raw_sig_result_unbal_int <- runModels(raw_cases,
                                      bump_hunter = T,
                                      bump_hunter_data = raw_unbal_int_sig_feat)
#
raw_sig_table_unbal_int <- extractResults(raw_sig_result_unbal_int,
                                          data_name = 'raw_unbal_sig_int')


###########
# rbind tables and save RDA file
###########
raw_table <- rbind(raw_table_10, raw_table_15,
                    raw_table_20,
                    raw_sig_table_10, raw_sig_table_15,
                    raw_sig_table_20,
                    raw_table_int, raw_sig_table_int,
                    raw_table_unbal_10, raw_table_unbal_15,
                    raw_table_unbal_20, 
                    raw_sig_table_unbal_10, raw_sig_table_unbal_15,
                    raw_sig_table_unbal_20, 
                    raw_table_unbal_int, 
                    raw_sig_table_unbal_int)



# # remove data 
rm(raw_table_10, raw_table_15,
   raw_table_20,
   raw_sig_table_10, raw_sig_table_15,
   raw_sig_table_20,
   raw_table_int, raw_sig_table_int,
   raw_table_unbal_10, raw_table_unbal_15,
   raw_table_unbal_20, 
   raw_sig_table_unbal_10, raw_sig_table_unbal_15,
   raw_sig_table_unbal_20, 
   raw_table_unbal_int, 
   raw_sig_table_unbal_int)

#save table 
saveRDS(raw_table, 
        file = paste0(raw_folder, '/raw_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(raw_result_10, 
        file = paste0(raw_folder, '/raw_result_10.rda'))

saveRDS(raw_result_15, 
        file = paste0(raw_folder, '/raw_result_15.rda'))

saveRDS(raw_result_20, 
        file = paste0(raw_folder, '/raw_result_20.rda'))

saveRDS(raw_result_unbal_10, 
        file = paste0(raw_folder, '/raw_result_unbal_10.rda'))

saveRDS(raw_result_unbal_15, 
        file = paste0(raw_folder, '/raw_result_unbal_15.rda'))

saveRDS(raw_result_unbal_20, 
        file = paste0(raw_folder, '/raw_result_unbal_20.rda'))
