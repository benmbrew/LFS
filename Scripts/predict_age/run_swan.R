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
# remove swan, swan, and swan objects
##########
rm(list = ls(pattern = "quan_*"))
rm(list = ls(pattern = "funnorm_*"))
rm(list = ls(pattern = "raw_*"))

###################################################################################################################################
## beta_swan

##########
# threshold 0.10, 0.15, 0.20, 0.30
##########

#NORMAL
# swan 10 
swan_result_10 <- runModels(swan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = swan_bh_10)

swan_table_10 <- extractResults(swan_result_10, 
                                   data_name = 'swan_10')


# swan 15
swan_result_15 <- runModels(swan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = swan_bh_15)

swan_table_15 <- extractResults(swan_result_15, 
                                   data_name = 'swan_15')

# swan 20
swan_result_20 <- runModels(swan_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = swan_bh_20)

swan_table_20 <- extractResults(swan_result_20, 
                                   data_name = 'swan_20')


#UNBAL
# swan 10 
swan_result_unbal_10 <- runModels(swan_cases, 
                                     bump_hunter = T, 
                                     bump_hunter_data = swan_unbal_bh_10)

swan_table_unbal_10 <- extractResults(swan_result_unbal_10, 
                                         data_name = 'swan_unbal_10')


# swan 15
swan_result_unbal_15 <- runModels(swan_cases, 
                                     bump_hunter = T, 
                                     bump_hunter_data = swan_unbal_bh_15)

swan_table_unbal_15 <- extractResults(swan_result_unbal_15, 
                                         data_name = 'swan_unbal_15')

# swan 20
swan_result_unbal_20 <- runModels(swan_cases, 
                                     bump_hunter = T, 
                                     bump_hunter_data = swan_unbal_bh_20)

swan_table_unbal_20 <- extractResults(swan_result_unbal_20, 
                                         data_name = 'swan_unbal_20')

# swan intersection
swan_result_int <- runModels(swan_cases, 
                                bump_hunter = T, 
                                bump_hunter_data = swan_int_feat)

swan_table_int <- extractResults(swan_result_int, 
                                    data_name = 'swan_int')

swan_result_unbal_int <- runModels(swan_cases, 
                                      bump_hunter = T, 
                                      bump_hunter_data = swan_unbal_int_feat)

swan_table_unbal_int <- extractResults(swan_result_unbal_int, 
                                          data_name = 'swan_unbal_int')



##########
# threshold 0.07 - 0.15 sig
##########

#NORMAL
# swan 10 
swan_sig_result_10 <- runModels(swan_cases, 
                                   bump_hunter = T, 
                                   bump_hunter_data = swan_bh_sig_10)

swan_sig_table_10 <- extractResults(swan_sig_result_10, 
                                       data_name = 'swan_sig_10')


# swan 15
swan_sig_result_15 <- runModels(swan_cases, 
                                   bump_hunter = T, 
                                   bump_hunter_data = swan_bh_sig_15)

swan_sig_table_15 <- extractResults(swan_sig_result_15, 
                                       data_name = 'swan_sig_15')

# swan 20
swan_sig_result_20 <- runModels(swan_cases, 
                                   bump_hunter = T, 
                                   bump_hunter_data = swan_bh_sig_20)

swan_sig_table_20 <- extractResults(swan_sig_result_20, 
                                       data_name = 'swan_sig_20')



# swan intersection
swan_sig_result_int <- runModels(swan_cases, 
                                    bump_hunter = T, 
                                    bump_hunter_data = swan_int_sig_feat)

swan_sig_table_int <- extractResults(swan_sig_result_int, 
                                        data_name = 'swan_sig_int')



#UNBAL
# swan 10 
swan_sig_result_unbal_10 <- runModels(swan_cases, 
                                         bump_hunter = T, 
                                         bump_hunter_data = swan_unbal_bh_sig_10)

swan_sig_table_unbal_10 <- extractResults(swan_sig_result_unbal_10, 
                                             data_name = 'swan_unbal_sig_10')


# swan 15
swan_sig_result_unbal_15 <- runModels(swan_cases, 
                                         bump_hunter = T, 
                                         bump_hunter_data = swan_unbal_bh_sig_15)

swan_sig_table_unbal_15 <- extractResults(swan_sig_result_unbal_15, 
                                             data_name = 'swan_unbal_sig_15')

# # swan 20
swan_sig_result_unbal_20 <- runModels(swan_cases,
                                         bump_hunter = T,
                                         bump_hunter_data = swan_unbal_bh_sig_20)

swan_sig_table_unbal_20 <- extractResults(swan_sig_result_unbal_20,
                                             data_name = 'swan_unbal_sig_20')


# # # swan intersection
swan_sig_result_unbal_int <- runModels(swan_cases,
                                          bump_hunter = T,
                                          bump_hunter_data = swan_unbal_int_sig_feat)
#
swan_sig_table_unbal_int <- extractResults(swan_sig_result_unbal_int,
                                              data_name = 'swan_unbal_sig_int')


###########
# rbind tables and save RDA file
###########
swan_table <- rbind(swan_table_10, swan_table_15,
                       swan_table_20,
                       swan_sig_table_10, swan_sig_table_15,
                       swan_sig_table_20,
                       swan_table_int, swan_sig_table_int,
                       swan_table_unbal_10, swan_table_unbal_15,
                       swan_table_unbal_20, 
                       swan_sig_table_unbal_10, swan_sig_table_unbal_15,
                       swan_sig_table_unbal_20, 
                       swan_table_unbal_int, 
                       swan_sig_table_unbal_int)



# # remove data 
rm(swan_table_10, swan_table_15,
   swan_table_20,
   swan_sig_table_10, swan_sig_table_15,
   swan_sig_table_20,
   swan_table_int, swan_sig_table_int,
   swan_table_unbal_10, swan_table_unbal_15,
   swan_table_unbal_20, 
   swan_sig_table_unbal_10, swan_sig_table_unbal_15,
   swan_sig_table_unbal_20, 
   swan_table_unbal_int, 
   swan_sig_table_unbal_int)

#save table 
saveRDS(swan_table, 
        file = paste0(swan_folder, '/swan_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(swan_result_10, 
        file = paste0(swan_folder, '/swan_result_10.rda'))

saveRDS(swan_result_15, 
        file = paste0(swan_folder, '/swan_result_15.rda'))

saveRDS(swan_result_20, 
        file = paste0(swan_folder, '/swan_result_20.rda'))

saveRDS(swan_result_unbal_10, 
        file = paste0(swan_folder, '/swan_result_unbal_10.rda'))

saveRDS(swan_result_unbal_15, 
        file = paste0(swan_folder, '/swan_result_unbal_15.rda'))

saveRDS(swan_result_unbal_20, 
        file = paste0(swan_folder, '/swan_result_unbal_20.rda'))

