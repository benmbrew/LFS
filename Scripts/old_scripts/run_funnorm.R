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
# remove funnorm, swan, and funnorm objects
##########
rm(list = ls(pattern = "quan_*"))
rm(list = ls(pattern = "swan_*"))
rm(list = ls(pattern = "raw_*"))

###################################################################################################################################
## beta_funnorm

##########
# threshold 0.10, 0.15, 0.20, 0.30
##########

#NORMAL
# funnorm 10 
funnorm_result_10 <- runModels(funnorm_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = funnorm_bh_10)

funnorm_table_10 <- extractResults(funnorm_result_10, 
                               data_name = 'funnorm_10')


# funnorm 15
funnorm_result_15 <- runModels(funnorm_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = funnorm_bh_15)

funnorm_table_15 <- extractResults(funnorm_result_15, 
                               data_name = 'funnorm_15')

# funnorm 20
funnorm_result_20 <- runModels(funnorm_cases, 
                           bump_hunter = T, 
                           bump_hunter_data = funnorm_bh_20)

funnorm_table_20 <- extractResults(funnorm_result_20, 
                               data_name = 'funnorm_20')


#UNBAL
# funnorm 10 
funnorm_result_unbal_10 <- runModels(funnorm_cases, 
                                 bump_hunter = T, 
                                 bump_hunter_data = funnorm_unbal_bh_10)

funnorm_table_unbal_10 <- extractResults(funnorm_result_unbal_10, 
                                     data_name = 'funnorm_unbal_10')


# funnorm 15
funnorm_result_unbal_15 <- runModels(funnorm_cases, 
                                 bump_hunter = T, 
                                 bump_hunter_data = funnorm_unbal_bh_15)

funnorm_table_unbal_15 <- extractResults(funnorm_result_unbal_15, 
                                     data_name = 'funnorm_unbal_15')

# funnorm 20
funnorm_result_unbal_20 <- runModels(funnorm_cases, 
                                 bump_hunter = T, 
                                 bump_hunter_data = funnorm_unbal_bh_20)

funnorm_table_unbal_20 <- extractResults(funnorm_result_unbal_20, 
                                     data_name = 'funnorm_unbal_20')

# funnorm intersection
funnorm_result_int <- runModels(funnorm_cases, 
                            bump_hunter = T, 
                            bump_hunter_data = funnorm_int_feat)

funnorm_table_int <- extractResults(funnorm_result_int, 
                                data_name = 'funnorm_int')

funnorm_result_unbal_int <- runModels(funnorm_cases, 
                                  bump_hunter = T, 
                                  bump_hunter_data = funnorm_unbal_int_feat)

funnorm_table_unbal_int <- extractResults(funnorm_result_unbal_int, 
                                      data_name = 'funnorm_unbal_int')



##########
# threshold 0.07 - 0.15 sig
##########

#NORMAL
# funnorm 10 
funnorm_sig_result_10 <- runModels(funnorm_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = funnorm_bh_sig_10)

funnorm_sig_table_10 <- extractResults(funnorm_sig_result_10, 
                                   data_name = 'funnorm_sig_10')


# funnorm 15
funnorm_sig_result_15 <- runModels(funnorm_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = funnorm_bh_sig_15)

funnorm_sig_table_15 <- extractResults(funnorm_sig_result_15, 
                                   data_name = 'funnorm_sig_15')

# funnorm 20
funnorm_sig_result_20 <- runModels(funnorm_cases, 
                               bump_hunter = T, 
                               bump_hunter_data = funnorm_bh_sig_20)

funnorm_sig_table_20 <- extractResults(funnorm_sig_result_20, 
                                   data_name = 'funnorm_sig_20')



# funnorm intersection
funnorm_sig_result_int <- runModels(funnorm_cases, 
                                bump_hunter = T, 
                                bump_hunter_data = funnorm_int_sig_feat)

funnorm_sig_table_int <- extractResults(funnorm_sig_result_int, 
                                    data_name = 'funnorm_sig_int')



#UNBAL
# funnorm 10 
funnorm_sig_result_unbal_10 <- runModels(funnorm_cases, 
                                     bump_hunter = T, 
                                     bump_hunter_data = funnorm_unbal_bh_sig_10)

funnorm_sig_table_unbal_10 <- extractResults(funnorm_sig_result_unbal_10, 
                                         data_name = 'funnorm_unbal_sig_10')


# funnorm 15
funnorm_sig_result_unbal_15 <- runModels(funnorm_cases, 
                                     bump_hunter = T, 
                                     bump_hunter_data = funnorm_unbal_bh_sig_15)

funnorm_sig_table_unbal_15 <- extractResults(funnorm_sig_result_unbal_15, 
                                         data_name = 'funnorm_unbal_sig_15')

# # funnorm 20
funnorm_sig_result_unbal_20 <- runModels(funnorm_cases,
                                     bump_hunter = T,
                                     bump_hunter_data = funnorm_unbal_bh_sig_20)

funnorm_sig_table_unbal_20 <- extractResults(funnorm_sig_result_unbal_20,
                                         data_name = 'funnorm_unbal_sig_20')


# # # funnorm intersection
funnorm_sig_result_unbal_int <- runModels(funnorm_cases,
                                      bump_hunter = T,
                                      bump_hunter_data = funnorm_unbal_int_sig_feat)
#
funnorm_sig_table_unbal_int <- extractResults(funnorm_sig_result_unbal_int,
                                          data_name = 'funnorm_unbal_sig_int')


###########
# rbind tables and save RDA file
###########
funnorm_table <- rbind(funnorm_table_10, funnorm_table_15,
                   funnorm_table_20,
                   funnorm_sig_table_10, funnorm_sig_table_15,
                   funnorm_sig_table_20,
                   funnorm_table_int, funnorm_sig_table_int,
                   funnorm_table_unbal_10, funnorm_table_unbal_15,
                   funnorm_table_unbal_20, 
                   funnorm_sig_table_unbal_10, funnorm_sig_table_unbal_15,
                   funnorm_sig_table_unbal_20, 
                   funnorm_table_unbal_int, 
                   funnorm_sig_table_unbal_int)



# # remove data 
rm(funnorm_table_10, funnorm_table_15,
   funnorm_table_20,
   funnorm_sig_table_10, funnorm_sig_table_15,
   funnorm_sig_table_20,
   funnorm_table_int, funnorm_sig_table_int,
   funnorm_table_unbal_10, funnorm_table_unbal_15,
   funnorm_table_unbal_20, 
   funnorm_sig_table_unbal_10, funnorm_sig_table_unbal_15,
   funnorm_sig_table_unbal_20, 
   funnorm_table_unbal_int, 
   funnorm_sig_table_unbal_int)

#save table 
saveRDS(funnorm_table, 
        file = paste0(funnorm_folder, '/funnorm_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(funnorm_result_10, 
        file = paste0(funnorm_folder, '/funnorm_result_10.rda'))

saveRDS(funnorm_result_15, 
        file = paste0(funnorm_folder, '/funnorm_result_15.rda'))

saveRDS(funnorm_result_20, 
        file = paste0(funnorm_folder, '/funnorm_result_20.rda'))

saveRDS(funnorm_result_unbal_10, 
        file = paste0(funnorm_folder, '/funnorm_result_unbal_10.rda'))

saveRDS(funnorm_result_unbal_15, 
        file = paste0(funnorm_folder, '/funnorm_result_unbal_15.rda'))

saveRDS(funnorm_result_unbal_20, 
        file = paste0(funnorm_folder, '/funnorm_result_unbal_20.rda'))

