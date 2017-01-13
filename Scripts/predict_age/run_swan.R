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
source(paste0(scripts_folder, '/predict_age/model_functions.R'))
source(paste0(scripts_folder, '/predict_age/run_models.R'))

##########
# remove raw, quan, and funnorm objects
##########
rm(list = ls(pattern = "beta_raw_*"))
rm(list = ls(pattern = "beta_quan_*"))
rm(list = ls(pattern = "beta_funnorm_*"))


###################################################################################################################################
## beta_swan

##########
# cancer
##########

# bal cancer
swan_bal_cancer_models <- runModels(beta_swan, 
                                    bump_hunter = T, 
                                    bump_hunter_data = beta_swan_bal_cancer_features)

swan_bal_cancer_table <- extractResults(swan_bal_cancer_models, 
                                        data_name = 'beta_swan_bal_cancer')

# bal cancer sig
swan_bal_cancer_sig_models <- runModels(beta_swan, 
                                        bump_hunter = T, 
                                        bump_hunter_data = beta_swan_bal_cancer_sig_features)

swan_bal_cancer_sig_table <- extractResults(swan_bal_cancer_sig_models, 
                                            data_name = 'swan_bal_cancer_sig')


# bal counts cancer
swan_bal_counts_cancer_models <- runModels(beta_swan, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_swan_bal_counts_cancer_features)

swan_bal_counts_cancer_table <- extractResults(swan_bal_counts_cancer_models, 
                                              data_name = 'swan_bal_counts_cancer')

# bal counts cancer sig
swan_bal_counts_cancer_sig_models <- runModels(beta_swan, 
                                              bump_hunter = T, 
                                              bump_hunter_data =  beta_swan_bal_counts_cancer_sig_features)

swan_bal_counts_cancer_sig_table <- extractResults(swan_bal_counts_cancer_sig_models, 
                                                  data_name = 'swan_bal_counts_cancer_sig')

# unbal cancer
swan_unbal_cancer_models <- runModels(beta_swan, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_swan_unbal_cancer_features)

swan_unbal_cancer_table <- extractResults(swan_unbal_cancer_models, 
                                         data_name = 'swan_unbal_cancer')


# unbal cancer sig
swan_unbal_cancer_sig_models <- runModels(beta_swan, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_swan_unbal_cancer_sig_features)

swan_unbal_cancer_sig_table <- extractResults(swan_unbal_cancer_sig_models, 
                                             data_name = 'swan_unbal_cancer_sig')

##########
# p53
##########

# bal p53
swan_bal_p53_models <- runModels(beta_swan, 
                                bump_hunter = T, 
                                bump_hunter_data = beta_swan_bal_p53_features)

swan_bal_p53_table <- extractResults(swan_bal_p53_models, 
                                    data_name = 'swan_bal_p53')

# bal p53 sig
swan_bal_p53_sig_models <- runModels(beta_swan, 
                                    bump_hunter = T, 
                                    bump_hunter_data = beta_swan_bal_p53_sig_features)

swan_bal_p53_sig_table <- extractResults(swan_bal_p53_sig_models, 
                                        data_name = 'swan_bal_p53_sig')

# bal counts p53
swan_bal_counts_p53_models <- runModels(beta_swan, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_swan_bal_counts_p53_features)

swan_bal_counts_p53_table <- extractResults(swan_bal_counts_p53_models, 
                                           data_name = 'swan_bal_counts_p53')

# bal counts p53 sig
swan_bal_counts_p53_sig_models <- runModels(beta_swan, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_swan_bal_counts_p53_sig_features)

swan_bal_counts_p53_sig_table <- extractResults(swan_bal_counts_p53_sig_models, 
                                               data_name = 'swan_bal_counts_p53_sig')

# unbal p53
swan_unbal_p53_models <- runModels(beta_swan, 
                                  bump_hunter = T, 
                                  bump_hunter_data = beta_swan_unbal_p53_features)

swan_unbal_p53_table <- extractResults(swan_unbal_p53_models, 
                                      data_name = 'swan_unbal_p53')


# unbal p53 sig
swan_unbal_p53_sig_models <- runModels(beta_swan, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_swan_unbal_p53_sig_features)

swan_unbal_p53_sig_table <- extractResults(swan_unbal_p53_sig_models, 
                                          data_name = 'swan_unbal_p53_sig')


##########
# cancer
##########

# cancer_intersection
swan_cancer_intersection_models <- runModels(beta_swan, 
                                            bump_hunter = T, 
                                            bump_hunter_data = beta_swan_cancer_intersection_features)

swan_cancer_intersection_table <- extractResults(swan_cancer_intersection_models, 
                                                data_name = 'swan_cancer_intersection')

# cancer_intersection sig
swan_cancer_intersection_sig_models <- runModels(beta_swan, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_swan_cancer_intersection_sig_features)

swan_cancer_intersection_sig_table <- extractResults(swan_cancer_intersection_sig_models, 
                                                    data_name = 'swan_cancer_intersection_sig')

# cancer_union
swan_cancer_union_models <- runModels(beta_swan, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_swan_cancer_union_features)

swan_cancer_union_table <- extractResults(swan_cancer_union_models, 
                                         data_name = 'swan_cancer_union')

# cancer_union sig
swan_cancer_union_sig_models <- runModels(beta_swan, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_swan_cancer_union_sig_features)

swan_cancer_union_sig_table <- extractResults(swan_cancer_union_sig_models, 
                                             data_name = 'swan_cancer_union_sig')
##########
# p53
##########
# p53_intersection
swan_p53_intersection_models <- runModels(beta_swan, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_swan_p53_intersection_features)

swan_p53_intersection_table <- extractResults(swan_p53_intersection_models, 
                                             data_name = 'swan_p53_intersection')

# p53_intersection sig
swan_p53_intersection_sig_models <- runModels(beta_swan, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_swan_p53_intersection_sig_features)

swan_p53_intersection_sig_table <- extractResults(swan_p53_intersection_sig_models, 
                                                 data_name = 'swan_p53_intersection_sig')

# p53_union
swan_p53_union_models <- runModels(beta_swan, 
                                  bump_hunter = T, 
                                  bump_hunter_data = beta_swan_p53_union_features)

swan_p53_union_table <- extractResults(swan_p53_union_models, 
                                      data_name = 'swan_p53_union')

# p53_union sig
swan_p53_union_sig_models <- runModels(beta_swan, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_swan_p53_union_sig_features)

swan_p53_union_sig_table <- extractResults(swan_p53_union_sig_models, 
                                          data_name = 'swan_p53_union_sig')


##########
# most balanced cancer 
##########

# intersection
swan_bal_counts_cancer_intersection_models <- runModels(beta_swan, 
                                                       bump_hunter = T, 
                                                       bump_hunter_data = beta_bal_counts_cancer_intersection_features)

swan_bal_counts_cancer_intersection_table <- extractResults(swan_bal_counts_cancer_intersection_models, 
                                                           data_name = 'swan_bal_counts_cancer_intersection')


swan_bal_counts_cancer_intersection_sig_models <- runModels(beta_swan, 
                                                           bump_hunter = T, 
                                                           bump_hunter_data = beta_bal_counts_cancer_intersection_sig_features)

swan_bal_counts_cancer_intersection_sig_table <- extractResults(swan_bal_counts_cancer_intersection_sig_models, 
                                                               data_name = 'swan_bal_counts_cancer_intersection_sig')

# union
swan_bal_counts_cancer_union_models <- runModels(beta_swan, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_bal_counts_cancer_union_features)

swan_bal_counts_cancer_union_table <- extractResults(swan_bal_counts_cancer_union_models, 
                                                    data_name = 'swan_bal_counts_cancer_union')


swan_bal_counts_cancer_union_sig_models <- runModels(beta_swan, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_cancer_union_sig_features)

swan_bal_counts_cancer_union_sig_table <- extractResults(swan_bal_counts_cancer_union_sig_models, 
                                                        data_name = 'swan_bal_counts_cancer_union_sig')

##########
# most balanced p53 
##########

# intersection
swan_bal_counts_p53_intersection_models <- runModels(beta_swan, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_p53_intersection_features)

swan_bal_counts_p53_intersection_table <- extractResults(swan_bal_counts_p53_intersection_models, 
                                                        data_name = 'swan_bal_counts_p53_intersection')


swan_bal_counts_p53_intersection_sig_models <- runModels(beta_swan, 
                                                        bump_hunter = T, 
                                                        bump_hunter_data = beta_bal_counts_p53_intersection_sig_features)

swan_bal_counts_p53_intersection_sig_table <- extractResults(swan_bal_counts_p53_intersection_sig_models, 
                                                            data_name = 'swan_bal_counts_p53_intersection_sig')


# union
swan_bal_counts_p53_union_models <- runModels(beta_swan, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_bal_counts_p53_union_features)

swan_bal_counts_p53_union_table <- extractResults(swan_bal_counts_p53_union_models, 
                                                 data_name = 'swan_bal_counts_p53_union')


swan_bal_counts_p53_union_sig_models <- runModels(beta_swan, 
                                                 bump_hunter = T, 
                                                 bump_hunter_data = beta_bal_counts_p53_union_sig_features)

swan_bal_counts_p53_union_sig_table <- extractResults(swan_bal_counts_p53_union_sig_models, 
                                                     data_name = 'swan_bal_counts_p53_union_sig')

###########
# complete cancer - intersection across each method intersection
###########

# complete cancer intersection
swan_complete_cancer_intersection_models <- runModels(beta_swan, 
                                                     bump_hunter = T, 
                                                     bump_hunter_data = beta_cancer_intersection_features)

# get table 
swan_complete_cancer_intersection_table <- extractResults(swan_complete_cancer_intersection_models, 
                                                         data_name = 'swan_complete_cancer_intersection')

# empty features set
# # complete cancer intersection sig
# swan_complete_cancer_intersection_sig_models <- runModels(beta_swan, 
#                                                          bump_hunter = T, 
#                                                          bump_hunter_data = beta_cancer_intersection_sig_features)
# 
# # get table 
# swan_complete_cancer_intersection_sig_table <- extractResults(swan_complete_cancer_intersection_sig_models, 
#                                                              data_name = 'swan_complete_cancer_intersection_sig')


# complete cancer union
swan_complete_cancer_union_models <- runModels(beta_swan, 
                                              bump_hunter = T, 
                                              bump_hunter_data = beta_cancer_union_features)

# get table 
swan_complete_cancer_union_table <- extractResults(swan_complete_cancer_union_models, 
                                                  data_name = 'swan_complete_cancer_union')

# complete cancer union sig
swan_complete_cancer_union_sig_models <- runModels(beta_swan, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_cancer_union_sig_features)

# get table 
swan_complete_cancer_union_sig_table <- extractResults(swan_complete_cancer_union_sig_models, 
                                                      data_name = 'swan_complete_cancer_union_sig')


###########
# complete p53 - intersection across each method intersection
###########

# complete p53 intersection
swan_complete_p53_intersection_models <- runModels(beta_swan, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_p53_intersection_features)

# get table 
swan_complete_p53_intersection_table <- extractResults(swan_complete_p53_intersection_models, 
                                                      data_name = 'swan_complete_p53_intersection')

# complete p53 intersection sig
swan_complete_p53_intersection_sig_models <- runModels(beta_swan, 
                                                      bump_hunter = T, 
                                                      bump_hunter_data = beta_p53_intersection_sig_features)

# get table 
swan_complete_p53_intersection_sig_table <- extractResults(swan_complete_p53_intersection_sig_models, 
                                                          data_name = 'swan_complete_p53_intersection_sig')


# complete p53 union
swan_complete_p53_union_models <- runModels(beta_swan, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_p53_union_features)

# get table 
swan_complete_p53_union_table <- extractResults(swan_complete_p53_union_models, 
                                               data_name = 'swan_complete_p53_union')

# complete p53 union sig
swan_complete_p53_union_sig_models <- runModels(beta_swan, 
                                               bump_hunter = T, 
                                               bump_hunter_data = beta_p53_union_sig_features)

# get table 
swan_complete_p53_union_sig_table <- extractResults(swan_complete_p53_union_sig_models, 
                                                   data_name = 'swan_complete_p53_union_sig')


###########
# rbind tables and save RDA file
###########
swan_table <- rbind(swan_bal_cancer_table, swan_bal_cancer_sig_table, swan_bal_counts_cancer_table, swan_bal_counts_cancer_sig_table,
                   swan_unbal_cancer_table, swan_unbal_cancer_sig_table, swan_bal_p53_table, swan_bal_p53_sig_table, swan_bal_counts_p53_table, 
                   swan_bal_counts_p53_sig_table, swan_unbal_p53_table, swan_unbal_p53_sig_table, swan_cancer_intersection_table, 
                   swan_cancer_intersection_sig_table, swan_cancer_union_table, swan_cancer_union_sig_table, swan_p53_intersection_table, 
                   swan_p53_intersection_sig_table, swan_p53_union_table, swan_p53_union_sig_table,
                   swan_bal_counts_cancer_intersection_table, swan_bal_counts_cancer_intersection_sig_table,
                   swan_bal_counts_cancer_union_table, swan_bal_counts_cancer_union_sig_table,
                   swan_bal_counts_p53_intersection_table, swan_bal_counts_p53_intersection_sig_table,
                   swan_bal_counts_p53_union_table, swan_bal_counts_p53_union_sig_table,
                   swan_complete_cancer_intersection_table,
                   swan_complete_cancer_union_table, swan_complete_cancer_union_sig_table,
                   swan_complete_p53_intersection_table, swan_complete_p53_intersection_sig_table,
                   swan_complete_p53_union_table, swan_complete_p53_union_sig_table)

# remove data 
rm(swan_bal_cancer_table, swan_bal_cancer_sig_table, swan_bal_counts_cancer_table, swan_bal_counts_cancer_sig_table,
   swan_unbal_cancer_table, swan_unbal_cancer_sig_table, swan_bal_p53_table, swan_bal_p53_sig_table, swan_bal_counts_p53_table, 
   swan_bal_counts_p53_sig_table, swan_unbal_p53_table, swan_unbal_p53_sig_table, swan_cancer_intersection_table, 
   swan_cancer_intersection_sig_table, swan_cancer_union_table, swan_cancer_union_sig_table, swan_p53_intersection_table, 
   swan_p53_intersection_sig_table, swan_p53_union_table, swan_p53_union_sig_table,
   swan_bal_counts_cancer_intersection_table, swan_bal_counts_cancer_intersection_sig_table,
   swan_bal_counts_cancer_union_table, swan_bal_counts_cancer_union_sig_table,
   swan_bal_counts_p53_intersection_table, swan_bal_counts_p53_intersection_sig_table,
   swan_bal_counts_p53_union_table, swan_bal_counts_p53_union_sig_table,
   swan_complete_cancer_intersection_table, 
   swan_complete_cancer_union_table, swan_complete_cancer_union_sig_table,
   swan_complete_p53_intersection_table, swan_complete_p53_intersection_sig_table,
   swan_complete_p53_union_table, swan_complete_p53_union_sig_table)


#save table 
saveRDS(swan_table, 
        file = paste0(swan_folder, '/swan_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(swan_bal_cancer_models, 
        file = paste0(swan_folder, '/swan_bal_cancer_models.rda'))
saveRDS(swan_bal_cancer_sig_models, 
        file = paste0(swan_folder, '/swan_bal_cancer_sig_models.rda'))
saveRDS(swan_bal_counts_cancer_models, 
        file = paste0(swan_folder, '/swan_bal_counts_cancer_models.rda'))
saveRDS(swan_bal_counts_cancer_sig_models, 
        file = paste0(swan_folder, '/swan_bal_counts_cancer_sig_models.rda'))
saveRDS(swan_unbal_cancer_models, 
        file = paste0(swan_folder, '/swan_unbal_cancer_models.rda'))
saveRDS(swan_unbal_cancer_sig_models, 
        file = paste0(swan_folder, '/swan_unbal_cancer_sig_models.rda'))
saveRDS(swan_bal_p53_models, 
        file = paste0(swan_folder, '/swan_bal_p53_models.rda'))
saveRDS(swan_bal_p53_sig_models, 
        file = paste0(swan_folder, '/swan_bal_p53_sig_models.rda'))
saveRDS(swan_bal_counts_p53_models, 
        file = paste0(swan_folder, '/swan_bal_counts_p53_models.rda'))
saveRDS(swan_bal_counts_p53_sig_models, 
        file = paste0(swan_folder, '/swan_bal_counts_p53_sig_models.rda'))
saveRDS(swan_unbal_p53_models, 
        file = paste0(swan_folder, '/swan_unbal_p53_models.rda'))
saveRDS(swan_unbal_p53_sig_models, 
        file = paste0(swan_folder, '/swan_unbal_p53_sig_models.rda'))
saveRDS(swan_cancer_intersection_models, 
        file = paste0(swan_folder, '/swan_cancer_intersection_models.rda'))
saveRDS(swan_cancer_intersection_sig_models, 
        file = paste0(swan_folder, '/swan_cancer_intersection_sig_models.rda'))
saveRDS(swan_cancer_union_models, 
        file = paste0(swan_folder, '/swan_cancer_union_models.rda'))
saveRDS(swan_cancer_union_sig_models, 
        file = paste0(swan_folder, '/swan_cancer_union_sig_models.rda'))
saveRDS(swan_p53_intersection_models, 
        file = paste0(swan_folder, '/swan_p53_intersection_models.rda'))
saveRDS(swan_p53_intersection_sig_models, 
        file = paste0(swan_folder, '/swan_p53_intersection_sig_models.rda'))
saveRDS(swan_p53_union_models, 
        file = paste0(swan_folder, '/swan_p53_union_models.rda'))
saveRDS(swan_p53_union_sig_models, 
        file = paste0(swan_folder, '/swan_p53_union_sig_models.rda'))
saveRDS(swan_bal_counts_cancer_intersection_models,
        file = paste0(swan_folder, '/swan_bal_counts_cancer_intersection_models.rda'))
saveRDS(swan_bal_counts_cancer_intersection_sig_models,
        file = paste0(swan_folder, '/swan_bal_counts_cancer_intersection_sig_models.rda'))
saveRDS(swan_bal_counts_cancer_union_models,
        file = paste0(swan_folder, '/swan_bal_counts_cancer_union_models.rda'))
saveRDS(swan_bal_counts_cancer_union_sig_models,
        file = paste0(swan_folder, '/swan_bal_counts_cancer_union_sig_models.rda'))
saveRDS(swan_bal_counts_p53_intersection_models,
        file = paste0(swan_folder, '/swan_bal_counts_p53_intersection_models.rda'))
saveRDS(swan_bal_counts_p53_intersection_sig_models,
        file = paste0(swan_folder, '/swan_bal_counts_p53_intersection_sig_models.rda'))
saveRDS(swan_bal_counts_p53_union_models,
        file = paste0(swan_folder, '/swan_bal_counts_p53_union_models.rda'))
saveRDS(swan_bal_counts_p53_union_sig_models,
        file = paste0(swan_folder, '/swan_bal_counts_p53_union_sig_models.rda'))
saveRDS(swan_complete_cancer_intersection_models,
        file = paste0(swan_folder, '/swan_complete_cancer_intersection_models.rda'))
# saveRDS(swan_complete_cancer_intersection_sig_models,
#         file = paste0(swan_folder, '/swan_complete_cancer_intersection_sig_models.rda'))
saveRDS(swan_complete_cancer_union_models,
        file = paste0(swan_folder, '/swan_complete_cancer_union_models.rda'))
saveRDS(swan_complete_cancer_union_sig_models,
        file = paste0(swan_folder, '/swan_complete_cancer_union_sig_models.rda'))
saveRDS(swan_complete_p53_intersection_models,
        file = paste0(swan_folder, '/swan_complete_p53_intersection_models.rda'))
saveRDS(swan_complete_p53_intersection_sig_models,
        file = paste0(swan_folder, '/swan_complete_p53_intersection_sig_models.rda'))
saveRDS(swan_complete_p53_union_models,
        file = paste0(swan_folder, '/swan_complete_p53_union_models.rda'))
saveRDS(swan_complete_p53_union_sig_models,
        file = paste0(swan_folder, '/swan_complete_p53_union_sig_models.rda'))



