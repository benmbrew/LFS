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
source(paste0(scripts_folder, '/predict_age/model_functions.R'))
source(paste0(scripts_folder, '/predict_age/run_models.R'))

##########
# remove raw, swan, and funnorm objects
##########
rm(list = ls(pattern = "beta_raw_*"))
rm(list = ls(pattern = "beta_swan_*"))
rm(list = ls(pattern = "beta_funnorm_*"))

###################################################################################################################################
## beta_quan

##########
# cancer
##########

# bal cancer
quan_bal_cancer_models <- runModels(beta_quan, 
                                    bump_hunter = T, 
                                    bump_hunter_data = beta_quan_bal_cancer_features)

quan_bal_cancer_table <- extractResults(quan_bal_cancer_models, 
                                        data_name = 'beta_quan_bal_cancer')

# bal cancer sig
quan_bal_cancer_sig_models <- runModels(beta_quan, 
                                        bump_hunter = T, 
                                        bump_hunter_data = beta_quan_bal_cancer_sig_features)

quan_bal_cancer_sig_table <- extractResults(quan_bal_cancer_sig_models, 
                                            data_name = 'quan_bal_cancer_sig')


# bal counts cancer
quan_bal_counts_cancer_models <- runModels(beta_quan, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_quan_bal_counts_cancer_features)

quan_bal_counts_cancer_table <- extractResults(quan_bal_counts_cancer_models, 
                                               data_name = 'quan_bal_counts_cancer')

# bal counts cancer sig
quan_bal_counts_cancer_sig_models <- runModels(beta_quan, 
                                               bump_hunter = T, 
                                               bump_hunter_data =  beta_quan_bal_counts_cancer_sig_features)

quan_bal_counts_cancer_sig_table <- extractResults(quan_bal_counts_cancer_sig_models, 
                                                   data_name = 'quan_bal_counts_cancer_sig')

# unbal cancer
quan_unbal_cancer_models <- runModels(beta_quan, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_quan_unbal_cancer_features)

quan_unbal_cancer_table <- extractResults(quan_unbal_cancer_models, 
                                          data_name = 'quan_unbal_cancer')


# unbal cancer sig
quan_unbal_cancer_sig_models <- runModels(beta_quan, 
                                          bump_hunter = T, 
                                          bump_hunter_data = beta_quan_unbal_cancer_sig_features)

quan_unbal_cancer_sig_table <- extractResults(quan_unbal_cancer_sig_models, 
                                              data_name = 'quan_unbal_cancer_sig')

##########
# p53
##########

# bal p53
quan_bal_p53_models <- runModels(beta_quan, 
                                 bump_hunter = T, 
                                 bump_hunter_data = beta_quan_bal_p53_features)

quan_bal_p53_table <- extractResults(quan_bal_p53_models, 
                                     data_name = 'quan_bal_p53')

# bal p53 sig
quan_bal_p53_sig_models <- runModels(beta_quan, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_quan_bal_p53_sig_features)

quan_bal_p53_sig_table <- extractResults(quan_bal_p53_sig_models, 
                                         data_name = 'quan_bal_p53_sig')

# bal counts p53
quan_bal_counts_p53_models <- runModels(beta_quan, 
                                        bump_hunter = T, 
                                        bump_hunter_data = beta_quan_bal_counts_p53_features)

quan_bal_counts_p53_table <- extractResults(quan_bal_counts_p53_models, 
                                            data_name = 'quan_bal_counts_p53')

# bal counts p53 sig
quan_bal_counts_p53_sig_models <- runModels(beta_quan, 
                                            bump_hunter = T, 
                                            bump_hunter_data = beta_quan_bal_counts_p53_sig_features)

quan_bal_counts_p53_sig_table <- extractResults(quan_bal_counts_p53_sig_models, 
                                                data_name = 'quan_bal_counts_p53_sig')

# unbal p53
quan_unbal_p53_models <- runModels(beta_quan, 
                                   bump_hunter = T, 
                                   bump_hunter_data = beta_quan_unbal_p53_features)

quan_unbal_p53_table <- extractResults(quan_unbal_p53_models, 
                                       data_name = 'quan_unbal_p53')


# unbal p53 sig
quan_unbal_p53_sig_models <- runModels(beta_quan, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_quan_unbal_p53_sig_features)

quan_unbal_p53_sig_table <- extractResults(quan_unbal_p53_sig_models, 
                                           data_name = 'quan_unbal_p53_sig')


##########
# cancer
##########

# cancer_intersection
quan_cancer_intersection_models <- runModels(beta_quan, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_quan_cancer_intersection_features)

quan_cancer_intersection_table <- extractResults(quan_cancer_intersection_models, 
                                                 data_name = 'quan_cancer_intersection')

# cancer_intersection sig
quan_cancer_intersection_sig_models <- runModels(beta_quan, 
                                                 bump_hunter = T, 
                                                 bump_hunter_data = beta_quan_cancer_intersection_sig_features)

quan_cancer_intersection_sig_table <- extractResults(quan_cancer_intersection_sig_models, 
                                                     data_name = 'quan_cancer_intersection_sig')

# cancer_union
quan_cancer_union_models <- runModels(beta_quan, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_quan_cancer_union_features)

quan_cancer_union_table <- extractResults(quan_cancer_union_models, 
                                          data_name = 'quan_cancer_union')

# cancer_union sig
quan_cancer_union_sig_models <- runModels(beta_quan, 
                                          bump_hunter = T, 
                                          bump_hunter_data = beta_quan_cancer_union_sig_features)

quan_cancer_union_sig_table <- extractResults(quan_cancer_union_sig_models, 
                                              data_name = 'quan_cancer_union_sig')
##########
# p53
##########
# p53_intersection
quan_p53_intersection_models <- runModels(beta_quan, 
                                          bump_hunter = T, 
                                          bump_hunter_data = beta_quan_p53_intersection_features)

quan_p53_intersection_table <- extractResults(quan_p53_intersection_models, 
                                             data_name = 'quan_p53_intersection')

# p53_intersection sig
quan_p53_intersection_sig_models <- runModels(beta_quan, 
                                              bump_hunter = T, 
                                              bump_hunter_data = beta_quan_p53_intersection_sig_features)

quan_p53_intersection_sig_table <- extractResults(quan_p53_intersection_sig_models, 
                                                  data_name = 'quan_p53_intersection_sig')

# p53_union
quan_p53_union_models <- runModels(beta_quan, 
                                   bump_hunter = T, 
                                   bump_hunter_data = beta_quan_p53_union_features)

quan_p53_union_table <- extractResults(quan_p53_union_models, 
                                       data_name = 'quan_p53_union')

# p53_union sig
quan_p53_union_sig_models <- runModels(beta_quan, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_quan_p53_union_sig_features)
 
quan_p53_union_sig_table <- extractResults(quan_p53_union_sig_models, 
                                           data_name = 'quan_p53_union_sig')


##########
# most balanced cancer 
##########

# intersection
quan_bal_counts_cancer_intersection_models <- runModels(beta_quan, 
                                                        bump_hunter = T, 
                                                        bump_hunter_data = beta_bal_counts_cancer_intersection_features)

quan_bal_counts_cancer_intersection_table <- extractResults(quan_bal_counts_cancer_intersection_models, 
                                                            data_name = 'quan_bal_counts_cancer_intersection')


quan_bal_counts_cancer_intersection_sig_models <- runModels(beta_quan, 
                                                            bump_hunter = T, 
                                                            bump_hunter_data = beta_bal_counts_cancer_intersection_sig_features)

quan_bal_counts_cancer_intersection_sig_table <- extractResults(quan_bal_counts_cancer_intersection_sig_models, 
                                                                data_name = 'quan_bal_counts_cancer_intersection_sig')

# union
quan_bal_counts_cancer_union_models <- runModels(beta_quan, 
                                                 bump_hunter = T, 
                                                 bump_hunter_data = beta_bal_counts_cancer_union_features)

quan_bal_counts_cancer_union_table <- extractResults(quan_bal_counts_cancer_union_models, 
                                                     data_name = 'quan_bal_counts_cancer_union')


quan_bal_counts_cancer_union_sig_models <- runModels(beta_quan, 
                                                     bump_hunter = T, 
                                                     bump_hunter_data = beta_bal_counts_cancer_union_sig_features)

quan_bal_counts_cancer_union_sig_table <- extractResults(quan_bal_counts_cancer_union_sig_models, 
                                                         data_name = 'quan_bal_counts_cancer_union_sig')

##########
# most balanced p53 
##########

# intersection
quan_bal_counts_p53_intersection_models <- runModels(beta_quan, 
                                                     bump_hunter = T, 
                                                     bump_hunter_data = beta_bal_counts_p53_intersection_features)

quan_bal_counts_p53_intersection_table <- extractResults(quan_bal_counts_p53_intersection_models, 
                                                         data_name = 'quan_bal_counts_p53_intersection')


quan_bal_counts_p53_intersection_sig_models <- runModels(beta_quan, 
                                                         bump_hunter = T, 
                                                         bump_hunter_data = beta_bal_counts_p53_intersection_sig_features)

quan_bal_counts_p53_intersection_sig_table <- extractResults(quan_bal_counts_p53_intersection_sig_models, 
                                                             data_name = 'quan_bal_counts_p53_intersection_sig')


# union
quan_bal_counts_p53_union_models <- runModels(beta_quan, 
                                              bump_hunter = T, 
                                              bump_hunter_data = beta_bal_counts_p53_union_features)

quan_bal_counts_p53_union_table <- extractResults(quan_bal_counts_p53_union_models, 
                                                  data_name = 'quan_bal_counts_p53_union')


quan_bal_counts_p53_union_sig_models <- runModels(beta_quan, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_bal_counts_p53_union_sig_features)

quan_bal_counts_p53_union_sig_table <- extractResults(quan_bal_counts_p53_union_sig_models, 
                                                      data_name = 'quan_bal_counts_p53_union_sig')

###########
# complete cancer - intersection across each method intersection
###########

# complete cancer intersection
quan_complete_cancer_intersection_models <- runModels(beta_quan, 
                                                      bump_hunter = T, 
                                                      bump_hunter_data = beta_cancer_intersection_features)

# get table 
quan_complete_cancer_intersection_table <- extractResults(quan_complete_cancer_intersection_models, 
                                                          data_name = 'quan_complete_cancer_intersection')

# empty features set
# # complete cancer intersection sig
# quan_complete_cancer_intersection_sig_models <- runModels(beta_quan, 
#                                                           bump_hunter = T, 
#                                                           bump_hunter_data = beta_cancer_intersection_sig_features)
# 
# # get table 
# quan_complete_cancer_intersection_sig_table <- extractResults(quan_complete_cancer_intersection_sig_models, 
#                                                               data_name = 'quan_complete_cancer_intersection_sig')


# complete cancer union
quan_complete_cancer_union_models <- runModels(beta_quan, 
                                               bump_hunter = T, 
                                               bump_hunter_data = beta_cancer_union_features)

# get table 
quan_complete_cancer_union_table <- extractResults(quan_complete_cancer_union_models, 
                                                   data_name = 'quan_complete_cancer_union')

# complete cancer union sig
quan_complete_cancer_union_sig_models <- runModels(beta_quan, 
                                                   bump_hunter = T, 
                                                   bump_hunter_data = beta_cancer_union_sig_features)

# get table 
quan_complete_cancer_union_sig_table <- extractResults(quan_complete_cancer_union_sig_models, 
                                                       data_name = 'quan_complete_cancer_union_sig')


###########
# complete p53 - intersection across each method intersection
###########

# complete p53 intersection
quan_complete_p53_intersection_models <- runModels(beta_quan, 
                                                   bump_hunter = T, 
                                                   bump_hunter_data = beta_p53_intersection_features)

# get table 
quan_complete_p53_intersection_table <- extractResults(quan_complete_p53_intersection_models, 
                                                       data_name = 'quan_complete_p53_intersection')

# complete p53 intersection sig
quan_complete_p53_intersection_sig_models <- runModels(beta_quan, 
                                                       bump_hunter = T, 
                                                       bump_hunter_data = beta_p53_intersection_sig_features)

# get table 
quan_complete_p53_intersection_sig_table <- extractResults(quan_complete_p53_intersection_sig_models, 
                                                           data_name = 'quan_complete_p53_intersection_sig')


# complete p53 union
quan_complete_p53_union_models <- runModels(beta_quan, 
                                            bump_hunter = T, 
                                            bump_hunter_data = beta_p53_union_features)

# get table 
quan_complete_p53_union_table <- extractResults(quan_complete_p53_union_models, 
                                                data_name = 'quan_complete_p53_union')

# complete p53 union sig
quan_complete_p53_union_sig_models <- runModels(beta_quan, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_p53_union_sig_features)

# get table 
quan_complete_p53_union_sig_table <- extractResults(quan_complete_p53_union_sig_models, 
                                                    data_name = 'quan_complete_p53_union_sig')


###########
# rbind tables and save RDA file
###########
quan_table <- rbind(quan_bal_cancer_table, quan_bal_cancer_sig_table, quan_bal_counts_cancer_table, quan_bal_counts_cancer_sig_table,
                   quan_unbal_cancer_table, quan_unbal_cancer_sig_table, quan_bal_p53_table, quan_bal_p53_sig_table, quan_bal_counts_p53_table, 
                   quan_bal_counts_p53_sig_table, quan_unbal_p53_table, quan_unbal_p53_sig_table, quan_cancer_intersection_table, 
                   quan_cancer_intersection_sig_table, quan_cancer_union_table, quan_cancer_union_sig_table, quan_p53_intersection_table, 
                   quan_p53_intersection_sig_table, quan_p53_union_table, quan_p53_union_sig_table,
                   quan_bal_counts_cancer_intersection_table, quan_bal_counts_cancer_intersection_sig_table,
                   quan_bal_counts_cancer_union_table, quan_bal_counts_cancer_union_sig_table,
                   quan_bal_counts_p53_intersection_table, quan_bal_counts_p53_intersection_sig_table,
                   quan_bal_counts_p53_union_table, quan_bal_counts_p53_union_sig_table,
                   quan_complete_cancer_intersection_table, #quan_complete_cancer_intersection_sig_table,
                   quan_complete_cancer_union_table, quan_complete_cancer_union_sig_table,
                   quan_complete_p53_intersection_table, quan_complete_p53_intersection_sig_table,
                   quan_complete_p53_union_table, quan_complete_p53_union_sig_table)

# remove data 
rm(quan_bal_cancer_table, quan_bal_cancer_sig_table, quan_bal_counts_cancer_table, quan_bal_counts_cancer_sig_table,
   quan_unbal_cancer_table, quan_unbal_cancer_sig_table, quan_bal_p53_table, quan_bal_p53_sig_table, quan_bal_counts_p53_table, 
   quan_bal_counts_p53_sig_table, quan_unbal_p53_table, quan_unbal_p53_sig_table, quan_cancer_intersection_table, 
   quan_cancer_intersection_sig_table, quan_cancer_union_table, quan_cancer_union_sig_table, quan_p53_intersection_table, 
   quan_p53_intersection_sig_table, quan_p53_union_table, quan_p53_union_sig_table,
   quan_bal_counts_cancer_intersection_table, quan_bal_counts_cancer_intersection_sig_table,
   quan_bal_counts_cancer_union_table, quan_bal_counts_cancer_union_sig_table,
   quan_bal_counts_p53_intersection_table, quan_bal_counts_p53_intersection_sig_table,
   quan_bal_counts_p53_union_table, quan_bal_counts_p53_union_sig_table,
   quan_complete_cancer_intersection_table, #quan_complete_cancer_intersection_sig_table,
   quan_complete_cancer_union_table, quan_complete_cancer_union_sig_table,
   quan_complete_p53_intersection_table, quan_complete_p53_intersection_sig_table,
   quan_complete_p53_union_table, quan_complete_p53_union_sig_table)


#save table 
saveRDS(quan_table, 
        file = paste0(quan_folder, '/quan_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(quan_bal_cancer_models, 
        file = paste0(quan_folder, '/quan_bal_cancer_models.rda'))
saveRDS(quan_bal_cancer_sig_models, 
        file = paste0(quan_folder, '/quan_bal_cancer_sig_models.rda'))
saveRDS(quan_bal_counts_cancer_models, 
        file = paste0(quan_folder, '/quan_bal_counts_cancer_models.rda'))
saveRDS(quan_bal_counts_cancer_sig_models, 
        file = paste0(quan_folder, '/quan_bal_counts_cancer_sig_models.rda'))
saveRDS(quan_unbal_cancer_models, 
        file = paste0(quan_folder, '/quan_unbal_cancer_models.rda'))
saveRDS(quan_unbal_cancer_sig_models, 
        file = paste0(quan_folder, '/quan_unbal_cancer_sig_models.rda'))
saveRDS(quan_bal_p53_models, 
        file = paste0(quan_folder, '/quan_bal_p53_models.rda'))
saveRDS(quan_bal_p53_sig_models, 
        file = paste0(quan_folder, '/quan_bal_p53_sig_models.rda'))
saveRDS(quan_bal_counts_p53_models, 
        file = paste0(quan_folder, '/quan_bal_counts_p53_models.rda'))
saveRDS(quan_bal_counts_p53_sig_models, 
        file = paste0(quan_folder, '/quan_bal_counts_p53_sig_models.rda'))
saveRDS(quan_unbal_p53_models, 
        file = paste0(quan_folder, '/quan_unbal_p53_models.rda'))
saveRDS(quan_unbal_p53_sig_models, 
        file = paste0(quan_folder, '/quan_unbal_p53_sig_models.rda'))
saveRDS(quan_cancer_intersection_models, 
        file = paste0(quan_folder, '/quan_cancer_intersection_models.rda'))
saveRDS(quan_cancer_intersection_sig_models, 
        file = paste0(quan_folder, '/quan_cancer_intersection_sig_models.rda'))
saveRDS(quan_cancer_union_models, 
        file = paste0(quan_folder, '/quan_cancer_union_models.rda'))
saveRDS(quan_cancer_union_sig_models, 
        file = paste0(quan_folder, '/quan_cancer_union_sig_models.rda'))
saveRDS(quan_p53_intersection_models, 
        file = paste0(quan_folder, '/quan_p53_intersection_models.rda'))
saveRDS(quan_p53_intersection_sig_models, 
        file = paste0(quan_folder, '/quan_p53_intersection_sig_models.rda'))
saveRDS(quan_p53_union_models, 
        file = paste0(quan_folder, '/quan_p53_union_models.rda'))
saveRDS(quan_p53_union_sig_models, 
        file = paste0(quan_folder, '/quan_p53_union_sig_models.rda'))
saveRDS(quan_bal_counts_cancer_intersection_models,
        file = paste0(quan_folder, '/quan_bal_counts_cancer_intersection_models.rda'))
saveRDS(quan_bal_counts_cancer_intersection_sig_models,
        file = paste0(quan_folder, '/quan_bal_counts_cancer_intersection_sig_models.rda'))
saveRDS(quan_bal_counts_cancer_union_models,
        file = paste0(quan_folder, '/quan_bal_counts_cancer_union_models.rda'))
saveRDS(quan_bal_counts_cancer_union_sig_models,
        file = paste0(quan_folder, '/quan_bal_counts_cancer_union_sig_models.rda'))
saveRDS(quan_bal_counts_p53_intersection_models,
        file = paste0(quan_folder, '/quan_bal_counts_p53_intersection_models.rda'))
saveRDS(quan_bal_counts_p53_intersection_sig_models,
        file = paste0(quan_folder, '/quan_bal_counts_p53_intersection_sig_models.rda'))
saveRDS(quan_bal_counts_p53_union_models,
        file = paste0(quan_folder, '/quan_bal_counts_p53_union_models.rda'))
saveRDS(quan_bal_counts_p53_union_sig_models,
        file = paste0(quan_folder, '/quan_bal_counts_p53_union_sig_models.rda'))
saveRDS(quan_complete_cancer_intersection_models,
        file = paste0(quan_folder, '/quan_complete_cancer_intersection_models.rda'))
# saveRDS(quan_complete_cancer_intersection_sig_models,
#         file = paste0(quan_folder, '/quan_complete_cancer_intersection_sig_models.rda'))
saveRDS(quan_complete_cancer_union_models,
        file = paste0(quan_folder, '/quan_complete_cancer_union_models.rda'))
saveRDS(quan_complete_cancer_union_sig_models,
        file = paste0(quan_folder, '/quan_complete_cancer_union_sig_models.rda'))
saveRDS(quan_complete_p53_intersection_models,
        file = paste0(quan_folder, '/quan_complete_p53_intersection_models.rda'))
saveRDS(quan_complete_p53_intersection_sig_models,
        file = paste0(quan_folder, '/quan_complete_p53_intersection_sig_models.rda'))
saveRDS(quan_complete_p53_union_models,
        file = paste0(quan_folder, '/quan_complete_p53_union_models.rda'))
saveRDS(quan_complete_p53_union_sig_models,
        file = paste0(quan_folder, '/quan_complete_p53_union_sig_models.rda'))



