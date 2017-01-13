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
source(paste0(scripts_folder, '/predict_age/model_functions.R'))
source(paste0(scripts_folder, '/predict_age/run_models.R'))

##########
# remove quan, swan, and funnorm objects
##########
rm(list = ls(pattern = "beta_quan_*"))
rm(list = ls(pattern = "beta_swan_*"))
rm(list = ls(pattern = "beta_funnorm_*"))


###################################################################################################################################
## beta_raw

##########
# cancer
##########

# bal cancer
raw_bal_cancer_models <- runModels(beta_raw, 
                                   bump_hunter = T, 
                                   bump_hunter_data = beta_raw_bal_cancer_features)

raw_bal_cancer_table <- extractResults(raw_bal_cancer_models, 
                                       data_name = 'beta_raw_bal_cancer')

# bal cancer sig
raw_bal_cancer_sig_models <- runModels(beta_raw, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_raw_bal_cancer_sig_features)

raw_bal_cancer_sig_table <- extractResults(raw_bal_cancer_sig_models, 
                                           data_name = 'raw_bal_cancer_sig')


# bal counts cancer
raw_bal_counts_cancer_models <- runModels(beta_raw, 
                                          bump_hunter = T, 
                                          bump_hunter_data = beta_raw_bal_counts_cancer_features)

raw_bal_counts_cancer_table <- extractResults(raw_bal_counts_cancer_models, 
                                              data_name = 'raw_bal_counts_cancer')

# bal counts cancer sig
raw_bal_counts_cancer_sig_models <- runModels(beta_raw, 
                                              bump_hunter = T, 
                                              bump_hunter_data =  beta_raw_bal_counts_cancer_sig_features)

raw_bal_counts_cancer_sig_table <- extractResults(raw_bal_counts_cancer_sig_models, 
                                                  data_name = 'raw_bal_counts_cancer_sig')

# unbal cancer
raw_unbal_cancer_models <- runModels(beta_raw, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_raw_unbal_cancer_features)

raw_unbal_cancer_table <- extractResults(raw_unbal_cancer_models, 
                                         data_name = 'raw_unbal_cancer')


# unbal cancer sig
raw_unbal_cancer_sig_models <- runModels(beta_raw, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_raw_unbal_cancer_sig_features)

raw_unbal_cancer_sig_table <- extractResults(raw_unbal_cancer_sig_models, 
                                             data_name = 'raw_unbal_cancer_sig')

##########
# p53
##########

# bal p53
raw_bal_p53_models <- runModels(beta_raw, 
                                bump_hunter = T, 
                                bump_hunter_data = beta_raw_bal_p53_features)

raw_bal_p53_table <- extractResults(raw_bal_p53_models, 
                                    data_name = 'raw_bal_p53')

# bal p53 sig
raw_bal_p53_sig_models <- runModels(beta_raw, 
                                    bump_hunter = T, 
                                    bump_hunter_data = beta_raw_bal_p53_sig_features)

raw_bal_p53_sig_table <- extractResults(raw_bal_p53_sig_models, 
                                        data_name = 'raw_bal_p53_sig')

# bal counts p53
raw_bal_counts_p53_models <- runModels(beta_raw, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_raw_bal_counts_p53_features)

raw_bal_counts_p53_table <- extractResults(raw_bal_counts_p53_models, 
                                           data_name = 'raw_bal_counts_p53')

# bal counts p53 sig
raw_bal_counts_p53_sig_models <- runModels(beta_raw, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_raw_bal_counts_p53_sig_features)

raw_bal_counts_p53_sig_table <- extractResults(raw_bal_counts_p53_sig_models, 
                                               data_name = 'raw_bal_counts_p53_sig')

# unbal p53
raw_unbal_p53_models <- runModels(beta_raw, 
                                  bump_hunter = T, 
                                  bump_hunter_data = beta_raw_unbal_p53_features)

raw_unbal_p53_table <- extractResults(raw_unbal_p53_models, 
                                      data_name = 'raw_unbal_p53')


# unbal p53 sig
raw_unbal_p53_sig_models <- runModels(beta_raw, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_raw_unbal_p53_sig_features)

raw_unbal_p53_sig_table <- extractResults(raw_unbal_p53_sig_models, 
                                          data_name = 'raw_unbal_p53_sig')


##########
# cancer
##########

# cancer_intersection
raw_cancer_intersection_models <- runModels(beta_raw, 
                                            bump_hunter = T, 
                                            bump_hunter_data = beta_raw_cancer_intersection_features)

raw_cancer_intersection_table <- extractResults(raw_cancer_intersection_models, 
                                                data_name = 'raw_cancer_intersection')

# cancer_intersection sig
raw_cancer_intersection_sig_models <- runModels(beta_raw, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_raw_cancer_intersection_sig_features)

raw_cancer_intersection_sig_table <- extractResults(raw_cancer_intersection_sig_models, 
                                                    data_name = 'raw_cancer_intersection_sig')

# cancer_union
raw_cancer_union_models <- runModels(beta_raw, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_raw_cancer_union_features)

raw_cancer_union_table <- extractResults(raw_cancer_union_models, 
                                         data_name = 'raw_cancer_union')

# cancer_union sig
raw_cancer_union_sig_models <- runModels(beta_raw, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_raw_cancer_union_sig_features)

raw_cancer_union_sig_table <- extractResults(raw_cancer_union_sig_models, 
                                             data_name = 'raw_cancer_union_sig')
##########
# p53
##########
# p53_intersection
raw_p53_intersection_models <- runModels(beta_raw, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_raw_p53_intersection_features)

raw_p53_intersection_table <- extractResults(raw_p53_intersection_models, 
                                             data_name = 'raw_p53_intersection')

# p53_intersection sig
raw_p53_intersection_sig_models <- runModels(beta_raw, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_raw_p53_intersection_sig_features)

raw_p53_intersection_sig_table <- extractResults(raw_p53_intersection_sig_models, 
                                                 data_name = 'raw_p53_intersection_sig')

# p53_union
raw_p53_union_models <- runModels(beta_raw, 
                                  bump_hunter = T, 
                                  bump_hunter_data = beta_raw_p53_union_features)

raw_p53_union_table <- extractResults(raw_p53_union_models, 
                                      data_name = 'raw_p53_union')

# p53_union sig
raw_p53_union_sig_models <- runModels(beta_raw, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_raw_p53_union_sig_features)

raw_p53_union_sig_table <- extractResults(raw_p53_union_sig_models, 
                                          data_name = 'raw_p53_union_sig')


##########
# most balanced cancer 
##########

# intersection
raw_bal_counts_cancer_intersection_models <- runModels(beta_raw, 
                                                       bump_hunter = T, 
                                                       bump_hunter_data = beta_bal_counts_cancer_intersection_features)

raw_bal_counts_cancer_intersection_table <- extractResults(raw_bal_counts_cancer_intersection_models, 
                                                           data_name = 'raw_bal_counts_cancer_intersection')


raw_bal_counts_cancer_intersection_sig_models <- runModels(beta_raw, 
                                                           bump_hunter = T, 
                                                           bump_hunter_data = beta_bal_counts_cancer_intersection_sig_features)

raw_bal_counts_cancer_intersection_sig_table <- extractResults(raw_bal_counts_cancer_intersection_sig_models, 
                                                               data_name = 'raw_bal_counts_cancer_intersection_sig')

# union
raw_bal_counts_cancer_union_models <- runModels(beta_raw, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_bal_counts_cancer_union_features)

raw_bal_counts_cancer_union_table <- extractResults(raw_bal_counts_cancer_union_models, 
                                                    data_name = 'raw_bal_counts_cancer_union')


raw_bal_counts_cancer_union_sig_models <- runModels(beta_raw, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_cancer_union_sig_features)

raw_bal_counts_cancer_union_sig_table <- extractResults(raw_bal_counts_cancer_union_sig_models, 
                                                        data_name = 'raw_bal_counts_cancer_union_sig')

##########
# most balanced p53 
##########

# intersection
raw_bal_counts_p53_intersection_models <- runModels(beta_raw, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_p53_intersection_features)

raw_bal_counts_p53_intersection_table <- extractResults(raw_bal_counts_p53_intersection_models, 
                                                        data_name = 'raw_bal_counts_p53_intersection')


raw_bal_counts_p53_intersection_sig_models <- runModels(beta_raw, 
                                                        bump_hunter = T, 
                                                        bump_hunter_data = beta_bal_counts_p53_intersection_sig_features)

raw_bal_counts_p53_intersection_sig_table <- extractResults(raw_bal_counts_p53_intersection_sig_models, 
                                                            data_name = 'raw_bal_counts_p53_intersection_sig')


# union
raw_bal_counts_p53_union_models <- runModels(beta_raw, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_bal_counts_p53_union_features)

raw_bal_counts_p53_union_table <- extractResults(raw_bal_counts_p53_union_models, 
                                                 data_name = 'raw_bal_counts_p53_union')


raw_bal_counts_p53_union_sig_models <- runModels(beta_raw, 
                                                 bump_hunter = T, 
                                                 bump_hunter_data = beta_bal_counts_p53_union_sig_features)

raw_bal_counts_p53_union_sig_table <- extractResults(raw_bal_counts_p53_union_sig_models, 
                                                     data_name = 'raw_bal_counts_p53_union_sig')

###########
# complete cancer - intersection across each method intersection
###########

# complete cancer intersection
raw_complete_cancer_intersection_models <- runModels(beta_raw, 
                                                     bump_hunter = T, 
                                                     bump_hunter_data = beta_cancer_intersection_features)

# get table 
raw_complete_cancer_intersection_table <- extractResults(raw_complete_cancer_intersection_models, 
                                                         data_name = 'raw_complete_cancer_intersection')

# empty features set
# # complete cancer intersection sig
# raw_complete_cancer_intersection_sig_models <- runModels(beta_raw, 
#                                                          bump_hunter = T, 
#                                                          bump_hunter_data = beta_cancer_intersection_sig_features)
# 
# # get table 
# raw_complete_cancer_intersection_sig_table <- extractResults(raw_complete_cancer_intersection_sig_models, 
#                                                              data_name = 'raw_complete_cancer_intersection_sig')


# complete cancer union
raw_complete_cancer_union_models <- runModels(beta_raw, 
                                              bump_hunter = T, 
                                              bump_hunter_data = beta_cancer_union_features)

# get table 
raw_complete_cancer_union_table <- extractResults(raw_complete_cancer_union_models, 
                                                  data_name = 'raw_complete_cancer_union')

# complete cancer union sig
raw_complete_cancer_union_sig_models <- runModels(beta_raw, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_cancer_union_sig_features)

# get table 
raw_complete_cancer_union_sig_table <- extractResults(raw_complete_cancer_union_sig_models, 
                                                      data_name = 'raw_complete_cancer_union_sig')


###########
# complete p53 - intersection across each method intersection
###########

# complete p53 intersection
raw_complete_p53_intersection_models <- runModels(beta_raw, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_p53_intersection_features)

# get table 
raw_complete_p53_intersection_table <- extractResults(raw_complete_p53_intersection_models, 
                                                      data_name = 'raw_complete_p53_intersection')

# complete p53 intersection sig
raw_complete_p53_intersection_sig_models <- runModels(beta_raw, 
                                                      bump_hunter = T, 
                                                      bump_hunter_data = beta_p53_intersection_sig_features)

# get table 
raw_complete_p53_intersection_sig_table <- extractResults(raw_complete_p53_intersection_sig_models, 
                                                          data_name = 'raw_complete_p53_intersection_sig')


# complete p53 union
raw_complete_p53_union_models <- runModels(beta_raw, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_p53_union_features)

# get table 
raw_complete_p53_union_table <- extractResults(raw_complete_p53_union_models, 
                                               data_name = 'raw_complete_p53_union')

# complete p53 union sig
raw_complete_p53_union_sig_models <- runModels(beta_raw, 
                                               bump_hunter = T, 
                                               bump_hunter_data = beta_p53_union_sig_features)

# get table 
raw_complete_p53_union_sig_table <- extractResults(raw_complete_p53_union_sig_models, 
                                                   data_name = 'raw_complete_p53_union_sig')


###########
# rbind tables and save RDA file
###########
raw_table <- rbind(raw_bal_cancer_table, raw_bal_cancer_sig_table, raw_bal_counts_cancer_table, raw_bal_counts_cancer_sig_table,
                   raw_unbal_cancer_table, raw_unbal_cancer_sig_table, raw_bal_p53_table, raw_bal_p53_sig_table, raw_bal_counts_p53_table, 
                   raw_bal_counts_p53_sig_table, raw_unbal_p53_table, raw_unbal_p53_sig_table, raw_cancer_intersection_table, 
                   raw_cancer_intersection_sig_table, raw_cancer_union_table, raw_cancer_union_sig_table, raw_p53_intersection_table, 
                   raw_p53_intersection_sig_table, raw_p53_union_table, raw_p53_union_sig_table,
                   raw_bal_counts_cancer_intersection_table, raw_bal_counts_cancer_intersection_sig_table,
                   raw_bal_counts_cancer_union_table, raw_bal_counts_cancer_union_sig_table,
                   raw_bal_counts_p53_intersection_table, raw_bal_counts_p53_intersection_sig_table,
                   raw_bal_counts_p53_union_table, raw_bal_counts_p53_union_sig_table,
                   raw_complete_cancer_intersection_table, raw_complete_cancer_intersection_sig_table,
                   raw_complete_cancer_union_table, raw_complete_cancer_union_sig_table,
                   raw_complete_p53_intersection_table, raw_complete_p53_intersection_sig_table,
                   raw_complete_p53_union_table, raw_complete_p53_union_sig_table)

# remove data 
rm(raw_bal_cancer_table, raw_bal_cancer_sig_table, raw_bal_counts_cancer_table, raw_bal_counts_cancer_sig_table,
   raw_unbal_cancer_table, raw_unbal_cancer_sig_table, raw_bal_p53_table, raw_bal_p53_sig_table, raw_bal_counts_p53_table, 
   raw_bal_counts_p53_sig_table, raw_unbal_p53_table, raw_unbal_p53_sig_table, raw_cancer_intersection_table, 
   raw_cancer_intersection_sig_table, raw_cancer_union_table, raw_cancer_union_sig_table, raw_p53_intersection_table, 
   raw_p53_intersection_sig_table, raw_p53_union_table, raw_p53_union_sig_table,
   raw_bal_counts_cancer_intersection_table, raw_bal_counts_cancer_intersection_sig_table,
   raw_bal_counts_cancer_union_table, raw_bal_counts_cancer_union_sig_table,
   raw_bal_counts_p53_intersection_table, raw_bal_counts_p53_intersection_sig_table,
   raw_bal_counts_p53_union_table, raw_bal_counts_p53_union_sig_table,
   raw_complete_cancer_intersection_table, raw_complete_cancer_intersection_sig_table,
   raw_complete_cancer_union_table, raw_complete_cancer_union_sig_table,
   raw_complete_p53_intersection_table, raw_complete_p53_intersection_sig_table,
   raw_complete_p53_union_table, raw_complete_p53_union_sig_table)


#save table 
saveRDS(raw_table, 
        file = paste0(raw_folder, '/raw_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(raw_bal_cancer_models, 
        file = paste0(raw_folder, '/raw_bal_cancer_models.rda'))
saveRDS(raw_bal_cancer_sig_models, 
        file = paste0(raw_folder, '/raw_bal_cancer_sig_models.rda'))
saveRDS(raw_bal_counts_cancer_models, 
        file = paste0(raw_folder, '/raw_bal_counts_cancer_models.rda'))
saveRDS(raw_bal_counts_cancer_sig_models, 
        file = paste0(raw_folder, '/raw_bal_counts_cancer_sig_models.rda'))
saveRDS(raw_unbal_cancer_models, 
        file = paste0(raw_folder, '/raw_unbal_cancer_models.rda'))
saveRDS(raw_unbal_cancer_sig_models, 
        file = paste0(raw_folder, '/raw_unbal_cancer_sig_models.rda'))
saveRDS(raw_bal_p53_models, 
        file = paste0(raw_folder, '/raw_bal_p53_models.rda'))
saveRDS(raw_bal_p53_sig_models, 
        file = paste0(raw_folder, '/raw_bal_p53_sig_models.rda'))
saveRDS(raw_bal_counts_p53_models, 
        file = paste0(raw_folder, '/raw_bal_counts_p53_models.rda'))
saveRDS(raw_bal_counts_p53_sig_models, 
        file = paste0(raw_folder, '/raw_bal_counts_p53_sig_models.rda'))
saveRDS(raw_unbal_p53_models, 
        file = paste0(raw_folder, '/raw_unbal_p53_models.rda'))
saveRDS(raw_unbal_p53_sig_models, 
        file = paste0(raw_folder, '/raw_unbal_p53_sig_models.rda'))
saveRDS(raw_cancer_intersection_models, 
        file = paste0(raw_folder, '/raw_cancer_intersection_models.rda'))
saveRDS(raw_cancer_intersection_sig_models, 
        file = paste0(raw_folder, '/raw_cancer_intersection_sig_models.rda'))
saveRDS(raw_cancer_union_models, 
        file = paste0(raw_folder, '/raw_cancer_union_models.rda'))
saveRDS(raw_cancer_union_sig_models, 
        file = paste0(raw_folder, '/raw_cancer_union_sig_models.rda'))
saveRDS(raw_p53_intersection_models, 
        file = paste0(raw_folder, '/raw_p53_intersection_models.rda'))
saveRDS(raw_p53_intersection_sig_models, 
        file = paste0(raw_folder, '/raw_p53_intersection_sig_models.rda'))
saveRDS(raw_p53_union_models, 
        file = paste0(raw_folder, '/raw_p53_union_models.rda'))
saveRDS(raw_p53_union_sig_models, 
        file = paste0(raw_folder, '/raw_p53_union_sig_models.rda'))
saveRDS(raw_bal_counts_cancer_intersection_models,
        file = paste0(raw_folder, '/raw_bal_counts_cancer_intersection_models.rda'))
saveRDS(raw_bal_counts_cancer_intersection_sig_models,
        file = paste0(raw_folder, '/raw_bal_counts_cancer_intersection_sig_models.rda'))
saveRDS(raw_bal_counts_cancer_union_models,
        file = paste0(raw_folder, '/raw_bal_counts_cancer_union_models.rda'))
saveRDS(raw_bal_counts_cancer_union_sig_models,
        file = paste0(raw_folder, '/raw_bal_counts_cancer_union_sig_models.rda'))
saveRDS(raw_bal_counts_p53_intersection_models,
        file = paste0(raw_folder, '/raw_bal_counts_p53_intersection_models.rda'))
saveRDS(raw_bal_counts_p53_intersection_sig_models,
        file = paste0(raw_folder, '/raw_bal_counts_p53_intersection_sig_models.rda'))
saveRDS(raw_bal_counts_p53_union_models,
        file = paste0(raw_folder, '/raw_bal_counts_p53_union_models.rda'))
saveRDS(raw_bal_counts_p53_union_sig_models,
        file = paste0(raw_folder, '/raw_bal_counts_p53_union_sig_models.rda'))
saveRDS(raw_complete_cancer_intersection_models,
        file = paste0(raw_folder, '/raw_complete_cancer_intersection_models.rda'))
# saveRDS(raw_complete_cancer_intersection_sig_models,
#         file = paste0(raw_folder, '/raw_complete_cancer_intersection_sig_models.rda'))
saveRDS(raw_complete_cancer_union_models,
        file = paste0(raw_folder, '/raw_complete_cancer_union_models.rda'))
saveRDS(raw_complete_cancer_union_sig_models,
        file = paste0(raw_folder, '/raw_complete_cancer_union_sig_models.rda'))
saveRDS(raw_complete_p53_intersection_models,
        file = paste0(raw_folder, '/raw_complete_p53_intersection_models.rda'))
saveRDS(raw_complete_p53_intersection_sig_models,
        file = paste0(raw_folder, '/raw_complete_p53_intersection_sig_models.rda'))
saveRDS(raw_complete_p53_union_models,
        file = paste0(raw_folder, '/raw_complete_p53_union_models.rda'))
saveRDS(raw_complete_p53_union_sig_models,
        file = paste0(raw_folder, '/raw_complete_p53_union_sig_models.rda'))



