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
source(paste0(scripts_folder, '/predict_age/model_functions.R'))
source(paste0(scripts_folder, '/predict_age/run_models.R'))

##########
# remove raw, swan, and quan objects
##########
rm(list = ls(pattern = "beta_raw_*"))
rm(list = ls(pattern = "beta_swan_*"))
rm(list = ls(pattern = "beta_quan_*"))


###################################################################################################################################
## beta_funnorm

##########
# cancer
##########

# bal cancer
funnorm_bal_cancer_models <- runModels(beta_funnorm, 
                                   bump_hunter = T, 
                                   bump_hunter_data = beta_funnorm_bal_cancer_features)

funnorm_bal_cancer_table <- extractResults(funnorm_bal_cancer_models, 
                                           data_name = 'beta_funnorm_bal_cancer')

# bal cancer sig
funnorm_bal_cancer_sig_models <- runModels(beta_funnorm, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_funnorm_bal_cancer_sig_features)

funnorm_bal_cancer_sig_table <- extractResults(funnorm_bal_cancer_sig_models, 
                                               data_name = 'funnorm_bal_cancer_sig')


# bal counts cancer
funnorm_bal_counts_cancer_models <- runModels(beta_funnorm, 
                                              bump_hunter = T, 
                                              bump_hunter_data = beta_funnorm_bal_counts_cancer_features)

funnorm_bal_counts_cancer_table <- extractResults(funnorm_bal_counts_cancer_models, 
                                              data_name = 'funnorm_bal_counts_cancer')

# bal counts cancer sig
funnorm_bal_counts_cancer_sig_models <- runModels(beta_funnorm, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data =  beta_funnorm_bal_counts_cancer_sig_features)

funnorm_bal_counts_cancer_sig_table <- extractResults(funnorm_bal_counts_cancer_sig_models, 
                                                      data_name = 'funnorm_bal_counts_cancer_sig')

# unbal cancer
funnorm_unbal_cancer_models <- runModels(beta_funnorm, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_funnorm_unbal_cancer_features)

funnorm_unbal_cancer_table <- extractResults(funnorm_unbal_cancer_models, 
                                         data_name = 'funnorm_unbal_cancer')


# unbal cancer sig
funnorm_unbal_cancer_sig_models <- runModels(beta_funnorm, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_funnorm_unbal_cancer_sig_features)

funnorm_unbal_cancer_sig_table <- extractResults(funnorm_unbal_cancer_sig_models, 
                                             data_name = 'funnorm_unbal_cancer_sig')

##########
# p53
##########

# bal p53
funnorm_bal_p53_models <- runModels(beta_funnorm, 
                                bump_hunter = T, 
                                bump_hunter_data = beta_funnorm_bal_p53_features)

funnorm_bal_p53_table <- extractResults(funnorm_bal_p53_models, 
                                    data_name = 'funnorm_bal_p53')

# bal p53 sig
funnorm_bal_p53_sig_models <- runModels(beta_funnorm, 
                                    bump_hunter = T, 
                                    bump_hunter_data = beta_funnorm_bal_p53_sig_features)

funnorm_bal_p53_sig_table <- extractResults(funnorm_bal_p53_sig_models, 
                                        data_name = 'funnorm_bal_p53_sig')

# bal counts p53
funnorm_bal_counts_p53_models <- runModels(beta_funnorm, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_funnorm_bal_counts_p53_features)

funnorm_bal_counts_p53_table <- extractResults(funnorm_bal_counts_p53_models, 
                                           data_name = 'funnorm_bal_counts_p53')

# bal counts p53 sig
funnorm_bal_counts_p53_sig_models <- runModels(beta_funnorm, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_funnorm_bal_counts_p53_sig_features)

funnorm_bal_counts_p53_sig_table <- extractResults(funnorm_bal_counts_p53_sig_models, 
                                               data_name = 'funnorm_bal_counts_p53_sig')

# unbal p53
funnorm_unbal_p53_models <- runModels(beta_funnorm, 
                                  bump_hunter = T, 
                                  bump_hunter_data = beta_funnorm_unbal_p53_features)

funnorm_unbal_p53_table <- extractResults(funnorm_unbal_p53_models, 
                                      data_name = 'funnorm_unbal_p53')


# unbal p53 sig
funnorm_unbal_p53_sig_models <- runModels(beta_funnorm, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_funnorm_unbal_p53_sig_features)

funnorm_unbal_p53_sig_table <- extractResults(funnorm_unbal_p53_sig_models, 
                                          data_name = 'funnorm_unbal_p53_sig')


##########
# cancer
##########

# cancer_intersection
funnorm_cancer_intersection_models <- runModels(beta_funnorm, 
                                            bump_hunter = T, 
                                            bump_hunter_data = beta_funnorm_cancer_intersection_features)

funnorm_cancer_intersection_table <- extractResults(funnorm_cancer_intersection_models, 
                                                data_name = 'funnorm_cancer_intersection')

# cancer_intersection sig
funnorm_cancer_intersection_sig_models <- runModels(beta_funnorm, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_funnorm_cancer_intersection_sig_features)

funnorm_cancer_intersection_sig_table <- extractResults(funnorm_cancer_intersection_sig_models, 
                                                    data_name = 'funnorm_cancer_intersection_sig')

# cancer_union
funnorm_cancer_union_models <- runModels(beta_funnorm, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_funnorm_cancer_union_features)

funnorm_cancer_union_table <- extractResults(funnorm_cancer_union_models, 
                                         data_name = 'funnorm_cancer_union')

# cancer_union sig
funnorm_cancer_union_sig_models <- runModels(beta_funnorm, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_funnorm_cancer_union_sig_features)

funnorm_cancer_union_sig_table <- extractResults(funnorm_cancer_union_sig_models, 
                                             data_name = 'funnorm_cancer_union_sig')
##########
# p53
##########
# p53_intersection
funnorm_p53_intersection_models <- runModels(beta_funnorm, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_funnorm_p53_intersection_features)

funnorm_p53_intersection_table <- extractResults(funnorm_p53_intersection_models, 
                                             data_name = 'funnorm_p53_intersection')

# p53_intersection sig
funnorm_p53_intersection_sig_models <- runModels(beta_funnorm, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_funnorm_p53_intersection_sig_features)

funnorm_p53_intersection_sig_table <- extractResults(funnorm_p53_intersection_sig_models, 
                                                 data_name = 'funnorm_p53_intersection_sig')

# p53_union
funnorm_p53_union_models <- runModels(beta_funnorm, 
                                  bump_hunter = T, 
                                  bump_hunter_data = beta_funnorm_p53_union_features)

funnorm_p53_union_table <- extractResults(funnorm_p53_union_models, 
                                      data_name = 'funnorm_p53_union')

# p53_union sig
funnorm_p53_union_sig_models <- runModels(beta_funnorm, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_funnorm_p53_union_sig_features)

funnorm_p53_union_sig_table <- extractResults(funnorm_p53_union_sig_models, 
                                          data_name = 'funnorm_p53_union_sig')


##########
# most balanced cancer 
##########

# intersection
funnorm_bal_counts_cancer_intersection_models <- runModels(beta_funnorm, 
                                                       bump_hunter = T, 
                                                       bump_hunter_data = beta_bal_counts_cancer_intersection_features)

funnorm_bal_counts_cancer_intersection_table <- extractResults(funnorm_bal_counts_cancer_intersection_models, 
                                                           data_name = 'funnorm_bal_counts_cancer_intersection')


funnorm_bal_counts_cancer_intersection_sig_models <- runModels(beta_funnorm, 
                                                           bump_hunter = T, 
                                                           bump_hunter_data = beta_bal_counts_cancer_intersection_sig_features)

funnorm_bal_counts_cancer_intersection_sig_table <- extractResults(funnorm_bal_counts_cancer_intersection_sig_models, 
                                                               data_name = 'funnorm_bal_counts_cancer_intersection_sig')

# union
funnorm_bal_counts_cancer_union_models <- runModels(beta_funnorm, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_bal_counts_cancer_union_features)

funnorm_bal_counts_cancer_union_table <- extractResults(funnorm_bal_counts_cancer_union_models, 
                                                    data_name = 'funnorm_bal_counts_cancer_union')


funnorm_bal_counts_cancer_union_sig_models <- runModels(beta_funnorm, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_cancer_union_sig_features)

funnorm_bal_counts_cancer_union_sig_table <- extractResults(funnorm_bal_counts_cancer_union_sig_models, 
                                                        data_name = 'funnorm_bal_counts_cancer_union_sig')

##########
# most balanced p53 
##########

# intersection
funnorm_bal_counts_p53_intersection_models <- runModels(beta_funnorm, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_p53_intersection_features)

funnorm_bal_counts_p53_intersection_table <- extractResults(funnorm_bal_counts_p53_intersection_models, 
                                                        data_name = 'funnorm_bal_counts_p53_intersection')


funnorm_bal_counts_p53_intersection_sig_models <- runModels(beta_funnorm, 
                                                        bump_hunter = T, 
                                                        bump_hunter_data = beta_bal_counts_p53_intersection_sig_features)

funnorm_bal_counts_p53_intersection_sig_table <- extractResults(funnorm_bal_counts_p53_intersection_sig_models, 
                                                            data_name = 'funnorm_bal_counts_p53_intersection_sig')


# union
funnorm_bal_counts_p53_union_models <- runModels(beta_funnorm, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_bal_counts_p53_union_features)

funnorm_bal_counts_p53_union_table <- extractResults(funnorm_bal_counts_p53_union_models, 
                                                 data_name = 'funnorm_bal_counts_p53_union')


funnorm_bal_counts_p53_union_sig_models <- runModels(beta_funnorm, 
                                                 bump_hunter = T, 
                                                 bump_hunter_data = beta_bal_counts_p53_union_sig_features)

funnorm_bal_counts_p53_union_sig_table <- extractResults(funnorm_bal_counts_p53_union_sig_models, 
                                                     data_name = 'funnorm_bal_counts_p53_union_sig')

###########
# complete cancer - intersection across each method intersection
###########

# complete cancer intersection
funnorm_complete_cancer_intersection_models <- runModels(beta_funnorm, 
                                                     bump_hunter = T, 
                                                     bump_hunter_data = beta_cancer_intersection_features)

# get table 
funnorm_complete_cancer_intersection_table <- extractResults(funnorm_complete_cancer_intersection_models, 
                                                         data_name = 'funnorm_complete_cancer_intersection')

# empty features set
# # complete cancer intersection sig
# funnorm_complete_cancer_intersection_sig_models <- runModels(beta_funnorm, 
#                                                          bump_hunter = T, 
#                                                          bump_hunter_data = beta_cancer_intersection_sig_features)
# 
# # get table 
# funnorm_complete_cancer_intersection_sig_table <- extractResults(funnorm_complete_cancer_intersection_sig_models, 
#                                                              data_name = 'funnorm_complete_cancer_intersection_sig')


# complete cancer union
funnorm_complete_cancer_union_models <- runModels(beta_funnorm, 
                                              bump_hunter = T, 
                                              bump_hunter_data = beta_cancer_union_features)

# get table 
funnorm_complete_cancer_union_table <- extractResults(funnorm_complete_cancer_union_models, 
                                                  data_name = 'funnorm_complete_cancer_union')

# complete cancer union sig
funnorm_complete_cancer_union_sig_models <- runModels(beta_funnorm, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_cancer_union_sig_features)

# get table 
funnorm_complete_cancer_union_sig_table <- extractResults(funnorm_complete_cancer_union_sig_models, 
                                                      data_name = 'funnorm_complete_cancer_union_sig')


###########
# complete p53 - intersection across each method intersection
###########

# complete p53 intersection
funnorm_complete_p53_intersection_models <- runModels(beta_funnorm, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_p53_intersection_features)

# get table 
funnorm_complete_p53_intersection_table <- extractResults(funnorm_complete_p53_intersection_models, 
                                                      data_name = 'funnorm_complete_p53_intersection')

# complete p53 intersection sig
funnorm_complete_p53_intersection_sig_models <- runModels(beta_funnorm, 
                                                      bump_hunter = T, 
                                                      bump_hunter_data = beta_p53_intersection_sig_features)

# get table 
funnorm_complete_p53_intersection_sig_table <- extractResults(funnorm_complete_p53_intersection_sig_models, 
                                                          data_name = 'funnorm_complete_p53_intersection_sig')


# complete p53 union
funnorm_complete_p53_union_models <- runModels(beta_funnorm, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_p53_union_features)

# get table 
funnorm_complete_p53_union_table <- extractResults(funnorm_complete_p53_union_models, 
                                               data_name = 'funnorm_complete_p53_union')

# complete p53 union sig
funnorm_complete_p53_union_sig_models <- runModels(beta_funnorm, 
                                               bump_hunter = T, 
                                               bump_hunter_data = beta_p53_union_sig_features)

# get table 
funnorm_complete_p53_union_sig_table <- extractResults(funnorm_complete_p53_union_sig_models, 
                                                   data_name = 'funnorm_complete_p53_union_sig')


###########
# rbind tables and save RDA file
###########
funnorm_table <- rbind(funnorm_bal_cancer_table, funnorm_bal_cancer_sig_table, funnorm_bal_counts_cancer_table, funnorm_bal_counts_cancer_sig_table,
                   funnorm_unbal_cancer_table, funnorm_unbal_cancer_sig_table, funnorm_bal_p53_table, funnorm_bal_p53_sig_table, funnorm_bal_counts_p53_table, 
                   funnorm_bal_counts_p53_sig_table, funnorm_unbal_p53_table, funnorm_unbal_p53_sig_table, funnorm_cancer_intersection_table, 
                   funnorm_cancer_intersection_sig_table, funnorm_cancer_union_table, funnorm_cancer_union_sig_table, funnorm_p53_intersection_table, 
                   funnorm_p53_intersection_sig_table, funnorm_p53_union_table, funnorm_p53_union_sig_table,
                   funnorm_bal_counts_cancer_intersection_table, funnorm_bal_counts_cancer_intersection_sig_table,
                   funnorm_bal_counts_cancer_union_table, funnorm_bal_counts_cancer_union_sig_table,
                   funnorm_bal_counts_p53_intersection_table, funnorm_bal_counts_p53_intersection_sig_table,
                   funnorm_bal_counts_p53_union_table, funnorm_bal_counts_p53_union_sig_table,
                   funnorm_complete_cancer_intersection_table, #funnorm_complete_cancer_intersection_sig_table,
                   funnorm_complete_cancer_union_table, funnorm_complete_cancer_union_sig_table,
                   funnorm_complete_p53_intersection_table, funnorm_complete_p53_intersection_sig_table,
                   funnorm_complete_p53_union_table, funnorm_complete_p53_union_sig_table)

# remove data 
rm(funnorm_bal_cancer_table, funnorm_bal_cancer_sig_table, funnorm_bal_counts_cancer_table, funnorm_bal_counts_cancer_sig_table,
   funnorm_unbal_cancer_table, funnorm_unbal_cancer_sig_table, funnorm_bal_p53_table, funnorm_bal_p53_sig_table, funnorm_bal_counts_p53_table, 
   funnorm_bal_counts_p53_sig_table, funnorm_unbal_p53_table, funnorm_unbal_p53_sig_table, funnorm_cancer_intersection_table, 
   funnorm_cancer_intersection_sig_table, funnorm_cancer_union_table, funnorm_cancer_union_sig_table, funnorm_p53_intersection_table, 
   funnorm_p53_intersection_sig_table, funnorm_p53_union_table, funnorm_p53_union_sig_table,
   funnorm_bal_counts_cancer_intersection_table, funnorm_bal_counts_cancer_intersection_sig_table,
   funnorm_bal_counts_cancer_union_table, funnorm_bal_counts_cancer_union_sig_table,
   funnorm_bal_counts_p53_intersection_table, funnorm_bal_counts_p53_intersection_sig_table,
   funnorm_bal_counts_p53_union_table, funnorm_bal_counts_p53_union_sig_table,
   funnorm_complete_cancer_intersection_table, #funnorm_complete_cancer_intersection_sig_table,
   funnorm_complete_cancer_union_table, funnorm_complete_cancer_union_sig_table,
   funnorm_complete_p53_intersection_table, funnorm_complete_p53_intersection_sig_table,
   funnorm_complete_p53_union_table, funnorm_complete_p53_union_sig_table)


#save table 
saveRDS(funnorm_table, 
        file = paste0(funnorm_folder, '/funnorm_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(funnorm_bal_cancer_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_cancer_models.rda'))
saveRDS(funnorm_bal_cancer_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_cancer_sig_models.rda'))
saveRDS(funnorm_bal_counts_cancer_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_counts_cancer_models.rda'))
saveRDS(funnorm_bal_counts_cancer_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_counts_cancer_sig_models.rda'))
saveRDS(funnorm_unbal_cancer_models, 
        file = paste0(funnorm_folder, '/funnorm_unbal_cancer_models.rda'))
saveRDS(funnorm_unbal_cancer_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_unbal_cancer_sig_models.rda'))
saveRDS(funnorm_bal_p53_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_p53_models.rda'))
saveRDS(funnorm_bal_p53_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_p53_sig_models.rda'))
saveRDS(funnorm_bal_counts_p53_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_counts_p53_models.rda'))
saveRDS(funnorm_bal_counts_p53_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_counts_p53_sig_models.rda'))
saveRDS(funnorm_unbal_p53_models, 
        file = paste0(funnorm_folder, '/funnorm_unbal_p53_models.rda'))
saveRDS(funnorm_unbal_p53_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_unbal_p53_sig_models.rda'))
saveRDS(funnorm_cancer_intersection_models, 
        file = paste0(funnorm_folder, '/funnorm_cancer_intersection_models.rda'))
saveRDS(funnorm_cancer_intersection_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_cancer_intersection_sig_models.rda'))
saveRDS(funnorm_cancer_union_models, 
        file = paste0(funnorm_folder, '/funnorm_cancer_union_models.rda'))
saveRDS(funnorm_cancer_union_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_cancer_union_sig_models.rda'))
saveRDS(funnorm_p53_intersection_models, 
        file = paste0(funnorm_folder, '/funnorm_p53_intersection_models.rda'))
saveRDS(funnorm_p53_intersection_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_p53_intersection_sig_models.rda'))
saveRDS(funnorm_p53_union_models, 
        file = paste0(funnorm_folder, '/funnorm_p53_union_models.rda'))
saveRDS(funnorm_p53_union_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_p53_union_sig_models.rda'))
saveRDS(funnorm_bal_counts_cancer_intersection_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_cancer_intersection_models.rda'))
saveRDS(funnorm_bal_counts_cancer_intersection_sig_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_cancer_intersection_sig_models.rda'))
saveRDS(funnorm_bal_counts_cancer_union_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_cancer_union_models.rda'))
saveRDS(funnorm_bal_counts_cancer_union_sig_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_cancer_union_sig_models.rda'))
saveRDS(funnorm_bal_counts_p53_intersection_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_p53_intersection_models.rda'))
saveRDS(funnorm_bal_counts_p53_intersection_sig_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_p53_intersection_sig_models.rda'))
saveRDS(funnorm_bal_counts_p53_union_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_p53_union_models.rda'))
saveRDS(funnorm_bal_counts_p53_union_sig_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_p53_union_sig_models.rda'))
saveRDS(funnorm_complete_cancer_intersection_models,
        file = paste0(funnorm_folder, '/funnorm_complete_cancer_intersection_models.rda'))
# saveRDS(funnorm_complete_cancer_intersection_sig_models,
#         file = paste0(funnorm_folder, '/funnorm_complete_cancer_intersection_sig_models.rda'))
saveRDS(funnorm_complete_cancer_union_models,
        file = paste0(funnorm_folder, '/funnorm_complete_cancer_union_models.rda'))
saveRDS(funnorm_complete_cancer_union_sig_models,
        file = paste0(funnorm_folder, '/funnorm_complete_cancer_union_sig_models.rda'))
saveRDS(funnorm_complete_p53_intersection_models,
        file = paste0(funnorm_folder, '/funnorm_complete_p53_intersection_models.rda'))
saveRDS(funnorm_complete_p53_intersection_sig_models,
        file = paste0(funnorm_folder, '/funnorm_complete_p53_intersection_sig_models.rda'))
saveRDS(funnorm_complete_p53_union_models,
        file = paste0(funnorm_folder, '/funnorm_complete_p53_union_models.rda'))
saveRDS(funnorm_complete_p53_union_sig_models,
        file = paste0(funnorm_folder, '/funnorm_complete_p53_union_sig_models.rda'))



