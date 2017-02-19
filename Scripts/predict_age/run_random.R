####### This script will predict age of onset with raw preprocessed data
# this is part of 7th step in pipleline


##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')
scripts_folder <- paste0(project_folder, '/Scripts')

##########
# source model_functions.R and run_models.R
##########
source(paste0(scripts_folder, '/predict_age/model_functions_short.R'))
source(paste0(scripts_folder, '/predict_age/run_models_short.R'))



#########################################################################################################################
# Random features - 100, 200, 500, 1000, 2000

###########
# function for random
###########
# data <- beta_raw
# rand_num <- 2
# num_feat <- 100
# data_name <- 'raw_rand_100'
# exlcude <- T
# union_features <-  beta_raw_union_features
getRand <- function(data, rand_num, num_feat, exclude, union_features, data_name) 
{
  
  if (exclude) {
    
    features <- colnames(data[, 8:ncol(data)])
    keep_these <- features[!features %in% union_features]
    data <- data[, c('id', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                     'age_sample_collection', keep_these)]
  }
  
  model_holder = list()
  table_holder = list()
  
  for (rand_num in 1:rand_num) {
    
    model_holder[[rand_num]] <- runModels(data, 
                                          bump_hunter = F,
                                          num_feat = num_feat,
                                          random = T,
                                          seed_num = rand_num)
    
    table_holder[[rand_num]] <- extractResults(model_holder[[rand_num]],
                                               data_name = data_name)
    
  }
  
  full_table <- do.call(rbind, table_holder)
  
  return(full_table)
}

###########
# get all features
###########
union_features <- Reduce(union, list(raw_union_feat, 
                                     swan_union_feat,
                                     quan_union_feat,
                                     funnorm_union_feat))
#############################################################################################################
###########
# raw
###########

# raw 07 features
raw_rand_07_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_07), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "raw_rand_07")

# raw 08 features
raw_rand_08_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_08), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "raw_rand_08")

# raw 09 features
raw_rand_09_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_09), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "raw_rand_09")


# raw 10 features
raw_rand_10_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_10), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "raw_rand_10")

# raw 11 features
raw_rand_11_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_11), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "raw_rand_11")

# raw 12 features
raw_rand_12_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_12), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "raw_rand_12")

# raw 13 features
raw_rand_13_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_13), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "raw_rand_13")

# raw 14 features
raw_rand_14_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_14), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "raw_rand_14")

# raw 15 features
raw_rand_15_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_15), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "raw_rand_15")

###########
# raw sig
###########

# raw 07 features
raw_rand_07_sig_table <- getRand(raw_cases, 
                                 rand_num = 5, 
                                 num_feat = length(raw_bh_sig_07), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "raw_rand_sig_07")

# raw 08 features
raw_rand_08_sig_table <- getRand(raw_cases, 
                                 rand_num = 5, 
                                 num_feat = length(raw_bh_sig_08), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "raw_rand_sig_08")

# raw 09 features
raw_rand_09_sig_table <- getRand(raw_cases, 
                                 rand_num = 5, 
                                 num_feat = length(raw_bh_sig_09), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "raw_rand_sig_09")


# raw 10 features
raw_rand_10_sig_table <- getRand(raw_cases, 
                                 rand_num = 5, 
                                 num_feat = length(raw_bh_sig_10), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "raw_rand_sig_10")

# raw 11 features
raw_rand_11_sig_table <- getRand(raw_cases, 
                                 rand_num = 5, 
                                 num_feat = length(raw_bh_sig_11), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "raw_rand_sig_11")

# raw 12 features
raw_rand_12_sig_table <- getRand(raw_cases, 
                                 rand_num = 5, 
                                 num_feat = length(raw_bh_sig_12), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "raw_rand_sig_12")

# raw 13 features
raw_rand_13_sig_table <- getRand(raw_cases, 
                                 rand_num = 5, 
                                 num_feat = length(raw_bh_sig_13), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "raw_rand_sig_13")

# raw 14 features
raw_rand_14_sig_table <- getRand(raw_cases, 
                                 rand_num = 5, 
                                 num_feat = length(raw_bh_sig_14), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "raw_rand_sig_14")

# raw 15 features
raw_rand_15_sig_table <- getRand(raw_cases, 
                                 rand_num = 5, 
                                 num_feat = length(raw_bh_sig_15), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "raw_rand_sig_15")

##########
# intersection and union
##########

# raw int features
raw_rand_int_table <- getRand(raw_cases, 
                              rand_num = 5, 
                              num_feat = length(raw_int_feat), 
                              exclude = T,
                              union_features = union_features,
                              data_name = "raw_rand_int")

# raw int sig features
raw_rand_int_sig_table <- getRand(raw_cases, 
                                  rand_num = 5, 
                                  num_feat = length(raw_int_sig_feat), 
                                  exclude = T,
                                  union_features = union_features,
                                  data_name = "raw_rand_int_sig")

# raw union features
raw_rand_union_table <- getRand(raw_cases, 
                                rand_num = 5, 
                                num_feat = length(raw_union_feat), 
                                exclude = T,
                                union_features = union_features,
                                data_name = "raw_rand_union")

# raw union sig features
raw_rand_union_sig_table <- getRand(raw_cases, 
                                    rand_num = 5, 
                                    num_feat = length(raw_union_sig_feat), 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "raw_rand_union_sig")



###########
# rbind tables and save RDA file
###########
raw_rand_table <- rbind(raw_rand_07_table, raw_rand_07_sig_table,
                        raw_rand_08_table, raw_rand_08_sig_table,
                        raw_rand_09_table, raw_rand_09_sig_table,
                        raw_rand_10_table, raw_rand_10_sig_table,
                        raw_rand_11_table, raw_rand_11_sig_table,
                        raw_rand_12_table, raw_rand_12_sig_table,
                        raw_rand_13_table, raw_rand_13_sig_table,
                        raw_rand_14_table, raw_rand_14_sig_table,
                        raw_rand_15_table, raw_rand_15_sig_table,
                        raw_rand_int_table, raw_rand_int_sig_table,
                        raw_rand_union_table, raw_rand_union_sig_table)

# remove data 
rm(raw_rand_07_table, raw_rand_07_sig_table,
   raw_rand_08_table, raw_rand_08_sig_table,
   raw_rand_09_table, raw_rand_09_sig_table,
   raw_rand_10_table, raw_rand_10_sig_table,
   raw_rand_11_table, raw_rand_11_sig_table,
   raw_rand_12_table, raw_rand_12_sig_table,
   raw_rand_13_table, raw_rand_13_sig_table,
   raw_rand_14_table, raw_rand_14_sig_table,
   raw_rand_15_table, raw_rand_15_sig_table,
   raw_rand_int_table, raw_rand_int_sig_table,
   raw_rand_union_table, raw_rand_union_sig_table)


#save table 
saveRDS(raw_rand_table, 
        file = paste0(rand_folder, '/raw_rand_table.rda'))

rm(list = ls(pattern = "raw_*"))


#############################################################################################################
###########
# quan
###########

# quan 07 features
quan_rand_07_table <- getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_bh_07), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "quan_rand_07")

# quan 08 features
quan_rand_08_table <- getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_bh_08), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "quan_rand_08")

# quan 09 features
quan_rand_09_table <- getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_bh_09), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "quan_rand_09")


# quan 10 features
quan_rand_10_table <- getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_bh_10), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "quan_rand_10")

# quan 11 features
quan_rand_11_table <- getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_bh_11), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "quan_rand_11")

# quan 12 features
quan_rand_12_table <- getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_bh_12), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "quan_rand_12")

# quan 13 features
quan_rand_13_table <- getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_bh_13), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "quan_rand_13")

# quan 14 features
quan_rand_14_table <- getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_bh_14), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "quan_rand_14")

# quan 15 features
quan_rand_15_table <- getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_bh_15), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "quan_rand_15")

###########
# quan sig
###########

# quan 07 features
quan_rand_07_sig_table <- getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_bh_sig_07), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_rand_sig_07")

# quan 08 features
quan_rand_08_sig_table <- getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_bh_sig_08), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_rand_sig_08")

# quan 09 features
quan_rand_09_sig_table <- getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_bh_sig_09), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_rand_sig_09")


# quan 10 features
quan_rand_10_sig_table <- getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_bh_sig_10), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_rand_sig_10")

# quan 11 features
quan_rand_11_sig_table <- getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_bh_sig_11), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_rand_sig_11")

# quan 12 features
quan_rand_12_sig_table <- getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_bh_sig_12), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_rand_sig_12")

# quan 13 features
quan_rand_13_sig_table <- getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_bh_sig_13), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_rand_sig_13")

# quan 14 features
quan_rand_14_sig_table <- getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_bh_sig_14), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_rand_sig_14")

# quan 15 features
quan_rand_15_sig_table <- getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_bh_sig_15), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_rand_sig_15")

##########
# intersection and union
##########

# quan int features
quan_rand_int_table <- getRand(quan_cases, 
                              rand_num = 5, 
                              num_feat = length(quan_int_feat), 
                              exclude = T,
                              union_features = union_features,
                              data_name = "quan_rand_int")

# quan int sig features
quan_rand_int_sig_table <- getRand(quan_cases, 
                                  rand_num = 5, 
                                  num_feat = length(quan_int_sig_feat), 
                                  exclude = T,
                                  union_features = union_features,
                                  data_name = "quan_rand_int_sig")

# quan union features
quan_rand_union_table <- getRand(quan_cases, 
                                rand_num = 5, 
                                num_feat = length(quan_union_feat), 
                                exclude = T,
                                union_features = union_features,
                                data_name = "quan_rand_union")

# quan union sig features
quan_rand_union_sig_table <- getRand(quan_cases, 
                                    rand_num = 5, 
                                    num_feat = length(quan_union_sig_feat), 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_rand_union_sig")



###########
# rbind tables and save RDA file
###########
quan_rand_table <- rbind(quan_rand_07_table, quan_rand_07_sig_table,
                        quan_rand_08_table, quan_rand_08_sig_table,
                        quan_rand_09_table, quan_rand_09_sig_table,
                        quan_rand_10_table, quan_rand_10_sig_table,
                        quan_rand_11_table, quan_rand_11_sig_table,
                        quan_rand_12_table, quan_rand_12_sig_table,
                        quan_rand_13_table, quan_rand_13_sig_table,
                        quan_rand_14_table, quan_rand_14_sig_table,
                        quan_rand_15_table, quan_rand_15_sig_table,
                        quan_rand_int_table, quan_rand_int_sig_table,
                        quan_rand_union_table, quan_rand_union_sig_table)

# remove data 
rm(quan_rand_07_table, quan_rand_07_sig_table,
   quan_rand_08_table, quan_rand_08_sig_table,
   quan_rand_09_table, quan_rand_09_sig_table,
   quan_rand_10_table, quan_rand_10_sig_table,
   quan_rand_11_table, quan_rand_11_sig_table,
   quan_rand_12_table, quan_rand_12_sig_table,
   quan_rand_13_table, quan_rand_13_sig_table,
   quan_rand_14_table, quan_rand_14_sig_table,
   quan_rand_15_table, quan_rand_15_sig_table,
   quan_rand_int_table, quan_rand_int_sig_table,
   quan_rand_union_table, quan_rand_union_sig_table)


#save table 
saveRDS(quan_rand_table, 
        file = paste0(rand_folder, '/quan_rand_table.rda'))

rm(list = ls(pattern = "quan_*"))


#############################################################################################################
###########
# swan
###########

# swan 07 features
swan_rand_07_table <- getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_bh_07), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "swan_rand_07")

# swan 08 features
swan_rand_08_table <- getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_bh_08), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "swan_rand_08")

# swan 09 features
swan_rand_09_table <- getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_bh_09), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "swan_rand_09")


# swan 10 features
swan_rand_10_table <- getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_bh_10), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "swan_rand_10")

# swan 11 features
swan_rand_11_table <- getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_bh_11), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "swan_rand_11")

# swan 12 features
swan_rand_12_table <- getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_bh_12), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "swan_rand_12")

# swan 13 features
swan_rand_13_table <- getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_bh_13), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "swan_rand_13")

# swan 14 features
swan_rand_14_table <- getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_bh_14), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "swan_rand_14")

# swan 15 features
swan_rand_15_table <- getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_bh_15), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "swan_rand_15")

###########
# swan sig
###########

# swan 07 features
swan_rand_07_sig_table <- getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_bh_sig_07), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "swan_rand_sig_07")

# swan 08 features
swan_rand_08_sig_table <- getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_bh_sig_08), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "swan_rand_sig_08")

# swan 09 features
swan_rand_09_sig_table <- getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_bh_sig_09), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "swan_rand_sig_09")


# swan 10 features
swan_rand_10_sig_table <- getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_bh_sig_10), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "swan_rand_sig_10")

# swan 11 features
swan_rand_11_sig_table <- getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_bh_sig_11), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "swan_rand_sig_11")

# swan 12 features
swan_rand_12_sig_table <- getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_bh_sig_12), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "swan_rand_sig_12")

# swan 13 features
swan_rand_13_sig_table <- getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_bh_sig_13), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "swan_rand_sig_13")

# swan 14 features
swan_rand_14_sig_table <- getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_bh_sig_14), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "swan_rand_sig_14")

# swan 15 features
swan_rand_15_sig_table <- getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_bh_sig_15), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "swan_rand_sig_15")

##########
# intersection and union
##########

# swan int features
swan_rand_int_table <- getRand(swan_cases, 
                              rand_num = 5, 
                              num_feat = length(swan_int_feat), 
                              exclude = T,
                              union_features = union_features,
                              data_name = "swan_rand_int")

# swan int sig features
swan_rand_int_sig_table <- getRand(swan_cases, 
                                  rand_num = 5, 
                                  num_feat = length(swan_int_sig_feat), 
                                  exclude = T,
                                  union_features = union_features,
                                  data_name = "swan_rand_int_sig")

# swan union features
swan_rand_union_table <- getRand(swan_cases, 
                                rand_num = 5, 
                                num_feat = length(swan_union_feat), 
                                exclude = T,
                                union_features = union_features,
                                data_name = "swan_rand_union")

# swan union sig features
swan_rand_union_sig_table <- getRand(swan_cases, 
                                    rand_num = 5, 
                                    num_feat = length(swan_union_sig_feat), 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "swan_rand_union_sig")



###########
# rbind tables and save RDA file
###########
swan_rand_table <- rbind(swan_rand_07_table, swan_rand_07_sig_table,
                        swan_rand_08_table, swan_rand_08_sig_table,
                        swan_rand_09_table, swan_rand_09_sig_table,
                        swan_rand_10_table, swan_rand_10_sig_table,
                        swan_rand_11_table, swan_rand_11_sig_table,
                        swan_rand_12_table, swan_rand_12_sig_table,
                        swan_rand_13_table, swan_rand_13_sig_table,
                        swan_rand_14_table, swan_rand_14_sig_table,
                        swan_rand_15_table, swan_rand_15_sig_table,
                        swan_rand_int_table, swan_rand_int_sig_table,
                        swan_rand_union_table, swan_rand_union_sig_table)

# remove data 
rm(swan_rand_07_table, swan_rand_07_sig_table,
   swan_rand_08_table, swan_rand_08_sig_table,
   swan_rand_09_table, swan_rand_09_sig_table,
   swan_rand_10_table, swan_rand_10_sig_table,
   swan_rand_11_table, swan_rand_11_sig_table,
   swan_rand_12_table, swan_rand_12_sig_table,
   swan_rand_13_table, swan_rand_13_sig_table,
   swan_rand_14_table, swan_rand_14_sig_table,
   swan_rand_15_table, swan_rand_15_sig_table,
   swan_rand_int_table, swan_rand_int_sig_table,
   swan_rand_union_table, swan_rand_union_sig_table)


#save table 
saveRDS(swan_rand_table, 
        file = paste0(rand_folder, '/swan_rand_table.rda'))

rm(list = ls(pattern = "swan_*"))



#############################################################################################################
###########
# funnorm
###########

# funnorm 07 features
funnorm_rand_07_table <- getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_bh_07), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "funnorm_rand_07")

# funnorm 08 features
funnorm_rand_08_table <- getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_bh_08), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "funnorm_rand_08")

# funnorm 09 features
funnorm_rand_09_table <- getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_bh_09), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "funnorm_rand_09")


# funnorm 10 features
funnorm_rand_10_table <- getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_bh_10), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "funnorm_rand_10")

# funnorm 11 features
funnorm_rand_11_table <- getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_bh_11), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "funnorm_rand_11")

# funnorm 12 features
funnorm_rand_12_table <- getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_bh_12), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "funnorm_rand_12")

# funnorm 13 features
funnorm_rand_13_table <- getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_bh_13), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "funnorm_rand_13")

# funnorm 14 features
funnorm_rand_14_table <- getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_bh_14), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "funnorm_rand_14")

# funnorm 15 features
funnorm_rand_15_table <- getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_bh_15), 
                             exclude = T,
                             union_features = union_features,
                             data_name = "funnorm_rand_15")

###########
# funnorm sig
###########

# funnorm 07 features
funnorm_rand_07_sig_table <- getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_bh_sig_07), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "funnorm_rand_sig_07")

# funnorm 08 features
funnorm_rand_08_sig_table <- getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_bh_sig_08), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "funnorm_rand_sig_08")

# funnorm 09 features
funnorm_rand_09_sig_table <- getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_bh_sig_09), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "funnorm_rand_sig_09")


# funnorm 10 features
funnorm_rand_10_sig_table <- getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_bh_sig_10), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "funnorm_rand_sig_10")

# funnorm 11 features
funnorm_rand_11_sig_table <- getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_bh_sig_11), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "funnorm_rand_sig_11")

# funnorm 12 features
funnorm_rand_12_sig_table <- getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_bh_sig_12), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "funnorm_rand_sig_12")

# funnorm 13 features
funnorm_rand_13_sig_table <- getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_bh_sig_13), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "funnorm_rand_sig_13")

# funnorm 14 features
funnorm_rand_14_sig_table <- getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_bh_sig_14), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "funnorm_rand_sig_14")

# funnorm 15 features
funnorm_rand_15_sig_table <- getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_bh_sig_15), 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "funnorm_rand_sig_15")

##########
# intersection and union
##########

# funnorm int features
funnorm_rand_int_table <- getRand(funnorm_cases, 
                              rand_num = 5, 
                              num_feat = length(funnorm_int_feat), 
                              exclude = T,
                              union_features = union_features,
                              data_name = "funnorm_rand_int")

# funnorm int sig features
funnorm_rand_int_sig_table <- getRand(funnorm_cases, 
                                  rand_num = 5, 
                                  num_feat = length(funnorm_int_sig_feat), 
                                  exclude = T,
                                  union_features = union_features,
                                  data_name = "funnorm_rand_int_sig")

# funnorm union features
funnorm_rand_union_table <- getRand(funnorm_cases, 
                                rand_num = 5, 
                                num_feat = length(funnorm_union_feat), 
                                exclude = T,
                                union_features = union_features,
                                data_name = "funnorm_rand_union")

# funnorm union sig features
funnorm_rand_union_sig_table <- getRand(funnorm_cases, 
                                    rand_num = 5, 
                                    num_feat = length(funnorm_union_sig_feat), 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "funnorm_rand_union_sig")



###########
# rbind tables and save RDA file
###########
funnorm_rand_table <- rbind(funnorm_rand_07_table, funnorm_rand_07_sig_table,
                        funnorm_rand_08_table, funnorm_rand_08_sig_table,
                        funnorm_rand_09_table, funnorm_rand_09_sig_table,
                        funnorm_rand_10_table, funnorm_rand_10_sig_table,
                        funnorm_rand_11_table, funnorm_rand_11_sig_table,
                        funnorm_rand_12_table, funnorm_rand_12_sig_table,
                        funnorm_rand_13_table, funnorm_rand_13_sig_table,
                        funnorm_rand_14_table, funnorm_rand_14_sig_table,
                        funnorm_rand_15_table, funnorm_rand_15_sig_table,
                        funnorm_rand_int_table, funnorm_rand_int_sig_table,
                        funnorm_rand_union_table, funnorm_rand_union_sig_table)

# remove data 
rm(funnorm_rand_07_table, funnorm_rand_07_sig_table,
   funnorm_rand_08_table, funnorm_rand_08_sig_table,
   funnorm_rand_09_table, funnorm_rand_09_sig_table,
   funnorm_rand_10_table, funnorm_rand_10_sig_table,
   funnorm_rand_11_table, funnorm_rand_11_sig_table,
   funnorm_rand_12_table, funnorm_rand_12_sig_table,
   funnorm_rand_13_table, funnorm_rand_13_sig_table,
   funnorm_rand_14_table, funnorm_rand_14_sig_table,
   funnorm_rand_15_table, funnorm_rand_15_sig_table,
   funnorm_rand_int_table, funnorm_rand_int_sig_table,
   funnorm_rand_union_table, funnorm_rand_union_sig_table)


#save table 
saveRDS(funnorm_rand_table, 
        file = paste0(rand_folder, '/funnorm_rand_table.rda'))

rm(list = ls(pattern = "funnorm_*"))


