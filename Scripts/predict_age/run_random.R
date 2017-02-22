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
# Random features - same number as each bh feature subset

###########
# function for random
###########
# data <- beta_raw
# rand_num <- 2
# num_feat <- 100
# data_name <- 'raw_rand_100'
# exlcude <- T
# union_features <-  union_features
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

# NORMAL
# raw 10 features
raw_rand_10_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_10), 
                             exclude = T,
                             union_features = raw_union_feat,
                             data_name = "raw_rand_10")


# raw 15 features
raw_rand_15_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_15), 
                             exclude = T,
                             union_features = raw_union_feat,
                             data_name = "raw_rand_15")

# raw 20 features
raw_rand_20_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_20), 
                             exclude = T,
                             union_features = raw_union_feat,
                             data_name = "raw_rand_20")

# UNBALL
# raw 10 features
raw_rand_unbal_10_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_unbal_bh_10), 
                             exclude = T,
                             union_features = raw_union_feat,
                             data_name = "raw_rand_unbal_10")


# raw 15 features
raw_rand_unbal_15_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_unbal_bh_15), 
                             exclude = T,
                             union_features = raw_union_feat,
                             data_name = "raw_rand_unbal_15")

# raw 20 features
raw_rand_unbal_20_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_unbal_bh_20), 
                             exclude = T,
                             union_features = raw_union_feat,
                             data_name = "raw_rand_unbal_20")

raw_rand_int_table <-getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_int_feat), 
                             exclude = T,
                             union_features = raw_union_features,
                             data_name = "raw_rand_int")

raw_rand_unbal_int_table <-getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_unbal_int_feat), 
                             exclude = T,
                             union_features = raw_union_feat,
                             data_name = "raw_rand_unbal_int")


###########
# raw sig
###########
# NORMAL
# raw 10 features
raw_rand_sig_10_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_sig_10), 
                             exclude = T,
                             union_features = raw_union_feat,
                             data_name = "raw_rand_sig_10")


# raw 15 features
raw_rand_sig_15_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_sig_15), 
                             exclude = T,
                             union_features = raw_union_feat,
                             data_name = "raw_rand_sig_15")

# raw 20 features
raw_rand_sig_20_table <- getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_bh_sig_20), 
                             exclude = T,
                             union_features = raw_union_feat,
                             data_name = "raw_rand_sig_20")

# UNBALL
# raw 10 features
raw_rand_unbal_sig_10_table <- getRand(raw_cases, 
                                   rand_num = 5, 
                                   num_feat = length(raw_unbal_bh_sig_10), 
                                   exclude = T,
                                   union_features = raw_union_feat,
                                   data_name = "raw_rand_unbal_sig_10")


# raw 15 features
raw_rand_unbal_sig_15_table <- getRand(raw_cases, 
                                   rand_num = 5, 
                                   num_feat = length(raw_unbal_bh_sig_15), 
                                   exclude = T,
                                   union_features = raw_union_feat,
                                   data_name = "raw_rand_unbal_sig_15")

# raw 20 features
raw_rand_unbal_sig_20_table <- getRand(raw_cases, 
                                   rand_num = 5, 
                                   num_feat = length(raw_unbal_bh_20), 
                                   exclude = T,
                                   union_features = raw_union_feat,
                                   data_name = "raw_rand_unbal_20")

raw_rand_sig_int_table <-getRand(raw_cases, 
                             rand_num = 5, 
                             num_feat = length(raw_int_sig_feat), 
                             exclude = T,
                             union_features = raw_union_feat,
                             data_name = "raw_rand_sig_int")

raw_rand_unbal_sig_int_table <-getRand(raw_cases, 
                                   rand_num = 5, 
                                   num_feat = length(raw_unbal_int_sig_feat), 
                                   exclude = T,
                                   union_features = raw_union_feat,
                                   data_name = "raw_rand_unbal_int_sig")




###########
# rbind tables and save RDA file
###########
raw_rand_table <- rbind(raw_rand_10_table, raw_rand_15_table, 
                        raw_rand_20_table,  raw_rand_unbal_10_table, 
                        raw_rand_unbal_15_table, raw_rand_unbal_20_table,
                        raw_rand_unbal_int_table, raw_rand_int_table,
                        raw_rand_sig_10_table, raw_rand_sig_15_table, 
                        raw_rand_sig_20_table,  raw_rand_unbal_sig_10_table, 
                        raw_rand_unbal_sig_15_table, raw_rand_unbal_sig_20_table,
                        raw_rand_unbal_sig_int_table, raw_rand_sig_int_table)

# remove data 
rm(raw_rand_10_table, raw_rand_15_table, 
   raw_rand_20_table,  raw_rand_unbal_10_table, 
   raw_rand_unbal_15_table, raw_rand_unbal_20_table,
   raw_rand_unbal_int_table, raw_rand_int_table,
   raw_rand_sig_10_table, raw_rand_sig_15_table, 
   raw_rand_sig_20_table,  raw_rand_unbal_sig_10_table, 
   raw_rand_unbal_sig_15_table, raw_rand_unbal_sig_20_table,
   raw_rand_unbal_sig_int_table, raw_rand_sig_int_table)


#save table 
saveRDS(raw_rand_table, 
        file = paste0(rand_folder, '/raw_rand_table.rda'))

rm(list = ls(pattern = "raw_*"))


#############################################################################################################
###########
# quan
###########

# NORMAL
# quan 10 features
quan_rand_10_table <- getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_bh_10), 
                             exclude = T,
                             union_features = quan_union_feat,
                             data_name = "quan_rand_10")


# quan 15 features
quan_rand_15_table <- getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_bh_15), 
                             exclude = T,
                             union_features = quan_union_feat,
                             data_name = "quan_rand_15")

# quan 20 features
quan_rand_20_table <- getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_bh_20), 
                             exclude = T,
                             union_features = quan_union_feat,
                             data_name = "quan_rand_20")

# UNBALL
# quan 10 features
quan_rand_unbal_10_table <- getRand(quan_cases, 
                                   rand_num = 5, 
                                   num_feat = length(quan_unbal_bh_10), 
                                   exclude = T,
                                   union_features = quan_union_feat,
                                   data_name = "quan_rand_unbal_10")


# quan 15 features
quan_rand_unbal_15_table <- getRand(quan_cases, 
                                   rand_num = 5, 
                                   num_feat = length(quan_unbal_bh_15), 
                                   exclude = T,
                                   union_features = quan_union_feat,
                                   data_name = "quan_rand_unbal_15")

# quan 20 features
quan_rand_unbal_20_table <- getRand(quan_cases, 
                                   rand_num = 5, 
                                   num_feat = length(quan_unbal_bh_20), 
                                   exclude = T,
                                   union_features = quan_union_feat,
                                   data_name = "quan_rand_unbal_20")

quan_rand_int_table <-getRand(quan_cases, 
                             rand_num = 5, 
                             num_feat = length(quan_int_feat), 
                             exclude = T,
                             union_features = quan_union_features,
                             data_name = "quan_rand_int")

quan_rand_unbal_int_table <-getRand(quan_cases, 
                                   rand_num = 5, 
                                   num_feat = length(quan_unbal_int_feat), 
                                   exclude = T,
                                   union_features = quan_union_feat,
                                   data_name = "quan_rand_unbal_int")


###########
# quan sig
###########
# NORMAL
# quan 10 features
quan_rand_sig_10_table <- getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_bh_sig_10), 
                                 exclude = T,
                                 union_features = quan_union_feat,
                                 data_name = "quan_rand_sig_10")


# quan 15 features
quan_rand_sig_15_table <- getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_bh_sig_15), 
                                 exclude = T,
                                 union_features = quan_union_feat,
                                 data_name = "quan_rand_sig_15")

# quan 20 features
quan_rand_sig_20_table <- getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_bh_sig_20), 
                                 exclude = T,
                                 union_features = quan_union_feat,
                                 data_name = "quan_rand_sig_20")

# UNBALL
# quan 10 features
quan_rand_unbal_sig_10_table <- getRand(quan_cases, 
                                       rand_num = 5, 
                                       num_feat = length(quan_unbal_bh_sig_10), 
                                       exclude = T,
                                       union_features = quan_union_feat,
                                       data_name = "quan_rand_unbal_sig_10")


# quan 15 features
quan_rand_unbal_sig_15_table <- getRand(quan_cases, 
                                       rand_num = 5, 
                                       num_feat = length(quan_unbal_bh_sig_15), 
                                       exclude = T,
                                       union_features = quan_union_feat,
                                       data_name = "quan_rand_unbal_sig_15")

# quan 20 features
quan_rand_unbal_sig_20_table <- getRand(quan_cases, 
                                       rand_num = 5, 
                                       num_feat = length(quan_unbal_bh_20), 
                                       exclude = T,
                                       union_features = quan_union_feat,
                                       data_name = "quan_rand_unbal_20")

quan_rand_sig_int_table <-getRand(quan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(quan_int_sig_feat), 
                                 exclude = T,
                                 union_features = quan_union_feat,
                                 data_name = "quan_rand_sig_int")

quan_rand_unbal_sig_int_table <-getRand(quan_cases, 
                                       rand_num = 5, 
                                       num_feat = length(quan_unbal_int_sig_feat), 
                                       exclude = T,
                                       union_features = quan_union_feat,
                                       data_name = "quan_rand_unbal_int_sig")




###########
# rbind tables and save RDA file
###########
quan_rand_table <- rbind(quan_rand_10_table, quan_rand_15_table, 
                        quan_rand_20_table,  quan_rand_unbal_10_table, 
                        quan_rand_unbal_15_table, quan_rand_unbal_20_table,
                        quan_rand_unbal_int_table, quan_rand_int_table,
                        quan_rand_sig_10_table, quan_rand_sig_15_table, 
                        quan_rand_sig_20_table,  quan_rand_unbal_sig_10_table, 
                        quan_rand_unbal_sig_15_table, quan_rand_unbal_sig_20_table,
                        quan_rand_unbal_sig_int_table, quan_rand_sig_int_table)

# remove data 
rm(quan_rand_10_table, quan_rand_15_table, 
   quan_rand_20_table,  quan_rand_unbal_10_table, 
   quan_rand_unbal_15_table, quan_rand_unbal_20_table,
   quan_rand_unbal_int_table, quan_rand_int_table,
   quan_rand_sig_10_table, quan_rand_sig_15_table, 
   quan_rand_sig_20_table,  quan_rand_unbal_sig_10_table, 
   quan_rand_unbal_sig_15_table, quan_rand_unbal_sig_20_table,
   quan_rand_unbal_sig_int_table, quan_rand_sig_int_table)


#save table 
saveRDS(quan_rand_table, 
        file = paste0(rand_folder, '/quan_rand_table.rda'))

rm(list = ls(pattern = "quan_*"))


###########
# swan
###########

# NORMAL
# swan 10 features
swan_rand_10_table <- getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_bh_10), 
                             exclude = T,
                             union_features = swan_union_feat,
                             data_name = "swan_rand_10")


# swan 15 features
swan_rand_15_table <- getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_bh_15), 
                             exclude = T,
                             union_features = swan_union_feat,
                             data_name = "swan_rand_15")

# swan 20 features
swan_rand_20_table <- getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_bh_20), 
                             exclude = T,
                             union_features = swan_union_feat,
                             data_name = "swan_rand_20")

# UNBALL
# swan 10 features
swan_rand_unbal_10_table <- getRand(swan_cases, 
                                   rand_num = 5, 
                                   num_feat = length(swan_unbal_bh_10), 
                                   exclude = T,
                                   union_features = swan_union_feat,
                                   data_name = "swan_rand_unbal_10")


# swan 15 features
swan_rand_unbal_15_table <- getRand(swan_cases, 
                                   rand_num = 5, 
                                   num_feat = length(swan_unbal_bh_15), 
                                   exclude = T,
                                   union_features = swan_union_feat,
                                   data_name = "swan_rand_unbal_15")

# swan 20 features
swan_rand_unbal_20_table <- getRand(swan_cases, 
                                   rand_num = 5, 
                                   num_feat = length(swan_unbal_bh_20), 
                                   exclude = T,
                                   union_features = swan_union_feat,
                                   data_name = "swan_rand_unbal_20")

swan_rand_int_table <-getRand(swan_cases, 
                             rand_num = 5, 
                             num_feat = length(swan_int_feat), 
                             exclude = T,
                             union_features = swan_union_features,
                             data_name = "swan_rand_int")

swan_rand_unbal_int_table <-getRand(swan_cases, 
                                   rand_num = 5, 
                                   num_feat = length(swan_unbal_int_feat), 
                                   exclude = T,
                                   union_features = swan_union_feat,
                                   data_name = "swan_rand_unbal_int")


###########
# swan sig
###########
# NORMAL
# swan 10 features
swan_rand_sig_10_table <- getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_bh_sig_10), 
                                 exclude = T,
                                 union_features = swan_union_feat,
                                 data_name = "swan_rand_sig_10")


# swan 15 features
swan_rand_sig_15_table <- getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_bh_sig_15), 
                                 exclude = T,
                                 union_features = swan_union_feat,
                                 data_name = "swan_rand_sig_15")

# swan 20 features
swan_rand_sig_20_table <- getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_bh_sig_20), 
                                 exclude = T,
                                 union_features = swan_union_feat,
                                 data_name = "swan_rand_sig_20")

# UNBALL
# swan 10 features
swan_rand_unbal_sig_10_table <- getRand(swan_cases, 
                                       rand_num = 5, 
                                       num_feat = length(swan_unbal_bh_sig_10), 
                                       exclude = T,
                                       union_features = swan_union_feat,
                                       data_name = "swan_rand_unbal_sig_10")


# swan 15 features
swan_rand_unbal_sig_15_table <- getRand(swan_cases, 
                                       rand_num = 5, 
                                       num_feat = length(swan_unbal_bh_sig_15), 
                                       exclude = T,
                                       union_features = swan_union_feat,
                                       data_name = "swan_rand_unbal_sig_15")

# swan 20 features
swan_rand_unbal_sig_20_table <- getRand(swan_cases, 
                                       rand_num = 5, 
                                       num_feat = length(swan_unbal_bh_20), 
                                       exclude = T,
                                       union_features = swan_union_feat,
                                       data_name = "swan_rand_unbal_20")

swan_rand_sig_int_table <-getRand(swan_cases, 
                                 rand_num = 5, 
                                 num_feat = length(swan_int_sig_feat), 
                                 exclude = T,
                                 union_features = swan_union_feat,
                                 data_name = "swan_rand_sig_int")

swan_rand_unbal_sig_int_table <-getRand(swan_cases, 
                                       rand_num = 5, 
                                       num_feat = length(swan_unbal_int_sig_feat), 
                                       exclude = T,
                                       union_features = swan_union_feat,
                                       data_name = "swan_rand_unbal_int_sig")




###########
# rbind tables and save RDA file
###########
swan_rand_table <- rbind(swan_rand_10_table, swan_rand_15_table, 
                        swan_rand_20_table,  swan_rand_unbal_10_table, 
                        swan_rand_unbal_15_table, swan_rand_unbal_20_table,
                        swan_rand_unbal_int_table, swan_rand_int_table,
                        swan_rand_sig_10_table, swan_rand_sig_15_table, 
                        swan_rand_sig_20_table,  swan_rand_unbal_sig_10_table, 
                        swan_rand_unbal_sig_15_table, swan_rand_unbal_sig_20_table,
                        swan_rand_unbal_sig_int_table, swan_rand_sig_int_table)

# remove data 
rm(swan_rand_10_table, swan_rand_15_table, 
   swan_rand_20_table,  swan_rand_unbal_10_table, 
   swan_rand_unbal_15_table, swan_rand_unbal_20_table,
   swan_rand_unbal_int_table, swan_rand_int_table,
   swan_rand_sig_10_table, swan_rand_sig_15_table, 
   swan_rand_sig_20_table,  swan_rand_unbal_sig_10_table, 
   swan_rand_unbal_sig_15_table, swan_rand_unbal_sig_20_table,
   swan_rand_unbal_sig_int_table, swan_rand_sig_int_table)


#save table 
saveRDS(swan_rand_table, 
        file = paste0(rand_folder, '/swan_rand_table.rda'))

rm(list = ls(pattern = "swan_*"))

###########
# funnorm
###########

# NORMAL
# funnorm 10 features
funnorm_rand_10_table <- getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_bh_10), 
                             exclude = T,
                             union_features = funnorm_union_feat,
                             data_name = "funnorm_rand_10")


# funnorm 15 features
funnorm_rand_15_table <- getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_bh_15), 
                             exclude = T,
                             union_features = funnorm_union_feat,
                             data_name = "funnorm_rand_15")

# funnorm 20 features
funnorm_rand_20_table <- getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_bh_20), 
                             exclude = T,
                             union_features = funnorm_union_feat,
                             data_name = "funnorm_rand_20")

# UNBALL
# funnorm 10 features
funnorm_rand_unbal_10_table <- getRand(funnorm_cases, 
                                   rand_num = 5, 
                                   num_feat = length(funnorm_unbal_bh_10), 
                                   exclude = T,
                                   union_features = funnorm_union_feat,
                                   data_name = "funnorm_rand_unbal_10")


# funnorm 15 features
funnorm_rand_unbal_15_table <- getRand(funnorm_cases, 
                                   rand_num = 5, 
                                   num_feat = length(funnorm_unbal_bh_15), 
                                   exclude = T,
                                   union_features = funnorm_union_feat,
                                   data_name = "funnorm_rand_unbal_15")

# funnorm 20 features
funnorm_rand_unbal_20_table <- getRand(funnorm_cases, 
                                   rand_num = 5, 
                                   num_feat = length(funnorm_unbal_bh_20), 
                                   exclude = T,
                                   union_features = funnorm_union_feat,
                                   data_name = "funnorm_rand_unbal_20")

funnorm_rand_int_table <-getRand(funnorm_cases, 
                             rand_num = 5, 
                             num_feat = length(funnorm_int_feat), 
                             exclude = T,
                             union_features = funnorm_union_features,
                             data_name = "funnorm_rand_int")

funnorm_rand_unbal_int_table <-getRand(funnorm_cases, 
                                   rand_num = 5, 
                                   num_feat = length(funnorm_unbal_int_feat), 
                                   exclude = T,
                                   union_features = funnorm_union_feat,
                                   data_name = "funnorm_rand_unbal_int")


###########
# funnorm sig
###########
# NORMAL
# funnorm 10 features
funnorm_rand_sig_10_table <- getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_bh_sig_10), 
                                 exclude = T,
                                 union_features = funnorm_union_feat,
                                 data_name = "funnorm_rand_sig_10")


# funnorm 15 features
funnorm_rand_sig_15_table <- getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_bh_sig_15), 
                                 exclude = T,
                                 union_features = funnorm_union_feat,
                                 data_name = "funnorm_rand_sig_15")

# funnorm 20 features
funnorm_rand_sig_20_table <- getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_bh_sig_20), 
                                 exclude = T,
                                 union_features = funnorm_union_feat,
                                 data_name = "funnorm_rand_sig_20")

# UNBALL
# funnorm 10 features
funnorm_rand_unbal_sig_10_table <- getRand(funnorm_cases, 
                                       rand_num = 5, 
                                       num_feat = length(funnorm_unbal_bh_sig_10), 
                                       exclude = T,
                                       union_features = funnorm_union_feat,
                                       data_name = "funnorm_rand_unbal_sig_10")


# funnorm 15 features
funnorm_rand_unbal_sig_15_table <- getRand(funnorm_cases, 
                                       rand_num = 5, 
                                       num_feat = length(funnorm_unbal_bh_sig_15), 
                                       exclude = T,
                                       union_features = funnorm_union_feat,
                                       data_name = "funnorm_rand_unbal_sig_15")

# funnorm 20 features
funnorm_rand_unbal_sig_20_table <- getRand(funnorm_cases, 
                                       rand_num = 5, 
                                       num_feat = length(funnorm_unbal_bh_20), 
                                       exclude = T,
                                       union_features = funnorm_union_feat,
                                       data_name = "funnorm_rand_unbal_20")

funnorm_rand_sig_int_table <-getRand(funnorm_cases, 
                                 rand_num = 5, 
                                 num_feat = length(funnorm_int_sig_feat), 
                                 exclude = T,
                                 union_features = funnorm_union_feat,
                                 data_name = "funnorm_rand_sig_int")

funnorm_rand_unbal_sig_int_table <-getRand(funnorm_cases, 
                                       rand_num = 5, 
                                       num_feat = length(funnorm_unbal_int_sig_feat), 
                                       exclude = T,
                                       union_features = funnorm_union_feat,
                                       data_name = "funnorm_rand_unbal_int_sig")




###########
# rbind tables and save RDA file
###########
funnorm_rand_table <- rbind(funnorm_rand_10_table, funnorm_rand_15_table, 
                        funnorm_rand_20_table,  funnorm_rand_unbal_10_table, 
                        funnorm_rand_unbal_15_table, funnorm_rand_unbal_20_table,
                        funnorm_rand_unbal_int_table, funnorm_rand_int_table,
                        funnorm_rand_sig_10_table, funnorm_rand_sig_15_table, 
                        funnorm_rand_sig_20_table,  funnorm_rand_unbal_sig_10_table, 
                        funnorm_rand_unbal_sig_15_table, funnorm_rand_unbal_sig_20_table,
                        funnorm_rand_unbal_sig_int_table, funnorm_rand_sig_int_table)

# remove data 
rm(funnorm_rand_10_table, funnorm_rand_15_table, 
   funnorm_rand_20_table,  funnorm_rand_unbal_10_table, 
   funnorm_rand_unbal_15_table, funnorm_rand_unbal_20_table,
   funnorm_rand_unbal_int_table, funnorm_rand_int_table,
   funnorm_rand_sig_10_table, funnorm_rand_sig_15_table, 
   funnorm_rand_sig_20_table,  funnorm_rand_unbal_sig_10_table, 
   funnorm_rand_unbal_sig_15_table, funnorm_rand_unbal_sig_20_table,
   funnorm_rand_unbal_sig_int_table, funnorm_rand_sig_int_table)


#save table 
saveRDS(funnorm_rand_table, 
        file = paste0(rand_folder, '/funnorm_rand_table.rda'))

rm(list = ls(pattern = "funnorm_*"))


