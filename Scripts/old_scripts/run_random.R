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

##### REmove one type of data
rm(list=ls(pattern="raw"))

#########################################################################################################################
# Random features - same number as each bh feature subset

###########
# function for random
###########
# data <- quan_cases
# rand_num <- 2
# num_feat <- 10
# data_name <- 'raw_rand_10'
# exlcude <- T
# union_features <-  union_features
getRand <- function(data, model, rand_num, num_feat, exclude, union_features, data_name) 
{
  
  if (exclude) {
    
    features <- colnames(data[, 8:ncol(data)])
    keep_these <- features[!features %in% union_features]
    data <- data[, c('id', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                     'age_sample_collection', keep_these)]
  }
  
  
  if (model == 'rf') {
    
    rf_holder = list()
    rf_table = list()
    
    for (rand_num in 1:rand_num) {
      
      rf_holder[[rand_num]] <- runModels(data, 
                                         model = 'rf',
                                         bump_hunter = F,
                                         num_feat = num_feat,
                                         random = T,
                                         seed_num = rand_num)
      
      rf_table[[rand_num]] <- extractResults(rf_holder[[rand_num]],
                                             data_name = data_name,
                                             regularize = F)
      
    }
    
    full_table <- do.call(rbind, rf_table)
  }
  
  if (model == 'enet') {
    
    enet_holder = list()
    enet_table = list()
    
    for (rand_num in 1:rand_num) {
      
      enet_holder[[rand_num]] <- runModels(data, 
                                           model = 'enet',
                                           bump_hunter = F,
                                           num_feat = num_feat,
                                           random = T,
                                           seed_num = rand_num)
      
      enet_table[[rand_num]] <- extractResults(enet_holder[[rand_num]],
                                               data_name = data_name,
                                               regularize = F)
      
    }
    
    full_table <- do.call(rbind, enet_table)
  }
  
  if (model == 'lasso') {
    
    lasso_holder = list()
    lasso_table = list()
    
    for (rand_num in 1:rand_num) {
      
      lasso_holder[[rand_num]] <- runModels(data, 
                                            model = 'lasso',
                                            bump_hunter = F,
                                            num_feat = num_feat,
                                            random = T,
                                            seed_num = rand_num)
      
      lasso_table[[rand_num]] <- extractResults(lasso_holder[[rand_num]],
                                                data_name = data_name,
                                                regularize = F)
      
    }
    
    full_table <- do.call(rbind, lasso_table)
  }
  
  
  return(full_table)
}

##########
# get union 
##########

# get list of features
df_names <- ls()[sapply(mget(ls(), .GlobalEnv), is.data.frame)]
dfs <-  sapply( df_names, function(x)  get( x )  )
dfs[grepl('cases|controls', names(dfs))] <- NULL

# combine and remove duplicates
feats <- do.call(rbind, dfs)

# union 
union_features <- as.character(feats[!duplicated(feats),])

#############################################################################################################

###########
# rf
###########

# quan 10 features
quan_unbatch_rand_rf_10 <- getRand(quan_cases, 
                                 model = 'rf',
                                 rand_num = 10, 
                                 num_feat = 10, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_unbatch_rand_rf_10")

# quan 20 features
quan_unbatch_rand_rf_20 <- getRand(quan_cases, 
                                 model = 'rf',
                                 rand_num = 10, 
                                 num_feat = 20, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_unbatch_rand_rf_20")


# quan 30 features
quan_unbatch_rand_rf_30 <- getRand(quan_cases, 
                                 model = 'rf',
                                 rand_num = 10, 
                                 num_feat = 30, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_unbatch_rand_rf_30")

# quan 40 features
quan_unbatch_rand_rf_40 <- getRand(quan_cases, 
                                 model = 'rf',
                                 rand_num = 10, 
                                 num_feat = 40, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_unbatch_rand_rf_40")

# quan 50 features
quan_unbatch_rand_rf_50 <- getRand(quan_cases, 
                                 model = 'rf',
                                 rand_num = 10, 
                                 num_feat = 50, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_unbatch_rand_rf_50")

# quan 60 features
quan_unbatch_rand_rf_60 <- getRand(quan_cases, 
                                 model = 'rf',
                                 rand_num = 10, 
                                 num_feat = 60, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_unbatch_rand_60")

# quan 70 features
quan_unbatch_rand_rf_70 <- getRand(quan_cases, 
                                 model = 'rf',
                                 rand_num = 10, 
                                 num_feat = 70, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_unbatch_rand_rf_70")

# quan 80 features
quan_unbatch_rand_rf_80 <- getRand(quan_cases, 
                                 model = 'rf',
                                 rand_num = 10, 
                                 num_feat = 80, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_unbatch_rand_rf_80")

# quan 90 features
quan_unbatch_rand_rf_90 <- getRand(quan_cases, 
                                 model = 'rf',
                                 rand_num = 10, 
                                 num_feat = 90, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_unbatch_rand_rf_90")

# quan 100 features
quan_unbatch_rand_rf_100 <- getRand(quan_cases, 
                                  model = 'rf',
                                  rand_num = 10, 
                                  num_feat = 100, 
                                  exclude = T,
                                  union_features = union_features,
                                  data_name = "quan_unbatch_rand_rf_100")

# quan 500 features
quan_unbatch_rand_rf_500 <- getRand(quan_cases, 
                                  model = 'rf',
                                  rand_num = 10, 
                                  num_feat = 500, 
                                  exclude = T,
                                  union_features = union_features,
                                  data_name = "quan_unbatch_rand_rf_500")

# quan 1000 features
quan_unbatch_rand_rf_1000 <- getRand(quan_cases, 
                                   model = 'rf',
                                   rand_num = 10, 
                                   num_feat = 1000, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_unbatch_rand_rf_1000")

# quan 5000 features
quan_unbatch_rand_rf_5000 <- getRand(quan_cases, 
                                   model = 'rf',
                                   rand_num = 10, 
                                   num_feat = 5000, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_unbatch_rand_rf_5000")

# quan 1000 features
quan_unbatch_rand_rf_10000 <- getRand(quan_cases, 
                                    model = 'rf',
                                    rand_num = 10, 
                                    num_feat = 10000, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_unbatch_rand_rf_10000")



###########
# rbind tables and save RDA file
###########
quan_unbatch_rand_rf_table <- rbind(quan_unbatch_rand_rf_10, quan_unbatch_rand_rf_20, 
                                  quan_unbatch_rand_rf_30, quan_unbatch_rand_rf_40,
                                  quan_unbatch_rand_rf_50, quan_unbatch_rand_rf_60,
                                  quan_unbatch_rand_rf_70, quan_unbatch_rand_rf_80,
                                  quan_unbatch_rand_rf_90, quan_unbatch_rand_rf_100,
                                  quan_unbatch_rand_rf_500, quan_unbatch_rand_rf_1000,
                                  quan_unbatch_rand_rf_5000, quan_unbatch_rand_rf_10000)

# remove data 
rm(quan_unbatch_rand_rf_10, quan_unbatch_rand_rf_20, 
   quan_unbatch_rand_rf_30, quan_unbatch_rand_rf_40,
   quan_unbatch_rand_rf_50, quan_unbatch_rand_rf_60,
   quan_unbatch_rand_rf_70, quan_unbatch_rand_rf_80,
   quan_unbatch_rand_rf_90, quan_unbatch_rand_rf_100,
   quan_unbatch_rand_rf_500, quan_unbatch_rand_rf_1000,
   quan_unbatch_rand_rf_5000, quan_unbatch_rand_rf_10000)


#save table 
saveRDS(quan_unbatch_rand_rf_table, 
        file = paste0(rand_folder, '/quan_unbatch_rand_rf_table.rda'))


#######################################################################################################

###########
# enet
###########

# quan 10 features
quan_unbatch_rand_enet_10 <- getRand(quan_cases, 
                                   model = 'enet',
                                   rand_num = 10, 
                                   num_feat = 10, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_unbatch_rand_enet_10")

# quan 20 features
quan_unbatch_rand_enet_20 <- getRand(quan_cases, 
                                   model = 'enet',
                                   rand_num = 10, 
                                   num_feat = 20, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_unbatch_rand_enet_20")


# quan 30 features
quan_unbatch_rand_enet_30 <- getRand(quan_cases, 
                                   model = 'enet',
                                   rand_num = 10, 
                                   num_feat = 30, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_unbatch_rand_enet_30")

# quan 40 features
quan_unbatch_rand_enet_40 <- getRand(quan_cases, 
                                   model = 'enet',
                                   rand_num = 10, 
                                   num_feat = 40, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_unbatch_rand_enet_40")

# quan 50 features
quan_unbatch_rand_enet_50 <- getRand(quan_cases, 
                                   model = 'enet',
                                   rand_num = 10, 
                                   num_feat = 50, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_unbatch_rand_enet_50")

# quan 60 features
quan_unbatch_rand_enet_60 <- getRand(quan_cases, 
                                   model = 'enet',
                                   rand_num = 10, 
                                   num_feat = 60, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_unbatch_rand_enet_60")

# quan 70 features
quan_unbatch_rand_enet_70 <- getRand(quan_cases, 
                                   model = 'enet',
                                   rand_num = 10, 
                                   num_feat = 70, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_unbatch_rand_enet_70")

# quan 80 features
quan_unbatch_rand_enet_80 <- getRand(quan_cases, 
                                   model = 'enet',
                                   rand_num = 10, 
                                   num_feat = 80, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_unbatch_rand_enet_80")

# quan 90 features
quan_unbatch_rand_enet_90 <- getRand(quan_cases, 
                                   model = 'enet',
                                   rand_num = 10, 
                                   num_feat = 90, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_unbatch_rand_enet_90")

# quan 100 features
quan_unbatch_rand_enet_100 <- getRand(quan_cases, 
                                    model = 'enet',
                                    rand_num = 10, 
                                    num_feat = 100, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_unbatch_rand_enet_100")

# quan 500 features
quan_unbatch_rand_enet_500 <- getRand(quan_cases, 
                                    model = 'enet',
                                    rand_num = 10, 
                                    num_feat = 500, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_unbatch_rand_enet_500")

# quan 1000 features
quan_unbatch_rand_enet_1000 <- getRand(quan_cases, 
                                     model = 'enet',
                                     rand_num = 10, 
                                     num_feat = 1000, 
                                     exclude = T,
                                     union_features = union_features,
                                     data_name = "quan_unbatch_rand_enet_1000")

# quan 5000 features
quan_unbatch_rand_enet_5000 <- getRand(quan_cases, 
                                     model = 'enet',
                                     rand_num = 10, 
                                     num_feat = 5000, 
                                     exclude = T,
                                     union_features = union_features,
                                     data_name = "quan_unbatch_rand_enet_5000")

# quan 1000 features
quan_unbatch_rand_enet_10000 <- getRand(quan_cases, 
                                      model = 'enet',
                                      rand_num = 10, 
                                      num_feat = 10000, 
                                      exclude = T,
                                      union_features = union_features,
                                      data_name = "quan_unbatch_rand_enet_10000")



###########
# rbind tables and save RDA file
###########
quan_unbatch_rand_enet_table <- rbind(quan_unbatch_rand_enet_10, quan_unbatch_rand_enet_20, 
                                    quan_unbatch_rand_enet_30, quan_unbatch_rand_enet_40,
                                    quan_unbatch_rand_enet_50, quan_unbatch_rand_enet_60,
                                    quan_unbatch_rand_enet_70, quan_unbatch_rand_enet_80,
                                    quan_unbatch_rand_enet_90, quan_unbatch_rand_enet_100,
                                    quan_unbatch_rand_enet_500, quan_unbatch_rand_enet_1000,
                                    quan_unbatch_rand_enet_5000, quan_unbatch_rand_enet_10000)

# remove data 
rm(quan_unbatch_rand_enet_10, quan_unbatch_rand_enet_20, 
   quan_unbatch_rand_enet_30, quan_unbatch_rand_enet_40,
   quan_unbatch_rand_enet_50, quan_unbatch_rand_enet_60,
   quan_unbatch_rand_enet_70, quan_unbatch_rand_enet_80,
   quan_unbatch_rand_enet_90, quan_unbatch_rand_enet_100,
   quan_unbatch_rand_enet_500, quan_unbatch_rand_enet_1000,
   quan_unbatch_rand_enet_5000, quan_unbatch_rand_enet_10000)


#save table 
saveRDS(quan_unbatch_rand_enet_table, 
        file = paste0(rand_folder, '/quan_unbatch_rand_enet_table.rda'))



###########
# lasso
###########

# quan 10 features
quan_unbatch_rand_lasso_10 <- getRand(quan_cases, 
                                    model = 'lasso',
                                    rand_num = 10, 
                                    num_feat = 10, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_unbatch_rand_lasso_10")

# quan 20 features
quan_unbatch_rand_lasso_20 <- getRand(quan_cases, 
                                    model = 'lasso',
                                    rand_num = 10, 
                                    num_feat = 20, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_unbatch_rand_lasso_20")


# quan 30 features
quan_unbatch_rand_lasso_30 <- getRand(quan_cases, 
                                    model = 'lasso',
                                    rand_num = 10, 
                                    num_feat = 30, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_unbatch_rand_lasso_30")

# quan 40 features
quan_unbatch_rand_lasso_40 <- getRand(quan_cases, 
                                    model = 'lasso',
                                    rand_num = 10, 
                                    num_feat = 40, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_unbatch_rand_lasso_40")

# quan 50 features
quan_unbatch_rand_lasso_50 <- getRand(quan_cases, 
                                    model = 'lasso',
                                    rand_num = 10, 
                                    num_feat = 50, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_unbatch_rand_lasso_50")

# quan 60 features
quan_unbatch_rand_lasso_60 <- getRand(quan_cases, 
                                    model = 'lasso',
                                    rand_num = 10, 
                                    num_feat = 60, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_unbatch_rand_lasso_60")

# quan 70 features
quan_unbatch_rand_lasso_70 <- getRand(quan_cases, 
                                    model = 'lasso',
                                    rand_num = 10, 
                                    num_feat = 70, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_unbatch_rand_lasso_70")

# quan 80 features
quan_unbatch_rand_lasso_80 <- getRand(quan_cases, 
                                    model = 'lasso',
                                    rand_num = 10, 
                                    num_feat = 80, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_unbatch_rand_lasso_80")

# quan 90 features
quan_unbatch_rand_lasso_90 <- getRand(quan_cases, 
                                    model = 'lasso',
                                    rand_num = 10, 
                                    num_feat = 90, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_unbatch_rand_lasso_90")

# quan 100 features
quan_unbatch_rand_lasso_100 <- getRand(quan_cases, 
                                     model = 'lasso',
                                     rand_num = 10, 
                                     num_feat = 100, 
                                     exclude = T,
                                     union_features = union_features,
                                     data_name = "quan_unbatch_rand_lasso_100")

# quan 500 features
quan_unbatch_rand_lasso_500 <- getRand(quan_cases, 
                                     model = 'lasso',
                                     rand_num = 10, 
                                     num_feat = 500, 
                                     exclude = T,
                                     union_features = union_features,
                                     data_name = "quan_unbatch_rand_lasso_500")

# quan 1000 features
quan_unbatch_rand_lasso_1000 <- getRand(quan_cases, 
                                      model = 'lasso',
                                      rand_num = 10, 
                                      num_feat = 1000, 
                                      exclude = T,
                                      union_features = union_features,
                                      data_name = "quan_unbatch_rand_lasso_1000")

# quan 5000 features
quan_unbatch_rand_lasso_5000 <- getRand(quan_cases, 
                                      model = 'lasso',
                                      rand_num = 10, 
                                      num_feat = 5000, 
                                      exclude = T,
                                      union_features = union_features,
                                      data_name = "quan_unbatch_rand_lasso_5000")

# quan 1000 features
quan_unbatch_rand_lasso_10000 <- getRand(quan_cases, 
                                       model = 'lasso',
                                       rand_num = 10, 
                                       num_feat = 10000, 
                                       exclude = T,
                                       union_features = union_features,
                                       data_name = "quan_unbatch_rand_lasso_10000")



###########
# rbind tables and save RDA file
###########
quan_unbatch_rand_lasso_table <- rbind(quan_unbatch_rand_lasso_10, quan_unbatch_rand_lasso_20, 
                                     quan_unbatch_rand_lasso_30, quan_unbatch_rand_lasso_40,
                                     quan_unbatch_rand_lasso_50, quan_unbatch_rand_lasso_60,
                                     quan_unbatch_rand_lasso_70, quan_unbatch_rand_lasso_80,
                                     quan_unbatch_rand_lasso_90, quan_unbatch_rand_lasso_100,
                                     quan_unbatch_rand_lasso_500, quan_unbatch_rand_lasso_1000,
                                     quan_unbatch_rand_lasso_5000, quan_unbatch_rand_lasso_10000)

# remove data 
rm(quan_unbatch_rand_lasso_10, quan_unbatch_rand_lasso_20, 
   quan_unbatch_rand_lasso_30, quan_unbatch_rand_lasso_40,
   quan_unbatch_rand_lasso_50, quan_unbatch_rand_lasso_60,
   quan_unbatch_rand_lasso_70, quan_unbatch_rand_lasso_80,
   quan_unbatch_rand_lasso_90, quan_unbatch_rand_lasso_100,
   quan_unbatch_rand_lasso_500, quan_unbatch_rand_lasso_1000,
   quan_unbatch_rand_lasso_5000, quan_unbatch_rand_lasso_10000)


#save table 
saveRDS(quan_unbatch_rand_lasso_table, 
        file = paste0(rand_folder, '/quan_unbatch_rand_lasso_table.rda'))

#############################################################################################################


###########
# rf
###########

# quan 10 features
quan_batch_rand_rf_10 <- getRand(quan_cases_batch, 
                           model = 'rf',
                           rand_num = 10, 
                           num_feat = 10, 
                           exclude = T,
                           union_features = union_features,
                           data_name = "quan_batch_rand_rf_10")

# quan 20 features
quan_batch_rand_rf_20 <- getRand(quan_cases_batch, 
                           model = 'rf',
                           rand_num = 10, 
                           num_feat = 20, 
                           exclude = T,
                           union_features = union_features,
                           data_name = "quan_batch_rand_rf_20")


# quan 30 features
quan_batch_rand_rf_30 <- getRand(quan_cases_batch, 
                           model = 'rf',
                           rand_num = 10, 
                           num_feat = 30, 
                           exclude = T,
                           union_features = union_features,
                           data_name = "quan_batch_rand_rf_30")

# quan 40 features
quan_batch_rand_rf_40 <- getRand(quan_cases_batch, 
                           model = 'rf',
                           rand_num = 10, 
                           num_feat = 40, 
                           exclude = T,
                           union_features = union_features,
                           data_name = "quan_batch_rand_rf_40")

# quan 50 features
quan_batch_rand_rf_50 <- getRand(quan_cases_batch, 
                           model = 'rf',
                           rand_num = 10, 
                           num_feat = 50, 
                           exclude = T,
                           union_features = union_features,
                           data_name = "quan_batch_rand_rf_50")

# quan 60 features
quan_batch_rand_rf_60 <- getRand(quan_cases_batch, 
                           model = 'rf',
                           rand_num = 10, 
                           num_feat = 60, 
                           exclude = T,
                           union_features = union_features,
                           data_name = "quan_batch_rand_rf_60")

# quan 70 features
quan_batch_rand_rf_70 <- getRand(quan_cases_batch, 
                           model = 'rf',
                           rand_num = 10, 
                           num_feat = 70, 
                           exclude = T,
                           union_features = union_features,
                           data_name = "quan_batch_rand_rf_70")

# quan 80 features
quan_batch_rand_rf_80 <- getRand(quan_cases_batch, 
                           model = 'rf',
                           rand_num = 10, 
                           num_feat = 80, 
                           exclude = T,
                           union_features = union_features,
                           data_name = "quan_batch_rand_rf_80")

# quan 90 features
quan_batch_rand_rf_90 <- getRand(quan_cases_batch, 
                           model = 'rf',
                           rand_num = 10, 
                           num_feat = 90, 
                           exclude = T,
                           union_features = union_features,
                           data_name = "quan_batch_rand_rf_90")

# quan 100 features
quan_batch_rand_rf_100 <- getRand(quan_cases_batch, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 100, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "quan_batch_rand_rf_100")

# quan 500 features
quan_batch_rand_rf_500 <- getRand(quan_cases_batch, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 500, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "quan_batch_rand_rf_500")

# quan 1000 features
quan_batch_rand_rf_1000 <- getRand(quan_cases_batch, 
                             model = 'rf',
                             rand_num = 10, 
                             num_feat = 1000, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "quan_batch_rand_rf_1000")

# quan 5000 features
quan_batch_rand_rf_5000 <- getRand(quan_cases_batch, 
                             model = 'rf',
                             rand_num = 10, 
                             num_feat = 5000, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "quan_batch_rand_rf_5000")

# quan 1000 features
quan_batch_rand_rf_10000 <- getRand(quan_cases_batch, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 10000, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "quan_batch_rand_rf_10000")



###########
# rbind tables and save RDA file
###########
quan_batch_rand_rf_table <- rbind(quan_batch_rand_rf_10, quan_batch_rand_rf_20, 
                            quan_batch_rand_rf_30, quan_batch_rand_rf_40,
                            quan_batch_rand_rf_50, quan_batch_rand_rf_60,
                            quan_batch_rand_rf_70, quan_batch_rand_rf_80,
                            quan_batch_rand_rf_90, quan_batch_rand_rf_100,
                            quan_batch_rand_rf_500, quan_batch_rand_rf_1000,
                            quan_batch_rand_rf_5000, quan_batch_rand_rf_10000)

# remove data 
rm(quan_batch_rand_rf_10, quan_batch_rand_rf_20, 
   quan_batch_rand_rf_30, quan_batch_rand_rf_40,
   quan_batch_rand_rf_50, quan_batch_rand_rf_60,
   quan_batch_rand_rf_70, quan_batch_rand_rf_80,
   quan_batch_rand_rf_90, quan_batch_rand_rf_100,
   quan_batch_rand_rf_500, quan_batch_rand_rf_1000,
   quan_batch_rand_rf_5000, quan_batch_rand_rf_10000)


#save table 
saveRDS(quan_batch_rand_rf_table, 
        file = paste0(rand_folder, '/quan_batch_rand_rf_table.rda'))


#######################################################################################################

###########
# enet
###########

# quan 10 features
quan_batch_rand_enet_10 <- getRand(quan_cases_batch, 
                                 model = 'enet',
                                 rand_num = 10, 
                                 num_feat = 10, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_enet_10")

# quan 20 features
quan_batch_rand_enet_20 <- getRand(quan_cases_batch, 
                                 model = 'enet',
                                 rand_num = 10, 
                                 num_feat = 20, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_enet_20")


# quan 30 features
quan_batch_rand_enet_30 <- getRand(quan_cases_batch, 
                                 model = 'enet',
                                 rand_num = 10, 
                                 num_feat = 30, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_enet_30")

# quan 40 features
quan_batch_rand_enet_40 <- getRand(quan_cases_batch, 
                                 model = 'enet',
                                 rand_num = 10, 
                                 num_feat = 40, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_enet_40")

# quan 50 features
quan_batch_rand_enet_50 <- getRand(quan_cases_batch, 
                                 model = 'enet',
                                 rand_num = 10, 
                                 num_feat = 50, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_enet_50")

# quan 60 features
quan_batch_rand_enet_60 <- getRand(quan_cases_batch, 
                                 model = 'enet',
                                 rand_num = 10, 
                                 num_feat = 60, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_enet_60")

# quan 70 features
quan_batch_rand_enet_70 <- getRand(quan_cases_batch, 
                                 model = 'enet',
                                 rand_num = 10, 
                                 num_feat = 70, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_enet_70")

# quan 80 features
quan_batch_rand_enet_80 <- getRand(quan_cases_batch, 
                                 model = 'enet',
                                 rand_num = 10, 
                                 num_feat = 80, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_enet_80")

# quan 90 features
quan_batch_rand_enet_90 <- getRand(quan_cases_batch, 
                                 model = 'enet',
                                 rand_num = 10, 
                                 num_feat = 90, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_enet_90")

# quan 100 features
quan_batch_rand_enet_100 <- getRand(quan_cases_batch, 
                                  model = 'enet',
                                  rand_num = 10, 
                                  num_feat = 100, 
                                  exclude = T,
                                  union_features = union_features,
                                  data_name = "quan_batch_rand_enet_100")

# quan 500 features
quan_batch_rand_enet_500 <- getRand(quan_cases_batch, 
                                  model = 'enet',
                                  rand_num = 10, 
                                  num_feat = 500, 
                                  exclude = T,
                                  union_features = union_features,
                                  data_name = "quan_batch_rand_enet_500")

# quan 1000 features
quan_batch_rand_enet_1000 <- getRand(quan_cases_batch, 
                                   model = 'enet',
                                   rand_num = 10, 
                                   num_feat = 1000, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_batch_rand_enet_1000")

# quan 5000 features
quan_batch_rand_enet_5000 <- getRand(quan_cases_batch, 
                                   model = 'enet',
                                   rand_num = 10, 
                                   num_feat = 5000, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_batch_rand_enet_5000")

# quan 1000 features
quan_batch_rand_enet_10000 <- getRand(quan_cases_batch, 
                                    model = 'enet',
                                    rand_num = 10, 
                                    num_feat = 10000, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_batch_rand_enet_10000")



###########
# rbind tables and save RDA file
###########
quan_batch_rand_enet_table <- rbind(quan_batch_rand_enet_10, quan_batch_rand_enet_20, 
                                  quan_batch_rand_enet_30, quan_batch_rand_enet_40,
                                  quan_batch_rand_enet_50, quan_batch_rand_enet_60,
                                  quan_batch_rand_enet_70, quan_batch_rand_enet_80,
                                  quan_batch_rand_enet_90, quan_batch_rand_enet_100,
                                  quan_batch_rand_enet_500, quan_batch_rand_enet_1000,
                                  quan_batch_rand_enet_5000, quan_batch_rand_enet_10000)

# remove data 
rm(quan_batch_rand_enet_10, quan_batch_rand_enet_20, 
   quan_batch_rand_enet_30, quan_batch_rand_enet_40,
   quan_batch_rand_enet_50, quan_batch_rand_enet_60,
   quan_batch_rand_enet_70, quan_batch_rand_enet_80,
   quan_batch_rand_enet_90, quan_batch_rand_enet_100,
   quan_batch_rand_enet_500, quan_batch_rand_enet_1000,
   quan_batch_rand_enet_5000, quan_batch_rand_enet_10000)


#save table 
saveRDS(quan_batch_rand_enet_table, 
        file = paste0(rand_folder, '/quan_batch_rand_enet_table.rda'))



###########
# lasso
###########

# quan 10 features
quan_batch_rand_lasso_10 <- getRand(quan_cases_batch, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 10, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_lasso_10")

# quan 20 features
quan_batch_rand_lasso_20 <- getRand(quan_cases_batch, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 20, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_lasso_20")


# quan 30 features
quan_batch_rand_lasso_30 <- getRand(quan_cases_batch, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 30, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_lasso_30")

# quan 40 features
quan_batch_rand_lasso_40 <- getRand(quan_cases_batch, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 40, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_lasso_40")

# quan 50 features
quan_batch_rand_lasso_50 <- getRand(quan_cases_batch, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 50, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_lasso_50")

# quan 60 features
quan_batch_rand_lasso_60 <- getRand(quan_cases_batch, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 60, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_lasso_60")

# quan 70 features
quan_batch_rand_lasso_70 <- getRand(quan_cases_batch, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 70, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_lasso_70")

# quan 80 features
quan_batch_rand_lasso_80 <- getRand(quan_cases_batch, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 80, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_lasso_80")

# quan 90 features
quan_batch_rand_lasso_90 <- getRand(quan_cases_batch, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 90, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "quan_batch_rand_lasso_90")

# quan 100 features
quan_batch_rand_lasso_100 <- getRand(quan_cases_batch, 
                                  model = 'lasso',
                                  rand_num = 10, 
                                  num_feat = 100, 
                                  exclude = T,
                                  union_features = union_features,
                                  data_name = "quan_batch_rand_lasso_100")

# quan 500 features
quan_batch_rand_lasso_500 <- getRand(quan_cases_batch, 
                                  model = 'lasso',
                                  rand_num = 10, 
                                  num_feat = 500, 
                                  exclude = T,
                                  union_features = union_features,
                                  data_name = "quan_batch_rand_lasso_500")

# quan 1000 features
quan_batch_rand_lasso_1000 <- getRand(quan_cases_batch, 
                                   model = 'lasso',
                                   rand_num = 10, 
                                   num_feat = 1000, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_batch_rand_lasso_1000")

# quan 5000 features
quan_batch_rand_lasso_5000 <- getRand(quan_cases_batch, 
                                   model = 'lasso',
                                   rand_num = 10, 
                                   num_feat = 5000, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "quan_batch_rand_lasso_5000")

# quan 1000 features
quan_batch_rand_lasso_10000 <- getRand(quan_cases_batch, 
                                    model = 'lasso',
                                    rand_num = 10, 
                                    num_feat = 10000, 
                                    exclude = T,
                                    union_features = union_features,
                                    data_name = "quan_batch_rand_lasso_10000")



###########
# rbind tables and save RDA file
###########
quan_batch_rand_lasso_table <- rbind(quan_batch_rand_lasso_10, quan_batch_rand_lasso_20, 
                                  quan_batch_rand_lasso_30, quan_batch_rand_lasso_40,
                                  quan_batch_rand_lasso_50, quan_batch_rand_lasso_60,
                                  quan_batch_rand_lasso_70, quan_batch_rand_lasso_80,
                                  quan_batch_rand_lasso_90, quan_batch_rand_lasso_100,
                                  quan_batch_rand_lasso_500, quan_batch_rand_lasso_1000,
                                  quan_batch_rand_lasso_5000, quan_batch_rand_lasso_10000)

# remove data 
rm(quan_batch_rand_lasso_10, quan_batch_rand_lasso_20, 
   quan_batch_rand_lasso_30, quan_batch_rand_lasso_40,
   quan_batch_rand_lasso_50, quan_batch_rand_lasso_60,
   quan_batch_rand_lasso_70, quan_batch_rand_lasso_80,
   quan_batch_rand_lasso_90, quan_batch_rand_lasso_100,
   quan_batch_rand_lasso_500, quan_batch_rand_lasso_1000,
   quan_batch_rand_lasso_5000, quan_batch_rand_lasso_10000)


#save table 
saveRDS(quan_batch_rand_lasso_table, 
        file = paste0(rand_folder, '/quan_batch_rand_lasso_table.rda'))

