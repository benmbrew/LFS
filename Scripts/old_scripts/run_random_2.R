####### Script will run random models

##########
# initialize libraries
##########
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(doParallel)
library(nnet)
library(dplyr)

registerDoParallel(1)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')
results_folder <- paste0(project_folder, '/Results')
quan_folder <- paste0(results_folder, '/quan')
rand_folder <- paste0(results_folder, '/rand')
scripts_folder <- paste0(project_folder, '/Scripts')


##########
# load data batch data
##########
# read cases
quan_cases <- readRDS(paste0(model_data, '/quan_cases.rda'))

# read batch corrected data for gender
quan_cases_gen <- readRDS(paste0(model_data, '/quan_cases_gen.rda'))

# read batch corrected data for sentrix id and SAM
quan_cases_sen <- readRDS(paste0(model_data, '/quan_cases_sen.rda'))

quan_cases_sam <- readRDS(paste0(model_data, '/quan_cases_sam.rda'))

# read batch corrected data for sentrix id and SAM and gender!
quan_cases_sen_gen <- readRDS(paste0(model_data, '/quan_cases_sen_gen.rda'))

quan_cases_sam_gen <- readRDS(paste0(model_data, '/quan_cases_sam_gen.rda'))

# load features
load(paste0(model_data, '/bh_feat.RData'))

source(paste0(scripts_folder, '/predict_age/functions.R'))

##########
# function to run models - subset, get residuals, get categorical, predict with regression and fac. 
##########
runModels <- function(data,
                      model,
                      random = F,
                      bump_hunter = F,
                      bump_hunter_data,
                      num_feat = NULL,
                      seed_num) 
{
  
  # get differenct variations of data
  if (random) {
    
    data <- subsetDat(data, 
                      random = T, 
                      num_feat,
                      seed_num)
  } else {
    
    data <- subsetDat(data, 
                      random = F, 
                      num_feat = NULL)
  }
  
  #########
  # get regression data and run regression
  #########
  # get data
  if (bump_hunter) {
    
    data <- bhSubset(data, bh_data = bump_hunter_data)
    
  }
  
  # data_resid <- getResidual(data)
  
  if (model == 'rf') {
    
    data_result <- rfPredictReg(data, cutoff = .7, iterations = 5)
    # data_resid_result <- rfPredictReg(data_resid, cutoff = .7, iterations = 5)
  }
  
  if (model == 'enet') {
    
    data_result <- enetPredReg(data, N_CV_REPEATS = 2, nfolds = 5,cutoff = .7, iterations = 5)
    # data_resid_result <- enetPredReg(data_resid, N_CV_REPEATS = 2, nfolds = 5,cutoff = .7, iterations = 5)
    
  }
  
  if (model == 'lasso') {
    
    data_result <- regPred(data, alpha = 1, nfolds = 5, cutoff = .7, iterations = 5)
    # data_resid_result <- regPred(data_resid, alpha = 1, nfolds = 5, cutoff = .7, iterations = 5)
    
  }
  
  # #########
  # # get classification data and run classifier
  # #########
  # data_fac <- makeFac(data, threshold = 48)
  # data_resid_fac <- makeFac(data_resid, threshold = 48)
  # 
  # data_fac_result  <- rfPredictFac(data_fac, cutoff = .7, iterations = 5)
  # data_resid_fac_result <- rfPredictFac(data_resid_fac, cutoff = .7, iterations = 5)
  
  return (data_result) 
  # data_fac_result, 
  # data_resid_fac_result))
  
  
}




###########
# function for random
###########
# data <- quan_cases
# rand_num <- 2
# num_feat <- 10
# data_name <- 'raw_rand_10'
# exlcude <- T
# union_features <-  union_features
getRand <- function(data, model, 
                    rand_num, 
                    num_feat, 
                    exclude, 
                    union_features, 
                    data_name) 
{
  
  if (exclude) {
    
    features <- colnames(data[, 10:ncol(data)])
    keep_these <- features[!features %in% union_features]
    data <- data[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
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
#10

###########
# rf
###########

# ungen 10 features
ungen_rand_rf_10 <- getRand(quan_cases, 
                             model = 'rf',
                             rand_num = 10, 
                             num_feat = 10, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "ungen_rand_rf_10")

# gen 10 features
gen_rand_rf_10 <- getRand(quan_cases_gen, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 10, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "gen_rand_rf_10")

# unsam 10 features
sam_rand_rf_10 <- getRand(quan_cases_sam, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 10, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sam_rand_rf_10")

# unsen 10 features
sen_rand_rf_10 <- getRand(quan_cases_sen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 10, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sen_rand_rf_10")

# unsam_gen 10 features
sam_gen_rand_rf_10 <- getRand(quan_cases_sam_gen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 10, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sam_gen_rand_rf_10")

# unsen_gen 10 features
sen_gen_rand_rf_10 <- getRand(quan_cases_sen_gen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 10, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sen_gen_rand_rf_10")




###########
# enet
###########

# ungen 10 features
ungen_rand_enet_10 <- getRand(quan_cases, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 10, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "ungen_rand_enet_10")

# gen 10 features
gen_rand_enet_10 <- getRand(quan_cases_gen, 
                          model = 'enet',
                          rand_num = 10, 
                          num_feat = 10, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "gen_rand_enet_10")

# unsam 10 features
sam_rand_enet_10 <- getRand(quan_cases_sam, 
                          model = 'enet',
                          rand_num = 10, 
                          num_feat = 10, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sam_rand_enet_10")

# unsen 10 features
sen_rand_enet_10 <- getRand(quan_cases_sen, 
                          model = 'enet',
                          rand_num = 10, 
                          num_feat = 10, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sen_rand_enet_10")

# unsam_gen 10 features
sam_gen_rand_enet_10 <- getRand(quan_cases_sam_gen, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 10, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sam_gen_rand_enet_10")

# unsen_gen 10 features
sen_gen_rand_enet_10 <- getRand(quan_cases_sen_gen, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 10, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sen_gen_rand_enet_10")


###########
# lasso
###########

# ungen 10 features
ungen_rand_lasso_10 <- getRand(quan_cases, 
                            model = 'lasso',
                            rand_num = 10, 
                            num_feat = 10, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "ungen_rand_lasso_10")

# gen 10 features
gen_rand_lasso_10 <- getRand(quan_cases_gen, 
                          model = 'lasso',
                          rand_num = 10, 
                          num_feat = 10, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "gen_rand_lasso_10")

# unsam 10 features
sam_rand_lasso_10 <- getRand(quan_cases_sam, 
                          model = 'lasso',
                          rand_num = 10, 
                          num_feat = 10, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sam_rand_lasso_10")

# unsen 10 features
sen_rand_lasso_10 <- getRand(quan_cases_sen, 
                          model = 'lasso',
                          rand_num = 10, 
                          num_feat = 10, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sen_rand_lasso_10")

# unsam_gen 10 features
sam_gen_rand_lasso_10 <- getRand(quan_cases_sam_gen, 
                              model = 'lasso',
                              rand_num = 10, 
                              num_feat = 10, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sam_gen_rand_lasso_10")

# unsen_gen 10 features
sen_gen_rand_lasso_10 <- getRand(quan_cases_sen_gen, 
                              model = 'lasso',
                              rand_num = 10, 
                              num_feat = 10, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sen_gen_rand_lasso_10")





rand_10 <- rbind(ungen_rand_rf_10,
                 gen_rand_rf_10,
                 sen_rand_rf_10,
                 sam_rand_rf_10,
                 sam_gen_rand_rf_10,
                 sen_gen_rand_rf_10,
                 ungen_rand_enet_10,
                 gen_rand_enet_10,
                 sen_rand_enet_10,
                 sam_rand_enet_10,
                 sam_gen_rand_enet_10,
                 sen_gen_rand_enet_10,
                 ungen_rand_lasso_10,
                 gen_rand_lasso_10,
                 sen_rand_lasso_10,
                 sam_rand_lasso_10,
                 sam_gen_rand_lasso_10,
                 sen_gen_rand_lasso_10)

rm(ungen_rand_rf_10,
   gen_rand_rf_10,
   sen_rand_rf_10,
   sam_rand_rf_10,
   sam_gen_rand_rf_10,
   sen_gen_rand_rf_10,
   ungen_rand_enet_10,
   gen_rand_enet_10,
   sen_rand_enet_10,
   sam_rand_enet_10,
   sam_gen_rand_enet_10,
   sen_gen_rand_enet_10,
   ungen_rand_lasso_10,
   gen_rand_lasso_10,
   sen_rand_lasso_10,
   sam_rand_lasso_10,
   sam_gen_rand_lasso_10,
   sen_gen_rand_lasso_10)

#############################################################################################################
#20

###########
# rf
###########

# ungen 20 features
ungen_rand_rf_20 <- getRand(quan_cases, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 20, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "ungen_rand_rf_20")

# gen 20 features
gen_rand_rf_20 <- getRand(quan_cases_gen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 20, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "gen_rand_rf_20")

# unsam 20 features
sam_rand_rf_20 <- getRand(quan_cases_sam, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 20, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sam_rand_rf_20")

# unsen 20 features
sen_rand_rf_20 <- getRand(quan_cases_sen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 20, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sen_rand_rf_20")

# unsam_gen 20 features
sam_gen_rand_rf_20 <- getRand(quan_cases_sam_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 20, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sam_gen_rand_rf_20")

# unsen_gen 20 features
sen_gen_rand_rf_20 <- getRand(quan_cases_sen_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 20, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sen_gen_rand_rf_20")




###########
# enet
###########

# ungen 20 features
ungen_rand_enet_20 <- getRand(quan_cases, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 20, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "ungen_rand_enet_20")

# gen 20 features
gen_rand_enet_20 <- getRand(quan_cases_gen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 20, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "gen_rand_enet_20")

# unsam 20 features
sam_rand_enet_20 <- getRand(quan_cases_sam, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 20, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sam_rand_enet_20")

# unsen 20 features
sen_rand_enet_20 <- getRand(quan_cases_sen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 20, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sen_rand_enet_20")

# unsam_gen 20 features
sam_gen_rand_enet_20 <- getRand(quan_cases_sam_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 20, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sam_gen_rand_enet_20")

# unsen_gen 20 features
sen_gen_rand_enet_20 <- getRand(quan_cases_sen_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 20, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sen_gen_rand_enet_20")


###########
# lasso
###########

# ungen 20 features
ungen_rand_lasso_20 <- getRand(quan_cases, 
                               model = 'lasso',
                               rand_num = 10, 
                               num_feat = 20, 
                               exclude = T,
                               union_features = union_features,
                               data_name = "ungen_rand_lasso_20")

# gen 20 features
gen_rand_lasso_20 <- getRand(quan_cases_gen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 20, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "gen_rand_lasso_20")

# unsam 20 features
sam_rand_lasso_20 <- getRand(quan_cases_sam, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 20, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sam_rand_lasso_20")

# unsen 20 features
sen_rand_lasso_20 <- getRand(quan_cases_sen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 20, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sen_rand_lasso_20")

# unsam_gen 20 features
sam_gen_rand_lasso_20 <- getRand(quan_cases_sam_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 20, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sam_gen_rand_lasso_20")

# unsen_gen 20 features
sen_gen_rand_lasso_20 <- getRand(quan_cases_sen_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 20, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sen_gen_rand_lasso_20")





rand_20 <- rbind(ungen_rand_rf_20,
                 gen_rand_rf_20,
                 sen_rand_rf_20,
                 sam_rand_rf_20,
                 sam_gen_rand_rf_20,
                 sen_gen_rand_rf_20,
                 ungen_rand_enet_20,
                 gen_rand_enet_20,
                 sen_rand_enet_20,
                 sam_rand_enet_20,
                 sam_gen_rand_enet_20,
                 sen_gen_rand_enet_20,
                 ungen_rand_lasso_20,
                 gen_rand_lasso_20,
                 sen_rand_lasso_20,
                 sam_rand_lasso_20,
                 sam_gen_rand_lasso_20,
                 sen_gen_rand_lasso_20)

rm(ungen_rand_rf_20,
   gen_rand_rf_20,
   sen_rand_rf_20,
   sam_rand_rf_20,
   sam_gen_rand_rf_20,
   sen_gen_rand_rf_20,
   ungen_rand_enet_20,
   gen_rand_enet_20,
   sen_rand_enet_20,
   sam_rand_enet_20,
   sam_gen_rand_enet_20,
   sen_gen_rand_enet_20,
   ungen_rand_lasso_20,
   gen_rand_lasso_20,
   sen_rand_lasso_20,
   sam_rand_lasso_20,
   sam_gen_rand_lasso_20,
   sen_gen_rand_lasso_20)

#############################################################################################################
#30

###########
# rf
###########

# ungen 30 features
ungen_rand_rf_30 <- getRand(quan_cases, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 30, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "ungen_rand_rf_30")

# gen 30 features
gen_rand_rf_30 <- getRand(quan_cases_gen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 30, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "gen_rand_rf_30")

# unsam 30 features
sam_rand_rf_30 <- getRand(quan_cases_sam, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 30, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sam_rand_rf_30")

# unsen 30 features
sen_rand_rf_30 <- getRand(quan_cases_sen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 30, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sen_rand_rf_30")

# unsam_gen 30 features
sam_gen_rand_rf_30 <- getRand(quan_cases_sam_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 30, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sam_gen_rand_rf_30")

# unsen_gen 30 features
sen_gen_rand_rf_30 <- getRand(quan_cases_sen_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 30, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sen_gen_rand_rf_30")




###########
# enet
###########

# ungen 30 features
ungen_rand_enet_30 <- getRand(quan_cases, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 30, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "ungen_rand_enet_30")

# gen 30 features
gen_rand_enet_30 <- getRand(quan_cases_gen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 30, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "gen_rand_enet_30")

# unsam 30 features
sam_rand_enet_30 <- getRand(quan_cases_sam, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 30, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sam_rand_enet_30")

# unsen 30 features
sen_rand_enet_30 <- getRand(quan_cases_sen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 30, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sen_rand_enet_30")

# unsam_gen 30 features
sam_gen_rand_enet_30 <- getRand(quan_cases_sam_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 30, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sam_gen_rand_enet_30")

# unsen_gen 30 features
sen_gen_rand_enet_30 <- getRand(quan_cases_sen_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 30, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sen_gen_rand_enet_30")


###########
# lasso
###########

# ungen 30 features
ungen_rand_lasso_30 <- getRand(quan_cases, 
                               model = 'lasso',
                               rand_num = 10, 
                               num_feat = 30, 
                               exclude = T,
                               union_features = union_features,
                               data_name = "ungen_rand_lasso_30")

# gen 30 features
gen_rand_lasso_30 <- getRand(quan_cases_gen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 30, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "gen_rand_lasso_30")

# unsam 30 features
sam_rand_lasso_30 <- getRand(quan_cases_sam, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 30, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sam_rand_lasso_30")

# unsen 30 features
sen_rand_lasso_30 <- getRand(quan_cases_sen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 30, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sen_rand_lasso_30")

# unsam_gen 30 features
sam_gen_rand_lasso_30 <- getRand(quan_cases_sam_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 30, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sam_gen_rand_lasso_30")

# unsen_gen 30 features
sen_gen_rand_lasso_30 <- getRand(quan_cases_sen_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 30, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sen_gen_rand_lasso_30")





rand_30 <- rbind(ungen_rand_rf_30,
                 gen_rand_rf_30,
                 sen_rand_rf_30,
                 sam_rand_rf_30,
                 sam_gen_rand_rf_30,
                 sen_gen_rand_rf_30,
                 ungen_rand_enet_30,
                 gen_rand_enet_30,
                 sen_rand_enet_30,
                 sam_rand_enet_30,
                 sam_gen_rand_enet_30,
                 sen_gen_rand_enet_30,
                 ungen_rand_lasso_30,
                 gen_rand_lasso_30,
                 sen_rand_lasso_30,
                 sam_rand_lasso_30,
                 sam_gen_rand_lasso_30,
                 sen_gen_rand_lasso_30)

rm(ungen_rand_rf_30,
   gen_rand_rf_30,
   sen_rand_rf_30,
   sam_rand_rf_30,
   sam_gen_rand_rf_30,
   sen_gen_rand_rf_30,
   ungen_rand_enet_30,
   gen_rand_enet_30,
   sen_rand_enet_30,
   sam_rand_enet_30,
   sam_gen_rand_enet_30,
   sen_gen_rand_enet_30,
   ungen_rand_lasso_30,
   gen_rand_lasso_30,
   sen_rand_lasso_30,
   sam_rand_lasso_30,
   sam_gen_rand_lasso_30,
   sen_gen_rand_lasso_30)


#############################################################################################################
#50

###########
# rf
###########

# ungen 50 features
ungen_rand_rf_50 <- getRand(quan_cases, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 50, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "ungen_rand_rf_50")

# gen 50 features
gen_rand_rf_50 <- getRand(quan_cases_gen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 50, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "gen_rand_rf_50")

# unsam 50 features
sam_rand_rf_50 <- getRand(quan_cases_sam, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 50, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sam_rand_rf_50")

# unsen 50 features
sen_rand_rf_50 <- getRand(quan_cases_sen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 50, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sen_rand_rf_50")

# unsam_gen 50 features
sam_gen_rand_rf_50 <- getRand(quan_cases_sam_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 50, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sam_gen_rand_rf_50")

# unsen_gen 50 features
sen_gen_rand_rf_50 <- getRand(quan_cases_sen_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 50, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sen_gen_rand_rf_50")




###########
# enet
###########

# ungen 50 features
ungen_rand_enet_50 <- getRand(quan_cases, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 50, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "ungen_rand_enet_50")

# gen 50 features
gen_rand_enet_50 <- getRand(quan_cases_gen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 50, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "gen_rand_enet_50")

# unsam 50 features
sam_rand_enet_50 <- getRand(quan_cases_sam, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 50, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sam_rand_enet_50")

# unsen 50 features
sen_rand_enet_50 <- getRand(quan_cases_sen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 50, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sen_rand_enet_50")

# unsam_gen 50 features
sam_gen_rand_enet_50 <- getRand(quan_cases_sam_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 50, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sam_gen_rand_enet_50")

# unsen_gen 50 features
sen_gen_rand_enet_50 <- getRand(quan_cases_sen_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 50, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sen_gen_rand_enet_50")


###########
# lasso
###########

# ungen 50 features
ungen_rand_lasso_50 <- getRand(quan_cases, 
                               model = 'lasso',
                               rand_num = 10, 
                               num_feat = 50, 
                               exclude = T,
                               union_features = union_features,
                               data_name = "ungen_rand_lasso_50")

# gen 50 features
gen_rand_lasso_50 <- getRand(quan_cases_gen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 50, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "gen_rand_lasso_50")

# unsam 50 features
sam_rand_lasso_50 <- getRand(quan_cases_sam, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 50, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sam_rand_lasso_50")

# unsen 50 features
sen_rand_lasso_50 <- getRand(quan_cases_sen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 50, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sen_rand_lasso_50")

# unsam_gen 50 features
sam_gen_rand_lasso_50 <- getRand(quan_cases_sam_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 50, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sam_gen_rand_lasso_50")

# unsen_gen 50 features
sen_gen_rand_lasso_50 <- getRand(quan_cases_sen_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 50, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sen_gen_rand_lasso_50")





rand_50 <- rbind(ungen_rand_rf_50,
                 gen_rand_rf_50,
                 sen_rand_rf_50,
                 sam_rand_rf_50,
                 sam_gen_rand_rf_50,
                 sen_gen_rand_rf_50,
                 ungen_rand_enet_50,
                 gen_rand_enet_50,
                 sen_rand_enet_50,
                 sam_rand_enet_50,
                 sam_gen_rand_enet_50,
                 sen_gen_rand_enet_50,
                 ungen_rand_lasso_50,
                 gen_rand_lasso_50,
                 sen_rand_lasso_50,
                 sam_rand_lasso_50,
                 sam_gen_rand_lasso_50,
                 sen_gen_rand_lasso_50)

rm(ungen_rand_rf_50,
   gen_rand_rf_50,
   sen_rand_rf_50,
   sam_rand_rf_50,
   sam_gen_rand_rf_50,
   sen_gen_rand_rf_50,
   ungen_rand_enet_50,
   gen_rand_enet_50,
   sen_rand_enet_50,
   sam_rand_enet_50,
   sam_gen_rand_enet_50,
   sen_gen_rand_enet_50,
   ungen_rand_lasso_50,
   gen_rand_lasso_50,
   sen_rand_lasso_50,
   sam_rand_lasso_50,
   sam_gen_rand_lasso_50,
   sen_gen_rand_lasso_50)


#############################################################################################################
#100

###########
# rf
###########

# ungen 100 features
ungen_rand_rf_100 <- getRand(quan_cases, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 100, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "ungen_rand_rf_100")

# gen 100 features
gen_rand_rf_100 <- getRand(quan_cases_gen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 100, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "gen_rand_rf_100")

# unsam 100 features
sam_rand_rf_100 <- getRand(quan_cases_sam, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 100, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sam_rand_rf_100")

# unsen 100 features
sen_rand_rf_100 <- getRand(quan_cases_sen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 100, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sen_rand_rf_100")

# unsam_gen 100 features
sam_gen_rand_rf_100 <- getRand(quan_cases_sam_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 100, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sam_gen_rand_rf_100")

# unsen_gen 100 features
sen_gen_rand_rf_100 <- getRand(quan_cases_sen_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 100, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sen_gen_rand_rf_100")




###########
# enet
###########

# ungen 100 features
ungen_rand_enet_100 <- getRand(quan_cases, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 100, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "ungen_rand_enet_100")

# gen 100 features
gen_rand_enet_100 <- getRand(quan_cases_gen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 100, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "gen_rand_enet_100")

# unsam 100 features
sam_rand_enet_100 <- getRand(quan_cases_sam, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 100, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sam_rand_enet_100")

# unsen 100 features
sen_rand_enet_100 <- getRand(quan_cases_sen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 100, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sen_rand_enet_100")

# unsam_gen 100 features
sam_gen_rand_enet_100 <- getRand(quan_cases_sam_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 100, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sam_gen_rand_enet_100")

# unsen_gen 100 features
sen_gen_rand_enet_100 <- getRand(quan_cases_sen_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 100, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sen_gen_rand_enet_100")


###########
# lasso
###########

# ungen 100 features
ungen_rand_lasso_100 <- getRand(quan_cases, 
                               model = 'lasso',
                               rand_num = 10, 
                               num_feat = 100, 
                               exclude = T,
                               union_features = union_features,
                               data_name = "ungen_rand_lasso_100")

# gen 100 features
gen_rand_lasso_100 <- getRand(quan_cases_gen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 100, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "gen_rand_lasso_100")

# unsam 100 features
sam_rand_lasso_100 <- getRand(quan_cases_sam, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 100, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sam_rand_lasso_100")

# unsen 100 features
sen_rand_lasso_100 <- getRand(quan_cases_sen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 100, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sen_rand_lasso_100")

# unsam_gen 100 features
sam_gen_rand_lasso_100 <- getRand(quan_cases_sam_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 100, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sam_gen_rand_lasso_100")

# unsen_gen 100 features
sen_gen_rand_lasso_100 <- getRand(quan_cases_sen_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 100, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sen_gen_rand_lasso_100")





rand_100 <- rbind(ungen_rand_rf_100,
                 gen_rand_rf_100,
                 sen_rand_rf_100,
                 sam_rand_rf_100,
                 sam_gen_rand_rf_100,
                 sen_gen_rand_rf_100,
                 ungen_rand_enet_100,
                 gen_rand_enet_100,
                 sen_rand_enet_100,
                 sam_rand_enet_100,
                 sam_gen_rand_enet_100,
                 sen_gen_rand_enet_100,
                 ungen_rand_lasso_100,
                 gen_rand_lasso_100,
                 sen_rand_lasso_100,
                 sam_rand_lasso_100,
                 sam_gen_rand_lasso_100,
                 sen_gen_rand_lasso_100)

rm(ungen_rand_rf_100,
   gen_rand_rf_100,
   sen_rand_rf_100,
   sam_rand_rf_100,
   sam_gen_rand_rf_100,
   sen_gen_rand_rf_100,
   ungen_rand_enet_100,
   gen_rand_enet_100,
   sen_rand_enet_100,
   sam_rand_enet_100,
   sam_gen_rand_enet_100,
   sen_gen_rand_enet_100,
   ungen_rand_lasso_100,
   gen_rand_lasso_100,
   sen_rand_lasso_100,
   sam_rand_lasso_100,
   sam_gen_rand_lasso_100,
   sen_gen_rand_lasso_100)


#############################################################################################################
#300

###########
# rf
###########

# ungen 300 features
ungen_rand_rf_300 <- getRand(quan_cases, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 300, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "ungen_rand_rf_300")

# gen 300 features
gen_rand_rf_300 <- getRand(quan_cases_gen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 300, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "gen_rand_rf_300")

# unsam 300 features
sam_rand_rf_300 <- getRand(quan_cases_sam, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 300, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sam_rand_rf_300")

# unsen 300 features
sen_rand_rf_300 <- getRand(quan_cases_sen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 300, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sen_rand_rf_300")

# unsam_gen 300 features
sam_gen_rand_rf_300 <- getRand(quan_cases_sam_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 300, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sam_gen_rand_rf_300")

# unsen_gen 300 features
sen_gen_rand_rf_300 <- getRand(quan_cases_sen_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 300, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sen_gen_rand_rf_300")




###########
# enet
###########

# ungen 300 features
ungen_rand_enet_300 <- getRand(quan_cases, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 300, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "ungen_rand_enet_300")

# gen 300 features
gen_rand_enet_300 <- getRand(quan_cases_gen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 300, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "gen_rand_enet_300")

# unsam 300 features
sam_rand_enet_300 <- getRand(quan_cases_sam, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 300, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sam_rand_enet_300")

# unsen 300 features
sen_rand_enet_300 <- getRand(quan_cases_sen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 300, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sen_rand_enet_300")

# unsam_gen 300 features
sam_gen_rand_enet_300 <- getRand(quan_cases_sam_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 300, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sam_gen_rand_enet_300")

# unsen_gen 300 features
sen_gen_rand_enet_300 <- getRand(quan_cases_sen_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 300, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sen_gen_rand_enet_300")


###########
# lasso
###########

# ungen 300 features
ungen_rand_lasso_300 <- getRand(quan_cases, 
                               model = 'lasso',
                               rand_num = 10, 
                               num_feat = 300, 
                               exclude = T,
                               union_features = union_features,
                               data_name = "ungen_rand_lasso_300")

# gen 300 features
gen_rand_lasso_300 <- getRand(quan_cases_gen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 300, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "gen_rand_lasso_300")

# unsam 300 features
sam_rand_lasso_300 <- getRand(quan_cases_sam, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 300, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sam_rand_lasso_300")

# unsen 300 features
sen_rand_lasso_300 <- getRand(quan_cases_sen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 300, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sen_rand_lasso_300")

# unsam_gen 300 features
sam_gen_rand_lasso_300 <- getRand(quan_cases_sam_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 300, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sam_gen_rand_lasso_300")

# unsen_gen 300 features
sen_gen_rand_lasso_300 <- getRand(quan_cases_sen_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 300, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sen_gen_rand_lasso_300")





rand_300 <- rbind(ungen_rand_rf_300,
                 gen_rand_rf_300,
                 sen_rand_rf_300,
                 sam_rand_rf_300,
                 sam_gen_rand_rf_300,
                 sen_gen_rand_rf_300,
                 ungen_rand_enet_300,
                 gen_rand_enet_300,
                 sen_rand_enet_300,
                 sam_rand_enet_300,
                 sam_gen_rand_enet_300,
                 sen_gen_rand_enet_300,
                 ungen_rand_lasso_300,
                 gen_rand_lasso_300,
                 sen_rand_lasso_300,
                 sam_rand_lasso_300,
                 sam_gen_rand_lasso_300,
                 sen_gen_rand_lasso_300)

rm(ungen_rand_rf_300,
   gen_rand_rf_300,
   sen_rand_rf_300,
   sam_rand_rf_300,
   sam_gen_rand_rf_300,
   sen_gen_rand_rf_300,
   ungen_rand_enet_300,
   gen_rand_enet_300,
   sen_rand_enet_300,
   sam_rand_enet_300,
   sam_gen_rand_enet_300,
   sen_gen_rand_enet_300,
   ungen_rand_lasso_300,
   gen_rand_lasso_300,
   sen_rand_lasso_300,
   sam_rand_lasso_300,
   sam_gen_rand_lasso_300,
   sen_gen_rand_lasso_300)

#############################################################################################################
#500

###########
# rf
###########

# ungen 500 features
ungen_rand_rf_500 <- getRand(quan_cases, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 500, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "ungen_rand_rf_500")

# gen 500 features
gen_rand_rf_500 <- getRand(quan_cases_gen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 500, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "gen_rand_rf_500")

# unsam 500 features
sam_rand_rf_500 <- getRand(quan_cases_sam, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 500, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sam_rand_rf_500")

# unsen 500 features
sen_rand_rf_500 <- getRand(quan_cases_sen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 500, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sen_rand_rf_500")

# unsam_gen 500 features
sam_gen_rand_rf_500 <- getRand(quan_cases_sam_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 500, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sam_gen_rand_rf_500")

# unsen_gen 500 features
sen_gen_rand_rf_500 <- getRand(quan_cases_sen_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 500, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sen_gen_rand_rf_500")




###########
# enet
###########

# ungen 500 features
ungen_rand_enet_500 <- getRand(quan_cases, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 500, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "ungen_rand_enet_500")

# gen 500 features
gen_rand_enet_500 <- getRand(quan_cases_gen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 500, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "gen_rand_enet_500")

# unsam 500 features
sam_rand_enet_500 <- getRand(quan_cases_sam, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 500, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sam_rand_enet_500")

# unsen 500 features
sen_rand_enet_500 <- getRand(quan_cases_sen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 500, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sen_rand_enet_500")

# unsam_gen 500 features
sam_gen_rand_enet_500 <- getRand(quan_cases_sam_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 500, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sam_gen_rand_enet_500")

# unsen_gen 500 features
sen_gen_rand_enet_500 <- getRand(quan_cases_sen_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 500, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sen_gen_rand_enet_500")


###########
# lasso
###########

# ungen 500 features
ungen_rand_lasso_500 <- getRand(quan_cases, 
                               model = 'lasso',
                               rand_num = 10, 
                               num_feat = 500, 
                               exclude = T,
                               union_features = union_features,
                               data_name = "ungen_rand_lasso_500")

# gen 500 features
gen_rand_lasso_500 <- getRand(quan_cases_gen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 500, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "gen_rand_lasso_500")

# unsam 500 features
sam_rand_lasso_500 <- getRand(quan_cases_sam, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 500, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sam_rand_lasso_500")

# unsen 500 features
sen_rand_lasso_500 <- getRand(quan_cases_sen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 500, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sen_rand_lasso_500")

# unsam_gen 500 features
sam_gen_rand_lasso_500 <- getRand(quan_cases_sam_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 500, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sam_gen_rand_lasso_500")

# unsen_gen 500 features
sen_gen_rand_lasso_500 <- getRand(quan_cases_sen_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 500, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sen_gen_rand_lasso_500")





rand_500 <- rbind(ungen_rand_rf_500,
                 gen_rand_rf_500,
                 sen_rand_rf_500,
                 sam_rand_rf_500,
                 sam_gen_rand_rf_500,
                 sen_gen_rand_rf_500,
                 ungen_rand_enet_500,
                 gen_rand_enet_500,
                 sen_rand_enet_500,
                 sam_rand_enet_500,
                 sam_gen_rand_enet_500,
                 sen_gen_rand_enet_500,
                 ungen_rand_lasso_500,
                 gen_rand_lasso_500,
                 sen_rand_lasso_500,
                 sam_rand_lasso_500,
                 sam_gen_rand_lasso_500,
                 sen_gen_rand_lasso_500)

rm(ungen_rand_rf_500,
   gen_rand_rf_500,
   sen_rand_rf_500,
   sam_rand_rf_500,
   sam_gen_rand_rf_500,
   sen_gen_rand_rf_500,
   ungen_rand_enet_500,
   gen_rand_enet_500,
   sen_rand_enet_500,
   sam_rand_enet_500,
   sam_gen_rand_enet_500,
   sen_gen_rand_enet_500,
   ungen_rand_lasso_500,
   gen_rand_lasso_500,
   sen_rand_lasso_500,
   sam_rand_lasso_500,
   sam_gen_rand_lasso_500,
   sen_gen_rand_lasso_500)


#############################################################################################################
#1000

###########
# rf
###########

# ungen 1000 features
ungen_rand_rf_1000 <- getRand(quan_cases, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 1000, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "ungen_rand_rf_1000")

# gen 1000 features
gen_rand_rf_1000 <- getRand(quan_cases_gen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 1000, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "gen_rand_rf_1000")

# unsam 1000 features
sam_rand_rf_1000 <- getRand(quan_cases_sam, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 1000, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sam_rand_rf_1000")

# unsen 1000 features
sen_rand_rf_1000 <- getRand(quan_cases_sen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 1000, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sen_rand_rf_1000")

# unsam_gen 1000 features
sam_gen_rand_rf_1000 <- getRand(quan_cases_sam_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 1000, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sam_gen_rand_rf_1000")

# unsen_gen 1000 features
sen_gen_rand_rf_1000 <- getRand(quan_cases_sen_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 1000, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sen_gen_rand_rf_1000")




###########
# enet
###########

# ungen 1000 features
ungen_rand_enet_1000 <- getRand(quan_cases, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 1000, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "ungen_rand_enet_1000")

# gen 1000 features
gen_rand_enet_1000 <- getRand(quan_cases_gen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 1000, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "gen_rand_enet_1000")

# unsam 1000 features
sam_rand_enet_1000 <- getRand(quan_cases_sam, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 1000, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sam_rand_enet_1000")

# unsen 1000 features
sen_rand_enet_1000 <- getRand(quan_cases_sen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 1000, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sen_rand_enet_1000")

# unsam_gen 1000 features
sam_gen_rand_enet_1000 <- getRand(quan_cases_sam_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 1000, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sam_gen_rand_enet_1000")

# unsen_gen 1000 features
sen_gen_rand_enet_1000 <- getRand(quan_cases_sen_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 1000, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sen_gen_rand_enet_1000")


###########
# lasso
###########

# ungen 1000 features
ungen_rand_lasso_1000 <- getRand(quan_cases, 
                               model = 'lasso',
                               rand_num = 10, 
                               num_feat = 1000, 
                               exclude = T,
                               union_features = union_features,
                               data_name = "ungen_rand_lasso_1000")

# gen 1000 features
gen_rand_lasso_1000 <- getRand(quan_cases_gen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 1000, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "gen_rand_lasso_1000")

# unsam 1000 features
sam_rand_lasso_1000 <- getRand(quan_cases_sam, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 1000, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sam_rand_lasso_1000")

# unsen 1000 features
sen_rand_lasso_1000 <- getRand(quan_cases_sen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 1000, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sen_rand_lasso_1000")

# unsam_gen 1000 features
sam_gen_rand_lasso_1000 <- getRand(quan_cases_sam_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 1000, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sam_gen_rand_lasso_1000")

# unsen_gen 1000 features
sen_gen_rand_lasso_1000 <- getRand(quan_cases_sen_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 1000, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sen_gen_rand_lasso_1000")





rand_1000 <- rbind(ungen_rand_rf_1000,
                 gen_rand_rf_1000,
                 sen_rand_rf_1000,
                 sam_rand_rf_1000,
                 sam_gen_rand_rf_1000,
                 sen_gen_rand_rf_1000,
                 ungen_rand_enet_1000,
                 gen_rand_enet_1000,
                 sen_rand_enet_1000,
                 sam_rand_enet_1000,
                 sam_gen_rand_enet_1000,
                 sen_gen_rand_enet_1000,
                 ungen_rand_lasso_1000,
                 gen_rand_lasso_1000,
                 sen_rand_lasso_1000,
                 sam_rand_lasso_1000,
                 sam_gen_rand_lasso_1000,
                 sen_gen_rand_lasso_1000)

rm(ungen_rand_rf_1000,
   gen_rand_rf_1000,
   sen_rand_rf_1000,
   sam_rand_rf_1000,
   sam_gen_rand_rf_1000,
   sen_gen_rand_rf_1000,
   ungen_rand_enet_1000,
   gen_rand_enet_1000,
   sen_rand_enet_1000,
   sam_rand_enet_1000,
   sam_gen_rand_enet_1000,
   sen_gen_rand_enet_1000,
   ungen_rand_lasso_1000,
   gen_rand_lasso_1000,
   sen_rand_lasso_1000,
   sam_rand_lasso_1000,
   sam_gen_rand_lasso_1000,
   sen_gen_rand_lasso_1000)



#############################################################################################################
#5000

###########
# rf
###########

# ungen 5000 features
ungen_rand_rf_5000 <- getRand(quan_cases, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 5000, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "ungen_rand_rf_5000")

# gen 5000 features
gen_rand_rf_5000 <- getRand(quan_cases_gen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 5000, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "gen_rand_rf_5000")

# unsam 5000 features
sam_rand_rf_5000 <- getRand(quan_cases_sam, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 5000, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sam_rand_rf_5000")

# unsen 5000 features
sen_rand_rf_5000 <- getRand(quan_cases_sen, 
                          model = 'rf',
                          rand_num = 10, 
                          num_feat = 5000, 
                          exclude = T,
                          union_features = union_features,
                          data_name = "sen_rand_rf_5000")

# unsam_gen 5000 features
sam_gen_rand_rf_5000 <- getRand(quan_cases_sam_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 5000, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sam_gen_rand_rf_5000")

# unsen_gen 5000 features
sen_gen_rand_rf_5000 <- getRand(quan_cases_sen_gen, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 5000, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sen_gen_rand_rf_5000")




###########
# enet
###########

# ungen 5000 features
ungen_rand_enet_5000 <- getRand(quan_cases, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 5000, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "ungen_rand_enet_5000")

# gen 5000 features
gen_rand_enet_5000 <- getRand(quan_cases_gen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 5000, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "gen_rand_enet_5000")

# unsam 5000 features
sam_rand_enet_5000 <- getRand(quan_cases_sam, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 5000, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sam_rand_enet_5000")

# unsen 5000 features
sen_rand_enet_5000 <- getRand(quan_cases_sen, 
                            model = 'enet',
                            rand_num = 10, 
                            num_feat = 5000, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sen_rand_enet_5000")

# unsam_gen 5000 features
sam_gen_rand_enet_5000 <- getRand(quan_cases_sam_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 5000, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sam_gen_rand_enet_5000")

# unsen_gen 5000 features
sen_gen_rand_enet_5000 <- getRand(quan_cases_sen_gen, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 5000, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sen_gen_rand_enet_5000")


###########
# lasso
###########

# ungen 5000 features
ungen_rand_lasso_5000 <- getRand(quan_cases, 
                               model = 'lasso',
                               rand_num = 10, 
                               num_feat = 5000, 
                               exclude = T,
                               union_features = union_features,
                               data_name = "ungen_rand_lasso_5000")

# gen 5000 features
gen_rand_lasso_5000 <- getRand(quan_cases_gen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 5000, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "gen_rand_lasso_5000")

# unsam 5000 features
sam_rand_lasso_5000 <- getRand(quan_cases_sam, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 5000, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sam_rand_lasso_5000")

# unsen 5000 features
sen_rand_lasso_5000 <- getRand(quan_cases_sen, 
                             model = 'lasso',
                             rand_num = 10, 
                             num_feat = 5000, 
                             exclude = T,
                             union_features = union_features,
                             data_name = "sen_rand_lasso_5000")

# unsam_gen 5000 features
sam_gen_rand_lasso_5000 <- getRand(quan_cases_sam_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 5000, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sam_gen_rand_lasso_5000")

# unsen_gen 5000 features
sen_gen_rand_lasso_5000 <- getRand(quan_cases_sen_gen, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 5000, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "sen_gen_rand_lasso_5000")





rand_5000 <- rbind(ungen_rand_rf_5000,
                 gen_rand_rf_5000,
                 sen_rand_rf_5000,
                 sam_rand_rf_5000,
                 sam_gen_rand_rf_5000,
                 sen_gen_rand_rf_5000,
                 ungen_rand_enet_5000,
                 gen_rand_enet_5000,
                 sen_rand_enet_5000,
                 sam_rand_enet_5000,
                 sam_gen_rand_enet_5000,
                 sen_gen_rand_enet_5000,
                 ungen_rand_lasso_5000,
                 gen_rand_lasso_5000,
                 sen_rand_lasso_5000,
                 sam_rand_lasso_5000,
                 sam_gen_rand_lasso_5000,
                 sen_gen_rand_lasso_5000)

rm(ungen_rand_rf_5000,
   gen_rand_rf_5000,
   sen_rand_rf_5000,
   sam_rand_rf_5000,
   sam_gen_rand_rf_5000,
   sen_gen_rand_rf_5000,
   ungen_rand_enet_5000,
   gen_rand_enet_5000,
   sen_rand_enet_5000,
   sam_rand_enet_5000,
   sam_gen_rand_enet_5000,
   sen_gen_rand_enet_5000,
   ungen_rand_lasso_5000,
   gen_rand_lasso_5000,
   sen_rand_lasso_5000,
   sam_rand_lasso_5000,
   sam_gen_rand_lasso_5000,
   sen_gen_rand_lasso_5000)




#############################################################################################################
#10000

###########
# rf
###########

# ungen 10000 features
ungen_rand_rf_10000 <- getRand(quan_cases, 
                              model = 'rf',
                              rand_num = 10, 
                              num_feat = 10000, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "ungen_rand_rf_10000")

# gen 10000 features
gen_rand_rf_10000 <- getRand(quan_cases_gen, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 10000, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "gen_rand_rf_10000")

# unsam 10000 features
sam_rand_rf_10000 <- getRand(quan_cases_sam, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 10000, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sam_rand_rf_10000")

# unsen 10000 features
sen_rand_rf_10000 <- getRand(quan_cases_sen, 
                            model = 'rf',
                            rand_num = 10, 
                            num_feat = 10000, 
                            exclude = T,
                            union_features = union_features,
                            data_name = "sen_rand_rf_10000")

# unsam_gen 10000 features
sam_gen_rand_rf_10000 <- getRand(quan_cases_sam_gen, 
                                model = 'rf',
                                rand_num = 10, 
                                num_feat = 10000, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sam_gen_rand_rf_10000")

# unsen_gen 10000 features
sen_gen_rand_rf_10000 <- getRand(quan_cases_sen_gen, 
                                model = 'rf',
                                rand_num = 10, 
                                num_feat = 10000, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "sen_gen_rand_rf_10000")




###########
# enet
###########

# ungen 10000 features
ungen_rand_enet_10000 <- getRand(quan_cases, 
                                model = 'enet',
                                rand_num = 10, 
                                num_feat = 10000, 
                                exclude = T,
                                union_features = union_features,
                                data_name = "ungen_rand_enet_10000")

# gen 10000 features
gen_rand_enet_10000 <- getRand(quan_cases_gen, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 10000, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "gen_rand_enet_10000")

# unsam 10000 features
sam_rand_enet_10000 <- getRand(quan_cases_sam, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 10000, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sam_rand_enet_10000")

# unsen 10000 features
sen_rand_enet_10000 <- getRand(quan_cases_sen, 
                              model = 'enet',
                              rand_num = 10, 
                              num_feat = 10000, 
                              exclude = T,
                              union_features = union_features,
                              data_name = "sen_rand_enet_10000")

# unsam_gen 10000 features
sam_gen_rand_enet_10000 <- getRand(quan_cases_sam_gen, 
                                  model = 'enet',
                                  rand_num = 10, 
                                  num_feat = 10000, 
                                  exclude = T,
                                  union_features = union_features,
                                  data_name = "sam_gen_rand_enet_10000")

# unsen_gen 10000 features
sen_gen_rand_enet_10000 <- getRand(quan_cases_sen_gen, 
                                  model = 'enet',
                                  rand_num = 10, 
                                  num_feat = 10000, 
                                  exclude = T,
                                  union_features = union_features,
                                  data_name = "sen_gen_rand_enet_10000")


###########
# lasso
###########

# ungen 10000 features
ungen_rand_lasso_10000 <- getRand(quan_cases, 
                                 model = 'lasso',
                                 rand_num = 10, 
                                 num_feat = 10000, 
                                 exclude = T,
                                 union_features = union_features,
                                 data_name = "ungen_rand_lasso_10000")

# gen 10000 features
gen_rand_lasso_10000 <- getRand(quan_cases_gen, 
                               model = 'lasso',
                               rand_num = 10, 
                               num_feat = 10000, 
                               exclude = T,
                               union_features = union_features,
                               data_name = "gen_rand_lasso_10000")

# unsam 10000 features
sam_rand_lasso_10000 <- getRand(quan_cases_sam, 
                               model = 'lasso',
                               rand_num = 10, 
                               num_feat = 10000, 
                               exclude = T,
                               union_features = union_features,
                               data_name = "sam_rand_lasso_10000")

# unsen 10000 features
sen_rand_lasso_10000 <- getRand(quan_cases_sen, 
                               model = 'lasso',
                               rand_num = 10, 
                               num_feat = 10000, 
                               exclude = T,
                               union_features = union_features,
                               data_name = "sen_rand_lasso_10000")

# unsam_gen 10000 features
sam_gen_rand_lasso_10000 <- getRand(quan_cases_sam_gen, 
                                   model = 'lasso',
                                   rand_num = 10, 
                                   num_feat = 10000, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "sam_gen_rand_lasso_10000")

# unsen_gen 10000 features
sen_gen_rand_lasso_10000 <- getRand(quan_cases_sen_gen, 
                                   model = 'lasso',
                                   rand_num = 10, 
                                   num_feat = 10000, 
                                   exclude = T,
                                   union_features = union_features,
                                   data_name = "sen_gen_rand_lasso_10000")


rand_10000 <- rbind(ungen_rand_rf_10000,
                   gen_rand_rf_10000,
                   sen_rand_rf_10000,
                   sam_rand_rf_10000,
                   sam_gen_rand_rf_10000,
                   sen_gen_rand_rf_10000,
                   ungen_rand_enet_10000,
                   gen_rand_enet_10000,
                   sen_rand_enet_10000,
                   sam_rand_enet_10000,
                   sam_gen_rand_enet_10000,
                   sen_gen_rand_enet_10000,
                   ungen_rand_lasso_10000,
                   gen_rand_lasso_10000,
                   sen_rand_lasso_10000,
                   sam_rand_lasso_10000,
                   sam_gen_rand_lasso_10000,
                   sen_gen_rand_lasso_10000)

rm(ungen_rand_rf_10000,
   gen_rand_rf_10000,
   sen_rand_rf_10000,
   sam_rand_rf_10000,
   sam_gen_rand_rf_10000,
   sen_gen_rand_rf_10000,
   ungen_rand_enet_10000,
   gen_rand_enet_10000,
   sen_rand_enet_10000,
   sam_rand_enet_10000,
   sam_gen_rand_enet_10000,
   sen_gen_rand_enet_10000,
   ungen_rand_lasso_10000,
   gen_rand_lasso_10000,
   sen_rand_lasso_10000,
   sam_rand_lasso_10000,
   sam_gen_rand_lasso_10000,
   sen_gen_rand_lasso_10000)



##########
# combine data
##########

rand <- rbind(rand_10, rand_20, rand_30, rand_50, rand_100, rand_300, rand_500,
              rand_1000, rand_5000, rand_10000)

saveRDS(rand, paste0(model_data, '/rand.rda'))
