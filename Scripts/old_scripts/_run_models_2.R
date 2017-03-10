####### Script will run models

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

##########
# source model_functions to get functions to run models 
##########
source(paste0(scripts_folder, '/predict_age/functions.R'))

# ##########
# # set parameters 
# ##########
# data_thresholds <- c(48, 60, 72)
# p53 <- c('Mut', 'WT')


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

########################################################################
# full data

even_full_rf <- runModels(quan_cases, 
                       model = 'rf',
                       rand = F,
                       bump_hunter = F,
                       bump_hunter_data)

even_full_table_rf <- extractResults(even_full_rf, data_name = 'even full rf', regularize = F)

even_full_enet <- runModels(quan_cases, 
                          model = 'enet',
                          rand = F,
                          bump_hunter = F,
                          bump_hunter_data)

even_full_table_enet <- extractResults(even_full_enet, data_name = 'even full enet', regularize = F)

even_full_lasso <- runModels(quan_cases, 
                          model = 'lasso',
                          rand = F,
                          bump_hunter = F,
                          bump_hunter_data)

even_full_table_lasso <- extractResults(even_full_lasso, data_name = 'even full lasso', regularize = F)

# get table 
even_full <- rbind(even_full_table_rf,
                   even_full_table_enet,
                   even_full_table_lasso)


rm(even_full_table_rf,
   even_full_table_enet,
   even_full_table_lasso)


# full data gen
even_full_gen_rf <- runModels(quan_cases_gen, 
                          model = 'rf',
                          rand = F,
                          bump_hunter = F,
                          bump_hunter_data)

even_full_gen_table_rf <- extractResults(even_full_gen_rf, data_name = 'even gen_full rf', regularize = F)

even_full_gen_enet <- runModels(quan_cases_gen, 
                            model = 'enet',
                            rand = F,
                            bump_hunter = F,
                            bump_hunter_data)

even_full_gen_table_enet <- extractResults(even_full_gen_enet, data_name = 'even gen_full enet', regularize = F)

even_full_gen_lasso <- runModels(quan_cases_gen, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = F,
                             bump_hunter_data)

even_full_gen_table_lasso <- extractResults(even_full_gen_lasso, data_name = 'even gen_full lasso', regularize = F)

# get table 
even_full_gen <- rbind(even_full_gen_table_rf,
                   even_full_gen_table_enet,
                   even_full_gen_table_lasso)


rm(even_full_gen_table_rf,
   even_full_gen_table_enet,
   even_full_gen_table_lasso)

# full data sam
even_full_sam_rf <- runModels(quan_cases_sam, 
                              model = 'rf',
                              rand = F,
                              bump_hunter = F,
                              bump_hunter_data)

even_full_sam_table_rf <- extractResults(even_full_sam_rf, data_name = 'even sam_full rf', regularize = F)

even_full_sam_enet <- runModels(quan_cases_sam, 
                                model = 'enet',
                                rand = F,
                                bump_hunter = F,
                                bump_hunter_data)

even_full_sam_table_enet <- extractResults(even_full_sam_enet, data_name = 'even sam_full enet', regularize = F)

even_full_sam_lasso <- runModels(quan_cases_sam, 
                                 model = 'lasso',
                                 rand = F,
                                 bump_hunter = F,
                                 bump_hunter_data)

even_full_sam_table_lasso <- extractResults(even_full_sam_lasso, data_name = 'even sam_full lasso', regularize = F)

#here
# get table 
even_full_sam <- rbind(even_full_sam_table_rf,
                       even_full_sam_table_enet,
                       even_full_sam_table_lasso)


rm(even_full_sam_table_rf,
   even_full_sam_table_enet,
   even_full_sam_table_lasso)


# full data sen
even_full_sen_rf <- runModels(quan_cases_sen, 
                              model = 'rf',
                              rand = F,
                              bump_hunter = F,
                              bump_hunter_data)

even_full_sen_table_rf <- extractResults(even_full_sen_rf, data_name = 'even sen_full rf', regularize = F)

even_full_sen_enet <- runModels(quan_cases_sen, 
                                model = 'enet',
                                rand = F,
                                bump_hunter = F,
                                bump_hunter_data)

even_full_sen_table_enet <- extractResults(even_full_sen_enet, data_name = 'even sen_full enet', regularize = F)

even_full_sen_lasso <- runModels(quan_cases_sen, 
                                 model = 'lasso',
                                 rand = F,
                                 bump_hunter = F,
                                 bump_hunter_data)

even_full_sen_table_lasso <- extractResults(even_full_sen_lasso, data_name = 'even sen_full lasso', regularize = F)

#here
# get table 
even_full_sen <- rbind(even_full_sen_table_rf,
                       even_full_sen_table_enet,
                       even_full_sen_table_lasso)


rm(even_full_sen_table_rf,
   even_full_sen_table_enet,
   even_full_sen_table_lasso)



# full data sam_gen
even_full_sam_gen_rf <- runModels(quan_cases_sam_gen, 
                              model = 'rf',
                              rand = F,
                              bump_hunter = F,
                              bump_hunter_data)

even_full_sam_gen_table_rf <- extractResults(even_full_sam_gen_rf, data_name = 'even sam_gen_full rf', regularize = F)

even_full_sam_gen_enet <- runModels(quan_cases_sam_gen, 
                                model = 'enet',
                                rand = F,
                                bump_hunter = F,
                                bump_hunter_data)

even_full_sam_gen_table_enet <- extractResults(even_full_sam_gen_enet, data_name = 'even sam_gen_full enet', regularize = F)

even_full_sam_gen_lasso <- runModels(quan_cases_sam_gen, 
                                 model = 'lasso',
                                 rand = F,
                                 bump_hunter = F,
                                 bump_hunter_data)

even_full_sam_gen_table_lasso <- extractResults(even_full_sam_gen_lasso, data_name = 'even sam_gen_full lasso', regularize = F)

# get table 
even_full_sam_gen <- rbind(even_full_sam_gen_table_rf,
                       even_full_sam_gen_table_enet,
                       even_full_sam_gen_table_lasso)


rm(even_full_sam_gen_table_rf,
   even_full_sam_gen_table_enet,
   even_full_sam_gen_table_lasso)



# full data sen_gen
even_full_sen_gen_rf <- runModels(quan_cases_sen_gen, 
                              model = 'rf',
                              rand = F,
                              bump_hunter = F,
                              bump_hunter_data)

even_full_sen_gen_table_rf <- extractResults(even_full_sen_gen_rf, data_name = 'even sen_gen_full rf', regularize = F)

even_full_sen_gen_enet <- runModels(quan_cases_sen_gen, 
                                model = 'enet',
                                rand = F,
                                bump_hunter = F,
                                bump_hunter_data)

even_full_sen_gen_table_enet <- extractResults(even_full_sen_gen_enet, data_name = 'even sen_gen_full enet', regularize = F)

even_full_sen_gen_lasso <- runModels(quan_cases_sen_gen, 
                                 model = 'lasso',
                                 rand = F,
                                 bump_hunter = F,
                                 bump_hunter_data)

even_full_sen_gen_table_lasso <- extractResults(even_full_sen_gen_lasso, data_name = 'even sen_gen_full lasso', regularize = F)

# get table 
even_full_sen_gen <- rbind(even_full_sen_gen_table_rf,
                       even_full_sen_gen_table_enet,
                       even_full_sen_gen_table_lasso)


rm(even_full_sen_gen_table_rf,
   even_full_sen_gen_table_enet,
   even_full_sen_gen_table_lasso)

#########
# put full together
#########

even_full <- rbind(even_full, even_full_gen, even_full_sam, even_full_sen,
                   even_full_sam_gen, even_full_sen_gen)

saveRDS(even_full, paste0(model_data, '/even_full.rds'))

########################################################################3
##########
# run models
##########

##################################
####### not significant

#####not gen

#### 10
even_10 <- runModels(quan_cases, 
                     model = 'rf',
                     rand = F,
                     bump_hunter = T,
                     bump_hunter_data = quan_even_10)

even_10_table <- extractResults(even_10, data_name = 'even_10 rf', regularize = F)

#### 20
even_20 <- runModels(quan_cases, 
                     model = 'rf',
                     rand = F,
                     bump_hunter = T,
                     bump_hunter_data = quan_even_20)

even_20_table <- extractResults(even_10, data_name = 'even_20 rf', regularize = F)


#### 30
even_30 <- runModels(quan_cases, 
                     model = 'rf',
                     rand = F,
                     bump_hunter = T,
                     bump_hunter_data = quan_even_30)

even_30_table <- extractResults(even_30, data_name = 'even_30 rf', regularize = F)


##### gen

#### 10
even_gen_10 <- runModels(quan_cases_gen, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_gen_10)

even_gen_10_table <- extractResults(even_gen_10, data_name = 'even gen_10 rf', regularize = F)


#### 20
even_gen_20 <- runModels(quan_cases_gen, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_gen_20)

even_gen_20_table <- extractResults(even_gen_20, data_name = 'even gen_20 rf', regularize = F)

#### 30
even_gen_30 <- runModels(quan_cases_gen, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_gen_30)

even_gen_30_table <- extractResults(even_gen_30, data_name = 'even gen_30 rf', regularize = F)

##### sen

#### 10
even_sen_10 <- runModels(quan_cases_sen, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_10)

even_sen_10_table <- extractResults(even_sen_10, data_name = 'even sen_10 rf', regularize = F)


#### 20
even_sen_20 <- runModels(quan_cases_sen, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_20)

even_sen_20_table <- extractResults(even_sen_20, data_name = 'even sen_20 rf', regularize = F)

#### 30
even_sen_30 <- runModels(quan_cases_sen, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_30)

even_sen_30_table <- extractResults(even_sen_30, data_name = 'even sen_30 rf', regularize = F)


##### sam

#### 10
even_sam_10 <- runModels(quan_cases_sam, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_10)

even_sam_10_table <- extractResults(even_sam_10, data_name = 'even sam_10 rf', regularize = F)


#### 20
even_sam_20 <- runModels(quan_cases_sam, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_20)

even_sam_20_table <- extractResults(even_sam_20, data_name = 'even sam_20 rf', regularize = F)

#### 30
even_sam_30 <- runModels(quan_cases_sam, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_30)

even_sam_30_table <- extractResults(even_sam_30, data_name = 'even sam_30 rf', regularize = F)


##### gen sam

#### 10
even_gen_sam_10 <- runModels(quan_cases_sam_gen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sam_10)

even_gen_sam_10_table <- extractResults(even_gen_sam_10, data_name = 'even gen sam_10 rf', regularize = F)


#### 20
even_gen_sam_20 <- runModels(quan_cases_sam_gen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sam_20)

even_gen_sam_20_table <- extractResults(even_gen_sam_20, data_name = 'even gen sam_20 rf', regularize = F)

#### 30
even_gen_sam_30 <- runModels(quan_cases_sam_gen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sam_30)

even_gen_sam_30_table <- extractResults(even_gen_sam_30, data_name = 'even gen sam_30 rf', regularize = F)



##### gen sen

#### 10
even_gen_sen_10 <- runModels(quan_cases_sen_gen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sen_10)

even_gen_sen_10_table <- extractResults(even_gen_sen_10, data_name = 'even gen sen_10 rf', regularize = F)


#### 20
even_gen_sen_20 <- runModels(quan_cases_sen_gen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sen_20)

even_gen_sen_20_table <- extractResults(even_gen_sen_20, data_name = 'even gen sen_20 rf', regularize = F)

#### 30
even_gen_sen_30 <- runModels(quan_cases_sen_gen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sen_30)

even_gen_sen_30_table <- extractResults(even_gen_sen_30, data_name = 'even gen sen_30 rf', regularize = F)

###########
# get table for even not sig
###########
even_not_sig_rf <- rbind(even_10_table, even_20_table, even_30_table,
                      even_gen_10_table, even_gen_20_table, even_gen_30_table,
                      even_sen_10_table, even_sen_20_table, even_sen_30_table,
                      even_sam_10_table, even_sam_20_table, even_sam_30_table,
                      even_gen_sen_10_table, even_gen_sen_20_table, even_gen_sen_30_table,
                      even_gen_sam_10_table, even_gen_sam_20_table, even_gen_sam_30_table)

rm(even_10_table, even_20_table, even_30_table,
   even_gen_10_table, even_gen_20_table, even_gen_30_table,
   even_sen_10_table, even_sen_20_table, even_sen_30_table,
   even_sam_10_table, even_sam_20_table, even_sam_30_table,
   even_gen_sen_10_table, even_gen_sen_20_table, even_gen_sen_30_table,
   even_gen_sam_10_table, even_gen_sam_20_table, even_gen_sam_30_table)


##################################
#######significant

#####not gen

#### 10
even_sig_10 <- runModels(quan_cases, 
                     model = 'rf',
                     rand = F,
                     bump_hunter = T,
                     bump_hunter_data = quan_even_sig_10)

even_sig_10_table <- extractResults(even_sig_10, data_name = 'even sig_10 rf', regularize = F)

#### 20
even_sig_20 <- runModels(quan_cases, 
                     model = 'rf',
                     rand = F,
                     bump_hunter = T,
                     bump_hunter_data = quan_even_sig_20)

even_sig_20_table <- extractResults(even_sig_20, data_name = 'even sig_20 rf', regularize = F)


# #### 30
# even_sig_30 <- runModels(quan_cases, 
#                      model = 'rf',
#                      rand = F,
#                      bump_hunter = T,
#                      bump_hunter_data = quan_even_sig_30)
# 
# even_sig_30_table <- extractResults(even_sig_30, data_name = 'even sig_30 rf', regularize = F)


##### gen

#### 10
even_gen_sig_10 <- runModels(quan_cases_gen, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_gen_sig_10)

even_gen_sig_10_table <- extractResults(even_gen_sig_10, data_name = 'even gen sig_10 rf', regularize = F)


#### 20
even_gen_sig_20 <- runModels(quan_cases_gen, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_gen_sig_20)

even_gen_20_sig_table <- extractResults(even_gen_sig_20, data_name = 'even gen sig_20 rf', regularize = F)

# #### 30
# even_gen_sig_30 <- runModels(quan_cases_gen, 
#                          model = 'rf',
#                          rand = F,
#                          bump_hunter = T,
#                          bump_hunter_data = quan_even_gen_sig_30)
# 
# even_gen_sig_30_table <- extractResults(even_gen_sig_30, data_name = 'even gen sig_30 rf', regularize = F)

##### sen

#### 10
even_sen_sig_10 <- runModels(quan_cases_sen, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_sig_10)

even_sen_sig_10_table <- extractResults(even_sen_sig_10, data_name = 'even sen sig_10 rf', regularize = F)


#### 20
even_sen_sig_20 <- runModels(quan_cases_sen, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_sig_20)

even_sen_sig_20_table <- extractResults(even_sen_sig_20, data_name = 'even sen sig_20 rf', regularize = F)

# #### 30
# even_sen_sig_30 <- runModels(quan_cases_sen, 
#                          model = 'rf',
#                          rand = F,
#                          bump_hunter = T,
#                          bump_hunter_data = quan_even_sig_30)
# 
# even_sen_sig_30_table <- extractResults(even_sen_sig_30, data_name = 'even sen sig_30 rf', regularize = F)


##### sam

#### 10
even_sam_sig_10 <- runModels(quan_cases_sam, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_sig_10)

even_sam_sig_10_table <- extractResults(even_sam_sig_10, data_name = 'even sam sig_10 rf', regularize = F)


#### 20
even_sam_sig_20 <- runModels(quan_cases_sam, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_sig_20)

even_sam_sig_20_table <- extractResults(even_sam_sig_20, data_name = 'even sam sig_20 rf', regularize = F)

# #### 30
# even_sam_sig_30 <- runModels(quan_cases_sam, 
#                          model = 'rf',
#                          rand = F,
#                          bump_hunter = T,
#                          bump_hunter_data = quan_even_sig_30)
# 
# even_sam_sig_30_table <- extractResults(even_sam_sig_30, data_name = 'even sam sig_30 rf', regularize = F)


##### gen sam

#### 10
even_gen_sam_sig_10 <- runModels(quan_cases_sam_gen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sam_sig_10)

even_gen_sam_10_sig_table <- extractResults(even_gen_sam_sig_10, data_name = 'even gen sam sig_10 rf', regularize = F)


#### 20
even_gen_sam_sig_20 <- runModels(quan_cases_sam_gen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sam_sig_20)

even_gen_sam_sig_20_table <- extractResults(even_gen_sam_sig_20, data_name = 'even gen sam sig_20 rf', regularize = F)

# #### 30
# even_gen_sam_sig_30 <- runModels(quan_cases_sam_gen, 
#                              model = 'rf',
#                              rand = F,
#                              bump_hunter = T,
#                              bump_hunter_data = quan_even_gen_sam_sig_30)
# 
# even_gen_sam_sig_30_table <- extractResults(even_gen_sam_sig_30, data_name = 'even gen sam sig_30 rf', regularize = F)



##### gen sen

#### 10
even_gen_sen_sig_10 <- runModels(quan_cases_sen_gen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sen_sig_10)

even_gen_sen_sig_10_table <- extractResults(even_gen_sen_sig_10, data_name = 'even gen sen sig_10 rf', regularize = F)


# #### 20
# even_gen_sen_sig_20 <- runModels(quan_cases_sen_gen, 
#                              model = 'rf',
#                              rand = F,
#                              bump_hunter = T,
#                              bump_hunter_data = quan_even_gen_sen_sig_20)
# 
# even_gen_sen_sig_20_table <- extractResults(even_gen_sen_sig_20, data_name = 'even gen sen sig_20 rf', regularize = F)

# #### 30
# even_gen_sen_sig_30 <- runModels(quan_cases_sen_gen, 
#                              model = 'rf',
#                              rand = F,
#                              bump_hunter = T,
#                              bump_hunter_data = quan_even_gen_sen_sig_30)
# 
# even_gen_sen_sig_30_table <- extractResults(even_gen_sen_sig_30, data_name = 'even gen sen sig_30 rf', regularize = F)

###########
# get table for even not sig
###########
even_sig_rf <- rbind(even_sig_10_table, even_sig_20_table,
                      even_gen_sig_10_table,
                      even_sen_sig_10_table, even_sen_sig_20_table,
                      even_sam_sig_10_table, even_sam_sig_20_table,
                      even_gen_sen_sig_10_table,even_gen_sam_sig_20_table)

rm(even_sig_10_table, even_sig_20_table,
   even_gen_sig_10_table,
   even_sen_sig_10_table, even_sen_sig_20_table,
   even_sam_sig_10_table, even_sam_sig_20_table,
   even_gen_sen_sig_10_table,even_gen_sam_sig_20_table)


##################################
#######fwer

#####not gen

#### 10
even_fwer_10 <- runModels(quan_cases, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_fwer_10)

even_fwer_10_table <- extractResults(even_fwer_10, data_name = 'even fwer 10 rf', regularize = F)

#### 20
even_fwer_20 <- runModels(quan_cases, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_fwer_20)

even_fwer_20_table <- extractResults(even_fwer_20, data_name = 'even fwer_20 rf', regularize = F)


#### 30
even_fwer_30 <- runModels(quan_cases, 
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_fwer_30)

even_fwer_30_table <- extractResults(even_fwer_30, data_name = 'even fwer_30 rf', regularize = F)


##### gen

#### 10
even_gen_fwer_10 <- runModels(quan_cases_gen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_fwer_10)

even_gen_fwer_10_table <- extractResults(even_gen_fwer_10, data_name = 'even gen fwer_10 rf', regularize = F)


#### 20
even_gen_fwer_20 <- runModels(quan_cases_gen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_fwer_20)

even_gen_fwer_20_table <- extractResults(even_gen_fwer_20, data_name = 'even gen fwer_20 rf', regularize = F)

#### 30
even_gen_fwer_30 <- runModels(quan_cases_gen,
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_gen_fwer_30)

even_gen_fwer_30_table <- extractResults(even_gen_fwer_30, data_name = 'even gen fwer_30 rf', regularize = F)

##### sen

#### 10
even_sen_fwer_10 <- runModels(quan_cases_sen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_fwer_10)

even_sen_fwer_10_table <- extractResults(even_sen_fwer_10, data_name = 'even sen fwer_10 rf', regularize = F)


#### 20
even_sen_fwer_20 <- runModels(quan_cases_sen, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_fwer_20)

even_sen_fwer_20_table <- extractResults(even_sen_fwer_20, data_name = 'even sen fwer_20 rf', regularize = F)

#### 30
even_sen_fwer_30 <- runModels(quan_cases_sen,
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_fwer_30)

even_sen_fwer_30_table <- extractResults(even_sen_fwer_30, data_name = 'even sen fwer_30 rf', regularize = F)


##### sam

#### 10
even_sam_fwer_10 <- runModels(quan_cases_sam, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_fwer_10)

even_sam_fwer_10_table <- extractResults(even_sam_fwer_10, data_name = 'even sam fwer_10 rf', regularize = F)


#### 20
even_sam_fwer_20 <- runModels(quan_cases_sam, 
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_fwer_20)

even_sam_fwer_20_table <- extractResults(even_sam_fwer_20, data_name = 'even sam fwer_20 rf', regularize = F)

#### 30
even_sam_fwer_30 <- runModels(quan_cases_sam,
                         model = 'rf',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_fwer_30)

even_sam_fwer_30_table <- extractResults(even_sam_fwer_30, data_name = 'even sam fwer_30 rf', regularize = F)


##### gen sam

#### 10
even_gen_sam_fwer_10 <- runModels(quan_cases_sam_gen, 
                                 model = 'rf',
                                 rand = F,
                                 bump_hunter = T,
                                 bump_hunter_data = quan_even_gen_sam_fwer_10)

even_gen_sam_fwer_10_table <- extractResults(even_gen_sam_fwer_10, data_name = 'even gen sam fwer_10 rf', regularize = F)


#### 20
even_gen_sam_fwer_20 <- runModels(quan_cases_sam_gen, 
                                 model = 'rf',
                                 rand = F,
                                 bump_hunter = T,
                                 bump_hunter_data = quan_even_gen_sam_fwer_20)

even_gen_sam_fwer_20_table <- extractResults(even_gen_sam_fwer_20, data_name = 'even gen sam fwer_20 rf', regularize = F)

#### 30
even_gen_sam_fwer_30 <- runModels(quan_cases_sam_gen,
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sam_fwer_30)

even_gen_sam_fwer_30_table <- extractResults(even_gen_sam_fwer_30, data_name = 'even gen sam fwer_30 rf', regularize = F)



##### gen sen

#### 10
even_gen_sen_fwer_10 <- runModels(quan_cases_sen_gen, 
                                 model = 'rf',
                                 rand = F,
                                 bump_hunter = T,
                                 bump_hunter_data = quan_even_gen_sen_fwer_10)

even_gen_sen_fwer_10_table <- extractResults(even_gen_sen_fwer_10, data_name = 'even gen sen fwer_10 rf', regularize = F)


#### 20
even_gen_sen_fwer_20 <- runModels(quan_cases_sen_gen,
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sen_fwer_20)

even_gen_sen_fwer_20_table <- extractResults(even_gen_sen_fwer_20, data_name = 'even gen sen fwer_20 rf', regularize = F)

#### 30
even_gen_sen_fwer_30 <- runModels(quan_cases_sen_gen,
                             model = 'rf',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sen_fwer_30)

even_gen_sen_fwer_30_table <- extractResults(even_gen_sen_fwer_30, data_name = 'even gen sen fwer_30 rf', regularize = F)

###########
# get table for even not fwer
###########
even_fwer_rf <- rbind(even_fwer_10_table, even_fwer_20_table, even_fwer_30_table,
                   even_gen_fwer_10_table, even_gen_fwer_20_table, even_gen_fwer_30_table,
                   even_sen_fwer_10_table, even_sen_fwer_20_table,even_sen_fwer_30_table,
                   even_sam_fwer_10_table, even_sam_fwer_20_table, even_sam_fwer_30_table,
                   even_gen_sen_fwer_10_table, even_gen_sen_fwer_20_table, even_gen_sen_fwer_30_table,
                   even_gen_sam_fwer_10_table, even_gen_sam_fwer_20_table, even_gen_sam_fwer_30_table)

rm(even_fwer_10_table, even_fwer_20_table, even_fwer_30_table,
   even_gen_fwer_10_table, even_gen_fwer_20_table, even_gen_fwer_30_table,
   even_sen_fwer_10_table, even_sen_fwer_20_table,even_sen_fwer_30_table,
   even_sam_fwer_10_table, even_sam_fwer_20_table, even_sam_fwer_30_table,
   even_gen_sen_fwer_10_table, even_gen_sen_fwer_20_table, even_gen_sen_fwer_30_table,
   even_gen_sam_fwer_10_table, even_gen_sam_fwer_20_table, even_gen_sam_fwer_30_table)



########################################################################3
##########
# run models
##########

##################################
####### not significant
#####not gen

#### 10
even_10 <- runModels(quan_cases, 
                     model = 'enet',
                     rand = F,
                     bump_hunter = T,
                     bump_hunter_data = quan_even_10)

even_10_table <- extractResults(even_10, data_name = 'even_10 enet', regularize = T)

#### 20
even_20 <- runModels(quan_cases, 
                     model = 'enet',
                     rand = F,
                     bump_hunter = T,
                     bump_hunter_data = quan_even_20)

even_20_table <- extractResults(even_10, data_name = 'even_20 enet', regularize = T)


#### 30
even_30 <- runModels(quan_cases, 
                     model = 'enet',
                     rand = F,
                     bump_hunter = T,
                     bump_hunter_data = quan_even_30)

even_30_table <- extractResults(even_30, data_name = 'even_30 enet', regularize = T)


##### gen

#### 10
even_gen_10 <- runModels(quan_cases_gen, 
                         model = 'enet',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_gen_10)

even_gen_10_table <- extractResults(even_gen_10, data_name = 'even gen_10 enet', regularize = T)


#### 20
even_gen_20 <- runModels(quan_cases_gen, 
                         model = 'enet',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_gen_20)

even_gen_20_table <- extractResults(even_gen_20, data_name = 'even gen_20 enet', regularize = T)

#### 30
even_gen_30 <- runModels(quan_cases_gen, 
                         model = 'enet',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_gen_30)

even_gen_30_table <- extractResults(even_gen_30, data_name = 'even gen_30 enet', regularize = T)

##### sen

#### 10
even_sen_10 <- runModels(quan_cases_sen, 
                         model = 'enet',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_10)

even_sen_10_table <- extractResults(even_sen_10, data_name = 'even sen_10 enet', regularize = T)


#### 20
even_sen_20 <- runModels(quan_cases_sen, 
                         model = 'enet',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_20)

even_sen_20_table <- extractResults(even_sen_20, data_name = 'even sen_20 enet', regularize = T)

#### 30
even_sen_30 <- runModels(quan_cases_sen, 
                         model = 'enet',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_30)

even_sen_30_table <- extractResults(even_sen_30, data_name = 'even sen_30 enet', regularize = T)


##### sam

#### 10
even_sam_10 <- runModels(quan_cases_sam, 
                         model = 'enet',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_10)

even_sam_10_table <- extractResults(even_sam_10, data_name = 'even sam_10 enet', regularize = T)


#### 20
even_sam_20 <- runModels(quan_cases_sam, 
                         model = 'enet',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_20)

even_sam_20_table <- extractResults(even_sam_20, data_name = 'even sam_20 enet', regularize = T)

#### 30
even_sam_30 <- runModels(quan_cases_sam, 
                         model = 'enet',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_30)

even_sam_30_table <- extractResults(even_sam_30, data_name = 'even sam_30 enet', regularize = T)


##### gen sam

#### 10
even_gen_sam_10 <- runModels(quan_cases_sam_gen, 
                             model = 'enet',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sam_10)

even_gen_sam_10_table <- extractResults(even_gen_sam_10, data_name = 'even gen sam_10 enet', regularize = T)


#### 20
even_gen_sam_20 <- runModels(quan_cases_sam_gen, 
                             model = 'enet',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sam_20)

even_gen_sam_20_table <- extractResults(even_gen_sam_20, data_name = 'even gen sam_20 enet', regularize = T)

#### 30
even_gen_sam_30 <- runModels(quan_cases_sam_gen, 
                             model = 'enet',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sam_30)

even_gen_sam_30_table <- extractResults(even_gen_sam_30, data_name = 'even gen sam_30 enet', regularize = T)



##### gen sen

#### 10
even_gen_sen_10 <- runModels(quan_cases_sen_gen, 
                             model = 'enet',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sen_10)

even_gen_sen_10_table <- extractResults(even_gen_sen_10, data_name = 'even gen sen_10 enet', regularize = T)


#### 20
even_gen_sen_20 <- runModels(quan_cases_sen_gen, 
                             model = 'enet',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sen_20)

even_gen_sen_20_table <- extractResults(even_gen_sen_20, data_name = 'even gen sen_20 enet', regularize = T)

#### 30
even_gen_sen_30 <- runModels(quan_cases_sen_gen, 
                             model = 'enet',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sen_30)

even_gen_sen_30_table <- extractResults(even_gen_sen_30, data_name = 'even gen sen_30 enet', regularize = T)

###########
# get table for even not sig
###########
even_not_sig_enet <- rbind(even_10_table, even_20_table, even_30_table,
                         even_gen_10_table, even_gen_20_table, even_gen_30_table,
                         even_sen_10_table, even_sen_20_table, even_sen_30_table,
                         even_sam_10_table, even_sam_20_table, even_sam_30_table,
                         even_gen_sen_10_table, even_gen_sen_20_table, even_gen_sen_30_table,
                         even_gen_sam_10_table, even_gen_sam_20_table, even_gen_sam_30_table)

rm(even_10_table, even_20_table, even_30_table,
   even_gen_10_table, even_gen_20_table, even_gen_30_table,
   even_sen_10_table, even_sen_20_table, even_sen_30_table,
   even_sam_10_table, even_sam_20_table, even_sam_30_table,
   even_gen_sen_10_table, even_gen_sen_20_table, even_gen_sen_30_table,
   even_gen_sam_10_table, even_gen_sam_20_table, even_gen_sam_30_table)


##################################
#######significant

#####not gen

#### 10
even_sig_10 <- runModels(quan_cases, 
                         model = 'enet',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_sig_10)

even_sig_10_table <- extractResults(even_sig_10, data_name = 'even sig_10 enet', regularize = T)

#### 20
even_sig_20 <- runModels(quan_cases, 
                         model = 'enet',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_sig_20)

even_sig_20_table <- extractResults(even_sig_20, data_name = 'even sig_20 enet', regularize = T)


# #### 30
# even_sig_30 <- runModels(quan_cases, 
#                          model = 'enet',
#                          rand = F,
#                          bump_hunter = T,
#                          bump_hunter_data = quan_even_sig_30)
# 
# even_sig_30_table <- extractResults(even_sig_30, data_name = 'even sig_30 enet', regularize = T)


##### gen

#### 10
even_gen_sig_10 <- runModels(quan_cases_gen, 
                             model = 'enet',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sig_10)

even_gen_sig_10_table <- extractResults(even_gen_sig_10, data_name = 'even gen sig_10 enet', regularize = T)


#### 20
even_gen_sig_20 <- runModels(quan_cases_gen, 
                             model = 'enet',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sig_20)

even_gen_20_sig_table <- extractResults(even_gen_sig_20, data_name = 'even gen sig_20 enet', regularize = T)

# #### 30
# even_gen_sig_30 <- runModels(quan_cases_gen, 
#                          model = 'enet',
#                          rand = F,
#                          bump_hunter = T,
#                          bump_hunter_data = quan_even_gen_sig_30)
# 
# even_gen_sig_30_table <- extractResults(even_gen_sig_30, data_name = 'even gen sig_30 enet', regularize = T)

##### sen

#### 10
even_sen_sig_10 <- runModels(quan_cases_sen, 
                             model = 'enet',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_sig_10)

even_sen_sig_10_table <- extractResults(even_sen_sig_10, data_name = 'even sen sig_10 enet', regularize = T)


#### 20
even_sen_sig_20 <- runModels(quan_cases_sen, 
                             model = 'enet',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_sig_20)

even_sen_sig_20_table <- extractResults(even_sen_sig_20, data_name = 'even sen sig_20 enet', regularize = T)

# #### 30
# even_sen_sig_30 <- runModels(quan_cases_sen, 
#                          model = 'enet',
#                          rand = F,
#                          bump_hunter = T,
#                          bump_hunter_data = quan_even_sig_30)
# 
# even_sen_sig_30_table <- extractResults(even_sen_sig_30, data_name = 'even sen sig_30 enet', regularize = T)


##### sam

#### 10
even_sam_sig_10 <- runModels(quan_cases_sam, 
                             model = 'enet',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_sig_10)

even_sam_sig_10_table <- extractResults(even_sam_sig_10, data_name = 'even sam sig_10 enet', regularize = T)


#### 20
even_sam_sig_20 <- runModels(quan_cases_sam, 
                             model = 'enet',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_sig_20)

even_sam_sig_20_table <- extractResults(even_sam_sig_20, data_name = 'even sam sig_20 enet', regularize = T)

# #### 30
# even_sam_sig_30 <- runModels(quan_cases_sam, 
#                          model = 'enet',
#                          rand = F,
#                          bump_hunter = T,
#                          bump_hunter_data = quan_even_sig_30)
# 
# even_sam_sig_30_table <- extractResults(even_sam_sig_30, data_name = 'even sam sig_30 enet', regularize = T)


##### gen sam

#### 10
even_gen_sam_sig_10 <- runModels(quan_cases_sam_gen, 
                                 model = 'enet',
                                 rand = F,
                                 bump_hunter = T,
                                 bump_hunter_data = quan_even_gen_sam_sig_10)

even_gen_sam_10_sig_table <- extractResults(even_gen_sam_sig_10, data_name = 'even gen sam sig_10 enet', regularize = T)


#### 20
even_gen_sam_sig_20 <- runModels(quan_cases_sam_gen, 
                                 model = 'enet',
                                 rand = F,
                                 bump_hunter = T,
                                 bump_hunter_data = quan_even_gen_sam_sig_20)

even_gen_sam_sig_20_table <- extractResults(even_gen_sam_sig_20, data_name = 'even gen sam sig_20 enet', regularize = T)

# #### 30
# even_gen_sam_sig_30 <- runModels(quan_cases_sam_gen, 
#                              model = 'enet',
#                              rand = F,
#                              bump_hunter = T,
#                              bump_hunter_data = quan_even_gen_sam_sig_30)
# 
# even_gen_sam_sig_30_table <- extractResults(even_gen_sam_sig_30, data_name = 'even gen sam sig_30 enet', regularize = T)



##### gen sen

#### 10
even_gen_sen_sig_10 <- runModels(quan_cases_sen_gen, 
                                 model = 'enet',
                                 rand = F,
                                 bump_hunter = T,
                                 bump_hunter_data = quan_even_gen_sen_sig_10)

even_gen_sen_sig_10_table <- extractResults(even_gen_sen_sig_10, data_name = 'even gen sen sig_10 enet', regularize = T)


# #### 20
# even_gen_sen_sig_20 <- runModels(quan_cases_sen_gen, 
#                              model = 'enet',
#                              rand = F,
#                              bump_hunter = T,
#                              bump_hunter_data = quan_even_gen_sen_sig_20)
# 
# even_gen_sen_sig_20_table <- extractResults(even_gen_sen_sig_20, data_name = 'even gen sen sig_20 enet', regularize = T)

# #### 30
# even_gen_sen_sig_30 <- runModels(quan_cases_sen_gen, 
#                              model = 'enet',
#                              rand = F,
#                              bump_hunter = T,
#                              bump_hunter_data = quan_even_gen_sen_sig_30)
# 
# even_gen_sen_sig_30_table <- extractResults(even_gen_sen_sig_30, data_name = 'even gen sen sig_30 enet', regularize = T)

###########
# get table for even not sig
###########
even_sig_enet <- rbind(even_sig_10_table, even_sig_20_table,
                     even_gen_sig_10_table,
                     even_sen_sig_10_table, even_sen_sig_20_table,
                     even_sam_sig_10_table, even_sam_sig_10_table,
                     even_gen_sen_sig_10_table)

rm(even_sig_10_table, even_sig_20_table,
   even_gen_sig_10_table,
   even_sen_sig_10_table, even_sen_sig_20_table,
   even_sam_sig_10_table,
   even_gen_sen_sig_10_table,even_gen_sam_sig_20_table)


##################################
#######fwer

#####not gen

#### 10
even_fwer_10 <- runModels(quan_cases, 
                          model = 'enet',
                          rand = F,
                          bump_hunter = T,
                          bump_hunter_data = quan_even_fwer_10)

even_fwer_10_table <- extractResults(even_fwer_10, data_name = 'even fwer 10 enet', regularize = T)

#### 20
even_fwer_20 <- runModels(quan_cases, 
                          model = 'enet',
                          rand = F,
                          bump_hunter = T,
                          bump_hunter_data = quan_even_fwer_20)

even_fwer_20_table <- extractResults(even_fwer_20, data_name = 'even fwer_20 enet', regularize = T)


#### 30
even_fwer_30 <- runModels(quan_cases, 
                          model = 'enet',
                          rand = F,
                          bump_hunter = T,
                          bump_hunter_data = quan_even_fwer_30)

even_fwer_30_table <- extractResults(even_fwer_30, data_name = 'even fwer_30 enet', regularize = T)


##### gen

#### 10
even_gen_fwer_10 <- runModels(quan_cases_gen, 
                              model = 'enet',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_gen_fwer_10)

even_gen_fwer_10_table <- extractResults(even_gen_fwer_10, data_name = 'even gen fwer_10 enet', regularize = T)


#### 20
even_gen_fwer_20 <- runModels(quan_cases_gen, 
                              model = 'enet',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_gen_fwer_20)

even_gen_fwer_20_table <- extractResults(even_gen_fwer_20, data_name = 'even gen fwer_20 enet', regularize = T)

#### 30
even_gen_fwer_30 <- runModels(quan_cases_gen,
                              model = 'enet',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_gen_fwer_30)

even_gen_fwer_30_table <- extractResults(even_gen_fwer_30, data_name = 'even gen fwer_30 enet', regularize = T)

##### sen

#### 10
even_sen_fwer_10 <- runModels(quan_cases_sen, 
                              model = 'enet',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_fwer_10)

even_sen_fwer_10_table <- extractResults(even_sen_fwer_10, data_name = 'even sen fwer_10 enet', regularize = T)


#### 20
even_sen_fwer_20 <- runModels(quan_cases_sen, 
                              model = 'enet',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_fwer_20)

even_sen_fwer_20_table <- extractResults(even_sen_fwer_20, data_name = 'even sen fwer_20 enet', regularize = T)

#### 30
even_sen_fwer_30 <- runModels(quan_cases_sen,
                              model = 'enet',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_fwer_30)

even_sen_fwer_30_table <- extractResults(even_sen_fwer_30, data_name = 'even sen fwer_30 enet', regularize = T)


##### sam

#### 10
even_sam_fwer_10 <- runModels(quan_cases_sam, 
                              model = 'enet',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_fwer_10)

even_sam_fwer_10_table <- extractResults(even_sam_fwer_10, data_name = 'even sam fwer_10 enet', regularize = T)


#### 20
even_sam_fwer_20 <- runModels(quan_cases_sam, 
                              model = 'enet',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_fwer_20)

even_sam_fwer_20_table <- extractResults(even_sam_fwer_20, data_name = 'even sam fwer_20 enet', regularize = T)

#### 30
even_sam_fwer_30 <- runModels(quan_cases_sam,
                              model = 'enet',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_fwer_30)

even_sam_fwer_30_table <- extractResults(even_sam_fwer_30, data_name = 'even sam fwer_30 enet', regularize = T)


##### gen sam

#### 10
even_gen_sam_fwer_10 <- runModels(quan_cases_sam_gen, 
                                  model = 'enet',
                                  rand = F,
                                  bump_hunter = T,
                                  bump_hunter_data = quan_even_gen_sam_fwer_10)

even_gen_sam_fwer_10_table <- extractResults(even_gen_sam_fwer_10, data_name = 'even gen sam fwer_10 enet', regularize = T)


#### 20
even_gen_sam_fwer_20 <- runModels(quan_cases_sam_gen, 
                                  model = 'enet',
                                  rand = F,
                                  bump_hunter = T,
                                  bump_hunter_data = quan_even_gen_sam_fwer_20)

even_gen_sam_fwer_20_table <- extractResults(even_gen_sam_fwer_20, data_name = 'even gen sam fwer_20 enet', regularize = T)

#### 30
even_gen_sam_fwer_30 <- runModels(quan_cases_sam_gen,
                                  model = 'enet',
                                  rand = F,
                                  bump_hunter = T,
                                  bump_hunter_data = quan_even_gen_sam_fwer_30)

even_gen_sam_fwer_30_table <- extractResults(even_gen_sam_fwer_30, data_name = 'even gen sam fwer_30 enet', regularize = T)



##### gen sen

#### 10
even_gen_sen_fwer_10 <- runModels(quan_cases_sen_gen, 
                                  model = 'enet',
                                  rand = F,
                                  bump_hunter = T,
                                  bump_hunter_data = quan_even_gen_sen_fwer_10)

even_gen_sen_fwer_10_table <- extractResults(even_gen_sen_fwer_10, data_name = 'even gen sen fwer_10 enet', regularize = T)


#### 20
even_gen_sen_fwer_20 <- runModels(quan_cases_sen_gen,
                                  model = 'enet',
                                  rand = F,
                                  bump_hunter = T,
                                  bump_hunter_data = quan_even_gen_sen_fwer_20)

even_gen_sen_fwer_20_table <- extractResults(even_gen_sen_fwer_20, data_name = 'even gen sen fwer_20 enet', regularize = T)

#### 30
even_gen_sen_fwer_30 <- runModels(quan_cases_sen_gen,
                                  model = 'enet',
                                  rand = F,
                                  bump_hunter = T,
                                  bump_hunter_data = quan_even_gen_sen_fwer_30)

even_gen_sen_fwer_30_table <- extractResults(even_gen_sen_fwer_30, data_name = 'even gen sen fwer_30 enet', regularize = T)

###########
# get table for even not fwer
###########


even_fwer_enet <- rbind(even_fwer_10_table, even_fwer_20_table, even_fwer_30_table,
                      even_gen_fwer_10_table, even_gen_fwer_20_table, even_gen_fwer_30_table,
                      even_sen_fwer_10_table, even_sen_fwer_20_table,even_sen_fwer_30_table,
                      even_sam_fwer_10_table, even_sam_fwer_20_table, even_sam_fwer_30_table,
                      even_gen_sen_fwer_10_table, even_gen_sen_fwer_20_table, even_gen_sen_fwer_30_table,
                      even_gen_sam_fwer_10_table, even_gen_sam_fwer_20_table, even_gen_sam_fwer_30_table)

rm(even_fwer_10_table, even_fwer_20_table, even_fwer_30_table,
   even_gen_fwer_10_table, even_gen_fwer_20_table, even_gen_fwer_30_table,
   even_sen_fwer_10_table, even_sen_fwer_20_table,even_sen_fwer_30_table,
   even_sam_fwer_10_table, even_sam_fwer_20_table, even_sam_fwer_30_table,
   even_gen_sen_fwer_10_table, even_gen_sen_fwer_20_table, even_gen_sen_fwer_30_table,
   even_gen_sam_fwer_10_table, even_gen_sam_fwer_20_table, even_gen_sam_fwer_30_table)



########################################################################3
##########
# run models
##########

##################################
####### not significant

#####not gen

#### 10
even_10 <- runModels(quan_cases, 
                     model = 'lasso',
                     rand = F,
                     bump_hunter = T,
                     bump_hunter_data = quan_even_10)

even_10_table <- extractResults(even_10, data_name = 'even_10 lasso', regularize = T)

#### 20
even_20 <- runModels(quan_cases, 
                     model = 'lasso',
                     rand = F,
                     bump_hunter = T,
                     bump_hunter_data = quan_even_20)

even_20_table <- extractResults(even_10, data_name = 'even_20 lasso', regularize = T)


#### 30
even_30 <- runModels(quan_cases, 
                     model = 'lasso',
                     rand = F,
                     bump_hunter = T,
                     bump_hunter_data = quan_even_30)

even_30_table <- extractResults(even_30, data_name = 'even_30 lasso', regularize = T)


##### gen

#### 10
even_gen_10 <- runModels(quan_cases_gen, 
                         model = 'lasso',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_gen_10)

even_gen_10_table <- extractResults(even_gen_10, data_name = 'even gen_10 lasso', regularize = T)


#### 20
even_gen_20 <- runModels(quan_cases_gen, 
                         model = 'lasso',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_gen_20)

even_gen_20_table <- extractResults(even_gen_20, data_name = 'even gen_20 lasso', regularize = T)

#### 30
even_gen_30 <- runModels(quan_cases_gen, 
                         model = 'lasso',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_gen_30)

even_gen_30_table <- extractResults(even_gen_30, data_name = 'even gen_30 lasso', regularize = T)

##### sen

#### 10
even_sen_10 <- runModels(quan_cases_sen, 
                         model = 'lasso',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_10)

even_sen_10_table <- extractResults(even_sen_10, data_name = 'even sen_10 lasso', regularize = T)


#### 20
even_sen_20 <- runModels(quan_cases_sen, 
                         model = 'lasso',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_20)

even_sen_20_table <- extractResults(even_sen_20, data_name = 'even sen_20 lasso', regularize = T)

#### 30
even_sen_30 <- runModels(quan_cases_sen, 
                         model = 'lasso',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_30)

even_sen_30_table <- extractResults(even_sen_30, data_name = 'even sen_30 lasso', regularize = T)


##### sam

#### 10
even_sam_10 <- runModels(quan_cases_sam, 
                         model = 'lasso',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_10)

even_sam_10_table <- extractResults(even_sam_10, data_name = 'even sam_10 lasso', regularize = T)


#### 20
even_sam_20 <- runModels(quan_cases_sam, 
                         model = 'lasso',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_20)

even_sam_20_table <- extractResults(even_sam_20, data_name = 'even sam_20 lasso', regularize = T)

#### 30
even_sam_30 <- runModels(quan_cases_sam, 
                         model = 'lasso',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_30)

even_sam_30_table <- extractResults(even_sam_30, data_name = 'even sam_30 lasso', regularize = T)


##### gen sam

#### 10
even_gen_sam_10 <- runModels(quan_cases_sam_gen, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sam_10)

even_gen_sam_10_table <- extractResults(even_gen_sam_10, data_name = 'even gen sam_10 lasso', regularize = T)


#### 20
even_gen_sam_20 <- runModels(quan_cases_sam_gen, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sam_20)

even_gen_sam_20_table <- extractResults(even_gen_sam_20, data_name = 'even gen sam_20 lasso', regularize = T)

#### 30
even_gen_sam_30 <- runModels(quan_cases_sam_gen, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sam_30)

even_gen_sam_30_table <- extractResults(even_gen_sam_30, data_name = 'even gen sam_30 lasso', regularize = T)



##### gen sen

#### 10
even_gen_sen_10 <- runModels(quan_cases_sen_gen, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sen_10)

even_gen_sen_10_table <- extractResults(even_gen_sen_10, data_name = 'even gen sen_10 lasso', regularize = T)


#### 20
even_gen_sen_20 <- runModels(quan_cases_sen_gen, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sen_20)

even_gen_sen_20_table <- extractResults(even_gen_sen_20, data_name = 'even gen sen_20 lasso', regularize = T)

#### 30
even_gen_sen_30 <- runModels(quan_cases_sen_gen, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sen_30)

even_gen_sen_30_table <- extractResults(even_gen_sen_30, data_name = 'even gen sen_30 lasso', regularize = T)

###########
# get table for even not sig
###########
even_not_sig_lasso <- rbind(even_10_table, even_20_table, even_30_table,
                         even_gen_10_table, even_gen_20_table, even_gen_30_table,
                         even_sen_10_table, even_sen_20_table, even_sen_30_table,
                         even_sam_10_table, even_sam_20_table, even_sam_30_table,
                         even_gen_sen_10_table, even_gen_sen_20_table, even_gen_sen_30_table,
                         even_gen_sam_10_table, even_gen_sam_20_table, even_gen_sam_30_table)

rm(even_10_table, even_20_table, even_30_table,
   even_gen_10_table, even_gen_20_table, even_gen_30_table,
   even_sen_10_table, even_sen_20_table, even_sen_30_table,
   even_sam_10_table, even_sam_20_table, even_sam_30_table,
   even_gen_sen_10_table, even_gen_sen_20_table, even_gen_sen_30_table,
   even_gen_sam_10_table, even_gen_sam_20_table, even_gen_sam_30_table)


##################################
#######significant

#####not gen

#### 10
even_sig_10 <- runModels(quan_cases, 
                         model = 'lasso',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_sig_10)

even_sig_10_table <- extractResults(even_sig_10, data_name = 'even sig_10 lasso', regularize = T)

#### 20
even_sig_20 <- runModels(quan_cases, 
                         model = 'lasso',
                         rand = F,
                         bump_hunter = T,
                         bump_hunter_data = quan_even_sig_20)

even_sig_20_table <- extractResults(even_sig_20, data_name = 'even sig_20 lasso', regularize = T)


# #### 30
# even_sig_30 <- runModels(quan_cases, 
#                          model = 'lasso',
#                          rand = F,
#                          bump_hunter = T,
#                          bump_hunter_data = quan_even_sig_30)
# 
# even_sig_30_table <- extractResults(even_sig_30, data_name = 'even sig_30 lasso', regularize = T)
# 

##### gen

#### 10
even_gen_sig_10 <- runModels(quan_cases_gen, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sig_10)

even_gen_sig_10_table <- extractResults(even_gen_sig_10, data_name = 'even gen sig_10 lasso', regularize = T)


#### 20
even_gen_sig_20 <- runModels(quan_cases_gen, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_gen_sig_20)

even_gen_20_sig_table <- extractResults(even_gen_sig_20, data_name = 'even gen sig_20 lasso', regularize = T)

# #### 30
# even_gen_sig_30 <- runModels(quan_cases_gen, 
#                          model = 'lasso',
#                          rand = F,
#                          bump_hunter = T,
#                          bump_hunter_data = quan_even_gen_sig_30)
# 
# even_gen_sig_30_table <- extractResults(even_gen_sig_30, data_name = 'even gen sig_30 lasso', regularize = T)

##### sen

#### 10
even_sen_sig_10 <- runModels(quan_cases_sen, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_sig_10)

even_sen_sig_10_table <- extractResults(even_sen_sig_10, data_name = 'even sen sig_10 lasso', regularize = T)


#### 20
even_sen_sig_20 <- runModels(quan_cases_sen, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_sig_20)

even_sen_sig_20_table <- extractResults(even_sen_sig_20, data_name = 'even sen sig_20 lasso', regularize = T)

# #### 30
# even_sen_sig_30 <- runModels(quan_cases_sen, 
#                          model = 'lasso',
#                          rand = F,
#                          bump_hunter = T,
#                          bump_hunter_data = quan_even_sig_30)
# 
# even_sen_sig_30_table <- extractResults(even_sen_sig_30, data_name = 'even sen sig_30 lasso', regularize = T)


##### sam

#### 10
even_sam_sig_10 <- runModels(quan_cases_sam, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_sig_10)

even_sam_sig_10_table <- extractResults(even_sam_sig_10, data_name = 'even sam sig_10 lasso', regularize = T)


#### 20
even_sam_sig_20 <- runModels(quan_cases_sam, 
                             model = 'lasso',
                             rand = F,
                             bump_hunter = T,
                             bump_hunter_data = quan_even_sig_20)

even_sam_sig_20_table <- extractResults(even_sam_sig_20, data_name = 'even sam sig_20 lasso', regularize = T)

# #### 30
# even_sam_sig_30 <- runModels(quan_cases_sam, 
#                          model = 'lasso',
#                          rand = F,
#                          bump_hunter = T,
#                          bump_hunter_data = quan_even_sig_30)
# 
# even_sam_sig_30_table <- extractResults(even_sam_sig_30, data_name = 'even sam sig_30 lasso', regularize = T)


##### gen sam

#### 10
even_gen_sam_sig_10 <- runModels(quan_cases_sam_gen, 
                                 model = 'lasso',
                                 rand = F,
                                 bump_hunter = T,
                                 bump_hunter_data = quan_even_gen_sam_sig_10)

even_gen_sam_10_sig_table <- extractResults(even_gen_sam_sig_10, data_name = 'even gen sam sig_10 lasso', regularize = T)


#### 20
even_gen_sam_sig_20 <- runModels(quan_cases_sam_gen, 
                                 model = 'lasso',
                                 rand = F,
                                 bump_hunter = T,
                                 bump_hunter_data = quan_even_gen_sam_sig_20)

even_gen_sam_sig_20_table <- extractResults(even_gen_sam_sig_20, data_name = 'even gen sam sig_20 lasso', regularize = T)

# #### 30
# even_gen_sam_sig_30 <- runModels(quan_cases_sam_gen, 
#                              model = 'lasso',
#                              rand = F,
#                              bump_hunter = T,
#                              bump_hunter_data = quan_even_gen_sam_sig_30)
# 
# even_gen_sam_sig_30_table <- extractResults(even_gen_sam_sig_30, data_name = 'even gen sam sig_30 lasso', regularize = T)



##### gen sen

#### 10
even_gen_sen_sig_10 <- runModels(quan_cases_sen_gen, 
                                 model = 'lasso',
                                 rand = F,
                                 bump_hunter = T,
                                 bump_hunter_data = quan_even_gen_sen_sig_10)

even_gen_sen_sig_10_table <- extractResults(even_gen_sen_sig_10, data_name = 'even gen sen sig_10 lasso', regularize = T)


# #### 20
# even_gen_sen_sig_20 <- runModels(quan_cases_sen_gen, 
#                              model = 'lasso',
#                              rand = F,
#                              bump_hunter = T,
#                              bump_hunter_data = quan_even_gen_sen_sig_20)
# 
# even_gen_sen_sig_20_table <- extractResults(even_gen_sen_sig_20, data_name = 'even gen sen sig_20 lasso', regularize = T)

# #### 30
# even_gen_sen_sig_30 <- runModels(quan_cases_sen_gen, 
#                              model = 'lasso',
#                              rand = F,
#                              bump_hunter = T,
#                              bump_hunter_data = quan_even_gen_sen_sig_30)
# 
# even_gen_sen_sig_30_table <- extractResults(even_gen_sen_sig_30, data_name = 'even gen sen sig_30 lasso', regularize = T)

###########
# get table for even not sig
###########
even_sig_lasso <- rbind(even_sig_10_table, even_sig_20_table,
                     even_gen_sig_10_table,
                     even_sen_sig_10_table, even_sen_sig_20_table,
                     even_sam_sig_10_table, 
                     even_gen_sen_sig_10_table,even_gen_sam_sig_20_table)

rm(even_sig_10_table, even_sig_20_table,
   even_gen_sig_10_table,
   even_sen_sig_10_table, even_sen_sig_20_table,
   even_sam_sig_10_table, 
   even_gen_sen_sig_10_table,even_gen_sam_sig_20_table)


##################################
#######fwer

#####not gen

#### 10
even_fwer_10 <- runModels(quan_cases, 
                          model = 'lasso',
                          rand = F,
                          bump_hunter = T,
                          bump_hunter_data = quan_even_fwer_10)

even_fwer_10_table <- extractResults(even_fwer_10, data_name = 'even fwer 10 lasso', regularize = T)

#### 20
even_fwer_20 <- runModels(quan_cases, 
                          model = 'lasso',
                          rand = F,
                          bump_hunter = T,
                          bump_hunter_data = quan_even_fwer_20)

even_fwer_20_table <- extractResults(even_fwer_20, data_name = 'even fwer_20 lasso', regularize = T)


#### 30
even_fwer_30 <- runModels(quan_cases, 
                          model = 'lasso',
                          rand = F,
                          bump_hunter = T,
                          bump_hunter_data = quan_even_fwer_30)

even_fwer_30_table <- extractResults(even_fwer_30, data_name = 'even fwer_30 lasso', regularize = T)


##### gen

#### 10
even_gen_fwer_10 <- runModels(quan_cases_gen, 
                              model = 'lasso',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_gen_fwer_10)

even_gen_fwer_10_table <- extractResults(even_gen_fwer_10, data_name = 'even gen fwer_10 lasso', regularize = T)


#### 20
even_gen_fwer_20 <- runModels(quan_cases_gen, 
                              model = 'lasso',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_gen_fwer_20)

even_gen_fwer_20_table <- extractResults(even_gen_fwer_20, data_name = 'even gen fwer_20 lasso', regularize = T)

#### 30
even_gen_fwer_30 <- runModels(quan_cases_gen,
                              model = 'lasso',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_gen_fwer_30)

even_gen_fwer_30_table <- extractResults(even_gen_fwer_30, data_name = 'even gen fwer_30 lasso', regularize = T)

##### sen

#### 10
even_sen_fwer_10 <- runModels(quan_cases_sen, 
                              model = 'lasso',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_fwer_10)

even_sen_fwer_10_table <- extractResults(even_sen_fwer_10, data_name = 'even sen fwer_10 lasso', regularize = T)


#### 20
even_sen_fwer_20 <- runModels(quan_cases_sen, 
                              model = 'lasso',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_fwer_20)

even_sen_fwer_20_table <- extractResults(even_sen_fwer_20, data_name = 'even sen fwer_20 lasso', regularize = T)

#### 30
even_sen_fwer_30 <- runModels(quan_cases_sen,
                              model = 'lasso',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_fwer_30)

even_sen_fwer_30_table <- extractResults(even_sen_fwer_30, data_name = 'even sen fwer_30 lasso', regularize = T)


##### sam

#### 10
even_sam_fwer_10 <- runModels(quan_cases_sam, 
                              model = 'lasso',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_fwer_10)

even_sam_fwer_10_table <- extractResults(even_sam_fwer_10, data_name = 'even sam fwer_10 lasso', regularize = T)


#### 20
even_sam_fwer_20 <- runModels(quan_cases_sam, 
                              model = 'lasso',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_fwer_20)

even_sam_fwer_20_table <- extractResults(even_sam_fwer_20, data_name = 'even sam fwer_20 lasso', regularize = T)

#### 30
even_sam_fwer_30 <- runModels(quan_cases_sam,
                              model = 'lasso',
                              rand = F,
                              bump_hunter = T,
                              bump_hunter_data = quan_even_fwer_30)

even_sam_fwer_30_table <- extractResults(even_sam_fwer_30, data_name = 'even sam fwer_30 lasso', regularize = T)


##### gen sam

#### 10
even_gen_sam_fwer_10 <- runModels(quan_cases_sam_gen, 
                                  model = 'lasso',
                                  rand = F,
                                  bump_hunter = T,
                                  bump_hunter_data = quan_even_gen_sam_fwer_10)

even_gen_sam_fwer_10_table <- extractResults(even_gen_sam_fwer_10, data_name = 'even gen sam fwer_10 lasso', regularize = T)


#### 20
even_gen_sam_fwer_20 <- runModels(quan_cases_sam_gen, 
                                  model = 'lasso',
                                  rand = F,
                                  bump_hunter = T,
                                  bump_hunter_data = quan_even_gen_sam_fwer_20)

even_gen_sam_fwer_20_table <- extractResults(even_gen_sam_fwer_20, data_name = 'even gen sam fwer_20 lasso', regularize = T)

#### 30
even_gen_sam_fwer_30 <- runModels(quan_cases_sam_gen,
                                  model = 'lasso',
                                  rand = F,
                                  bump_hunter = T,
                                  bump_hunter_data = quan_even_gen_sam_fwer_30)

even_gen_sam_fwer_30_table <- extractResults(even_gen_sam_fwer_30, data_name = 'even gen sam fwer_30 lasso', regularize = T)



##### gen sen

#### 10
even_gen_sen_fwer_10 <- runModels(quan_cases_sen_gen, 
                                  model = 'lasso',
                                  rand = F,
                                  bump_hunter = T,
                                  bump_hunter_data = quan_even_gen_sen_fwer_10)

even_gen_sen_fwer_10_table <- extractResults(even_gen_sen_fwer_10, data_name = 'even gen sen fwer_10 lasso', regularize = T)


#### 20
even_gen_sen_fwer_20 <- runModels(quan_cases_sen_gen,
                                  model = 'lasso',
                                  rand = F,
                                  bump_hunter = T,
                                  bump_hunter_data = quan_even_gen_sen_fwer_20)

even_gen_sen_fwer_20_table <- extractResults(even_gen_sen_fwer_20, data_name = 'even gen sen fwer_20 lasso', regularize = T)

#### 30
even_gen_sen_fwer_30 <- runModels(quan_cases_sen_gen,
                                  model = 'lasso',
                                  rand = F,
                                  bump_hunter = T,
                                  bump_hunter_data = quan_even_gen_sen_fwer_30)

even_gen_sen_fwer_30_table <- extractResults(even_gen_sen_fwer_30, data_name = 'even gen sen fwer_30 lasso', regularize = T)

###########
# get table for even not fwer
###########
even_fwer_lasso <- rbind(even_fwer_10_table, even_fwer_20_table, even_fwer_30_table,
                      even_gen_fwer_10_table, even_gen_fwer_20_table, even_gen_fwer_30_table,
                      even_sen_fwer_10_table, even_sen_fwer_20_table,even_sen_fwer_30_table,
                      even_sam_fwer_10_table, even_sam_fwer_20_table, even_sam_fwer_30_table,
                      even_gen_sen_fwer_10_table, even_gen_sen_fwer_20_table, even_gen_sen_fwer_30_table,
                      even_gen_sam_fwer_10_table, even_gen_sam_fwer_20_table, even_gen_sam_fwer_30_table)

rm(even_fwer_10_table, even_fwer_20_table, even_fwer_30_table,
   even_gen_fwer_10_table, even_gen_fwer_20_table, even_gen_fwer_30_table,
   even_sen_fwer_10_table, even_sen_fwer_20_table,even_sen_fwer_30_table,
   even_sam_fwer_10_table, even_sam_fwer_20_table, even_sam_fwer_30_table,
   even_gen_sen_fwer_10_table, even_gen_sen_fwer_20_table, even_gen_sen_fwer_30_table,
   even_gen_sam_fwer_10_table, even_gen_sam_fwer_20_table, even_gen_sam_fwer_30_table)


#########################################################################################
# save tables

# combine

result_table <- rbind(even_not_sig_rf,
                      even_sig_rf,
                      even_fwer_rf,
                      even_not_sig_enet,
                      even_sig_enet,
                      even_fwer_enet,
                      even_not_sig_lasso,
                      even_sig_lasso,
                      even_fwer_lasso)

# save 
saveRDS(result_table, paste0(model_data, '/result_table.rds'))
