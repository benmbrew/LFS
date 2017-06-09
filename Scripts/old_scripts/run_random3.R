###################################################################3
# run random for both cases and  controls

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
quan_cases_sam <- readRDS(paste0(model_data, '/quan_cases_sam.rda'))

quan_cases_sam_gen <- readRDS(paste0(model_data, '/quan_cases_sam_gen.rda'))

# read controls
quan_controls <- readRDS(paste0(model_data, '/quan_controls.rda'))

quan_controls_gen <- readRDS(paste0(model_data, '/quan_controls_gen.rda'))


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
                      seed_num,
                      control = control) 
{
  
  # get differenct variations of data
  if (random) {
    
    data <- subsetDat(data, 
                      random = T, 
                      num_feat,
                      seed_num = seed_num)
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
    
    data_result <- rfPredictReg(data, cutoff = .7, iterations = 5, control = control)
    # data_resid_result <- rfPredictReg(data_resid, cutoff = .7, iterations = 5)
  }
  
  if (model == 'enet') {
    
    data_result <- enetPredReg(data, N_CV_REPEATS = 2, nfolds = 5,cutoff = .7, iterations = 5, control = control)
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
# data <- quan_cases_sam
# rand_num <- 5
# num_feat <- 10
# data_name <- 'raw_rand_10'
# exlcude <- T
# union_features <-  union_features
getRand <- function(data, model,
                    control,
                    rand_num, 
                    num_feat, 
                    exclude, 
                    union_features, 
                    data_name) 
{
  
  if (exclude) {
    
    features <- colnames(data[, 10:ncol(data)])
    # keep_these <- features[!features %in% union_features]
    data <- data[, c('ids', 'type', 'gender', 'sentrix_id', 'batch', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                     'age_sample_collection', features)]
  }
  
  
  if (model == 'rf') {
    
    rf_holder = list()
    rf_table = list()
    
    for (rand_number in 1:rand_num) {
      
      rf_holder[[rand_number]] <- runModels(data, 
                                         model = 'rf',
                                         bump_hunter = F,
                                         num_feat = num_feat,
                                         random = T,
                                         seed_num = rand_number,
                                         control = control)
      
      rf_table[[rand_number]] <- extractResults(rf_holder[[rand_number]],
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
                                           seed_num = rand_num,
                                           control = control)
      
      enet_table[[rand_num]] <- extractResults(enet_holder[[rand_num]],
                                               data_name = data_name,
                                               regularize = T)
      
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

# # get list of features
# df_names <- ls()[sapply(mget(ls(), .GlobalEnv), is.data.frame)]
# dfs <-  sapply( df_names, function(x)  get( x )  )
# dfs[grepl('cases|controls', names(dfs))] <- NULL
#
# # combine and remove duplicates
# feats <- do.call(rbind, dfs)
#
# # union
# union_features <- as.character(feats[!duplicated(feats),])


###CASES

##########
# rf
##########

# cases_10
cases_10 <- getRand(quan_cases_sam_gen, 
                    model = 'rf',
                    control = F,
                    rand_num = 10, 
                    num_feat = 10, 
                    exclude = T,
                    union_features = union_features,
                    data_name = "cases_10")


# cases_100
cases_100 <- getRand(quan_cases_sam_gen, 
                    model = 'rf',
                    control = F,
                    rand_num = 10, 
                    num_feat = 100, 
                    exclude = T,
                    union_features = union_features,
                    data_name = "cases_100")

# cases_500
cases_500 <- getRand(quan_cases_sam_gen, 
                     model = 'rf',
                     control = F,
                     rand_num = 10, 
                     num_feat = 500, 
                     exclude = T,
                     union_features = union_features,
                     data_name = "cases_500")

# cases_1000
cases_1000 <- getRand(quan_cases_sam_gen, 
                     model = 'rf',
                     control = F,
                     rand_num = 10, 
                     num_feat = 1000, 
                     exclude = T,
                     union_features = union_features,
                     data_name = "cases_1000")

# cases_5000
cases_5000 <- getRand(quan_cases_sam_gen, 
                      model = 'rf',
                      control = F,
                      rand_num = 10, 
                      num_feat = 5000, 
                      exclude = T,
                      union_features = union_features,
                      data_name = "cases_5000")

# cases_10000
cases_10000 <- getRand(quan_cases_sam_gen, 
                      model = 'rf',
                      control = F,
                      rand_num = 10, 
                      num_feat = 10000, 
                      exclude = T,
                      union_features = union_features,
                      data_name = "cases_10000")


# cases_50000
cases_50000 <- getRand(quan_cases_sam_gen, 
                      model = 'rf',
                      control = F,
                      rand_num = 10, 
                      num_feat = 50000, 
                      exclude = T,
                      union_features = union_features,
                      data_name = "cases_50000")

# cases_1000
cases_100000 <- getRand(quan_cases_sam_gen, 
                      model = 'rf',
                      control = F,
                      rand_num = 10, 
                      num_feat = 100000, 
                      exclude = T,
                      union_features = union_features,
                      data_name = "cases_100000")



##########
# enet
##########

# cases_10
cases_10 <- getRand(quan_cases_sam_gen, 
                    model = 'enet',
                    control = F,
                    rand_num = 10, 
                    num_feat = 10, 
                    exclude = T,
                    union_features = union_features,
                    data_name = "cases_10")


# cases_100
cases_100 <- getRand(quan_cases_sam_gen, 
                     model = 'enet',
                     control = F,
                     rand_num = 10, 
                     num_feat = 100, 
                     exclude = T,
                     union_features = union_features,
                     data_name = "cases_100")

# cases_500
cases_500 <- getRand(quan_cases_sam_gen, 
                     model = 'enet',
                     control = F,
                     rand_num = 10, 
                     num_feat = 500, 
                     exclude = T,
                     union_features = union_features,
                     data_name = "cases_500")

# cases_1000
cases_1000 <- getRand(quan_cases_sam_gen, 
                      model = 'enet',
                      control = F,
                      rand_num = 10, 
                      num_feat = 1000, 
                      exclude = T,
                      union_features = union_features,
                      data_name = "cases_1000")

# cases_5000
cases_5000 <- getRand(quan_cases_sam_gen, 
                      model = 'enet',
                      control = F,
                      rand_num = 10, 
                      num_feat = 5000, 
                      exclude = T,
                      union_features = union_features,
                      data_name = "cases_5000")

# cases_10000
cases_10000 <- getRand(quan_cases_sam_gen, 
                       model = 'enet',
                       control = F,
                       rand_num = 10, 
                       num_feat = 10000, 
                       exclude = T,
                       union_features = union_features,
                       data_name = "cases_10000")


# cases_50000
cases_50000 <- getRand(quan_cases_sam_gen, 
                       model = 'enet',
                       control = F,
                       rand_num = 10, 
                       num_feat = 50000, 
                       exclude = T,
                       union_features = union_features,
                       data_name = "cases_50000")

# cases_1000
cases_100000 <- getRand(quan_cases_sam_gen, 
                        model = 'enet',
                        control = F,
                        rand_num = 10, 
                        num_feat = 100000, 
                        exclude = T,
                        union_features = union_features,
                        data_name = "cases_100000")



############################################################################################3
### Controls

##########
# rf
##########
# ungen 10 features
con_10 <- getRand(quan_controls_gen, 
                  model = 'rf',
                  control = T,
                  rand_num = 10, 
                  num_feat = 10, 
                  exclude = T,
                  union_features = union_features,
                  data_name = "con_10")

# ungen 100 features
con_100 <- getRand(quan_controls_gen, 
                  model = 'rf',
                  control = T,
                  rand_num = 10, 
                  num_feat = 100, 
                  exclude = T,
                  union_features = union_features,
                  data_name = "con_100")

# ungen 500 features
con_500 <- getRand(quan_controls_gen, 
                   model = 'rf',
                   control = T,
                   rand_num = 10, 
                   num_feat = 500, 
                   exclude = T,
                   union_features = union_features,
                   data_name = "con_500")

# ungen 1000 features
con_1000 <- getRand(quan_controls_gen, 
                   model = 'rf',
                   control = T,
                   rand_num = 10, 
                   num_feat = 1000, 
                   exclude = T,
                   union_features = union_features,
                   data_name = "con_1000")

# ungen 10000 features
con_10000 <- getRand(quan_controls_gen, 
                   model = 'rf',
                   control = T,
                   rand_num = 10, 
                   num_feat = 10000, 
                   exclude = T,
                   union_features = union_features,
                   data_name = "con_10000")

# ungen 100 features
con_50000 <- getRand(quan_controls_gen, 
                   model = 'rf',
                   control = T,
                   rand_num = 10, 
                   num_feat = 50000, 
                   exclude = T,
                   union_features = union_features,
                   data_name = "con_50000")

# ungen 100 features
con_100000 <- getRand(quan_controls_gen, 
                     model = 'rf',
                     control = T,
                     rand_num = 10, 
                     num_feat = 100000, 
                     exclude = T,
                     union_features = union_features,
                     data_name = "con_100000")


##########
#enet
##########

### Controls
# ungen 10 features
con_10 <- getRand(quan_controls_gen, 
                  model = 'enet',
                  control = T,
                  rand_num = 10, 
                  num_feat = 10, 
                  exclude = T,
                  union_features = union_features,
                  data_name = "con_10")

# ungen 100 features
con_100 <- getRand(quan_controls_gen, 
                   model = 'enet',
                   control = T,
                   rand_num = 10, 
                   num_feat = 100, 
                   exclude = T,
                   union_features = union_features,
                   data_name = "con_100")

# ungen 500 features
con_500 <- getRand(quan_controls_gen, 
                   model = 'enet',
                   control = T,
                   rand_num = 10, 
                   num_feat = 500, 
                   exclude = T,
                   union_features = union_features,
                   data_name = "con_500")

# ungen 1000 features
con_1000 <- getRand(quan_controls_gen, 
                    model = 'enet',
                    control = T,
                    rand_num = 10, 
                    num_feat = 1000, 
                    exclude = T,
                    union_features = union_features,
                    data_name = "con_1000")

# ungen 10000 features
con_10000 <- getRand(quan_controls_gen, 
                     model = 'enet',
                     control = T,
                     rand_num = 10, 
                     num_feat = 10000, 
                     exclude = T,
                     union_features = union_features,
                     data_name = "con_10000")

# ungen 100 features
con_50000 <- getRand(quan_controls_gen, 
                     model = 'enet',
                     control = T,
                     rand_num = 10, 
                     num_feat = 50000, 
                     exclude = T,
                     union_features = union_features,
                     data_name = "con_50000")

# ungen 100 features
con_100000 <- getRand(quan_controls_gen, 
                      model = 'enet',
                      control = T,
                      rand_num = 10, 
                      num_feat = 100000, 
                      exclude = T,
                      union_features = union_features,
                      data_name = "con_100000")
##########
# combine data
##########

rand_cases_enet <- rbind(cases_10, cases_100, cases_500, cases_1000,
                    cases_5000, cases_10000, cases_50000, cases_100000)

rand_controls_enet <- rbind(con_10, con_100, con_500, con_1000,
                    con_10000,  con_50000, con_100000)

rand_cases_enet <- rand_cases_enet[order(rand_cases_enet$score, decreasing = T),]
rand_controls_enet <- rand_controls_enet[order(rand_controls_enet$score, decreasing = T),]



save.image('/home/benbrew/Desktop/temp/rand_model_enet.RData')

save.image('/home/benbrew/Desktop/temp/rand_model.RData')

