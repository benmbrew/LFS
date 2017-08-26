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
library(bumphunter)
library(sqldf)

registerDoParallel(1)
##########
# initialize folders
##########
home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
results_folder <- paste0(project_folder, '/Scripts/predict_age/Results/test_reg_results_05')


##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

# load in alpha mzx
alpha_max <- readRDS(paste0(model_data, '/alpha_max.rda'))

##########
# fixed variables
##########
method = 'quan'
k = 5
type = 'no_transform'

##########
# load data
##########
if (type == 'original') {
  
  betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new.rda')))
  betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new.rda'))) #34 449936
  betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new.rda'))) #35 449783
  
} else if (type == 'transform') {
  
  betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new.rda')))
  betaControls <- readRDS(paste0(model_data, paste0('/controls_transform','_' , method, '.rda'))) # 34 435813
  betaValid <- readRDS(paste0(model_data, paste0('/valid_transform','_' , method, '.rda'))) # 45 400842
  
} else {
  
  betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new.rda')))
  betaControls <- readRDS(paste0(model_data, paste0('/controls_no_transform','_' , method, '.rda'))) # 34 435813
  betaValid <- readRDS(paste0(model_data, paste0('/valid_no_transform','_' , method, '.rda')))# 45 400842
  
}

###########
# make id into ids
###########
colnames(betaCases)[1] <- 'ids'
colnames(betaControls)[1] <- 'ids'
colnames(betaValid)[1] <- 'ids'

###########
# get extra controls from 450k 
###########
betaControlsOld <- getControls(betaCases, mut = T)

###########
# get model data
###########
betaCases <- getModData(betaCases)

# get rid of cancer samples in controls 
betaControls <- betaControls[grepl('Unaffected', betaControls$cancer_diagnosis_diagnoses),]


# load cg_locations
cg_locations <- read.csv(paste0(model_data, 
                                '/cg_locations.csv'))

cg_locations$X <- NULL
##########
# get intersecting colnames and prepare data for modeling
##########

intersect_names <- Reduce(intersect, list(colnames(betaCases)[8:ncol(betaCases)], 
                                          colnames(betaControls)[8:ncol(betaControls)], 
                                          colnames(betaValid)[8:ncol(betaValid)]))
# organize each data set accordling

# cases
betaCases <- betaCases[, c('ids',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           intersect_names)]
# controls
betaControls <- betaControls[, c('ids',
                                 'age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender', 
                                 intersect_names)]

# controls
betaControlsOld <- betaControlsOld[, c('ids', 
                                       'age_diagnosis', 
                                       'age_sample_collection', 
                                       'cancer_diagnosis_diagnoses', 
                                       'gender', 
                                       intersect_names)]

#validation
betaValid <- betaValid[, c('ids', 
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           intersect_names)]

betaControls$ids <- as.character(betaControls$ids) 
betaControls$cancer_diagnosis_diagnoses <- as.character(betaControls$cancer_diagnosis_diagnoses) 
betaControls$gender <- as.character(betaControls$gender) 
betaControls$age_diagnosis <- as.numeric(as.character(betaControls$age_diagnosis))
betaControls$age_sample_collection <- as.numeric(as.character(betaControls$age_sample_collection))


# get controls full
betaControlsFull <- rbind(betaControls,
                          betaControlsOld)

# remove duplicates from betaControlsFull
length(which(duplicated(betaControlsFull$ids)))
betaControlsFull <- betaControlsFull[!duplicated(betaControlsFull$ids),]

##########
# remove samples with no age data
##########
betaControls <- betaControls[!is.na(betaControls$age_sample_collection),]
betaControlsOld <- betaControlsOld[!is.na(betaControlsOld$age_sample_collection),]
betaControlsFull <- betaControlsFull[!is.na(betaControlsFull$age_sample_collection),]


##########
# test for distribution sameness 
########## 
testKS(x = betaCases$age_sample_collection, y = betaControls$age_sample_collection)

# 
# betaCases <- betaCases[, c(1:3000, ncol(betaCases))]
# betaControls <- betaControls[, c(1:3000, ncol(betaControls))]
# betaControlsOld <- betaControlsOld[, c(1:3000, ncol(betaControlsOld))]
# betaControlsFull <- betaControlsFull[, c(1:3000, ncol(betaControlsFull))]
# betaValid <- betaValid[, c(1:3000, ncol(betaValid))]
##########
# run bumphunter
##########
# bh_feat <- bumpHunterSurv(dat_cases = betaCases, dat_controls = betaControls)
# bh_feat_old <- bumpHunterSurv(dat_cases = betaCases, dat_controls = betaControlsOld)
# bh_feat_full <- bumpHunterSurv(dat_cases = betaCases, dat_controls = betaControlsFull)

##########
# get features
##########
# bh_feat_all <- getProbe(bh_feat)
# bh_feat_tot <- getRun(bh_feat_all[[1]], run_num = .05)
# 
# bh_feat_all_old <- getProbe(bh_feat_old)
# bh_feat_tot_old <- getRun(bh_feat_all_old[[1]], run_num = .05)
# 
# bh_feat_all_full <- getProbe(bh_feat_full)
# bh_feat_tot_full <- getRun(bh_feat_all_full[[1]], run_num = .05)


###########
# save bh_features
###########
# saveRDS(bh_feat_tot, paste0(model_data, '/bh_feat_all', '_', method, '_', type))
# saveRDS(bh_feat_tot_old, paste0(model_data, '/bh_feat_old', '_', method, '_', type))
# saveRDS(bh_feat_tot_full, paste0(model_data, '/bh_feat_full', '_', method, '_', type))

bh_feat <- readRDS(paste0(model_data, '/bh_feat_all', '_', method, '_', type))
bh_feat_old <- readRDS(paste0(model_data, '/bh_feat_old', '_', method, '_', type))
bh_feat_full <- readRDS(paste0(model_data, '/bh_feat_full', '_', method, '_', type))

# get gender dummy variable
betaCases <- cbind(as.data.frame(class.ind(betaCases$gender)), betaCases)
betaControls <- cbind(as.data.frame(class.ind(betaControls$gender)), betaControls)
betaControlsOld <- cbind(as.data.frame(class.ind(betaControlsOld$gender)), betaControlsOld)
betaControlsFull <- cbind(as.data.frame(class.ind(betaControlsFull$gender)), betaControlsFull)
betaValid <- cbind(as.data.frame(class.ind(betaValid$gender)), betaValid)

#subset valid
betaValid <- betaValid[!betaValid$ids %in% betaCases$ids,]


##########
# test model
# ##########
# cases_dat <- betaCases
# controls_dat <- betaControls
# controls_dat_old <- betaControlsOld
# controls_dat_full <- betaControlsFull
# valid_dat <- betaValid
# bh_features <- bh_feat_tot
# alpha = 0.4
testModel <- function(cases_dat,
                      controls_dat,
                      controls_dat_old,
                      controls_dat_full,
                      valid_dat,
                      bh_features,
                      alpha)
{
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  bh_features <- append('M', bh_features)
  bh_features <- append('F', bh_features)
  
  intersected_feats <- intersect(bh_features, colnames(cases_dat))
  
  # get y
  train_y <- as.numeric(cases_dat$age_diagnosis)
  test_y_controls <- as.numeric(controls_dat$age_sample_collection)
  test_y_controls_old <- as.numeric(controls_dat_old$age_sample_collection)
  test_y_controls_full <- as.numeric(controls_dat_full$age_sample_collection)
  test_y_valid <- as.numeric(valid_dat$age_diagnosis)
  
  # get bumphunter features
  cases_dat <- cases_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  controls_dat_old <- controls_dat_old[, intersected_feats]
  controls_dat_full <- controls_dat_full[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
  N_CV_REPEATS = 2
  nfolds = 3
  
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  elastic_net.cv_model = cv.glmnet(x = as.matrix(cases_dat)
                                   , y =  train_y
                                   , alpha = alpha
                                   , type.measure = 'deviance'
                                   , family = 'gaussian'
                                   , standardize=FALSE
                                   , nlambda = 100
                                   , nfolds = nfolds
                                   , parallel = TRUE
  )
    
    
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    print(temp.non_zero_coeff)  
  
  
    model  = glmnet(x = as.matrix(cases_dat)
                  , y =  train_y
                  ,alpha = alpha
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = 'gaussian')
  
  # get controls
  temp_test.predictions_controls <- predict(model,
                                            data.matrix(controls_dat),
                                            type = 'response')



  test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
  
  # get controls
  temp_test.predictions_controls_old <- predict(model,
                                            data.matrix(controls_dat_old),
                                            type = 'response')
  
  
  
  test.predictions_controls_old <- temp_test.predictions_controls_old[, temp.min_lambda_index]
  
  # get controls full
  temp_test.predictions_controls_full <- predict(model,
                                                data.matrix(controls_dat_full),
                                                type = 'response')
  
  
  
  test.predictions_controls_full <- temp_test.predictions_controls_full[, temp.min_lambda_index]
  


  # get validation
  temp_test.predictions_valid <- predict(model,
                                         data.matrix(valid_dat),
                                         type = 'response')



  test.predictions_valid  <- temp_test.predictions_valid[, temp.min_lambda_index]

  importance <- coef(model)
  
  lambda_value <- elastic_net.cv_model$lambda.min
  
  # # for each iteration, this should always be the same.
  controls_cor <- cor(test_y_controls, test.predictions_controls)
  
  controls_cor_old <- cor(test_y_controls_old, test.predictions_controls_old)
  
  
  controls_cor_full <- cor(test_y_controls_full, test.predictions_controls_full)
  
  valid_cor  <- cor(test_y_valid, test.predictions_valid)

  return(list(controls_cor, controls_cor_old, controls_cor_full, valid_cor))
  
  
}

runModel <- function(bh_features) {
 
  alphas <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  
  result_list <- list()
  
  for(i in 1:length(alphas)) {
    
    alpha <- alphas[i]
    
    temp.result  <- testModel(cases_dat = betaCases,
                              controls_dat = betaControls, 
                              controls_dat_old = betaControlsOld,
                              controls_dat_full = betaControlsFull,
                              valid_dat = betaValid,
                              bh_features = bh_features,
                              alpha = alpha)
    
    temp.result_dat <- as.data.frame(t(do.call(rbind, temp.result)))
    temp.result_dat$alpha <- alpha
    colnames(temp.result_dat) <- c('controls', 'controls_old', 'controls_full', 'valid', 'alpha')
  
    result_list[[i]] <- temp.result_dat
    
    print(i)
  }
  
  result_table <- do.call(rbind, result_list)
  
  return(result_table)
  
   
}


# get results
results_norm <- runModel(bh_features = bh_feat)
results_old <- runModel(bh_features = bh_feat_old)
results_full <- runModel(bh_features = bh_feat_full)

# add columna and combine
results_norm$feat <- 'normal'
results_old$feat <- 'old'
results_full$feat <- 'full'

results <- rbind(results_norm,
                 results_old,
                 results_full)


saveRDS(results, paste0(results_folder, '/results_test_pipeline', '_', method, '_', type, '_', '05', '_' ,'.rda'))

##########
# calssification
##########

testModelFac <- function(cases_dat,
                         controls_dat,
                         controls_dat_full,
                         valid_dat,
                         bh_features,
                         cutoff,
                         alpha)
{
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(cases_dat))
  
  intersected_feats <- append('gender', intersected_feats)
  cases_dat$gender <- as.numeric(as.factor(cases_dat$gender))
  controls_dat$gender <- as.numeric(as.factor(controls_dat$gender))
  controls_dat_full$gender <- as.numeric(as.factor(controls_dat_full$gender))
  valid_dat$gender <- as.numeric(as.factor(valid_dat$gender))
  
  # # get y
  train_y <- factor(ifelse(cases_dat$age_diagnosis <= cutoff, 'yes', 'no'), 
                    levels = c('yes', 'no'))
  test_y_controls <- factor(ifelse(controls_dat$age_sample_collection <= cutoff, 'yes', 'no'), 
                            levels = c('yes', 'no'))
  test_y_controls_full <- factor(ifelse(controls_dat_full$age_sample_collection <= cutoff, 'yes', 'no'), 
                                 levels = c('yes', 'no'))
  test_y_valid <- factor(ifelse(valid_dat$age_sample_collection <= cutoff, 'yes', 'no'), 
                         levels = c('yes', 'no'))
  
  # get bumphunter features
  cases_dat <- cases_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  controls_dat_full <- controls_dat_full[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
  N_CV_REPEATS = 2
  nfolds = 3
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  elastic_net.cv_model = cv.glmnet(x = as.matrix(cases_dat)
                                   , y =  train_y
                                   , alpha = alpha
                                   , type.measure = 'auc'
                                   , family = type_family
                                   , standardize=FALSE
                                   , nlambda = 100
                                   , nfolds = nfolds
                                   , parallel = TRUE
  )
  
  
  temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
  
  # # number of non zero coefficients at that lambda    
  temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
  print(temp.non_zero_coeff)  
  
  
  model  = glmnet(x = as.matrix(cases_dat)
                  , y =  train_y
                  ,alpha = alpha
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  # get controls
  temp_test.predictions_controls <- predict(model,
                                            data.matrix(controls_dat),
                                            type = 'class')
  
  
  
  test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
  
  test.predictions_controls <- factor(test.predictions_controls, levels = c('yes', 'no'))
  
  
  test_stats_controls <- confusionMatrix(test_y_controls, test.predictions_controls)
  
  
  # get controls full
  temp_test.predictions_controls_full <- predict(model,
                                                 data.matrix(controls_dat_full),
                                                 type = 'class')
  
  
  
  test.predictions_controls_full <- temp_test.predictions_controls_full[, temp.min_lambda_index]
  
  test.predictions_controls_full <- factor(test.predictions_controls_full, levels = c('yes', 'no'))
  
  
  test_stats_controls_full <- confusionMatrix(test_y_controls_full, test.predictions_controls_full)
  
  
  # get validation
  temp_test.predictions_valid <- predict(model,
                                         data.matrix(valid_dat),
                                         type = 'class')
  
  
  
  test.predictions_valid  <- temp_test.predictions_valid[, temp.min_lambda_index]
  
  test.predictions_valid <- factor(test.predictions_valid, levels = c('yes', 'no'))
  
  
  test_stats_valid <- confusionMatrix(test_y_valid, test.predictions_valid)
  
  
  importance <- coef(model)
  
  lambda_value <- elastic_net.cv_model$lambda.min
  
  
  
  return(list(test_stats_controls, test_stats_controls_full, test_stats_valid))
  
  
}


testModelFac(betaCases,
             betaControls,
             betaControlsFull,
             betaValid,
             bh_feat_sig,
             cutoff = 60,
             alpha = 0.9)
