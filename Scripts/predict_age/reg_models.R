# test models on this scripts

## This script will run all iterations of regression models 

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
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')

##########
# load cases 
##########

# load
cases_quan <- readRDS(paste0(model_data, '/cases_quan.rda'))

cases_sub_quan <- readRDS(paste0(model_data, '/cases_sub_quan.rda'))

valid_quan <- readRDS(paste0(model_data, '/valid_quan.rda'))


##########
# load controls 
##########

#load
controls_quan <- readRDS(paste0(model_data, '/controls_quan.rda'))

controls_full_quan <- readRDS(paste0(model_data, '/controls_full_quan.rda'))

##########
# check similarity of test data - names, structure
##########
# str(valid)
# str(cases)


##########
# load surveillance features or prediction features
##########
load(paste0(model_data, '/surv_10.RData'))
load(paste0(model_data, '/surv_20.RData'))
load(paste0(model_data, '/surv_30.RData'))

load(paste0(model_data, '/pred_10.RData'))
load(paste0(model_data, '/pred_20.RData'))
load(paste0(model_data, '/pred_30.RData'))


# model_dat <- cases_sub
# bh_features <- pred_bal_20_fwer
# classifier <- 'enet'
# gender <- T
# k <- 4
# seed_num <- 1
# i = 1
# valid_dat <- valid
# controls_dat <- controls

predAge <- function(model_dat,
                    controls_dat,
                    valid_dat,
                    classifier,
                    bh_features,
                    gender,
                    k,
                    seed_num)

{
  
  # create place list place holders
  model <- list()
  importance <- list()
  cases_cor <- list()
  controls_cor <- list()
  valid_cor <- list()
  alpha <- list()
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(model_dat))
  
  # get bumphunter features
  model_dat <- model_dat[, c('age_diagnosis', 'age_sample_collection', 'gender',  intersected_feats)]
  
  # 
  if(gender) {
    intersected_feats <- append('gender', intersected_feats)
    model_dat$gender <- as.factor(model_dat$gender)
    controls_dat$gender <- as.factor(controls_dat$gender)
    valid_dat$gender <- as.factor(valid_dat$gender)
    
  }

  set.seed(seed_num)
  # assign folds
  model_dat$folds <- sample(1:k, nrow(model_dat), replace = T)
  
  # now write forloop to 
  for (i in 1:k) { 
    
    # get x 
    train_index <- !grepl(i, model_dat$folds)
    test_index <- !train_index
    
    # get train and test data
    train_x <- model_dat[train_index, c('age_diagnosis', intersected_feats)]
    test_x <- model_dat[test_index, c('age_diagnosis', intersected_feats)]
    test_x_controls <- controls_dat[, c('age_sample_collection', intersected_feats)]
    test_x_valid <- valid_dat[, c('age_diagnosis', intersected_feats)]
    
    
    # get y
    train_y = as.numeric(train_x$age_diagnosis)
    test_y = as.numeric(test_x$age_diagnosis)
    test_y_controls <- as.numeric(test_x_controls$age_sample_collection)
    test_y_valid <- as.numeric(test_x_valid$age_diagnosis)
    
    if (classifier == 'rf'){
      
      # determines how you train the model.
      NFOLDS =  3
      fitControl <- trainControl( 
        method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
        number = min(10, NFOLDS),
        repeats = 3,
        allowParallel = TRUE
      )
      
      # mtry: Number of variables randomly sampled as candidates at each split.
      # ntree: Number of trees to grow.
      mtry <- sqrt(ncol(train_x))
      tunegrid <- expand.grid(.mtry=mtry)
      
      # train model
      model[[i]] <- train(x = train_x[, intersected_feats]
                          , y = train_y
                          , method = "rf"
                          , trControl = fitControl
                          , tuneGrid = tunegrid
                          , importance = T
                          , verbose = FALSE)
      
      temp <- varImp(model[[i]])[[1]]
      importance[[i]] <- cbind(rownames(temp), temp$Overall)
      
      
      test.predictions  <- predict(model[[i]] 
                                    , newdata = test_x[, intersected_feats])
      
      # get controls
      test.predictions_controls <- predict(model[[i]], 
                                           test_x_controls[, intersected_feats])
      
      # get valid
      test.predictions_valid <- predict(model[[i]], 
                                           test_x_valid[, intersected_feats])
      
      cases_cor[[i]] <- cor(test_y, test.predictions)
      
      # for each iteration, this should always be the same.
      controls_cor[[i]] <- cor(test_y_controls, test.predictions_controls)
      
      valid_cor[[i]] <- cor(test_y_valid, test.predictions_valid)
      
    }
    
    if (classifier == 'enet') {
      N_CV_REPEATS = 3
      nfolds = 3
      
      ###### ENET
      # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
      # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
      elastic_net.cv_error = vector()
      elastic_net.cv_model = list()
      elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
      
      # set parameters for training model
      type_family <- 'gaussian'
      
      if(gender){
        # make gender numeric
        train_x$gender <- as.numeric(train_x$gender)
        test_x$gender <- as.numeric(test_x$gender)
      }
      
      
      # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
      # or if you have a high number fo N_CV_REPEATS
      temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
        for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
        {      
          elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(train_x[, intersected_feats])
                                                    , y =  train_y
                                                    , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                    , type.measure = 'deviance'
                                                    , family = type_family
                                                    , standardize = FALSE 
                                                    , nfolds = nfolds 
                                                    , nlambda = 10
                                                    , parallel = TRUE
          )
          elastic_net.cv_error[alpha] = min(elastic_net.cv_model[[alpha]]$cvm)
        }
        elastic_net.cv_error # stores 9 errors    
      }
      
      if (N_CV_REPEATS == 1) {
        temp.cv_error_mean = temp.cv_error_matrix
      } else {
        temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean) # take the mean of the 5 iterations  
        # as your value for alpha
      }
      
      # stop if you did not recover error for any models 
      stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
      
      # get index of best alpha (lowest error) - alpha is values 0.1-0.9
      temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
      print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
      best_alpha <- elastic_net.ALPHA[temp.best_alpha_index]
      temp.non_zero_coeff = 0
      temp.loop_count = 0
      # loop runs initially because temp.non_zero coefficient <3 and then stops 
      # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
      # it they are never greater than 1, then the model does not converge. 
      while (temp.non_zero_coeff < 1) { 
        elastic_net.cv_model = cv.glmnet(x = as.matrix(train_x[, intersected_feats])
                                         , y =  train_y
                                         , alpha = elastic_net.ALPHA[temp.best_alpha_index]
                                         , type.measure = 'deviance'
                                         , family = type_family
                                         , standardize=FALSE
                                         , nlambda = 100
                                         , nfolds = nfolds
                                         , parallel = TRUE
        )
        
        # get optimal lambda - the tuning parameter for ridge and lasso
        # THIS IS IMPORTANT BECAUSE WHEN YOU TRAIN THE MODEL ON 100 SEPERATE VALUES OF LAMBDA
        # AND WHEN YOU TEST THE MODEL IT WILL RETURN PREDCITION FOR ALL THOSE VALUES (1-100). YOU NEED TO 
        # GRAB THE PREDICTION WITH SAME LAMBDA THAT YOU TRAINED ON. ITS ALL IN THE CODE, BUT JUST WANTED TO 
        # GIVE YOU REASONS
        temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
        
        # # number of non zero coefficients at that lambda    
        temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
        temp.loop_count = temp.loop_count + 1
        
        # set seed for next loop iteration
        as.numeric(Sys.time())-> t 
        set.seed((t - floor(t)) * 1e8 -> seed) 
        if (temp.loop_count > 10) {
          print("diverged")
          temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
          break
        }
      }# while loop ends 
      print(temp.non_zero_coeff)  
      
      model[[i]] = glmnet(x = as.matrix(train_x[, intersected_feats])
                          , y =  train_y
                          ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                          ,standardize=FALSE
                          ,nlambda = 100
                          ,family = type_family)
      
      # This returns 100 prediction with 1-100 lambdas
      temp_test.predictions <- predict(model[[i]], 
                                       data.matrix(test_x[, intersected_feats]),
                                       type = 'response')
      
      
      
      test.predictions <- temp_test.predictions[, temp.min_lambda_index]
      
      # get controls
      temp_test.predictions_controls <- predict(model[[i]], 
                                                data.matrix(test_x_controls[, intersected_feats]),
                                                type = 'response')
      
      
      
      test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
      
      
      # get validation
      temp_test.predictions_valid <- predict(model[[i]], 
                                             data.matrix(valid_dat[, intersected_feats]),
                                             type = 'response')
      
      
      
      test.predictions_valid  <- temp_test.predictions_valid[, temp.min_lambda_index]
      
      importance[[i]] <- coef(model[[i]])
      
      cases_cor[[i]] <- cor(test_y, test.predictions)
      
      # for each iteration, this should always be the same.
      controls_cor[[i]] <- cor(test_y_controls, test.predictions_controls)
      
      valid_cor[[i]] <- cor(test_y_valid, test.predictions_valid)
      
      alpha[[i]] <- best_alpha
    }
    
      
    
    print(i)

  }
  
  # get mean correlations in list
  mean_cor_cases <- mean(unlist(cases_cor))
  mean_cor_controls <- mean(unlist(controls_cor))
  mean_cor_valid <- mean(unlist(valid_cor))
  
 
  
  if(classifier == 'enet') {
    mean_alpha <- mean(unlist(alpha))
    
    result_table <- data.frame(mean_cor_cases = mean_cor_cases,
                               mean_cor_controls = mean_cor_controls,
                               mean_cor_valid = mean_cor_valid,
                               mean_alpha)

  } else {
    result_table <- data.frame(mean_cor_cases = mean_cor_cases,
                               mean_cor_controls = mean_cor_controls,
                               mean_cor_valid = mean_cor_valid)
  }
  
  return(list(result_table, importance, model))
}

############################################################
# surv
##########
# full, sub, 
##########

# full, sub, enet, 10, all, gender
f_s_surve_10_all <- predAge(model_dat = cases_quan, 
                            controls_dat = controls_full_quan, 
                            valid_dat = valid_quan, 
                            classifier = 'enet', 
                            bh_features = surv_ff_20_all, 
                            gender = F,
                            k = 4,
                            seed_num = 2)

f_s_surve_10_all[[1]]
#83 66 81
# full, sub, enet, 10, all, gender
f_s_surve_10_all_1 <- predAge(model_dat = cases_quan, 
                            controls_dat = controls_full_quan, 
                            valid_dat = valid_quan, 
                            classifier = 'enet', 
                            bh_features = surv_ff_20_sig, 
                            gender = F,
                            k = 4,
                            seed_num = 2)

f_s_surve_10_all_1[[1]]
# 85 68 82

# full, sub, enet, 10, all, gender
f_s_surve_10_all_2 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_full_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = surv_ff_30_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_2[[1]]

# 77 64 75
# full, sub, enet, 10, all, gender
f_s_surve_10_all_3 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_full_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = surv_ff_30_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_3[[1]]
# 84 74 79

# full, sub, enet, 10, all, gender
f_s_surve_10_all_4 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = surv_fs_10_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_4[[1]]
# 

# full, sub, enet, 10, all, gender
f_s_surve_10_all_5 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = surv_fs_10_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_5[[1]]

# full, sub, enet, 10, all, gender
f_s_surve_10_all_6 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = surv_fs_20_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_6[[1]]
# 80 63 75

# full, sub, enet, 10, all, gender
f_s_surve_10_all_7 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = surv_fs_20_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_7[[1]]
# 86 71 75

# full, sub, enet, 10, all, gender
f_s_surve_10_all_8 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = surv_fs_30_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_8[[1]]

# full, sub, enet, 10, all, gender
f_s_surve_10_all_9 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = surv_fs_30_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_9[[1]]



############################################################
# pred
##########
# full, sub, 
##########

# full, sub, enet, 10, all, gender
f_s_prede_10_all <- predAge(model_dat = cases_quan, 
                            controls_dat = controls_full_quan, 
                            valid_dat = valid_quan, 
                            classifier = 'enet', 
                            bh_features = pred_bal_full_10_sig, 
                            gender = F,
                            k = 4,
                            seed_num = 2)

f_s_prede_10_all[[1]]
# 88 83 88

# full, sub, enet, 10, all, gender
f_s_prede_10_all_1 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_full_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = pred_bal_full_20_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_1[[1]]
# 90 84 85

# full, sub, enet, 10, all, gender
f_s_prede_10_all_2 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_full_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = pred_bal_full_20_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_2[[1]]
# 88 82 86


# full, sub, enet, 10, all, gender
f_s_prede_10_all_3 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_full_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = pred_bal_full_30_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_3[[1]]
# 89 86 86

# full, sub, enet, 10, all, gender
f_s_prede_10_all_4 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = pred_bal_10_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_4[[1]]
# 90 83 87

# full, sub, enet, 10, all, gender
f_s_prede_10_all_5 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = pred_bal_20_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_5[[1]]
# 89 83 84

# full, sub, enet, 10, all, gender
f_s_prede_10_all_6 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = pred_bal_20_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_6[[1]]
# 91 83 88

# full, sub, enet, 10, all, gender
f_s_prede_10_all_7 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = pred_bal_30_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_7[[1]]
# 89 87 85

# full, sub, enet, 10, all, gender
f_s_prede_10_all_8 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = pred_bal_full_30_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_8[[1]]

# full, sub, enet, 10, all, gender
f_s_prede_10_all_9 <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'enet', 
                              bh_features = pred_bal_20_fwer, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_9[[1]]

########
# Random forest

############################################################
# surv
##########
# full, sub, 
##########

# full, sub, rf, 10, all, gender
f_s_surve_10_all_rf <- predAge(model_dat = cases_quan, 
                            controls_dat = controls_full_quan, 
                            valid_dat = valid_quan, 
                            classifier = 'rf', 
                            bh_features = surv_ff_20_all, 
                            gender = F,
                            k = 4,
                            seed_num = 2)

f_s_surve_10_all_rf[[1]]
#83 66 81
# full, sub, rf, 10, all, gender
f_s_surve_10_all_1_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_full_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = surv_ff_20_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_1_rf[[1]]
# 85 68 82

# full, sub, rf, 10, all, gender
f_s_surve_10_all_2_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_full_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = surv_ff_30_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_2_rf[[1]]

# 77 64 75
# full, sub, rf, 10, all, gender
f_s_surve_10_all_3_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_full_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = surv_ff_30_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_3_rf[[1]]
# 84 74 79

# full, sub, rf, 10, all, gender
f_s_surve_10_all_4_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = surv_fs_10_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_4_rf[[1]]
# 

# full, sub, rf, 10, all, gender
f_s_surve_10_all_5_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = surv_fs_10_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_5_rf[[1]]

# full, sub, rf, 10, all, gender
f_s_surve_10_all_6_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = surv_fs_20_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_6_rf[[1]]
# 80 63 75

# full, sub, rf, 10, all, gender
f_s_surve_10_all_7_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = surv_fs_20_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_7_rf[[1]]
# 86 71 75

# full, sub, rf, 10, all, gender
f_s_surve_10_all_8_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = surv_fs_30_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_8_rf[[1]]

# full, sub, rf, 10, all, gender
f_s_surve_10_all_9_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = surv_fs_30_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_surve_10_all_9_rf[[1]]



############################################################
# pred
##########
# full, sub, 
##########

# full, sub, rf, 10, all, gender
f_s_prede_10_all_rf <- predAge(model_dat = cases_quan, 
                            controls_dat = controls_full_quan, 
                            valid_dat = valid_quan, 
                            classifier = 'rf', 
                            bh_features = pred_bal_full_10_sig, 
                            gender = F,
                            k = 4,
                            seed_num = 2)

f_s_prede_10_all_rf[[1]]
# 88 83 88

# full, sub, rf, 10, all, gender
f_s_prede_10_all_1_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_full_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = pred_bal_full_20_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_1_rf[[1]]
# 90 84 85

# full, sub, rf, 10, all, gender
f_s_prede_10_all_2_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_full_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = pred_bal_full_20_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_2_rf[[1]]
# 88 82 86


# full, sub, rf, 10, all, gender
f_s_prede_10_all_3_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_full_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = pred_bal_full_30_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_3_rf[[1]]
# 89 86 86

# full, sub, rf, 10, all, gender
f_s_prede_10_all_4_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = pred_bal_10_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_4_rf[[1]]
# 90 83 87

# full, sub, rf, 10, all, gender
f_s_prede_10_all_5_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = pred_bal_20_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_5_rf[[1]]
# 89 83 84

# full, sub, rf, 10, all, gender
f_s_prede_10_all_6_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = pred_bal_20_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_6_rf[[1]]
# 91 83 88

# full, sub, rf, 10, all, gender
f_s_prede_10_all_7_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = pred_bal_30_all, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_7_rf[[1]]
# 89 87 85

# full, sub, rf, 10, all, gender
f_s_prede_10_all_8_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = pred_bal_full_30_sig, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_8_rf[[1]]

# full, sub, rf, 10, all, gender
f_s_prede_10_all_9_rf <- predAge(model_dat = cases_quan, 
                              controls_dat = controls_quan, 
                              valid_dat = valid_quan, 
                              classifier = 'rf', 
                              bh_features = pred_bal_20_fwer, 
                              gender = F,
                              k = 4,
                              seed_num = 2)

f_s_prede_10_all_9_rf[[1]]
