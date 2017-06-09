### 
# this script will take the final model and test it on controls and validation.

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
# remove valid ids that are also in cases
##########
intersected_ids <- paste(intersect(valid_quan$ids, cases_quan$ids), collapse = '|')

valid_quan <- valid_quan[!grepl(intersected_ids, valid_quan$ids),]


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


# ##########
# # load surveillance features or prediction features
# ##########
load(paste0(model_data, '/surv_10.RData'))
load(paste0(model_data, '/surv_20.RData'))
load(paste0(model_data, '/surv_30.RData'))


load(paste0(model_data, '/pred_10.RData'))
load(paste0(model_data, '/pred_20.RData'))
load(paste0(model_data, '/pred_30.RData'))



model_dat <- cases_quan
bh_features <- surv_ff_20_sig
classifier <- 'enet'
gender <- F
k <- 4
seed_num <- 1
alpha_val <- 0.9


predAgeValid <- function(model_dat,
                         alpha_val,
                         classifier,
                         bh_features,
                         gender,
                         k,
                         seed_num)
  
{

  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(model_dat))
  
  # get bumphunter features
  model_dat <- model_dat[, c('age_diagnosis', 'age_sample_collection', 'gender',  intersected_feats)]
  
  # y
  y <- model_dat$age_diagnosis
  # 
  if(gender) {
    intersected_feats <- append('gender', intersected_feats)
    model_dat$gender <- as.factor(model_dat$gender)
    controls_dat$gender <- as.factor(controls_dat$gender)
    valid_dat$gender <- as.factor(valid_dat$gender)
    
  }
  
  valid_quan$gender <- as.factor(valid_quan$gender)
  controls_quan$gender <- as.factor(controls_quan$gender)
  controls_full_quan$gender <- as.factor(controls_full_quan$gender)
  
  y_valid <- valid_quan$age_diagnosis
  y_controls <- controls_quan$age_sample_collection
  y_controls_full <- controls_full_quan$age_sample_collection
  
    
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
      mtry <- sqrt(ncol(model_dat))
      tunegrid <- expand.grid(.mtry=mtry)
      
      # train model
      model  <- train(x = model_dat[, intersected_feats]
                          , y = y
                          , method = "rf"
                          , trControl = fitControl
                          , tuneGrid = tunegrid
                          , importance = T
                          , verbose = FALSE)
      
      temp <- varImp(model)[[1]]
      importance  <- cbind(rownames(temp), temp$Overall)
      
      
      test.predictions_valid  <- predict(model
                                   , newdata = valid_quan[, intersected_feats])
      
      # get controls
      test.predictions_controls <- predict(model, 
                                           controls_quan[, intersected_feats])
      
      # get controls
      test.predictions_controls_full <- predict(model, 
                                           controls_full_quan[, intersected_feats])
      
  
      # for each iteration, this should always be the same.
      controls_cor <- cor(test.predictions_controls, y_controls)
      controls_full_cor <- cor(test.predictions_controls, y_controls_full)
      valid_cor <- cor(y_valid, test.predictions_valid)
      
    }
    
    if (classifier == 'enet') {
      N_CV_REPEATS = 3
      nfolds = 3
      
      ###### ENET
      # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
      # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
      temp.non_zero_coeff = 0
      temp.loop_count = 0
      # set parameters for training model
      type_family <- 'gaussian'
      
      
      if(gender){
        # make gender numeric
        train_x$gender <- as.numeric(train_x$gender)
        test_x$gender <- as.numeric(test_x$gender)
      }
      
      
      # loop runs initially because temp.non_zero coefficient <3 and then stops 
      # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
      # it they are never greater than 1, then the model does not converge. 
      while (temp.non_zero_coeff < 1) { 
        elastic_net.cv_model = cv.glmnet(x = as.matrix(model_dat[, intersected_feats])
                                         , y =  y
                                         , alpha = alpha_val
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
      
      model = glmnet(x = as.matrix(model_dat[, intersected_feats])
                          , y =  y
                          ,alpha = alpha_val
                          ,standardize=FALSE
                          ,nlambda = 100
                          ,family = type_family)
      
      # This returns 100 prediction with 1-100 lambdas
      temp_test.predictions <- predict(model, 
                                       data.matrix(model_dat[, intersected_feats]),
                                       type = 'response')
      
      
      
      test.predictions <- temp_test.predictions[, temp.min_lambda_index]
      
      # get controls
      temp_test.predictions_controls <- predict(model, 
                                                data.matrix(controls_quan[, intersected_feats]),
                                                type = 'response')
      
      
      
      test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
      
      # get controls full
      temp_test.predictions_controls_full <- predict(model, 
                                                data.matrix(controls_full_quan[, intersected_feats]),
                                                type = 'response')
      
      
      
      test.predictions_controls_full <- temp_test.predictions_controls_full[, temp.min_lambda_index]
      
      
      # get validation
      temp_test.predictions_valid <- predict(model, 
                                             data.matrix(valid_quan[, intersected_feats]),
                                             type = 'response')
      
      
      
      test.predictions_valid  <- temp_test.predictions_valid[, temp.min_lambda_index]
      
      importance <- coef(model)
      
      train_cor <- cor(test.predictions, y)
      
      # for each iteration, this should always be the same.
      controls_cor <- cor(y_controls, test.predictions_controls)
      controls_cor_full <- cor(y_controls_full, test.predictions_controls_full)
      valid_cor <- cor(y_valid, test.predictions_valid)
      
      
    }
  
  result_table <- data.frame(controls_cor = controls_cor, 
                             controls_cor_full = controls_cor_full,
                             valid_cor = valid_cor)
  
    
    return(list(result_table, model, importance))

}

surv_ff_20_all_mod <- predAgeValid(model_dat = cases_quan, 
                               alpha_val = 0.9, 
                               classifier = 'enet', 
                               bh_features = surv_ff_20_all, 
                               gender = F, 
                               seed_num = 1)

surv_ff_20_all_mod[[1]]
  

surv_ff_20_sig_mod <- predAgeValid(model_dat = cases_quan, 
                               alpha_val = 0.9, 
                               classifier = 'enet', 
                               bh_features = surv_sf_20_fwer, 
                               gender = F, 
                               seed_num = 1)

surv_ff_20_sig_mod[[1]]

surv_ff_30_all_mod <- predAgeValid(model_dat = cases_quan, 
                                   alpha_val = 0.9, 
                                   classifier = 'enet', 
                                   bh_features = surv_ff_30_all, 
                                   gender = F, 
                                   seed_num = 1)

surv_ff_30_all_mod[[1]]


surv_ff_30_sig_mod <- predAgeValid(model_dat = cases_quan, 
                                   alpha_val = 0.9, 
                                   classifier = 'enet', 
                                   bh_features = surv_ff_30_sig, 
                                   gender = F, 
                                   seed_num = 1)

surv_ff_30_sig_mod[[1]]


