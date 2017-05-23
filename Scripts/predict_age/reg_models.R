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
# load all data
##########
load(paste0(model_data, '/modal_feat_surv.RData'))

# remove everything but cases and controls
# remove cases and controls for now
rm(list=ls(pattern="even"))
rm(list=ls(pattern="uneven"))


##########
# load bh features for surveillance and pred data
##########
load(paste0(model_data, '/bh_feat_surv.RData'))
rm(list=ls(pattern="controls"))

# load(paste0(model_data, '/bh_feat_surv.RData'))

# model_dat <- quan_cases_full
# bh_features <- quan_even_full
# cases <- T
# gender <- T
# cv = 'fold' #other option here is "loo"
# k <- 4
# seed_num <- 1
# i = 1

predAge <- function(model_dat,
                    bh_features,
                    classifier,
                    gender,
                    k,
                    seed_num)

{
  
  # create place list place holders
  model <- list()
  importance <- list()
  test.predictions <- list()
  y.test <- list()
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(model_dat))
  
  # get bumphunter features
  model_dat <- model_dat[, c('age_diagnosis', 'age_sample_collection', 'gender',  intersected_feats)]
  # get features and remove unneeded columns
  
  if(gender) {
    intersected_feats <- append('gender', intersected_feats)
    model_dat$gender <- as.factor(model_dat$gender)
    
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
    
    # get y
    train_y = as.numeric(train_x$age_diagnosis)
    test_y = as.numeric(test_x$age_diagnosis)
    
    # determines how you train the model.
    NFOLDS =  2
    fitControl <- trainControl( 
      method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
      number = min(10, NFOLDS),
      repeats = 2,
      allowParallel = TRUE
    )
   
    
    
    if (classifier == 'rf'){
      
      
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
      
      
      test.predictions[[i]] <- predict(model[[i]] 
                                       , newdata = test_x[, intersected_feats])
      
      y.test[[i]] <- test_y
      
    }
    
    
    if(classifier == 'svm'){
     
      model[[i]] <- train(x = train_x[, intersected_feats]
                          , y = train_y
                          , method = 'svmLinear'
                          , trControl = fitControl                   
                          , verbose = FALSE)
      
      importance[[i]] <- 'svm_no_import'
      
      
      test.predictions[[i]] <- predict(model[[i]]
                                       , newdata = test_x[, intersected_feats])
      
      y.test[[i]] <- test_y
      
    }
    
    if(classifier == 'enet'){
      
      N_CV_REPEATS = 2
      nfolds = 3
      
      ###### ENET
      # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
      # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
      elastic_net.cv_error = vector()
      elastic_net.cv_model = list()
      elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
      
      # set parameters for training model
      type_family <- 'gaussian'
      
      # make gender numeric
      train_x$gender <- as.numeric(train_x$gender)
      test_x$gender <- as.numeric(test_x$gender)
      
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
                                      as.matrix(test_x[, intersected_feats]),
                                       type = 'response')
      
      test.predictions[[i]] <- temp_test.predictions[, temp.min_lambda_index]
      
      importance[[i]] <- coef(model[[i]])
      
      y.test[[i]] <- test_y
      
      
    }
    
    print(i)
    
    

  }
   return(list(y.test, test.predictions, importance, model))
  
}
# bh:  even_full, even_full_sig, even_full_fwer,
#      even_sub_bal, even_sub_bal_sig, even_sub_bal_fwer,
#      uneven_full, uneven_full_sig, uneven_full_fwer,
#      uneven_sub, uneven_sub_sig, uneven_sub_full 

# data: cases_full, cases_sub
##########
# apply model to raw
##########
# even, sub_bal
raw_sub_mod <- predAge(model_dat = raw_cases_full, 
                        bh_features = raw_even_sub_bal_sig, 
                        classifier = 'rf', 
                        gender = T, 
                        k = 3, 
                        seed_num = 1)

cor(unlist(raw_sub_mod[[1]]), unlist(raw_sub_mod[[2]]))

