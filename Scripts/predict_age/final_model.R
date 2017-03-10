##########################################
# Final model
# 1) uses data quan_cases, quan_cases_gen, quan_cases_sen, quan_cases_sam, quan_cases_sen_gen, quan_cases_sam_gen
# 2) Will have function that subsets data, gets residuals, does classification threshold for 48,60,72, and gender. 
# 3) on training set either select variables through bumphunter or SVa
# 4) run Reg or Fac
# 5) RF or ENET

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
mod_data_folder <- paste0(data_folder, '/model_data')
results_folder <- paste0(project_folder, '/Results')
quan_folder <- paste0(results_folder, '/quan')
rand_folder <- paste0(results_folder, '/rand')
scripts_folder <- paste0(project_folder, '/Scripts')


##########
# load data batch data
##########
# read cases
quan_cases <- readRDS(paste0(mod_data_folder, '/quan_cases.rda'))

# read batch corrected data for gender
quan_cases_gen <- readRDS(paste0(mod_data_folder, '/quan_cases_gen.rda'))

# read batch corrected data for sentrix id and SAM
quan_cases_sen <- readRDS(paste0(mod_data_folder, '/quan_cases_sen.rda'))

quan_cases_sam <- readRDS(paste0(mod_data_folder, '/quan_cases_sam.rda'))

# read batch corrected data for sentrix id and SAM and gender!
quan_cases_sen_gen <- readRDS(paste0(mod_data_folder, '/quan_cases_sen_gen.rda'))

quan_cases_sam_gen <- readRDS(paste0(mod_data_folder, '/quan_cases_sam_gen.rda'))

# read controls
quan_controls <- readRDS(paste0(mod_data_folder, '/quan_controls.rda'))

quan_controls_gen <- readRDS(paste0(mod_data_folder, '/quan_controls_gen.rda'))

# load gene locations
cg_locations <- read.csv(paste0(mod_data_folder, '/cg_locations.csv'))


##########
# change column names of cg_locations for merging
##########
cg_locations$X <- NULL

# get functions 
source(paste0(scripts_folder, '/predict_age/final_functions.R'))

##########
# function to run models - subset, get residuals, get categorical, predict with regression and fac. 
##########
runModels <- function(data,
                      data_controls,
                      resid,
                      model,
                      gender,
                      bump_hunter,
                      num_it) 
  
{
    
  if (resid) {
    
    # get residuals
    data <- getResidual(data)
    
  } else {
    # get subset 
    data <- getSubset(data)
  }
  
  controls <- getSubset(data_controls)
  
  # get features
  all_feat <- colnames(data)[5:ncol(data)]
  
  # get holders
  model <- list()
  best_features <- list()
  importance <- list()
  test.predictions <- list()
  
  
  # split train and test at 70/30 and loop through iterations
  for(i in 1:num_it) {
    
    # set seed and get training index
    set.seed(i)
    train_index <- sample(nrow(data), nrow(data) *.70, replace = F)
    # creat training data frame
    dat_train <- data[train_index, ]

    # create testing data frame 
    dat_test <- data[-train_index,]

    ##########
    # run bumphunter to get features
    ##########
    # get balanced data from train
    dat <- getBal(dat_train, controls)
    
    # get bh results
    bh_dat <- bumpHunter(dat)
    
    # extract probe
    probe_features <- getProbe(bh_dat, cg_locations)
    model_feat <- as.character(probe_features$probe)
    
    # intersect model_feat with all_feat
    int_feat <- intersect(model_feat, all_feat)
      
    ##########
    # get model training and testing data
    ##########
    # get training outcome 
    y_train <- dat_train$age_diagnosis# gets training outcome
    y_test <- dat_test$age_diagnosis# testing outcome (real y)
    
    # get model data 
    dat_train <- dat_train[, c(int_feat)]
    dat_test <- dat_test[, int_feat]
    
    dims <- dim(dat_train)
    
    ##########
    # models
    ##########
    if (model == 'rf') {
      
      fitControl <- trainControl( 
        method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
        number = 4,      
        repeats = 1,
        allowParallel = TRUE
      )
      
      # mtry: Number of variables randomly sampled as candidates at each split.
      # ntree: Number of trees to grow.
      mtry <- sqrt(ncol(dat_train))
      tunegrid <- expand.grid(.mtry=mtry)
      
      model[[i]] <- train(x = dat_train
                    , y = y_train
                    , method = "rf"
                    , trControl = fitControl
                    , tuneGrid = tunegrid
                    , importance = T
                    , verbose = FALSE)
      
      temp <- varImp(model[[i]])[[1]]
      importance <- cbind(rownames(temp), temp$Overall)
      
      
      test.predictions <- predict(model 
                                  , newdata = dat_test)
      
    }
    
    if (model == 'enet') {
      
      N_CV_REPEATS = 2
      nfolds = 5
      
      ###### ENET
      # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
      # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
      elastic_net.cv_error = vector()
      elastic_net.cv_model = list()
      elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
      
      # set parameters for training model
      type_family <- 'gaussian'
      
      # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
      # or if you have a high number fo N_CV_REPEATS
      temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
        for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
        {      
          elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(dat_train)
                                                    , y =  y_train
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
        elastic_net.cv_model = cv.glmnet(as.matrix(dat_train)
                                         , y_train
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
      
      model[[i]] = glmnet(x = as.matrix(dat_train)
                          ,y =  y_train
                          ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                          ,standardize=FALSE
                          ,nlambda = 100
                          ,family = type_family)
      
      # This returns 100 prediction with 1-100 lambdas
      temp_test.predictions <- predict(model[[i]], as.matrix(dat_test),
                                       type = 'response')
      
      test.predictions[[i]] <- temp_test.predictions[, temp.min_lambda_index]
      
      importance[[i]] <- coef(model[[i]])
    

      
      
    }

  
  }

  return (model,test.predictions, y_test,  importance, dims) 
  
}



# data, data_controls, resid, model, gender, bump_hunter, num_it

quan_rf_

