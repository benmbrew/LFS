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
cases <- readRDS(paste0(model_data, '/cases.rda'))

cases_sub <- readRDS(paste0(model_data, '/cases_sub.rda'))

##########
# load controls 
##########

#load
controls <- readRDS(paste0(model_data, '/controls.rda'))

controls_full <- readRDS(paste0(model_data, '/controls_full.rda'))


model_dat <- cases
num_feat <- 10
gender <- T
k <- 4
seed_num <- 1
i = 1

predAge <- function(model_dat,
                    num_feat,
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
  
  
  all_feat <- colnames(model_dat)[9:ncol(model_dat)]
  intersected_feats <- sample(all_feat, num_feat)
  
  # 
  if(gender) {
    intersected_feats <- append('gender', intersected_feats)
    model_dat$gender <- as.factor(model_dat$gender)

  }
  
  # get bumphunter features
  model_dat <- model_dat[, c('age_diagnosis', 'age_sample_collection',intersected_feats)]
  
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
    
    print(i)
    
  }
  correlation <- cor(unlist(test.predictions), unlist(y.test))
  
 return(correlation)
}

##########
# for loop that runs 3 times for each number of features
##########

getRand <- function(num_feat, iterations) 
{
  
  result <- list()
  
  for (seed_num in 1:iterations) {
    
    
    result[[seed_num]] <- predAge(model_dat = cases, 
                                  num_feat = num_feat,
                                  classifier = 'enet', 
                                  gender = T, 
                                  k = 3, 
                                  seed_num = seed_num)
  }
  
  mean_cor <- mean(unlist(result))
  sd_cor <- sd(unlist(result))
  n <- iterations
  
  # CI
  error <- qnorm(0.975)*sd_cor/sqrt(n)
  left <- mean_cor-error
  right <- mean_cor+error
  
  # make into data frame 
  result_data <- data.frame(mean_cor = mean_cor, upper = right, lower = left, n = n, num_feats = num_feat)
  
  
  return(result_data)
}


# 10 features
table_10 <- getRand(num_feat = 10, iterations = 20)

# 20 feautes
table_20 <- getRand(num_feat = 20, iterations = 20)

# 30 feautes
table_30 <- getRand(num_feat = 30, iterations = 20)

# 40 feautes
table_40 <- getRand(num_feat = 40, iterations = 20)

# 50 feautes
table_50 <- getRand(num_feat = 50, iterations = 20)

# 60 feautes
table_60 <- getRand(num_feat = 60, iterations = 20)

# 70 feautes
table_70 <- getRand(num_feat = 70, iterations = 20)

# 80 feautes
table_80 <- getRand(num_feat = 80, iterations = 20)

# 90 feautes
table_90 <- getRand(num_feat = 90, iterations = 20)

# 100 feautes
table_100 <- getRand(num_feat = 100, iterations = 20)

# 200 feautes
table_200 <- getRand(num_feat = 200, iterations = 20)

# 300 feautes
table_300 <- getRand(num_feat = 300, iterations = 20)

# 400 feautes
table_400 <- getRand(num_feat = 400, iterations = 20)

# 500 feautes
table_500 <- getRand(num_feat = 500, iterations = 20)

# 600 feautes
table_600 <- getRand(num_feat = 600, iterations = 20)

# 700 feautes
table_700 <- getRand(num_feat = 700, iterations = 20)

# 800 feautes
table_800 <- getRand(num_feat = 800, iterations = 20)

# 900 feautes
table_900 <- getRand(num_feat = 900, iterations = 20)

# 1000 feautes
table_1000 <- getRand(num_feat = 1000, iterations = 20)

# 2000 feautes
table_2000 <- getRand(num_feat = 2000, iterations = 20)

# 3000 feautes
table_3000 <- getRand(num_feat = 3000, iterations = 20)

# 4000 feautes
table_4000 <- getRand(num_feat = 4000, iterations = 20)

# 5000 feautes
table_5000 <- getRand(num_feat = 5000, iterations = 20)

# 6000 feautes
table_6000 <- getRand(num_feat = 6000, iterations = 20)

# 7000 feautes
table_7000 <- getRand(num_feat = 7000, iterations = 20)

# 8000 feautes
table_8000 <- getRand(num_feat = 8000, iterations = 20)

# 9000 feautes
table_9000 <- getRand(num_feat = 9000, iterations = 20)

# 10000 feautes
table_10000 <- getRand(num_feat = 10000, iterations = 20)

# 20000 feautes
table_20000 <- getRand(num_feat = 20000, iterations = 20)

# 30000 feautes
table_30000 <- getRand(num_feat = 30000, iterations = 20)

# 40000 feautes
table_40000 <- getRand(num_feat = 40000, iterations = 20)

# 50000 feautes
table_50000 <- getRand(num_feat = 50000, iterations = 20)

rand_results <- rbind(table_10, table_20, table_30, table_40, table_50,
                      table_60, table_70, table_80, table_90, table_100,
                      table_200, table_300, table_400, table_500, table_600,
                      table_700, table_800, table_900, table_1000)

# sor the data by score
rand_results <- rand_results[order(rand_results$mean_cor, decreasing = T),]

# make num_feat and mean_cor numeric
rand_results$mean_cor <- as.numeric(rand_results$mean_cor)
rand_results$num_feats <- as.numeric(rand_results$num_feats)


library(ggthemes)
ggplot(rand_results, aes(num_feats, mean_cor)) + 
  xlab('Number of Random Features') + ylab('Mean correlation') + ggtitle('Predictions and Real Ages') +
  geom_errorbar(aes(ymin = lower, ymax =upper), colour = 'black') + 
  geom_point(size = 2) + theme_calc()


saveRDS(rand_results, paste0(data_folder, '/rand_results.rda'))



