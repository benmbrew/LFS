#############################################################################################################
# This script will create functions to prepare data for modeling and functions for actual modeling.
##################################################################################################

# function that sunsets by age of diagnosis and mut
subsetDat <- function(model_data,
                      random = random,
                      num_feat = num_feat,
                      seed_num = seed_num) 
{
  
  feature_names <- colnames(model_data)[10:ncol(model_data)]
  
  
  if (random) {
    set.seed(seed_num)
    keep_index <- sample(feature_names, num_feat)
    model_data <- model_data[, c('age_diagnosis', 'age_sample_collection', keep_index) ]
    
  } else {
    model_data <- model_data[, c('age_diagnosis', 'age_sample_collection', feature_names)]
    
  }
  
  return(model_data)
  
}

runControl <- function(model_data, bh_data) 
{
  
  feature_names <- colnames(model_data)[8:ncol(model_data)]
  
  model_data <- model_data[, c('age_sample_collection', feature_names)]
  
  bh_features <- bh_data
  model_features <- colnames(model_data)
  bh_intersect <- intersect(bh_features, model_features)
  model_data <- model_data[, c('age_sample_collection', bh_intersect)]
  
  return(model_data)
}

# function that takes a full data set and subsets it by bumphunter features
bhSubset <- function(model_data, 
                     bh_data) 
{
  
  bh_features <- as.character(bh_data[,1])
  model_features <- colnames(model_data)
  bh_intersect <- intersect(bh_features, model_features)
  model_data <- model_data[, c('age_diagnosis', 'age_sample_collection', bh_intersect)]
  
  
  return(model_data)
}


# Function that takes model model_data and obstains residual features based on regressin each gene/probe on age of sample collection

getResidual <- function(model_data) 
{
  
  data_list <- list()
  
  sub_data <- model_data
  # subset by mut, and complete cases for age diagnosis and age sample collection
  sub_data <- sub_data[complete.cases(sub_data),]
  
  feature_names <- colnames(sub_data)[3:ncol(sub_data)]
  
  resid <- list()
  
  for (i in 3:ncol(sub_data)) {
    
    temp_response <- sub_data[, i]
    temp_var <- sub_data$age_sample_collection
    
    resid[[i]] <- lm(temp_response ~ temp_var)$residuals
    
    print(i)
    
  }
  
  resid_data <- as.data.frame(do.call('cbind', resid))
  sub_data <- cbind(sub_data$age_diagnosis, sub_data$age_sample_collection, resid_data)
  colnames(sub_data) <- c('age_diagnosis', 'age_sample_collection', feature_names)
  model_data <- sub_data
  
  return(model_data)
}

# Function that creates factor for a given threshold
makeFac <- function(model_data, 
                    threshold) 
{
  
  
  sub_data <- model_data
  feature_names <- colnames(sub_data)[3:ncol(sub_data)]
  
  sub_data$age_diagnosis_fac <- as.integer(ifelse(sub_data$age_diagnosis <= threshold, 1, 2))
  sub_data$age_sample_fac <- as.integer(ifelse(sub_data$age_sample_collection <= threshold, 1, 2))
  
  
  sub_data <- sub_data[, c('age_diagnosis_fac', 'age_sample_fac', feature_names)]
  
  model_data <- sub_data
  
  
  return(model_data)
  
}


# Function that takes results list from regression and plots predictions against ground truth

plotModel <- function(result_list,
                      main1,
                      main2,
                      xlim,
                      ylim) 
{
  
  # plot predictions against ground truth
  plot(unlist(result_list[[4]]), unlist(result_list[[6]]), 
       xlim = xlim,
       ylim = ylim,
       bty = 'n',
       col = adjustcolor('blue', alpha.f = 0.6),
       pch = 16,
       xlab = 'Predictions',
       ylab = 'Real Age of Diagnosis',
       main = main1)
  abline(0,1, lty = 3)
  corr <- round(cor(unlist(result_list[[4]]), unlist(result_list[[6]])), 2)
  legend("topleft", legend = paste0('correlation = ', corr), cex = 1, bty = 'n')
  legend("bottomright", legend = paste0('# samples = ', 
                                        result_list[[15]][1], 'and', result_list[[15]][2], 'features used' ), cex = 1, bty = 'n')
  
  
  # plot predictions against ground truth
  plot(unlist(result_list[[4]]), unlist(result_list[[8]]), 
       xlim = xlim,
       ylim = ylim,
       bty = 'n',
       col = adjustcolor('blue', alpha.f = 0.6),
       pch = 16,
       xlab = 'Predictions',
       ylab = 'Real Age of Sample Collection',
       main = main2)
  abline(0,1, lty = 3)
  corr <- round(cor(unlist(result_list[[4]]), unlist(result_list[[8]]), use = "complete.obs"), 2)
  legend("topleft", legend = paste0('correlation = ', corr), cex = 1, bty = 'n')
  legend("bottomright", legend = paste0('# dims = ', result_list[[15]]), cex = 1, bty = 'n')
  
}


conMatrix <- function(results) 
{
  
  # test acc for age of diagnosis
  acc_age <- mean(unlist(results[[7]]))
  
  # test acc for age of sample collection
  acc_samp <-mean(unlist(results[[9]]))
  
  # confustion matrix age of diagnosis 10
  iterations <- 10
  temp <- list()
  for (i in 1:10){
    temp[[i]] <- results[[8]][[i]]$table
  }
  mat <- unlist(temp)
  new_mat <- matrix(, 2, 2)
  
  mat_index <- seq(1, length(mat), 4)
  
  new_mat[1,1] <- sum(mat[mat_index])/iterations
  new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
  new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
  new_mat[2,2] <- sum(mat[mat_index + 3])/iterations
  
  return(list(new_mat, acc_age, acc_samp))
  
}


# function that runs random forest as a classification

rfPredictFac <- function(model_data,
                         cutoff,
                         iterations) 
{
  
  model <- list()
  best_features <- list()
  importance <- list()
  train.predictions <- list()
  test.predictions <- list()
  train.ground_truth <- list()
  test.ground_truth <- list()
  train.sample_collection <- list()
  test.sample_collection <- list()
  test_acc <- list()
  test_stats  <- list()
  test_acc_samp <- list()
  test_stats_samp <- list()
  
  
  dims <- dim(model_data)
  
  selected_features <- names(model_data[, 3:ncol(model_data)])
  
  for (i in 1:iterations) {
    
    set.seed(i)
    
    train_index <- sample(nrow(model_data), nrow(model_data) *cutoff, replace = F)
    
    # determines how you train the model.
    NFOLDS <- 2
    fitControl <- trainControl( 
      method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
      number = min(10, NFOLDS),
      classProbs = TRUE,
      repeats = 1,
      allowParallel = TRUE,
      summaryFunction = twoClassSummary
      
    )
    
    y = make.names(as.factor(model_data$age_diagnosis_fac[train_index]))
    
    # mtry: Number of variables randomly sampled as candidates at each split.
    # ntree: Number of trees to grow.
    
    mtry <- sqrt(ncol(model_data[train_index, selected_features]))
    tunegrid <- expand.grid(.mtry=mtry)
    
    model[[i]] <- train(x = model_data[train_index, selected_features]
                        , y = y
                        , method = "rf"
                        , trControl = fitControl
                        , tuneGrid = tunegrid
                        , importance = T
                        , verbose = FALSE)
    
    temp <- varImp(model[[i]])[[1]]
    importance[[i]] <- cbind(rownames(temp), temp$X1)
    
    test.predictions[[i]] <- predict(model[[i]] 
                                     , newdata = model_data[-train_index, selected_features]
                                     , type = "prob")
    
    train.predictions[[i]] <- predict(model[[i]] 
                                      , newdata = model_data[train_index, selected_features]
                                      ,type = "prob")
    
    
    
    
    
    train.ground_truth[[i]] <- as.factor(make.names(model_data$age_diagnosis_fac[train_index]))
    test.ground_truth[[i]] <- as.factor(make.names(model_data$age_diagnosis_fac[-train_index]))
    train.sample_collection[[i]] = as.factor(make.names(model_data$age_sample_fac[train_index]))
    train.sample_collection[[i]] <- factor(train.sample_collection[[i]], levels = c("X1", "X2"))
    # temp <- as.factor(make.names(model_data$age_sample_fac[-train_index]))
    test.sample_collection[[i]] = as.factor(make.names(model_data$age_sample_fac[-train_index]))
    test.sample_collection[[i]] <- factor(test.sample_collection[[i]], levels = c("X1", "X2"))
    
    
    # For age of diagnosis
    # Accuracy
    test_acc[[i]] <- sum(levels(test.ground_truth[[i]])[apply(test.predictions[[i]], 1, which.is.max)] == test.ground_truth[[i]]) / dim(test.predictions[[i]])[1]
    # Confustion Matrix
    test_stats[[i]] <- confusionMatrix(levels(test.ground_truth[[i]])[apply(test.predictions[[i]], 1, which.is.max)], test.ground_truth[[i]])
    # print(rf.test_stats)
    missing_ind <- !is.na(unlist(test.sample_collection[[i]]))
    test.sample_collection[[i]] <- unlist(test.sample_collection[[i]])[missing_ind]
    test.predictions[[i]] <- test.predictions[[i]][missing_ind,]
    # for age of collection
    # subset test.predictions by missing index in age.sample_collection, and remove the NAs in age.sample collection
    
    # Accuracy
    test_acc_samp[[i]] <- sum(levels(test.sample_collection[[i]])[apply(test.predictions[[i]], 1, which.is.max)] == test.sample_collection[[i]], na.rm = T) / dim(test.predictions[[i]])[1]
    # Confustion Matrix
    # subset to remove NAs
    test_stats_samp[[i]] <- confusionMatrix(levels(test.sample_collection[[i]])[apply(test.predictions[[i]], 1, which.is.max)], test.sample_collection[[i]])
    # print(rf.test_stats)
    
    
    
    print(i)
    
  }
  
  return(list(train.predictions, test.predictions, train.ground_truth, test.ground_truth, train.sample_collection,
              test.sample_collection, test_acc, test_stats, test_acc_samp, test_stats_samp, model, importance, dims))
  
}

# function that runs random forest as a regression

rfPredictReg <- function(model_data,
                         cutoff,
                         features,
                         iterations) 
{
  
  model <- list()
  best_features <- list()
  importance <- list()
  test.predictions <- list()
  test.ground_truth <- list()
  
  dims <- dim(model_data)
  selected_features <- names(model_data)[3:ncol(model_data)]
  
  
  for (i in 1:iterations) {
    
    set.seed(i)
    
    train_index <- sample(nrow(model_data), nrow(model_data) *cutoff, replace = F)
    
    NFOLDS <- 2
    fitControl <- trainControl( 
      method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
      number = min(10, NFOLDS),      
      repeats = 1,
      allowParallel = TRUE
    )
    
    y <- model_data$age_diagnosis[train_index]
    
    # mtry: Number of variables randomly sampled as candidates at each split.
    # ntree: Number of trees to grow.
    mtry <- sqrt(ncol(model_data[train_index, selected_features]))
    tunegrid <- expand.grid(.mtry=mtry)
    
    model[[i]] <- train(x = model_data[train_index, selected_features]
                        , y = y
                        , method = "rf"
                        , trControl = fitControl
                        , tuneGrid = tunegrid
                        , importance = T
                        , verbose = FALSE)
    
    temp <- varImp(model[[i]])[[1]]
    importance[[i]] <- cbind(rownames(temp), temp$Overall)
    
    
    test.predictions[[i]] <- predict(model[[i]] 
                                     , newdata = model_data[-train_index, selected_features])
    
    test.ground_truth[[i]] <- model_data$age_diagnosis[-train_index]
    
    
    
    
    print(i)
    
  }
  
  return(list(test.predictions, model,
              test.ground_truth,dims))
}


# now with my function
enetPredReg <- function(model_data,
                        N_CV_REPEATS,
                        nfolds,
                        cutoff,
                        iterations) {
  
  # get lists to store results
  model <- list()
  best_features <- list()
  importance <- list()
  test.predictions <- list()
  test.ground_truth <- list()
  
  # remove all rows with NAs - you can also impute
  # model_data <- model_data[complete.cases(model_data),]
  dims <- dim(model_data)
  
  for (i in 1:iterations){
    
    # set seed and get training index
    set.seed(i)
    train_index <- sample(nrow(model_data), nrow(model_data) *cutoff, replace = F)
    
    # get training outcome 
    y <- model_data$age_diagnosis[train_index] # gets training outcome
    selected_features <- names(model_data)[3:ncol(model_data)]
    
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
        elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(model_data[train_index, selected_features])
                                                  , y =  y
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
      elastic_net.cv_model = cv.glmnet(as.matrix(model_data[train_index, selected_features])
                                       , y
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
    
    model[[i]] = glmnet(x = as.matrix(model_data[train_index,  selected_features])
                        ,y =  y
                        ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                        ,standardize=FALSE
                        ,nlambda = 100
                        ,family = type_family)
    
    # This returns 100 prediction with 1-100 lambdas
    temp_test.predictions <- predict(model[[i]], as.matrix(model_data[-train_index, selected_features,]),
                                     type = 'response')
    
    test.predictions[[i]] <- temp_test.predictions[, temp.min_lambda_index]
    
    test.ground_truth[[i]] <- model_data$age_diagnosis[-train_index]
    
    print(i)
    
  }
  
  
  return(list(test.predictions, model,
              test.ground_truth, dims,
              temp.non_zero_coeff))
  
}


# variables to zero and therefore is a lot faster than random forest and other models
regPred <- function(model_data,
                    alpha,
                    nfolds,
                    cutoff,
                    iterations) {
  
  # get lists to store results
  model <- list()
  best_features <- list()
  importance <- list()
  test.predictions <- list()
  test.ground_truth <- list()
  
  # remove all rows with NAs - you can also impute
  dims <- dim(model_data)
  selected_features <- names(model_data)[3:ncol(model_data)]
  
  for (i in 1:iterations){
    
    # set seed and get training index
    set.seed(i)
    
    train_index <- sample(nrow(model_data), nrow(model_data)*cutoff, replace = F)
    
    # get training outcome 
    y <- model_data$age_diagnosis[train_index] # gets training outcome
    
    # set parameters for training model
    type_family <- 'gaussian'
    
    temp.non_zero_coeff = 0
    temp.loop_count = 0
    # loop runs initially because temp.non_zero coefficient <3 and then stops 
    # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
    # it they are never greater than 1, then the model does not converge. 
    while (temp.non_zero_coeff < 1) {      
      temp.cv_model = cv.glmnet(as.matrix(model_data[train_index, selected_features])
                                , y = y
                                , alpha = alpha
                                , type.measure = "deviance"
                                , family = type_family
                                , standardize = FALSE
                                , nlambda = 10
                                , nfolds = nfolds
                                , parallel = TRUE
      )
      temp.min_lambda_index = which(temp.cv_model$lambda == temp.cv_model$lambda.min)
      temp.non_zero_coeff = temp.cv_model$nzero[temp.min_lambda_index]
      temp.loop_count = temp.loop_count + 1
      as.numeric(Sys.time())-> t 
      set.seed((t - floor(t)) * 1e8 -> seed) 
      #print(paste0("seed: ", seed))
      if (temp.loop_count > 5) {
        print("diverged")
        temp.min_lambda_index = 50
        break
      }    
    }  
    print(temp.non_zero_coeff)  
    
    model[[i]] = glmnet(x = as.matrix(model_data[train_index,  selected_features])
                        ,y =  y
                        ,alpha = alpha # this means lasso, would be between 0-1 if enet
                        ,standardize=FALSE
                        ,nlambda = 100
                        ,family = type_family)
    
    # This returns 100 prediction with 1-100 lambdas
    temp_test.predictions <- predict(model[[i]], as.matrix(model_data[-train_index, selected_features,]),
                                     type = 'response')
    
    # get predictions with same lambda used in training 
    test.predictions[[i]] <- temp_test.predictions[, temp.min_lambda_index]
    
    test.ground_truth[[i]] <- model_data$age_diagnosis[-train_index]
    
    
    print(i)
    
  }
  
  
  return(list(test.predictions, model,
              test.ground_truth, dims, 
              temp.non_zero_coeff))
  
}



extractResults <- function (result_list,
                            data_name,
                            regularize) 
{
  
  # extract regression normal, correlation 
  temp.1 <- result_list
  reg_cor <- round(cor(unlist(temp.1[[1]]), unlist(temp.1[[3]])), 2)
  
  
  reg_cor <- as.data.frame(reg_cor)
  colnames(reg_cor) <- 'score'
  reg_cor$age <- 'regression'
  reg_cor$type <- 'normal'
  reg_cor$data <- data_name
  
  if(regularize) {
    reg_cor$features <- paste0(temp.1[[5]])
    
  } else {
    reg_cor$features <- as.numeric(strsplit(as.character(temp.1[[4]]), split = ' ')[[2]])
    
  }
  
  
  # # extract regression resid , resid_correlation 
  # 
  # temp.1 <- result_list[[2]]
  # reg_resid_cor <- round(cor(unlist(temp.1[[1]]), unlist(temp.1[[3]])), 2)
  # 
  # 
  # reg_resid_cor <- as.data.frame(reg_resid_cor)
  # colnames(reg_resid_cor) <- 'score'
  # reg_resid_cor$age <- 'regression'
  # reg_resid_cor$type <- 'resid'
  # reg_resid_cor$data <- data_name
  # reg_resid_cor$features <- paste0(result_list[[1]][[4]], collapse = '_')
  # 
  # # extract fac normal, acc for mut and all
  # temp.1 <- result_list[[3]]
  # temp.2 <- round(mean(unlist(temp.1[[9]])), 2)
  #   
  # 
  # 
  # fac_final <- as.data.frame(temp.2)
  # colnames(fac_final) <- 'score'
  # fac_final$age <- 48
  # fac_final$type <- 'normal'
  # fac_final$data <- data_name
  # fac_final$features <- paste0(unlist(result_list[[3]][13]), collapse = '_')
  # 
  # 
  # 
  # temp.1 <- result_list[[4]]
  # temp.2 <- round(mean(unlist(temp.1[[9]])),2)
  # 
  # fac_final_resid <- as.data.frame(unlist(temp.2))
  # colnames(fac_final_resid) <- 'score'
  # fac_final_resid$age <- 48
  # fac_final_resid$type <- 'resid'
  # fac_final_resid$data <- data_name
  # fac_final_resid$features <-paste0(unlist(result_list[[3]][13]), collapse = '_')
  # 
  # # combine all four result tables 
  # 
  # result_table <- rbind(reg_cor, reg_resid_cor)#, fac_final, fac_final_resid)
  
  return(reg_cor)
  
}


# make function that takes a result list and creates an object for plotting 
plotObject <- function (result_list, residual, p53_mut) 
{
  
  if (residual) {
    
    if (p53_mut) {
      plot_object <- result_list[[2]][[1]]
    } else {
      plot_object <- result_list[[2]][[2]]
    }
    
  } else {
    
    if (p53_mut) {
      plot_object <- result_list[[1]][[1]]
    } else {
      plot_object <- result_list[[1]][[2]]
    }
  }
  
  return(plot_object)
}

matObject <- function(result_list, age, residual, p53_mut) 
{
  
  if(residual) {
    
    if (age == 48) {
      if (p53_mut) {
        con_object <- result_list[[4]][[1]][[1]]
        
      } else {
        con_object <- result_list[[4]][[1]][[2]]
      }
    }
    
    if (age == 60) {
      if (p53_mut) {
        con_object <- result_list[[4]][[2]][[1]]
        
      } else {
        con_object <- result_list[[4]][[2]][[2]]
      }
    }
    
    if (age == 72) {
      if (p53_mut) {
        con_object <- result_list[[4]][[3]][[1]]
        
      } else {
        con_object <- result_list[[4]][[3]][[2]]
      }
    }
    
    if (age == 84) {
      if (p53_mut) {
        con_object <- result_list[[4]][[4]][[1]]
        
      } else {
        con_object <- result_list[[4]][[4]][[2]]
      }
    }
    
  } else {
    
    if (age == 48) {
      if (p53_mut) {
        con_object <- result_list[[3]][[1]][[1]]
        
      } else {
        con_object <- result_list[[3]][[1]][[2]]
      }
    }
    
    if (age == 60) {
      if (p53_mut) {
        con_object <- result_list[[3]][[2]][[1]]
        
      } else {
        con_object <- result_list[[3]][[2]][[2]]
      }
    }
    
    if (age == 72) {
      if (p53_mut) {
        con_object <- result_list[[3]][[3]][[1]]
        
      } else {
        con_object <- result_list[[3]][[3]][[2]]
      }
    }
    
    if (age == 84) {
      if (p53_mut) {
        con_object <- result_list[[3]][[4]][[1]]
        
      } else {
        con_object <- result_list[[3]][[4]][[2]]
      }
    }
    
  }
  
  
  return(con_object)
  
}


##########
# function takes the result table and creates a row and column variable from list "features" in model table"
##########
getDims <- function(results_table) 
{
  
  results_table$features <- as.character(results_table$features)
  
  for (i in 1:nrow(results_table)) {
    
    sub_dat <- results_table$features[[i]]
    if(grepl('c', sub_dat)) {
      sub_dat <- gsub('c(', '', sub_dat, fixed = T)
      sub_dat <- gsub(')', '', sub_dat, fixed = T)
      
      results_table$rows[i] <- strsplit(sub_dat, ',')[[1]][[1]]
      results_table$columns [i] <- strsplit(sub_dat, ',')[[1]][[2]]
      
    } else {
      results_table$rows[i] <- strsplit(sub_dat, ' ')[[1]][[1]]
      results_table$columns[i] <- strsplit(sub_dat, ' ')[[1]][[2]]
    }
    
  }
  
  results_table$features <- NULL
  return(results_table)
  
}
