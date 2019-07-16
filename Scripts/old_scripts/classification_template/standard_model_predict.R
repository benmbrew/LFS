standard_model_predict <- function(data, ground_truth, partition, selected_features = NULL, NFOLDS = 5, N_CV_REPEATS = 2, run_ind = 0, type_measure = "auc") {

  # uses all samples in CV, as opposed to Other_Model_Predict2 which uses only patients
  
  # Checking for input parameters ---------------------------------
  require("glmnet")
  require("rpart")
  require("foreach")
  require("doParallel")
  require("ROCR")
  require("caret")
  require("nnet")
  source(paste0(test, "/multiclass.R"))
  
  stopifnot(dim(data)[1] > 10) # dimensions of data must be greater than 10
  stopifnot(dim(data)[2] > 10)
  stopifnot(typeof(partition) == "list")  # partition must be list
  #stopifnot(length(ground_truth) == dim(data)[1])  
  
  if (is.null(selected_features)) { # Null, which is the default, then selected_features
    # is all columns
    selected_features = 1:dim(data)[2]
  }
  
  type_measure.glmnet = "auc" # evaluate glmnet with auc
  type_measure.other_models = "ROC" # evaluate other models with ROC
  type_measure.performance = "auc"   # evalute perfomance auc
  
  if (length(partition[[run_ind]]$train_index) < 30) { # length of training set less that 30
    type_measure.glmnet = "class"
  }
  
  if (type_measure == "auc") {
    ground_truth <- as.factor(ground_truth)  
    #stopifnot(length(levels(ground_truth)) == 2)  # changes false to a and true to b
    levels(ground_truth) = c("a", "b", "c", "d", "e", "f", "g", "h")[1:length(levels(ground_truth))]
  } else if (type_measure == "acc") {
    ground_truth <- as.numeric(ground_truth)  # default is auc, so this is probably not relevant
    type_measure.glmnet = "class"
    type_measure.other_models = "ROC"
    type_measure.performance = "acc"    
  }

  #if (length(levels(ground_truth)) == 2) {
  #  type_family <- "binomial"
  #} else {
    type_family <- "multinomial" # packages will treat binomial and multinomial the same. 
  #}
  
  # Setup ---------------------------------
  # 0) Decision Tree
  dtree.gt <- ground_truth[partition[[run_ind]]$train_index] # get ground truth for training
  dtree.x <- data[partition[[run_ind]]$train_index, selected_features] # gets training x matrix 
  dtree.model <- rpart(
                  dtree.gt ~ .
                  , data=data.frame(dtree.x) 
                  , method='class'
                  , parms = list(split = "information") # splitting criteria
                  , control=rpart.control(xval = NFOLDS, minbucket = 5, cp = 0) # cross validates, 
                  #minbucket = the minimum number of observations in any terminal <leaf> node. 
                  #cp = complexity parameter. Any split that does not decrease the overall lack of fit by a factor of 
                  # cp is not attempted. overall R squared must increase by cp at each step. 
  )
  printcp(dtree.model)

  temp.predictions <- predict(# predicts on new data.
                      dtree.model
                      , newdata = data.frame(data[partition[[run_ind]]$test_index, selected_features])
                      , type = "prob" # gives probability
  )
  
  print(dim(temp.predictions))
  dtree.predictions <- temp.predictions
  #print(head(dtree.predictions))

  test.ground_truth <- ground_truth[partition[[run_ind]]$test_index] #get the real results for the current test
  test.ground_truth_1inN <- class.ind(ground_truth[partition[[run_ind]]$test_index]) # creates matrix of 0 and 1 for T and F
  #print(table(levels(test.ground_truth)[apply(dtree.predictions, 1, which.is.max)], test.ground_truth))
  # AUC
  #create ROCR prediction object
  temp.predictions <- prediction(dtree.predictions, test.ground_truth_1inN)
  dtree.test_auc <- unlist(slot(performance(temp.predictions, type_measure.performance), "y.values"))          
  print(paste("Decision Tree test AUC:", dtree.test_auc))  
  # Accuracy - the accuracy is the proportion of true results (both true positives and true negatives) 
  # among the total number of cases examined
  dtree.test_acc <- sum(levels(test.ground_truth)[apply(dtree.predictions, 1, which.is.max)] == test.ground_truth) / dim(dtree.predictions)[1]
  # basically anything, whether true or false, that the model got right/total observations. 
  print(paste("Decision Tree test acc:", dtree.test_acc)) 
  # Compute Confusion Matrix and Statistics
  #confusionMatrix(pred, truth)
  dtree.test_stats <- confusionMatrix(levels(test.ground_truth)[apply(dtree.predictions, 1, which.is.max)], test.ground_truth)
  print(dtree.test_stats)

  dtree = list( # saves everything in list
    # model = dtree.model
    # , predict = dtree.predictions
    # , 
    test_auc = dtree.test_auc
    , test_acc = dtree.test_stats$overall["Accuracy"]
    , test_stats = dtree.test_stats
  )

  rm(list = ls(pattern="dtree."))
  rm(list = ls(pattern="temp"))

  # 1) Elastic Net Logistic Regression
  elastic_net.cv_error = vector()
  elastic_net.cv_model = list()
  elastic_net.ALPHA <- c(1:9) / 10 # creates decimals
  # alpha controls the relative weighting of the ridge and lasso constraints. 
  
  # create error matrix for default twice and combine with rows
  # for opitmal alpha
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha_index in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha_index]] = cv.glmnet(data[partition[[run_ind]]$train_index, selected_features],
                                          ground_truth[partition[[run_ind]]$train_index],
                                          alpha = elastic_net.ALPHA[alpha_index] # first time with 0.1 and so on
                                          , type.measure = type_measure.glmnet
                                          , family = type_family
                                          , standardize = FALSE 
                                          , nfolds = NFOLDS # five folds
                                          , nlambda = 100
                                          , parallel = TRUE
      )
      elastic_net.cv_error[alpha_index] = min(elastic_net.cv_model[[alpha_index]]$cvm)
    }
    elastic_net.cv_error # stores 9 errors    
  }
  
  if (N_CV_REPEATS == 1) {
    temp.cv_error_mean = temp.cv_error_matrix
  } else {
    temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean) # take the mean of all the iterations  
  }
  
  stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
  # get index of best alpha (lowest alpha)
  temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
  print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
  
  # get optimal lambda
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  while (temp.non_zero_coeff < 3) { # loop runs initially because temp.non_zero coefficient <3 and then stops 
    # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
    # it they are never greater than 3, then the model does not converge. 
    elastic_net.cv_model = cv.glmnet(
                          data[partition[[run_ind]]$train_index, selected_features]
                          , ground_truth[partition[[run_ind]]$train_index]
                          , alpha = elastic_net.ALPHA[temp.best_alpha_index]
                          , type.measure = type_measure.glmnet
                          , family = type_family
                          , standardize=FALSE, 
                          , nlambda = 100
                          , nfolds = NFOLDS
                          , parallel = TRUE
    )
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) # get the min lambda
    # after 100 folds of cross validation
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] # number of non zero coefficients at that lambda    
    temp.loop_count = temp.loop_count + 1
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) # floor is opposite of ceiling. This just sets seed
    #print(paste0("seed: ", seed))
    if (temp.loop_count > 5) {
      print("diverged")
      temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
      break
    }
  }# while loop ends 
  print(temp.non_zero_coeff)  
  
  # Now that optimal level of alpha is chosen and the model does not diverge, run the model.
  elasticNet.model = glmnet(data[partition[[run_ind]]$train_index, selected_features], 
                            ground_truth[partition[[run_ind]]$train_index], 
                            alpha = elastic_net.ALPHA[temp.best_alpha_index],
                            standardize=FALSE,
                            nlambda = 100,
                            family = type_family)
  
  # This returns 100 prediction with 1-100 lambdas
  temp.predictions <- predict(elasticNet.model, data[partition[[run_ind]]$test_index, selected_features], type = "response")
  #print(dim(temp.predictions))
  elasticNet.predictions <- temp.predictions[,,  temp.min_lambda_index]  # this grabs the opitmal lambda 
  temp.l <- min(length(elastic_net.cv_model$lambda), length(elasticNet.model$lambda)) 
  stopifnot(elastic_net.cv_model$lambda[1:temp.l] == elasticNet.model$lambda[1:temp.l])  
  
  test.ground_truth <- ground_truth[partition[[run_ind]]$test_index]
  test.ground_truth_1inN <- class.ind(ground_truth[partition[[run_ind]]$test_index])
  #print(table(levels(test.ground_truth)[apply(elasticNet.predictions, 1, which.is.max)], test.ground_truth))
  # AUC
  #create ROCR prediction object
  temp.predictions <- prediction(elasticNet.predictions, test.ground_truth_1inN)
  elasticNet.test_auc <- unlist(slot(performance(temp.predictions, type_measure.performance), "y.values"))          
  print(paste("Elastic Net test AUC:", elasticNet.test_auc))  
  # Accuracy
  elasticNet.test_acc <- sum(levels(test.ground_truth)[apply(elasticNet.predictions, 1, which.is.max)] == test.ground_truth) / dim(elasticNet.predictions)[1]
  print(paste("Elastic Net test acc:", elasticNet.test_acc))  
  # Compute Confusion Matrix and Statistics
  #confusionMatrix(pred, truth)
  elasticNet.test_stats <- confusionMatrix(levels(test.ground_truth)[apply(elasticNet.predictions, 1, which.is.max)], test.ground_truth)
  print(elasticNet.test_stats)

  elasticNet = list(
    # model = elasticNet.model
    # , predict = elasticNet.predictions  
    # , 
    test_auc = elasticNet.test_auc
    , test_acc = elasticNet.test_stats$overall["Accuracy"]
    , test_stats = elasticNet.test_stats
    , alpha = elastic_net.ALPHA[temp.best_alpha_index]        
  )
  
  rm(list = ls(pattern="elasticNet."))
  rm(list = ls(pattern="elastic_net."))
  rm(list = ls(pattern="temp"))
  rm(list = ls(pattern="alpha."))  
  
  # 2) Lasso Logistic Regression  
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  while (temp.non_zero_coeff < 3) {      
    temp.cv_model = cv.glmnet(data[partition[[run_ind]]$train_index, selected_features]
                         , ground_truth[partition[[run_ind]]$train_index]
                         , alpha = 1
                         , type.measure = type_measure.glmnet
                         , family = type_family
                         , standardize = FALSE
                         , nlambda = 100
                         , nfolds = NFOLDS
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
  
  lasso.model = glmnet(data[partition[[run_ind]]$train_index, selected_features], 
                       ground_truth[partition[[run_ind]]$train_index], 
                       alpha = 1,
                       standardize=FALSE,
                       nlambda = 100,
                       family = type_family)
  
  temp.predictions <- predict(lasso.model, data[partition[[run_ind]]$test_index, selected_features], type = "response")  
  lasso.predictions <- temp.predictions[, , temp.min_lambda_index]  
  
  test.ground_truth <- ground_truth[partition[[run_ind]]$test_index]
  test.ground_truth_1inN <- class.ind(ground_truth[partition[[run_ind]]$test_index])
  #print(table(levels(test.ground_truth)[apply(lasso.predictions, 1, which.is.max)], test.ground_truth))
  # AUC
  #create ROCR prediction object
  temp.predictions <- prediction(lasso.predictions, test.ground_truth_1inN)
  lasso.test_auc <- unlist(slot(performance(temp.predictions, type_measure.performance), "y.values"))          
  print(paste("lasso test AUC:", lasso.test_auc))  
  # Accuracy
  lasso.test_acc <- sum(levels(test.ground_truth)[apply(lasso.predictions, 1, which.is.max)] == test.ground_truth) / dim(lasso.predictions)[1]
  print(paste("lasso test acc:", lasso.test_acc))  
  # Compute Confusion Matrix and Statistics
  #confusionMatrix(pred, truth)
  lasso.test_stats <- confusionMatrix(levels(test.ground_truth)[apply(lasso.predictions, 1, which.is.max)], test.ground_truth)
  print(lasso.test_stats)

  lasso = list(
    # model = lasso.model
    # , predict = lasso.predictions
    # , 
    test_auc = lasso.test_auc
    , test_acc = lasso.test_stats$overall["Accuracy"]
    , test_stats = lasso.test_stats
  )
  
  rm(list = ls(pattern="lasso."))  
  rm(list = ls(pattern="temp"))
  
  # 3) Ridge Regularized Logistic Regression  
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  while (temp.non_zero_coeff < 3) {      
    temp.cv_model = cv.glmnet(data[partition[[run_ind]]$train_index, selected_features]
                              , ground_truth[partition[[run_ind]]$train_index]
                              , alpha = 0
                              , type.measure = type_measure.glmnet
                              , family = type_family
                              , standardize = FALSE
                              , nlambda = 100
                              , nfolds = NFOLDS
                              , parallel = FALSE
    )
    temp.min_lambda_index = which(temp.cv_model$lambda == temp.cv_model$lambda.min)
    temp.non_zero_coeff = temp.cv_model$nzero[temp.min_lambda_index]
    temp.loop_count = temp.loop_count + 1
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    print(paste0("seed: ", seed))
    print(temp.loop_count)
    if (temp.loop_count > 5) {
      print("diverged")
      temp.min_lambda_index = 50
      break
    }        
  } 
  print(temp.non_zero_coeff)  
  
  ridge.model = glmnet(data[partition[[run_ind]]$train_index, selected_features], 
                       ground_truth[partition[[run_ind]]$train_index], 
                       alpha = 0,
                       standardize=FALSE,
                       nlambda = 100,
                       family = type_family)
  
  temp.predictions <- predict(ridge.model, data[partition[[run_ind]]$test_index, selected_features], type = "response")  
  ridge.predictions <- temp.predictions[, , temp.min_lambda_index]  
  
  test.ground_truth <- ground_truth[partition[[run_ind]]$test_index]
  test.ground_truth_1inN <- class.ind(ground_truth[partition[[run_ind]]$test_index])
  #print(table(levels(test.ground_truth)[apply(ridge.predictions, 1, which.is.max)], test.ground_truth))
  # AUC
  #create ROCR prediction object
  temp.predictions <- prediction(ridge.predictions, test.ground_truth_1inN)
  ridge.test_auc <- unlist(slot(performance(temp.predictions, type_measure.performance), "y.values"))          
  print(paste("Ridge test AUC:", ridge.test_auc))  
  # Accuracy
  ridge.test_acc <- sum(levels(test.ground_truth)[apply(ridge.predictions, 1, which.is.max)] == test.ground_truth) / dim(ridge.predictions)[1]
  print(paste("Ridge test acc:", ridge.test_acc))
  # Compute Confusion Matrix and Statistics
  #confusionMatrix(pred, truth)
  ridge.test_stats <- confusionMatrix(levels(test.ground_truth)[apply(ridge.predictions, 1, which.is.max)], test.ground_truth)
  print(ridge.test_stats)

  ridge = list(    
    # model = ridge.model
    # , predict = ridge.predictions   
    # , 
    test_auc = ridge.test_auc
    , test_acc = ridge.test_stats$overall["Accuracy"]
    , test_stats = ridge.test_stats
  )
  
  rm(list = ls(pattern="ridge."))  
  rm(list = ls(pattern="temp"))  
  
  # 4) Random Forest 
  other_model.y = as.factor(ground_truth[partition[[run_ind]]$train_index])
  levels(other_model.y) <- c("a", "b", "c", "d", "e", "f", "g", "h")[1:length(levels(ground_truth))]
  
  if (length(levels(ground_truth)) == 2) {
    summaryFunc <- twoClassSummary
  } else {
    summaryFunc <- multiClassSummary
  }

  # determines how you train the model.
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = min(10, NFOLDS), 
    classProbs = TRUE,     
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = summaryFunc)

  rf.model <- train(x = data[partition[[run_ind]]$train_index, selected_features]
                    , y = other_model.y
                    , method = "rf"
                    , trControl = fitControl                   
                    , verbose = FALSE
                    , metric = "ROC")

  rf.predictions <- predict(rf.model 
                        , newdata = data[partition[[run_ind]]$test_index, selected_features]
                        , type = "prob")
  
  test.ground_truth <- ground_truth[partition[[run_ind]]$test_index]
  test.ground_truth_1inN <- class.ind(ground_truth[partition[[run_ind]]$test_index])
  #print(table(levels(test.ground_truth)[apply(rf.predictions, 1, which.is.max)], test.ground_truth))
  # AUC
  #create ROCR prediction object
  temp.predict <- prediction(rf.predictions, test.ground_truth_1inN)  
  rf.test_auc <- unlist(slot(performance(temp.predict, "auc"), "y.values"))  
  print(paste("RF test AUC:", rf.test_auc))  
  # Accuracy
  rf.test_acc <- sum(levels(test.ground_truth)[apply(rf.predictions, 1, which.is.max)] == test.ground_truth) / dim(rf.predictions)[1]
  print(paste("RF test acc:", rf.test_acc))  
  # Compute Confusion Matrix and Statistics
  #confusionMatrix(pred, truth)
  rf.test_stats <- confusionMatrix(levels(test.ground_truth)[apply(rf.predictions, 1, which.is.max)], test.ground_truth)
  print(rf.test_stats)

  rf = list(
    # model = rf.model
    # , predict = rf.predictions
    # , 
    test_auc = rf.test_auc
    , test_acc = rf.test_stats$overall["Accuracy"]
    , test_stats = rf.test_stats
  )
  
  rm(list = ls(pattern="rf."))
  rm(list = ls(pattern="temp."))
  
  # 5) SVM with linear kernel
  
  svmLinear.model <- train(x = data[partition[[run_ind]]$train_index, selected_features]
                           , y = other_model.y
                           , method = "svmLinear"
                           , trControl = fitControl                   
                           , verbose = FALSE
                           , metric = type_measure.other_models # ROC instead of AUC                
  )
  
  svmLinear.predictions <- predict(svmLinear.model 
                               , newdata = data[partition[[run_ind]]$test_index, selected_features]
                               , type = "prob")

  test.ground_truth <- ground_truth[partition[[run_ind]]$test_index]
  test.ground_truth_1inN <- class.ind(ground_truth[partition[[run_ind]]$test_index])
  #print(table(levels(test.ground_truth)[apply(svmLinear.predictions, 1, which.is.max)], test.ground_truth))
  # AUC
  #create ROCR prediction object
  temp.predict <- prediction(svmLinear.predictions, test.ground_truth_1inN)  
  svmLinear.test_auc <- unlist(slot(performance(temp.predict, "auc"), "y.values"))  
  print(paste("SVM linear test AUC:", svmLinear.test_auc))     
  # Accuracy
  svmLinear.test_acc <- sum(levels(test.ground_truth)[apply(svmLinear.predictions, 1, which.is.max)] == test.ground_truth) / dim(svmLinear.predictions)[1]
  print(paste("SVM linear test acc:", svmLinear.test_acc))    
  # Compute Confusion Matrix and Statistics
  #confusionMatrix(pred, truth)
  svmLinear.test_stats <- confusionMatrix(levels(test.ground_truth)[apply(svmLinear.predictions, 1, which.is.max)], test.ground_truth)
  print(svmLinear.test_stats)

  svmLinear = list(    
    # model = svmLinear.model
    # , predict = svmLinear.predictions
    # , 
    test_auc = svmLinear.test_auc
    , test_acc = svmLinear.test_stats$overall["Accuracy"]
    , test_stats = svmLinear.test_stats
  )
  
  rm(list = ls(pattern="svmLinear."))
  rm(list = ls(pattern="temp."))
  
  # 6) SVM with radial kernel
  svmRadial.model <- train(x = data[partition[[run_ind]]$train_index, selected_features]
                           , y = other_model.y
                           , method = "svmRadial"
                           , trControl = fitControl
                           , verbose = FALSE
                           , metric = type_measure.other_models
  )
  
  svmRadial.predictions <- predict(svmRadial.model 
                               , newdata = data[partition[[run_ind]]$test_index, selected_features]
                               , type = "prob")

  test.ground_truth <- ground_truth[partition[[run_ind]]$test_index]
  test.ground_truth_1inN <- class.ind(ground_truth[partition[[run_ind]]$test_index])
  #print(table(levels(test.ground_truth)[apply(svmRadial.predictions, 1, which.is.max)], test.ground_truth))
  # AUC
  #create ROCR prediction object
  temp.predict <- prediction(svmRadial.predictions, test.ground_truth_1inN)  
  svmRadial.test_auc <- unlist(slot(performance(temp.predict, "auc"), "y.values"))  
  print(paste("SVM radial test AUC:", svmRadial.test_auc))     
  # Accuracy
  svmRadial.test_acc <- sum(levels(test.ground_truth)[apply(svmRadial.predictions, 1, which.is.max)] == test.ground_truth) / dim(svmRadial.predictions)[1]
  print(paste("SVM radial test acc:", svmRadial.test_acc))
  # Compute Confusion Matrix and Statistics
  #confusionMatrix(pred, truth)
  svmRadial.test_stats <- confusionMatrix(levels(test.ground_truth)[apply(svmRadial.predictions, 1, which.is.max)], test.ground_truth)
  print(svmRadial.test_stats)

  svmRadial = list(    
    # model = svmRadial.model
    # , predict = svmRadial.predictions
    # , 
    test_auc = svmRadial.test_auc
    , test_acc = svmRadial.test_stats$overall["Accuracy"]
    , test_stats = svmRadial.test_stats
  )
  
  rm(list = ls(pattern="svmRadial."))
  rm(list = ls(pattern="temp."))
  
  return (list(
    dtree = dtree
    , elasticNet = elasticNet
    , lasso = lasso
    , ridge = ridge
    , rf = rf
    , svmLinear = svmLinear
    , svmRadial = svmRadial
  ))
}