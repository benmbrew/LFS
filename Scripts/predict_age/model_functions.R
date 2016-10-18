#############################################################################################################
# This script will create functions to prepare data for modeling and functions for actual modeling.
##################################################################################################
# function that sunsets by age of diagnosis and mut
subsetDat <- function(model_data) {
  
  feature_names <- colnames(model_data)[5:ncol(model_data)]
  
  model_data_all <- model_data[!is.na(model_data$age_diagnosis),]

  model_data_mut <- model_data_all[model_data_all$p53_germline == 'Mut',]
  
  model_data_mut <- model_data_mut[, c('age_diagnosis', 'age_sample_collection', feature_names)]
  model_data_all <- model_data_all[, c('age_diagnosis', 'age_sample_collection', feature_names)]
  
  return(list(model_data_mut, model_data_all))
  
}

# function that takes a full data set and subsets it by bumphunter features
bhSubset <- function(model_data, 
                     bh_data) {
  
  bh_features <- bh_data$probe
  model_features <- colnames(model_data[[1]])
  bh_intersect <- intersect(bh_features, model_features)
  model_data_mut <- model_data[[1]][, c('age_diagnosis', 'age_sample_collection', bh_intersect)]
  model_data_all <- model_data[[2]][, c('age_diagnosis', 'age_sample_collection', bh_intersect)]
  
  
  return(list(model_data_mut, model_data_all))
}

# function that selects random features for model
getRand <- function(model_data) {
  
  feature_names <- colnames(model_data[[1]])[3:ncol(model_data[[1]])]
  rand_features <- sample(feature_names, 1000, replace = F)
  
  model_data_mut <- model_data[[1]][, c('age_diagnosis', 'age_sample_collection', rand_features)]
  model_data_all <- model_data[[2]][, c('age_diagnosis', 'age_sample_collection', rand_features)]
  
  return(list(model_data_mut, model_data_all))
  
}

# Function that takes model model_data and obstains residual features based on regressin each gene/probe on age of sample collection

getResidual <- function(model_data) {
  
  data_list <- list()
  
  for (data in 1:length(model_data)) {
    
    sub_data <- model_data[[data]]
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
    data_list[[data]] <- sub_data
  }
  
  return(data_list)
}

# Function that creates factor for a given threshold
makeFac <- function(model_data, 
                    threshold) {
  
  data_list <- list()
  for(data in 1:length(model_data)) {
    
    sub_data <- model_data[[data]]
    feature_names <- colnames(sub_data)[3:ncol(sub_data)]
    
    sub_data$age_diagnosis_fac <- as.integer(ifelse(sub_data$age_diagnosis <= threshold, 1, 2))
    sub_data$age_sample_fac <- as.integer(ifelse(sub_data$age_sample_collection <= threshold, 1, 2))
    
    
    sub_data <- sub_data[, c('age_diagnosis_fac', 'age_sample_fac', feature_names)]
    
    data_list[[data]] <- sub_data
    
  }
 
  return(data_list)
  
}


# Function that takes results list from regression and plots predictions against ground truth

plotModel <- function(result_list,
                      main1,
                      main2,
                      xlim,
                      ylim) {
  
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
  legend("bottomright", legend = paste0('# dims = ', result_list[[15]]), cex = 1, bty = 'n')
  
  
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


conMatrix <- function(results) {
  
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
                         iterations) {
  
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
    importance[[i]] <- cbind(rownames(temp), temp$Overall)
    
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
                        iterations) {
  
  model <- list()
  best_features <- list()
  importance <- list()
  train.mse <- list()
  test.mse <- list()
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
      
      train.predictions[[i]] <- predict(model[[i]] 
                                        , newdata = model_data[train_index, selected_features])
      
      
      
      
   
      train.ground_truth[[i]] <- model_data$age_diagnosis[train_index]
      test.ground_truth[[i]] <- model_data$age_diagnosis[-train_index]
      train.sample_collection[[i]] = model_data$age_sample_collection[train_index]
      test.sample_collection[[i]] = model_data$age_sample_collection[-train_index]
      train.mse[[i]] <- rmse(unlist(train.predictions[[i]]), unlist(train.ground_truth[[i]]))
      test.mse[[i]] <- rmse(unlist(test.predictions[[i]]), unlist(test.ground_truth[[i]]))
      
      
      
    
    print(i)
    
  }
  
  return(list(train.mse, test.mse,  train.predictions, test.predictions, train.ground_truth, test.ground_truth, train.sample_collection,
              test.sample_collection, test_acc, test_stats, test_acc_samp, test_stats_samp, model, importance, dims))
  
}


