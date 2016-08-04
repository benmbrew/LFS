##################################################################################################
# this script will read in different versions of subsetted full data and run random forest. 
library(dplyr)
library(stringr)
library(impute)
library(mlbench)
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(dplyr)
library(Metrics)


# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/Analyze')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

# Read in 3 different data sets 
full_data <- read.csv(paste0(data_folder, '/full_data.csv'), stringsAsFactors = F)
full_data_cor <- read.csv(paste0(data_folder, '/full_data_cor.csv'), stringsAsFactors = F)
full_data_rf <- read.csv(paste0(data_folder, '/full_data_rf.csv'), stringsAsFactors = F)

full_data$X <- NULL
full_data_cor$X <- NULL
full_data_rf$X <- NULL


# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(model_name, 
                       data,
                       subset, 
                       selected_features,
                       iterations) {
  
  model <- list()
  predictions <- list()
  mse <- list()
  importance <- list()
  test.ground_truth <- list()
  
  genes <- colnames(data)[28:ncol(data)]

  data <- data[, c(subset, genes)]
  
  # Try the model with all different selection of features based on number of missinginess. 
  data <- data[complete.cases(data),]
  
  obs <- nrow(data)
  
  # convert characters to factors 
  for ( i in 1:ncol(data)){
    
    if(typeof(data[,i]) == 'character' || grepl('num', names(data[i]))) {
      data[,i] <- as.factor(data[,i])
      print(i)
    } 
  }
  
  for (i in 1:iterations){
    
    set.seed(i)
    train_index <- sample(nrow(data), nrow(data) *.7, replace = F)
    
    # determines how you train the model.
    NFOLDS <- 2
    fitControl <- trainControl( 
      method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
      number = min(10, NFOLDS),      
      repeats = 1,
      allowParallel = TRUE
      #summaryFunction = summaryFunc
    )
    
    # mtry: Number of variables randomly sampled as candidates at each split.
    # ntree: Number of trees to grow.
    if (model_name == 'rf') {
      
      mtry <- sqrt(ncol(data))
      tunegrid <- expand.grid(.mtry=mtry)
      
      rf_y = data$age_diagnosis[train_index]
      
      
      model[[i]] <- train(x = data[train_index, c(selected_features, genes)]
                          , y = rf_y
                          , method = "rf"
                          , trControl = fitControl
                          , tuneGrid = tunegrid
                          , importance = T
                          , verbose = FALSE)
      
      importance[[i]] <- varImp(model[[i]])
      
      
      predictions[[i]] <- predict(model[[i]] 
                                  , newdata = data[-train_index, c(selected_features, genes)])
      
      test.ground_truth[[i]] <- data$age_diagnosis[-train_index]
      mse[[i]] <- rmse(unlist(predictions[[i]]), unlist(test.ground_truth[[i]]))
      
    }
    
    if (model_name == 'enet') {
      
      enet_y = data$age_diagnosis[train_index]
      N_CV_REPEATS <- 2
      elastic_net.cv_error = vector()
      elastic_net.cv_model = list()
      elastic_net.ALPHA <- c(1:9) / 10 # c
      
      temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
        
        for (alpha_index in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
        {      
          elastic_net.cv_model[[alpha_index]] = cv.glmnet(data.matrix(data[train_index, c(selected_features, genes)])
                                                          , enet_y
                                                          , alpha = elastic_net.ALPHA[alpha_index] # first time with 0.1 and so on
                                                          , type.measure = 'deviance'
                                                          , family = 'gaussian'
                                                          , standardize = FALSE 
                                                          , nfolds = 5# five folds
                                                          , nlambda = 100
                                                          , parallel = T
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
          data.matrix(data[train_index, c(selected_features, genes)])
          , enet_y
          , alpha = elastic_net.ALPHA[temp.best_alpha_index]
          , type.measure = 'deviance'
          , family = 'gaussian'
          , standardize=FALSE
          , nlambda = 100
          , nfolds = 3
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
      model[[i]] = glmnet(data.matrix(data[train_index, c(selected_features, genes)]), 
                          enet_y, 
                          alpha = elastic_net.ALPHA[temp.best_alpha_index],
                          standardize=FALSE,
                          nlambda = 100,
                          family = 'gaussian')
      
      # This returns 100 prediction with 1-100 lambdas
      temp_predictions <- predict(model[[i]], data.matrix(data[-train_index, c(selected_features, genes)]))
      #print(dim(temp.predictions))
      predictions[[i]] <- temp_predictions[ , temp.min_lambda_index]  # this grabs the opitmal lambda 
      temp.l <- min(length(elastic_net.cv_model$lambda), length(model[[i]]$lambda)) 
      stopifnot(elastic_net.cv_model$lambda[1:temp.l] == model[[i]]$lambda[1:temp.l])  
      test.ground_truth[[i]] <- data$age_diagnosis[-train_index]
      mse[[i]] <- rmse(unlist(predictions[[i]]), unlist(test.ground_truth[[i]]))
      
    }
    
    if(model_name == "lasso"){
      
      # 2) Lasso Logistic Regression 
      lasso_y <- data$age_diagnosis[train_index]
      temp.non_zero_coeff = 0
      temp.loop_count = 0
      while (temp.non_zero_coeff < 3) {     
        temp.cv_model = cv.glmnet(data.matrix(data[train_index, c(selected_features, genes)])
                                  , lasso_y
                                  , alpha = 1 # first time with 0.1 and so on
                                  , type.measure = 'deviance'
                                  , family = 'gaussian'
                                  , standardize = FALSE 
                                  , nfolds = 5# five folds
                                  , nlambda = 100
                                  , parallel = T
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
      
      model[[i]] = glmnet(data.matrix(data[train_index, c(selected_features, genes)]), 
                          lasso_y, 
                          alpha = 1,
                          standardize=FALSE,
                          nlambda = 100,
                          family = 'gaussian')
      
      
      # This returns 100 prediction with 1-100 lambdas
      temp_predictions <- predict(model[[i]], data.matrix(data[-train_index, c(selected_features, genes)]))
      #print(dim(temp.predictions))
      predictions[[i]] <- temp_predictions[ , temp.min_lambda_index]  # this grabs the opitmal lambda 
      
      test.ground_truth[[i]] <- data$age_diagnosis[-train_index]
      mse[[i]] <- rmse(unlist(predictions[[i]]), unlist(test.ground_truth[[i]]))
      
    }
    
    if(model_name == "svm") {
      
      
      svm_y = data$age_diagnosis[train_index]
      
      model[[i]] <- train(x = data[train_index, c(selected_features, genes)]
                          , y = svm_y
                          , method = "svmLinear"
                          , trControl = fitControl  
                          , importance = T
                          , verbose = FALSE
                          
      )
      
      predictions[[i]] <- predict(model[[i]] 
                                  , newdata = data[-train_index, selected_features])
      
      importance[[i]] <- varImp(model[[i]])
      
      test.ground_truth[[i]] <- data$age_diagnosis[-train_index]
      mse[[i]] <- rmse(unlist(predictions[[i]]), unlist(test.ground_truth[[i]]))
      
    }
    
    print(i)
    
  }
  
  return(list(mse, predictions, model, importance, test.ground_truth, obs))
  
}


########################################################################################################gdna.exon.intron

# variables missing
# gender 0
# gdna.base.change 164
# gdna.codon 164
# protein.codon.change 177
# gdna.exon.intron 492
# codon72.npro 517
# splice.delins.snv 519
# protein.codon.num 549
# mdm2.nG 652
##################################################################################################################3
# Random forest
#############################
# Full Data

# just methylation
rf_methyl <- predictAll(model_name = 'rf', 
                     data = full_data,
                     subset <- c("age_diagnosis"), 
                     selected_features = NULL, iterations = 50)

plot(unlist(rf_methyl[[2]]), unlist(rf_methyl[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# gender and gdna.base.change
rf_mut <- predictAll(model_name = 'rf', 
                     data = full_data,
                     subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                     selected_features = c("gender", "gdna.base.change"), iterations = 50)

plot(unlist(rf_mut[[2]]), unlist(rf_mut[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# add gdna.codon
rf_mut1 <- predictAll(model_name = 'rf', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 50)

plot(unlist(rf_mut1[[2]]), unlist(rf_mut1[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.change
rf_mut2 <- predictAll(model_name = 'rf', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), iterations = 50)

plot(unlist(rf_mut2[[2]]), unlist(rf_mut2[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add gdna.exon.intron
rf_mut3 <- predictAll(model_name = 'rf', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron"), iterations = 50)

plot(unlist(rf_mut3[[2]]), unlist(rf_mut3[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
rf_mut4 <- predictAll(model_name = 'rf', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro"), iterations = 50)

plot(unlist(rf_mut4[[2]]), unlist(rf_mut4[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add splice.delins.snv
rf_mut5 <- predictAll(model_name = 'rf', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), iterations = 50)

plot(unlist(rf_mut5[[2]]), unlist(rf_mut5[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add protein.codon.num
rf_mut6 <- predictAll(model_name = 'rf', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), iterations = 50)

plot(unlist(rf_mut6[[2]]), unlist(rf_mut6[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add mdm2.nG
rf_mut7 <- predictAll(model_name = 'rf', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), iterations = 50)

plot(unlist(rf_mut7[[2]]), unlist(rf_mut7[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


#############################
# Corelation data

# just methylation
rf_methyl_cor <- predictAll(model_name = 'rf', 
                        data = full_data_cor,
                        subset <- c("age_diagnosis"), 
                        selected_features = NULL, iterations = 50)

plot(unlist(rf_methyl_cor[[2]]), unlist(rf_methyl_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# gender and gdna.base.change
rf_mut_cor <- predictAll(model_name = 'rf', 
                     data = full_data_cor,
                     subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                     selected_features = c("gender", "gdna.base.change"), iterations = 50)

plot(unlist(rf_mut_cor[[2]]), unlist(rf_mut_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# add gdna.codon
rf_mut1_cor <- predictAll(model_name = 'rf', 
                      data = full_data_cor,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 50)

plot(unlist(rf_mut1_cor[[2]]), unlist(rf_mut1_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.change
rf_mut2_cor <- predictAll(model_name = 'rf', 
                      data = full_data_cor,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), iterations = 50)

plot(unlist(rf_mut2_cor[[2]]), unlist(rf_mut2_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add gdna.exon.intron
rf_mut3_cor <- predictAll(model_name = 'rf', 
                      data = full_data_cor,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron"), iterations = 50)

plot(unlist(rf_mut3_cor[[2]]), unlist(rf_mut3_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
rf_mut4_cor <- predictAll(model_name = 'rf', 
                      data = full_data_cor,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro"), iterations = 50)

plot(unlist(rf_mut4_cor[[2]]), unlist(rf_mut4_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add splice.delins.snv
rf_mut5_cor <- predictAll(model_name = 'rf', 
                      data = full_data_cor,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), iterations = 50)

plot(unlist(rf_mut5_cor[[2]]), unlist(rf_mut5_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add protein.codon.num
rf_mut6_cor <- predictAll(model_name = 'rf', 
                      data = full_data_cor,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), iterations = 50)

plot(unlist(rf_mut6_cor[[2]]), unlist(rf_mut6_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add mdm2.nG
rf_mut7_cor <- predictAll(model_name = 'rf', 
                      data = full_data_cor,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), iterations = 50)

plot(unlist(rf_mut7_cor[[2]]), unlist(rf_mut7_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)



#############################
# recursive elimination data

# just methylation
rf_methyl_reg <- predictAll(model_name = 'rf', 
                            data = full_data_rf,
                            subset <- c("age_diagnosis"), 
                            selected_features = NULL, iterations = 50)

plot(unlist(rf_methyl_reg[[2]]), unlist(rf_methyl_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# gender and gdna.base.change
rf_mut_reg <- predictAll(model_name = 'rf', 
                         data = full_data_rf,
                         subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                         selected_features = c("gender", "gdna.base.change"), iterations = 50)

plot(unlist(rf_mut_reg[[2]]), unlist(rf_mut_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# add gdna.codon
rf_mut1_reg <- predictAll(model_name = 'rf', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 50)

plot(unlist(rf_mut1_reg[[2]]), unlist(rf_mut1_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.change
rf_mut2_reg <- predictAll(model_name = 'rf', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), iterations = 50)

plot(unlist(rf_mut2_reg[[2]]), unlist(rf_mut2_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add gdna.exon.intron
rf_mut3_reg <- predictAll(model_name = 'rf', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron"), iterations = 50)

plot(unlist(rf_mut3_reg[[2]]), unlist(rf_mut3_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
rf_mut4_reg <- predictAll(model_name = 'rf', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro"), iterations = 50)

plot(unlist(rf_mut4_reg[[2]]), unlist(rf_mut4_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add splice.delins.snv
rf_mut5_reg <- predictAll(model_name = 'rf', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), iterations = 50)

plot(unlist(rf_mut5_reg[[2]]), unlist(rf_mut5_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add protein.codon.num
rf_mut6_reg <- predictAll(model_name = 'rf', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), iterations = 50)

plot(unlist(rf_mut6_reg[[2]]), unlist(rf_mut6_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add mdm2.nG
rf_mut7_reg <- predictAll(model_name = 'rf', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), iterations = 50)

plot(unlist(rf_mut7_reg[[2]]), unlist(rf_mut7_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# combine all tests 

rf_all <- rbind (
  append('just_methyl', c(mean(unlist(rf_methyl[[1]])), rf_methyl[[6]], 'all_data')),
  append('gender_and_gdna.base.change', c(mean(unlist(rf_mut[[1]])), rf_mut[[6]], 'all_data')),
  append('gdna.codon', c(mean(unlist(rf_mut1[[1]])), rf_mut1[[6]], 'all_data')),
  append('protein.codon.change', c(mean(unlist(rf_mut2[[1]])), rf_mut2[[6]], 'all_data')),
  append('gdna.exon.intron', c(mean(unlist(rf_mut3[[1]])), rf_mut3[[6]], 'all_data')),
  append('codon72.npro', c(mean(unlist(rf_mut4[[1]])), rf_mut4[[6]], 'all_data')),
  append('splice.delins.snv', c(mean(unlist(rf_mut5[[1]])), rf_mut5[[6]], 'all_data')),
  append('protein.codon.num', c(mean(unlist(rf_mut6[[1]])), rf_mut6[[6]], 'all_data')),
  append('mdm2.nG', c(mean(unlist(rf_mut7[[1]])), rf_mut7[[6]], 'all_data')),
  
  append('just_methyl', c(mean(unlist(rf_methyl_cor[[1]])), rf_methyl_cor[[6]], 'corr')),
  append('gender_and_gdna.base.change', c(mean(unlist(rf_mut_cor[[1]])), rf_mut_cor[[6]], 'corr')),
  append('gdna.codon', c(mean(unlist(rf_mut1_cor[[1]])), rf_mut1_cor[[6]], 'corr')),
  append('protein.codon.change', c(mean(unlist(rf_mut2_cor[[1]])), rf_mut2_cor[[6]], 'corr')),
  append('gdna.exon.intron', c(mean(unlist(rf_mut3_cor[[1]])), rf_mut3_cor[[6]], 'corr')),
  append('codon72.npro', c(mean(unlist(rf_mut4_cor[[1]])), rf_mut4_cor[[6]], 'corr')),
  append('splice.delins.snv', c(mean(unlist(rf_mut5_cor[[1]])), rf_mut5_cor[[6]], 'corr')),
  append('protein.codon.num', c(mean(unlist(rf_mut6_cor[[1]])), rf_mut6_cor[[6]], 'corr')),
  append('mdm2.nG', c(mean(unlist(rf_mut7_cor[[1]])), rf_mut7_cor[[6]], 'corr')),
  
  append('just_methyl', c(mean(unlist(rf_methyl_reg[[1]])), rf_methyl_reg[[6]], 'recursive')),
  append('gender_and_gdna.base.change', c(mean(unlist(rf_mut_reg[[1]])), rf_mut_reg[[6]], 'recursive')),
  append('gdna.codon', c(mean(unlist(rf_mut1_reg[[1]])), rf_mut1_reg[[6]], 'recursive')),
  append('protein.codon.change', c(mean(unlist(rf_mut2_reg[[1]])), rf_mut2_reg[[6]], 'recursive')),
  append('gdna.exon.intron', c(mean(unlist(rf_mut3_reg[[1]])), rf_mut3_reg[[6]], 'recursive')),
  append('codon72.npro', c(mean(unlist(rf_mut4_reg[[1]])), rf_mut4_reg[[6]], 'recursive')),
  append('splice.delins.snv', c(mean(unlist(rf_mut5_reg[[1]])), rf_mut5_reg[[6]], 'recursive')),
  append('protein.codon.num', c(mean(unlist(rf_mut6_reg[[1]])), rf_mut6_reg[[6]], 'recursive')),
  append('mdm2.nG', c(mean(unlist(rf_mut7_reg[[1]])), rf_mut7_reg[[6]], 'recursive'))
)



##################################################################################################################3
# Elastic net
#############################
# Full Data

# just methylation
enet_methyl <- predictAll(model_name = 'enet', 
                        data = full_data,
                        subset <- c("age_diagnosis"), 
                        selected_features = NULL, iterations = 50)

plot(unlist(enet_methyl[[2]]), unlist(enet_methyl[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# gender and gdna.base.change
enet_mut <- predictAll(model_name = 'enet', 
                     data = full_data,
                     subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                     selected_features = c("gender", "gdna.base.change"), iterations = 50)

plot(unlist(enet_mut[[2]]), unlist(enet_mut[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# add gdna.codon
enet_mut1 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 50)

plot(unlist(enet_mut1[[2]]), unlist(enet_mut1[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.change
enet_mut2 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), iterations = 50)

plot(unlist(enet_mut2[[2]]), unlist(enet_mut2[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add gdna.exon.intron
enet_mut3 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron"), iterations = 50)

plot(unlist(enet_mut3[[2]]), unlist(enet_mut3[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
enet_mut4 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro"), iterations = 50)

plot(unlist(enet_mut4[[2]]), unlist(enet_mut4[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add splice.delins.snv
enet_mut5 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), iterations = 50)

plot(unlist(enet_mut5[[2]]), unlist(enet_mut5[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add protein.codon.num
enet_mut6 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), iterations = 50)

plot(unlist(enet_mut6[[2]]), unlist(enet_mut6[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add mdm2.nG
enet_mut7 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), iterations = 50)

plot(unlist(enet_mut7[[2]]), unlist(enet_mut7[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


#############################
# Corelation data

# just methylation
enet_methyl_cor <- predictAll(model_name = 'enet', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis"), 
                            selected_features = NULL, iterations = 50)

plot(unlist(enet_methyl_cor[[2]]), unlist(enet_methyl_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# gender and gdna.base.change
enet_mut_cor <- predictAll(model_name = 'enet', 
                         data = full_data_cor,
                         subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                         selected_features = c("gender", "gdna.base.change"), iterations = 50)

plot(unlist(enet_mut_cor[[2]]), unlist(enet_mut_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# add gdna.codon
enet_mut1_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 50)

plot(unlist(enet_mut1_cor[[2]]), unlist(enet_mut1_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.change
enet_mut2_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), iterations = 50)

plot(unlist(enet_mut2_cor[[2]]), unlist(enet_mut2_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add gdna.exon.intron
enet_mut3_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron"), iterations = 50)

plot(unlist(enet_mut3_cor[[2]]), unlist(enet_mut3_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
enet_mut4_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro"), iterations = 50)

plot(unlist(enet_mut4_cor[[2]]), unlist(enet_mut4_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add splice.delins.snv
enet_mut5_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), iterations = 50)

plot(unlist(enet_mut5_cor[[2]]), unlist(enet_mut5_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add protein.codon.num
enet_mut6_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), iterations = 50)

plot(unlist(enet_mut6_cor[[2]]), unlist(enet_mut6_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add mdm2.nG
enet_mut7_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), iterations = 50)

plot(unlist(enet_mut7_cor[[2]]), unlist(enet_mut7_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)



#############################
# recursive elimination data

# just methylation
enet_methyl_reg <- predictAll(model_name = 'enet', 
                            data = full_data_rf,
                            subset <- c("age_diagnosis"), 
                            selected_features = NULL, iterations = 50)

plot(unlist(enet_methyl_reg[[2]]), unlist(enet_methyl_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# gender and gdna.base.change
enet_mut_reg <- predictAll(model_name = 'enet', 
                         data = full_data_rf,
                         subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                         selected_features = c("gender", "gdna.base.change"), iterations = 50)

plot(unlist(enet_mut_reg[[2]]), unlist(enet_mut_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# add gdna.codon
enet_mut1_reg <- predictAll(model_name = 'enet', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 50)

plot(unlist(enet_mut1_reg[[2]]), unlist(enet_mut1_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.change
enet_mut2_reg <- predictAll(model_name = 'enet', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), iterations = 50)

plot(unlist(enet_mut2_reg[[2]]), unlist(enet_mut2_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add gdna.exon.intron
enet_mut3_reg <- predictAll(model_name = 'enet', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron"), iterations = 50)

plot(unlist(enet_mut3_reg[[2]]), unlist(enet_mut3_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
enet_mut4_reg <- predictAll(model_name = 'enet', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro"), iterations = 50)

plot(unlist(enet_mut4_reg[[2]]), unlist(enet_mut4_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add splice.delins.snv
enet_mut5_reg <- predictAll(model_name = 'enet', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), iterations = 50)

plot(unlist(enet_mut5_reg[[2]]), unlist(enet_mut5_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add protein.codon.num
enet_mut6_reg <- predictAll(model_name = 'enet', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), iterations = 50)

plot(unlist(enet_mut6_reg[[2]]), unlist(enet_mut6_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add mdm2.nG
enet_mut7_reg <- predictAll(model_name = 'enet', 
                          data = full_data_rf,
                          subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), iterations = 50)

plot(unlist(enet_mut7_reg[[2]]), unlist(enet_mut7_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)



# combine all tests 

enet_all <- rbind (
  append('just_methyl', c(mean(unlist(enet_methyl[[1]])), enet_methyl[[6]], 'all_data')),
  append('gender_and_gdna.base.change', c(mean(unlist(enet_mut[[1]])), enet_mut[[6]], 'all_data')),
  append('gdna.codon', c(mean(unlist(enet_mut1[[1]])), enet_mut1[[6]], 'all_data')),
  append('protein.codon.change', c(mean(unlist(enet_mut2[[1]])), enet_mut2[[6]], 'all_data')),
  append('gdna.exon.intron', c(mean(unlist(enet_mut3[[1]])), enet_mut3[[6]], 'all_data')),
  append('codon72.npro', c(mean(unlist(enet_mut4[[1]])), enet_mut4[[6]], 'all_data')),
  append('splice.delins.snv', c(mean(unlist(enet_mut5[[1]])), enet_mut5[[6]], 'all_data')),
  append('protein.codon.num', c(mean(unlist(enet_mut6[[1]])), enet_mut6[[6]], 'all_data')),
  append('mdm2.nG', c(mean(unlist(enet_mut7[[1]])), enet_mut7[[6]], 'all_data')),
  
  append('just_methyl', c(mean(unlist(enet_methyl_cor[[1]])), enet_methyl_cor[[6]], 'corr')),
  append('gender_and_gdna.base.change', c(mean(unlist(enet_mut_cor[[1]])), enet_mut_cor[[6]], 'corr')),
  append('gdna.codon', c(mean(unlist(enet_mut1_cor[[1]])), enet_mut1_cor[[6]], 'corr')),
  append('protein.codon.change', c(mean(unlist(enet_mut2_cor[[1]])), enet_mut2_cor[[6]], 'corr')),
  append('gdna.exon.intron', c(mean(unlist(enet_mut3_cor[[1]])), enet_mut3_cor[[6]], 'corr')),
  append('codon72.npro', c(mean(unlist(enet_mut4_cor[[1]])), enet_mut4_cor[[6]], 'corr')),
  append('splice.delins.snv', c(mean(unlist(enet_mut5_cor[[1]])), enet_mut5_cor[[6]], 'corr')),
  append('protein.codon.num', c(mean(unlist(enet_mut6_cor[[1]])), enet_mut6_cor[[6]], 'corr')),
  append('mdm2.nG', c(mean(unlist(enet_mut7_cor[[1]])), enet_mut7_cor[[6]], 'corr')),
  
  append('just_methyl', c(mean(unlist(enet_methyl_reg[[1]])), enet_methyl_reg[[6]], 'recursive')),
  append('gender_and_gdna.base.change', c(mean(unlist(enet_mut_reg[[1]])), enet_mut_reg[[6]], 'recursive')),
  append('gdna.codon', c(mean(unlist(enet_mut1_reg[[1]])), enet_mut1_reg[[6]], 'recursive')),
  append('protein.codon.change', c(mean(unlist(enet_mut2_reg[[1]])), enet_mut2_reg[[6]], 'recursive')),
  append('gdna.exon.intron', c(mean(unlist(enet_mut3_reg[[1]])), enet_mut3_reg[[6]], 'recursive')),
  append('codon72.npro', c(mean(unlist(enet_mut4_reg[[1]])), enet_mut4_reg[[6]], 'recursive')),
  append('splice.delins.snv', c(mean(unlist(enet_mut5_reg[[1]])), enet_mut5_reg[[6]], 'recursive')),
  append('protein.codon.num', c(mean(unlist(enet_mut6_reg[[1]])), enet_mut6_reg[[6]], 'recursive')),
  append('mdm2.nG', c(mean(unlist(enet_mut7_reg[[1]])), enet_mut7_reg[[6]], 'recursive'))
)



##################################################################################################################3
# Lasso
#############################
# Full Data

# just methylation
lasso_methyl <- predictAll(model_name = 'lasso', 
                          data = full_data,
                          subset <- c("age_diagnosis"), 
                          selected_features = NULL, iterations = 50)

plot(unlist(lasso_methyl[[2]]), unlist(lasso_methyl[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# gender and gdna.base.change
lasso_mut <- predictAll(model_name = 'lasso', 
                       data = full_data,
                       subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                       selected_features = c("gender", "gdna.base.change"), iterations = 50)

plot(unlist(lasso_mut[[2]]), unlist(lasso_mut[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# add gdna.codon
lasso_mut1 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 50)

plot(unlist(lasso_mut1[[2]]), unlist(lasso_mut1[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.change
lasso_mut2 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), iterations = 50)

plot(unlist(lasso_mut2[[2]]), unlist(lasso_mut2[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add gdna.exon.intron
lasso_mut3 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                    "gdna.exon.intron"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron"), iterations = 50)

plot(unlist(lasso_mut3[[2]]), unlist(lasso_mut3[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
lasso_mut4 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                    "gdna.exon.intron", "codon72.npro"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron", "codon72.npro"), iterations = 50)

plot(unlist(lasso_mut4[[2]]), unlist(lasso_mut4[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add splice.delins.snv
lasso_mut5 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                    "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), iterations = 50)

plot(unlist(lasso_mut5[[2]]), unlist(lasso_mut5[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add protein.codon.num
lasso_mut6 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                    "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), iterations = 50)

plot(unlist(lasso_mut6[[2]]), unlist(lasso_mut6[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add mdm2.nG
lasso_mut7 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                    "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), iterations = 50)

plot(unlist(lasso_mut7[[2]]), unlist(lasso_mut7[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


#############################
# Corelation data

# just methylation
lasso_methyl_cor <- predictAll(model_name = 'lasso', 
                              data = full_data_cor,
                              subset <- c("age_diagnosis"), 
                              selected_features = NULL, iterations = 50)

plot(unlist(lasso_methyl_cor[[2]]), unlist(lasso_methyl_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# gender and gdna.base.change
lasso_mut_cor <- predictAll(model_name = 'lasso', 
                           data = full_data_cor,
                           subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                           selected_features = c("gender", "gdna.base.change"), iterations = 50)

plot(unlist(lasso_mut_cor[[2]]), unlist(lasso_mut_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# add gdna.codon
lasso_mut1_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 50)

plot(unlist(lasso_mut1_cor[[2]]), unlist(lasso_mut1_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.change
lasso_mut2_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), iterations = 50)

plot(unlist(lasso_mut2_cor[[2]]), unlist(lasso_mut2_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add gdna.exon.intron
lasso_mut3_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron"), iterations = 50)

plot(unlist(lasso_mut3_cor[[2]]), unlist(lasso_mut3_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
lasso_mut4_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro"), iterations = 50)

plot(unlist(lasso_mut4_cor[[2]]), unlist(lasso_mut4_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add splice.delins.snv
lasso_mut5_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), iterations = 50)

plot(unlist(lasso_mut5_cor[[2]]), unlist(lasso_mut5_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add protein.codon.num
lasso_mut6_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), iterations = 50)

plot(unlist(lasso_mut6_cor[[2]]), unlist(lasso_mut6_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add mdm2.nG
lasso_mut7_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), iterations = 50)

plot(unlist(lasso_mut7_cor[[2]]), unlist(lasso_mut7_cor[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)



#############################
# recursive elimination data

# just methylation
lasso_methyl_reg <- predictAll(model_name = 'lasso', 
                              data = full_data_rf,
                              subset <- c("age_diagnosis"), 
                              selected_features = NULL, iterations = 50)

plot(unlist(lasso_methyl_reg[[2]]), unlist(lasso_methyl_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# gender and gdna.base.change
lasso_mut_reg <- predictAll(model_name = 'lasso', 
                           data = full_data_rf,
                           subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                           selected_features = c("gender", "gdna.base.change"), iterations = 50)

plot(unlist(lasso_mut_reg[[2]]), unlist(lasso_mut_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# add gdna.codon
lasso_mut1_reg <- predictAll(model_name = 'lasso', 
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 50)

plot(unlist(lasso_mut1_reg[[2]]), unlist(lasso_mut1_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.change
lasso_mut2_reg <- predictAll(model_name = 'lasso', 
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), iterations = 50)

plot(unlist(lasso_mut2_reg[[2]]), unlist(lasso_mut2_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add gdna.exon.intron
lasso_mut3_reg <- predictAll(model_name = 'lasso', 
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron"), iterations = 50)

plot(unlist(lasso_mut3_reg[[2]]), unlist(lasso_mut3_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
lasso_mut4_reg <- predictAll(model_name = 'lasso', 
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro"), iterations = 50)

plot(unlist(lasso_mut4_reg[[2]]), unlist(lasso_mut4_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add splice.delins.snv
lasso_mut5_reg <- predictAll(model_name = 'lasso', 
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), iterations = 50)

plot(unlist(lasso_mut5_reg[[2]]), unlist(lasso_mut5_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add protein.codon.num
lasso_mut6_reg <- predictAll(model_name = 'lasso', 
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), iterations = 50)

plot(unlist(lasso_mut6_reg[[2]]), unlist(lasso_mut6_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add mdm2.nG
lasso_mut7_reg <- predictAll(model_name = 'lasso', 
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), iterations = 50)

plot(unlist(lasso_mut7_reg[[2]]), unlist(lasso_mut7_reg[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)



# combine all tests 

lasso_all <- rbind (
  append('just_methyl', c(mean(unlist(lasso_methyl[[1]])), lasso_methyl[[6]], 'all_data')),
  append('gender_and_gdna.base.change', c(mean(unlist(lasso_mut[[1]])), lasso_mut[[6]], 'all_data')),
  append('gdna.codon', c(mean(unlist(lasso_mut1[[1]])), lasso_mut1[[6]], 'all_data')),
  append('protein.codon.change', c(mean(unlist(lasso_mut2[[1]])), lasso_mut2[[6]], 'all_data')),
  append('gdna.exon.intron', c(mean(unlist(lasso_mut3[[1]])), lasso_mut3[[6]], 'all_data')),
  append('codon72.npro', c(mean(unlist(lasso_mut4[[1]])), lasso_mut4[[6]], 'all_data')),
  append('splice.delins.snv', c(mean(unlist(lasso_mut5[[1]])), lasso_mut5[[6]], 'all_data')),
  append('protein.codon.num', c(mean(unlist(lasso_mut6[[1]])), lasso_mut6[[6]], 'all_data')),
  append('mdm2.nG', c(mean(unlist(lasso_mut7[[1]])), lasso_mut7[[6]], 'all_data')),
  
  append('just_methyl', c(mean(unlist(lasso_methyl_cor[[1]])), lasso_methyl_cor[[6]], 'corr')),
  append('gender_and_gdna.base.change', c(mean(unlist(lasso_mut_cor[[1]])), lasso_mut_cor[[6]], 'corr')),
  append('gdna.codon', c(mean(unlist(lasso_mut1_cor[[1]])), lasso_mut1_cor[[6]], 'corr')),
  append('protein.codon.change', c(mean(unlist(lasso_mut2_cor[[1]])), lasso_mut2_cor[[6]], 'corr')),
  append('gdna.exon.intron', c(mean(unlist(lasso_mut3_cor[[1]])), lasso_mut3_cor[[6]], 'corr')),
  append('codon72.npro', c(mean(unlist(lasso_mut4_cor[[1]])), lasso_mut4_cor[[6]], 'corr')),
  append('splice.delins.snv', c(mean(unlist(lasso_mut5_cor[[1]])), lasso_mut5_cor[[6]], 'corr')),
  append('protein.codon.num', c(mean(unlist(lasso_mut6_cor[[1]])), lasso_mut6_cor[[6]], 'corr')),
  append('mdm2.nG', c(mean(unlist(lasso_mut7_cor[[1]])), lasso_mut7_cor[[6]], 'corr')),
  
  append('just_methyl', c(mean(unlist(lasso_methyl_reg[[1]])), lasso_methyl_reg[[6]], 'recursive')),
  append('gender_and_gdna.base.change', c(mean(unlist(lasso_mut_reg[[1]])), lasso_mut_reg[[6]], 'recursive')),
  append('gdna.codon', c(mean(unlist(lasso_mut1_reg[[1]])), lasso_mut1_reg[[6]], 'recursive')),
  append('protein.codon.change', c(mean(unlist(lasso_mut2_reg[[1]])), lasso_mut2_reg[[6]], 'recursive')),
  append('gdna.exon.intron', c(mean(unlist(lasso_mut3_reg[[1]])), lasso_mut3_reg[[6]], 'recursive')),
  append('codon72.npro', c(mean(unlist(lasso_mut4_reg[[1]])), lasso_mut4_reg[[6]], 'recursive')),
  append('splice.delins.snv', c(mean(unlist(lasso_mut5_reg[[1]])), lasso_mut5_reg[[6]], 'recursive')),
  append('protein.codon.num', c(mean(unlist(lasso_mut6_reg[[1]])), lasso_mut6_reg[[6]], 'recursive')),
  append('mdm2.nG', c(mean(unlist(lasso_mut7_reg[[1]])), lasso_mut7_reg[[6]], 'recursive'))
)

# write.csv(rf_all, paste0(data_folder, '/rf_all.csv'))
# write.csv(enet_all, paste0(data_folder, '/enet_all.csv'))
# write.csv(lasso_all, paste0(data_folder, '/lasso_all.csv'))

##########################################################################################
# combine ran_forest, enet, lasso
ran_forest <- as.data.frame(rf_all)
ran_forest$model <- 'rand_forest'
names(ran_forest) <- c('variables', 'rmse', 'observations', 'data', 'model')

enet <- as.data.frame(enet_all)
enet$model <- 'enet'
names(enet) <- c('variables', 'rmse', 'observations', 'data', 'model')

lasso <- as.data.frame(lasso_all)
lasso$model <- 'lasso'
names(lasso) <- c('variables', 'rmse', 'observations', 'data', 'model')

models <- rbind(ran_forest, enet, lasso)
models$rmse <- as.numeric(as.character(models$rmse))


# write.csv(models, paste0(data_folder, '/clin_methyl_model_results.csv'))






