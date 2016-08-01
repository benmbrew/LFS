### 
# this script will run the models

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/regression_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


# Load Libraries and Dependencies ---------------------------------
library("igraph")
library("foreach")
library("ROCR")
library("doParallel")
library("caret")
library("glmnet")
library("randomForest")
library("kernlab")
library("pROC")
library(dplyr)
library(Metrics)

# Load in clinical data
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

#############################################################################################
# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(model_name, 
                       subset, 
                       selected_features,
                       iterations) {
  
  model <- list()
  predictions <- list()
  mse <- list()
  importance <- list()
  test.ground_truth <- list()
  
  clin <- clin[, subset]
  
  # Try the model with all different selection of features based on number of missinginess. 
  clin <- clin[complete.cases(clin),]
  
  
  # convert characters to factors 
  for ( i in 1:ncol(clin)){
    
  if(!grepl('age', colnames(clin[i]))) {
      clin[,i] <- as.factor(clin[,i])
      
    } 
  }
  
  for (i in 1:iterations){
    
    set.seed(i)
    train_index <- sample(nrow(clin), nrow(clin) *.7, replace = F)
    
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
      
      mtry <- sqrt(ncol(clin))
      tunegrid <- expand.grid(.mtry=mtry)
      
      rf_y = clin$age_diagnosis[train_index]
      
      
      model[[i]] <- train(x = clin[train_index, selected_features]
                          , y = rf_y
                          , method = "rf"
                          , trControl = fitControl
                          , tuneGrid = tunegrid
                          , importance = T
                          , verbose = FALSE)
      
      importance[[i]] <- varImp(model[[i]])
      
      
      predictions[[i]] <- predict(model[[i]] 
                                  , newdata = clin[-train_index, selected_features])
      
      test.ground_truth[[i]] <- clin$age_diagnosis[-train_index]
      mse[[i]] <- rmse(unlist(predictions[[i]]), unlist(test.ground_truth[[i]]))
      
    }
    
    if (model_name == 'enet') {
      
      enet_y = clin$age_diagnosis[train_index]
      N_CV_REPEATS <- 2
      elastic_net.cv_error = vector()
      elastic_net.cv_model = list()
      elastic_net.ALPHA <- c(1:9) / 10 # c
      
      temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
        
        for (alpha_index in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
        {      
          elastic_net.cv_model[[alpha_index]] = cv.glmnet(data.matrix(clin[train_index, selected_features])
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
          data.matrix(clin[train_index, selected_features])
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
      model[[i]] = glmnet(data.matrix(clin[train_index, selected_features]), 
                          enet_y, 
                          alpha = elastic_net.ALPHA[temp.best_alpha_index],
                          standardize=FALSE,
                          nlambda = 100,
                          family = 'gaussian')
      
      # This returns 100 prediction with 1-100 lambdas
      temp_predictions <- predict(model[[i]], data.matrix(clin[-train_index, selected_features]))
      #print(dim(temp.predictions))
      predictions[[i]] <- temp_predictions[ , temp.min_lambda_index]  # this grabs the opitmal lambda 
      temp.l <- min(length(elastic_net.cv_model$lambda), length(model[[i]]$lambda)) 
      stopifnot(elastic_net.cv_model$lambda[1:temp.l] == model[[i]]$lambda[1:temp.l])  
      test.ground_truth[[i]] <- clin$age_diagnosis[-train_index]
      mse[[i]] <- rmse(unlist(predictions[[i]]), unlist(test.ground_truth[[i]]))
      
    }
    
    if(model_name == "lasso"){
      
      # 2) Lasso Logistic Regression 
      lasso_y <- clin$age_diagnosis[train_index]
      temp.non_zero_coeff = 0
      temp.loop_count = 0
      while (temp.non_zero_coeff < 3) {     
          temp.cv_model = cv.glmnet(data.matrix(clin[train_index, selected_features])
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
      
      model[[i]] = glmnet(data.matrix(clin[train_index, selected_features]), 
                           lasso_y, 
                           alpha = 1,
                           standardize=FALSE,
                           nlambda = 100,
                           family = 'gaussian')
      
    
      # This returns 100 prediction with 1-100 lambdas
      temp_predictions <- predict(model[[i]], data.matrix(clin[-train_index, selected_features]))
      #print(dim(temp.predictions))
      predictions[[i]] <- temp_predictions[ , temp.min_lambda_index]  # this grabs the opitmal lambda 
     
      test.ground_truth[[i]] <- clin$age_diagnosis[-train_index]
      mse[[i]] <- rmse(unlist(predictions[[i]]), unlist(test.ground_truth[[i]]))

    }
    
    if(model_name == "svm") {
      
      
      svm_y = clin$age_diagnosis[train_index]
      
      model[[i]] <- train(x = clin[train_index, selected_features]
                          , y = svm_y
                          , method = "svmLinear"
                          , trControl = fitControl  
                          , importance = T
                          , verbose = FALSE
                          
      )
      
      predictions[[i]] <- predict(model[[i]] 
                                  , newdata = clin[-train_index, selected_features])
      
      importance[[i]] <- varImp(model[[i]])
      
      test.ground_truth[[i]] <- clin$age_diagnosis[-train_index]
      mse[[i]] <- rmse(unlist(predictions[[i]]), unlist(test.ground_truth[[i]]))
      
    }
    
    print(i)
    
  }
  
  return(list(mse, predictions, model, importance, test.ground_truth))
  
}


########################################################################################################gdna.exon.intron

# variables missing
# gdna.exon.intron 492
# gdna.base.change 164
# gdna.codon 164
# protein.codon.change 177
# protein.codon.num 549
# splice.delins.snv 519
# codon72.npro 517
# mdm2.nG 652
##################################################################################################################3
# Random forest
# all variables
rf_mut <- predictAll(model_name = 'rf', 
                     subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                 "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                 "mdm2.nG"), 
                     selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                           "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                           "mdm2.nG"), iterations = 50)

plot(unlist(rf_mut[[2]]), unlist(rf_mut[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# remove mdm2.nG
rf_mut1 <- predictAll(model_name = 'rf', 
                     subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                 "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro"), 
                     selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                           "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro"), iterations = 50)

plot(unlist(rf_mut1[[2]]), unlist(rf_mut1[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# remove gdna.exon.intron
rf_mut2 <- predictAll(model_name = 'rf', 
                     subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.base.change",
                                 "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                 "mdm2.nG"), 
                     selected_features = c("gender", "gdna.base.change",
                                           "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                           "mdm2.nG"), iterations = 50)

plot(unlist(rf_mut2[[2]]), unlist(rf_mut2[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# remove protein.codon.num
rf_mut3 <- predictAll(model_name = 'rf', 
                     subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                 "gdna.codon", "protein.codon.change", "splice.delins.snv", "codon72.npro",
                                 "mdm2.nG"), 
                     selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                           "gdna.codon", "protein.codon.change", "splice.delins.snv", "codon72.npro",
                                           "mdm2.nG"), iterations = 50)

plot(unlist(rf_mut3[[2]]), unlist(rf_mut3[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# remove splice.delins.snv
rf_mut4 <- predictAll(model_name = 'rf', 
                     subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                 "gdna.codon", "protein.codon.change", "protein.codon.num", "codon72.npro",
                                 "mdm2.nG"), 
                     selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                           "gdna.codon", "protein.codon.change", "protein.codon.num", "codon72.npro",
                                           "mdm2.nG"), iterations = 50)

plot(unlist(rf_mut4[[2]]), unlist(rf_mut4[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# remove gdna.base.change
rf_mut5 <- predictAll(model_name = 'rf', 
                     subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron",
                                 "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                 "mdm2.nG"), 
                     selected_features = c("gender", "gdna.exon.intron",
                                           "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                           "mdm2.nG"), iterations = 50)

plot(unlist(rf_mut5[[2]]), unlist(rf_mut5[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# remove gdna.codon
rf_mut6 <- predictAll(model_name = 'rf', 
                     subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change", 
                                 "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                 "mdm2.nG"), 
                     selected_features = c("gender", "gdna.exon.intron", "gdna.base.change", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                           "mdm2.nG"), iterations = 50)

plot(unlist(rf_mut6[[2]]), unlist(rf_mut6[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# remove protein.codon.change 
rf_mut7 <- predictAll(model_name = 'rf', 
                     subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                 "gdna.codon", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                 "mdm2.nG"), 
                     selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                           "gdna.codon", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                           "mdm2.nG"), iterations = 50)

plot(unlist(rf_mut7[[2]]), unlist(rf_mut7[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


rf <- rbind (
  append('all_variables', mean(unlist(rf_mut[[1]]))),
  append('no_mdm2.nG', mean(unlist(rf_mut1[[1]]))),
  append('no_gdna.exon.intron', mean(unlist(rf_mut2[[1]]))),
  append('no_protein.codon.num',mean(unlist(rf_mut3[[1]]))),
  append('no_splice.delins.snv', mean(unlist(rf_mut4[[1]]))),
  append('no_gdna.base.change', mean(unlist(rf_mut5[[1]]))),
  append('no_gdna.codon', mean(unlist(rf_mut6[[1]]))),
  append('no_proteinp.codon.change', mean(unlist(rf_mut7[[1]])))
)
         
##################################################################################################################3
# Elastic Net
# all variables
enet_mut <- predictAll(model_name = 'enet', 
                     subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                 "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                 "mdm2.nG"), 
                     selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                           "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                           "mdm2.nG"), iterations = 50)

plot(unlist(enet_mut[[2]]), unlist(enet_mut[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# remove mdm2.nG
enet_mut1 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                  "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro"), 
                      selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                            "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro"), iterations = 50)

plot(unlist(enet_mut1[[2]]), unlist(enet_mut1[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# remove gdna.exon.intron
enet_mut2 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.base.change",
                                  "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                  "mdm2.nG"), 
                      selected_features = c("gender", "gdna.base.change",
                                            "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                            "mdm2.nG"), iterations = 50)

plot(unlist(enet_mut2[[2]]), unlist(enet_mut2[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# remove protein.codon.num
enet_mut3 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                  "gdna.codon", "protein.codon.change", "splice.delins.snv", "codon72.npro",
                                  "mdm2.nG"), 
                      selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                            "gdna.codon", "protein.codon.change", "splice.delins.snv", "codon72.npro",
                                            "mdm2.nG"), iterations = 50)

plot(unlist(enet_mut3[[2]]), unlist(enet_mut3[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# remove splice.delins.snv
enet_mut4 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                  "gdna.codon", "protein.codon.change", "protein.codon.num", "codon72.npro",
                                  "mdm2.nG"), 
                      selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                            "gdna.codon", "protein.codon.change", "protein.codon.num", "codon72.npro",
                                            "mdm2.nG"), iterations = 50)

plot(unlist(enet_mut4[[2]]), unlist(enet_mut4[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# remove gdna.base.change
enet_mut5 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron",
                                  "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                  "mdm2.nG"), 
                      selected_features = c("gender", "gdna.exon.intron",
                                            "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                            "mdm2.nG"), iterations = 50)

plot(unlist(enet_mut5[[2]]), unlist(enet_mut5[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# remove gdna.codon
enet_mut6 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change", 
                                  "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                  "mdm2.nG"), 
                      selected_features = c("gender", "gdna.exon.intron", "gdna.base.change", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                            "mdm2.nG"), iterations = 50)

plot(unlist(enet_mut6[[2]]), unlist(enet_mut6[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# remove protein.codon.change 
enet_mut7 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                  "gdna.codon", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                  "mdm2.nG"), 
                      selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                            "gdna.codon", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                            "mdm2.nG"), iterations = 50)

plot(unlist(enet_mut7[[2]]), unlist(enet_mut7[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


enet <- rbind (
  append('all_variables', mean(unlist(enet_mut[[1]]))),
  append('no_mdm2.nG', mean(unlist(enet_mut1[[1]]))),
  append('no_gdna.exon.intron', mean(unlist(enet_mut2[[1]]))),
  append('no_protein.codon.num',mean(unlist(enet_mut3[[1]]))),
  append('no_splice.delins.snv', mean(unlist(enet_mut4[[1]]))),
  append('no_gdna.base.change', mean(unlist(enet_mut5[[1]]))),
  append('no_gdna.codon', mean(unlist(enet_mut6[[1]]))),
  append('no_proteinp.codon.change', mean(unlist(enet_mut7[[1]])))
)


##################################################################################################################3
# Lasso
# all variables
lasso_mut <- predictAll(model_name = 'lasso', 
                     subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                 "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                 "mdm2.nG"), 
                     selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                           "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                           "mdm2.nG"), iterations = 50)

plot(unlist(lasso_mut[[2]]), unlist(lasso_mut[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# remove mdm2.nG
lasso_mut1 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                  "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro"), 
                      selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                            "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro"), iterations = 50)

plot(unlist(lasso_mut1[[2]]), unlist(lasso_mut1[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# remove gdna.exon.intron
lasso_mut2 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.base.change",
                                  "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                  "mdm2.nG"), 
                      selected_features = c("gender", "gdna.base.change",
                                            "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                            "mdm2.nG"), iterations = 50)

plot(unlist(lasso_mut2[[2]]), unlist(lasso_mut2[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# remove protein.codon.num
lasso_mut3 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                  "gdna.codon", "protein.codon.change", "splice.delins.snv", "codon72.npro",
                                  "mdm2.nG"), 
                      selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                            "gdna.codon", "protein.codon.change", "splice.delins.snv", "codon72.npro",
                                            "mdm2.nG"), iterations = 50)

plot(unlist(lasso_mut3[[2]]), unlist(lasso_mut3[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# remove splice.delins.snv
lasso_mut4 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                  "gdna.codon", "protein.codon.change", "protein.codon.num", "codon72.npro",
                                  "mdm2.nG"), 
                      selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                            "gdna.codon", "protein.codon.change", "protein.codon.num", "codon72.npro",
                                            "mdm2.nG"), iterations = 50)

plot(unlist(lasso_mut4[[2]]), unlist(lasso_mut4[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# remove gdna.base.change
lasso_mut5 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron",
                                  "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                  "mdm2.nG"), 
                      selected_features = c("gender", "gdna.exon.intron",
                                            "gdna.codon", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                            "mdm2.nG"), iterations = 50)

plot(unlist(lasso_mut5[[2]]), unlist(lasso_mut5[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# remove gdna.codon
lasso_mut6 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change", 
                                  "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                  "mdm2.nG"), 
                      selected_features = c("gender", "gdna.exon.intron", "gdna.base.change", "protein.codon.change", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                            "mdm2.nG"), iterations = 50)

plot(unlist(lasso_mut6[[2]]), unlist(lasso_mut6[[5]]), xlab = 'Predictions', ylab = 'Actual Age')

# remove protein.codon.change 
lasso_mut7 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", 'methyl_indicator', "gdna.exon.intron", "gdna.base.change",
                                  "gdna.codon", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                  "mdm2.nG"), 
                      selected_features = c("gender", "gdna.exon.intron", "gdna.base.change",
                                            "gdna.codon", "protein.codon.num", "splice.delins.snv", "codon72.npro",
                                            "mdm2.nG"), iterations = 50)

plot(unlist(lasso_mut7[[2]]), unlist(lasso_mut7[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


lasso <- rbind (
  append('all_variables', mean(unlist(lasso_mut[[1]]))),
  append('no_mdm2.nG', mean(unlist(lasso_mut1[[1]]))),
  append('no_gdna.exon.intron', mean(unlist(lasso_mut2[[1]]))),
  append('no_protein.codon.num',mean(unlist(lasso_mut3[[1]]))),
  append('no_splice.delins.snv', mean(unlist(lasso_mut4[[1]]))),
  append('no_gdna.base.change', mean(unlist(lasso_mut5[[1]]))),
  append('no_gdna.codon', mean(unlist(lasso_mut6[[1]]))),
  append('no_proteinp.codon.change', mean(unlist(lasso_mut7[[1]])))
)

##########################################################################################
# combine ran_forest, enet, lasso
ran_forest <- as.data.frame(rf)
ran_forest$model <- 'rand_forest'
names(ran_forest) <- c('variables', 'mse', 'model')

enet <- as.data.frame(enet)
enet$model <- 'enet'
names(enet) <- c('variables', 'mse', 'model')

lasso <- as.data.frame(lasso)
lasso$model <- 'lasso'
names(lasso) <- c('variables', 'mse', 'model')

models <- rbind(ran_forest, enet, lasso)
models$mse <- as.numeric(as.character(models$mse))

ggplot(models, aes(model, mse, group = variables, fill = variables)) + 
  geom_bar(stat = 'identity', 
           position = 'dodge')





