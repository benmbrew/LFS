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
library(doParallel)

registerDoParallel(1)

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
                       log,
                       iterations) {
  
  model <- list()
  importance <- list()
  train.mse <- list()
  test.mse <- list()
  train.predictions <- list()
  test.predictions <- list()
  train.ground_truth <- list()
  test.ground_truth <- list()
  train.sample_collection <- list()
  test.sample_collection <- list()

  
  # set log transformation
  if(log) {
    
    data[,c(6,8, 27:ncol(data))]  <- log(data[,c(6,8,27:ncol(data))])
  }
  genes <- colnames(data)[27:ncol(data)]
  
  data <- data[, c(subset, genes)]
  
  # remove rows where age of diagnosis meissin
  data <- data[!(is.na(data$age_diagnos)),]
  
  for (i in 3:ncol(data)) {
    if(any(is.na(data[, i]))) {
      data <- data[!is.na(data[,i]),]
    }
  }
  
  
  obs <- nrow(data)
  
  # convert characters to factors
  for ( i in 3:ncol(data)){

    if(typeof(data[,i]) == 'character' || grepl('int', names(data[i]))) {
      data[,i] <- as.factor(data[,i])
      print(i)
    }
  }
  
  for (i in 1:iterations){
    
    set.seed(i)
    train_index <- sample(nrow(data), nrow(data) *cutoff, replace = F)
    
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
      
      temp <- varImp(model[[i]])[[1]]
      importance[[i]] <- cbind(rownames(temp), temp$Overall)
      
      
      test.predictions[[i]] <- predict(model[[i]] 
                                  , newdata = data[-train_index, c(selected_features, genes)])
      
      train.predictions[[i]] <- predict(model[[i]] 
                                  , newdata = data[train_index, c(selected_features, genes)])
      
      train.ground_truth[[i]] <- data$age_diagnosis[train_index]
      test.ground_truth[[i]] <- data$age_diagnosis[-train_index]
      train.sample_collection[[i]] = data$age_sample_collection[train_index]
      test.sample_collection[[i]] = data$age_sample_collection[-train_index]
      train.mse[[i]] <- rmse(unlist(train.predictions[[i]]), unlist(train.ground_truth[[i]]))
      test.mse[[i]] <- rmse(unlist(test.predictions[[i]]), unlist(test.ground_truth[[i]]))
      
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
      real_y[[i]] = data$age_sample_collection[-train_index]
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
      real_y[[i]] = data$age_sample_collection[-train_index]
      mse[[i]] <- rmse(unlist(predictions[[i]]), unlist(test.ground_truth[[i]]))
      
    }
    
    if(model_name == "lm") {
      
      
      lm_y = data$age_diagnosis[train_index]
      
      model[[i]] <- train(x = data.matrix(data[train_index, c(selected_features, genes)])
                          , y = lm_y
                          , method = "lm"
                          , trControl = fitControl  

      )
      
      predictions[[i]] <- predict(model[[i]] 
                                  , newdata = data[-train_index, selected_features])
      
      importance[[i]] <- varImp(model[[i]])
      
      test.ground_truth[[i]] <- data$age_diagnosis[-train_index]
      mse[[i]] <- rmse(unlist(predictions[[i]]), unlist(test.ground_truth[[i]]))
      
    }
    
    print(i)
    
  }
  
  return(list(train.mse, test.mse,  train.predictions, test.predictions, train.ground_truth, test.ground_truth, train.sample_collection,
              test.sample_collection, model, importance, obs))
  
}


########################################################################################################gdna.exon.intron
# save.image(paste0(data_folder, '/clin_methyl_models.RData'))
# load(paste0(data_folder, '/clin_methyl_models.RData'))

# histogram of age of diagnosis
hist(full_data_rf$age_diagnosis, xlab = 'Age of Diagnosis', main = 'Distribution of Age of Diagnosis', col = 'lightblue')
hist(full_data_rf$age_sample_collection, xlab = 'Age of Sample Collection', col = 'lightblue', main = 'Distribution of Sample Collection')


# plot of age of onset vs age of diagnosis with r squared
plot(full_data$age_diagnosis, full_data$age_sample_collection, xlab = 'Age of Diagnosis',
     ylab = 'Age of Sample Collection')
abline(0,1)
r_squared <- round(summary(lm(full_data$age_diagnosis ~ full_data$age_sample_collection))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))


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
pdf('/home/benbrew/Desktop/rf_methylation.pdf')

# Random forest
#############################
# Full Data
# 
# # just methylation
# rf_methyl <- predictAll(model_name = 'rf',
#                      data = full_data,
#                      subset <- c("age_diagnosis", "age_sample_collection"),
#                      selected_features = NULL, iterations = 10)
# 
# 
# plot(unlist(rf_methyl[[2]]), unlist(rf_methyl[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = 'Just methylation')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_methyl[[2]]) ~ unlist(rf_methyl[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_methyl[[6]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# 
# 
# plot(unlist(rf_methyl[[2]]), unlist(rf_methyl[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = 'Just methylation')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_methyl[[6]]))
# 
# # gender and gdna.base.change
# rf_mut <- predictAll(model_name = 'rf', 
#                      data = full_data,
#                      subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change"), 
#                      selected_features = c("gender", "gdna.base.change"), iterations = 10)
# 
# plot(unlist(rf_mut[[2]]), unlist(rf_mut[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut[[6]]))
# 
# plot(unlist(rf_mut[[2]]), unlist(rf_mut[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut[[6]]))
# 
# 
# 
# 
# # add gdna.codon
# rf_mut1 <- predictAll(model_name = 'rf', 
#                       data = full_data,
#                       subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.base.change", "gdna.codon"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 10)
# 
# plot(unlist(rf_mut1[[2]]), unlist(rf_mut1[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut1[[6]]))
# 
# plot(unlist(rf_mut1[[2]]), unlist(rf_mut1[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut1[[6]]))
# 
# 
# # add protein.codon.change
# rf_mut2 <- predictAll(model_name = 'rf', 
#                       data = full_data,
#                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
#                       iterations = 10)
# 
# plot(unlist(rf_mut2[[2]]), unlist(rf_mut2[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut2[[6]]))
# 
# plot(unlist(rf_mut2[[2]]), unlist(rf_mut2[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut2[[6]]))
# 
# 
# 
# # add gdna.exon.intron
# rf_mut3 <- predictAll(model_name = 'rf', 
#                       data = full_data,
#                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                   "gdna.exon.intron"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                             "gdna.exon.intron"), iterations = 10)
# 
# plot(unlist(rf_mut3[[2]]), unlist(rf_mut3[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut3[[6]]))
# 
# plot(unlist(rf_mut3[[2]]), unlist(rf_mut3[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut3[[6]]))
# 
# # add codon72.npro
# rf_mut4 <- predictAll(model_name = 'rf', 
#                       data = full_data,
#                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                   "gdna.exon.intron", "codon72.npro"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                             "gdna.exon.intron", "codon72.npro"), iterations = 10)
# 
# plot(unlist(rf_mut4[[2]]), unlist(rf_mut4[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut4[[6]]))
# 
# plot(unlist(rf_mut4[[2]]), unlist(rf_mut4[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut4[[6]]))
# 
# 
# 
# # add splice.delins.snv
# rf_mut5 <- predictAll(model_name = 'rf', 
#                       data = full_data,
#                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                   "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                             "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
#                       iterations = 10)
# 
# plot(unlist(rf_mut5[[2]]), unlist(rf_mut5[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut5[[6]]))
# 
# plot(unlist(rf_mut5[[2]]), unlist(rf_mut5[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut5[[6]]))
# 
# 
# # add protein.codon.num
# rf_mut6 <- predictAll(model_name = 'rf', 
#                       data = full_data,
#                       subset <- c("age_diagnosis","age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                   "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                             "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
#                       iterations = 10)
# 
# plot(unlist(rf_mut6[[2]]), unlist(rf_mut6[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut6[[6]]))
# 
# plot(unlist(rf_mut6[[2]]), unlist(rf_mut6[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut6[[6]]))
# 
# 
# # add mdm2.nG
# rf_mut7 <- predictAll(model_name = 'rf', 
#                       data = full_data,
#                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                   "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                             "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                       iterations = 10)
# 
# plot(unlist(rf_mut7[[2]]), unlist(rf_mut7[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
#      mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut7[[6]]))
# 
# plot(unlist(rf_mut7[[2]]), unlist(rf_mut7[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
#      mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut7[[6]]))



####################################################################
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
# Run with #gdna.codon, protein.codon.num, mdm2.nG, protein.codon.change
# return(list(train.mse, test.mse,  train.predictions, test.predictions, train.ground_truth, 
#             test.ground_truth, train.sample_collection,
#             test.sample_collection, model, importance, obs))
##################################################################################################
# just methylation
rf_methyl <- predictAll(model_name = 'rf',
                        data = full_data_rf,
                        subset <- c("age_diagnosis", "age_sample_collection"),
                        selected_features = NULL, 
                        log = FALSE,
                        iterations = 10)

rf_methyl[[10]]
top <- as.data.frame(do.call('rbind', rf_methyl[[10]]))
top$V2 <- as.numeric(as.character(top$V2))

mean_top <- top %>%
  group_by(V1) %>%
  summarise(mean_score = mean(V2, na.rm = T))

# rank them 
mean_top <- mean_top[order(mean_top$mean_score, decreasing = T), ]
###############################
# # plot train age diagnosis
# plot(unlist(rf_methyl[[3]]), unlist(rf_methyl[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just methylation no log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_methyl[[3]]) ~ unlist(rf_methyl[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_methyl[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_methyl[[4]]), unlist(rf_methyl[[6]]), 
     xlab = 'Test Predictions', 
     ylab = 'Test Age of Diagnosis',
     xlim= c(0, 900),
     ylim= c(0, 900),
     main = 'Just methylation')
abline(lm(unlist(rf_methyl[[4]]) ~ unlist(rf_methyl[[6]])))
r_squared <- round(summary(lm(unlist(rf_methyl[[4]]) ~ unlist(rf_methyl[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_methyl[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_methyl[[3]]), unlist(rf_methyl[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just methylation')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_methyl[[3]]) ~ unlist(rf_methyl[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_methyl[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_methyl[[4]]), unlist(rf_methyl[[8]]), 
     xlab = 'Test Predictions', 
     ylab = 'Test Age Sample Collection',
     xlim= c(0, 900),
     ylim= c(0, 900),
     main = 'Just methylation')
abline(lm(unlist(rf_methyl[[4]]) ~ unlist(rf_methyl[[8]])))
r_squared <- round(summary(lm(unlist(rf_methyl[[4]]) ~ unlist(rf_methyl[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_methyl[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))



# just methylation with log
rf_methyl_log <- predictAll(model_name = 'rf',
                        data = full_data_rf,
                        subset <- c("age_diagnosis", "age_sample_collection"),
                        selected_features = NULL, 
                        log = TRUE,
                        iterations = 10)


###############################
# # plot train age diagnosis
# plot(unlist(rf_methyl_log[[3]]), unlist(rf_methyl_log[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just methylation with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_methyl_log[[3]]) ~ unlist(rf_methyl_log[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_methyl_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_methyl_log[[4]]), unlist(rf_methyl_log[[6]]), 
     xlab = 'Test Predictions', 
     ylab = 'Test Age of Diagnosis',
     xlim = c(0, 8),
     ylim = c(0, 8),
     main = 'Just methylation with log')
abline(lm(unlist(rf_methyl_log[[4]]) ~ unlist(rf_methyl_log[[6]])))
r_squared <- round(summary(lm(unlist(rf_methyl_log[[4]]) ~ unlist(rf_methyl_log[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_methyl_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_methyl_log[[3]]), unlist(rf_methyl_log[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just methylation with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_methyl_log[[3]]) ~ unlist(rf_methyl_log[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_methyl_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_methyl_log[[4]]), unlist(rf_methyl_log[[8]]), 
     xlab = 'Test Predictions', 
     ylab = 'Test Age Sample Collection',
     xlim = c(0, 8),
     ylim = c(0, 8),
     
     main = 'Just methylation with log')
abline(lm(unlist(rf_methyl_log[[4]]) ~ unlist(rf_methyl_log[[8]])))
r_squared <- round(summary(lm(unlist(rf_methyl_log[[4]]) ~ unlist(rf_methyl_log[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_methyl_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


##################################################################################################
# just gdna.base.change
rf_gdna.base.change <- predictAll(model_name = 'rf',
                        data = full_data_rf,
                        subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change"),
                        selected_features = c("gender", "gdna.base.change"),
                        log = FALSE,
                        iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_gdna.base.change[[3]]), unlist(rf_gdna.base.change[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just gdna.base.change no log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_gdna.base.change[[3]]) ~ unlist(rf_gdna.base.change[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.base.change[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_gdna.base.change[[4]]), unlist(rf_gdna.base.change[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+gdna.base.change')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_gdna.base.change[[4]]) ~ unlist(rf_gdna.base.change[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.base.change[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_gdna.base.change[[3]]), unlist(rf_gdna.base.change[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just gdna.base.change')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_gdna.base.change[[3]]) ~ unlist(rf_gdna.base.change[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.base.change[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_gdna.base.change[[4]]), unlist(rf_gdna.base.change[[8]]), xlab = 'Test Predictions', ylab = 'Test Age Sample Collection',
     main = 'gender+gdna.base.change ')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_gdna.base.change[[4]]) ~ unlist(rf_gdna.base.change[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.base.change[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))



# just methylation with log
rf_gdna.base.change_log <- predictAll(model_name = 'rf',
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change"),
                            selected_features = c("gender", "gdna.base.change"), 
                            log = TRUE,
                            iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_gdna.base.change_log[[3]]), unlist(rf_gdna.base.change_log[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just gdna.base.change with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_gdna.base.change_log[[3]]) ~ unlist(rf_gdna.base.change_log[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.base.change_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_gdna.base.change_log[[4]]), unlist(rf_gdna.base.change_log[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+gdna.base.change with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_gdna.base.change_log[[4]]) ~ unlist(rf_gdna.base.change_log[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.base.change_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


# ###############################
# # plot train age sample collection
# plot(unlist(rf_gdna.base.change_log[[3]]), unlist(rf_gdna.base.change_log[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just gdna.base.change with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_gdna.base.change_log[[3]]) ~ unlist(rf_gdna.base.change_log[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.base.change_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_gdna.base.change_log[[4]]), unlist(rf_gdna.base.change_log[[8]]), xlab = 'Test Predictions', ylab = 'Test Age Sample Collection',
     main = 'gender+gdna.base.change with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_gdna.base.change_log[[4]]) ~ unlist(rf_gdna.base.change_log[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.base.change_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


##################################################################################################
# just gdna.codon
rf_gdna.codon <- predictAll(model_name = 'rf',
                                  data = full_data_rf,
                                  subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon"),
                                  selected_features = c("gender", "gdna.codon"),
                                  log = FALSE,
                                  iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_gdna.codon[[3]]), unlist(rf_gdna.codon[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just gdna.codon no log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_gdna.codon[[3]]) ~ unlist(rf_gdna.codon[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_gdna.codon[[4]]), unlist(rf_gdna.codon[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+gdna.codon')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_gdna.codon[[4]]) ~ unlist(rf_gdna.codon[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_gdna.codon[[3]]), unlist(rf_gdna.codon[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just gdna.codon')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_gdna.codon[[3]]) ~ unlist(rf_gdna.codon[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_gdna.codon[[4]]), unlist(rf_gdna.codon[[8]]), xlab = 'Test Predictions', ylab = 'Test Age Sample Collection',
     main = 'Just gdna.codon')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_gdna.codon[[4]]) ~ unlist(rf_gdna.codon[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))



# just methylation with log
rf_gdna.codon_log <- predictAll(model_name = 'rf',
                                      data = full_data_rf,
                                      subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon"),
                                      selected_features = c("gender", "gdna.codon"),
                                      log = TRUE,
                                      iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_gdna.codon_log[[3]]), unlist(rf_gdna.codon_log[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just gdna.codon with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_gdna.codon_log[[3]]) ~ unlist(rf_gdna.codon_log[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_gdna.codon_log[[4]]), unlist(rf_gdna.codon_log[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+gdna.codon with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_gdna.codon_log[[4]]) ~ unlist(rf_gdna.codon_log[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_gdna.codon_log[[3]]), unlist(rf_gdna.codon_log[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just gdna.codon with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_gdna.codon_log[[3]]) ~ unlist(rf_gdna.codon_log[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_gdna.codon_log[[4]]), unlist(rf_gdna.codon_log[[8]]), xlab = 'Test Predictions', ylab = 'Test Age Sample Collection',
     main = 'gender+gdna.codon with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_gdna.codon_log[[4]]) ~ unlist(rf_gdna.codon_log[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))



##################################################################################################
# just protein.codon.change
rf_protein.codon.change <- predictAll(model_name = 'rf',
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "age_sample_collection", "gender", "protein.codon.change"),
                            selected_features = c("gender", "protein.codon.change"),
                            log = FALSE,
                            iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_protein.codon.change[[3]]), unlist(rf_protein.codon.change[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'gender+protein.codon.change no log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_protein.codon.change[[3]]) ~ unlist(rf_protein.codon.change[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.change[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_protein.codon.change[[4]]), unlist(rf_protein.codon.change[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+protein.codon.change')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_protein.codon.change[[4]]) ~ unlist(rf_protein.codon.change[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.change[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_protein.codon.change[[3]]), unlist(rf_protein.codon.change[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'gender+protein.codon.change')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_protein.codon.change[[3]]) ~ unlist(rf_protein.codon.change[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.change[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_protein.codon.change[[4]]), unlist(rf_protein.codon.change[[8]]), xlab = 'Test Predictions', ylab = 'Test Age Sample Collection',
     main = 'gender+protein.codon.change')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_protein.codon.change[[4]]) ~ unlist(rf_protein.codon.change[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.change[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))



# just methylation with log
rf_protein.codon.change_log <- predictAll(model_name = 'rf',
                                data = full_data_rf,
                                subset <- c("age_diagnosis", "age_sample_collection", "gender", "protein.codon.change"),
                                selected_features = c("gender", "protein.codon.change"),
                                log = TRUE,
                                iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_protein.codon.change_log[[3]]), unlist(rf_protein.codon.change_log[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just protein.codon.change with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_protein.codon.change_log[[3]]) ~ unlist(rf_protein.codon.change_log[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.change_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_protein.codon.change_log[[4]]), unlist(rf_protein.codon.change_log[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+protein.codon.change with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_protein.codon.change_log[[4]]) ~ unlist(rf_protein.codon.change_log[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.change_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_protein.codon.change_log[[3]]), unlist(rf_protein.codon.change_log[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just protein.codon.change with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_protein.codon.change_log[[3]]) ~ unlist(rf_protein.codon.change_log[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.change_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_protein.codon.change_log[[4]]), unlist(rf_protein.codon.change_log[[8]]), xlab = 'Test Predictions', ylab = 'Test Age Sample Collection',
     main = 'gender+protein.codon.change with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_protein.codon.change_log[[4]]) ~ unlist(rf_protein.codon.change_log[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.change_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


##################################################################################################
# just gdna.exon.intron
rf_gdna.exon.intron <- predictAll(model_name = 'rf',
                                      data = full_data_rf,
                                      subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.exon.intron"),
                                      selected_features = c("gender", "gdna.exon.intron"),
                                      log = FALSE,
                                      iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_gdna.exon.intron[[3]]), unlist(rf_gdna.exon.intron[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just gdna.exon.intron no log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_gdna.exon.intron[[3]]) ~ unlist(rf_gdna.exon.intron[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.exon.intron[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

#test 
plot(unlist(rf_gdna.exon.intron[[4]]), unlist(rf_gdna.exon.intron[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+gdna.exon.intron')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_gdna.exon.intron[[4]]) ~ unlist(rf_gdna.exon.intron[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.exon.intron[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_gdna.exon.intron[[3]]), unlist(rf_gdna.exon.intron[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just gdna.exon.intron')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_gdna.exon.intron[[3]]) ~ unlist(rf_gdna.exon.intron[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.exon.intron[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_gdna.exon.intron[[4]]), unlist(rf_gdna.exon.intron[[8]]), xlab = 'Test Predictions', ylab = 'Test Age Sample Collection',
     main = 'gender+gdna.exon.intron')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_gdna.exon.intron[[4]]) ~ unlist(rf_gdna.exon.intron[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.exon.intron[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))



# just methylation with log
rf_gdna.exon.intron_log <- predictAll(model_name = 'rf',
                                          data = full_data_rf,
                                          subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.exon.intron"),
                                          selected_features = c("gender", "gdna.exon.intron"),
                                          log = TRUE,
                                          iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_gdna.exon.intron_log[[3]]), unlist(rf_gdna.exon.intron_log[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just gdna.exon.intron with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_gdna.exon.intron_log[[3]]) ~ unlist(rf_gdna.exon.intron_log[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.exon.intron_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_gdna.exon.intron_log[[4]]), unlist(rf_gdna.exon.intron_log[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+gdna.exon.intron with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_gdna.exon.intron_log[[4]]) ~ unlist(rf_gdna.exon.intron_log[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.exon.intron_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_gdna.exon.intron_log[[3]]), unlist(rf_gdna.exon.intron_log[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just gdna.exon.intron with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_gdna.exon.intron_log[[3]]) ~ unlist(rf_gdna.exon.intron_log[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.exon.intron_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_gdna.exon.intron_log[[4]]), unlist(rf_gdna.exon.intron_log[[8]]), xlab = 'Test Predictions', ylab = 'Test Age Sample Collection',
     main = 'gender+gdna.exon.intron with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_gdna.exon.intron_log[[4]]) ~ unlist(rf_gdna.exon.intron_log[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.exon.intron_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


##################################################################################################
# just codon72.npro
rf_codon72.npro <- predictAll(model_name = 'rf',
                                  data = full_data_rf,
                                  subset <- c("age_diagnosis", "age_sample_collection", "gender", "codon72.npro"),
                                  selected_features = c("gender", "codon72.npro"),
                                  log = FALSE,
                                  iterations = 20)


###############################
# plot train age diagnosis
# plot(unlist(rf_codon72.npro[[3]]), unlist(rf_codon72.npro[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just codon72.npro no log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_codon72.npro[[3]]) ~ unlist(rf_codon72.npro[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_codon72.npro[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_codon72.npro[[4]]), unlist(rf_codon72.npro[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+codon72.npro')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_codon72.npro[[4]]) ~ unlist(rf_codon72.npro[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_codon72.npro[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_codon72.npro[[3]]), unlist(rf_codon72.npro[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just codon72.npro')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_codon72.npro[[3]]) ~ unlist(rf_codon72.npro[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_codon72.npro[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_codon72.npro[[4]]), unlist(rf_codon72.npro[[8]]), xlab = 'Test Predictions', ylab = 'Test Age Sample Collection',
     main = 'gender+codon72.npro')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_codon72.npro[[4]]) ~ unlist(rf_codon72.npro[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_codon72.npro[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))



# just methylation with log
rf_codon72.npro_log <- predictAll(model_name = 'rf',
                                      data = full_data_rf,
                                      subset <- c("age_diagnosis", "age_sample_collection", "gender", "codon72.npro"),
                                      selected_features = c("gender", "codon72.npro"),
                                      log = TRUE,
                                      iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_codon72.npro_log[[3]]), unlist(rf_codon72.npro_log[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just codon72.npro with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_codon72.npro_log[[3]]) ~ unlist(rf_codon72.npro_log[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_codon72.npro_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_codon72.npro_log[[4]]), unlist(rf_codon72.npro_log[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+codon72.npro with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_codon72.npro_log[[4]]) ~ unlist(rf_codon72.npro_log[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_codon72.npro_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_codon72.npro_log[[3]]), unlist(rf_codon72.npro_log[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just codon72.npro with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_codon72.npro_log[[3]]) ~ unlist(rf_codon72.npro_log[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_codon72.npro_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_codon72.npro_log[[4]]), unlist(rf_codon72.npro_log[[8]]), xlab = 'Test Predictions', ylab = 'Test Age Sample Collection',
     main = 'gender+codon72.npro with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_codon72.npro_log[[4]]) ~ unlist(rf_codon72.npro_log[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_codon72.npro_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


##################################################################################################
# just splice.delins.snv
rf_splice.delins.snv <- predictAll(model_name = 'rf',
                              data = full_data_rf,
                              subset <- c("age_diagnosis", "age_sample_collection", "gender", "splice.delins.snv"),
                              selected_features = c("gender", "splice.delins.snv"),
                              log = FALSE,
                              iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_splice.delins.snv[[3]]), unlist(rf_splice.delins.snv[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just splice.delins.snv no log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_splice.delins.snv[[3]]) ~ unlist(rf_splice.delins.snv[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_splice.delins.snv[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_splice.delins.snv[[4]]), unlist(rf_splice.delins.snv[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+splice.delins.snv')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_splice.delins.snv[[4]]) ~ unlist(rf_splice.delins.snv[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_splice.delins.snv[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_splice.delins.snv[[3]]), unlist(rf_splice.delins.snv[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just splice.delins.snv')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_splice.delins.snv[[3]]) ~ unlist(rf_splice.delins.snv[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_splice.delins.snv[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_splice.delins.snv[[4]]), unlist(rf_splice.delins.snv[[8]]), xlab = 'Test Predictions', ylab = 'Test Age of Sample Collection',
     main = 'gender+splice.delins.snv')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_splice.delins.snv[[4]]) ~ unlist(rf_splice.delins.snv[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_splice.delins.snv[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))



# just methylation with log
rf_splice.delins.snv_log <- predictAll(model_name = 'rf',
                                  data = full_data_rf,
                                  subset <- c("age_diagnosis", "age_sample_collection", "gender", "splice.delins.snv"),
                                  selected_features = c("gender", "splice.delins.snv"),
                                  log = TRUE,
                                  iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_splice.delins.snv_log[[3]]), unlist(rf_splice.delins.snv_log[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just splice.delins.snv with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_splice.delins.snv_log[[3]]) ~ unlist(rf_splice.delins.snv_log[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_splice.delins.snv_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_splice.delins.snv_log[[4]]), unlist(rf_splice.delins.snv_log[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+splice.delins.snv with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_splice.delins.snv_log[[4]]) ~ unlist(rf_splice.delins.snv_log[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_splice.delins.snv_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_splice.delins.snv_log[[3]]), unlist(rf_splice.delins.snv_log[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just splice.delins.snv with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_splice.delins.snv_log[[3]]) ~ unlist(rf_splice.delins.snv_log[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_splice.delins.snv_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_splice.delins.snv_log[[4]]), unlist(rf_splice.delins.snv_log[[8]]), xlab = 'Test Predictions', ylab = 'Test Age of Sample Collection',
     main = 'Just splice.delins.snv with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_splice.delins.snv_log[[4]]) ~ unlist(rf_splice.delins.snv_log[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_splice.delins.snv_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


##################################################################################################
# just protein.codon.num
rf_protein.codon.num <- predictAll(model_name = 'rf',
                                   data = full_data_rf,
                                   subset <- c("age_diagnosis", "age_sample_collection", "gender", "protein.codon.num"),
                                   selected_features = c("gender", "protein.codon.num"),
                                   log = FALSE,
                                   iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_protein.codon.num[[3]]), unlist(rf_protein.codon.num[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just protein.codon.num no log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_protein.codon.num[[3]]) ~ unlist(rf_protein.codon.num[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.num[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_protein.codon.num[[4]]), unlist(rf_protein.codon.num[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+protein.codon.num')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_protein.codon.num[[4]]) ~ unlist(rf_protein.codon.num[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.num[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_protein.codon.num[[3]]), unlist(rf_protein.codon.num[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just protein.codon.num')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_protein.codon.num[[3]]) ~ unlist(rf_protein.codon.num[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.num[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_protein.codon.num[[4]]), unlist(rf_protein.codon.num[[8]]), xlab = 'Test Predictions', ylab = 'Test Age of Sample Collection',
     main = 'gender+protein.codon.num')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_protein.codon.num[[4]]) ~ unlist(rf_protein.codon.num[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.num[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))



# just methylation with log
rf_protein.codon.num_log <- predictAll(model_name = 'rf',
                                       data = full_data_rf,
                                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "protein.codon.num"),
                                       selected_features = c("gender", "protein.codon.num"),
                                       log = TRUE,
                                       iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_protein.codon.num_log[[3]]), unlist(rf_protein.codon.num_log[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just protein.codon.num with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_protein.codon.num_log[[3]]) ~ unlist(rf_protein.codon.num_log[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.num_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_protein.codon.num_log[[4]]), unlist(rf_protein.codon.num_log[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+protein.codon.num with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_protein.codon.num_log[[4]]) ~ unlist(rf_protein.codon.num_log[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.num_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_protein.codon.num_log[[3]]), unlist(rf_protein.codon.num_log[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just protein.codon.num with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_protein.codon.num_log[[3]]) ~ unlist(rf_protein.codon.num_log[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.num_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_protein.codon.num_log[[4]]), unlist(rf_protein.codon.num_log[[8]]), xlab = 'Test Predictions', ylab = 'Test Age of Sample Collection',
     main = 'Just protein.codon.num with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_protein.codon.num_log[[4]]) ~ unlist(rf_protein.codon.num_log[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.num_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


##################################################################################################
# just mdm2.nG
rf_mdm2.nG <- predictAll(model_name = 'rf',
                                   data = full_data_rf,
                                   subset <- c("age_diagnosis", "age_sample_collection", "gender", "mdm2.nG"),
                                   selected_features = c("gender", "mdm2.nG"),
                                   log = FALSE,
                                   iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_mdm2.nG[[3]]), unlist(rf_mdm2.nG[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just mdm2.nG no log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_mdm2.nG[[3]]) ~ unlist(rf_mdm2.nG[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_mdm2.nG[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_mdm2.nG[[4]]), unlist(rf_mdm2.nG[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+mdm2.nG')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_mdm2.nG[[4]]) ~ unlist(rf_mdm2.nG[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_mdm2.nG[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# plot train age sample collection
# plot(unlist(rf_mdm2.nG[[3]]), unlist(rf_mdm2.nG[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just mdm2.nG')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_mdm2.nG[[3]]) ~ unlist(rf_mdm2.nG[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_mdm2.nG[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_mdm2.nG[[4]]), unlist(rf_mdm2.nG[[8]]), xlab = 'Test Predictions', ylab = 'Test Age of Sample Collection',
     main = 'gender+mdm2.nG')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_mdm2.nG[[4]]) ~ unlist(rf_mdm2.nG[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_mdm2.nG[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))



# just methylation with log
rf_mdm2.nG_log <- predictAll(model_name = 'rf',
                                       data = full_data_rf,
                                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "mdm2.nG"),
                                       selected_features = c("gender", "mdm2.nG"),
                                       log = TRUE,
                                       iterations = 20)


###############################
# # plot train age diagnosis
# plot(unlist(rf_mdm2.nG_log[[3]]), unlist(rf_mdm2.nG_log[[5]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just mdm2.nG with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_mdm2.nG_log[[3]]) ~ unlist(rf_mdm2.nG_log[[5]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_mdm2.nG_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age diagnosis
plot(unlist(rf_mdm2.nG_log[[4]]), unlist(rf_mdm2.nG_log[[6]]), xlab = 'Test Predictions', ylab = 'Test Age of Diagnosis',
     main = 'gender+mdm2.nG with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_mdm2.nG_log[[4]]) ~ unlist(rf_mdm2.nG_log[[6]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_mdm2.nG_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


###############################
# # plot train age sample collection
# plot(unlist(rf_mdm2.nG_log[[3]]), unlist(rf_mdm2.nG_log[[7]]), xlab = 'Train Predictions', ylab = 'Train Age of Diagnosis',
#      main = 'Just mdm2.nG with log')
# abline(0,1)
# r_squared <- round(summary(lm(unlist(rf_mdm2.nG_log[[3]]) ~ unlist(rf_mdm2.nG_log[[7]])))$adj.r.squared, 2)
# legend("bottomright", legend = paste0('# obs = ', rf_mdm2.nG_log[[11]]))
# legend("topleft", legend = paste0('r_squared = ', r_squared))

# plot test age sample collection
plot(unlist(rf_mdm2.nG_log[[4]]), unlist(rf_mdm2.nG_log[[8]]), xlab = 'Test Predictions', ylab = 'Test Age of Sample Collection',
     main = 'Just mdm2.nG with log')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_mdm2.nG_log[[4]]) ~ unlist(rf_mdm2.nG_log[[8]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_mdm2.nG_log[[11]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))



####################################################################
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
# Run with #gdna.codon, protein.codon.num, mdm2.nG, protein.codon.change
# return(list(train.mse, test.mse,  train.predictions, test.predictions, train.ground_truth, 
#             test.ground_truth, train.sample_collection,
#             test.sample_collection, model, importance, obs))

# #########################
# # combinations
# # add gdna.codon and protein.codon.change
# rf_gdna.codon_protein.codon.change <- predictAll(model_name = 'rf', 
#                             data = full_data,
#                             subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.codon", "protein.codon.change"), 
#                             selected_features = c("gender", "gdna.codon", "protein.codon.change"), 
#                             iterations = 10)
# 
# 
# plot(unlist(rf_gdna.codon_protein.codon.change[[2]]), unlist(rf_gdna.codon_protein.codon.change[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.codon+protein.codon.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_protein.codon.change[[6]]))
# 
# plot(unlist(rf_gdna.codon_protein.codon.change[[2]]), unlist(rf_gdna.codon_protein.codon.change[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.codon+protein.codon.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_protein.codon.change[[6]]))
# 
# 
# # add protein.codon.num
# rf_add_protein.codon.num <- predictAll(model_name = 'rf', 
#                                    data = full_data,
#                                    subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num"), 
#                                    selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num"), 
#                                    iterations = 10)
# 
# plot(unlist(rf_add_protein.codon.num[[2]]), unlist(rf_add_protein.codon.num[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_add_protein.codon.num[[6]]))
# 
# plot(unlist(rf_add_protein.codon.num[[2]]), unlist(rf_add_protein.codon.num[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_add_protein.codon.num[[6]]))
# 
# # add mdm2.nG
# rf_add_mdm2.nG <- predictAll(model_name = 'rf', 
#                                    data = full_data,
#                                    subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num", "mdm2.nG"), 
#                                    selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num", "mdm2.nG"), 
#                                    iterations = 10)
# 
# plot(unlist(rf_add_mdm2.nG[[2]]), unlist(rf_add_mdm2.nG[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_add_mdm2.nG[[6]]))
# 
# plot(unlist(rf_add_mdm2.nG[[2]]), unlist(rf_add_mdm2.nG[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_add_mdm2.nG[[6]]))
# 
# dev.off()
# ########################################################################################
# Corelation data
# pdf('/home/benbrew/Desktop/rf_methylation_cor.pdf')
# 
# # just methylation
# rf_methyl_cor <- predictAll(model_name = 'rf', 
#                         data = full_data_cor,
#                         subset <- c("age_diagnosis", "age_sample_collection"), 
#                         selected_features = NULL, iterations = 10)
# 
# 
# plot(unlist(rf_methyl_cor[[2]]), unlist(rf_methyl_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = 'Just methylation')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_methyl_cor[[6]]))
# 
# plot(unlist(rf_methyl_cor[[2]]), unlist(rf_methyl_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = 'Just methylation')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_methyl_cor[[6]]))
# 
# # gender and gdna.base.change
# rf_mut_cor <- predictAll(model_name = 'rf', 
#                      data = full_data_cor,
#                      subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change"), 
#                      selected_features = c("gender", "gdna.base.change"), iterations = 10)
# 
# plot(unlist(rf_mut_cor[[2]]), unlist(rf_mut_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut_cor[[6]]))
# 
# plot(unlist(rf_mut_cor[[2]]), unlist(rf_mut_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut_cor[[6]]))
# 
# 
# 
# # add gdna.codon
# rf_mut1_cor <- predictAll(model_name = 'rf', 
#                       data = full_data_cor,
#                       subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.base.change", "gdna.codon"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 10)
# 
# plot(unlist(rf_mut1_cor[[2]]), unlist(rf_mut1_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut1_cor[[6]]))
# 
# plot(unlist(rf_mut1_cor[[2]]), unlist(rf_mut1_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut1_cor[[6]]))
# 
# 
# # add protein.codon.change
# rf_mut2_cor <- predictAll(model_name = 'rf', 
#                       data = full_data_cor,
#                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
#                       iterations = 10)
# 
# plot(unlist(rf_mut2_cor[[2]]), unlist(rf_mut2_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut2_cor[[6]]))
# 
# plot(unlist(rf_mut2_cor[[2]]), unlist(rf_mut2_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut2_cor[[6]]))
# 
# 
# 
# # add gdna.exon.intron
# rf_mut3_cor <- predictAll(model_name = 'rf', 
#                       data = full_data_cor,
#                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                   "gdna.exon.intron"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                             "gdna.exon.intron"), iterations = 10)
# 
# plot(unlist(rf_mut3_cor[[2]]), unlist(rf_mut3_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut3_cor[[6]]))
# 
# plot(unlist(rf_mut3_cor[[2]]), unlist(rf_mut3_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut3_cor[[6]]))
# 
# # add codon72.npro
# rf_mut4_cor <- predictAll(model_name = 'rf', 
#                       data = full_data_cor,
#                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                   "gdna.exon.intron", "codon72.npro"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                             "gdna.exon.intron", "codon72.npro"), iterations = 10)
# 
# plot(unlist(rf_mut4_cor[[2]]), unlist(rf_mut4_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut4_cor[[6]]))
# 
# plot(unlist(rf_mut4_cor[[2]]), unlist(rf_mut4_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut4_cor[[6]]))
# 
# 
# 
# # add splice.delins.snv
# rf_mut5_cor <- predictAll(model_name = 'rf', 
#                       data = full_data_cor,
#                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                   "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                             "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
#                       iterations = 10)
# 
# plot(unlist(rf_mut5_cor[[2]]), unlist(rf_mut5_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut5_cor[[6]]))
# 
# plot(unlist(rf_mut5_cor[[2]]), unlist(rf_mut5_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut5_cor[[6]]))
# 
# 
# # add protein.codon.num
# rf_mut6_cor <- predictAll(model_name = 'rf', 
#                       data = full_data_cor,
#                       subset <- c("age_diagnosis","age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                   "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                             "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
#                       iterations = 10)
# 
# plot(unlist(rf_mut6_cor[[2]]), unlist(rf_mut6_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut6_cor[[6]]))
# 
# plot(unlist(rf_mut6_cor[[2]]), unlist(rf_mut6_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut6_cor[[6]]))
# 
# 
# # add mdm2.nG
# rf_mut7_cor <- predictAll(model_name = 'rf', 
#                       data = full_data_cor,
#                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                   "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                             "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                       iterations = 10)
# 
# plot(unlist(rf_mut7_cor[[2]]), unlist(rf_mut7_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
#      mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut7_cor[[6]]))
# 
# plot(unlist(rf_mut7_cor[[2]]), unlist(rf_mut7_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
#      mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut7_cor[[6]]))
# 
# 
# 
# ####################################################################
# # variables missing
# # gender 0
# # gdna.base.change 164
# # gdna.codon 164
# # protein.codon.change 177
# # gdna.exon.intron 492
# # codon72.npro 517
# # splice.delins.snv 519
# # protein.codon.num 549
# # mdm2.nG 652
# # Run with #gdna.codon, protein.codon.num, mdm2.nG, protein.codon.change
# # add mdm2.nG
# rf_mdm2.nG_cor <- predictAll(model_name = 'rf', 
#                          data = full_data_cor,
#                          subset <- c("age_diagnosis", "age_sample_collection", "gender", "mdm2.nG"), 
#                          selected_features = c("gender", "mdm2.nG"), iterations = 10)
# 
# plot(unlist(rf_mdm2.nG_cor[[2]]), unlist(rf_mdm2.nG_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mdm2.nG_cor[[6]]))
# 
# plot(unlist(rf_mdm2.nG_cor[[2]]), unlist(rf_mdm2.nG_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mdm2.nG_cor[[6]]))
# 
# # add gdna.codon
# rf_gdna.codon_cor <- predictAll(model_name = 'rf', 
#                             data = full_data_cor,
#                             subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon"), 
#                             selected_features = c("gender", "gdna.codon"), iterations = 10)
# 
# plot(unlist(rf_gdna.codon_cor[[2]]), unlist(rf_gdna.codon_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.codon')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_cor[[6]]))
# 
# plot(unlist(rf_gdna.codon_cor[[2]]), unlist(rf_gdna.codon_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.codon')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_cor[[6]]))
# 
# # add protein.codon.num
# rf_protein.codon.num_cor <- predictAll(model_name = 'rf', 
#                                    data = full_data_cor,
#                                    subset <- c("age_diagnosis", "age_sample_collection", "gender", "protein.codon.num"), 
#                                    selected_features = c("gender", "protein.codon.num"), iterations = 10)
# 
# plot(unlist(rf_protein.codon.num_cor[[2]]), unlist(rf_protein.codon.num_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.num_cor[[6]]))
# 
# plot(unlist(rf_protein.codon.num_cor[[2]]), unlist(rf_protein.codon.num_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.num_cor[[6]]))
# 
# # add protein.codon.change
# rf_protein.codon.change_cor <- predictAll(model_name = 'rf', 
#                                       data = full_data_cor,
#                                       subset <- c("age_diagnosis","age_sample_collection",  "gender", "protein.codon.change"), 
#                                       selected_features = c("gender", "protein.codon.change"), 
#                                       iterations = 10)
# 
# plot(unlist(rf_protein.codon.change_cor[[2]]), unlist(rf_protein.codon.change_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.change_cor[[6]]))
# 
# plot(unlist(rf_protein.codon.change_cor[[2]]), unlist(rf_protein.codon.change_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.change_cor[[6]]))
# 
# #########################
# # combinations
# # add gdna.codon and protein.codon.change
# rf_gdna.codon_protein.codon.change_cor <- predictAll(model_name = 'rf', 
#                                                  data = full_data_cor,
#                                                  subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.codon", "protein.codon.change"), 
#                                                  selected_features = c("gender", "gdna.codon", "protein.codon.change"), 
#                                                  iterations = 10)
# 
# 
# plot(unlist(rf_gdna.codon_protein.codon.change_cor[[2]]), unlist(rf_gdna.codon_protein.codon.change_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.codon+protein.codon.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_protein.codon.change_cor[[6]]))
# 
# plot(unlist(rf_gdna.codon_protein.codon.change_cor[[2]]), unlist(rf_gdna.codon_protein.codon.change_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.codon+protein.codon.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_protein.codon.change_cor[[6]]))
# 
# 
# # add protein.codon.num
# rf_add_protein.codon.num_cor <- predictAll(model_name = 'rf', 
#                                        data = full_data_cor,
#                                        subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num"), 
#                                        selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num"), 
#                                        iterations = 10)
# 
# plot(unlist(rf_add_protein.codon.num_cor[[2]]), unlist(rf_add_protein.codon.num_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_add_protein.codon.num_cor[[6]]))
# 
# plot(unlist(rf_add_protein.codon.num_cor[[2]]), unlist(rf_add_protein.codon.num_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_add_protein.codon.num_cor[[6]]))
# 
# # add mdm2.nG
# rf_add_mdm2.nG_cor <- predictAll(model_name = 'rf', 
#                              data = full_data_cor,
#                              subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num", "mdm2.nG"), 
#                              selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num", "mdm2.nG"), 
#                              iterations = 10)
# 
# plot(unlist(rf_add_mdm2.nG_cor[[2]]), unlist(rf_add_mdm2.nG_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_add_mdm2.nG_cor[[6]]))
# 
# plot(unlist(rf_add_mdm2.nG_cor[[2]]), unlist(rf_add_mdm2.nG_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_add_mdm2.nG_cor[[6]]))
# 
# dev.off()
################################################################################################
# rf
pdf('/home/benbrew/Desktop/rf_methylation_reg.pdf')
# just methylation
rf_methyl_rf <- predictAll(model_name = 'rf',
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "age_sample_collection"),
                            selected_features = NULL, iterations = 10)


plot(unlist(rf_methyl_rf[[2]]), unlist(rf_methyl_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = 'Just methylation')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_methyl_rf[[2]]) ~ unlist(rf_methyl_rf[[5]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_methyl_rf[[6]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


plot(unlist(rf_methyl_rf[[2]]), unlist(rf_methyl_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = 'Just methylation')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_methyl_rf[[2]]) ~ unlist(rf_methyl_rf[[7]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_methyl_rf[[6]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))
# # gender and gdna.base.change
# rf_mut_rf <- predictAll(model_name = 'rf', 
#                          data = full_data_rf,
#                          subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change"), 
#                          selected_features = c("gender", "gdna.base.change"), iterations = 10)
# 
# plot(unlist(rf_mut_rf[[2]]), unlist(rf_mut_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut_rf[[6]]))
# 
# plot(unlist(rf_mut_rf[[2]]), unlist(rf_mut_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut_rf[[6]]))
# 
# 
# 
# # add gdna.codon
# rf_mut1_rf <- predictAll(model_name = 'rf', 
#                           data = full_data_rf,
#                           subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.base.change", "gdna.codon"), 
#                           selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 10)
# 
# plot(unlist(rf_mut1_rf[[2]]), unlist(rf_mut1_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut1_rf[[6]]))
# 
# plot(unlist(rf_mut1_rf[[2]]), unlist(rf_mut1_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut1_rf[[6]]))
# 
# 
# # add protein.codon.change
# rf_mut2_rf <- predictAll(model_name = 'rf', 
#                       data = full_data_rf,
#                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
#                       iterations = 10)
# 
# plot(unlist(rf_mut2_rf[[2]]), unlist(rf_mut2_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut2_rf[[6]]))
# 
# plot(unlist(rf_mut2_rf[[2]]), unlist(rf_mut2_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut2_rf[[6]]))
# 
# 
# 
# # add gdna.exon.intron
# rf_mut3_rf <- predictAll(model_name = 'rf', 
#                       data = full_data_rf,
#                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                   "gdna.exon.intron"), 
#                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                             "gdna.exon.intron"), iterations = 10)
# 
# plot(unlist(rf_mut3_rf[[2]]), unlist(rf_mut3_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut3_rf[[6]]))
# 
# plot(unlist(rf_mut3_rf[[2]]), unlist(rf_mut3_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut3_rf[[6]]))
# 
# # add codon72.npro
# rf_mut4_rf <- predictAll(model_name = 'rf', 
#                           data = full_data_rf,
#                           subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                       "gdna.exon.intron", "codon72.npro"), 
#                           selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                                 "gdna.exon.intron", "codon72.npro"), iterations = 10)
# 
# plot(unlist(rf_mut4_rf[[2]]), unlist(rf_mut4_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut4_rf[[6]]))
# 
# plot(unlist(rf_mut4_rf[[2]]), unlist(rf_mut4_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut4_rf[[6]]))
# 
# 
# 
# # add splice.delins.snv
# rf_mut5_rf <- predictAll(model_name = 'rf', 
#                           data = full_data_rf,
#                           subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                       "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
#                           selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                                 "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
#                           iterations = 10)
# 
# plot(unlist(rf_mut5_rf[[2]]), unlist(rf_mut5_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut5_rf[[6]]))
# 
# plot(unlist(rf_mut5_rf[[2]]), unlist(rf_mut5_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut5_rf[[6]]))
# 
# 
# # add protein.codon.num
# rf_mut6_rf <- predictAll(model_name = 'rf', 
#                           data = full_data_rf,
#                           subset <- c("age_diagnosis","age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                       "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
#                           selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                                 "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
#                           iterations = 10)
# 
# plot(unlist(rf_mut6_rf[[2]]), unlist(rf_mut6_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut6_rf[[6]]))
# 
# plot(unlist(rf_mut6_rf[[2]]), unlist(rf_mut6_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut6_rf[[6]]))
# 
# 
# # add mdm2.nG
# rf_mut7_rf <- predictAll(model_name = 'rf', 
#                           data = full_data_rf,
#                           subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                       "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                           selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                                 "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                           iterations = 10)
# 
# plot(unlist(rf_mut7_rf[[2]]), unlist(rf_mut7_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
#      mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut7_rf[[6]]))
# 
# plot(unlist(rf_mut7_rf[[2]]), unlist(rf_mut7_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
#      gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
#      mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_mut7_rf[[6]]))
# 
# 

####################################################################
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
# Run with #gdna.codon, protein.codon.num, mdm2.nG, protein.codon.change
# add mdm2.nG
# add mdm2.nG
rf_mdm2.nG_rf <- predictAll(model_name = 'rf', 
                         data = full_data_rf,
                         subset <- c("age_diagnosis", "age_sample_collection", "gender", "mdm2.nG"), 
                         selected_features = c("gender", "mdm2.nG"), iterations = 10)

plot(unlist(rf_mdm2.nG_rf[[2]]), unlist(rf_mdm2.nG_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+mdm2.nG')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_mdm2.nG_rf[[2]]) ~ unlist(rf_mdm2.nG_rf[[5]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))
legend("bottomright", legend = paste0('# obs = ', rf_mdm2.nG_rf[[6]]))



plot(unlist(rf_mdm2.nG_rf[[2]]), unlist(rf_mdm2.nG_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+mdm2.nG')
abline(0,1)
r_squared <- round(summary(lm(unlist(rf_mdm2.nG_rf[[2]]) ~ unlist(rf_mdm2.nG_rf[[7]])))$adj.r.squared, 2)
legend("bottomright", legend = paste0('# obs = ', rf_mdm2.nG_rf[[6]]))
legend("topleft", legend = paste0('r_squared = ', r_squared))


# add gdna.codon
rf_gdna.codon_rf <- predictAll(model_name = 'rf', 
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon"), 
                            selected_features = c("gender", "gdna.codon"), iterations = 10)

plot(unlist(rf_gdna.codon_rf[[2]]), unlist(rf_gdna.codon_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_gdna.codon_rf[[2]]) ~ unlist(rf_gdna.codon_rf[[5]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_gdna.codon_rf[[2]]), unlist(rf_gdna.codon_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_gdna.codon_rf[[2]]) ~ unlist(rf_gdna.codon_rf[[7]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))

# add protein.codon.num
rf_protein.codon.num_rf <- predictAll(model_name = 'rf', 
                                   data = full_data_rf,
                                   subset <- c("age_diagnosis", "age_sample_collection", "gender", "protein.codon.num"), 
                                   selected_features = c("gender", "protein.codon.num"), iterations = 10)

plot(unlist(rf_protein.codon.num_rf[[2]]), unlist(rf_protein.codon.num_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.num_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_protein.codon.num_rf[[2]]) ~ unlist(rf_protein.codon.num_rf[[5]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_protein.codon.num_rf[[2]]), unlist(rf_protein.codon.num_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.num_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_protein.codon.num_rf[[2]]) ~ unlist(rf_protein.codon.num_rf[[7]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))


# add protein.codon.change
rf_protein.codon.change_rf <- predictAll(model_name = 'rf', 
                                      data = full_data_rf,
                                      subset <- c("age_diagnosis","age_sample_collection",  "gender", "protein.codon.change"), 
                                      selected_features = c("gender", "protein.codon.change"), 
                                      iterations = 10)

plot(unlist(rf_protein.codon.change_rf[[2]]), unlist(rf_protein.codon.change_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.change_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_protein.codon.change_rf[[2]]) ~ unlist(rf_protein.codon.change_rf[[5]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_protein.codon.change[[2]]), unlist(rf_protein.codon.change[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_protein.codon.change_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_protein.codon.change_rf[[2]]) ~ unlist(rf_protein.codon.change_rf[[7]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))


# add gdna.base.change
rf_gdna.base.change_rf <- predictAll(model_name = 'rf', 
                                  data = full_data_rf,
                                  subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.base.change"), 
                                  selected_features = c("gender", "gdna.base.change"), 
                                  iterations = 10)

plot(unlist(rf_gdna.base.change_rf[[2]]), unlist(rf_gdna.base.change_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.base.change_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_gdna.base.change_rf[[2]]) ~ unlist(rf_gdna.base.change_rf[[5]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_gdna.base.change_rf[[2]]), unlist(rf_gdna.base.change_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.base.change_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_gdna.base.change_rf[[2]]) ~ unlist(rf_gdna.base.change_rf[[7]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))


#splice.delins.snv
rf_splice.delins.snv_rf <- predictAll(model_name = 'rf', 
                                   data = full_data_rf,
                                   subset <- c("age_diagnosis","age_sample_collection",  "gender", "splice.delins.snv"), 
                                   selected_features = c("gender", "splice.delins.snv"), 
                                   iterations = 10)

plot(unlist(rf_splice.delins.snv_rf[[2]]), unlist(rf_splice.delins.snv_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_splice.delins.snv_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_splice.delins.snv_rf[[2]]) ~ unlist(rf_splice.delins.snv_rf[[5]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_splice.delins.snv_rf[[2]]), unlist(rf_splice.delins.snv_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_splice.delins.snv_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_splice.delins.snv_rf[[2]]) ~ unlist(rf_splice.delins.snv_rf[[7]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))


#codon72.npro
rf_codon72.npro_rf <- predictAll(model_name = 'rf', 
                              data = full_data_rf,
                              subset <- c("age_diagnosis","age_sample_collection",  "gender", "codon72.npro"), 
                              selected_features = c("gender", "codon72.npro"), 
                              iterations = 10)

plot(unlist(rf_codon72.npro_rf[[2]]), unlist(rf_codon72.npro_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_codon72.npro_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_codon72.npro_rf[[2]]) ~ unlist(rf_codon72.npro_rf[[5]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_codon72.npro_rf[[2]]), unlist(rf_codon72.npro_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_codon72.npro_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_codon72.npro_rf[[2]]) ~ unlist(rf_codon72.npro_rf[[7]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))

# gdna.exon.intron

rf_gdna.exon.intron_rf <- predictAll(model_name = 'rf', 
                                  data = full_data_rf,
                                  subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.exon.intron"), 
                                  selected_features = c("gender", "gdna.exon.intron"), 
                                  iterations = 10)

plot(unlist(rf_gdna.exon.intron_rf[[2]]), unlist(rf_gdna.exon.intron_rf[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.exon.intron_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_gdna.exon.intron_rf[[2]]) ~ unlist(rf_gdna.exon.intron_rf[[5]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))

plot(unlist(rf_gdna.exon.intron_rf[[2]]), unlist(rf_gdna.exon.intron_rf[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', rf_gdna.exon.intron_rf[[6]]))
r_squared <- round(summary(lm(unlist(rf_gdna.exon.intron_rf[[2]]) ~ unlist(rf_gdna.exon.intron_rf[[7]])))$adj.r.squared, 2)
legend("topleft", legend = paste0('r_squared = ', r_squared))

# #########################
# # combinations
# # add gdna.codon and protein.codon.change
# rf_gdna.codon_protein.codon.change <- predictAll(model_name = 'rf', 
#                             data = full_data,
#                             subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.codon", "protein.codon.change"), 
#                             selected_features = c("gender", "gdna.codon", "protein.codon.change"), 
#                             iterations = 10)
# 
# 
# plot(unlist(rf_gdna.codon_protein.codon.change[[2]]), unlist(rf_gdna.codon_protein.codon.change[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.codon+protein.codon.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_protein.codon.change[[6]]))
# 
# plot(unlist(rf_gdna.codon_protein.codon.change[[2]]), unlist(rf_gdna.codon_protein.codon.change[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.codon+protein.codon.change')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_gdna.codon_protein.codon.change[[6]]))
# 
# 
# # add protein.codon.num
# rf_add_protein.codon.num <- predictAll(model_name = 'rf', 
#                                    data = full_data,
#                                    subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num"), 
#                                    selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num"), 
#                                    iterations = 10)
# 
# plot(unlist(rf_add_protein.codon.num[[2]]), unlist(rf_add_protein.codon.num[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_add_protein.codon.num[[6]]))
# 
# plot(unlist(rf_add_protein.codon.num[[2]]), unlist(rf_add_protein.codon.num[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_add_protein.codon.num[[6]]))
# 
# # add mdm2.nG
# rf_add_mdm2.nG <- predictAll(model_name = 'rf', 
#                                    data = full_data,
#                                    subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num", "mdm2.nG"), 
#                                    selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num", "mdm2.nG"), 
#                                    iterations = 10)
# 
# plot(unlist(rf_add_mdm2.nG[[2]]), unlist(rf_add_mdm2.nG[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
#      main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_add_mdm2.nG[[6]]))
# 
# plot(unlist(rf_add_mdm2.nG[[2]]), unlist(rf_add_mdm2.nG[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
#      main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
# abline(0,1)
# legend("bottomright", legend = paste0('# obs = ', rf_add_mdm2.nG[[6]]))
# 
# dev.off()
dev.off()
# # combine all tests 
# 
# rf_all <- rbind (
#   append('just_methyl', c(mean(unlist(rf_methyl[[1]])), rf_methyl[[6]], 'all_data')),
#   append('gender_and_gdna.base.change', c(mean(unlist(rf_mut[[1]])), rf_mut[[6]], 'all_data')),
#   append('gdna.codon', c(mean(unlist(rf_mut1[[1]])), rf_mut1[[6]], 'all_data')),
#   append('protein.codon.change', c(mean(unlist(rf_mut2[[1]])), rf_mut2[[6]], 'all_data')),
#   append('gdna.exon.intron', c(mean(unlist(rf_mut3[[1]])), rf_mut3[[6]], 'all_data')),
#   append('codon72.npro', c(mean(unlist(rf_mut4[[1]])), rf_mut4[[6]], 'all_data')),
#   append('splice.delins.snv', c(mean(unlist(rf_mut5[[1]])), rf_mut5[[6]], 'all_data')),
#   append('protein.codon.num', c(mean(unlist(rf_mut6[[1]])), rf_mut6[[6]], 'all_data')),
#   append('mdm2.nG', c(mean(unlist(rf_mut7[[1]])), rf_mut7[[6]], 'all_data')),
#   
#   append('just_methyl', c(mean(unlist(rf_methyl_cor[[1]])), rf_methyl_cor[[6]], 'corr')),
#   append('gender_and_gdna.base.change', c(mean(unlist(rf_mut_cor[[1]])), rf_mut_cor[[6]], 'corr')),
#   append('gdna.codon', c(mean(unlist(rf_mut1_cor[[1]])), rf_mut1_cor[[6]], 'corr')),
#   append('protein.codon.change', c(mean(unlist(rf_mut2_cor[[1]])), rf_mut2_cor[[6]], 'corr')),
#   append('gdna.exon.intron', c(mean(unlist(rf_mut3_cor[[1]])), rf_mut3_cor[[6]], 'corr')),
#   append('codon72.npro', c(mean(unlist(rf_mut4_cor[[1]])), rf_mut4_cor[[6]], 'corr')),
#   append('splice.delins.snv', c(mean(unlist(rf_mut5_cor[[1]])), rf_mut5_cor[[6]], 'corr')),
#   append('protein.codon.num', c(mean(unlist(rf_mut6_cor[[1]])), rf_mut6_cor[[6]], 'corr')),
#   append('mdm2.nG', c(mean(unlist(rf_mut7_cor[[1]])), rf_mut7_cor[[6]], 'corr')),
#   
#   append('just_methyl', c(mean(unlist(rf_methyl_reg[[1]])), rf_methyl_reg[[6]], 'recursive')),
#   append('gender_and_gdna.base.change', c(mean(unlist(rf_mut_reg[[1]])), rf_mut_reg[[6]], 'recursive')),
#   append('gdna.codon', c(mean(unlist(rf_mut1_reg[[1]])), rf_mut1_reg[[6]], 'recursive')),
#   append('protein.codon.change', c(mean(unlist(rf_mut2_reg[[1]])), rf_mut2_reg[[6]], 'recursive')),
#   append('gdna.exon.intron', c(mean(unlist(rf_mut3_reg[[1]])), rf_mut3_reg[[6]], 'recursive')),
#   append('codon72.npro', c(mean(unlist(rf_mut4_reg[[1]])), rf_mut4_reg[[6]], 'recursive')),
#   append('splice.delins.snv', c(mean(unlist(rf_mut5_reg[[1]])), rf_mut5_reg[[6]], 'recursive')),
#   append('protein.codon.num', c(mean(unlist(rf_mut6_reg[[1]])), rf_mut6_reg[[6]], 'recursive')),
#   append('mdm2.nG', c(mean(unlist(rf_mut7_reg[[1]])), rf_mut7_reg[[6]], 'recursive'))
# )




##################################################################################################################3
##################################################################################################################3
pdf('/home/benbrew/Desktop/enet_methylation.pdf')
# enet forest
#############################
# Full Data

# just methylation
enet_methyl <- predictAll(model_name = 'enet', 
                        data = full_data,
                        subset <- c("age_diagnosis", "age_sample_collection"), 
                        selected_features = NULL, iterations = 10)


plot(unlist(enet_methyl[[2]]), unlist(enet_methyl[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = 'Just methylation')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_methyl[[6]]))

plot(unlist(enet_methyl[[2]]), unlist(enet_methyl[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = 'Just methylation')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_methyl[[6]]))

# gender and gdna.base.change
enet_mut <- predictAll(model_name = 'enet', 
                     data = full_data,
                     subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change"), 
                     selected_features = c("gender", "gdna.base.change"), iterations = 10)

plot(unlist(enet_mut[[2]]), unlist(enet_mut[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut[[6]]))

plot(unlist(enet_mut[[2]]), unlist(enet_mut[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut[[6]]))



# add gdna.codon
enet_mut1 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.base.change", "gdna.codon"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 10)

plot(unlist(enet_mut1[[2]]), unlist(enet_mut1[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut1[[6]]))

plot(unlist(enet_mut1[[2]]), unlist(enet_mut1[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut1[[6]]))


# add protein.codon.change
enet_mut2 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                      iterations = 10)

plot(unlist(enet_mut2[[2]]), unlist(enet_mut2[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut2[[6]]))

plot(unlist(enet_mut2[[2]]), unlist(enet_mut2[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut2[[6]]))



# add gdna.exon.intron
enet_mut3 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron"), iterations = 10)

plot(unlist(enet_mut3[[2]]), unlist(enet_mut3[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut3[[6]]))

plot(unlist(enet_mut3[[2]]), unlist(enet_mut3[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut3[[6]]))

# add codon72.npro
enet_mut4 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro"), iterations = 10)

plot(unlist(enet_mut4[[2]]), unlist(enet_mut4[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut4[[6]]))

plot(unlist(enet_mut4[[2]]), unlist(enet_mut4[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut4[[6]]))



# add splice.delins.snv
enet_mut5 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                      iterations = 10)

plot(unlist(enet_mut5[[2]]), unlist(enet_mut5[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut5[[6]]))

plot(unlist(enet_mut5[[2]]), unlist(enet_mut5[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut5[[6]]))


# add protein.codon.num
enet_mut6 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis","age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                      iterations = 10)

plot(unlist(enet_mut6[[2]]), unlist(enet_mut6[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut6[[6]]))

plot(unlist(enet_mut6[[2]]), unlist(enet_mut6[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut6[[6]]))


# add mdm2.nG
enet_mut7 <- predictAll(model_name = 'enet', 
                      data = full_data,
                      subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                      iterations = 10)

plot(unlist(enet_mut7[[2]]), unlist(enet_mut7[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
     mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut7[[6]]))

plot(unlist(enet_mut7[[2]]), unlist(enet_mut7[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
     mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut7[[6]]))



####################################################################
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
# Run with #gdna.codon, protein.codon.num, mdm2.nG, protein.codon.change
# add mdm2.nG
enet_mdm2.nG <- predictAll(model_name = 'enet', 
                         data = full_data,
                         subset <- c("age_diagnosis", "age_sample_collection", "gender", "mdm2.nG"), 
                         selected_features = c("gender", "mdm2.nG"), iterations = 10)

plot(unlist(enet_mdm2.nG[[2]]), unlist(enet_mdm2.nG[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mdm2.nG[[6]]))

plot(unlist(enet_mdm2.nG[[2]]), unlist(enet_mdm2.nG[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mdm2.nG[[6]]))

# add gdna.codon
enet_gdna.codon <- predictAll(model_name = 'enet', 
                            data = full_data,
                            subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon"), 
                            selected_features = c("gender", "gdna.codon"), iterations = 10)

plot(unlist(enet_gdna.codon[[2]]), unlist(enet_gdna.codon[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_gdna.codon[[6]]))

plot(unlist(enet_gdna.codon[[2]]), unlist(enet_gdna.codon[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_gdna.codon[[6]]))

# add protein.codon.num
enet_protein.codon.num <- predictAll(model_name = 'enet', 
                                   data = full_data,
                                   subset <- c("age_diagnosis", "age_sample_collection", "gender", "protein.codon.num"), 
                                   selected_features = c("gender", "protein.codon.num"), iterations = 10)

plot(unlist(enet_protein.codon.num[[2]]), unlist(enet_protein.codon.num[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_protein.codon.num[[6]]))

plot(unlist(enet_protein.codon.num[[2]]), unlist(enet_protein.codon.num[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_protein.codon.num[[6]]))

# add protein.codon.change
enet_protein.codon.change <- predictAll(model_name = 'enet', 
                                      data = full_data,
                                      subset <- c("age_diagnosis","age_sample_collection",  "gender", "protein.codon.change"), 
                                      selected_features = c("gender", "protein.codon.change"), 
                                      iterations = 10)

plot(unlist(enet_protein.codon.change[[2]]), unlist(enet_protein.codon.change[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_protein.codon.change[[6]]))

plot(unlist(enet_protein.codon.change[[2]]), unlist(enet_protein.codon.change[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_protein.codon.change[[6]]))

#########################
# combinations
# add gdna.codon and protein.codon.change
enet_gdna.codon_protein.codon.change <- predictAll(model_name = 'enet', 
                                                 data = full_data,
                                                 subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.codon", "protein.codon.change"), 
                                                 selected_features = c("gender", "gdna.codon", "protein.codon.change"), 
                                                 iterations = 10)


plot(unlist(enet_gdna.codon_protein.codon.change[[2]]), unlist(enet_gdna.codon_protein.codon.change[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_gdna.codon_protein.codon.change[[6]]))

plot(unlist(enet_gdna.codon_protein.codon.change[[2]]), unlist(enet_gdna.codon_protein.codon.change[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_gdna.codon_protein.codon.change[[6]]))


# add protein.codon.num
enet_add_protein.codon.num <- predictAll(model_name = 'enet', 
                                       data = full_data,
                                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num"), 
                                       selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num"), 
                                       iterations = 10)

plot(unlist(enet_add_protein.codon.num[[2]]), unlist(enet_add_protein.codon.num[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_add_protein.codon.num[[6]]))

plot(unlist(enet_add_protein.codon.num[[2]]), unlist(enet_add_protein.codon.num[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_add_protein.codon.num[[6]]))

# add mdm2.nG
enet_add_mdm2.nG <- predictAll(model_name = 'enet', 
                             data = full_data,
                             subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num", "mdm2.nG"), 
                             selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num", "mdm2.nG"), 
                             iterations = 10)

plot(unlist(enet_add_mdm2.nG[[2]]), unlist(enet_add_mdm2.nG[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_add_mdm2.nG[[6]]))

plot(unlist(enet_add_mdm2.nG[[2]]), unlist(enet_add_mdm2.nG[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_add_mdm2.nG[[6]]))

dev.off()
########################################################################################
# Corelation data
pdf('/home/benbrew/Desktop/enet_methylation_cor.pdf')

# just methylation
enet_methyl_cor <- predictAll(model_name = 'enet', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "age_sample_collection"), 
                            selected_features = NULL, iterations = 10)


plot(unlist(enet_methyl_cor[[2]]), unlist(enet_methyl_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = 'Just methylation')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_methyl_cor[[6]]))

plot(unlist(enet_methyl_cor[[2]]), unlist(enet_methyl_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = 'Just methylation')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_methyl_cor[[6]]))

# gender and gdna.base.change
enet_mut_cor <- predictAll(model_name = 'enet', 
                         data = full_data_cor,
                         subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change"), 
                         selected_features = c("gender", "gdna.base.change"), iterations = 10)

plot(unlist(enet_mut_cor[[2]]), unlist(enet_mut_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut_cor[[6]]))

plot(unlist(enet_mut_cor[[2]]), unlist(enet_mut_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut_cor[[6]]))



# add gdna.codon
enet_mut1_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.base.change", "gdna.codon"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 10)

plot(unlist(enet_mut1_cor[[2]]), unlist(enet_mut1_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut1_cor[[6]]))

plot(unlist(enet_mut1_cor[[2]]), unlist(enet_mut1_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut1_cor[[6]]))


# add protein.codon.change
enet_mut2_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                          iterations = 10)

plot(unlist(enet_mut2_cor[[2]]), unlist(enet_mut2_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut2_cor[[6]]))

plot(unlist(enet_mut2_cor[[2]]), unlist(enet_mut2_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut2_cor[[6]]))



# add gdna.exon.intron
enet_mut3_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron"), iterations = 10)

plot(unlist(enet_mut3_cor[[2]]), unlist(enet_mut3_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut3_cor[[6]]))

plot(unlist(enet_mut3_cor[[2]]), unlist(enet_mut3_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut3_cor[[6]]))

# add codon72.npro
enet_mut4_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro"), iterations = 10)

plot(unlist(enet_mut4_cor[[2]]), unlist(enet_mut4_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut4_cor[[6]]))

plot(unlist(enet_mut4_cor[[2]]), unlist(enet_mut4_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut4_cor[[6]]))



# add splice.delins.snv
enet_mut5_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                          iterations = 10)

plot(unlist(enet_mut5_cor[[2]]), unlist(enet_mut5_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut5_cor[[6]]))

plot(unlist(enet_mut5_cor[[2]]), unlist(enet_mut5_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut5_cor[[6]]))


# add protein.codon.num
enet_mut6_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis","age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                          iterations = 10)

plot(unlist(enet_mut6_cor[[2]]), unlist(enet_mut6_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut6_cor[[6]]))

plot(unlist(enet_mut6_cor[[2]]), unlist(enet_mut6_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut6_cor[[6]]))


# add mdm2.nG
enet_mut7_cor <- predictAll(model_name = 'enet', 
                          data = full_data_cor,
                          subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                          selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                          iterations = 10)

plot(unlist(enet_mut7_cor[[2]]), unlist(enet_mut7_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
     mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut7_cor[[6]]))

plot(unlist(enet_mut7_cor[[2]]), unlist(enet_mut7_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
     mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut7_cor[[6]]))



####################################################################
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
# Run with #gdna.codon, protein.codon.num, mdm2.nG, protein.codon.change
# add mdm2.nG
enet_mdm2.nG_cor <- predictAll(model_name = 'enet', 
                             data = full_data_cor,
                             subset <- c("age_diagnosis", "age_sample_collection", "gender", "mdm2.nG"), 
                             selected_features = c("gender", "mdm2.nG"), iterations = 10)

plot(unlist(enet_mdm2.nG_cor[[2]]), unlist(enet_mdm2.nG_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mdm2.nG_cor[[6]]))

plot(unlist(enet_mdm2.nG_cor[[2]]), unlist(enet_mdm2.nG_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mdm2.nG_cor[[6]]))

# add gdna.codon
enet_gdna.codon_cor <- predictAll(model_name = 'enet', 
                                data = full_data_cor,
                                subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon"), 
                                selected_features = c("gender", "gdna.codon"), iterations = 10)

plot(unlist(enet_gdna.codon_cor[[2]]), unlist(enet_gdna.codon_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_gdna.codon_cor[[6]]))

plot(unlist(enet_gdna.codon_cor[[2]]), unlist(enet_gdna.codon_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_gdna.codon_cor[[6]]))

# add protein.codon.num
enet_protein.codon.num_cor <- predictAll(model_name = 'enet', 
                                       data = full_data_cor,
                                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "protein.codon.num"), 
                                       selected_features = c("gender", "protein.codon.num"), iterations = 10)

plot(unlist(enet_protein.codon.num_cor[[2]]), unlist(enet_protein.codon.num_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_protein.codon.num_cor[[6]]))

plot(unlist(enet_protein.codon.num_cor[[2]]), unlist(enet_protein.codon.num_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_protein.codon.num_cor[[6]]))

# add protein.codon.change
enet_protein.codon.change_cor <- predictAll(model_name = 'enet', 
                                          data = full_data_cor,
                                          subset <- c("age_diagnosis","age_sample_collection",  "gender", "protein.codon.change"), 
                                          selected_features = c("gender", "protein.codon.change"), 
                                          iterations = 10)

plot(unlist(enet_protein.codon.change_cor[[2]]), unlist(enet_protein.codon.change_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_protein.codon.change_cor[[6]]))

plot(unlist(enet_protein.codon.change_cor[[2]]), unlist(enet_protein.codon.change_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_protein.codon.change_cor[[6]]))

#########################
# combinations
# add gdna.codon and protein.codon.change
enet_gdna.codon_protein.codon.change_cor <- predictAll(model_name = 'enet', 
                                                     data = full_data_cor,
                                                     subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.codon", "protein.codon.change"), 
                                                     selected_features = c("gender", "gdna.codon", "protein.codon.change"), 
                                                     iterations = 10)


plot(unlist(enet_gdna.codon_protein.codon.change_cor[[2]]), unlist(enet_gdna.codon_protein.codon.change_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_gdna.codon_protein.codon.change_cor[[6]]))

plot(unlist(enet_gdna.codon_protein.codon.change_cor[[2]]), unlist(enet_gdna.codon_protein.codon.change_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_gdna.codon_protein.codon.change_cor[[6]]))


# add protein.codon.num
enet_add_protein.codon.num_cor <- predictAll(model_name = 'enet', 
                                           data = full_data_cor,
                                           subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num"), 
                                           selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num"), 
                                           iterations = 10)

plot(unlist(enet_add_protein.codon.num_cor[[2]]), unlist(enet_add_protein.codon.num_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_add_protein.codon.num_cor[[6]]))

plot(unlist(enet_add_protein.codon.num_cor[[2]]), unlist(enet_add_protein.codon.num_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_add_protein.codon.num_cor[[6]]))

# add mdm2.nG
enet_add_mdm2.nG_cor <- predictAll(model_name = 'enet', 
                                 data = full_data_cor,
                                 subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num", "mdm2.nG"), 
                                 selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num", "mdm2.nG"), 
                                 iterations = 10)

plot(unlist(enet_add_mdm2.nG_cor[[2]]), unlist(enet_add_mdm2.nG_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_add_mdm2.nG_cor[[6]]))

plot(unlist(enet_add_mdm2.nG_cor[[2]]), unlist(enet_add_mdm2.nG_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_add_mdm2.nG_cor[[6]]))

dev.off()
################################################################################################
# enet
pdf('/home/benbrew/Desktop/enet_methylation_reg.pdf')
# just methylation
# just methylation
enet_methyl_enet <- predictAll(model_name = 'enet', 
                           data = full_data_rf,
                           subset <- c("age_diagnosis", "age_sample_collection"), 
                           selected_features = NULL, iterations = 10)


plot(unlist(enet_methyl_enet[[2]]), unlist(enet_methyl_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = 'Just methylation')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_methyl_enet[[6]]))

plot(unlist(enet_methyl_enet[[2]]), unlist(enet_methyl_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = 'Just methylation')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_methyl_enet[[6]]))

# gender and gdna.base.change
enet_mut_enet <- predictAll(model_name = 'enet', 
                        data = full_data_rf,
                        subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change"), 
                        selected_features = c("gender", "gdna.base.change"), iterations = 10)

plot(unlist(enet_mut_enet[[2]]), unlist(enet_mut_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut_enet[[6]]))

plot(unlist(enet_mut_enet[[2]]), unlist(enet_mut_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut_enet[[6]]))



# add gdna.codon
enet_mut1_enet <- predictAll(model_name = 'enet', 
                         data = full_data_rf,
                         subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.base.change", "gdna.codon"), 
                         selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 10)

plot(unlist(enet_mut1_enet[[2]]), unlist(enet_mut1_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut1_enet[[6]]))

plot(unlist(enet_mut1_enet[[2]]), unlist(enet_mut1_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut1_enet[[6]]))


# add protein.codon.change
enet_mut2_enet <- predictAll(model_name = 'enet', 
                         data = full_data_rf,
                         subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                         selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                         iterations = 10)

plot(unlist(enet_mut2_enet[[2]]), unlist(enet_mut2_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut2_enet[[6]]))

plot(unlist(enet_mut2_enet[[2]]), unlist(enet_mut2_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut2_enet[[6]]))



# add gdna.exon.intron
enet_mut3_enet <- predictAll(model_name = 'enet', 
                         data = full_data_rf,
                         subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                     "gdna.exon.intron"), 
                         selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                               "gdna.exon.intron"), iterations = 10)

plot(unlist(enet_mut3_enet[[2]]), unlist(enet_mut3_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut3_enet[[6]]))

plot(unlist(enet_mut3_enet[[2]]), unlist(enet_mut3_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut3_enet[[6]]))

# add codon72.npro
enet_mut4_enet <- predictAll(model_name = 'enet', 
                         data = full_data_rf,
                         subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                     "gdna.exon.intron", "codon72.npro"), 
                         selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                               "gdna.exon.intron", "codon72.npro"), iterations = 10)

plot(unlist(enet_mut4_enet[[2]]), unlist(enet_mut4_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut4_enet[[6]]))

plot(unlist(enet_mut4_enet[[2]]), unlist(enet_mut4_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut4_enet[[6]]))



# add splice.delins.snv
enet_mut5_enet <- predictAll(model_name = 'enet', 
                         data = full_data_rf,
                         subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                     "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                         selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                               "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                         iterations = 10)

plot(unlist(enet_mut5_enet[[2]]), unlist(enet_mut5_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut5_enet[[6]]))

plot(unlist(enet_mut5_enet[[2]]), unlist(enet_mut5_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut5_enet[[6]]))


# add protein.codon.num
enet_mut6_enet <- predictAll(model_name = 'enet', 
                         data = full_data_rf,
                         subset <- c("age_diagnosis","age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                     "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                         selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                               "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                         iterations = 10)

plot(unlist(enet_mut6_enet[[2]]), unlist(enet_mut6_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut6_enet[[6]]))

plot(unlist(enet_mut6_enet[[2]]), unlist(enet_mut6_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut6_enet[[6]]))


# add mdm2.nG
enet_mut7_enet <- predictAll(model_name = 'enet', 
                         data = full_data_rf,
                         subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                     "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                         selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                               "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                         iterations = 10)

plot(unlist(enet_mut7_enet[[2]]), unlist(enet_mut7_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
     mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut7_enet[[6]]))

plot(unlist(enet_mut7_enet[[2]]), unlist(enet_mut7_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
     mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mut7_enet[[6]]))



####################################################################
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
# Run with #gdna.codon, protein.codon.num, mdm2.nG, protein.codon.change
# add mdm2.nG
enet_mdm2.nG_enet <- predictAll(model_name = 'enet', 
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "age_sample_collection", "gender", "mdm2.nG"), 
                            selected_features = c("gender", "mdm2.nG"), iterations = 10)

plot(unlist(enet_mdm2.nG_enet[[2]]), unlist(enet_mdm2.nG_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mdm2.nG_enet[[6]]))

plot(unlist(enet_mdm2.nG_enet[[2]]), unlist(enet_mdm2.nG_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_mdm2.nG_enet[[6]]))

# add gdna.codon
enet_gdna.codon_enet <- predictAll(model_name = 'enet', 
                               data = full_data_rf,
                               subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon"), 
                               selected_features = c("gender", "gdna.codon"), iterations = 10)

plot(unlist(enet_gdna.codon_enet[[2]]), unlist(enet_gdna.codon_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_gdna.codon_enet[[6]]))

plot(unlist(enet_gdna.codon_enet[[2]]), unlist(enet_gdna.codon_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_gdna.codon_enet[[6]]))

# add protein.codon.num
enet_protein.codon.num_enet <- predictAll(model_name = 'enet', 
                                      data = full_data_rf,
                                      subset <- c("age_diagnosis", "age_sample_collection", "gender", "protein.codon.num"), 
                                      selected_features = c("gender", "protein.codon.num"), iterations = 10)

plot(unlist(enet_protein.codon.num_enet[[2]]), unlist(enet_protein.codon.num_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_protein.codon.num_enet[[6]]))

plot(unlist(enet_protein.codon.num_enet[[2]]), unlist(enet_protein.codon.num_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_protein.codon.num_enet[[6]]))

# add protein.codon.change
enet_protein.codon.change_enet <- predictAll(model_name = 'enet', 
                                         data = full_data_rf,
                                         subset <- c("age_diagnosis","age_sample_collection",  "gender", "protein.codon.change"), 
                                         selected_features = c("gender", "protein.codon.change"), 
                                         iterations = 10)

plot(unlist(enet_protein.codon.change_enet[[2]]), unlist(enet_protein.codon.change_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_protein.codon.change_enet[[6]]))

plot(unlist(enet_protein.codon.change_enet[[2]]), unlist(enet_protein.codon.change_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_protein.codon.change_enet[[6]]))

#########################
# combinations
# add gdna.codon and protein.codon.change
enet_gdna.codon_protein.codon.change_enet <- predictAll(model_name = 'enet', 
                                                    data = full_data_rf,
                                                    subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.codon", "protein.codon.change"), 
                                                    selected_features = c("gender", "gdna.codon", "protein.codon.change"), 
                                                    iterations = 10)


plot(unlist(enet_gdna.codon_protein.codon.change_enet[[2]]), unlist(enet_gdna.codon_protein.codon.change_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_gdna.codon_protein.codon.change_enet[[6]]))

plot(unlist(enet_gdna.codon_protein.codon.change_enet[[2]]), unlist(enet_gdna.codon_protein.codon.change_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_gdna.codon_protein.codon.change_enet[[6]]))


# add protein.codon.num
enet_add_protein.codon.num_enet <- predictAll(model_name = 'enet', 
                                          data = full_data_rf,
                                          subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num"), 
                                          selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num"), 
                                          iterations = 10)

plot(unlist(enet_add_protein.codon.num_enet[[2]]), unlist(enet_add_protein.codon.num_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_add_protein.codon.num_enet[[6]]))

plot(unlist(enet_add_protein.codon.num_enet[[2]]), unlist(enet_add_protein.codon.num_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_add_protein.codon.num_enet[[6]]))

# add mdm2.nG
enet_add_mdm2.nG_enet <- predictAll(model_name = 'enet', 
                                data = full_data_rf,
                                subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num", "mdm2.nG"), 
                                selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num", "mdm2.nG"), 
                                iterations = 10)

plot(unlist(enet_add_mdm2.nG_enet[[2]]), unlist(enet_add_mdm2.nG_enet[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_add_mdm2.nG_enet[[6]]))

plot(unlist(enet_add_mdm2.nG_enet[[2]]), unlist(enet_add_mdm2.nG_enet[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', enet_add_mdm2.nG_enet[[6]]))

dev.off()
# # combine all tests 
# 
# enet_all <- rbind (
#   append('just_methyl', c(mean(unlist(enet_methyl[[1]])), enet_methyl[[6]], 'all_data')),
#   append('gender_and_gdna.base.change', c(mean(unlist(enet_mut[[1]])), enet_mut[[6]], 'all_data')),
#   append('gdna.codon', c(mean(unlist(enet_mut1[[1]])), enet_mut1[[6]], 'all_data')),
#   append('protein.codon.change', c(mean(unlist(enet_mut2[[1]])), enet_mut2[[6]], 'all_data')),
#   append('gdna.exon.intron', c(mean(unlist(enet_mut3[[1]])), enet_mut3[[6]], 'all_data')),
#   append('codon72.npro', c(mean(unlist(enet_mut4[[1]])), enet_mut4[[6]], 'all_data')),
#   append('splice.delins.snv', c(mean(unlist(enet_mut5[[1]])), enet_mut5[[6]], 'all_data')),
#   append('protein.codon.num', c(mean(unlist(enet_mut6[[1]])), enet_mut6[[6]], 'all_data')),
#   append('mdm2.nG', c(mean(unlist(enet_mut7[[1]])), enet_mut7[[6]], 'all_data')),
#   
#   append('just_methyl', c(mean(unlist(enet_methyl_cor[[1]])), enet_methyl_cor[[6]], 'corr')),
#   append('gender_and_gdna.base.change', c(mean(unlist(enet_mut_cor[[1]])), enet_mut_cor[[6]], 'corr')),
#   append('gdna.codon', c(mean(unlist(enet_mut1_cor[[1]])), enet_mut1_cor[[6]], 'corr')),
#   append('protein.codon.change', c(mean(unlist(enet_mut2_cor[[1]])), enet_mut2_cor[[6]], 'corr')),
#   append('gdna.exon.intron', c(mean(unlist(enet_mut3_cor[[1]])), enet_mut3_cor[[6]], 'corr')),
#   append('codon72.npro', c(mean(unlist(enet_mut4_cor[[1]])), enet_mut4_cor[[6]], 'corr')),
#   append('splice.delins.snv', c(mean(unlist(enet_mut5_cor[[1]])), enet_mut5_cor[[6]], 'corr')),
#   append('protein.codon.num', c(mean(unlist(enet_mut6_cor[[1]])), enet_mut6_cor[[6]], 'corr')),
#   append('mdm2.nG', c(mean(unlist(enet_mut7_cor[[1]])), enet_mut7_cor[[6]], 'corr')),
#   
#   append('just_methyl', c(mean(unlist(enet_methyl_reg[[1]])), enet_methyl_reg[[6]], 'recursive')),
#   append('gender_and_gdna.base.change', c(mean(unlist(enet_mut_reg[[1]])), enet_mut_reg[[6]], 'recursive')),
#   append('gdna.codon', c(mean(unlist(enet_mut1_reg[[1]])), enet_mut1_reg[[6]], 'recursive')),
#   append('protein.codon.change', c(mean(unlist(enet_mut2_reg[[1]])), enet_mut2_reg[[6]], 'recursive')),
#   append('gdna.exon.intron', c(mean(unlist(enet_mut3_reg[[1]])), enet_mut3_reg[[6]], 'recursive')),
#   append('codon72.npro', c(mean(unlist(enet_mut4_reg[[1]])), enet_mut4_reg[[6]], 'recursive')),
#   append('splice.delins.snv', c(mean(unlist(enet_mut5_reg[[1]])), enet_mut5_reg[[6]], 'recursive')),
#   append('protein.codon.num', c(mean(unlist(enet_mut6_reg[[1]])), enet_mut6_reg[[6]], 'recursive')),
#   append('mdm2.nG', c(mean(unlist(enet_mut7_reg[[1]])), enet_mut7_reg[[6]], 'recursive'))
# )




##################################################################################################################3
# Lasso
##################################################################################################################3
pdf('/home/benbrew/Desktop/lasso_methylation.pdf')

# just methylation
lasso_methyl <- predictAll(model_name = 'lasso', 
                          data = full_data,
                          subset <- c("age_diagnosis", "age_sample_collection"), 
                          selected_features = NULL, iterations = 10)


plot(unlist(lasso_methyl[[2]]), unlist(lasso_methyl[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = 'Just methylation')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_methyl[[6]]))

plot(unlist(lasso_methyl[[2]]), unlist(lasso_methyl[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = 'Just methylation')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_methyl[[6]]))

# gender and gdna.base.change
lasso_mut <- predictAll(model_name = 'lasso', 
                       data = full_data,
                       subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change"), 
                       selected_features = c("gender", "gdna.base.change"), iterations = 10)

plot(unlist(lasso_mut[[2]]), unlist(lasso_mut[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut[[6]]))

plot(unlist(lasso_mut[[2]]), unlist(lasso_mut[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut[[6]]))



# add gdna.codon
lasso_mut1 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.base.change", "gdna.codon"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 10)

plot(unlist(lasso_mut1[[2]]), unlist(lasso_mut1[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut1[[6]]))

plot(unlist(lasso_mut1[[2]]), unlist(lasso_mut1[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut1[[6]]))


# add protein.codon.change
lasso_mut2 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                        iterations = 10)

plot(unlist(lasso_mut2[[2]]), unlist(lasso_mut2[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut2[[6]]))

plot(unlist(lasso_mut2[[2]]), unlist(lasso_mut2[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut2[[6]]))



# add gdna.exon.intron
lasso_mut3 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                    "gdna.exon.intron"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron"), iterations = 10)

plot(unlist(lasso_mut3[[2]]), unlist(lasso_mut3[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut3[[6]]))

plot(unlist(lasso_mut3[[2]]), unlist(lasso_mut3[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut3[[6]]))

# add codon72.npro
lasso_mut4 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                    "gdna.exon.intron", "codon72.npro"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron", "codon72.npro"), iterations = 10)

plot(unlist(lasso_mut4[[2]]), unlist(lasso_mut4[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut4[[6]]))

plot(unlist(lasso_mut4[[2]]), unlist(lasso_mut4[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut4[[6]]))



# add splice.delins.snv
lasso_mut5 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                    "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                        iterations = 10)

plot(unlist(lasso_mut5[[2]]), unlist(lasso_mut5[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut5[[6]]))

plot(unlist(lasso_mut5[[2]]), unlist(lasso_mut5[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut5[[6]]))


# add protein.codon.num
lasso_mut6 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis","age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                    "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                        iterations = 10)

plot(unlist(lasso_mut6[[2]]), unlist(lasso_mut6[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut6[[6]]))

plot(unlist(lasso_mut6[[2]]), unlist(lasso_mut6[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut6[[6]]))


# add mdm2.nG
lasso_mut7 <- predictAll(model_name = 'lasso', 
                        data = full_data,
                        subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                    "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                        selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                        iterations = 10)

plot(unlist(lasso_mut7[[2]]), unlist(lasso_mut7[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
     mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut7[[6]]))

plot(unlist(lasso_mut7[[2]]), unlist(lasso_mut7[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
     mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut7[[6]]))



####################################################################
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
# Run with #gdna.codon, protein.codon.num, mdm2.nG, protein.codon.change
# add mdm2.nG
lasso_mdm2.nG <- predictAll(model_name = 'lasso', 
                           data = full_data,
                           subset <- c("age_diagnosis", "age_sample_collection", "gender", "mdm2.nG"), 
                           selected_features = c("gender", "mdm2.nG"), iterations = 10)

plot(unlist(lasso_mdm2.nG[[2]]), unlist(lasso_mdm2.nG[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mdm2.nG[[6]]))

plot(unlist(lasso_mdm2.nG[[2]]), unlist(lasso_mdm2.nG[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mdm2.nG[[6]]))

# add gdna.codon
lasso_gdna.codon <- predictAll(model_name = 'lasso', 
                              data = full_data,
                              subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon"), 
                              selected_features = c("gender", "gdna.codon"), iterations = 10)

plot(unlist(lasso_gdna.codon[[2]]), unlist(lasso_gdna.codon[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_gdna.codon[[6]]))

plot(unlist(lasso_gdna.codon[[2]]), unlist(lasso_gdna.codon[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_gdna.codon[[6]]))

# add protein.codon.num
lasso_protein.codon.num <- predictAll(model_name = 'lasso', 
                                     data = full_data,
                                     subset <- c("age_diagnosis", "age_sample_collection", "gender", "protein.codon.num"), 
                                     selected_features = c("gender", "protein.codon.num"), iterations = 10)

plot(unlist(lasso_protein.codon.num[[2]]), unlist(lasso_protein.codon.num[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_protein.codon.num[[6]]))

plot(unlist(lasso_protein.codon.num[[2]]), unlist(lasso_protein.codon.num[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_protein.codon.num[[6]]))

# add protein.codon.change
lasso_protein.codon.change <- predictAll(model_name = 'lasso', 
                                        data = full_data,
                                        subset <- c("age_diagnosis","age_sample_collection",  "gender", "protein.codon.change"), 
                                        selected_features = c("gender", "protein.codon.change"), 
                                        iterations = 10)

plot(unlist(lasso_protein.codon.change[[2]]), unlist(lasso_protein.codon.change[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_protein.codon.change[[6]]))

plot(unlist(lasso_protein.codon.change[[2]]), unlist(lasso_protein.codon.change[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_protein.codon.change[[6]]))

#########################
# combinations
# add gdna.codon and protein.codon.change
lasso_gdna.codon_protein.codon.change <- predictAll(model_name = 'lasso', 
                                                   data = full_data,
                                                   subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.codon", "protein.codon.change"), 
                                                   selected_features = c("gender", "gdna.codon", "protein.codon.change"), 
                                                   iterations = 10)


plot(unlist(lasso_gdna.codon_protein.codon.change[[2]]), unlist(lasso_gdna.codon_protein.codon.change[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_gdna.codon_protein.codon.change[[6]]))

plot(unlist(lasso_gdna.codon_protein.codon.change[[2]]), unlist(lasso_gdna.codon_protein.codon.change[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_gdna.codon_protein.codon.change[[6]]))


# add protein.codon.num
lasso_add_protein.codon.num <- predictAll(model_name = 'lasso', 
                                         data = full_data,
                                         subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num"), 
                                         selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num"), 
                                         iterations = 10)

plot(unlist(lasso_add_protein.codon.num[[2]]), unlist(lasso_add_protein.codon.num[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_add_protein.codon.num[[6]]))

plot(unlist(lasso_add_protein.codon.num[[2]]), unlist(lasso_add_protein.codon.num[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_add_protein.codon.num[[6]]))

# add mdm2.nG
lasso_add_mdm2.nG <- predictAll(model_name = 'lasso', 
                               data = full_data,
                               subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num", "mdm2.nG"), 
                               selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num", "mdm2.nG"), 
                               iterations = 10)

plot(unlist(lasso_add_mdm2.nG[[2]]), unlist(lasso_add_mdm2.nG[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_add_mdm2.nG[[6]]))

plot(unlist(lasso_add_mdm2.nG[[2]]), unlist(lasso_add_mdm2.nG[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_add_mdm2.nG[[6]]))

dev.off()
########################################################################################
# Corelation data
pdf('/home/benbrew/Desktop/lasso_methylation_cor.pdf')

# just methylation
lasso_methyl_cor <- predictAll(model_name = 'lasso', 
                              data = full_data_cor,
                              subset <- c("age_diagnosis", "age_sample_collection"), 
                              selected_features = NULL, iterations = 10)


plot(unlist(lasso_methyl_cor[[2]]), unlist(lasso_methyl_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = 'Just methylation')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_methyl_cor[[6]]))

plot(unlist(lasso_methyl_cor[[2]]), unlist(lasso_methyl_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = 'Just methylation')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_methyl_cor[[6]]))

# gender and gdna.base.change
lasso_mut_cor <- predictAll(model_name = 'lasso', 
                           data = full_data_cor,
                           subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change"), 
                           selected_features = c("gender", "gdna.base.change"), iterations = 10)

plot(unlist(lasso_mut_cor[[2]]), unlist(lasso_mut_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut_cor[[6]]))

plot(unlist(lasso_mut_cor[[2]]), unlist(lasso_mut_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut_cor[[6]]))



# add gdna.codon
lasso_mut1_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.base.change", "gdna.codon"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 10)

plot(unlist(lasso_mut1_cor[[2]]), unlist(lasso_mut1_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut1_cor[[6]]))

plot(unlist(lasso_mut1_cor[[2]]), unlist(lasso_mut1_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut1_cor[[6]]))


# add protein.codon.change
lasso_mut2_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                            iterations = 10)

plot(unlist(lasso_mut2_cor[[2]]), unlist(lasso_mut2_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut2_cor[[6]]))

plot(unlist(lasso_mut2_cor[[2]]), unlist(lasso_mut2_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut2_cor[[6]]))



# add gdna.exon.intron
lasso_mut3_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron"), iterations = 10)

plot(unlist(lasso_mut3_cor[[2]]), unlist(lasso_mut3_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut3_cor[[6]]))

plot(unlist(lasso_mut3_cor[[2]]), unlist(lasso_mut3_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut3_cor[[6]]))

# add codon72.npro
lasso_mut4_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro"), iterations = 10)

plot(unlist(lasso_mut4_cor[[2]]), unlist(lasso_mut4_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut4_cor[[6]]))

plot(unlist(lasso_mut4_cor[[2]]), unlist(lasso_mut4_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut4_cor[[6]]))



# add splice.delins.snv
lasso_mut5_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                            iterations = 10)

plot(unlist(lasso_mut5_cor[[2]]), unlist(lasso_mut5_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut5_cor[[6]]))

plot(unlist(lasso_mut5_cor[[2]]), unlist(lasso_mut5_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut5_cor[[6]]))


# add protein.codon.num
lasso_mut6_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis","age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                            iterations = 10)

plot(unlist(lasso_mut6_cor[[2]]), unlist(lasso_mut6_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut6_cor[[6]]))

plot(unlist(lasso_mut6_cor[[2]]), unlist(lasso_mut6_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut6_cor[[6]]))


# add mdm2.nG
lasso_mut7_cor <- predictAll(model_name = 'lasso', 
                            data = full_data_cor,
                            subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                            iterations = 10)

plot(unlist(lasso_mut7_cor[[2]]), unlist(lasso_mut7_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
     mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut7_cor[[6]]))

plot(unlist(lasso_mut7_cor[[2]]), unlist(lasso_mut7_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
     mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut7_cor[[6]]))



####################################################################
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
# Run with #gdna.codon, protein.codon.num, mdm2.nG, protein.codon.change
# add mdm2.nG
lasso_mdm2.nG_cor <- predictAll(model_name = 'lasso', 
                               data = full_data_cor,
                               subset <- c("age_diagnosis", "age_sample_collection", "gender", "mdm2.nG"), 
                               selected_features = c("gender", "mdm2.nG"), iterations = 10)

plot(unlist(lasso_mdm2.nG_cor[[2]]), unlist(lasso_mdm2.nG_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mdm2.nG_cor[[6]]))

plot(unlist(lasso_mdm2.nG_cor[[2]]), unlist(lasso_mdm2.nG_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mdm2.nG_cor[[6]]))

# add gdna.codon
lasso_gdna.codon_cor <- predictAll(model_name = 'lasso', 
                                  data = full_data_cor,
                                  subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon"), 
                                  selected_features = c("gender", "gdna.codon"), iterations = 10)

plot(unlist(lasso_gdna.codon_cor[[2]]), unlist(lasso_gdna.codon_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_gdna.codon_cor[[6]]))

plot(unlist(lasso_gdna.codon_cor[[2]]), unlist(lasso_gdna.codon_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_gdna.codon_cor[[6]]))

# add protein.codon.num
lasso_protein.codon.num_cor <- predictAll(model_name = 'lasso', 
                                         data = full_data_cor,
                                         subset <- c("age_diagnosis", "age_sample_collection", "gender", "protein.codon.num"), 
                                         selected_features = c("gender", "protein.codon.num"), iterations = 10)

plot(unlist(lasso_protein.codon.num_cor[[2]]), unlist(lasso_protein.codon.num_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_protein.codon.num_cor[[6]]))

plot(unlist(lasso_protein.codon.num_cor[[2]]), unlist(lasso_protein.codon.num_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_protein.codon.num_cor[[6]]))

# add protein.codon.change
lasso_protein.codon.change_cor <- predictAll(model_name = 'lasso', 
                                            data = full_data_cor,
                                            subset <- c("age_diagnosis","age_sample_collection",  "gender", "protein.codon.change"), 
                                            selected_features = c("gender", "protein.codon.change"), 
                                            iterations = 10)

plot(unlist(lasso_protein.codon.change_cor[[2]]), unlist(lasso_protein.codon.change_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_protein.codon.change_cor[[6]]))

plot(unlist(lasso_protein.codon.change_cor[[2]]), unlist(lasso_protein.codon.change_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_protein.codon.change_cor[[6]]))

#########################
# combinations
# add gdna.codon and protein.codon.change
lasso_gdna.codon_protein.codon.change_cor <- predictAll(model_name = 'lasso', 
                                                       data = full_data_cor,
                                                       subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.codon", "protein.codon.change"), 
                                                       selected_features = c("gender", "gdna.codon", "protein.codon.change"), 
                                                       iterations = 10)


plot(unlist(lasso_gdna.codon_protein.codon.change_cor[[2]]), unlist(lasso_gdna.codon_protein.codon.change_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_gdna.codon_protein.codon.change_cor[[6]]))

plot(unlist(lasso_gdna.codon_protein.codon.change_cor[[2]]), unlist(lasso_gdna.codon_protein.codon.change_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_gdna.codon_protein.codon.change_cor[[6]]))


# add protein.codon.num
lasso_add_protein.codon.num_cor <- predictAll(model_name = 'lasso', 
                                             data = full_data_cor,
                                             subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num"), 
                                             selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num"), 
                                             iterations = 10)

plot(unlist(lasso_add_protein.codon.num_cor[[2]]), unlist(lasso_add_protein.codon.num_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_add_protein.codon.num_cor[[6]]))

plot(unlist(lasso_add_protein.codon.num_cor[[2]]), unlist(lasso_add_protein.codon.num_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_add_protein.codon.num_cor[[6]]))

# add mdm2.nG
lasso_add_mdm2.nG_cor <- predictAll(model_name = 'lasso', 
                                   data = full_data_cor,
                                   subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num", "mdm2.nG"), 
                                   selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num", "mdm2.nG"), 
                                   iterations = 10)

plot(unlist(lasso_add_mdm2.nG_cor[[2]]), unlist(lasso_add_mdm2.nG_cor[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_add_mdm2.nG_cor[[6]]))

plot(unlist(lasso_add_mdm2.nG_cor[[2]]), unlist(lasso_add_mdm2.nG_cor[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_add_mdm2.nG_cor[[6]]))

dev.off()
################################################################################################
# lasso
pdf('/home/benbrew/Desktop/lasso_methylation_reg.pdf')
# just methylation
# just methylation
lasso_methyl_lasso <- predictAll(model_name = 'lasso', 
                               data = full_data_rf,
                               subset <- c("age_diagnosis", "age_sample_collection"), 
                               selected_features = NULL, iterations = 10)


plot(unlist(lasso_methyl_lasso[[2]]), unlist(lasso_methyl_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = 'Just methylation')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_methyl_lasso[[6]]))

plot(unlist(lasso_methyl_lasso[[2]]), unlist(lasso_methyl_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = 'Just methylation')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_methyl_lasso[[6]]))

# gender and gdna.base.change
lasso_mut_lasso <- predictAll(model_name = 'lasso', 
                            data = full_data_rf,
                            subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change"), 
                            selected_features = c("gender", "gdna.base.change"), iterations = 10)

plot(unlist(lasso_mut_lasso[[2]]), unlist(lasso_mut_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut_lasso[[6]]))

plot(unlist(lasso_mut_lasso[[2]]), unlist(lasso_mut_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut_lasso[[6]]))



# add gdna.codon
lasso_mut1_lasso <- predictAll(model_name = 'lasso', 
                             data = full_data_rf,
                             subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.base.change", "gdna.codon"), 
                             selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 10)

plot(unlist(lasso_mut1_lasso[[2]]), unlist(lasso_mut1_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut1_lasso[[6]]))

plot(unlist(lasso_mut1_lasso[[2]]), unlist(lasso_mut1_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut1_lasso[[6]]))


# add protein.codon.change
lasso_mut2_lasso <- predictAll(model_name = 'lasso', 
                             data = full_data_rf,
                             subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                             selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                             iterations = 10)

plot(unlist(lasso_mut2_lasso[[2]]), unlist(lasso_mut2_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut2_lasso[[6]]))

plot(unlist(lasso_mut2_lasso[[2]]), unlist(lasso_mut2_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut2_lasso[[6]]))



# add gdna.exon.intron
lasso_mut3_lasso <- predictAll(model_name = 'lasso', 
                             data = full_data_rf,
                             subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                         "gdna.exon.intron"), 
                             selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                   "gdna.exon.intron"), iterations = 10)

plot(unlist(lasso_mut3_lasso[[2]]), unlist(lasso_mut3_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut3_lasso[[6]]))

plot(unlist(lasso_mut3_lasso[[2]]), unlist(lasso_mut3_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+gdna.exon.intron')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut3_lasso[[6]]))

# add codon72.npro
lasso_mut4_lasso <- predictAll(model_name = 'lasso', 
                             data = full_data_rf,
                             subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                         "gdna.exon.intron", "codon72.npro"), 
                             selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                   "gdna.exon.intron", "codon72.npro"), iterations = 10)

plot(unlist(lasso_mut4_lasso[[2]]), unlist(lasso_mut4_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut4_lasso[[6]]))

plot(unlist(lasso_mut4_lasso[[2]]), unlist(lasso_mut4_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut4_lasso[[6]]))



# add splice.delins.snv
lasso_mut5_lasso <- predictAll(model_name = 'lasso', 
                             data = full_data_rf,
                             subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                         "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                             selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                   "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                             iterations = 10)

plot(unlist(lasso_mut5_lasso[[2]]), unlist(lasso_mut5_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut5_lasso[[6]]))

plot(unlist(lasso_mut5_lasso[[2]]), unlist(lasso_mut5_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut5_lasso[[6]]))


# add protein.codon.num
lasso_mut6_lasso <- predictAll(model_name = 'lasso', 
                             data = full_data_rf,
                             subset <- c("age_diagnosis","age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                         "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                             selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                   "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                             iterations = 10)

plot(unlist(lasso_mut6_lasso[[2]]), unlist(lasso_mut6_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut6_lasso[[6]]))

plot(unlist(lasso_mut6_lasso[[2]]), unlist(lasso_mut6_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut6_lasso[[6]]))


# add mdm2.nG
lasso_mut7_lasso <- predictAll(model_name = 'lasso', 
                             data = full_data_rf,
                             subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                         "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                             selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                   "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                             iterations = 10)

plot(unlist(lasso_mut7_lasso[[2]]), unlist(lasso_mut7_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
     mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut7_lasso[[6]]))

plot(unlist(lasso_mut7_lasso[[2]]), unlist(lasso_mut7_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.base.change+gdna.codon+protein.codon.change+
     gdna.exon.intron+codon72.npro+splice.delins.snv+protein.codon.num+
     mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mut7_lasso[[6]]))



####################################################################
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
# Run with #gdna.codon, protein.codon.num, mdm2.nG, protein.codon.change
# add mdm2.nG
lasso_mdm2.nG_lasso <- predictAll(model_name = 'lasso', 
                                data = full_data_rf,
                                subset <- c("age_diagnosis", "age_sample_collection", "gender", "mdm2.nG"), 
                                selected_features = c("gender", "mdm2.nG"), iterations = 10)

plot(unlist(lasso_mdm2.nG_lasso[[2]]), unlist(lasso_mdm2.nG_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mdm2.nG_lasso[[6]]))

plot(unlist(lasso_mdm2.nG_lasso[[2]]), unlist(lasso_mdm2.nG_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_mdm2.nG_lasso[[6]]))

# add gdna.codon
lasso_gdna.codon_lasso <- predictAll(model_name = 'lasso', 
                                   data = full_data_rf,
                                   subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon"), 
                                   selected_features = c("gender", "gdna.codon"), iterations = 10)

plot(unlist(lasso_gdna.codon_lasso[[2]]), unlist(lasso_gdna.codon_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_gdna.codon_lasso[[6]]))

plot(unlist(lasso_gdna.codon_lasso[[2]]), unlist(lasso_gdna.codon_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_gdna.codon_lasso[[6]]))

# add protein.codon.num
lasso_protein.codon.num_lasso <- predictAll(model_name = 'lasso', 
                                          data = full_data_rf,
                                          subset <- c("age_diagnosis", "age_sample_collection", "gender", "protein.codon.num"), 
                                          selected_features = c("gender", "protein.codon.num"), iterations = 10)

plot(unlist(lasso_protein.codon.num_lasso[[2]]), unlist(lasso_protein.codon.num_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_protein.codon.num_lasso[[6]]))

plot(unlist(lasso_protein.codon.num_lasso[[2]]), unlist(lasso_protein.codon.num_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_protein.codon.num_lasso[[6]]))

# add protein.codon.change
lasso_protein.codon.change_lasso <- predictAll(model_name = 'lasso', 
                                             data = full_data_rf,
                                             subset <- c("age_diagnosis","age_sample_collection",  "gender", "protein.codon.change"), 
                                             selected_features = c("gender", "protein.codon.change"), 
                                             iterations = 10)

plot(unlist(lasso_protein.codon.change_lasso[[2]]), unlist(lasso_protein.codon.change_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_protein.codon.change_lasso[[6]]))

plot(unlist(lasso_protein.codon.change_lasso[[2]]), unlist(lasso_protein.codon.change_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_protein.codon.change_lasso[[6]]))

#########################
# combinations
# add gdna.codon and protein.codon.change
lasso_gdna.codon_protein.codon.change_lasso <- predictAll(model_name = 'lasso', 
                                                        data = full_data_rf,
                                                        subset <- c("age_diagnosis","age_sample_collection",  "gender", "gdna.codon", "protein.codon.change"), 
                                                        selected_features = c("gender", "gdna.codon", "protein.codon.change"), 
                                                        iterations = 10)


plot(unlist(lasso_gdna.codon_protein.codon.change_lasso[[2]]), unlist(lasso_gdna.codon_protein.codon.change_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_gdna.codon_protein.codon.change_lasso[[6]]))

plot(unlist(lasso_gdna.codon_protein.codon.change_lasso[[2]]), unlist(lasso_gdna.codon_protein.codon.change_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_gdna.codon_protein.codon.change_lasso[[6]]))


# add protein.codon.num
lasso_add_protein.codon.num_lasso <- predictAll(model_name = 'lasso', 
                                              data = full_data_rf,
                                              subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num"), 
                                              selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num"), 
                                              iterations = 10)

plot(unlist(lasso_add_protein.codon.num_lasso[[2]]), unlist(lasso_add_protein.codon.num_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_add_protein.codon.num_lasso[[6]]))

plot(unlist(lasso_add_protein.codon.num_lasso[[2]]), unlist(lasso_add_protein.codon.num_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_add_protein.codon.num_lasso[[6]]))

# add mdm2.nG
lasso_add_mdm2.nG_lasso <- predictAll(model_name = 'lasso', 
                                    data = full_data_rf,
                                    subset <- c("age_diagnosis", "age_sample_collection", "gender", "gdna.codon", "protein.codon.change","protein.codon.num", "mdm2.nG"), 
                                    selected_features = c("gender", "gdna.codon", "protein.codon.change", "protein.codon.num", "mdm2.nG"), 
                                    iterations = 10)

plot(unlist(lasso_add_mdm2.nG_lasso[[2]]), unlist(lasso_add_mdm2.nG_lasso[[5]]), xlab = 'Predictions', ylab = 'Age of Diagnosis',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_add_mdm2.nG_lasso[[6]]))

plot(unlist(lasso_add_mdm2.nG_lasso[[2]]), unlist(lasso_add_mdm2.nG_lasso[[7]]), xlab = 'Predictions', ylab = 'Age of Sample Collection',
     main = '+gender+gdna.codon+protein.codon.change+protein.codon.num+mdm2.nG')
abline(0,1)
legend("bottomright", legend = paste0('# obs = ', lasso_add_mdm2.nG_lasso[[6]]))

dev.off()
# # combine all tests 
# 
# lasso_all <- rbind (
#   append('just_methyl', c(mean(unlist(lasso_methyl[[1]])), lasso_methyl[[6]], 'all_data')),
#   append('gender_and_gdna.base.change', c(mean(unlist(lasso_mut[[1]])), lasso_mut[[6]], 'all_data')),
#   append('gdna.codon', c(mean(unlist(lasso_mut1[[1]])), lasso_mut1[[6]], 'all_data')),
#   append('protein.codon.change', c(mean(unlist(lasso_mut2[[1]])), lasso_mut2[[6]], 'all_data')),
#   append('gdna.exon.intron', c(mean(unlist(lasso_mut3[[1]])), lasso_mut3[[6]], 'all_data')),
#   append('codon72.npro', c(mean(unlist(lasso_mut4[[1]])), lasso_mut4[[6]], 'all_data')),
#   append('splice.delins.snv', c(mean(unlist(lasso_mut5[[1]])), lasso_mut5[[6]], 'all_data')),
#   append('protein.codon.num', c(mean(unlist(lasso_mut6[[1]])), lasso_mut6[[6]], 'all_data')),
#   append('mdm2.nG', c(mean(unlist(lasso_mut7[[1]])), lasso_mut7[[6]], 'all_data')),
#   
#   append('just_methyl', c(mean(unlist(lasso_methyl_cor[[1]])), lasso_methyl_cor[[6]], 'corr')),
#   append('gender_and_gdna.base.change', c(mean(unlist(lasso_mut_cor[[1]])), lasso_mut_cor[[6]], 'corr')),
#   append('gdna.codon', c(mean(unlist(lasso_mut1_cor[[1]])), lasso_mut1_cor[[6]], 'corr')),
#   append('protein.codon.change', c(mean(unlist(lasso_mut2_cor[[1]])), lasso_mut2_cor[[6]], 'corr')),
#   append('gdna.exon.intron', c(mean(unlist(lasso_mut3_cor[[1]])), lasso_mut3_cor[[6]], 'corr')),
#   append('codon72.npro', c(mean(unlist(lasso_mut4_cor[[1]])), lasso_mut4_cor[[6]], 'corr')),
#   append('splice.delins.snv', c(mean(unlist(lasso_mut5_cor[[1]])), lasso_mut5_cor[[6]], 'corr')),
#   append('protein.codon.num', c(mean(unlist(lasso_mut6_cor[[1]])), lasso_mut6_cor[[6]], 'corr')),
#   append('mdm2.nG', c(mean(unlist(lasso_mut7_cor[[1]])), lasso_mut7_cor[[6]], 'corr')),
#   
#   append('just_methyl', c(mean(unlist(lasso_methyl_reg[[1]])), lasso_methyl_reg[[6]], 'recursive')),
#   append('gender_and_gdna.base.change', c(mean(unlist(lasso_mut_reg[[1]])), lasso_mut_reg[[6]], 'recursive')),
#   append('gdna.codon', c(mean(unlist(lasso_mut1_reg[[1]])), lasso_mut1_reg[[6]], 'recursive')),
#   append('protein.codon.change', c(mean(unlist(lasso_mut2_reg[[1]])), lasso_mut2_reg[[6]], 'recursive')),
#   append('gdna.exon.intron', c(mean(unlist(lasso_mut3_reg[[1]])), lasso_mut3_reg[[6]], 'recursive')),
#   append('codon72.npro', c(mean(unlist(lasso_mut4_reg[[1]])), lasso_mut4_reg[[6]], 'recursive')),
#   append('splice.delins.snv', c(mean(unlist(lasso_mut5_reg[[1]])), lasso_mut5_reg[[6]], 'recursive')),
#   append('protein.codon.num', c(mean(unlist(lasso_mut6_reg[[1]])), lasso_mut6_reg[[6]], 'recursive')),
#   append('mdm2.nG', c(mean(unlist(lasso_mut7_reg[[1]])), lasso_mut7_reg[[6]], 'recursive'))
# )


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
