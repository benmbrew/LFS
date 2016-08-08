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
  
  variables <- nrow(clin)
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
      
      temp <- varImp(model[[i]])[[1]]
      importance[[i]] <- cbind(rownames(temp), temp$Overall)

      
      
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
  
  return(list(mse, predictions, model, importance, test.ground_truth, variables))
  
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
# gender and gdna.base.change
rf_mut <- predictAll(model_name = 'rf', 
                     subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                     selected_features = c("gender", "gdna.base.change"), iterations = 50)

# add gdna.codon
rf_mut1 <- predictAll(model_name = 'rf', 
                     subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                     selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 50)

# add protein.codon.change
rf_mut2 <- predictAll(model_name = 'rf', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), iterations = 50)

# add gdna.exon.intron
rf_mut3 <- predictAll(model_name = 'rf', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron"), iterations = 50)

# add codon72.npro
rf_mut4 <- predictAll(model_name = 'rf', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro"), iterations = 50)

# add splice.delins.snv
rf_mut5 <- predictAll(model_name = 'rf', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), iterations = 50)

# add protein.codon.num
rf_mut6 <- predictAll(model_name = 'rf', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), iterations = 50)

# add mdm2.nG
rf_mut7 <- predictAll(model_name = 'rf', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), iterations = 50)

# combine all tests 

rf <- rbind (
  append('gender_and_gdna.base.change', c(mean(unlist(rf_mut[[1]])), rf_mut[[6]])),
  append('gdna.codon', c(mean(unlist(rf_mut1[[1]])), rf_mut1[[6]])),
  append('protein.codon.change', c(mean(unlist(rf_mut2[[1]])), rf_mut2[[6]])),
  append('gdna.exon.intron', c(mean(unlist(rf_mut3[[1]])), rf_mut3[[6]])),
  append('codon72.npro', c(mean(unlist(rf_mut4[[1]])), rf_mut4[[6]])),
  append('splice.delins.snv', c(mean(unlist(rf_mut5[[1]])), rf_mut5[[6]])),
  append('protein.codon.num', c(mean(unlist(rf_mut6[[1]])), rf_mut6[[6]])),
  append('mdm2.nG', c(mean(unlist(rf_mut7[[1]])), rf_mut7[[6]]))
)

####
# check variable importance 
rf_mut7[[4]]
temp <- as.data.frame(do.call('rbind', rf_mut7[[4]]))
temp$V2 <- as.numeric(as.character(temp$V2))
import <- temp %>%
  group_by(V1) %>%
  summarise(mean_import = mean(V2))
#gdna.codon, protein.codon.num, mdm2.nG,protein.codon.change
  
########################################################################################
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
# Random forest Predict backwards
# gender and mdm2.nG

rf_mut_back <- predictAll(model_name = 'rf', 
                     subset <- c("age_diagnosis", "gender", "mdm2.nG"), 
                     selected_features = c("gender", "mdm2.nG"), iterations = 50)


# add protein.codon.num
rf_mut1_back <- predictAll(model_name = 'rf', 
                      subset <- c("age_diagnosis", "gender", "mdm2.nG", "protein.codon.num"), 
                      selected_features = c("gender", "mdm2.nG", "protein.codon.num"), iterations = 50)


# add splice.delins.snv
rf_mut2_back <- predictAll(model_name = 'rf', 
                      subset <- c("age_diagnosis", "gender", "mdm2.nG", "protein.codon.num", "splice.delins.snv"), 
                      selected_features = c("gender", "mdm2.nG", "protein.codon.num", "splice.delins.snv"), iterations = 50)


# add codon72.npro
rf_mut3_back <- predictAll(model_name = 'rf', 
                      subset <- c("age_diagnosis", "gender", "mdm2.nG", "protein.codon.num", "splice.delins.snv",
                                  "codon72.npro"), 
                      selected_features = c("gender", "mdm2.nG", "protein.codon.num", "splice.delins.snv",
                                            "codon72.npro"), iterations = 50)

# add gdna.exon.intron
rf_mut4_back <- predictAll(model_name = 'rf', 
                      subset <- c("age_diagnosis", "gender", "mdm2.nG", "protein.codon.num", "splice.delins.snv",
                                  "codon72.npro", "gdna.exon.intron"), 
                      selected_features = c("gender", "mdm2.nG", "protein.codon.num", "splice.delins.snv",
                                            "codon72.npro", "gdna.exon.intron"), iterations = 50)

# add protein.codon.change
rf_mut5_back <- predictAll(model_name = 'rf', 
                      subset <- c("age_diagnosis", "gender", "mdm2.nG", "protein.codon.num", "splice.delins.snv",
                                  "codon72.npro", "gdna.exon.intron", "protein.codon.change"), 
                      selected_features = c("gender", "mdm2.nG", "protein.codon.num", "splice.delins.snv",
                                            "codon72.npro", "gdna.exon.intron", "protein.codon.change"), iterations = 50)


# add gdna.codon
rf_mut6_back <- predictAll(model_name = 'rf', 
                     subset <- c("age_diagnosis", "gender", "mdm2.nG", "protein.codon.num", "splice.delins.snv",
                                 "codon72.npro", "gdna.exon.intron", "protein.codon.change", "gdna.codon"), 
                     selected_features = c("gender", "mdm2.nG", "protein.codon.num", "splice.delins.snv",
                                           "codon72.npro", "gdna.exon.intron", "protein.codon.change", "gdna.codon"), iterations = 50)



# add gdna.base.change
rf_mut7_back <- predictAll(model_name = 'rf', 
                      subset <- c("age_diagnosis", "gender", "mdm2.nG", "protein.codon.num", "splice.delins.snv",
                                  "codon72.npro", "gdna.exon.intron", "protein.codon.change", "gdna.codon", "gdna.base.change"), 
                      selected_features = c("gender", "mdm2.nG", "protein.codon.num", "splice.delins.snv",
                                            "codon72.npro", "gdna.exon.intron", "protein.codon.change", "gdna.codon", "gdna.base.change"), iterations = 50)

# combine all tests 

rf_back <- rbind (
  append('gender_and_mdm2.nG', c(mean(unlist(rf_mut_back[[1]])), rf_mut_back[[6]])),
  append('protein.codon.num', c(mean(unlist(rf_mut1_back[[1]])), rf_mut1_back[[6]])),
  append('splice.delins.snv', c(mean(unlist(rf_mut2_back[[1]])), rf_mut2_back[[6]])),
  append('codon72.npro', c(mean(unlist(rf_mut3_back[[1]])), rf_mut3_back[[6]])),
  append('gdna.exon.intron', c(mean(unlist(rf_mut4_back[[1]])), rf_mut4_back[[6]])),
  append('protein.codon.change', c(mean(unlist(rf_mut5_back[[1]])), rf_mut5_back[[6]])),
  append('gdna.codon', c(mean(unlist(rf_mut6_back[[1]])), rf_mut6_back[[6]])),
  append('gdna.base.change', c(mean(unlist(rf_mut7_back[[1]])), rf_mut7_back[[6]]))
)


##################################################################################################################3
# Random forest Predict each variable with gender

rf_mut_mdm2 <- predictAll(model_name = 'rf', 
                          subset <- c("age_diagnosis", "gender", "mdm2.nG"), 
                          selected_features = c("gender", "mdm2.nG"), iterations = 50)

rf_mut_protein_codon_num <- predictAll(model_name = 'rf', 
                                       subset <- c("age_diagnosis", "gender", "protein.codon.num"), 
                                       selected_features = c("gender", "protein.codon.num"), iterations = 50)

rf_mut_splice_delins_snv <- predictAll(model_name = 'rf', 
                                       subset <- c("age_diagnosis", "gender", "splice.delins.snv"), 
                                       selected_features = c("gender", "splice.delins.snv"), iterations = 50)

rf_mut_codon72.npro <- predictAll(model_name = 'rf', 
                                  subset <- c("age_diagnosis", "gender", "codon72.npro"), 
                                  selected_features = c("gender", "codon72.npro"), iterations = 50)

rf_mut_gdna_exon_intron <- predictAll(model_name = 'rf', 
                                      subset <- c("age_diagnosis", "gender", "gdna.exon.intron"), 
                                      selected_features = c("gender", "gdna.exon.intron"), iterations = 50)

rf_mut_protein_codon_change <- predictAll(model_name = 'rf', 
                                          subset <- c("age_diagnosis", "gender", "protein.codon.change"), 
                                          selected_features = c("gender", "protein.codon.change"), iterations = 50)

rf_mut_gdna_codon <- predictAll(model_name = 'rf', 
                                subset <- c("age_diagnosis", "gender", "gdna.codon"), 
                                selected_features = c("gender", "gdna.codon"), iterations = 50)

rf_mut_gdna_base_change <- predictAll(model_name = 'rf', 
                                      subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                                      selected_features = c("gender", "gdna.base.change"), iterations = 50)
# combine all tests 

rf_single <- rbind (
  append('gender_and_mdm2.nG', c(mean(unlist(rf_mut_mdm2[[1]])), rf_mut_mdm2[[6]])),
  append('protein.codon.num', c(mean(unlist(rf_mut_protein_codon_num[[1]])), rf_mut_protein_codon_num[[6]])),
  append('splice.delins.snv', c(mean(unlist(rf_mut_splice_delins_snv[[1]])), rf_mut_splice_delins_snv[[6]])),
  append('codon72.npro', c(mean(unlist(rf_mut_codon72.npro[[1]])), rf_mut_codon72.npro[[6]])),
  append('gdna.exon.intron', c(mean(unlist(rf_mut_gdna_exon_intron[[1]])), rf_mut_gdna_exon_intron[[6]])),
  append('protein.codon.change', c(mean(unlist(rf_mut_protein_codon_change[[1]])), rf_mut_protein_codon_change[[6]])),
  append('gdna.codon', c(mean(unlist(rf_mut_gdna_codon[[1]])), rf_mut_gdna_codon[[6]])),
  append('gdna.base.change', c(mean(unlist(rf_mut_gdna_base_change[[1]])), rf_mut_gdna_base_change[[6]]))
)


########################################################################3
# Random forest Predict backwards
# gender and mdm2.nG - combine with every variable once 

rf_mut_mdm2 <- predictAll(model_name = 'rf', 
                          subset <- c("age_diagnosis", "gender", "mdm2.nG"), 
                          selected_features = c("gender", "mdm2.nG"), iterations = 50)

rf_mut_protein_codon_num <- predictAll(model_name = 'rf', 
                          subset <- c("age_diagnosis", "gender", "mdm2.nG", "protein.codon.num"), 
                          selected_features = c("gender", "mdm2.nG", "protein.codon.num"), iterations = 50)

rf_mut_splice_delins_snv <- predictAll(model_name = 'rf', 
                          subset <- c("age_diagnosis", "gender", "mdm2.nG", "splice.delins.snv"), 
                          selected_features = c("gender", "mdm2.nG", "splice.delins.snv"), iterations = 50)

rf_mut_codon72.npro <- predictAll(model_name = 'rf', 
                          subset <- c("age_diagnosis", "gender", "mdm2.nG", "codon72.npro"), 
                          selected_features = c("gender", "mdm2.nG", "codon72.npro"), iterations = 50)

rf_mut_gdna_exon_intron <- predictAll(model_name = 'rf', 
                          subset <- c("age_diagnosis", "gender", "mdm2.nG", "gdna.exon.intron"), 
                          selected_features = c("gender", "mdm2.nG", "gdna.exon.intron"), iterations = 50)

rf_mut_protein_codon_change <- predictAll(model_name = 'rf', 
                          subset <- c("age_diagnosis", "gender", "mdm2.nG", "protein.codon.change"), 
                          selected_features = c("gender", "mdm2.nG", "protein.codon.change"), iterations = 50)

rf_mut_gdna_codon <- predictAll(model_name = 'rf', 
                          subset <- c("age_diagnosis", "gender", "mdm2.nG", "gdna.codon"), 
                          selected_features = c("gender", "mdm2.nG", "gdna.codon"), iterations = 50)

rf_mut_gdna_base_change <- predictAll(model_name = 'rf', 
                          subset <- c("age_diagnosis", "gender", "mdm2.nG", "gdna.base.change"), 
                          selected_features = c("gender", "mdm2.nG", "gdna.base.change"), iterations = 50)
# combine all tests 

rf_mdm2 <- rbind (
  append('gender_and_mdm2.nG', c(mean(unlist(rf_mut_mdm2[[1]])), rf_mut_mdm2[[6]])),
  append('protein.codon.num', c(mean(unlist(rf_mut_protein_codon_num[[1]])), rf_mut_protein_codon_num[[6]])),
  append('splice.delins.snv', c(mean(unlist(rf_mut_splice_delins_snv[[1]])), rf_mut_splice_delins_snv[[6]])),
  append('codon72.npro', c(mean(unlist(rf_mut_codon72.npro[[1]])), rf_mut_codon72.npro[[6]])),
  append('gdna.exon.intron', c(mean(unlist(rf_mut_gdna_exon_intron[[1]])), rf_mut_gdna_exon_intron[[6]])),
  append('protein.codon.change', c(mean(unlist(rf_mut_protein_codon_change[[1]])), rf_mut_protein_codon_change[[6]])),
  append('gdna.codon', c(mean(unlist(rf_mut_gdna_codon[[1]])), rf_mut_gdna_codon[[6]])),
  append('gdna.base.change', c(mean(unlist(rf_mut_gdna_base_change[[1]])), rf_mut_gdna_base_change[[6]]))
)



##################################################################################################################3
# Elastic Net
# gender and gdna.base.change
enet_mut <- predictAll(model_name = 'enet', 
                     subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                     selected_features = c("gender", "gdna.base.change"), iterations = 50)

plot(unlist(enet_mut[[2]]), unlist(enet_mut[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# add gdna.codon
enet_mut1 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 50)

plot(unlist(enet_mut1[[2]]), unlist(enet_mut1[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.change
enet_mut2 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), iterations = 50)

plot(unlist(enet_mut2[[2]]), unlist(enet_mut2[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add gdna.exon.intron
enet_mut3 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron"), iterations = 50)

plot(unlist(enet_mut3[[2]]), unlist(enet_mut3[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
enet_mut4 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro"), iterations = 50)

plot(unlist(enet_mut4[[2]]), unlist(enet_mut4[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add splice.delins.snv
enet_mut5 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), iterations = 50)

plot(unlist(enet_mut5[[2]]), unlist(enet_mut5[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add protein.codon.num
enet_mut6 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), iterations = 50)

plot(unlist(enet_mut6[[2]]), unlist(enet_mut6[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add mdm2.nG
enet_mut7 <- predictAll(model_name = 'enet', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), iterations = 50)

plot(unlist(enet_mut7[[2]]), unlist(enet_mut7[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# combine all tests 

enet <- rbind (
  append('gender_and_gdna.base.change', c(mean(unlist(enet_mut[[1]])), enet_mut[[6]])),
  append('gdna.codon', c(mean(unlist(enet_mut1[[1]])), enet_mut1[[6]])),
  append('protein.codon.change', c(mean(unlist(enet_mut2[[1]])), enet_mut2[[6]])),
  append('gdna.exon.intron', c(mean(unlist(enet_mut3[[1]])), enet_mut3[[6]])),
  append('codon72.npro', c(mean(unlist(enet_mut4[[1]])), enet_mut4[[6]])),
  append('splice.delins.snv', c(mean(unlist(enet_mut5[[1]])), enet_mut5[[6]])),
  append('protein.codon.num', c(mean(unlist(enet_mut6[[1]])), enet_mut6[[6]])),
  append('mdm2.nG', c(mean(unlist(enet_mut7[[1]])), enet_mut7[[6]]))
)
         

##################################################################################################################3
# Lasso
# gender and gdna.base.change
lasso_mut <- predictAll(model_name = 'lasso', 
                     subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                     selected_features = c("gender", "gdna.base.change"), iterations = 50)

plot(unlist(lasso_mut[[2]]), unlist(lasso_mut[[5]]), xlab = 'Predictions', ylab = 'Actual Age')


# add gdna.codon
lasso_mut1 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon"), iterations = 50)

plot(unlist(lasso_mut1[[2]]), unlist(lasso_mut1[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.change
lasso_mut2 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), iterations = 50)

plot(unlist(lasso_mut2[[2]]), unlist(lasso_mut2[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add gdna.exon.intron
lasso_mut3 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron"), iterations = 50)

plot(unlist(lasso_mut3[[2]]), unlist(lasso_mut3[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
lasso_mut4 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro"), iterations = 50)

plot(unlist(lasso_mut4[[2]]), unlist(lasso_mut4[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add splice.delins.snv
lasso_mut5 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), iterations = 50)

plot(unlist(lasso_mut5[[2]]), unlist(lasso_mut5[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add protein.codon.num
lasso_mut6 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), iterations = 50)

plot(unlist(lasso_mut6[[2]]), unlist(lasso_mut6[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add mdm2.nG
lasso_mut7 <- predictAll(model_name = 'lasso', 
                      subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                      selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), iterations = 50)

plot(unlist(lasso_mut7[[2]]), unlist(lasso_mut7[[5]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# combine all tests 

lasso <- rbind (
  append('gender_and_gdna.base.change', c(mean(unlist(lasso_mut[[1]])), lasso_mut[[6]])),
  append('gdna.codon', c(mean(unlist(lasso_mut1[[1]])), lasso_mut1[[6]])),
  append('protein.codon.change', c(mean(unlist(lasso_mut2[[1]])), lasso_mut2[[6]])),
  append('gdna.exon.intron', c(mean(unlist(lasso_mut3[[1]])), lasso_mut3[[6]])),
  append('codon72.npro', c(mean(unlist(lasso_mut4[[1]])), lasso_mut4[[6]])),
  append('splice.delins.snv', c(mean(unlist(lasso_mut5[[1]])), lasso_mut5[[6]])),
  append('protein.codon.num', c(mean(unlist(lasso_mut6[[1]])), lasso_mut6[[6]])),
  append('mdm2.nG', c(mean(unlist(lasso_mut7[[1]])), lasso_mut7[[6]]))
)


##########################################################################################
# combine ran_forest, enet, lasso
ran_forest <- as.data.frame(rf)
ran_forest$model <- 'rand_forest'
names(ran_forest) <- c('variables', 'rmse', 'observations', 'model')

enet <- as.data.frame(enet)
enet$model <- 'enet'
names(enet) <- c('variables', 'rmse', 'observations', 'model')

lasso <- as.data.frame(lasso)
lasso$model <- 'lasso'
names(lasso) <- c('variables', 'rmse', 'observations', 'model')

models <- rbind(ran_forest, enet, lasso)
models$rmse <- as.numeric(as.character(models$rmse))

ggplot(models, aes(model, mse, group = variables, fill = variables)) + 
  geom_bar(stat = 'identity', 
           position = 'dodge')


# write.csv(models, paste0(data_folder, '/clin_model_results.csv'))





