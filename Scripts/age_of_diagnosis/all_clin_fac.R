### 
# this script will modles on just clinical predicting binomial and multinomial
# This script will predict a multinomal age model ibrary(dplyr)
library(doParallel)
library(randomForest)
library(caret)
library(nnet)

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

# Create binary variables for age of diagnosis and age of sample collection on both data sets
clin$age_diagnosis_fac <- ifelse(clin$age_diagnosis <= 48, 1, 2)

clin$age_sample_fac<- ifelse(clin$age_sample_collection <= 48, 1, 2)


# Create multinomial variables for age of diagnosis and age of sample collection
clin$age_diagnosis_multi <- ifelse(clin$age_diagnosis <= 30, 1, 
                                        ifelse(clin$age_diagnosis > 30 & clin$age_diagnosis <= 50, 2,
                                               ifelse(clin$age_diagnosis > 50 & clin$age_diagnosis <= 300, 3, 4)))

clin$age_sample_multi <- ifelse(clin$age_sample_collection <= 30, 1, 
                                     ifelse(clin$age_sample_collection > 30 & clin$age_sample_collection <= 50, 2,
                                            ifelse(clin$age_sample_collection > 50 & clin$age_sample_collection <= 300, 3, 4)))



# Create multinomal variables for age of diagnosis and age of sample collection
clin$age_diagnosis_year <- ifelse(clin$age_diagnosis <= 24, 1, 
                                       ifelse(clin$age_diagnosis > 24 & clin$age_diagnosis <= 48, 2,
                                              ifelse(clin$age_diagnosis > 48 & clin$age_diagnosis <= 180, 3,
                                                     ifelse(clin$age_diagnosis > 180 & clin$age_diagnosis <= 360, 4, 5))))

clin$age_sample_year <- ifelse(clin$age_diagnosis <= 24, 1, 
                                    ifelse(clin$age_diagnosis > 24 & clin$age_diagnosis <= 48, 2,
                                           ifelse(clin$age_diagnosis > 48 & clin$age_diagnosis <= 180, 3,
                                                  ifelse(clin$age_diagnosis > 180 & clin$age_diagnosis <= 360, 4, 5))))


#############################################################################################
# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       subset, 
                       selected_features,
                       binary,
                       multi,
                       iterations) {
  
  model <- list()
  predictions <- list()
  test.ground_truth <- list()
  rf.test_acc <- list()
  rf.test_stats <- list()
  
  data <- data[, subset]
  
  # Try the model with all different selection of features based on number of missinginess. 
  data <- data[complete.cases(data),]
  
  obs <- nrow(data)
  
  # convert characters to factors 
  for ( i in 1:ncol(data)){
    
    if(typeof(data[,i]) == 'character') {
      data[,i] <- as.factor(data[,i])
      print(i)
    } 
  }
  
  for (i in 1:iterations){
    
    set.seed(i)
    temp_levels <- 1
    loop_count <- i
    
    variable <- setdiff(subset, selected_features)
    
    if (binary){
      level_count <- 2
    } else if (multi) {
      level_count <- 4
    } else {
      level_count <- 5
    }
    
    while (length(temp_levels) < level_count) {
      set.seed(loop_count)
      train_index <- sample(nrow(data), nrow(data) *.7, replace = F)
      temp_levels <- levels(as.factor(data[, variable][-train_index]))
      # if(length(temp_levels) < level_count) {
      # loop_count <- i + 1
      # }
      loop_count <- loop_count + 1
    }
    
    rf_y = as.factor(make.names(data[, variable])[train_index])
    
    # 4) Random Forest 
    
    if (length(levels(rf_y)) == 2) {
      summaryFunc <- twoClassSummary
    } else {
      summaryFunc <- multiClassSummary
    }
    
    # determines how you train the model.
    fitControl <- trainControl( 
      # method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
      number = 1, 
      classProbs = TRUE,  
      repeats = 1,
      allowParallel = TRUE,
      summaryFunction = summaryFunc)
    
    # mtry <- sqrt(ncol(data))
    # # tunegrid <- expand.grid(.mtry=mtry)
    # model[[i]] <- randomForest(x = data[train_index, selected_features]
    #              , y = rf_y)

    model[[i]] <- train(x = data[train_index, selected_features]
                        , y = rf_y
                        , method = "rf"
                        , trControl = fitControl
                        , verbose = FALSE
                        , metric = "ROC")
    
    predictions[[i]] <- predict(model[[i]] 
                                , newdata = data[-train_index, selected_features]
                                , type = "prob")
    
    test.ground_truth[[i]] <- as.factor(data[, variable][-train_index])
    # test.ground_truth_1inN <- as.factor(class.ind(data[, variable][-train_index]))
    #print(table(levels(test.ground_truth)[apply(rf.predictions, 1, which.is.max)], test.ground_truth))
    # AUC
    # #create ROCR prediction object
    # temp.predict <- prediction(unlist(predictions[[i]]), test.ground_truth_1inN)  
    # 
    # if(auc){
    #   rf.test_auc <- unlist(slot(performance(temp.predict, "auc"), "y.values"))  
    #   print(paste("RF test AUC:", rf.test_auc))  
    # }
    
    # Accuracy
    rf.test_acc[[i]] <- sum(levels(test.ground_truth[[i]])[apply(predictions[[i]], 1, which.is.max)] == test.ground_truth[[i]]) / dim(predictions[[i]])[1]
    # print(paste("RF test acc:", rf.test_acc))  
    # Compute Confusion Matrix and Statistics
    #confusionMatrix(pred, truth)
    rf.test_stats[[i]] <- confusionMatrix(levels(test.ground_truth[[i]])[apply(predictions[[i]], 1, which.is.max)], test.ground_truth[[i]])
    # print(rf.test_stats)
    
    print(i)
    
  }
  
  return(list(predictions, test.ground_truth, rf.test_acc, model, rf.test_stats, obs))
  
}


##################################################################################################################3
# Random forest Predict each variable with gender
# pdf('/home/benbrew/Desktop/rf_clin.pdf')
#binomial
rf_mdm2 <- predictAll(data = clin, 
                      subset <- c("age_diagnosis_fac", "gender", "mdm2.nG"), 
                      selected_features = c("gender", "mdm2.nG"), 
                      binary = T,
                      multi = F,
                      iterations = 10)


mean(unlist(rf_mdm2[[3]]), na.rm = T)


#multiomial
rf_mdm2_multi <- predictAll(data = clin, 
                      subset <- c("age_diagnosis_multi", "gender", "mdm2.nG"), 
                      selected_features = c("gender", "mdm2.nG"),
                      binary = F,
                      multi = T,
                      iterations = 10)

mean(unlist(rf_mdm2_multi[[3]]), na.rm = T)

#multiomial year
rf_mdm2_year <- predictAll(data = clin, 
                      subset <- c("age_diagnosis_year", "gender", "mdm2.nG"), 
                      selected_features = c("gender", "mdm2.nG"), 
                      binary = F,
                      multi = F,
                      iterations = 10)

mean(unlist(rf_mdm2_year[[3]]), na.rm = T)

#######################
#binomial
rf_protein_codon_num <- predictAll(data = clin, 
                                   subset <- c("age_diagnosis_fac","gender", "protein.codon.num"),
                                   binary = TRUE,
                                   selected_features = c("gender", "protein.codon.num"), 
                                   iterations = 10)

mean(unlist(rf_protein_codon_num[[3]]))
#multinomial
rf_protein_codon_num_multi <- predictAll(data = clin, 
                                   subset <- c("age_diagnosis_multi","gender", "protein.codon.num"),
                                   binary = TRUE,
                                   selected_features = c("gender", "protein.codon.num"), 
                                   iterations = 10)

mean(unlist(rf_protein_codon_num_multi[[3]]))
#multinomial year
rf_protein_codon_num_year <- predictAll(data = clin, 
                                   subset <- c("age_diagnosis_year","gender", "protein.codon.num"),
                                   binary = TRUE,
                                   selected_features = c("gender", "protein.codon.num"), 
                                   iterations = 10)

mean(unlist(rf_protein_codon_num_year[[3]]))

#########################
#binomial
rf_splice_delins_snv <- predictAll(data = clin, 
                                   subset <- c("age_diagnosis_fac","gender", "splice.delins.snv"), 
                                   selected_features = c("gender", "splice.delins.snv"), 
                                   binary = T,
                                   iterations = 10)


mean(unlist(rf_splice_delins_snv[[3]]))
#multinomial
rf_splice_delins_snv_multi <- predictAll(data = clin, 
                                   subset <- c("age_diagnosis_multi","gender", "splice.delins.snv"), 
                                   selected_features = c("gender", "splice.delins.snv"),
                                   binary = F,
                                   multi = T,
                                   iterations = 10)

mean(unlist(rf_splice_delins_snv_multi[[3]]))

#multinomial year
rf_splice_delins_snv_year <- predictAll(data = clin, 
                                   subset <- c("age_diagnosis_year","gender", "splice.delins.snv"), 
                                   selected_features = c("gender", "splice.delins.snv"), 
                                   binary = F,
                                   multi = F,
                                   iterations = 10)
mean(unlist(rf_splice_delins_snv_year[[3]]))

##########################
#binomial
rf_codon72.npro <- predictAll(data = clin, 
                              subset <- c("age_diagnosis_fac", "gender", "codon72.npro"), 
                              selected_features = c("gender", "codon72.npro"), 
                              binary = T,
                              iterations = 10)

mean(unlist(rf_codon72.npro[[3]]))
#multinomial
rf_codon72.npro_multi <- predictAll(data = clin, 
                              subset <- c("age_diagnosis_multi", "gender", "codon72.npro"), 
                              selected_features = c("gender", "codon72.npro"), 
                              binary = F,
                              multi = T,
                              iterations = 10)
mean(unlist(rf_codon72.npro_multi[[3]]))

#multinomial year
rf_codon72.npro_year <- predictAll(data = clin, 
                              subset <- c("age_diagnosis_year", "gender", "codon72.npro"), 
                              selected_features = c("gender", "codon72.npro"), 
                              binary = F,
                              multi = F,
                              iterations = 10)
mean(unlist(rf_codon72.npro_year[[3]]))

######################################
#binomial
rf_gdna_exon_intron <- predictAll(data = clin, 
                                  subset <- c("age_diagnosis_fac", "gender", "gdna.exon.intron"), 
                                  selected_features = c("gender", "gdna.exon.intron"), 
                                  binary = T,
                                  iterations = 10)

mean(unlist(rf_gdna_exon_intron[[3]]))
#best at 62
# look at confustion matrices
temp <- list()
for (i in 1:10){
  temp[[i]] <- rf_gdna_exon_intron[[5]][[i]]$table
  
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)
spot_1<- (mat[1] + mat[5] + mat[9] + mat[13] + mat[17] + mat[21] +mat[25] + mat[29] + mat[33] + mat[37])/10
spot_2 <- (mat[2] + mat[6] + mat[10] + mat[14] + mat[18] + mat[22] +mat[26] + mat[30] + mat[34] + mat[38])/10
spot_3<- (mat[3] + mat[7] + mat[11] + mat[15] + mat[19] + mat[23] +mat[27] + mat[31] + mat[35] + mat[39])/10
spot_4 <- (mat[4] + mat[8] + mat[12] + mat[16] + mat[20] + mat[24] +mat[28] + mat[32] + mat[36] + mat[40])/10

new_mat[1,1] <- spot_1
new_mat[2,1] <- spot_2
new_mat[1,2] <- spot_3
new_mat[2,2] <- spot_4


#multinomial
rf_gdna_exon_intron_multi <- predictAll(data = clin, 
                                  subset <- c("age_diagnosis_multi", "gender", "gdna.exon.intron"), 
                                  selected_features = c("gender", "gdna.exon.intron"), 
                                  binary = F,
                                  multi = T,
                                  iterations = 10)

mean(unlist(rf_gdna_exon_intron_multi[[3]]))

#multinomial year
rf_gdna_exon_intron_year <- predictAll(data = clin, 
                                  subset <- c("age_diagnosis_year", "gender", "gdna.exon.intron"), 
                                  selected_features = c("gender", "gdna.exon.intron"), 
                                  binary = F,
                                  multi = F,
                                  iterations = 10)
mean(unlist(rf_gdna_exon_intron_year[[3]]))



######################################
#binomial
rf_protein_codon_change <- predictAll(data = clin, 
                                      subset <- c("age_diagnosis_fac", "gender", "protein.codon.change"), 
                                      selected_features = c("gender", "protein.codon.change"), 
                                      binary = T,
                                      iterations = 10)


mean(unlist(rf_protein_codon_change[[3]]))
#multi
rf_protein_codon_change_multi <- predictAll(data = clin, 
                                      subset <- c("age_diagnosis_multi", "gender", "protein.codon.change"), 
                                      selected_features = c("gender", "protein.codon.change"), 
                                      binary = F,
                                      multi = T,
                                      iterations = 10)

mean(unlist(rf_protein_codon_change_multi[[3]]))

#multi year 
rf_protein_codon_change_year <- predictAll(data = clin, 
                                      subset <- c("age_diagnosis_year", "gender", "protein.codon.change"), 
                                      selected_features = c("gender", "protein.codon.change"), 
                                      binary = F,
                                      multi = F,
                                      iterations = 10)

mean(unlist(rf_protein_codon_change_year[[3]]))
#here
######################################
#binomial
rf_gdna_codon <- predictAll(data = clin, 
                            subset <- c("age_diagnosis_fac","gender", "gdna.codon"), 
                            selected_features = c("gender", "gdna.codon"), 
                            binary = T,
                            iterations = 10)

mean(unlist(rf_gdna_codon[[3]]))
#multi
rf_gdna_codon_multi <- predictAll(data = clin, 
                            subset <- c("age_diagnosis_multi","gender", "gdna.codon"), 
                            selected_features = c("gender", "gdna.codon"), 
                            binary = F,
                            multi = T,
                            iterations = 10)

mean(unlist(rf_gdna_codon_multi[[3]]))

#multi year
rf_gdna_codon_year <- predictAll(data = clin, 
                            subset <- c("age_diagnosis_year","gender", "gdna.codon"), 
                            selected_features = c("gender", "gdna.codon"), 
                            binary = F,
                            multi = F,
                            iterations = 10)

mean(unlist(rf_gdna_codon_year[[3]]))

################################
#binomial
rf_gdna_base_change <- predictAll(data = clin, 
                                  subset <- c("age_diagnosis_fac", "gender", "gdna.base.change"), 
                                  selected_features = c("gender", "gdna.base.change"), 
                                  binary = T,
                                  iterations = 10)

mean(unlist(rf_gdna_base_change[[3]]))
# maybe 2nd best at 56

#multi
rf_gdna_base_change_multi <- predictAll(data = clin, 
                                  subset <- c("age_diagnosis_multi", "gender", "gdna.base.change"), 
                                  selected_features = c("gender", "gdna.base.change"), 
                                  binary = F,
                                  multi = T,
                                  iterations = 10)


mean(unlist(rf_gdna_base_change_multi[[3]]))

#multi year
rf_gdna_base_change_year <- predictAll(data = clin, 
                                  subset <- c("age_diagnosis_year", "gender", "gdna.base.change"), 
                                  selected_features = c("gender", "gdna.base.change"), 
                                  binary = F,
                                  multi = F,
                                  iterations = 10)

mean(unlist(rf_gdna_base_change_year[[3]]))

##############################
# Combinations 
rf_gdna_exon_intron_gdna_base_change <- predictAll(data = clin, 
                                  subset <- c("age_diagnosis_fac", "gender", "gdna.exon.intron", "gdna.base.change"), 
                                  selected_features = c("gender", "gdna.exon.intron", "gdna.base.change"), 
                                  binary = T,
                                  iterations = 10)


mean(unlist(rf_gdna_exon_intron_gdna_base_change[[3]]))

#########################33
rf_gdna_exon_intron_gdna_base_change_multi <- predictAll(data = clin, 
                                                   subset <- c("age_diagnosis_fac", "gender", "gdna.exon.intron", "gdna.base.change"), 
                                                   selected_features = c("gender", "gdna.exon.intron", "gdna.base.change"), 
                                                   binary = F,
                                                   multi = T,
                                                   iterations = 10)


mean(unlist(rf_gdna_exon_intron_gdna_base_change_multi[[3]]))

#####################################3
rf_gdna_exon_intron_gdna_base_change_year <- predictAll(data = clin, 
                                                   subset <- c("age_diagnosis_fac", "gender", "gdna.exon.intron", "gdna.base.change"), 
                                                   selected_features = c("gender", "gdna.exon.intron", "gdna.base.change"), 
                                                   binary = F,
                                                   multi = F,
                                                   iterations = 10)


mean(unlist(rf_gdna_exon_intron_gdna_base_change_year[[3]]))


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