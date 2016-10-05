#############################################################################################################
# This script will be a clean version of running all the dmr data and its residuals on both probe and gene level. 
# required data: methylation (imputed with knn, clean clinical data, different dmr variations, residuals)
# this will also check random features 

##################################################################################################
# this script will read in different versions of subsetted full data and run random forest. 
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(doParallel)
library(nnet)
library(dplyr)

registerDoParallel(1)

###### Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
clin_data <- paste0(data_folder, '/clin_data')

###### load methylation knn imputation
load(paste0(data_folder, '/methyl_knn.RData'))
methyl_knn_probe <- methylation
rm(tab, methylation)

###### load methylation lsa imutation (genes)
methyl <- read.csv(paste0(data_folder, '/methyl_impute_raw.csv'))
methyl$X <- NULL
methyl_lsa_gene <- methyl
rm(methyl)

###### load methylation knn imutation (genes)
methyl <- read.csv(paste0(data_folder, '/methyl_impute_knn.csv'))
methyl$X <- NULL
methyl_knn_gene <- methyl
rm(methyl)

###### Load in clinical data
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)
# make clinical id column and take away special characters 
clin$id <- clin$blood_dna_malkin_lab_ <- gsub('_|A|B', '', clin$blood_dna_malkin_lab_)

###### read in dmr variations

# cancer and global
bh_cancer <- read.csv(paste0(data_folder, '/bh_cancer_full.csv'))
bh_global <- read.csv(paste0(data_folder, '/bh_global_full.csv'))

# subselection
bh_cancer_sub <- read.csv(paste0(data_folder, '/bh_cancer_sub_full.csv'))
bh_global_sub<- read.csv(paste0(data_folder, '/bh_global_sub_full.csv'))

# balanced
bh_cancer_bal <- read.csv(paste0(data_folder, '/bh_cancer_bal_full.csv'))
bh_global_bal <- read.csv(paste0(data_folder, '/bh_global_bal_full.csv'))

# read in union and intersection
bh_union <- read.csv(paste0(data_folder, '/bh_union.csv'))
bh_intersection <- read.csv(paste0(data_folder, '/bh_intersection.csv'))

# read in bumphunter results from cancer indicator
bh_cancer_ind <-  read.csv(paste0(data_folder, '/bh_cancer_ind_full.csv'))

# remove X column from data frames 
bh_cancer$X <- bh_global$X <- bh_cancer_sub$X <- bh_global_sub$X <- bh_cancer_bal$X <- bh_global_bal$X <- bh_union$X <-
  bh_intersection$X <- bh_cancer_ind$X <- NULL


######################################
# function will have to take full_data

# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       threshold,
                       lsa_gene,
                       knn_gene,
                       knn_probe,
                       fac,
                       log,
                       cutoff,
                       residual,
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
  
  # if (random) {
  #   
  #   # remove extra characters 
  #   methyl_knn_probe$id <- gsub('_|A|B', '', methyl_knn_probe$id)
  #   methyl_knn_probe$id <- substr(methyl_knn_probe$id, 1, 4)
  #   
  #   # select probes 
  #   genes <- sample(colnames(methyl_knn_probe), threshold, replace = F)
  #   methyl_knn_probe <- methyl_knn_probe[, c(2:genes)]
  #   
  #   #join methyl_gene and clin
  #   model_data <- inner_join(methyl_knn_probe, clin, by = 'id')
  #   model_data <- model_data[!is.na(model_data$age_diagnosis),]
  #   
  #   # remove duplicates
  #   model_data <- model_data[!duplicated(model_data$id),]
  #   # keep only biological data and age data
  #   model_data <- model_data[c('age_diagnosis', 'age_sample_collection', probe_intersection)]
  #   
  #   
  # }
  
  
  if (lsa_gene) {
    
    # get vector of gene names from bumphunter with cutoff 
    genes <- data$nearestGeneSymbol[1:threshold]
    # get intersection of genes and colnames of methly_gene
    gene_intersection <- intersect(genes, colnames(methyl_lsa_gene))
    # subset methyl gene by bumphunter data 
    methyl_lsa_gene <- cbind(id = methyl_lsa_gene$id, methyl_lsa_gene[, gene_intersection])
    
    # remove extra characters 
    methyl_lsa_gene$id <- gsub('_|A|B', '', methyl_lsa_gene$id)
    methyl_lsa_gene$id <- substr(methyl_lsa_gene$id, 1, 4)
    
    #join methyl_gene and clin
    model_data <- inner_join(methyl_lsa_gene, clin, by = 'id')
    model_data <- model_data[!is.na(model_data$age_diagnosis),]
    
    # remove duplicates
    model_data <- model_data[!duplicated(model_data$id),]
    
    # keep only biological data and age data
    model_data <- model_data[c('age_diagnosis', 'age_sample_collection', gene_intersection)]
    
    if (fac) {
      
      model_data$age_diagnosis_fac <- as.integer(ifelse(model_data$age_diagnosis <= 48, 1, 2))
      model_data$age_sample_fac <- as.integer(ifelse(model_data$age_sample_collection <= 48, 1, 2))
      
      # keep only biological data and age data
      model_data <- model_data[c('age_diagnosis_fac', 'age_sample_fac', gene_intersection)]
      
    }
    
  }
  
  if (knn_gene) {
    
    # get vector of gene names from bumphunter with cutoff 
    genes <- data$nearestGeneSymbol[1:threshold]
    # get intersection of genes and colnames of methly_gene
    gene_intersection <- intersect(genes, colnames(methyl_knn_gene))
    # subset methyl gene by bumphunter data 
    methyl_knn_gene <- cbind(id = methyl_knn_gene$id, methyl_knn_gene[, gene_intersection])
    
    # remove extra characters 
    methyl_knn_gene$id <- gsub('_|A|B', '', methyl_knn_gene$id)
    methyl_knn_gene$id <- substr(methyl_knn_gene$id, 1, 4)
    
    #join methyl_gene and clin
    model_data <- inner_join(methyl_knn_gene, clin, by = 'id')
    model_data <- model_data[!is.na(model_data$age_diagnosis),]
    
    # remove duplicates
    model_data <- model_data[!duplicated(model_data$id),]
    
    # keep only biological data and age data
    model_data <- model_data[c('age_diagnosis', 'age_sample_collection', gene_intersection)]
    
    if (fac) {
      
      model_data$age_diagnosis_fac <- as.integer(ifelse(model_data$age_diagnosis <= 48, 1, 2))
      model_data$age_sample_fac <- as.integer(ifelse(model_data$age_sample_collection <= 48, 1, 2))
      
      # keep only biological data and age data
      model_data <- model_data[c('age_diagnosis_fac', 'age_sample_fac', gene_intersection)]
      
    }
    
  }
  
  if (knn_probe) {
    
    # get vector of gene names from bumphunter with cutoff 
    probes <- data$probe[1:threshold]
    # get intersection of genes and colnames of methly_gene
    probe_intersection <- intersect(probes, colnames(methyl_knn_probe))
    # subset methyl gene by bumphunter data 
    methyl_knn_probe <- cbind(id = methyl_knn_probe$id, methyl_knn_probe[, probe_intersection])
    # remove extra characters 
    methyl_knn_probe$id <- gsub('_|A|B', '', methyl_knn_probe$id)
    methyl_knn_probe$id <- substr(methyl_knn_probe$id, 1, 4)
    #join methyl_gene and clin
    model_data <- inner_join(methyl_knn_probe, clin, by = 'id')
    model_data <- model_data[!is.na(model_data$age_diagnosis),]
    
    # remove duplicates
    model_data <- model_data[!duplicated(model_data$id),]
    # keep only biological data and age data
    model_data <- model_data[c('age_diagnosis', 'age_sample_collection', probe_intersection)]
    
    if (fac) {
      
      model_data$age_diagnosis_fac <- as.integer(ifelse(model_data$age_diagnosis <= 48, 1, 2))
      model_data$age_sample_fac <- as.integer(ifelse(model_data$age_sample_collection <= 48, 1, 2))
      
      # keep only biological data and age data
      model_data <- model_data[c('age_diagnosis_fac', 'age_sample_fac', probe_intersection)]
      
    }
    
  }
  
  
  if (residual & !fac) {
    
    feature_names <- colnames(model_data)[3:ncol(model_data)]
    
    # first subset data to complete cases- that is keep only samples where age of diagnosis and age of sample collection 
    # are not missing 
    
    model_data <- model_data[complete.cases(model_data),]
    
    resid <- list()
    
    for (i in 3:ncol(model_data)){
      
      resid[[i]] <- lm(model_data[, i] ~ model_data$age_sample_collection, data = model_data)$residuals
      
      print(i)
      
    }
    
    resid_data <- as.data.frame(do.call('cbind', resid))
    model_data <- cbind(model_data$age_diagnosis, model_data$age_sample_collection, resid_data)
    colnames(model_data) <- c('age_diagnosis', 'age_sample_collection', feature_names)
    
    
    if (log) {
      
      # this condition is only observed if residual is also true -> this squares the negatives so log can be taken.
      
      # get min of model data and take absolute value then add that + 1 to model data
      model_data[, 3:ncol(model_data)] <-  model_data[, 3:ncol(model_data)] + abs(min(model_data[, 3:ncol(model_data)])) + .1
      model_data <- log(model_data)
      
    }
    
  }
  
  if (residual & fac) {
    
    feature_names <- colnames(model_data)[3:ncol(model_data)]
    
    # first subset data to complete cases- that is keep only samples where age of diagnosis and age of sample collection 
    # are not missing 
    
    model_data <- model_data[complete.cases(model_data),]
    
    resid <- list()
    
    model_data$age_sample_fac <- as.factor(model_data$age_sample_fac)
    
    for (i in 3:ncol(model_data)){
      
      resid[[i]] <- lm(model_data[, i] ~ model_data$age_sample_fac, data = model_data)$residuals
      
      print(i)
      
    }
    
    resid_data <- as.data.frame(do.call('cbind', resid))
    model_data <- cbind(model_data$age_diagnosis_fac, model_data$age_sample_fac, resid_data)
    colnames(model_data) <- c('age_diagnosis_fac', 'age_sample_fac', feature_names)
    
  }
  
  
  if (log & !residual) {
    
    model_data <- log(model_data)
    
  }
  
  # subset by removing rows where age of diagnosis is missing
  # model_data <- model_data[complete.cases(model_data$age_diagnosis),]
  
  obs <- nrow(model_data)
  
  selected_features <- colnames(model_data)[3:ncol(model_data)]
  
  for (i in 1:iterations) {
    
    set.seed(i)
    train_index <- sample(nrow(model_data), nrow(model_data) *cutoff, replace = F)
    
    if(fac) {
      
      # determines how you train the model.
      NFOLDS <- 2
      fitControl <- trainControl( 
        method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
        number = min(10, NFOLDS),
        classProbs = TRUE,
        repeats = 1,
        allowParallel = TRUE,
        summaryFunction = twoClassSummary)
      
      y = make.names(as.factor(model_data$age_diagnosis_fac[train_index]))
      
    } else {
      
      NFOLDS <- 2
      fitControl <- trainControl( 
        method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
        number = min(10, NFOLDS),      
        repeats = 1,
        allowParallel = TRUE)
      
      y <- model_data$age_diagnosis[train_index]
      
    }
    
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
    
    if(fac){
      test.predictions[[i]] <- predict(model[[i]] 
                                       , newdata = model_data[-train_index, selected_features]
                                       , type = "prob")
      
      train.predictions[[i]] <- predict(model[[i]] 
                                        , newdata = model_data[train_index, selected_features]
                                        ,type = "prob")
      
    } else {
      test.predictions[[i]] <- predict(model[[i]] 
                                       , newdata = model_data[-train_index, selected_features])
      
      train.predictions[[i]] <- predict(model[[i]] 
                                        , newdata = model_data[train_index, selected_features])
      
    }
    
    
    if (fac) {
      
      train.ground_truth[[i]] <- as.factor(make.names(model_data$age_diagnosis_fac[train_index]))
      test.ground_truth[[i]] <- as.factor(make.names(model_data$age_diagnosis_fac[-train_index]))
      train.sample_collection[[i]] = as.factor(make.names(model_data$age_sample_fac[train_index]))
      train.sample_collection[[i]] <- factor(train.sample_collection[[i]], levels = c("X1", "X2"))
      # temp <- as.factor(make.names(data$age_sample_fac[-train_index]))
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
    } else {
      
      train.ground_truth[[i]] <- model_data$age_diagnosis[train_index]
      test.ground_truth[[i]] <- model_data$age_diagnosis[-train_index]
      train.sample_collection[[i]] = model_data$age_sample_collection[train_index]
      test.sample_collection[[i]] = model_data$age_sample_collection[-train_index]
      train.mse[[i]] <- rmse(unlist(train.predictions[[i]]), unlist(train.ground_truth[[i]]))
      test.mse[[i]] <- rmse(unlist(test.predictions[[i]]), unlist(test.ground_truth[[i]]))
      
      test.ground_truth
      
    }
    
    
    print(i)
    
  }
  
  return(list(train.mse, test.mse,  train.predictions, test.predictions, train.ground_truth, test.ground_truth, train.sample_collection,
              test.sample_collection, test_acc, test_stats, test_acc_samp, test_stats_samp, model, importance, obs))
  
}

# plot function 
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
  legend("bottomright", legend = paste0('# obs = ', result_list[[15]]), cex = 1, bty = 'n')
  
  
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
  legend("bottomright", legend = paste0('# obs = ', result_list[[15]]), cex = 1, bty = 'n')
  
  
  
}

###### read in dmr variations
# 
# # cancer and global
# bh_cancer, bh_global 
# # subselection
# bh_cancer_sub, bh_global_sub
# # balanced
# bh_cancer_bal, bh_global_bal
# # read in union and intersection
# bh_union, bh_intersection
# # read in bumphunter results from cancer indicator
# bh_cancer_ind 


###################
# Regression, not log, no residual
###################

# bh_global_bal with lsa_gene
bh_global_bal_lsa_gene <- predictAll(data = bh_global_bal,
                                     lsa_gene = T,
                                     knn_gene = F,
                                     knn_probe = F,
                                     threshold = nrow(bh_global_bal),
                                     fac = F,
                                     log = F,
                                     cutoff = .7,
                                     residual = F,
                                     iterations = 10)


plotModel(bh_global_bal_lsa_gene,
          xlim = c(0, 1000),
          ylim = c(0, 1000),
          main1 = 'bh_global_bal_lsa_gene',
          main2 = 'bh_global_bal_lsa_gene')


# bh_global_bal with knn_gene
bh_global_bal_knn_gene <- predictAll(data = bh_global_bal,
                                     lsa_gene = F,
                                     knn_gene = T,
                                     knn_probe = F,
                                     threshold = nrow(bh_global_bal),
                                     fac = F,
                                     log = F,
                                     cutoff = .7,
                                     residual = F,
                                     iterations = 10)


plotModel(bh_global_bal_knn_gene,
          xlim = c(0, 1000),
          ylim = c(0, 1000),
          main1 = 'bh_global_bal_knn_gene',
          main2 = 'bh_global_bal_knn_gene')


# bh_global_bal with knn_probe
bh_global_bal_knn_probe <- predictAll(data = bh_global_bal,
                                      lsa_gene = F,
                                      knn_gene = F,
                                      knn_probe = T,
                                      threshold = nrow(bh_global_bal),
                                      fac = F,
                                      log = F,
                                      cutoff = .7,
                                      residual = F,
                                      iterations = 10)


plotModel(bh_global_bal_knn_probe,
          xlim = c(0, 1000),
          ylim = c(0, 1000),
          main1 = 'Bump Hunter Features',
          main2 = 'Bump Hunter Features Sample Collection')
# 480 features
# get top features 

temp <- as.data.frame(do.call(rbind, bh_global_bal_knn_probe[[14]]))
colnames(temp) <- c('probe', 'importance')
temp$importance <- as.numeric(as.character(temp$importance))
temp <- temp[order(temp$importance, decreasing = T),]
temp <- temp[!duplicated(temp$probe),]

important_vars <- temp

# join with bh_global_bal to get gene info
location_info <- inner_join(important_vars, bh_global_bal, by = 'probe')

location_info <- location_info[with(location_info, order(-importance, p.value)), ]

# write.csv(location_info, '/home/benbrew/Desktop/location_info.csv')
# write.csv(location_info_resid, '/home/benbrew/Desktop/location_info_resid.csv')


###################
# Regression, not log, residual
###################

# bh_global_bal with lsa_gene
bh_global_bal_lsa_gene_resid <- predictAll(data = bh_global_bal,
                                           lsa_gene = T,
                                           knn_gene = F,
                                           knn_probe = F,
                                           threshold = nrow(bh_global_bal),
                                           fac = F,
                                           log = F,
                                           cutoff = .7,
                                           residual = T,
                                           iterations = 10)


plotModel(bh_global_bal_lsa_gene_resid,
          xlim = c(0, 1000),
          ylim = c(0, 1000),
          main1 = 'bh_global_bal_lsa_gene_resid',
          main2 = 'bh_global_bal_lsa_gene_resid')


# bh_global_bal with knn_gene
bh_global_bal_knn_gene_resid <- predictAll(data = bh_global_bal,
                                           lsa_gene = F,
                                           knn_gene = T,
                                           knn_probe = F,
                                           threshold = nrow(bh_global_bal),
                                           fac = F,
                                           log = F,
                                           cutoff = .7,
                                           residual = T,
                                           iterations = 10)


plotModel(bh_global_bal_knn_gene_resid,
          xlim = c(0, 1000),
          ylim = c(0, 1000),
          main1 = 'bh_global_bal_knn_gene_resid',
          main2 = 'bh_global_bal_knn_gene_resid')


# bh_global_bal with knn_probe
bh_global_bal_knn_probe_resid <- predictAll(data = bh_global_bal,
                                            lsa_gene = F,
                                            knn_gene = F,
                                            knn_probe = T,
                                            threshold = nrow(bh_global_bal),
                                            fac = F,
                                            log = F,
                                            cutoff = .7,
                                            residual = T,
                                            iterations = 10)


plotModel(bh_global_bal_knn_probe_resid,
          xlim = c(0, 1000),
          ylim = c(0, 1000),
          main1 = 'Bump Hunter Residual Features',
          main2 = 'bh_global_bal_knn_probe_resid')

temp <- as.data.frame(do.call(rbind, bh_global_bal_knn_probe_resid[[14]]))
colnames(temp) <- c('probe', 'importance')
temp$importance <- as.numeric(as.character(temp$importance))
temp <- temp[order(temp$importance, decreasing = T),]
temp <- temp[!duplicated(temp$probe),]

important_vars <- temp

# join with bh_global_bal to get gene info
location_info <- inner_join(important_vars, bh_global_bal, by = 'probe')

location_info_resid <- location_info

location_info_resid <- location_info_resid[with(location_info_resid, order(-importance, p.value)), ]

###################
# Regression, log, residual
###################

# bh_global_bal with lsa_gene
bh_global_bal_lsa_gene_resid_log <- predictAll(data = bh_global_bal,
                                               lsa_gene = T,
                                               knn_gene = F,
                                               knn_probe = F,
                                               threshold = nrow(bh_global_bal),
                                               fac = T,
                                               log = F,
                                               cutoff = .7,
                                               residual = T,
                                               iterations = 10)


plotModel(bh_global_bal_lsa_gene_resid_log,
          xlim = c(0, 10),
          ylim = c(0, 10),
          main1 = 'bh_global_bal_lsa_gene_resid_log',
          main2 = 'bh_global_bal_lsa_gene_resid_log')


# bh_global_bal with knn_gene
bh_global_bal_knn_gene_resid_log <- predictAll(data = bh_global_bal,
                                               lsa_gene = F,
                                               knn_gene = T,
                                               knn_probe = F,
                                               threshold = nrow(bh_global_bal),
                                               fac = F,
                                               log = T,
                                               cutoff = .7,
                                               residual = T,
                                               iterations = 10)


plotModel(bh_global_bal_knn_gene_resid_log,
          xlim = c(0, 10),
          ylim = c(0, 10),
          main1 = 'bh_global_bal_knn_gene_resid_log',
          main2 = 'bh_global_bal_knn_gene_resid_log')


# bh_global_bal with knn_probe
bh_global_bal_knn_probe_resid_log <- predictAll(data = bh_global_bal,
                                                lsa_gene = F,
                                                knn_gene = F,
                                                knn_probe = T,
                                                threshold = nrow(bh_global_bal),
                                                fac = F,
                                                log = T,
                                                cutoff = .7,
                                                residual = T,
                                                iterations = 10)


plotModel(bh_global_bal_knn_probe_resid,
          xlim = c(0, 10),
          ylim = c(0, 10),
          main1 = 'Bump Hunter Residual Features (Log)',
          main2 = 'bh_global_bal_knn_probe_resid_log')



############################################################################
# Same with factors as outcome
###########################################################################3

# bh_global with lsa_gene
bh_global_lsa_gene <- predictAll(data = bh_global,
                                     lsa_gene = T,
                                     knn_gene = F,
                                     knn_probe = F,
                                     threshold = nrow(bh_global),
                                     fac = T,
                                     log = F,
                                     cutoff = .7,
                                     residual = F,
                                     iterations = 10)

# test acc for age of diagnosis
mean(unlist(bh_global_lsa_gene[[9]]))

# test acc for age of sample collection
mean(unlist(bh_global_lsa_gene[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- bh_global_lsa_gene[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations


# bh_global with knn_gene
bh_global_knn_gene <- predictAll(data = bh_global,
                                     lsa_gene = F,
                                     knn_gene = T,
                                     knn_probe = F,
                                     threshold = nrow(bh_global),
                                     fac = T,
                                     log = F,
                                     cutoff = .7,
                                     residual = F,
                                     iterations = 10)


# test acc for age of diagnosis
mean(unlist(bh_global_knn_gene[[9]]))

# test acc for age of sample collection
mean(unlist(bh_global_knn_gene[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- bh_global_knn_gene[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations


# bh_global with knn_probe
bh_global_knn_probe <- predictAll(data = bh_global,
                                      lsa_gene = F,
                                      knn_gene = F,
                                      knn_probe = T,
                                      threshold = nrow(bh_global),
                                      fac = T,
                                      log = F,
                                      cutoff = .7,
                                      residual = F,
                                      iterations = 10)

# test acc for age of diagnosis
mean(unlist(bh_global_knn_probe[[9]]))

# test acc for age of sample collection
mean(unlist(bh_global_knn_probe[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- bh_global_knn_probe[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations



###################
# classification, not log, residual
###################

# bh_global with lsa_gene
bh_global_lsa_gene_resid <- predictAll(data = bh_global,
                                           lsa_gene = T,
                                           knn_gene = F,
                                           knn_probe = F,
                                           threshold = nrow(bh_global),
                                           fac = T,
                                           log = F,
                                           cutoff = .7,
                                           residual = T,
                                           iterations = 10)

# test acc for age of diagnosis
mean(unlist(bh_global_lsa_gene_resid[[9]]))

# test acc for age of sample collection
mean(unlist(bh_global_lsa_gene_resid[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- bh_global_lsa_gene_resid[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations


# bh_global with knn_gene
bh_global_knn_gene_resid <- predictAll(data = bh_global,
                                           lsa_gene = F,
                                           knn_gene = T,
                                           knn_probe = F,
                                           threshold = nrow(bh_global),
                                           fac = T,
                                           log = F,
                                           cutoff = .7,
                                           residual = T,
                                           iterations = 10)

# test acc for age of diagnosis
mean(unlist(bh_global_knn_gene_resid[[9]]))

# test acc for age of sample collection
mean(unlist(bh_global_knn_gene_resid[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- bh_global_knn_gene_resid[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations


# bh_global with knn_probe
bh_global_knn_probe_resid <- predictAll(data = bh_global,
                                            lsa_gene = F,
                                            knn_gene = F,
                                            knn_probe = T,
                                            threshold = nrow(bh_global),
                                            fac = T,
                                            log = F,
                                            cutoff = .7,
                                            residual = T,
                                            iterations = 10)

# test acc for age of diagnosis
mean(unlist(bh_global_knn_probe_resid[[9]]))

# test acc for age of sample collection
mean(unlist(bh_global_knn_probe_resid[[11]]))

# confustion matrix age of diagnosis 10
iterations <- 10
temp <- list()
for (i in 1:10){
  temp[[i]] <- bh_global_knn_probe_resid[[10]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)

mat_index <- seq(1, length(mat), 4)

new_mat[1,1] <- sum(mat[mat_index])/iterations
new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
new_mat[2,2] <- sum(mat[mat_index + 3])/iterations


