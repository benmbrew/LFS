
#!/hpf/tools/centos6/R/3.2.3/bin/Rscript

# Script to evaluate the how each imputation method affects the
# performance of the clustering methods

# argv <- as.numeric(commandArgs(T))

##########
# This script will predict cancer with cases and controls, both 450k and 850k

# predict cancer = methylation (p53 mutants)

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
library(e1071)

registerDoParallel(1)
##########
# initialize folders
##########
home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'raw'
k = 5
seed_num = 1
##########
# load data
##########
# read in full m value data 
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new_m.rda'))) #34 449936
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new_m.rda')))
#35 449783

###########
# make id into ids
###########
colnames(betaCases)[1] <- 'ids'
colnames(betaControls)[1] <- 'ids'
colnames(betaValid)[1] <- 'ids'

##########
# remove inf
##########
betaCases <- removeInf(betaCases, probe_start = 8)
betaControls <- removeInf(betaControls, probe_start = 8)
betaValid<- removeInf(betaValid, probe_start = 8)


# get old controls - Mut and 'Unaffected'
betaControlsOld <- subset(betaCases, p53_germline == 'Mut' & 
                            cancer_diagnosis_diagnoses == 'Unaffected')

# get p53, not 'Unaffected'
betaCases <- subset(betaCases, p53_germline == 'Mut' & 
                      cancer_diagnosis_diagnoses != 'Unaffected')

# get rid of cancer samples in controls 
betaControls <- betaControls[grepl('Unaffected', betaControls$cancer_diagnosis_diagnoses),]

#subset valid
betaValid <- betaValid[!betaValid$ids %in% betaCases$ids,]

##########
# get intersecting colnames and prepare data for modeling
##########

intersect_names <- Reduce(intersect, list(colnames(betaCases)[8:ncol(betaCases)], 
                                          colnames(betaControls)[8:ncol(betaControls)], 
                                          colnames(betaValid)[8:ncol(betaValid)]))
# assign dataframe identifier
betaCases$type <- '450k'
betaControls$type <- '850k'
betaControlsOld$type <- '450k'
betaValid$type <- '850k'


# cases
betaCases <- betaCases[, c('ids',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender',
                           'type',
                           intersect_names)]
# controls
betaControls <- betaControls[, c('ids',
                                 'age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender', 
                                 'type',
                                 intersect_names)]

# controls
betaControlsOld <- betaControlsOld[, c('ids',
                                 'age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender', 
                                 'type',
                                 intersect_names)]

#validation
betaValid <- betaValid[, c('ids', 
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           'type',
                           intersect_names)]



###########
# get full data
###########
beta_full <- rbind(betaCases,
                   betaControls,
                   betaControlsOld,
                   betaValid)

rm(betaCases,
   betaControls,
   betaControlsOld,
   betaValid)

##########
# remove duplicates 
##########
length(which(duplicated(beta_full$ids)))

# remove 
beta_full <- beta_full[!duplicated(beta_full$ids),]

# summary of data types 
summary(as.factor(beta_full$type))

# get a column for each dataset indicating the fold
beta_full <- getFolds(beta_full, seed_number = seed_num, k = k)


# lets have this model either do 450 only, 850 only, or both, no combinations
trainTest <- function(dat,
                      methyl_tech,
                      k) 
{
  
  
  if(methyl_tech == '450k'){
    
    beta_dat <- dat[grepl('450k', dat$type),]
    
  } else if (methyl_tech =='850k') {
    
    beta_dat <- dat[grepl('850k', dat$type),]
    
  } else if(methyl_tech == '450_850') {
    
    beta_dat <- subset(dat, type == '450k' |  (type == '850k' & cancer_diagnosis_diagnoses == 'Unaffected'))
    
  } else if (methyl_tech =='850_450') {
    
    beta_dat <- subset(dat, type == '850k' |  (type == '450k' & cancer_diagnosis_diagnoses == 'Unaffected'))
    
  } else {
    
    beta_dat <- dat
    
  }
  
  
  # list to store results
  bh_feat <- list()
  model_results <- list()
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, beta_dat$folds)
    test_index <- !train_index
    
    mod_feats <- colnames(beta_dat)[7:(ncol(beta_dat) -1)]
    
    mod_feats <- sample(mod_feats, 5000, replace = T)
    
    mod_result <- predCancer(training_dat = beta_dat[train_index,], 
                             test_dat = beta_dat[test_index,], 
                             bh_features = mod_feats)
    
    
    model_results[[i]] <- getResultsCancer(mod_result)
    
    
  }
  
  return(model_results)
  
}
set.seed(4)
mod_results <- trainTest(dat = beta_full, methyl_tech = '850', k = 5)


getClassResutls <- function(temp.result) {
  
  # creat matrix
  mat <- matrix(NA, nrow = 2, ncol = 2)

  for (j in 1:4) {
    
    if(j > 1) {
      
      mat <- mat + temp.result[[j]][[4]]$table[1:2,1:2]
      acc_temp <- acc_temp + temp.result[[j]][[4]]$overall[1]
      tpr_temp <- tpr_temp + temp.result[[j]][[4]]$byClass[1]
      tnr_temp <- tnr_temp + temp.result[[j]][[4]]$byClass[2]
      alpha_temp <- alpha_temp + temp.result[[j]][[1]]
      lambda_temp <- lambda_temp + temp.result[[j]][[1]]
      
      
    } else{
      
      mat <- temp.result[[j]][[4]]$table[1:2,1:2]
      acc_temp <- as.numeric(temp.result[[j]][[4]]$overall[1])
      tpr_temp <- as.numeric(temp.result[[j]][[4]]$byClass[1])
      tnr_temp <- as.numeric(temp.result[[j]][[4]]$byClass[2])
      alpha_temp <- as.numeric(temp.result[[j]][[1]])
      lambda_temp <- as.numeric(temp.result[[j]][[1]])
      
    }
    
  }
  
  # get average from j
  confusion_table <- mat/4
  acc_table <- acc_temp/4
  tnr_table <- tnr_temp/4
  tpr_table <- tpr_temp/4
  alpha_table <- alpha_temp/4
  lambda_table <- lambda_temp/4
  
  # store at i 
  table_results <- confusion_table
  
  table_results <- as.data.frame(table_results)
  
  table_scores <- as.data.frame(cbind(alpha_table, lambda_table, acc_table, tnr_table, tpr_table))
  
  # remove NAs - bad predictions dont yield a value for precision, etc
  table_scores <- table_scores[complete.cases(table_scores),]
  
  # fpr, fnr
  table_scores$fpr <- 1 - table_scores$tnr
  table_scores$fnr <- 1 - table_scores$tpr
  

  return(list(table_results, table_scores))
}

# get results for age 48
result_48 <- getClassResutls(mod_results)
conMat_48 <- result_48[[1]]
resultTab_48 <- result_48[[2]]
