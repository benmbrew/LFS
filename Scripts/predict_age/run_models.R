##################################################################3
# this script will source model_functions.R and run models
##################################################################

library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(doParallel)
library(nnet)

registerDoParallel(1)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
scripts_folder <- paste0(project_folder, '/Scripts')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')
model_data <- paste0(data_folder, '/model_data')


# source model_functions to get functions to run models 
source(paste0(scripts_folder, '/predict_age/model_functions.R'))


###########################################
# Read in data- gene_knn, gene_lsa, 
# probe_knn, probe_lsa, and bh_features
###########################################

load(paste0(model_data, '/model_data.RData'))
load(paste0(model_data, '/bh_features.RData'))
rm(cg_locations)

# Data types: 
# 1) full_data : gene_knn, gene_lsa, probe_knn, probe_lsa 
# 2) bh features: bh_probe_knn_cancer_features, bh_probe_lsa_cancer_features,
#    bh_probe_knn_global_features, bh_probe_lsa_global_features

# function to run models - subset, get residuals, get categorical, predict with regression and fac. 
runModels <- function(data,
                      bump_hunter) {
  
  # get differenct variations of data
  data <- subsetDat(data)
  data_resid <- getResidual(data)
  
  if (bump_hunter) {
    
    data <- bhSubset(data)
    data_resid <- bhSubset(data_resid)
    
  }
  
  data_fac <- makeFac(data, threshold = 48)
  data_resid_fac <- makeFac(data_resid, threshold = 48)
  
  # run models 
  data_result <- rfPredictReg(data, cutoff = .7, iterations = 10)
  data_resid_result <- rfPredictReg(data_resid, cutoff = .7, iterations = 10)
  data_fac_result <- rfPredictFac(data_fac, cutoff = .7, iterations = 10)
  data_resid_fac_result <- rfPredictFac(data_resid_fac, cutoff = .7, iterations = 10)
  
  return(list(data_result, data_resid_result, data_fac_result, data_resid_fac_result))
  
}

###################################
# first run each gene data - gene_knn, gene_lsa, with fac
###################################

###### GENE KNN
gene_knn_models <- runModels(gene_knn, bump_hunter = F)
# plot models 
plotModel(gene_knn_models[[1]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Gene KNN Diagnosis ',
          main2 = 'Gene KNN Sample' )

plotModel(gene_knn_models[[2]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Gene KNN Diagnosis (Resid) ',
          main2 = 'Gene KNN Sample (Resid)' )

# get confusion matrix 
conMatrix(gene_knn_models[[3]])
conMatrix(gene_knn_models[[4]])

###### GENE lsa
gene_lsa_models <- runModels(gene_lsa, bump_hunter = F)
# plot models 
plotModel(gene_lsa_models[[1]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Gene lsa Diagnosis ',
          main2 = 'Gene lsa Sample' )

plotModel(gene_lsa_models[[2]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Gene lsa Diagnosis (Resid) ',
          main2 = 'Gene lsa Sample (Resid)' )

# get confusion matrix 
conMatrix(gene_lsa_models[[3]])
conMatrix(gene_lsa_models[[4]])






























# get differenct variations of data
gene_knn <- subsetDat(gene_knn)
gene_knn_resid <- getResidual(gene_knn)

gene_knn_fac <- makeFac(gene_knn, threshold = 48)
gene_knn_resid_fac <- makeFac(gene_knn_resid, threshold = 48)

# run models 
gene_knn_result <- rfPredictReg(gene_knn, cutoff = .7, iterations = 10)
gene_knn_resid_result <- rfPredictReg(gene_knn_resid, cutoff = .7, iterations = 10)
gene_knn_fac_result <- rfPredictFac(gene_knn_fac, cutoff = .7, iterations = 10)
gene_knn_resid_fac_result <- rfPredictFac(gene_knn_resid_fac, cutoff = .7, iterations = 10)

# plot models 
plotModel(gene_knn_result, xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Gene KNN Diagnosis ',
          main2 = 'Gene KNN Sample' )

plotModel(gene_knn_resid_result, xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Gene KNN Diagnosis (Resid) ',
          main2 = 'Gene KNN Sample (Resid)' )

# get confusion matrix 
conMatrix(gene_knn_fac_result)
conMatrix(gene_knn_resid_fac_result)



