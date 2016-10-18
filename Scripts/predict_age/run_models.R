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

# set parameters 
data_thresholds <- c(48, 60, 72, 84)

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
                      random = F,
                      bump_hunter = F,
                      bump_hunter_data) {
  
  # get differenct variations of data
  data <- subsetDat(data)
  
  if (bump_hunter) {
    
    data <- bhSubset(data, bh_data = bump_hunter_data)

  }
  
  if (random) {
    
    data <- getRand(data)
    
  }
  
  data_resid <- getResidual(data)
  
  data_fac <- list()
  data_resid_fac <- list()
  
  for (thresh in 1:length(data_thresholds)) {
    
    data_fac[[thresh]] <- makeFac(data, threshold = data_thresholds[thresh])
    data_resid_fac[[thresh]] <- makeFac(data_resid, threshold = data_thresholds[thresh])
    
  }
  
  # run regressions
  data_result <- list()
  data_resid_result <- list()
 
  for (dat in 1:length(data)) {
    
    sub_dat <- data[[dat]]
    sub_data_resid <- data_resid[[dat]]
    data_result[[dat]] <- rfPredictReg(sub_dat, cutoff = .7, iterations = 10)
    data_resid_result[[dat]] <- rfPredictReg(sub_data_resid, cutoff = .7, iterations = 10)
  }
  
  # run classification
  data_fac_result <- list()
  data_resid_fac_result <- list()
  
  for(dat in 1:length(data_fac)) {
    
    sub_dat_fac <- data_fac[[dat]]
    sub_dat_resid_fac <- data_resid_fac[[dat]]
    
    temp.data_fac_result <- list()
    temp.data_resid_fac_result <- list()
    
    for(sub_dat in 1:length(sub_dat_fac)) {
      
      temp.sub_dat_fac <- sub_dat_fac[[sub_dat]]
      temp.sub_dat_resid_fac <- sub_dat_resid_fac[[sub_dat]]
      
      temp.data_fac_result[[sub_dat]] <- rfPredictFac(temp.sub_dat_fac, cutoff = .7, iterations = 10)
      temp.data_resid_fac_result[[sub_dat]] <- rfPredictFac(temp.sub_dat_resid_fac, cutoff = .7, iterations = 10)
      
    }
    
    data_fac_result[[dat]] <- temp.data_fac_result
    data_resid_fac_result[[dat]] <- temp.data_fac_result
    
  }
  
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


###################################
# second run each probe data - probe_knn, probe_lsa, with fac
###################################

###### probe KNN
probe_knn_models <- runModels(probe_knn, bump_hunter = F)
# plot models 
plotModel(probe_knn_models[[1]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe KNN Diagnosis ',
          main2 = 'Probe KNN Sample' )

plotModel(probe_knn_models[[2]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe KNN Diagnosis (Resid) ',
          main2 = 'Probe KNN Sample (Resid)' )

# get confusion matrix 
conMatrix(probe_knn_models[[3]])
conMatrix(probe_knn_models[[4]])

###### probe lsa
probe_lsa_models <- runModels(probe_lsa, bump_hunter = F)
# plot models 
plotModel(probe_lsa_models[[1]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe lsa Diagnosis ',
          main2 = 'Probe lsa Sample' )

plotModel(probe_lsa_models[[2]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe lsa Diagnosis (Resid) ',
          main2 = 'Probe lsa Sample (Resid)' )

# get confusion matrix 
conMatrix(probe_lsa_models[[3]])
conMatrix(probe_lsa_models[[4]])

###################################
# Third, run probe_knn and probe_lsa with all different bumphunter features
###################################

###### probe knn with global features
probe_knn_global_models <- runModels(probe_knn, random = F, bump_hunter = T, 
                                     bump_hunter_data = bh_probe_knn_global_features)
# plot models 
plotModel(probe_knn_global_models[[1]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe knn global Diagnosis ',
          main2 = 'Probe knn global Sample' )

plotModel(probe_knn_global_models[[2]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe lsa Diagnosis (Resid) ',
          main2 = 'Probe lsa Sample (Resid)' )

# get confusion matrix 
conMatrix(probe_knn_global_models[[3]])
conMatrix(probe_knn_global_models[[4]])

# dims 
probe_knn_global_models[[1]][15]
probe_knn_global_models[[2]][15]
probe_knn_global_models[[3]][13]
probe_knn_global_models[[4]][13]

###### probe knn with cancer features
probe_knn_cancer_models <- runModels(probe_knn, bump_hunter = T, 
                                     bump_hunter_data = bh_probe_knn_cancer_features)
# plot models 
plotModel(probe_knn_cancer_models[[1]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe knn cancer Diagnosis ',
          main2 = 'Probe knn cancer Sample' )

plotModel(probe_knn_cancer_models[[2]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe lsa Diagnosis (Resid) ',
          main2 = 'Probe lsa Sample (Resid)' )

# get confusion matrix 
conMatrix(probe_knn_cancer_models[[3]])
conMatrix(probe_knn_cancer_models[[4]])

# dims 
probe_knn_cancer_models[[1]][15]
probe_knn_cancer_models[[2]][15]
probe_knn_cancer_models[[3]][13]
probe_knn_cancer_models[[4]][13]


###### probe lsa with global features
probe_lsa_global_models <- runModels(probe_lsa, bump_hunter = T, 
                                     bump_hunter_data = bh_probe_lsa_global_features)
# plot models 
plotModel(probe_lsa_global_models[[1]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe lsa global Diagnosis ',
          main2 = 'Probe lsa global Sample' )

plotModel(probe_lsa_global_models[[2]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe lsa global Diagnosis (Resid) ',
          main2 = 'Probe lsa global Sample (Resid)' )

# get confusion matrix 
conMatrix(probe_lsa_global_models[[3]])
conMatrix(probe_lsa_global_models[[4]])

# dims 
probe_lsa_global_models[[1]][15]
probe_lsa_global_models[[2]][15]
probe_lsa_global_models[[3]][13]
probe_lsa_global_models[[4]][13]

###### probe lsa with cancer features
probe_lsa_cancer_models <- runModels(probe_lsa, bump_hunter = T, 
                                     bump_hunter_data = bh_probe_lsa_cancer_features)
# plot models 
plotModel(probe_lsa_cancer_models[[1]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe lsa cancer Diagnosis ',
          main2 = 'Probe lsa cancer Sample' )

plotModel(probe_lsa_cancer_models[[2]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe lsa cancer Diagnosis (Resid) ',
          main2 = 'Probe lsa cancer Sample (Resid)' )

# get confusion matrix 
conMatrix(probe_lsa_cancer_models[[3]])
conMatrix(probe_lsa_cancer_models[[4]])

# dims 
probe_lsa_cancer_models[[1]][15]
probe_lsa_cancer_models[[2]][15]
probe_lsa_cancer_models[[3]][13]
probe_lsa_cancer_models[[4]][13]


###################################
# Finally run each gene and probe data with random features
###################################

###### GENE KNN
gene_knn_rand <- runModels(gene_knn, random = T)
# plot rand 
plotModel(gene_knn_rand[[1]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Gene KNN Rand Diagnosis ',
          main2 = 'Gene KNN Rand Sample' )

plotModel(gene_knn_rand[[2]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Gene KNN Rand Diagnosis (Resid) ',
          main2 = 'Gene KNN Rand Sample (Resid)' )

# get confusion matrix 
conMatrix(gene_knn_rand[[3]])
conMatrix(gene_knn_rand[[4]])

# dims 
gene_knn_rand[[1]][15]
gene_knn_rand[[2]][15]
gene_knn_rand[[3]][13]
gene_knn_rand[[4]][13]



###### GENE lsa
gene_lsa_rand <- runModels(gene_lsa, random = T)
# plot rand 
plotModel(gene_lsa_rand[[1]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Gene lsa Rand Diagnosis ',
          main2 = 'Gene lsa Rand Sample' )

plotModel(gene_lsa_rand[[2]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Gene lsa Rand Diagnosis (Resid) ',
          main2 = 'Gene lsa Rand Sample (Resid)' )

# get confusion matrix 
conMatrix(gene_lsa_rand[[3]])
conMatrix(gene_lsa_rand[[4]])

# dims 
gene_lsa_rand[[1]][15]
gene_lsa_rand[[2]][15]
gene_lsa_rand[[3]][13]
gene_lsa_rand[[4]][13]


###### probe KNN
probe_knn_rand <- runModels(probe_knn, random = T)
# plot rand 
plotModel(probe_knn_rand[[1]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe KNN Rand Diagnosis ',
          main2 = 'Probe KNN Rand  Sample' )

plotModel(probe_knn_rand[[2]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe KNN Rand Diagnosis (Resid) ',
          main2 = 'Probe KNN Rand Sample (Resid)' )

# get confusion matrix 
conMatrix(probe_knn_rand[[3]])
conMatrix(probe_knn_rand[[4]])

# dims 
probe_knn_rand[[1]][15]
probe_knn_rand[[2]][15]
probe_knn_rand[[3]][13]
probe_knn_rand[[4]][13]


###### probe lsa
probe_lsa_rand <- runModels(probe_lsa, random = T)
# plot rand 
plotModel(probe_lsa_rand[[1]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe lsa Rand Diagnosis ',
          main2 = 'Probe lsa Rand Sample' )

plotModel(probe_lsa_rand[[2]], xlim = c(0, 1000), ylim = c(0,1000), main1 = 'Probe lsa Rand Diagnosis (Resid) ',
          main2 = 'Probe lsa Rand Sample (Resid)' )

# get confusion matrix 
conMatrix(probe_lsa_rand[[3]])
conMatrix(probe_lsa_rand[[4]])

# dims 
probe_lsa_rand[[1]][15]
probe_lsa_rand[[2]][15]
probe_lsa_rand[[3]][13]
probe_lsa_rand[[4]][13]

# # Save main model data (use of all features)
# save.image(paste0(model_data, '/model_results_all_features.RData'))
# load(paste0(model_data, '/model_results_all_features.RData'))




