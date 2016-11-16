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
idat_data <- paste0(methyl_data, '/raw_files')
clin_data <- paste0(data_folder, '/clin_data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')
model_data <- paste0(data_folder, '/model_data')


# source model_functions to get functions to run models 
source(paste0(scripts_folder, '/predict_age/model_functions.R'))

# set parameters 
data_thresholds <- c(48, 60, 72, 84)
p53 <- c('Mut', 'WT')

###########################################
# Read in data- gene_knn, gene_lsa, 
# probe_knn, probe_lsa, and bh_features
###########################################

load(paste0(idat_data, '/imputed_idat_betas_final.RData'))
load(paste0(model_data, '/bh_features_idat.RData'))

# load(paste0(model_data, '/model_data.RData'))
# load(paste0(model_data, '/bh_features.RData'))

# Data types: 
# 1) full_data : gene_knn, gene_lsa, probe_knn, probe_lsa 
# 2) bh features_balanced: bh_probe_knn_cancer_features, bh_probe_lsa_cancer_features,
#    bh_probe_knn_global_features, bh_probe_lsa_global_features
# 3) bh features_unbalanced: bh_probe_knn_cancer_unbal_features, bh_probe_lsa_cancer_unbal_features,
#    bh_probe_knn_global_unbal_features, bh_probe_lsa_global_unbal_features
# 4) random with each one

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
    data_resid_fac_result[[dat]] <- temp.data_resid_fac_result
    
  }
  
  return(list(data_result, data_resid_result, data_fac_result, data_resid_fac_result))
  
}

##################################
# beta_raw_
##################################
# cancer balanced
beta_raw_cancer_bal_models <- runModels(beta_raw, random = F, bump_hunter = T, 
                              bump_hunter_data = beta_raw_cancer_bal_features)

beta_raw_cancer_bal_table <- extractResults(beta_raw_cancer_bal_models, data_name = 'beta_raw_cancer_bal')

#global balanced
beta_raw_global_bal_models <- runModels(beta_raw, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_raw_global_bal_features)

beta_raw_global_bal_table <- extractResults(beta_raw_cancer_bal_models, data_name = 'beta_raw_global_bal')

# cancer unbalanced
beta_raw_cancer_unbal_models <- runModels(beta_raw, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_raw_cancer_unbal_features)

beta_raw_cancer_unbal_table <- extractResults(beta_raw_cancer_unbal_models, data_name = 'beta_raw_cancer_unbal')

#global unbalanced
beta_raw_global_unbal_models <- runModels(beta_raw, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_raw_global_unbal_features)

beta_raw_global_unbal_table <- extractResults(beta_raw_cancer_unbal_models, data_name = 'beta_raw_global_unbal')


##################################
# beta_swan_
##################################
# cancer balanced
beta_swan_cancer_bal_models <- runModels(beta_swan, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_swan_cancer_bal_features)

beta_swan_cancer_bal_table <- extractResults(beta_swan_cancer_bal_models, data_name = 'beta_swan_cancer_bal')

#global balanced
beta_swan_global_bal_models <- runModels(beta_swan, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_swan_global_bal_features)

beta_swan_global_bal_table <- extractResults(beta_swan_cancer_bal_models, data_name = 'beta_swan_global_bal')

# cancer unbalanced
beta_swan_cancer_unbal_models <- runModels(beta_swan, random = F, bump_hunter = T, 
                                          bump_hunter_data = beta_swan_cancer_unbal_features)

beta_swan_cancer_unbal_table <- extractResults(beta_swan_cancer_unbal_models, data_name = 'beta_swan_cancer_unbal')

#global unbalanced
beta_swan_global_unbal_models <- runModels(beta_swan, random = F, bump_hunter = T, 
                                          bump_hunter_data = beta_swan_global_unbal_features)

beta_swan_global_unbal_table <- extractResults(beta_swan_cancer_unbal_models, data_name = 'beta_swan_global_unbal')


##################################
# beta_quan_
##################################
# cancer balanced
beta_quan_cancer_bal_models <- runModels(beta_quan, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_quan_cancer_bal_features)

beta_quan_cancer_bal_table <- extractResults(beta_quan_cancer_bal_models, data_name = 'beta_quan_cancer_bal')

#global balanced
beta_quan_global_bal_models <- runModels(beta_quan, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_quan_global_bal_features)

beta_quan_global_bal_table <- extractResults(beta_quan_cancer_bal_models, data_name = 'beta_quan_global_bal')

# cancer unbalanced
beta_quan_cancer_unbal_models <- runModels(beta_quan, random = F, bump_hunter = T, 
                                          bump_hunter_data = beta_quan_cancer_unbal_features)

beta_quan_cancer_unbal_table <- extractResults(beta_quan_cancer_unbal_models, data_name = 'beta_quan_cancer_unbal')

#global unbalanced
beta_quan_global_unbal_models <- runModels(beta_quan, random = F, bump_hunter = T, 
                                          bump_hunter_data = beta_quan_global_unbal_features)

beta_quan_global_unbal_table <- extractResults(beta_quan_cancer_unbal_models, data_name = 'beta_quan_global_unbal')


##################################
# beta_funnorm_
##################################
# cancer balanced
beta_funnorm_cancer_bal_models <- runModels(beta_funnorm, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_funnorm_cancer_bal_features)

beta_funnorm_cancer_bal_table <- extractResults(beta_funnorm_cancer_bal_models, data_name = 'beta_funnorm_cancer_bal')

#global balanced
beta_funnorm_global_bal_models <- runModels(beta_funnorm, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_funnorm_global_bal_features)

beta_funnorm_global_bal_table <- extractResults(beta_funnorm_cancer_bal_models, data_name = 'beta_funnorm_global_bal')

# cancer unbalanced
beta_funnorm_cancer_unbal_models <- runModels(beta_funnorm, random = F, bump_hunter = T, 
                                          bump_hunter_data = beta_funnorm_cancer_unbal_features)

beta_funnorm_cancer_unbal_table <- extractResults(beta_funnorm_cancer_unbal_models, data_name = 'beta_funnorm_cancer_unbal')

#global unbalanced
beta_funnorm_global_unbal_models <- runModels(beta_funnorm, random = F, bump_hunter = T, 
                                          bump_hunter_data = beta_funnorm_global_unbal_features)

beta_funnorm_global_unbal_table <- extractResults(beta_funnorm_cancer_unbal_models, data_name = 'beta_funnorm_global_unbal')

#####################################
# Now run each beta with random features
####################################
beta_raw_rand <- runModels(beta_raw, random = T, bump_hunter = F)
beta_raw_rand_table <- extractResults(beta_raw_rand, data_name = 'beta_raw_rand')

beta_swan_rand <- runModels(beta_swan, random = T, bump_hunter = F)
beta_swan_rand_table <- extractResults(beta_swan_rand, data_name = 'beta_swan_rand')

beta_quan_rand <- runModels(beta_quan, random = T, bump_hunter = F)
beta_quan_rand_table <- extractResults(beta_quan_rand, data_name = 'beta_quan_rand')

beta_funnorm_rand <- runModels(beta_quan, random = T, bump_hunter = F)
beta_funnorm_rand_table <- extractResults(beta_funnorm_rand, data_name = 'beta_funnorm_rand')

#######################################
# rbind final results
#######################################
beta_idat_results <- rbind(beta_raw_cancer_bal_table, beta_raw_cancer_unbal_table, beta_raw_global_bal_table, beta_raw_global_unbal_table,
                           beta_swan_cancer_bal_table, beta_swan_cancer_unbal_table, beta_swan_global_bal_table, beta_swan_global_unbal_table,
                           beta_quan_cancer_bal_table, beta_quan_cancer_unbal_table, beta_quan_global_bal_table, beta_quan_global_unbal_table,
                           beta_funnorm_cancer_bal_table, beta_funnorm_cancer_unbal_table, beta_funnorm_global_bal_table, beta_funnorm_global_unbal_table,
                           beta_raw_rand_table, beta_swan_rand_table, beta_quan_rand_table, beta_funnorm_rand_table)


rm(beta_raw_cancer_bal_features, beta_raw_cancer_bal_models, beta_raw_cancer_unbal_features, beta_raw_cancer_unbal_models,
   beta_raw_global_bal_features, beta_raw_global_bal_models, beta_raw_global_unbal_features, beta_raw_global_unbal_models,
   beta_swan_cancer_bal_features, beta_swan_cancer_bal_models, beta_swan_cancer_unbal_features, beta_swan_cancer_unbal_models,
   beta_swan_global_bal_features, beta_swan_global_bal_models, beta_swan_global_unbal_features, beta_swan_global_unbal_models,
   beta_quan_cancer_bal_features, beta_quan_cancer_bal_models, beta_quan_cancer_unbal_features, beta_quan_cancer_unbal_models,
   beta_quan_global_bal_features, beta_quan_global_bal_models, beta_quan_global_unbal_features, beta_quan_global_unbal_models,
   beta_funnorm_cancer_bal_features, beta_funnorm_cancer_bal_models, beta_funnorm_cancer_unbal_features, beta_funnorm_cancer_unbal_models,
   beta_funnorm_global_bal_features, beta_funnorm_global_bal_models, beta_funnorm_global_unbal_features, beta_funnorm_global_unbal_models,
   beta_raw, beta_quan, beta_swan, beta_funnorm, beta_funnorm_rand, beta_raw_rand, beta_quan_rand, beta_swan_rand)

load(paste0(model_data, '/idat_beta_table_results.RData'))

#####################################
# # first run each gene data - gene_knn, gene_lsa, with fac
# ###################################
# 
# ###### GENE KNN
# gene_knn_models <- runModels(gene_knn, bump_hunter = F)
# # get result table 
# gene_knn_table <- extractResults(gene_knn_models, data_name = 'gene knn all features')
# 
# ###### GENE lsa
# gene_lsa_models <- runModels(gene_lsa, bump_hunter = F)
# # get result table 
# gene_lsa_table <- extractResults(gene_lsa_models, data_name = 'gene lsa all features')
# 
# 
# # # Save main model data (use of all features)
# # save.image(paste0(model_data, '/model_results_all_features.RData'))
# # load(paste0(model_data, '/model_results_all_features.RData'))
# ###################################
# # second run each probe data - probe_knn, probe_lsa, with fac
# ###################################
# 
# ###### probe KNN
# probe_knn_models <- runModels(probe_knn, bump_hunter = F)
# # get result table 
# probe_knn_table <- extractResults(probe_knn_models, data_name = 'probe knn all features')
# 
# ###### probe lsa
# probe_lsa_models <- runModels(probe_lsa, bump_hunter = F)
# # get result table 
# probe_lsa_table <- extractResults(probe_lsa_models, data_name = 'probe lsa all features')
# 
# ###################################
# # Third, run probe_knn and probe_lsa with all different bumphunter features
# ###################################
# 
# ###### probe knn with global features
# probe_knn_global_models <- runModels(probe_knn, random = F, bump_hunter = T, 
#                                      bump_hunter_data = bh_probe_knn_global_features)
# # get result table 
# probe_knn_global_table <- extractResults(probe_knn_global_models, data_name = 'probe knn global')
# 
# ###### probe knn with cancer features
# probe_knn_cancer_models <- runModels(probe_knn, bump_hunter = T, 
#                                      bump_hunter_data = bh_probe_knn_cancer_features)
# # get result table 
# probe_knn_cancer_table <- extractResults(probe_knn_cancer_models, data_name = 'probe knn cancer')
# 
# ###### probe lsa with global features
# probe_lsa_global_models <- runModels(probe_lsa, bump_hunter = T, 
#                                      bump_hunter_data = bh_probe_lsa_global_features)
# # get result table 
# probe_lsa_global_table <- extractResults(probe_lsa_global_models, data_name = 'probe lsa global')
# 
# ###### probe lsa with cancer features
# probe_lsa_cancer_models <- runModels(probe_lsa, bump_hunter = T, 
#                                      bump_hunter_data = bh_probe_lsa_cancer_features)
# # get result table 
# probe_lsa_cancer_table <- extractResults(probe_lsa_cancer_models, data_name = 'probe lsa cancer')
# 
# ###### probe knn union 
# probe_knn_union_models <- runModels(probe_knn, random = F, bump_hunter = T, 
#                                     bump_hunter_data = bh_union_bal)
# # get result table 
# probe_knn_union_table <- extractResults(probe_knn_union_models, data_name = 'probe knn union')
# 
# ###### probe knn intersection
# probe_knn_int_models <- runModels(probe_knn, random = F, bump_hunter = T, 
#                                     bump_hunter_data = bh_intersection_bal)
# # get result table 
# probe_knn_int_table <- extractResults(probe_knn_int_models, data_name = 'probe knn intersection')
# 
# ###### probe lsa union 
# probe_lsa_union_models <- runModels(probe_lsa, random = F, bump_hunter = T, 
#                                     bump_hunter_data = bh_union_bal)
# # get result table 
# probe_lsa_union_table <- extractResults(probe_lsa_union_models, data_name = 'probe lsa union')
# 
# ###### probe lsa intersection
# probe_lsa_int_models <- runModels(probe_lsa, random = F, bump_hunter = T, 
#                                   bump_hunter_data = bh_intersection_bal)
# # get result table 
# probe_lsa_int_table <- extractResults(probe_lsa_int_models, data_name = 'probe lsa intersection')
# 
# ###################################
# # Fourth, run probe_knn and probe_lsa with all different bumphunter features that are unbalanced 
# ###################################
# 
# ###### probe knn with global features
# probe_knn_global_unbal_models <- runModels(probe_knn, random = F, bump_hunter = T, 
#                                      bump_hunter_data = bh_probe_knn_global_unbal_features)
# # get result table 
# probe_knn_global_unbal_table <- extractResults(probe_knn_global_unbal_models, data_name = 'probe knn global unbalanced')
# 
# ###### probe knn with cancer features
# probe_knn_cancer_unbal_models <- runModels(probe_knn, bump_hunter = T, 
#                                      bump_hunter_data = bh_probe_knn_cancer_unbal_features)
# # get result table 
# probe_knn_cancer_unbal_table <- extractResults(probe_knn_cancer_unbal_models, data_name = 'probe knn cancer unbalanced')
# 
# ###### probe lsa with global features
# probe_lsa_global_unbal_models <- runModels(probe_lsa, bump_hunter = T, 
#                                      bump_hunter_data = bh_probe_lsa_global_unbal_features)
# # get result table 
# probe_lsa_global_unbal_table <- extractResults(probe_lsa_global_unbal_models, data_name = 'probe lsa global unbalanced')
# 
# ###### probe lsa with cancer features
# probe_lsa_cancer_unbal_models <- runModels(probe_lsa, bump_hunter = T, 
#                                      bump_hunter_data = bh_probe_lsa_cancer_unbal_features)
# # get result table 
# probe_lsa_cancer_unbal_table <- extractResults(probe_lsa_cancer_unbal_models, data_name = 'probe lsa cancer unbalanced')
# 
# ###### probe knn union 
# probe_knn_union_unbal_models <- runModels(probe_knn, random = F, bump_hunter = T, 
#                                     bump_hunter_data = bh_union_unbal)
# # get result table 
# probe_knn_union_unbal_table <- extractResults(probe_knn_union_unbal_models, data_name = 'probe knn union unbalanced')
# 
# ###### probe knn intersection
# probe_knn_int_unbal_models <- runModels(probe_knn, random = F, bump_hunter = T,  
#                                   bump_hunter_data = bh_intersection_unbal)
# # get result table 
# probe_knn_int_unbal_table <- extractResults(probe_knn_int_unbal_models, data_name = 'probe knn intersection unbalanced')
# 
# ###### probe lsa union 
# probe_lsa_union_unbal_models <- runModels(probe_lsa, random = F, bump_hunter = T, 
#                                     bump_hunter_data = bh_union_unbal)
# # get result table 
# probe_lsa_union_unbal_table <- extractResults(probe_lsa_union_unbal_models, data_name = 'probe lsa union unbalanced')
# 
# ###### probe lsa intersection
# probe_lsa_int_unbal_models <- runModels(probe_lsa, random = F, bump_hunter = T, 
#                                   bump_hunter_data = bh_intersection_unbal)
# # get result table 
# probe_lsa_int_unbal_table <- extractResults(probe_lsa_int_unbal_models, data_name = 'probe lsa intersection unbalanced')
# 
# 
# ###################################
# # Finally run each gene and probe data with random features
# ###################################
# 
# ###### GENE KNN
# gene_knn_rand_models <- runModels(gene_knn, random = T)
# 
# # get table
# gene_knn_rand_table <- extractResults(gene_knn_rand_models, data_name = 'gene knn random')
# 
# ###### GENE lsa
# gene_lsa_rand_models <- runModels(gene_lsa, random = T)
# 
# # get table
# gene_lsa_rand_table <- extractResults(gene_lsa_rand_models, data_name = 'gene lsa random')
# 
# ###### probe KNN
# probe_knn_rand_models <- runModels(probe_knn, random = T)
# 
# # get table
# probe_knn_rand_table <- extractResults(probe_knn_rand_models, data_name = 'probe knn random')
# 
# ###### probe lsa
# probe_lsa_rand_models <- runModels(probe_lsa, random = T)
# 
# # get table
# probe_lsa_rand_table <- extractResults(probe_lsa_rand_models, data_name = 'probe lsa random')
# 
# ###################################
# # Aggregate results tables 
# ###################################
# final_results <- rbind(probe_knn_global_table, 
#                        probe_knn_cancer_table,
#                        probe_lsa_global_table,
#                        probe_lsa_cancer_table,
#                        probe_knn_union_table, 
#                        probe_knn_int_table,
#                        probe_lsa_union_table,
#                        probe_lsa_int_table, 
#                        probe_knn_global_unbal_table,
#                        probe_knn_cancer_unbal_table,
#                        probe_lsa_global_unbal_table,
#                        probe_lsa_cancer_unbal_table, 
#                        probe_knn_union_unbal_table,
#                        probe_knn_int_unbal_table,
#                        probe_lsa_union_unbal_table,
#                        probe_lsa_int_unbal_table,
#                        gene_knn_rand_table, 
#                        gene_lsa_rand_table, 
#                        probe_knn_rand_table,
#                        probe_lsa_rand_table)
# 
# 
# # # Save main model data (use of all features)
# # save.image(paste0(model_data, '/partial_bh_results.RData'))
# # load(paste0(model_data, '/partial_bh_results.RData'))
# # 
# # # get plot objects
# # plot_mut <- plotObject(gene_knn_models, residual = F, p53_mut = T)
# # plot_all <- plotObject(gene_knn_models, residual = F, p53_mut = F)
# # plot_resid_mut <- plotObject(gene_knn_models, residual = T, p53_mut = T)
# # plot_resid_all <- plotObject(gene_knn_models, residual = T, p53_mut = F)
# # 
# # # get confusion matrix objects
# # con_48_mut <- matObject(gene_knn_models, age = 48, residual = F, p53_mut = T)
# # con_48_all <- matObject(gene_knn_models, age = 48, residual = F, p53_mut = F)
# # con_60_mut <- matObject(gene_knn_models, age = 60, residual = F, p53_mut = T)
# # con_60_all <- matObject(gene_knn_models, age = 60, residual = F, p53_mut = F)
# # con_72_mut <- matObject(gene_knn_models, age = 72, residual = F, p53_mut = T)
# # con_72_all <- matObject(gene_knn_models, age = 72, residual = F, p53_mut = F)
# # con_84_mut <- matObject(gene_knn_models, age = 84, residual = F, p53_mut = T)
# # con_84_all <- matObject(gene_knn_models, age = 84, residual = F, p53_mut = F)
# # 
# # con_48_mut <- matObject(gene_knn_models, age = 48, residual = F, p53_mut = T)
# # con_48_all <- matObject(gene_knn_models, age = 48, residual = F, p53_mut = F)
# # con_60_mut <- matObject(gene_knn_models, age = 60, residual = F, p53_mut = T)
# # con_60_all <- matObject(gene_knn_models, age = 60, residual = F, p53_mut = F)
# # con_72_mut <- matObject(gene_knn_models, age = 72, residual = F, p53_mut = T)
# # con_72_all <- matObject(gene_knn_models, age = 72, residual = F, p53_mut = F)
# # con_84_mut <- matObject(gene_knn_models, age = 84, residual = F, p53_mut = T)
# # con_84_all <- matObject(gene_knn_models, age = 84, residual = F, p53_mut = F)
# 
# 
# 
# 
# 
# 
# 
# 
