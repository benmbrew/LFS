####### This script will predict age of onset and save results
# this is 7th step in pipeline

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

registerDoParallel(1)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
scripts_folder <- paste0(project_folder, '/Scripts')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
idat_data <- paste0(methyl_data, '/raw_files')
clin_data <- paste0(data_folder, '/clin_data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')
model_data <- paste0(data_folder, '/model_data')
results_folder <- paste0(project_folder, '/Results')

##########
# source model_functions to get functions to run models 
##########
source(paste0(scripts_folder, '/predict_age/model_functions.R'))

##########
# set parameters 
##########
data_thresholds <- c(48, 60, 72, 84)
p53 <- c('Mut', 'WT')

##########
# Read in cases and bumphunter features
##########
load(paste0(model_data, '/model_data_cases.RData'))
load(paste0(model_data, '/bh_features.RData'))

##########
# Read in controls (if needed)
##########
# load(paste0(model_data, '/model_data_controls.RData'))

##########
# function to run models - subset, get residuals, get categorical, predict with regression and fac. 
##########
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

##########
# get stats of data used in model
##########
# get data info from 
data_stats <- dataStats(beta_raw)

saveRDS(data_stats, file = paste0(results_folder, '/model_data_stats.rda'))

# remove object
rm(data_stats)

##########
# beta_raw
##########
# cancer balanced 
beta_raw_bal_cancer_models <- runModels(beta_raw, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_raw_bal_cancer_features)

beta_raw_bal_cancer_table <- extractResults(beta_raw_bal_cancer_models, data_name = 'beta_raw_bal_cancer')


# cancer balanced counts
beta_raw_bal_counts_cancer_models <- runModels(beta_raw, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_raw_bal_counts_cancer_features)

beta_raw_bal_counts_cancer_table <- extractResults(beta_raw_bal_counts_cancer_models, data_name = 'beta_raw_bal_counts_cancer')

# cancer unbalanced 
beta_raw_unbal_cancer_models <- runModels(beta_raw, random = F, bump_hunter = T, 
                                               bump_hunter_data = beta_raw_unbal_cancer_features)

beta_raw_unbal_cancer_table <- extractResults(beta_raw_unbal_cancer_models, data_name = 'beta_raw_unbal_cancer')

# p53 balanced 
beta_raw_bal_p53_models <- runModels(beta_raw, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_raw_bal_p53_features)

beta_raw_bal_p53_table <- extractResults(beta_raw_bal_p53_models, data_name = 'beta_raw_bal_p53')


# p53 balanced counts
beta_raw_bal_counts_p53_models <- runModels(beta_raw, random = F, bump_hunter = T, 
                                               bump_hunter_data = beta_raw_bal_counts_p53_features)

beta_raw_bal_counts_p53_table <- extractResults(beta_raw_bal_counts_p53_models, data_name = 'beta_raw_bal_counts_p53')

# p53 unbalanced 
beta_raw_unbal_p53_models <- runModels(beta_raw, random = F, bump_hunter = T, 
                                          bump_hunter_data = beta_raw_unbal_p53_features)

beta_raw_unbal_p53_table <- extractResults(beta_raw_unbal_p53_models, data_name = 'beta_raw_unbal_p53')

# load('/home/benbrew/Desktop/temp_raw_model_results.RData')

##########
# rbind results
##########
beta_raw_results <- as.data.frame(rbind(beta_raw_bal_cancer_table, beta_raw_bal_counts_cancer_table, beta_raw_unbal_cancer_table,
                                    beta_raw_bal_p53_table, beta_raw_bal_counts_p53_table, beta_raw_unbal_p53_table))


# get rows and columns variable 
beta_raw_results <- getDims(beta_raw_results)

##########
# save results table and models for beta raw
##########
# first save results table or beta_raw
saveRDS(beta_raw_results, file = paste0(results_folder, '/beta_raw_model_results.rda'))

# save models for beta_raw cancer
saveRDS(beta_raw_unbal_cancer_models, file = paste0(results_folder, '/beta_raw_unbal_cancer_models.rda'))
saveRDS(beta_raw_bal_cancer_models, file = paste0(results_folder, '/beta_raw_bal_cancer_models.rda'))
saveRDS(beta_raw_bal_counts_cancer_models, file = paste0(results_folder, '/beta_raw_bal_counts_cancer_models.rda'))

# save models for beta_raw p53
saveRDS(beta_raw_unbal_p53_models, file = paste0(results_folder, '/beta_raw_unbal_p53_models.rda'))
saveRDS(beta_raw_bal_p53_models, file = paste0(results_folder, '/beta_raw_bal_p53_models.rda'))
saveRDS(beta_raw_bal_counts_p53_models, file = paste0(results_folder, '/beta_raw_bal_counts_p53_models.rda'))

##########
# remove unneeded objects
##########
rm(beta_raw_bal_cancer_features, beta_raw_bal_counts_cancer_features, beta_raw_unbal_cancer_features,
   beta_raw_bal_p53_features, beta_raw_bal_counts_p53_features, beta_raw_unbal_p53_features,
   beta_raw_bal_cancer_models, beta_raw_bal_counts_cancer_models, beta_raw_unbal_cancer_models,
   beta_raw_bal_p53_models, beta_raw_bal_counts_p53_models, beta_raw_unbal_p53_models,
   beta_raw_bal_cancer_table, beta_raw_bal_counts_cancer_table, beta_raw_unbal_cancer_table,
   beta_raw_bal_p53_table, beta_raw_bal_counts_p53_table, beta_raw_unbal_p53_table, beta_raw_results)

######################################################################################################################

##########
# beta_swan
##########
# cancer balanced 
beta_swan_bal_cancer_models <- runModels(beta_swan, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_swan_bal_cancer_features)

beta_swan_bal_cancer_table <- extractResults(beta_swan_bal_cancer_models, data_name = 'beta_swan_bal_cancer')


# cancer balanced counts
beta_swan_bal_counts_cancer_models <- runModels(beta_swan, random = F, bump_hunter = T, 
                                               bump_hunter_data = beta_swan_bal_counts_cancer_features)

beta_swan_bal_counts_cancer_table <- extractResults(beta_swan_bal_counts_cancer_models, data_name = 'beta_swan_bal_counts_cancer')

# cancer unbalanced 
beta_swan_unbal_cancer_models <- runModels(beta_swan, random = F, bump_hunter = T, 
                                          bump_hunter_data = beta_swan_unbal_cancer_features)

beta_swan_unbal_cancer_table <- extractResults(beta_swan_unbal_cancer_models, data_name = 'beta_swan_unbal_cancer')

# p53 balanced 
beta_swan_bal_p53_models <- runModels(beta_swan, random = F, bump_hunter = T, 
                                     bump_hunter_data = beta_swan_bal_p53_features)

beta_swan_bal_p53_table <- extractResults(beta_swan_bal_p53_models, data_name = 'beta_swan_bal_p53')


# p53 balanced counts
beta_swan_bal_counts_p53_models <- runModels(beta_swan, random = F, bump_hunter = T, 
                                            bump_hunter_data = beta_swan_bal_counts_p53_features)

beta_swan_bal_counts_p53_table <- extractResults(beta_swan_bal_counts_p53_models, data_name = 'beta_swan_bal_counts_p53')

# p53 unbalanced 
beta_swan_unbal_p53_models <- runModels(beta_swan, random = F, bump_hunter = T, 
                                       bump_hunter_data = beta_swan_unbal_p53_features)

beta_swan_unbal_p53_table <- extractResults(beta_swan_unbal_p53_models, data_name = 'beta_swan_unbal_p53')

##########
# rbind results
##########
beta_swan_results <- as.data.frame(rbind(beta_swan_bal_cancer_table, beta_swan_bal_counts_cancer_table, beta_swan_unbal_cancer_table,
                                        beta_swan_bal_p53_table, beta_swan_bal_counts_p53_table, beta_swan_unbal_p53_table))


# get rows and columns variable 
beta_swan_results <- getDims(beta_swan_results)

##########
# save results table and models for beta swan
##########
# first save results table or beta_swan
saveRDS(beta_swan_results, file = paste0(results_folder, '/beta_swan_model_results.rda'))

# save models for beta_swan cancer
saveRDS(beta_swan_unbal_cancer_models, file = paste0(results_folder, '/beta_swan_unbal_cancer_models.rda'))
saveRDS(beta_swan_bal_cancer_models, file = paste0(results_folder, '/beta_swan_bal_cancer_models.rda'))
saveRDS(beta_swan_bal_counts_cancer_models, file = paste0(results_folder, '/beta_swan_bal_counts_cancer_models.rda'))

# save models for beta_swan p53
saveRDS(beta_swan_unbal_p53_models, file = paste0(results_folder, '/beta_swan_unbal_p53_models.rda'))
saveRDS(beta_swan_bal_p53_models, file = paste0(results_folder, '/beta_swan_bal_p53_models.rda'))
saveRDS(beta_swan_bal_counts_p53_models, file = paste0(results_folder, '/beta_swan_bal_counts_p53_models.rda'))

##########
# remove unneeded objects
##########
rm(beta_swan_bal_cancer_features, beta_swan_bal_counts_cancer_features, beta_swan_unbal_cancer_features,
   beta_swan_bal_p53_features, beta_swan_bal_counts_p53_features, beta_swan_unbal_p53_features,
   beta_swan_bal_cancer_models, beta_swan_bal_counts_cancer_models, beta_swan_unbal_cancer_models,
   beta_swan_bal_p53_models, beta_swan_bal_counts_p53_models, beta_swan_unbal_p53_models,
   beta_swan_bal_cancer_table, beta_swan_bal_counts_cancer_table, beta_swan_unbal_cancer_table,
   beta_swan_bal_p53_table, beta_swan_bal_counts_p53_table, beta_swan_unbal_p53_table, beta_swan_results)


######################################################################################################################

##########
# beta_quan
##########
# cancer balanced 
beta_quan_bal_cancer_models <- runModels(beta_quan, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_quan_bal_cancer_features)

beta_quan_bal_cancer_table <- extractResults(beta_quan_bal_cancer_models, data_name = 'beta_quan_bal_cancer')


# cancer balanced counts
beta_quan_bal_counts_cancer_models <- runModels(beta_quan, random = F, bump_hunter = T, 
                                               bump_hunter_data = beta_quan_bal_counts_cancer_features)

beta_quan_bal_counts_cancer_table <- extractResults(beta_quan_bal_counts_cancer_models, data_name = 'beta_quan_bal_counts_cancer')

# cancer unbalanced 
beta_quan_unbal_cancer_models <- runModels(beta_quan, random = F, bump_hunter = T, 
                                          bump_hunter_data = beta_quan_unbal_cancer_features)

beta_quan_unbal_cancer_table <- extractResults(beta_quan_unbal_cancer_models, data_name = 'beta_quan_unbal_cancer')

# p53 balanced 
beta_quan_bal_p53_models <- runModels(beta_quan, random = F, bump_hunter = T, 
                                     bump_hunter_data = beta_quan_bal_p53_features)

beta_quan_bal_p53_table <- extractResults(beta_quan_bal_p53_models, data_name = 'beta_quan_bal_p53')


# p53 balanced counts
beta_quan_bal_counts_p53_models <- runModels(beta_quan, random = F, bump_hunter = T, 
                                            bump_hunter_data = beta_quan_bal_counts_p53_features)

beta_quan_bal_counts_p53_table <- extractResults(beta_quan_bal_counts_p53_models, data_name = 'beta_quan_bal_counts_p53')

# p53 unbalanced 
beta_quan_unbal_p53_models <- runModels(beta_quan, random = F, bump_hunter = T, 
                                       bump_hunter_data = beta_quan_unbal_p53_features)

beta_quan_unbal_p53_table <- extractResults(beta_quan_unbal_p53_models, data_name = 'beta_quan_unbal_p53')

##########
# rbind results
##########
beta_quan_results <- as.data.frame(rbind(beta_quan_bal_cancer_table, beta_quan_bal_counts_cancer_table, beta_quan_unbal_cancer_table,
                                        beta_quan_bal_p53_table, beta_quan_bal_counts_p53_table, beta_quan_unbal_p53_table))

# get rows and columns variable 
beta_quan_results <- getDims(beta_quan_results)

##########
# save results table and models for beta quan
##########
# first save results table or beta_quan
saveRDS(beta_quan_results, file = paste0(results_folder, '/beta_quan_model_results.rda'))

# save models for beta_quan cancer
saveRDS(beta_quan_unbal_cancer_models, file = paste0(results_folder, '/beta_quan_unbal_cancer_models.rda'))
saveRDS(beta_quan_bal_cancer_models, file = paste0(results_folder, '/beta_quan_bal_cancer_models.rda'))
saveRDS(beta_quan_bal_counts_cancer_models, file = paste0(results_folder, '/beta_quan_bal_counts_cancer_models.rda'))

# save models for beta_quan p53
saveRDS(beta_quan_unbal_p53_models, file = paste0(results_folder, '/beta_quan_unbal_p53_models.rda'))
saveRDS(beta_quan_bal_p53_models, file = paste0(results_folder, '/beta_quan_bal_p53_models.rda'))
saveRDS(beta_quan_bal_counts_p53_models, file = paste0(results_folder, '/beta_quan_bal_counts_p53_models.rda'))

##########
# remove unneeded objects
##########
rm(beta_quan_bal_cancer_features, beta_quan_bal_counts_cancer_features, beta_quan_unbal_cancer_features,
   beta_quan_bal_p53_features, beta_quan_bal_counts_p53_features, beta_quan_unbal_p53_features,
   beta_quan_bal_cancer_models, beta_quan_bal_counts_cancer_models, beta_quan_unbal_cancer_models,
   beta_quan_bal_p53_models, beta_quan_bal_counts_p53_models, beta_quan_unbal_p53_models,
   beta_quan_bal_cancer_table, beta_quan_bal_counts_cancer_table, beta_quan_unbal_cancer_table,
   beta_quan_bal_p53_table, beta_quan_bal_counts_p53_table, beta_quan_unbal_p53_table, beta_quan_results)

#########################################################################################################################

##########
# beta_funnorm
##########
# cancer balanced 
beta_funnorm_bal_cancer_models <- runModels(beta_funnorm, random = F, bump_hunter = T, 
                                        bump_hunter_data = beta_funnorm_bal_cancer_features)

beta_funnorm_bal_cancer_table <- extractResults(beta_funnorm_bal_cancer_models, data_name = 'beta_funnorm_bal_cancer')


# cancer balanced counts
beta_funnorm_bal_counts_cancer_models <- runModels(beta_funnorm, random = F, bump_hunter = T, 
                                               bump_hunter_data = beta_funnorm_bal_counts_cancer_features)

beta_funnorm_bal_counts_cancer_table <- extractResults(beta_funnorm_bal_counts_cancer_models, data_name = 'beta_funnorm_bal_counts_cancer')

# cancer unbalanced 
beta_funnorm_unbal_cancer_models <- runModels(beta_funnorm, random = F, bump_hunter = T, 
                                          bump_hunter_data = beta_funnorm_unbal_cancer_features)

beta_funnorm_unbal_cancer_table <- extractResults(beta_funnorm_unbal_cancer_models, data_name = 'beta_funnorm_unbal_cancer')

# p53 balanced 
beta_funnorm_bal_p53_models <- runModels(beta_funnorm, random = F, bump_hunter = T, 
                                     bump_hunter_data = beta_funnorm_bal_p53_features)

beta_funnorm_bal_p53_table <- extractResults(beta_funnorm_bal_p53_models, data_name = 'beta_funnorm_bal_p53')


# p53 balanced counts
beta_funnorm_bal_counts_p53_models <- runModels(beta_funnorm, random = F, bump_hunter = T, 
                                            bump_hunter_data = beta_funnorm_bal_counts_p53_features)

beta_funnorm_bal_counts_p53_table <- extractResults(beta_funnorm_bal_counts_p53_models, data_name = 'beta_funnorm_bal_counts_p53')

# p53 unbalanced 
beta_funnorm_unbal_p53_models <- runModels(beta_funnorm, random = F, bump_hunter = T, 
                                       bump_hunter_data = beta_funnorm_unbal_p53_features)

beta_funnorm_unbal_p53_table <- extractResults(beta_funnorm_unbal_p53_models, data_name = 'beta_funnorm_unbal_p53')

##########
# rbind results
##########
beta_funnorm_results <- as.data.frame(rbind(beta_funnorm_bal_cancer_table, beta_funnorm_bal_counts_cancer_table, beta_funnorm_unbal_cancer_table,
                                        beta_funnorm_bal_p53_table, beta_funnorm_bal_counts_p53_table, beta_funnorm_unbal_p53_table))


# get rows and columns variable 
beta_funnorm_results <- getDims(beta_funnorm_results)

##########
# save results table and models for beta funnorm
##########
# first save results table or beta_funnorm
saveRDS(beta_funnorm_results, file = paste0(results_folder, '/beta_funnorm_model_results.rda'))

# save models for beta_funnorm cancer
saveRDS(beta_funnorm_unbal_cancer_models, file = paste0(results_folder, '/beta_funnorm_unbal_cancer_models.rda'))
saveRDS(beta_funnorm_bal_cancer_models, file = paste0(results_folder, '/beta_funnorm_bal_cancer_models.rda'))
saveRDS(beta_funnorm_bal_counts_cancer_models, file = paste0(results_folder, '/beta_funnorm_bal_counts_cancer_models.rda'))

# save models for beta_funnorm p53
saveRDS(beta_funnorm_unbal_p53_models, file = paste0(results_folder, '/beta_funnorm_unbal_p53_models.rda'))
saveRDS(beta_funnorm_bal_p53_models, file = paste0(results_folder, '/beta_funnorm_bal_p53_models.rda'))
saveRDS(beta_funnorm_bal_counts_p53_models, file = paste0(results_folder, '/beta_funnorm_bal_counts_p53_models.rda'))

##########
# remove unneeded objects
##########
rm(beta_funnorm_bal_cancer_features, beta_funnorm_bal_counts_cancer_features, beta_funnorm_unbal_cancer_features,
   beta_funnorm_bal_p53_features, beta_funnorm_bal_counts_p53_features, beta_funnorm_unbal_p53_features,
   beta_funnorm_bal_cancer_models, beta_funnorm_bal_counts_cancer_models, beta_funnorm_unbal_cancer_models,
   beta_funnorm_bal_p53_models, beta_funnorm_bal_counts_p53_models, beta_funnorm_unbal_p53_models,
   beta_funnorm_bal_cancer_table, beta_funnorm_bal_counts_cancer_table, beta_funnorm_unbal_cancer_table,
   beta_funnorm_bal_p53_table, beta_funnorm_bal_counts_p53_table, beta_funnorm_unbal_p53_table, beta_funnorm_results)

######################################################################################################################

##########
# Now run each beta with random features
##########
beta_raw_rand_models <- runModels(beta_raw, random = T, bump_hunter = F)
beta_raw_rand_table <- extractResults(beta_raw_rand_models, data_name = 'beta_raw_rand')

beta_swan_rand_models <- runModels(beta_swan, random = T, bump_hunter = F)
beta_swan_rand_table <- extractResults(beta_swan_rand_models, data_name = 'beta_swan_rand')

beta_quan_rand_models <- runModels(beta_quan, random = T, bump_hunter = F)
beta_quan_rand_table <- extractResults(beta_quan_rand_models, data_name = 'beta_quan_rand')

beta_funnorm_rand_models <- runModels(beta_quan, random = T, bump_hunter = F)
beta_funnorm_rand_table <- extractResults(beta_funnorm_rand_models, data_name = 'beta_funnorm_rand')

##########
# rbind results
##########
beta_rand_results <- as.data.frame(rbind(beta_raw_rand_table, 
                                         beta_swan_rand_table, 
                                         beta_quan_rand_table,
                                         beta_funnorm_rand_table))


# get rows and columns variable 
beta_rand_results <- getDims(beta_rand_results)

##########
# save results table and models for beta raw
##########
# first save results table or beta_raw
saveRDS(beta_rand_results, file = paste0(results_folder, '/beta_rand_model_results.rda'))

# save models for beta_raw cancer
saveRDS(beta_raw_rand_models, file = paste0(results_folder, '/beta_raw_rand_models.rda'))
saveRDS(beta_swan_rand_models, file = paste0(results_folder, '/beta_swan_rand_models.rda'))
saveRDS(beta_quan_rand_models, file = paste0(results_folder, '/beta_quan_rand_models.rda'))
saveRDS(beta_funnorm_rand_models, file = paste0(results_folder, '/beta_funnorm_rand_models.rda'))

##########
# remove unneeded objects
##########
rm(beta_raw_rand_models, beta_swan_rand_models, beta_quan_rand_models, beta_funnorm_rand_models,
   beta_raw_rand_table, beta_swan_rand_table, beta_quan_rand_table, beta_funnorm_rand_table, beta_rand_results)
