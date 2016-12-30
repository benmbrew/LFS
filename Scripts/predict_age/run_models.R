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
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')
results_folder <- paste0(project_folder, '/Results')
scripts_folder <- paste0(project_folder, '/Scripts')


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
# cases - raw, swan, quan, funnorm, and clin
load(paste0(model_data, '/model_data_cases.RData'))
# bh_features
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
                      bump_hunter_data,
                      num_features = NULL) 
{
  
  # get differenct variations of data
  data <- subsetDat(data)
  
  if (bump_hunter) {
    
    data <- bhSubset(data, bh_data = bump_hunter_data)

  }
  
  if (random) {
    
  data <- getRand(data, num_features)
    
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
  
  for (dat in 1:length(data_fac)) {
    
    sub_dat_fac <- data_fac[[dat]]
    sub_dat_resid_fac <- data_resid_fac[[dat]]
    
    temp.data_fac_result <- list()
    temp.data_resid_fac_result <- list()
    
    for (sub_dat in 1:length(sub_dat_fac)) {
      
      temp.sub_dat_fac <- sub_dat_fac[[sub_dat]]
      temp.sub_dat_resid_fac <- sub_dat_resid_fac[[sub_dat]]
      
      temp.data_fac_result[[sub_dat]] <- rfPredictFac(temp.sub_dat_fac, cutoff = .7, iterations = 10)
      temp.data_resid_fac_result[[sub_dat]] <- rfPredictFac(temp.sub_dat_resid_fac, cutoff = .7, iterations = 10)
      
    }
    
    data_fac_result[[dat]] <- temp.data_fac_result
    data_resid_fac_result[[dat]] <- temp.data_resid_fac_result
    
  }
  
  return (list(data_result, 
              data_resid_result, 
              data_fac_result, 
              data_resid_fac_result))
  
}

##########
# get stats of data used in model
##########
# get data info from 
data_stats <- dataStats(beta_raw)

saveRDS(data_stats, 
        file = paste0(results_folder, 
                                  '/model_data_stats.rda'))

# remove object
rm(data_stats)

###################################################################################################################################
## beta_raw

##########
# cancer
##########

# bal cancer
raw_bal_cancer_models <- runModels(beta_raw, 
                                   random = F, 
                                   bump_hunter = T, 
                                   bump_hunter_data = beta_raw_bal_cancer_features)

raw_bal_cancer_table <- extractResults(raw_bal_cancer_models, 
                                       data_name = 'beta_raw_bal_cancer')

# bal cancer sig
raw_bal_cancer_sig_models <- runModels(beta_raw, 
                                       random = F, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_raw_bal_cancer_sig_features)

raw_bal_cancer_sig_table <- extractResults(raw_bal_cancer_sig_models, 
                                           data_name = 'raw_bal_cancer_sig')


# bal counts cancer
raw_bal_counts_cancer_models <- runModels(beta_raw, 
                                          random = F, 
                                          bump_hunter = T, 
                                          bump_hunter_data = beta_raw_bal_counts_cancer_features)

raw_bal_counts_cancer_table <- extractResults(raw_bal_counts_cancer_models, 
                                              data_name = 'raw_bal_counts_cancer')

# bal counts cancer sig
raw_bal_counts_cancer_sig_models <- runModels(beta_raw, 
                                              random = F, 
                                              bump_hunter = T, 
                                              bump_hunter_data =  beta_raw_bal_counts_cancer_sig_features)

raw_bal_counts_cancer_sig_table <- extractResults(raw_bal_counts_cancer_sig_models, 
                                                  data_name = 'raw_bal_counts_cancer_sig')

# unbal cancer
raw_unbal_cancer_models <- runModels(beta_raw, 
                                     random = F, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_raw_unbal_cancer_features)

raw_unbal_cancer_table <- extractResults(raw_unbal_cancer_models, 
                                         data_name = 'raw_unbal_cancer')


# unbal cancer sig
raw_unbal_cancer_sig_models <- runModels(beta_raw, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_raw_unbal_cancer_sig_features)

raw_unbal_cancer_sig_table <- extractResults(raw_unbal_cancer_sig_models, 
                                             data_name = 'raw_unbal_cancer_sig')

##########
# p53
##########

# bal p53
raw_bal_p53_models <- runModels(beta_raw, 
                                random = F, 
                                bump_hunter = T, 
                                bump_hunter_data = beta_raw_bal_p53_features)

raw_bal_p53_table <- extractResults(raw_bal_p53_models, 
                                    data_name = 'raw_bal_p53')

# bal p53 sig
raw_bal_p53_sig_models <- runModels(beta_raw, 
                                    random = F, 
                                    bump_hunter = T, 
                                    bump_hunter_data = beta_raw_bal_p53_sig_features)

raw_bal_p53_sig_table <- extractResults(raw_bal_p53_sig_models, 
                                        data_name = 'raw_bal_p53_sig')

# bal counts p53
raw_bal_counts_p53_models <- runModels(beta_raw, 
                                       random = F, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_raw_bal_counts_p53_features)

raw_bal_counts_p53_table <- extractResults(raw_bal_counts_p53_models, 
                                           data_name = 'raw_bal_counts_p53')

# bal counts p53 sig
raw_bal_counts_p53_sig_models <- runModels(beta_raw, 
                                           random = F, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_raw_bal_counts_p53_sig_features)

raw_bal_counts_p53_sig_table <- extractResults(raw_bal_counts_p53_sig_models, 
                                               data_name = 'raw_bal_counts_p53_sig')

# unbal p53
raw_unbal_p53_models <- runModels(beta_raw, 
                                     random = F, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_raw_unbal_p53_features)

raw_unbal_p53_table <- extractResults(raw_unbal_p53_models, 
                                         data_name = 'raw_unbal_p53')


# unbal p53 sig
raw_unbal_p53_sig_models <- runModels(beta_raw, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_raw_unbal_p53_sig_features)

raw_unbal_p53_sig_table <- extractResults(raw_unbal_p53_sig_models, 
                                             data_name = 'raw_unbal_p53_sig')


##########
# cancer
##########

# cancer_intersection
raw_cancer_intersection_models <- runModels(beta_raw, 
                                        random = F, 
                                        bump_hunter = T, 
                                        bump_hunter_data = beta_raw_cancer_intersection_features)

raw_cancer_intersection_table <- extractResults(raw_cancer_intersection_models, 
                                            data_name = 'raw_cancer_intersection')

# cancer_intersection sig
raw_cancer_intersection_sig_models <- runModels(beta_raw, 
                                            random = F, 
                                            bump_hunter = T, 
                                            bump_hunter_data = beta_raw_cancer_intersection_sig_features)

raw_cancer_intersection_sig_table <- extractResults(raw_cancer_intersection_sig_models, 
                                            data_name = 'raw_cancer_intersection_sig')

# cancer_union
raw_cancer_union_models <- runModels(beta_raw, 
                                 random = F, 
                                 bump_hunter = T, 
                                 bump_hunter_data = beta_raw_cancer_union_features)

raw_cancer_union_table <- extractResults(raw_cancer_union_models, 
                                     data_name = 'raw_cancer_union')

# cancer_union sig
raw_cancer_union_sig_models <- runModels(beta_raw, 
                                     random = F, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_raw_cancer_union_sig_features)

raw_cancer_union_sig_table <- extractResults(raw_cancer_union_sig_models, 
                                         data_name = 'raw_cancer_union_sig')
##########
# p53
##########
# p53_intersection
raw_p53_intersection_models <- runModels(beta_raw, 
                                     random = F, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_raw_p53_intersection_features)

raw_p53_intersection_table <- extractResults(raw_p53_intersection_models, 
                                            data_name = 'raw_p53_intersection')

# p53_intersection sig
raw_p53_intersection_sig_models <- runModels(beta_raw, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_raw_p53_intersection_sig_features)

raw_p53_intersection_sig_table <- extractResults(raw_p53_intersection_sig_models, 
                                             data_name = 'raw_p53_intersection_sig')

# p53_union
raw_p53_union_models <- runModels(beta_raw, 
                              random = F, 
                              bump_hunter = T, 
                              bump_hunter_data = beta_raw_p53_union_features)

raw_p53_union_table <- extractResults(raw_p53_union_models, 
                                      data_name = 'raw_p53_union')

# p53_union sig
raw_p53_union_sig_models <- runModels(beta_raw, 
                                      random = F, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_raw_p53_union_sig_features)

raw_p53_union_sig_table <- extractResults(raw_p53_union_sig_models, 
                                      data_name = 'raw_p53_union_sig')


##########
# most balanced cancer 
##########

# intersection
raw_bal_counts_cancer_intersection_models <- runModels(beta_raw, 
                                                       random = F, 
                                                       bump_hunter = T, 
                                                       bump_hunter_data = beta_bal_counts_cancer_intersection_features)

raw_bal_counts_cancer_intersection_table <- extractResults(raw_bal_counts_cancer_intersection_models, 
                                                           data_name = 'raw_bal_counts_cancer_intersection')


raw_bal_counts_cancer_intersection_sig_models <- runModels(beta_raw, 
                                                           random = F, 
                                                           bump_hunter = T, 
                                                           bump_hunter_data = beta_bal_counts_cancer_intersection_sig_features)

raw_bal_counts_cancer_intersection_sig_table <- extractResults(raw_bal_counts_cancer_intersection_sig_models, 
                                                               data_name = 'raw_bal_counts_cancer_intersection_sig')

# union
raw_bal_counts_cancer_union_models <- runModels(beta_raw, 
                                                random = F, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_bal_counts_cancer_union_features)

raw_bal_counts_cancer_union_table <- extractResults(raw_bal_counts_cancer_union_models, 
                                                    data_name = 'raw_bal_counts_cancer_union')


raw_bal_counts_cancer_union_sig_models <- runModels(beta_raw, 
                                                    random = F, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_cancer_union_sig_features)

raw_bal_counts_cancer_union_sig_table <- extractResults(raw_bal_counts_cancer_union_sig_models, 
                                                        data_name = 'raw_bal_counts_cancer_union_sig')

##########
# most balanced p53 
##########

# intersection
raw_bal_counts_p53_intersection_models <- runModels(beta_raw, 
                                                    random = F, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_p53_intersection_features)

raw_bal_counts_p53_intersection_table <- extractResults(raw_bal_counts_p53_intersection_models, 
                                                        data_name = 'raw_bal_counts_p53_intersection')


raw_bal_counts_p53_intersection_sig_models <- runModels(beta_raw, 
                                                        random = F, 
                                                        bump_hunter = T, 
                                                        bump_hunter_data = beta_bal_counts_p53_intersection_sig_features)

raw_bal_counts_p53_intersection_sig_table <- extractResults(raw_bal_counts_p53_intersection_sig_models, 
                                                            data_name = 'raw_bal_counts_p53_intersection_sig')


# union
raw_bal_counts_p53_union_models <- runModels(beta_raw, 
                                             random = F, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_bal_counts_p53_union_features)

raw_bal_counts_p53_union_table <- extractResults(raw_bal_counts_p53_union_models, 
                                                 data_name = 'raw_bal_counts_p53_union')


raw_bal_counts_p53_union_sig_models <- runModels(beta_raw, 
                                                 random = F, 
                                                 bump_hunter = T, 
                                                 bump_hunter_data = beta_bal_counts_p53_union_sig_features)

raw_bal_counts_p53_union_sig_table <- extractResults(raw_bal_counts_p53_union_sig_models, 
                                                     data_name = 'raw_bal_counts_p53_union_sig')

###########
# complete cancer - intersection across each method intersection
###########

# complete cancer intersection
raw_complete_cancer_intersection_models <- runModels(beta_raw, 
                                                     random = F, 
                                                     bump_hunter = T, 
                                                     bump_hunter_data = beta_cancer_intersection_features)

# get table 
raw_complete_cancer_intersection_table <- extractResults(raw_complete_cancer_intersection_models, 
                                                         data_name = 'raw_complete_cancer_intersection')

# complete cancer intersection sig
raw_complete_cancer_intersection_sig_models <- runModels(beta_raw, 
                                                         random = F, 
                                                         bump_hunter = T, 
                                                         bump_hunter_data = beta_cancer_intersection_sig_features)

# get table 
raw_complete_cancer_intersection_sig_table <- extractResults(raw_complete_cancer_intersection_sig_models, 
                                                             data_name = 'raw_complete_cancer_intersection_sig')


# complete cancer union
raw_complete_cancer_union_models <- runModels(beta_raw, 
                                                     random = F, 
                                                     bump_hunter = T, 
                                                     bump_hunter_data = beta_cancer_union_features)

# get table 
raw_complete_cancer_union_table <- extractResults(raw_complete_cancer_union_models, 
                                                         data_name = 'raw_complete_cancer_union')

# complete cancer union sig
raw_complete_cancer_union_sig_models <- runModels(beta_raw, 
                                                         random = F, 
                                                         bump_hunter = T, 
                                                         bump_hunter_data = beta_cancer_union_sig_features)

# get table 
raw_complete_cancer_union_sig_table <- extractResults(raw_complete_cancer_union_sig_models, 
                                                             data_name = 'raw_complete_cancer_union_sig')


###########
# complete p53 - intersection across each method intersection
###########

# complete p53 intersection
raw_complete_p53_intersection_models <- runModels(beta_raw, 
                                                     random = F, 
                                                     bump_hunter = T, 
                                                     bump_hunter_data = beta_p53_intersection_features)

# get table 
raw_complete_p53_intersection_table <- extractResults(raw_complete_p53_intersection_models, 
                                                         data_name = 'raw_complete_p53_intersection')

# complete p53 intersection sig
raw_complete_p53_intersection_sig_models <- runModels(beta_raw, 
                                                         random = F, 
                                                         bump_hunter = T, 
                                                         bump_hunter_data = beta_p53_intersection_sig_features)

# get table 
raw_complete_p53_intersection_sig_table <- extractResults(raw_complete_p53_intersection_sig_models, 
                                                             data_name = 'raw_complete_p53_intersection_sig')


# complete p53 union
raw_complete_p53_union_models <- runModels(beta_raw, 
                                              random = F, 
                                              bump_hunter = T, 
                                              bump_hunter_data = beta_p53_union_features)

# get table 
raw_complete_p53_union_table <- extractResults(raw_complete_p53_union_models, 
                                                  data_name = 'raw_complete_p53_union')

# complete p53 union sig
raw_complete_p53_union_sig_models <- runModels(beta_raw, 
                                                  random = F, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_p53_union_sig_features)

# get table 
raw_complete_p53_union_sig_table <- extractResults(raw_complete_p53_union_sig_models, 
                                                      data_name = 'raw_complete_p53_union_sig')


###########
# rbind tables and save RDA file
###########
raw_table <- rbind(raw_bal_cancer_table, raw_bal_cancer_sig_table, raw_bal_counts_cancer_table, raw_bal_counts_cancer_sig_table,
                   raw_unbal_cancer_table, raw_unbal_cancer_sig_table, raw_bal_p53_table, raw_bal_p53_sig_table, raw_bal_counts_p53_table, 
                   raw_bal_counts_p53_sig_table, raw_unbal_p53_table, raw_unbal_p53_sig_table, raw_cancer_intersection_table, 
                   raw_cancer_intersection_sig_table, raw_cancer_union_table, raw_cancer_union_sig_table, raw_p53_intersection_table, 
                   raw_p53_intersection_sig_table, raw_p53_union_table, raw_p53_union_sig_table,
                   raw_bal_counts_cancer_intersection_table, raw_bal_counts_cancer_intersection_sig_table,
                   raw_bal_counts_cancer_union_table, raw_bal_counts_cancer_union_sig_table,
                   raw_bal_counts_p53_intersection_table, raw_bal_counts_p53_intersection_sig_table,
                   raw_bal_counts_p53_union_table, raw_bal_counts_p53_union_sig_table,
                   raw_complete_cancer_intersection_table, raw_complete_cancer_intersection_sig_table,
                   raw_complete_cancer_union_table, raw_complete_cancer_union_sig_table,
                   raw_complete_p53_intersection_table, raw_complete_p53_intersection_sig_table,
                   raw_complete_p53_union_table, raw_complete_p53_union_sig_table)

# remove data 
rm(raw_bal_cancer_table, raw_bal_cancer_sig_table, raw_bal_counts_cancer_table, raw_bal_counts_cancer_sig_table,
   raw_unbal_cancer_table, raw_unbal_cancer_sig_table, raw_bal_p53_table, raw_bal_p53_sig_table, raw_bal_counts_p53_table, 
   raw_bal_counts_p53_sig_table, raw_unbal_p53_table, raw_unbal_p53_sig_table, raw_cancer_intersection_table, 
   raw_cancer_intersection_sig_table, raw_cancer_union_table, raw_cancer_union_sig_table, raw_p53_intersection_table, 
   raw_p53_intersection_sig_table, raw_p53_union_table, raw_p53_union_sig_table,
   raw_bal_counts_cancer_intersection_table, raw_bal_counts_cancer_intersection_sig_table,
   raw_bal_counts_cancer_union_table, raw_bal_counts_cancer_union_sig_table,
   raw_bal_counts_p53_intersection_table, raw_bal_counts_p53_intersection_sig_table,
   raw_bal_counts_p53_union_table, raw_bal_counts_p53_union_sig_table,
   raw_complete_cancer_intersection_table, raw_complete_cancer_intersection_sig_table,
   raw_complete_cancer_union_table, raw_complete_cancer_union_sig_table,
   raw_complete_p53_intersection_table, raw_complete_p53_intersection_sig_table,
   raw_complete_p53_union_table, raw_complete_p53_union_sig_table)


#save table 
saveRDS(raw_table, 
        file = paste0(results_folder, '/raw_table.rda'))
# run 8 new models
#LOAD raw table and assign to temp_new
#rbind with new models 
# save raw table
# remove tables
# save models
temp_new <- readRDS(file = paste0(results_folder, '/raw_table.rda'))
raw_table <- rbind(temp_new, raw_complete_cancer_intersection_table, raw_complete_cancer_intersection_sig_table,
                   raw_complete_cancer_union_table, raw_complete_cancer_union_sig_table,
                   raw_complete_p53_intersection_table, raw_complete_p53_intersection_sig_table,
                   raw_complete_p53_union_table, raw_complete_p53_union_sig_table)

rm(raw_complete_cancer_intersection_table, raw_complete_cancer_intersection_sig_table,
   raw_complete_cancer_union_table, raw_complete_cancer_union_sig_table,
   raw_complete_p53_intersection_table, raw_complete_p53_intersection_sig_table,
   raw_complete_p53_union_table, raw_complete_p53_union_sig_table)

###########
# save all model objects as RDA file
###########
saveRDS(raw_bal_cancer_models, 
        file = paste0(results_folder, '/raw_bal_cancer_models.rda'))
saveRDS(raw_bal_cancer_sig_models, 
        file = paste0(results_folder, '/raw_bal_cancer_sig_models.rda'))
saveRDS(raw_bal_counts_cancer_models, 
        file = paste0(results_folder, '/raw_bal_counts_cancer_models.rda'))
saveRDS(raw_bal_counts_cancer_sig_models, 
        file = paste0(results_folder, '/raw_bal_counts_cancer_sig_models.rda'))
saveRDS(raw_unbal_cancer_models, 
        file = paste0(results_folder, '/raw_unbal_cancer_models.rda'))
saveRDS(raw_unbal_cancer_sig_models, 
        file = paste0(results_folder, '/raw_unbal_cancer_sig_models.rda'))
saveRDS(raw_bal_p53_models, 
        file = paste0(results_folder, '/raw_bal_p53_models.rda'))
saveRDS(raw_bal_p53_sig_models, 
        file = paste0(results_folder, '/raw_bal_p53_sig_models.rda'))
saveRDS(raw_bal_counts_p53_models, 
        file = paste0(results_folder, '/raw_bal_counts_p53_models.rda'))
saveRDS(raw_bal_counts_p53_sig_models, 
        file = paste0(results_folder, '/raw_bal_counts_p53_sig_models.rda'))
saveRDS(raw_unbal_p53_models, 
        file = paste0(results_folder, '/raw_unbal_p53_models.rda'))
saveRDS(raw_unbal_p53_sig_models, 
        file = paste0(results_folder, '/raw_unbal_p53_sig_models.rda'))
saveRDS(raw_cancer_intersection_models, 
        file = paste0(results_folder, '/raw_cancer_intersection_models.rda'))
saveRDS(raw_cancer_intersection_sig_models, 
        file = paste0(results_folder, '/raw_cancer_intersection_sig_models.rda'))
saveRDS(raw_cancer_union_models, 
        file = paste0(results_folder, '/raw_cancer_union_models.rda'))
saveRDS(raw_cancer_union_sig_models, 
        file = paste0(results_folder, '/raw_cancer_union_sig_models.rda'))
saveRDS(raw_p53_intersection_models, 
        file = paste0(results_folder, '/raw_p53_intersection_models.rda'))
saveRDS(raw_p53_intersection_sig_models, 
        file = paste0(results_folder, '/raw_p53_intersection_sig_models.rda'))
saveRDS(raw_p53_union_models, 
        file = paste0(results_folder, '/raw_p53_union_models.rda'))
saveRDS(raw_p53_union_sig_models, 
        file = paste0(results_folder, '/raw_p53_union_sig_models.rda'))
saveRDS(raw_bal_counts_cancer_intersection_models,
        file = paste0(results_folder, '/raw_bal_counts_cancer_intersection_models.rda'))
saveRDS(raw_bal_counts_cancer_intersection_sig_models,
        file = paste0(results_folder, '/raw_bal_counts_cancer_intersection_sig_models.rda'))
saveRDS(raw_bal_counts_cancer_union_models,
        file = paste0(results_folder, '/raw_bal_counts_cancer_union_models.rda'))
saveRDS(raw_bal_counts_cancer_union_sig_models,
        file = paste0(results_folder, '/raw_bal_counts_cancer_union_sig_models.rda'))
saveRDS(raw_bal_counts_p53_intersection_models,
        file = paste0(results_folder, '/raw_bal_counts_p53_intersection_models.rda'))
saveRDS(raw_bal_counts_p53_intersection_sig_models,
        file = paste0(results_folder, '/raw_bal_counts_p53_intersection_sig_models.rda'))
saveRDS(raw_bal_counts_p53_union_models,
        file = paste0(results_folder, '/raw_bal_counts_p53_union_models.rda'))
saveRDS(raw_bal_counts_p53_union_sig_models,
        file = paste0(results_folder, '/raw_bal_counts_p53_union_sig_models.rda'))
saveRDS(raw_complete_cancer_intersection_models,
        file = paste0(results_folder, '/raw_complete_cancer_intersection_models.rda'))
saveRDS(raw_complete_cancer_intersection_sig_models,
        file = paste0(results_folder, '/raw_complete_cancer_intersection_sig_models.rda'))
saveRDS(raw_complete_cancer_intersection_models,
        file = paste0(results_folder, '/raw_complete_cancer_union_models.rda'))
saveRDS(raw_complete_cancer_intersection_sig_models,
        file = paste0(results_folder, '/raw_complete_cancer_union_sig_models.rda'))
saveRDS(raw_complete_p53_intersection_models,
        file = paste0(results_folder, '/raw_complete_p53_intersection_models.rda'))
saveRDS(raw_complete_p53_intersection_sig_models,
        file = paste0(results_folder, '/raw_complete_p53_intersection_sig_models.rda'))
saveRDS(raw_complete_p53_intersection_models,
        file = paste0(results_folder, '/raw_complete_p53_union_models.rda'))
saveRDS(raw_complete_p53_intersection_sig_models,
        file = paste0(results_folder, '/raw_complete_p53_union_sig_models.rda'))





save.image('/home/benbrew/Desktop/temp_Results.RData')
load('/home/benbrew/Desktop/temp_Results.RData')





















