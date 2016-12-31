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
raw_folder <- paste0(results_folder, '/raw')
swan_folder <- paste0(results_folder, '/swan')
quan_folder <- paste0(results_folder, '/quan')
funnorm_folder <- paste0(results_folder, '/funnorm')
rand_folder <- paste0(results_folder, '/rand')
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

# empty features set
# # complete cancer intersection sig
# raw_complete_cancer_intersection_sig_models <- runModels(beta_raw, 
#                                                          random = F, 
#                                                          bump_hunter = T, 
#                                                          bump_hunter_data = beta_cancer_intersection_sig_features)
# 
# # get table 
# raw_complete_cancer_intersection_sig_table <- extractResults(raw_complete_cancer_intersection_sig_models, 
#                                                              data_name = 'raw_complete_cancer_intersection_sig')


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
        file = paste0(raw_folder, '/raw_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(raw_bal_cancer_models, 
        file = paste0(raw_folder, '/raw_bal_cancer_models.rda'))
saveRDS(raw_bal_cancer_sig_models, 
        file = paste0(raw_folder, '/raw_bal_cancer_sig_models.rda'))
saveRDS(raw_bal_counts_cancer_models, 
        file = paste0(raw_folder, '/raw_bal_counts_cancer_models.rda'))
saveRDS(raw_bal_counts_cancer_sig_models, 
        file = paste0(raw_folder, '/raw_bal_counts_cancer_sig_models.rda'))
saveRDS(raw_unbal_cancer_models, 
        file = paste0(raw_folder, '/raw_unbal_cancer_models.rda'))
saveRDS(raw_unbal_cancer_sig_models, 
        file = paste0(raw_folder, '/raw_unbal_cancer_sig_models.rda'))
saveRDS(raw_bal_p53_models, 
        file = paste0(raw_folder, '/raw_bal_p53_models.rda'))
saveRDS(raw_bal_p53_sig_models, 
        file = paste0(raw_folder, '/raw_bal_p53_sig_models.rda'))
saveRDS(raw_bal_counts_p53_models, 
        file = paste0(raw_folder, '/raw_bal_counts_p53_models.rda'))
saveRDS(raw_bal_counts_p53_sig_models, 
        file = paste0(raw_folder, '/raw_bal_counts_p53_sig_models.rda'))
saveRDS(raw_unbal_p53_models, 
        file = paste0(raw_folder, '/raw_unbal_p53_models.rda'))
saveRDS(raw_unbal_p53_sig_models, 
        file = paste0(raw_folder, '/raw_unbal_p53_sig_models.rda'))
saveRDS(raw_cancer_intersection_models, 
        file = paste0(raw_folder, '/raw_cancer_intersection_models.rda'))
saveRDS(raw_cancer_intersection_sig_models, 
        file = paste0(raw_folder, '/raw_cancer_intersection_sig_models.rda'))
saveRDS(raw_cancer_union_models, 
        file = paste0(raw_folder, '/raw_cancer_union_models.rda'))
saveRDS(raw_cancer_union_sig_models, 
        file = paste0(raw_folder, '/raw_cancer_union_sig_models.rda'))
saveRDS(raw_p53_intersection_models, 
        file = paste0(raw_folder, '/raw_p53_intersection_models.rda'))
saveRDS(raw_p53_intersection_sig_models, 
        file = paste0(raw_folder, '/raw_p53_intersection_sig_models.rda'))
saveRDS(raw_p53_union_models, 
        file = paste0(raw_folder, '/raw_p53_union_models.rda'))
saveRDS(raw_p53_union_sig_models, 
        file = paste0(raw_folder, '/raw_p53_union_sig_models.rda'))
saveRDS(raw_bal_counts_cancer_intersection_models,
        file = paste0(raw_folder, '/raw_bal_counts_cancer_intersection_models.rda'))
saveRDS(raw_bal_counts_cancer_intersection_sig_models,
        file = paste0(raw_folder, '/raw_bal_counts_cancer_intersection_sig_models.rda'))
saveRDS(raw_bal_counts_cancer_union_models,
        file = paste0(raw_folder, '/raw_bal_counts_cancer_union_models.rda'))
saveRDS(raw_bal_counts_cancer_union_sig_models,
        file = paste0(raw_folder, '/raw_bal_counts_cancer_union_sig_models.rda'))
saveRDS(raw_bal_counts_p53_intersection_models,
        file = paste0(raw_folder, '/raw_bal_counts_p53_intersection_models.rda'))
saveRDS(raw_bal_counts_p53_intersection_sig_models,
        file = paste0(raw_folder, '/raw_bal_counts_p53_intersection_sig_models.rda'))
saveRDS(raw_bal_counts_p53_union_models,
        file = paste0(raw_folder, '/raw_bal_counts_p53_union_models.rda'))
saveRDS(raw_bal_counts_p53_union_sig_models,
        file = paste0(raw_folder, '/raw_bal_counts_p53_union_sig_models.rda'))
saveRDS(raw_complete_cancer_intersection_models,
        file = paste0(raw_folder, '/raw_complete_cancer_intersection_models.rda'))
# saveRDS(raw_complete_cancer_intersection_sig_models,
#         file = paste0(raw_folder, '/raw_complete_cancer_intersection_sig_models.rda'))
saveRDS(raw_complete_cancer_union_models,
        file = paste0(raw_folder, '/raw_complete_cancer_union_models.rda'))
saveRDS(raw_complete_cancer_union_sig_models,
        file = paste0(raw_folder, '/raw_complete_cancer_union_sig_models.rda'))
saveRDS(raw_complete_p53_intersection_models,
        file = paste0(raw_folder, '/raw_complete_p53_intersection_models.rda'))
saveRDS(raw_complete_p53_intersection_sig_models,
        file = paste0(raw_folder, '/raw_complete_p53_intersection_sig_models.rda'))
saveRDS(raw_complete_p53_union_models,
        file = paste0(raw_folder, '/raw_complete_p53_union_models.rda'))
saveRDS(raw_complete_p53_union_sig_models,
        file = paste0(raw_folder, '/raw_complete_p53_union_sig_models.rda'))

###################################################################################################################################
## beta_swan

##########
# cancer
##########

# bal cancer
swan_bal_cancer_models <- runModels(beta_swan, 
                                   random = F, 
                                   bump_hunter = T, 
                                   bump_hunter_data = beta_swan_bal_cancer_features)

swan_bal_cancer_table <- extractResults(swan_bal_cancer_models, 
                                       data_name = 'beta_swan_bal_cancer')

# bal cancer sig
swan_bal_cancer_sig_models <- runModels(beta_swan, 
                                       random = F, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_swan_bal_cancer_sig_features)

swan_bal_cancer_sig_table <- extractResults(swan_bal_cancer_sig_models, 
                                           data_name = 'swan_bal_cancer_sig')


# bal counts cancer
swan_bal_counts_cancer_models <- runModels(beta_swan, 
                                          random = F, 
                                          bump_hunter = T, 
                                          bump_hunter_data = beta_swan_bal_counts_cancer_features)

swan_bal_counts_cancer_table <- extractResults(swan_bal_counts_cancer_models, 
                                              data_name = 'swan_bal_counts_cancer')

# bal counts cancer sig
swan_bal_counts_cancer_sig_models <- runModels(beta_swan, 
                                              random = F, 
                                              bump_hunter = T, 
                                              bump_hunter_data =  beta_swan_bal_counts_cancer_sig_features)

swan_bal_counts_cancer_sig_table <- extractResults(swan_bal_counts_cancer_sig_models, 
                                                  data_name = 'swan_bal_counts_cancer_sig')

# unbal cancer
swan_unbal_cancer_models <- runModels(beta_swan, 
                                     random = F, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_swan_unbal_cancer_features)

swan_unbal_cancer_table <- extractResults(swan_unbal_cancer_models, 
                                         data_name = 'swan_unbal_cancer')


# unbal cancer sig
swan_unbal_cancer_sig_models <- runModels(beta_swan, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_swan_unbal_cancer_sig_features)

swan_unbal_cancer_sig_table <- extractResults(swan_unbal_cancer_sig_models, 
                                             data_name = 'swan_unbal_cancer_sig')

##########
# p53
##########

# bal p53
swan_bal_p53_models <- runModels(beta_swan, 
                                random = F, 
                                bump_hunter = T, 
                                bump_hunter_data = beta_swan_bal_p53_features)

swan_bal_p53_table <- extractResults(swan_bal_p53_models, 
                                    data_name = 'swan_bal_p53')

# bal p53 sig
swan_bal_p53_sig_models <- runModels(beta_swan, 
                                    random = F, 
                                    bump_hunter = T, 
                                    bump_hunter_data = beta_swan_bal_p53_sig_features)

swan_bal_p53_sig_table <- extractResults(swan_bal_p53_sig_models, 
                                        data_name = 'swan_bal_p53_sig')

# bal counts p53
swan_bal_counts_p53_models <- runModels(beta_swan, 
                                       random = F, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_swan_bal_counts_p53_features)

swan_bal_counts_p53_table <- extractResults(swan_bal_counts_p53_models, 
                                           data_name = 'swan_bal_counts_p53')

# bal counts p53 sig
swan_bal_counts_p53_sig_models <- runModels(beta_swan, 
                                           random = F, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_swan_bal_counts_p53_sig_features)

swan_bal_counts_p53_sig_table <- extractResults(swan_bal_counts_p53_sig_models, 
                                               data_name = 'swan_bal_counts_p53_sig')

# unbal p53
swan_unbal_p53_models <- runModels(beta_swan, 
                                  random = F, 
                                  bump_hunter = T, 
                                  bump_hunter_data = beta_swan_unbal_p53_features)

swan_unbal_p53_table <- extractResults(swan_unbal_p53_models, 
                                      data_name = 'swan_unbal_p53')


# unbal p53 sig
swan_unbal_p53_sig_models <- runModels(beta_swan, 
                                      random = F, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_swan_unbal_p53_sig_features)

swan_unbal_p53_sig_table <- extractResults(swan_unbal_p53_sig_models, 
                                          data_name = 'swan_unbal_p53_sig')


##########
# cancer
##########

# cancer_intersection
swan_cancer_intersection_models <- runModels(beta_swan, 
                                            random = F, 
                                            bump_hunter = T, 
                                            bump_hunter_data = beta_swan_cancer_intersection_features)

swan_cancer_intersection_table <- extractResults(swan_cancer_intersection_models, 
                                                data_name = 'swan_cancer_intersection')

# cancer_intersection sig
swan_cancer_intersection_sig_models <- runModels(beta_swan, 
                                                random = F, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_swan_cancer_intersection_sig_features)

swan_cancer_intersection_sig_table <- extractResults(swan_cancer_intersection_sig_models, 
                                                    data_name = 'swan_cancer_intersection_sig')

# cancer_union
swan_cancer_union_models <- runModels(beta_swan, 
                                     random = F, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_swan_cancer_union_features)

swan_cancer_union_table <- extractResults(swan_cancer_union_models, 
                                         data_name = 'swan_cancer_union')

# cancer_union sig
swan_cancer_union_sig_models <- runModels(beta_swan, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_swan_cancer_union_sig_features)

swan_cancer_union_sig_table <- extractResults(swan_cancer_union_sig_models, 
                                             data_name = 'swan_cancer_union_sig')
##########
# p53
##########
# p53_intersection
swan_p53_intersection_models <- runModels(beta_swan, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_swan_p53_intersection_features)

swan_p53_intersection_table <- extractResults(swan_p53_intersection_models, 
                                             data_name = 'swan_p53_intersection')

# p53_intersection sig
swan_p53_intersection_sig_models <- runModels(beta_swan, 
                                             random = F, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_swan_p53_intersection_sig_features)

swan_p53_intersection_sig_table <- extractResults(swan_p53_intersection_sig_models, 
                                                 data_name = 'swan_p53_intersection_sig')

# p53_union
swan_p53_union_models <- runModels(beta_swan, 
                                  random = F, 
                                  bump_hunter = T, 
                                  bump_hunter_data = beta_swan_p53_union_features)

swan_p53_union_table <- extractResults(swan_p53_union_models, 
                                      data_name = 'swan_p53_union')

# p53_union sig
swan_p53_union_sig_models <- runModels(beta_swan, 
                                      random = F, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_swan_p53_union_sig_features)

swan_p53_union_sig_table <- extractResults(swan_p53_union_sig_models, 
                                          data_name = 'swan_p53_union_sig')


##########
# most balanced cancer 
##########

# intersection
swan_bal_counts_cancer_intersection_models <- runModels(beta_swan, 
                                                       random = F, 
                                                       bump_hunter = T, 
                                                       bump_hunter_data = beta_bal_counts_cancer_intersection_features)

swan_bal_counts_cancer_intersection_table <- extractResults(swan_bal_counts_cancer_intersection_models, 
                                                           data_name = 'swan_bal_counts_cancer_intersection')


swan_bal_counts_cancer_intersection_sig_models <- runModels(beta_swan, 
                                                           random = F, 
                                                           bump_hunter = T, 
                                                           bump_hunter_data = beta_bal_counts_cancer_intersection_sig_features)

swan_bal_counts_cancer_intersection_sig_table <- extractResults(swan_bal_counts_cancer_intersection_sig_models, 
                                                               data_name = 'swan_bal_counts_cancer_intersection_sig')

# union
swan_bal_counts_cancer_union_models <- runModels(beta_swan, 
                                                random = F, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_bal_counts_cancer_union_features)

swan_bal_counts_cancer_union_table <- extractResults(swan_bal_counts_cancer_union_models, 
                                                    data_name = 'swan_bal_counts_cancer_union')


swan_bal_counts_cancer_union_sig_models <- runModels(beta_swan, 
                                                    random = F, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_cancer_union_sig_features)

swan_bal_counts_cancer_union_sig_table <- extractResults(swan_bal_counts_cancer_union_sig_models, 
                                                        data_name = 'swan_bal_counts_cancer_union_sig')

##########
# most balanced p53 
##########

# intersection
swan_bal_counts_p53_intersection_models <- runModels(beta_swan, 
                                                    random = F, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_p53_intersection_features)

swan_bal_counts_p53_intersection_table <- extractResults(swan_bal_counts_p53_intersection_models, 
                                                        data_name = 'swan_bal_counts_p53_intersection')


swan_bal_counts_p53_intersection_sig_models <- runModels(beta_swan, 
                                                        random = F, 
                                                        bump_hunter = T, 
                                                        bump_hunter_data = beta_bal_counts_p53_intersection_sig_features)

swan_bal_counts_p53_intersection_sig_table <- extractResults(swan_bal_counts_p53_intersection_sig_models, 
                                                            data_name = 'swan_bal_counts_p53_intersection_sig')


# union
swan_bal_counts_p53_union_models <- runModels(beta_swan, 
                                             random = F, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_bal_counts_p53_union_features)

swan_bal_counts_p53_union_table <- extractResults(swan_bal_counts_p53_union_models, 
                                                 data_name = 'swan_bal_counts_p53_union')


swan_bal_counts_p53_union_sig_models <- runModels(beta_swan, 
                                                 random = F, 
                                                 bump_hunter = T, 
                                                 bump_hunter_data = beta_bal_counts_p53_union_sig_features)

swan_bal_counts_p53_union_sig_table <- extractResults(swan_bal_counts_p53_union_sig_models, 
                                                     data_name = 'swan_bal_counts_p53_union_sig')

###########
# complete cancer - intersection across each method intersection
###########

# complete cancer intersection
swan_complete_cancer_intersection_models <- runModels(beta_swan, 
                                                     random = F, 
                                                     bump_hunter = T, 
                                                     bump_hunter_data = beta_cancer_intersection_features)

# get table 
swan_complete_cancer_intersection_table <- extractResults(swan_complete_cancer_intersection_models, 
                                                         data_name = 'swan_complete_cancer_intersection')

# empty features set
# # complete cancer intersection sig
# swan_complete_cancer_intersection_sig_models <- runModels(beta_swan, 
#                                                          random = F, 
#                                                          bump_hunter = T, 
#                                                          bump_hunter_data = beta_cancer_intersection_sig_features)
# 
# # get table 
# swan_complete_cancer_intersection_sig_table <- extractResults(swan_complete_cancer_intersection_sig_models, 
#                                                              data_name = 'swan_complete_cancer_intersection_sig')


# complete cancer union
swan_complete_cancer_union_models <- runModels(beta_swan, 
                                              random = F, 
                                              bump_hunter = T, 
                                              bump_hunter_data = beta_cancer_union_features)

# get table 
swan_complete_cancer_union_table <- extractResults(swan_complete_cancer_union_models, 
                                                  data_name = 'swan_complete_cancer_union')

# complete cancer union sig
swan_complete_cancer_union_sig_models <- runModels(beta_swan, 
                                                  random = F, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_cancer_union_sig_features)

# get table 
swan_complete_cancer_union_sig_table <- extractResults(swan_complete_cancer_union_sig_models, 
                                                      data_name = 'swan_complete_cancer_union_sig')


###########
# complete p53 - intersection across each method intersection
###########

# complete p53 intersection
swan_complete_p53_intersection_models <- runModels(beta_swan, 
                                                  random = F, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_p53_intersection_features)

# get table 
swan_complete_p53_intersection_table <- extractResults(swan_complete_p53_intersection_models, 
                                                      data_name = 'swan_complete_p53_intersection')

# complete p53 intersection sig
swan_complete_p53_intersection_sig_models <- runModels(beta_swan, 
                                                      random = F, 
                                                      bump_hunter = T, 
                                                      bump_hunter_data = beta_p53_intersection_sig_features)

# get table 
swan_complete_p53_intersection_sig_table <- extractResults(swan_complete_p53_intersection_sig_models, 
                                                          data_name = 'swan_complete_p53_intersection_sig')


# complete p53 union
swan_complete_p53_union_models <- runModels(beta_swan, 
                                           random = F, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_p53_union_features)

# get table 
swan_complete_p53_union_table <- extractResults(swan_complete_p53_union_models, 
                                               data_name = 'swan_complete_p53_union')

# complete p53 union sig
swan_complete_p53_union_sig_models <- runModels(beta_swan, 
                                               random = F, 
                                               bump_hunter = T, 
                                               bump_hunter_data = beta_p53_union_sig_features)

# get table 
swan_complete_p53_union_sig_table <- extractResults(swan_complete_p53_union_sig_models, 
                                                   data_name = 'swan_complete_p53_union_sig')


###########
# rbind tables and save RDA file
###########
swan_table <- rbind(swan_bal_cancer_table, swan_bal_cancer_sig_table, swan_bal_counts_cancer_table, swan_bal_counts_cancer_sig_table,
                   swan_unbal_cancer_table, swan_unbal_cancer_sig_table, swan_bal_p53_table, swan_bal_p53_sig_table, swan_bal_counts_p53_table, 
                   swan_bal_counts_p53_sig_table, swan_unbal_p53_table, swan_unbal_p53_sig_table, swan_cancer_intersection_table, 
                   swan_cancer_intersection_sig_table, swan_cancer_union_table, swan_cancer_union_sig_table, swan_p53_intersection_table, 
                   swan_p53_intersection_sig_table, swan_p53_union_table, swan_p53_union_sig_table,
                   swan_bal_counts_cancer_intersection_table, swan_bal_counts_cancer_intersection_sig_table,
                   swan_bal_counts_cancer_union_table, swan_bal_counts_cancer_union_sig_table,
                   swan_bal_counts_p53_intersection_table, swan_bal_counts_p53_intersection_sig_table,
                   swan_bal_counts_p53_union_table, swan_bal_counts_p53_union_sig_table,
                   swan_complete_cancer_intersection_table, swan_complete_cancer_intersection_sig_table,
                   swan_complete_cancer_union_table, swan_complete_cancer_union_sig_table,
                   swan_complete_p53_intersection_table, swan_complete_p53_intersection_sig_table,
                   swan_complete_p53_union_table, swan_complete_p53_union_sig_table)

# remove data 
rm(swan_bal_cancer_table, swan_bal_cancer_sig_table, swan_bal_counts_cancer_table, swan_bal_counts_cancer_sig_table,
   swan_unbal_cancer_table, swan_unbal_cancer_sig_table, swan_bal_p53_table, swan_bal_p53_sig_table, swan_bal_counts_p53_table, 
   swan_bal_counts_p53_sig_table, swan_unbal_p53_table, swan_unbal_p53_sig_table, swan_cancer_intersection_table, 
   swan_cancer_intersection_sig_table, swan_cancer_union_table, swan_cancer_union_sig_table, swan_p53_intersection_table, 
   swan_p53_intersection_sig_table, swan_p53_union_table, swan_p53_union_sig_table,
   swan_bal_counts_cancer_intersection_table, swan_bal_counts_cancer_intersection_sig_table,
   swan_bal_counts_cancer_union_table, swan_bal_counts_cancer_union_sig_table,
   swan_bal_counts_p53_intersection_table, swan_bal_counts_p53_intersection_sig_table,
   swan_bal_counts_p53_union_table, swan_bal_counts_p53_union_sig_table,
   swan_complete_cancer_intersection_table, swan_complete_cancer_intersection_sig_table,
   swan_complete_cancer_union_table, swan_complete_cancer_union_sig_table,
   swan_complete_p53_intersection_table, swan_complete_p53_intersection_sig_table,
   swan_complete_p53_union_table, swan_complete_p53_union_sig_table)


#save table 
saveRDS(swan_table, 
        file = paste0(swan_folder, '/swan_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(swan_bal_cancer_models, 
        file = paste0(swan_folder, '/swan_bal_cancer_models.rda'))
saveRDS(swan_bal_cancer_sig_models, 
        file = paste0(swan_folder, '/swan_bal_cancer_sig_models.rda'))
saveRDS(swan_bal_counts_cancer_models, 
        file = paste0(swan_folder, '/swan_bal_counts_cancer_models.rda'))
saveRDS(swan_bal_counts_cancer_sig_models, 
        file = paste0(swan_folder, '/swan_bal_counts_cancer_sig_models.rda'))
saveRDS(swan_unbal_cancer_models, 
        file = paste0(swan_folder, '/swan_unbal_cancer_models.rda'))
saveRDS(swan_unbal_cancer_sig_models, 
        file = paste0(swan_folder, '/swan_unbal_cancer_sig_models.rda'))
saveRDS(swan_bal_p53_models, 
        file = paste0(swan_folder, '/swan_bal_p53_models.rda'))
saveRDS(swan_bal_p53_sig_models, 
        file = paste0(swan_folder, '/swan_bal_p53_sig_models.rda'))
saveRDS(swan_bal_counts_p53_models, 
        file = paste0(swan_folder, '/swan_bal_counts_p53_models.rda'))
saveRDS(swan_bal_counts_p53_sig_models, 
        file = paste0(swan_folder, '/swan_bal_counts_p53_sig_models.rda'))
saveRDS(swan_unbal_p53_models, 
        file = paste0(swan_folder, '/swan_unbal_p53_models.rda'))
saveRDS(swan_unbal_p53_sig_models, 
        file = paste0(swan_folder, '/swan_unbal_p53_sig_models.rda'))
saveRDS(swan_cancer_intersection_models, 
        file = paste0(swan_folder, '/swan_cancer_intersection_models.rda'))
saveRDS(swan_cancer_intersection_sig_models, 
        file = paste0(swan_folder, '/swan_cancer_intersection_sig_models.rda'))
saveRDS(swan_cancer_union_models, 
        file = paste0(swan_folder, '/swan_cancer_union_models.rda'))
saveRDS(swan_cancer_union_sig_models, 
        file = paste0(swan_folder, '/swan_cancer_union_sig_models.rda'))
saveRDS(swan_p53_intersection_models, 
        file = paste0(swan_folder, '/swan_p53_intersection_models.rda'))
saveRDS(swan_p53_intersection_sig_models, 
        file = paste0(swan_folder, '/swan_p53_intersection_sig_models.rda'))
saveRDS(swan_p53_union_models, 
        file = paste0(swan_folder, '/swan_p53_union_models.rda'))
saveRDS(swan_p53_union_sig_models, 
        file = paste0(swan_folder, '/swan_p53_union_sig_models.rda'))
saveRDS(swan_bal_counts_cancer_intersection_models,
        file = paste0(swan_folder, '/swan_bal_counts_cancer_intersection_models.rda'))
saveRDS(swan_bal_counts_cancer_intersection_sig_models,
        file = paste0(swan_folder, '/swan_bal_counts_cancer_intersection_sig_models.rda'))
saveRDS(swan_bal_counts_cancer_union_models,
        file = paste0(swan_folder, '/swan_bal_counts_cancer_union_models.rda'))
saveRDS(swan_bal_counts_cancer_union_sig_models,
        file = paste0(swan_folder, '/swan_bal_counts_cancer_union_sig_models.rda'))
saveRDS(swan_bal_counts_p53_intersection_models,
        file = paste0(swan_folder, '/swan_bal_counts_p53_intersection_models.rda'))
saveRDS(swan_bal_counts_p53_intersection_sig_models,
        file = paste0(swan_folder, '/swan_bal_counts_p53_intersection_sig_models.rda'))
saveRDS(swan_bal_counts_p53_union_models,
        file = paste0(swan_folder, '/swan_bal_counts_p53_union_models.rda'))
saveRDS(swan_bal_counts_p53_union_sig_models,
        file = paste0(swan_folder, '/swan_bal_counts_p53_union_sig_models.rda'))
saveRDS(swan_complete_cancer_intersection_models,
        file = paste0(swan_folder, '/swan_complete_cancer_intersection_models.rda'))
# saveRDS(swan_complete_cancer_intersection_sig_models,
#         file = paste0(swan_folder, '/swan_complete_cancer_intersection_sig_models.rda'))
saveRDS(swan_complete_cancer_union_models,
        file = paste0(swan_folder, '/swan_complete_cancer_union_models.rda'))
saveRDS(swan_complete_cancer_union_sig_models,
        file = paste0(swan_folder, '/swan_complete_cancer_union_sig_models.rda'))
saveRDS(swan_complete_p53_intersection_models,
        file = paste0(swan_folder, '/swan_complete_p53_intersection_models.rda'))
saveRDS(swan_complete_p53_intersection_sig_models,
        file = paste0(swan_folder, '/swan_complete_p53_intersection_sig_models.rda'))
saveRDS(swan_complete_p53_union_models,
        file = paste0(swan_folder, '/swan_complete_p53_union_models.rda'))
saveRDS(swan_complete_p53_union_sig_models,
        file = paste0(swan_folder, '/swan_complete_p53_union_sig_models.rda'))


###################################################################################################################################
## beta_quan

##########
# cancer
##########

# bal cancer
quan_bal_cancer_models <- runModels(beta_quan, 
                                   random = F, 
                                   bump_hunter = T, 
                                   bump_hunter_data = beta_quan_bal_cancer_features)

quan_bal_cancer_table <- extractResults(quan_bal_cancer_models, 
                                       data_name = 'beta_quan_bal_cancer')

# bal cancer sig
quan_bal_cancer_sig_models <- runModels(beta_quan, 
                                       random = F, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_quan_bal_cancer_sig_features)

quan_bal_cancer_sig_table <- extractResults(quan_bal_cancer_sig_models, 
                                           data_name = 'quan_bal_cancer_sig')


# bal counts cancer
quan_bal_counts_cancer_models <- runModels(beta_quan, 
                                          random = F, 
                                          bump_hunter = T, 
                                          bump_hunter_data = beta_quan_bal_counts_cancer_features)

quan_bal_counts_cancer_table <- extractResults(quan_bal_counts_cancer_models, 
                                              data_name = 'quan_bal_counts_cancer')

# bal counts cancer sig
quan_bal_counts_cancer_sig_models <- runModels(beta_quan, 
                                              random = F, 
                                              bump_hunter = T, 
                                              bump_hunter_data =  beta_quan_bal_counts_cancer_sig_features)

quan_bal_counts_cancer_sig_table <- extractResults(quan_bal_counts_cancer_sig_models, 
                                                  data_name = 'quan_bal_counts_cancer_sig')

# unbal cancer
quan_unbal_cancer_models <- runModels(beta_quan, 
                                     random = F, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_quan_unbal_cancer_features)

quan_unbal_cancer_table <- extractResults(quan_unbal_cancer_models, 
                                         data_name = 'quan_unbal_cancer')


# unbal cancer sig
quan_unbal_cancer_sig_models <- runModels(beta_quan, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_quan_unbal_cancer_sig_features)

quan_unbal_cancer_sig_table <- extractResults(quan_unbal_cancer_sig_models, 
                                             data_name = 'quan_unbal_cancer_sig')

##########
# p53
##########

# bal p53
quan_bal_p53_models <- runModels(beta_quan, 
                                random = F, 
                                bump_hunter = T, 
                                bump_hunter_data = beta_quan_bal_p53_features)

quan_bal_p53_table <- extractResults(quan_bal_p53_models, 
                                    data_name = 'quan_bal_p53')

# bal p53 sig
quan_bal_p53_sig_models <- runModels(beta_quan, 
                                    random = F, 
                                    bump_hunter = T, 
                                    bump_hunter_data = beta_quan_bal_p53_sig_features)

quan_bal_p53_sig_table <- extractResults(quan_bal_p53_sig_models, 
                                        data_name = 'quan_bal_p53_sig')

# bal counts p53
quan_bal_counts_p53_models <- runModels(beta_quan, 
                                       random = F, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_quan_bal_counts_p53_features)

quan_bal_counts_p53_table <- extractResults(quan_bal_counts_p53_models, 
                                           data_name = 'quan_bal_counts_p53')

# bal counts p53 sig
quan_bal_counts_p53_sig_models <- runModels(beta_quan, 
                                           random = F, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_quan_bal_counts_p53_sig_features)

quan_bal_counts_p53_sig_table <- extractResults(quan_bal_counts_p53_sig_models, 
                                               data_name = 'quan_bal_counts_p53_sig')

# unbal p53
quan_unbal_p53_models <- runModels(beta_quan, 
                                  random = F, 
                                  bump_hunter = T, 
                                  bump_hunter_data = beta_quan_unbal_p53_features)

quan_unbal_p53_table <- extractResults(quan_unbal_p53_models, 
                                      data_name = 'quan_unbal_p53')


# unbal p53 sig
quan_unbal_p53_sig_models <- runModels(beta_quan, 
                                      random = F, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_quan_unbal_p53_sig_features)

quan_unbal_p53_sig_table <- extractResults(quan_unbal_p53_sig_models, 
                                          data_name = 'quan_unbal_p53_sig')


##########
# cancer
##########

# cancer_intersection
quan_cancer_intersection_models <- runModels(beta_quan, 
                                            random = F, 
                                            bump_hunter = T, 
                                            bump_hunter_data = beta_quan_cancer_intersection_features)

quan_cancer_intersection_table <- extractResults(quan_cancer_intersection_models, 
                                                data_name = 'quan_cancer_intersection')

# cancer_intersection sig
quan_cancer_intersection_sig_models <- runModels(beta_quan, 
                                                random = F, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_quan_cancer_intersection_sig_features)

quan_cancer_intersection_sig_table <- extractResults(quan_cancer_intersection_sig_models, 
                                                    data_name = 'quan_cancer_intersection_sig')

# cancer_union
quan_cancer_union_models <- runModels(beta_quan, 
                                     random = F, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_quan_cancer_union_features)

quan_cancer_union_table <- extractResults(quan_cancer_union_models, 
                                         data_name = 'quan_cancer_union')

# cancer_union sig
quan_cancer_union_sig_models <- runModels(beta_quan, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_quan_cancer_union_sig_features)

quan_cancer_union_sig_table <- extractResults(quan_cancer_union_sig_models, 
                                             data_name = 'quan_cancer_union_sig')
##########
# p53
##########
# p53_intersection
quan_p53_intersection_models <- runModels(beta_quan, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_quan_p53_intersection_features)

quan_p53_intersection_table <- extractResults(quan_p53_intersection_models, 
                                             data_name = 'quan_p53_intersection')

# p53_intersection sig
quan_p53_intersection_sig_models <- runModels(beta_quan, 
                                             random = F, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_quan_p53_intersection_sig_features)

quan_p53_intersection_sig_table <- extractResults(quan_p53_intersection_sig_models, 
                                                 data_name = 'quan_p53_intersection_sig')

# p53_union
quan_p53_union_models <- runModels(beta_quan, 
                                  random = F, 
                                  bump_hunter = T, 
                                  bump_hunter_data = beta_quan_p53_union_features)

quan_p53_union_table <- extractResults(quan_p53_union_models, 
                                      data_name = 'quan_p53_union')

# p53_union sig
quan_p53_union_sig_models <- runModels(beta_quan, 
                                      random = F, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_quan_p53_union_sig_features)

quan_p53_union_sig_table <- extractResults(quan_p53_union_sig_models, 
                                          data_name = 'quan_p53_union_sig')


##########
# most balanced cancer 
##########

# intersection
quan_bal_counts_cancer_intersection_models <- runModels(beta_quan, 
                                                       random = F, 
                                                       bump_hunter = T, 
                                                       bump_hunter_data = beta_bal_counts_cancer_intersection_features)

quan_bal_counts_cancer_intersection_table <- extractResults(quan_bal_counts_cancer_intersection_models, 
                                                           data_name = 'quan_bal_counts_cancer_intersection')


quan_bal_counts_cancer_intersection_sig_models <- runModels(beta_quan, 
                                                           random = F, 
                                                           bump_hunter = T, 
                                                           bump_hunter_data = beta_bal_counts_cancer_intersection_sig_features)

quan_bal_counts_cancer_intersection_sig_table <- extractResults(quan_bal_counts_cancer_intersection_sig_models, 
                                                               data_name = 'quan_bal_counts_cancer_intersection_sig')

# union
quan_bal_counts_cancer_union_models <- runModels(beta_quan, 
                                                random = F, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_bal_counts_cancer_union_features)

quan_bal_counts_cancer_union_table <- extractResults(quan_bal_counts_cancer_union_models, 
                                                    data_name = 'quan_bal_counts_cancer_union')


quan_bal_counts_cancer_union_sig_models <- runModels(beta_quan, 
                                                    random = F, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_cancer_union_sig_features)

quan_bal_counts_cancer_union_sig_table <- extractResults(quan_bal_counts_cancer_union_sig_models, 
                                                        data_name = 'quan_bal_counts_cancer_union_sig')

##########
# most balanced p53 
##########

# intersection
quan_bal_counts_p53_intersection_models <- runModels(beta_quan, 
                                                    random = F, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_p53_intersection_features)

quan_bal_counts_p53_intersection_table <- extractResults(quan_bal_counts_p53_intersection_models, 
                                                        data_name = 'quan_bal_counts_p53_intersection')


quan_bal_counts_p53_intersection_sig_models <- runModels(beta_quan, 
                                                        random = F, 
                                                        bump_hunter = T, 
                                                        bump_hunter_data = beta_bal_counts_p53_intersection_sig_features)

quan_bal_counts_p53_intersection_sig_table <- extractResults(quan_bal_counts_p53_intersection_sig_models, 
                                                            data_name = 'quan_bal_counts_p53_intersection_sig')


# union
quan_bal_counts_p53_union_models <- runModels(beta_quan, 
                                             random = F, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_bal_counts_p53_union_features)

quan_bal_counts_p53_union_table <- extractResults(quan_bal_counts_p53_union_models, 
                                                 data_name = 'quan_bal_counts_p53_union')


quan_bal_counts_p53_union_sig_models <- runModels(beta_quan, 
                                                 random = F, 
                                                 bump_hunter = T, 
                                                 bump_hunter_data = beta_bal_counts_p53_union_sig_features)

quan_bal_counts_p53_union_sig_table <- extractResults(quan_bal_counts_p53_union_sig_models, 
                                                     data_name = 'quan_bal_counts_p53_union_sig')

###########
# complete cancer - intersection across each method intersection
###########

# complete cancer intersection
quan_complete_cancer_intersection_models <- runModels(beta_quan, 
                                                     random = F, 
                                                     bump_hunter = T, 
                                                     bump_hunter_data = beta_cancer_intersection_features)

# get table 
quan_complete_cancer_intersection_table <- extractResults(quan_complete_cancer_intersection_models, 
                                                         data_name = 'quan_complete_cancer_intersection')

# empty features set
# # complete cancer intersection sig
# quan_complete_cancer_intersection_sig_models <- runModels(beta_quan, 
#                                                          random = F, 
#                                                          bump_hunter = T, 
#                                                          bump_hunter_data = beta_cancer_intersection_sig_features)
# 
# # get table 
# quan_complete_cancer_intersection_sig_table <- extractResults(quan_complete_cancer_intersection_sig_models, 
#                                                              data_name = 'quan_complete_cancer_intersection_sig')


# complete cancer union
quan_complete_cancer_union_models <- runModels(beta_quan, 
                                              random = F, 
                                              bump_hunter = T, 
                                              bump_hunter_data = beta_cancer_union_features)

# get table 
quan_complete_cancer_union_table <- extractResults(quan_complete_cancer_union_models, 
                                                  data_name = 'quan_complete_cancer_union')

# complete cancer union sig
quan_complete_cancer_union_sig_models <- runModels(beta_quan, 
                                                  random = F, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_cancer_union_sig_features)

# get table 
quan_complete_cancer_union_sig_table <- extractResults(quan_complete_cancer_union_sig_models, 
                                                      data_name = 'quan_complete_cancer_union_sig')


###########
# complete p53 - intersection across each method intersection
###########

# complete p53 intersection
quan_complete_p53_intersection_models <- runModels(beta_quan, 
                                                  random = F, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_p53_intersection_features)

# get table 
quan_complete_p53_intersection_table <- extractResults(quan_complete_p53_intersection_models, 
                                                      data_name = 'quan_complete_p53_intersection')

# complete p53 intersection sig
quan_complete_p53_intersection_sig_models <- runModels(beta_quan, 
                                                      random = F, 
                                                      bump_hunter = T, 
                                                      bump_hunter_data = beta_p53_intersection_sig_features)

# get table 
quan_complete_p53_intersection_sig_table <- extractResults(quan_complete_p53_intersection_sig_models, 
                                                          data_name = 'quan_complete_p53_intersection_sig')


# complete p53 union
quan_complete_p53_union_models <- runModels(beta_quan, 
                                           random = F, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_p53_union_features)

# get table 
quan_complete_p53_union_table <- extractResults(quan_complete_p53_union_models, 
                                               data_name = 'quan_complete_p53_union')

# complete p53 union sig
quan_complete_p53_union_sig_models <- runModels(beta_quan, 
                                               random = F, 
                                               bump_hunter = T, 
                                               bump_hunter_data = beta_p53_union_sig_features)

# get table 
quan_complete_p53_union_sig_table <- extractResults(quan_complete_p53_union_sig_models, 
                                                   data_name = 'quan_complete_p53_union_sig')


###########
# rbind tables and save RDA file
###########
quan_table <- rbind(quan_bal_cancer_table, quan_bal_cancer_sig_table, quan_bal_counts_cancer_table, quan_bal_counts_cancer_sig_table,
                   quan_unbal_cancer_table, quan_unbal_cancer_sig_table, quan_bal_p53_table, quan_bal_p53_sig_table, quan_bal_counts_p53_table, 
                   quan_bal_counts_p53_sig_table, quan_unbal_p53_table, quan_unbal_p53_sig_table, quan_cancer_intersection_table, 
                   quan_cancer_intersection_sig_table, quan_cancer_union_table, quan_cancer_union_sig_table, quan_p53_intersection_table, 
                   quan_p53_intersection_sig_table, quan_p53_union_table, quan_p53_union_sig_table,
                   quan_bal_counts_cancer_intersection_table, quan_bal_counts_cancer_intersection_sig_table,
                   quan_bal_counts_cancer_union_table, quan_bal_counts_cancer_union_sig_table,
                   quan_bal_counts_p53_intersection_table, quan_bal_counts_p53_intersection_sig_table,
                   quan_bal_counts_p53_union_table, quan_bal_counts_p53_union_sig_table,
                   quan_complete_cancer_intersection_table, quan_complete_cancer_intersection_sig_table,
                   quan_complete_cancer_union_table, quan_complete_cancer_union_sig_table,
                   quan_complete_p53_intersection_table, quan_complete_p53_intersection_sig_table,
                   quan_complete_p53_union_table, quan_complete_p53_union_sig_table)

# remove data 
rm(quan_bal_cancer_table, quan_bal_cancer_sig_table, quan_bal_counts_cancer_table, quan_bal_counts_cancer_sig_table,
   quan_unbal_cancer_table, quan_unbal_cancer_sig_table, quan_bal_p53_table, quan_bal_p53_sig_table, quan_bal_counts_p53_table, 
   quan_bal_counts_p53_sig_table, quan_unbal_p53_table, quan_unbal_p53_sig_table, quan_cancer_intersection_table, 
   quan_cancer_intersection_sig_table, quan_cancer_union_table, quan_cancer_union_sig_table, quan_p53_intersection_table, 
   quan_p53_intersection_sig_table, quan_p53_union_table, quan_p53_union_sig_table,
   quan_bal_counts_cancer_intersection_table, quan_bal_counts_cancer_intersection_sig_table,
   quan_bal_counts_cancer_union_table, quan_bal_counts_cancer_union_sig_table,
   quan_bal_counts_p53_intersection_table, quan_bal_counts_p53_intersection_sig_table,
   quan_bal_counts_p53_union_table, quan_bal_counts_p53_union_sig_table,
   quan_complete_cancer_intersection_table, quan_complete_cancer_intersection_sig_table,
   quan_complete_cancer_union_table, quan_complete_cancer_union_sig_table,
   quan_complete_p53_intersection_table, quan_complete_p53_intersection_sig_table,
   quan_complete_p53_union_table, quan_complete_p53_union_sig_table)


#save table 
saveRDS(quan_table, 
        file = paste0(quan_folder, '/quan_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(quan_bal_cancer_models, 
        file = paste0(quan_folder, '/quan_bal_cancer_models.rda'))
saveRDS(quan_bal_cancer_sig_models, 
        file = paste0(quan_folder, '/quan_bal_cancer_sig_models.rda'))
saveRDS(quan_bal_counts_cancer_models, 
        file = paste0(quan_folder, '/quan_bal_counts_cancer_models.rda'))
saveRDS(quan_bal_counts_cancer_sig_models, 
        file = paste0(quan_folder, '/quan_bal_counts_cancer_sig_models.rda'))
saveRDS(quan_unbal_cancer_models, 
        file = paste0(quan_folder, '/quan_unbal_cancer_models.rda'))
saveRDS(quan_unbal_cancer_sig_models, 
        file = paste0(quan_folder, '/quan_unbal_cancer_sig_models.rda'))
saveRDS(quan_bal_p53_models, 
        file = paste0(quan_folder, '/quan_bal_p53_models.rda'))
saveRDS(quan_bal_p53_sig_models, 
        file = paste0(quan_folder, '/quan_bal_p53_sig_models.rda'))
saveRDS(quan_bal_counts_p53_models, 
        file = paste0(quan_folder, '/quan_bal_counts_p53_models.rda'))
saveRDS(quan_bal_counts_p53_sig_models, 
        file = paste0(quan_folder, '/quan_bal_counts_p53_sig_models.rda'))
saveRDS(quan_unbal_p53_models, 
        file = paste0(quan_folder, '/quan_unbal_p53_models.rda'))
saveRDS(quan_unbal_p53_sig_models, 
        file = paste0(quan_folder, '/quan_unbal_p53_sig_models.rda'))
saveRDS(quan_cancer_intersection_models, 
        file = paste0(quan_folder, '/quan_cancer_intersection_models.rda'))
saveRDS(quan_cancer_intersection_sig_models, 
        file = paste0(quan_folder, '/quan_cancer_intersection_sig_models.rda'))
saveRDS(quan_cancer_union_models, 
        file = paste0(quan_folder, '/quan_cancer_union_models.rda'))
saveRDS(quan_cancer_union_sig_models, 
        file = paste0(quan_folder, '/quan_cancer_union_sig_models.rda'))
saveRDS(quan_p53_intersection_models, 
        file = paste0(quan_folder, '/quan_p53_intersection_models.rda'))
saveRDS(quan_p53_intersection_sig_models, 
        file = paste0(quan_folder, '/quan_p53_intersection_sig_models.rda'))
saveRDS(quan_p53_union_models, 
        file = paste0(quan_folder, '/quan_p53_union_models.rda'))
saveRDS(quan_p53_union_sig_models, 
        file = paste0(quan_folder, '/quan_p53_union_sig_models.rda'))
saveRDS(quan_bal_counts_cancer_intersection_models,
        file = paste0(quan_folder, '/quan_bal_counts_cancer_intersection_models.rda'))
saveRDS(quan_bal_counts_cancer_intersection_sig_models,
        file = paste0(quan_folder, '/quan_bal_counts_cancer_intersection_sig_models.rda'))
saveRDS(quan_bal_counts_cancer_union_models,
        file = paste0(quan_folder, '/quan_bal_counts_cancer_union_models.rda'))
saveRDS(quan_bal_counts_cancer_union_sig_models,
        file = paste0(quan_folder, '/quan_bal_counts_cancer_union_sig_models.rda'))
saveRDS(quan_bal_counts_p53_intersection_models,
        file = paste0(quan_folder, '/quan_bal_counts_p53_intersection_models.rda'))
saveRDS(quan_bal_counts_p53_intersection_sig_models,
        file = paste0(quan_folder, '/quan_bal_counts_p53_intersection_sig_models.rda'))
saveRDS(quan_bal_counts_p53_union_models,
        file = paste0(quan_folder, '/quan_bal_counts_p53_union_models.rda'))
saveRDS(quan_bal_counts_p53_union_sig_models,
        file = paste0(quan_folder, '/quan_bal_counts_p53_union_sig_models.rda'))
saveRDS(quan_complete_cancer_intersection_models,
        file = paste0(quan_folder, '/quan_complete_cancer_intersection_models.rda'))
# saveRDS(quan_complete_cancer_intersection_sig_models,
#         file = paste0(quan_folder, '/quan_complete_cancer_intersection_sig_models.rda'))
saveRDS(quan_complete_cancer_union_models,
        file = paste0(quan_folder, '/quan_complete_cancer_union_models.rda'))
saveRDS(quan_complete_cancer_union_sig_models,
        file = paste0(quan_folder, '/quan_complete_cancer_union_sig_models.rda'))
saveRDS(quan_complete_p53_intersection_models,
        file = paste0(quan_folder, '/quan_complete_p53_intersection_models.rda'))
saveRDS(quan_complete_p53_intersection_sig_models,
        file = paste0(quan_folder, '/quan_complete_p53_intersection_sig_models.rda'))
saveRDS(quan_complete_p53_union_models,
        file = paste0(quan_folder, '/quan_complete_p53_union_models.rda'))
saveRDS(quan_complete_p53_union_sig_models,
        file = paste0(quan_folder, '/quan_complete_p53_union_sig_models.rda'))


###################################################################################################################################
## beta_funnorm

##########
# cancer
##########

# bal cancer
funnorm_bal_cancer_models <- runModels(beta_funnorm, 
                                   random = F, 
                                   bump_hunter = T, 
                                   bump_hunter_data = beta_funnorm_bal_cancer_features)

funnorm_bal_cancer_table <- extractResults(funnorm_bal_cancer_models, 
                                       data_name = 'beta_funnorm_bal_cancer')

# bal cancer sig
funnorm_bal_cancer_sig_models <- runModels(beta_funnorm, 
                                       random = F, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_funnorm_bal_cancer_sig_features)

funnorm_bal_cancer_sig_table <- extractResults(funnorm_bal_cancer_sig_models, 
                                           data_name = 'funnorm_bal_cancer_sig')


# bal counts cancer
funnorm_bal_counts_cancer_models <- runModels(beta_funnorm, 
                                          random = F, 
                                          bump_hunter = T, 
                                          bump_hunter_data = beta_funnorm_bal_counts_cancer_features)

funnorm_bal_counts_cancer_table <- extractResults(funnorm_bal_counts_cancer_models, 
                                              data_name = 'funnorm_bal_counts_cancer')

# bal counts cancer sig
funnorm_bal_counts_cancer_sig_models <- runModels(beta_funnorm, 
                                              random = F, 
                                              bump_hunter = T, 
                                              bump_hunter_data =  beta_funnorm_bal_counts_cancer_sig_features)

funnorm_bal_counts_cancer_sig_table <- extractResults(funnorm_bal_counts_cancer_sig_models, 
                                                  data_name = 'funnorm_bal_counts_cancer_sig')

# unbal cancer
funnorm_unbal_cancer_models <- runModels(beta_funnorm, 
                                     random = F, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_funnorm_unbal_cancer_features)

funnorm_unbal_cancer_table <- extractResults(funnorm_unbal_cancer_models, 
                                         data_name = 'funnorm_unbal_cancer')


# unbal cancer sig
funnorm_unbal_cancer_sig_models <- runModels(beta_funnorm, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_funnorm_unbal_cancer_sig_features)

funnorm_unbal_cancer_sig_table <- extractResults(funnorm_unbal_cancer_sig_models, 
                                             data_name = 'funnorm_unbal_cancer_sig')

##########
# p53
##########

# bal p53
funnorm_bal_p53_models <- runModels(beta_funnorm, 
                                random = F, 
                                bump_hunter = T, 
                                bump_hunter_data = beta_funnorm_bal_p53_features)

funnorm_bal_p53_table <- extractResults(funnorm_bal_p53_models, 
                                    data_name = 'funnorm_bal_p53')

# bal p53 sig
funnorm_bal_p53_sig_models <- runModels(beta_funnorm, 
                                    random = F, 
                                    bump_hunter = T, 
                                    bump_hunter_data = beta_funnorm_bal_p53_sig_features)

funnorm_bal_p53_sig_table <- extractResults(funnorm_bal_p53_sig_models, 
                                        data_name = 'funnorm_bal_p53_sig')

# bal counts p53
funnorm_bal_counts_p53_models <- runModels(beta_funnorm, 
                                       random = F, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_funnorm_bal_counts_p53_features)

funnorm_bal_counts_p53_table <- extractResults(funnorm_bal_counts_p53_models, 
                                           data_name = 'funnorm_bal_counts_p53')

# bal counts p53 sig
funnorm_bal_counts_p53_sig_models <- runModels(beta_funnorm, 
                                           random = F, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_funnorm_bal_counts_p53_sig_features)

funnorm_bal_counts_p53_sig_table <- extractResults(funnorm_bal_counts_p53_sig_models, 
                                               data_name = 'funnorm_bal_counts_p53_sig')

# unbal p53
funnorm_unbal_p53_models <- runModels(beta_funnorm, 
                                  random = F, 
                                  bump_hunter = T, 
                                  bump_hunter_data = beta_funnorm_unbal_p53_features)

funnorm_unbal_p53_table <- extractResults(funnorm_unbal_p53_models, 
                                      data_name = 'funnorm_unbal_p53')


# unbal p53 sig
funnorm_unbal_p53_sig_models <- runModels(beta_funnorm, 
                                      random = F, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_funnorm_unbal_p53_sig_features)

funnorm_unbal_p53_sig_table <- extractResults(funnorm_unbal_p53_sig_models, 
                                          data_name = 'funnorm_unbal_p53_sig')


##########
# cancer
##########

# cancer_intersection
funnorm_cancer_intersection_models <- runModels(beta_funnorm, 
                                            random = F, 
                                            bump_hunter = T, 
                                            bump_hunter_data = beta_funnorm_cancer_intersection_features)

funnorm_cancer_intersection_table <- extractResults(funnorm_cancer_intersection_models, 
                                                data_name = 'funnorm_cancer_intersection')

# cancer_intersection sig
funnorm_cancer_intersection_sig_models <- runModels(beta_funnorm, 
                                                random = F, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_funnorm_cancer_intersection_sig_features)

funnorm_cancer_intersection_sig_table <- extractResults(funnorm_cancer_intersection_sig_models, 
                                                    data_name = 'funnorm_cancer_intersection_sig')

# cancer_union
funnorm_cancer_union_models <- runModels(beta_funnorm, 
                                     random = F, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_funnorm_cancer_union_features)

funnorm_cancer_union_table <- extractResults(funnorm_cancer_union_models, 
                                         data_name = 'funnorm_cancer_union')

# cancer_union sig
funnorm_cancer_union_sig_models <- runModels(beta_funnorm, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_funnorm_cancer_union_sig_features)

funnorm_cancer_union_sig_table <- extractResults(funnorm_cancer_union_sig_models, 
                                             data_name = 'funnorm_cancer_union_sig')
##########
# p53
##########
# p53_intersection
funnorm_p53_intersection_models <- runModels(beta_funnorm, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_funnorm_p53_intersection_features)

funnorm_p53_intersection_table <- extractResults(funnorm_p53_intersection_models, 
                                             data_name = 'funnorm_p53_intersection')

# p53_intersection sig
funnorm_p53_intersection_sig_models <- runModels(beta_funnorm, 
                                             random = F, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_funnorm_p53_intersection_sig_features)

funnorm_p53_intersection_sig_table <- extractResults(funnorm_p53_intersection_sig_models, 
                                                 data_name = 'funnorm_p53_intersection_sig')

# p53_union
funnorm_p53_union_models <- runModels(beta_funnorm, 
                                  random = F, 
                                  bump_hunter = T, 
                                  bump_hunter_data = beta_funnorm_p53_union_features)

funnorm_p53_union_table <- extractResults(funnorm_p53_union_models, 
                                      data_name = 'funnorm_p53_union')

# p53_union sig
funnorm_p53_union_sig_models <- runModels(beta_funnorm, 
                                      random = F, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_funnorm_p53_union_sig_features)

funnorm_p53_union_sig_table <- extractResults(funnorm_p53_union_sig_models, 
                                          data_name = 'funnorm_p53_union_sig')


##########
# most balanced cancer 
##########

# intersection
funnorm_bal_counts_cancer_intersection_models <- runModels(beta_funnorm, 
                                                       random = F, 
                                                       bump_hunter = T, 
                                                       bump_hunter_data = beta_bal_counts_cancer_intersection_features)

funnorm_bal_counts_cancer_intersection_table <- extractResults(funnorm_bal_counts_cancer_intersection_models, 
                                                           data_name = 'funnorm_bal_counts_cancer_intersection')


funnorm_bal_counts_cancer_intersection_sig_models <- runModels(beta_funnorm, 
                                                           random = F, 
                                                           bump_hunter = T, 
                                                           bump_hunter_data = beta_bal_counts_cancer_intersection_sig_features)

funnorm_bal_counts_cancer_intersection_sig_table <- extractResults(funnorm_bal_counts_cancer_intersection_sig_models, 
                                                               data_name = 'funnorm_bal_counts_cancer_intersection_sig')

# union
funnorm_bal_counts_cancer_union_models <- runModels(beta_funnorm, 
                                                random = F, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_bal_counts_cancer_union_features)

funnorm_bal_counts_cancer_union_table <- extractResults(funnorm_bal_counts_cancer_union_models, 
                                                    data_name = 'funnorm_bal_counts_cancer_union')


funnorm_bal_counts_cancer_union_sig_models <- runModels(beta_funnorm, 
                                                    random = F, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_cancer_union_sig_features)

funnorm_bal_counts_cancer_union_sig_table <- extractResults(funnorm_bal_counts_cancer_union_sig_models, 
                                                        data_name = 'funnorm_bal_counts_cancer_union_sig')

##########
# most balanced p53 
##########

# intersection
funnorm_bal_counts_p53_intersection_models <- runModels(beta_funnorm, 
                                                    random = F, 
                                                    bump_hunter = T, 
                                                    bump_hunter_data = beta_bal_counts_p53_intersection_features)

funnorm_bal_counts_p53_intersection_table <- extractResults(funnorm_bal_counts_p53_intersection_models, 
                                                        data_name = 'funnorm_bal_counts_p53_intersection')


funnorm_bal_counts_p53_intersection_sig_models <- runModels(beta_funnorm, 
                                                        random = F, 
                                                        bump_hunter = T, 
                                                        bump_hunter_data = beta_bal_counts_p53_intersection_sig_features)

funnorm_bal_counts_p53_intersection_sig_table <- extractResults(funnorm_bal_counts_p53_intersection_sig_models, 
                                                            data_name = 'funnorm_bal_counts_p53_intersection_sig')


# union
funnorm_bal_counts_p53_union_models <- runModels(beta_funnorm, 
                                             random = F, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_bal_counts_p53_union_features)

funnorm_bal_counts_p53_union_table <- extractResults(funnorm_bal_counts_p53_union_models, 
                                                 data_name = 'funnorm_bal_counts_p53_union')


funnorm_bal_counts_p53_union_sig_models <- runModels(beta_funnorm, 
                                                 random = F, 
                                                 bump_hunter = T, 
                                                 bump_hunter_data = beta_bal_counts_p53_union_sig_features)

funnorm_bal_counts_p53_union_sig_table <- extractResults(funnorm_bal_counts_p53_union_sig_models, 
                                                     data_name = 'funnorm_bal_counts_p53_union_sig')

###########
# complete cancer - intersection across each method intersection
###########

# complete cancer intersection
funnorm_complete_cancer_intersection_models <- runModels(beta_funnorm, 
                                                     random = F, 
                                                     bump_hunter = T, 
                                                     bump_hunter_data = beta_cancer_intersection_features)

# get table 
funnorm_complete_cancer_intersection_table <- extractResults(funnorm_complete_cancer_intersection_models, 
                                                         data_name = 'funnorm_complete_cancer_intersection')

# empty features set
# # complete cancer intersection sig
# funnorm_complete_cancer_intersection_sig_models <- runModels(beta_funnorm, 
#                                                          random = F, 
#                                                          bump_hunter = T, 
#                                                          bump_hunter_data = beta_cancer_intersection_sig_features)
# 
# # get table 
# funnorm_complete_cancer_intersection_sig_table <- extractResults(funnorm_complete_cancer_intersection_sig_models, 
#                                                              data_name = 'funnorm_complete_cancer_intersection_sig')


# complete cancer union
funnorm_complete_cancer_union_models <- runModels(beta_funnorm, 
                                              random = F, 
                                              bump_hunter = T, 
                                              bump_hunter_data = beta_cancer_union_features)

# get table 
funnorm_complete_cancer_union_table <- extractResults(funnorm_complete_cancer_union_models, 
                                                  data_name = 'funnorm_complete_cancer_union')

# complete cancer union sig
funnorm_complete_cancer_union_sig_models <- runModels(beta_funnorm, 
                                                  random = F, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_cancer_union_sig_features)

# get table 
funnorm_complete_cancer_union_sig_table <- extractResults(funnorm_complete_cancer_union_sig_models, 
                                                      data_name = 'funnorm_complete_cancer_union_sig')


###########
# complete p53 - intersection across each method intersection
###########

# complete p53 intersection
funnorm_complete_p53_intersection_models <- runModels(beta_funnorm, 
                                                  random = F, 
                                                  bump_hunter = T, 
                                                  bump_hunter_data = beta_p53_intersection_features)

# get table 
funnorm_complete_p53_intersection_table <- extractResults(funnorm_complete_p53_intersection_models, 
                                                      data_name = 'funnorm_complete_p53_intersection')

# complete p53 intersection sig
funnorm_complete_p53_intersection_sig_models <- runModels(beta_funnorm, 
                                                      random = F, 
                                                      bump_hunter = T, 
                                                      bump_hunter_data = beta_p53_intersection_sig_features)

# get table 
funnorm_complete_p53_intersection_sig_table <- extractResults(funnorm_complete_p53_intersection_sig_models, 
                                                          data_name = 'funnorm_complete_p53_intersection_sig')


# complete p53 union
funnorm_complete_p53_union_models <- runModels(beta_funnorm, 
                                           random = F, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_p53_union_features)

# get table 
funnorm_complete_p53_union_table <- extractResults(funnorm_complete_p53_union_models, 
                                               data_name = 'funnorm_complete_p53_union')

# complete p53 union sig
funnorm_complete_p53_union_sig_models <- runModels(beta_funnorm, 
                                               random = F, 
                                               bump_hunter = T, 
                                               bump_hunter_data = beta_p53_union_sig_features)

# get table 
funnorm_complete_p53_union_sig_table <- extractResults(funnorm_complete_p53_union_sig_models, 
                                                   data_name = 'funnorm_complete_p53_union_sig')


###########
# rbind tables and save RDA file
###########
funnorm_table <- rbind(funnorm_bal_cancer_table, funnorm_bal_cancer_sig_table, funnorm_bal_counts_cancer_table, funnorm_bal_counts_cancer_sig_table,
                   funnorm_unbal_cancer_table, funnorm_unbal_cancer_sig_table, funnorm_bal_p53_table, funnorm_bal_p53_sig_table, funnorm_bal_counts_p53_table, 
                   funnorm_bal_counts_p53_sig_table, funnorm_unbal_p53_table, funnorm_unbal_p53_sig_table, funnorm_cancer_intersection_table, 
                   funnorm_cancer_intersection_sig_table, funnorm_cancer_union_table, funnorm_cancer_union_sig_table, funnorm_p53_intersection_table, 
                   funnorm_p53_intersection_sig_table, funnorm_p53_union_table, funnorm_p53_union_sig_table,
                   funnorm_bal_counts_cancer_intersection_table, funnorm_bal_counts_cancer_intersection_sig_table,
                   funnorm_bal_counts_cancer_union_table, funnorm_bal_counts_cancer_union_sig_table,
                   funnorm_bal_counts_p53_intersection_table, funnorm_bal_counts_p53_intersection_sig_table,
                   funnorm_bal_counts_p53_union_table, funnorm_bal_counts_p53_union_sig_table,
                   funnorm_complete_cancer_intersection_table, funnorm_complete_cancer_intersection_sig_table,
                   funnorm_complete_cancer_union_table, funnorm_complete_cancer_union_sig_table,
                   funnorm_complete_p53_intersection_table, funnorm_complete_p53_intersection_sig_table,
                   funnorm_complete_p53_union_table, funnorm_complete_p53_union_sig_table)

# remove data 
rm(funnorm_bal_cancer_table, funnorm_bal_cancer_sig_table, funnorm_bal_counts_cancer_table, funnorm_bal_counts_cancer_sig_table,
   funnorm_unbal_cancer_table, funnorm_unbal_cancer_sig_table, funnorm_bal_p53_table, funnorm_bal_p53_sig_table, funnorm_bal_counts_p53_table, 
   funnorm_bal_counts_p53_sig_table, funnorm_unbal_p53_table, funnorm_unbal_p53_sig_table, funnorm_cancer_intersection_table, 
   funnorm_cancer_intersection_sig_table, funnorm_cancer_union_table, funnorm_cancer_union_sig_table, funnorm_p53_intersection_table, 
   funnorm_p53_intersection_sig_table, funnorm_p53_union_table, funnorm_p53_union_sig_table,
   funnorm_bal_counts_cancer_intersection_table, funnorm_bal_counts_cancer_intersection_sig_table,
   funnorm_bal_counts_cancer_union_table, funnorm_bal_counts_cancer_union_sig_table,
   funnorm_bal_counts_p53_intersection_table, funnorm_bal_counts_p53_intersection_sig_table,
   funnorm_bal_counts_p53_union_table, funnorm_bal_counts_p53_union_sig_table,
   funnorm_complete_cancer_intersection_table, funnorm_complete_cancer_intersection_sig_table,
   funnorm_complete_cancer_union_table, funnorm_complete_cancer_union_sig_table,
   funnorm_complete_p53_intersection_table, funnorm_complete_p53_intersection_sig_table,
   funnorm_complete_p53_union_table, funnorm_complete_p53_union_sig_table)


#save table 
saveRDS(funnorm_table, 
        file = paste0(funnorm_folder, '/funnorm_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(funnorm_bal_cancer_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_cancer_models.rda'))
saveRDS(funnorm_bal_cancer_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_cancer_sig_models.rda'))
saveRDS(funnorm_bal_counts_cancer_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_counts_cancer_models.rda'))
saveRDS(funnorm_bal_counts_cancer_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_counts_cancer_sig_models.rda'))
saveRDS(funnorm_unbal_cancer_models, 
        file = paste0(funnorm_folder, '/funnorm_unbal_cancer_models.rda'))
saveRDS(funnorm_unbal_cancer_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_unbal_cancer_sig_models.rda'))
saveRDS(funnorm_bal_p53_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_p53_models.rda'))
saveRDS(funnorm_bal_p53_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_p53_sig_models.rda'))
saveRDS(funnorm_bal_counts_p53_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_counts_p53_models.rda'))
saveRDS(funnorm_bal_counts_p53_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_bal_counts_p53_sig_models.rda'))
saveRDS(funnorm_unbal_p53_models, 
        file = paste0(funnorm_folder, '/funnorm_unbal_p53_models.rda'))
saveRDS(funnorm_unbal_p53_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_unbal_p53_sig_models.rda'))
saveRDS(funnorm_cancer_intersection_models, 
        file = paste0(funnorm_folder, '/funnorm_cancer_intersection_models.rda'))
saveRDS(funnorm_cancer_intersection_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_cancer_intersection_sig_models.rda'))
saveRDS(funnorm_cancer_union_models, 
        file = paste0(funnorm_folder, '/funnorm_cancer_union_models.rda'))
saveRDS(funnorm_cancer_union_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_cancer_union_sig_models.rda'))
saveRDS(funnorm_p53_intersection_models, 
        file = paste0(funnorm_folder, '/funnorm_p53_intersection_models.rda'))
saveRDS(funnorm_p53_intersection_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_p53_intersection_sig_models.rda'))
saveRDS(funnorm_p53_union_models, 
        file = paste0(funnorm_folder, '/funnorm_p53_union_models.rda'))
saveRDS(funnorm_p53_union_sig_models, 
        file = paste0(funnorm_folder, '/funnorm_p53_union_sig_models.rda'))
saveRDS(funnorm_bal_counts_cancer_intersection_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_cancer_intersection_models.rda'))
saveRDS(funnorm_bal_counts_cancer_intersection_sig_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_cancer_intersection_sig_models.rda'))
saveRDS(funnorm_bal_counts_cancer_union_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_cancer_union_models.rda'))
saveRDS(funnorm_bal_counts_cancer_union_sig_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_cancer_union_sig_models.rda'))
saveRDS(funnorm_bal_counts_p53_intersection_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_p53_intersection_models.rda'))
saveRDS(funnorm_bal_counts_p53_intersection_sig_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_p53_intersection_sig_models.rda'))
saveRDS(funnorm_bal_counts_p53_union_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_p53_union_models.rda'))
saveRDS(funnorm_bal_counts_p53_union_sig_models,
        file = paste0(funnorm_folder, '/funnorm_bal_counts_p53_union_sig_models.rda'))
saveRDS(funnorm_complete_cancer_intersection_models,
        file = paste0(funnorm_folder, '/funnorm_complete_cancer_intersection_models.rda'))
# saveRDS(funnorm_complete_cancer_intersection_sig_models,
#         file = paste0(funnorm_folder, '/funnorm_complete_cancer_intersection_sig_models.rda'))
saveRDS(funnorm_complete_cancer_union_models,
        file = paste0(funnorm_folder, '/funnorm_complete_cancer_union_models.rda'))
saveRDS(funnorm_complete_cancer_union_sig_models,
        file = paste0(funnorm_folder, '/funnorm_complete_cancer_union_sig_models.rda'))
saveRDS(funnorm_complete_p53_intersection_models,
        file = paste0(funnorm_folder, '/funnorm_complete_p53_intersection_models.rda'))
saveRDS(funnorm_complete_p53_intersection_sig_models,
        file = paste0(funnorm_folder, '/funnorm_complete_p53_intersection_sig_models.rda'))
saveRDS(funnorm_complete_p53_union_models,
        file = paste0(funnorm_folder, '/funnorm_complete_p53_union_models.rda'))
saveRDS(funnorm_complete_p53_union_sig_models,
        file = paste0(funnorm_folder, '/funnorm_complete_p53_union_sig_models.rda'))



#########################################################################################################################
# Random features - 100, 200, 500, 1000, 2000, 10000

###########
# raw
###########

# raw_100
raw_rand_100_models <- runModels(beta_raw, 
                          random = T,
                          num_features = 100,
                          bump_hunter = F)

# get table
raw_rand_100_table <- extractResults(raw_rand_100_models, 
                                     data_name = 'raw_rand_100')


# raw_200
raw_rand_200_models <- runModels(beta_raw, 
                                 random = T, 
                                 num_features = 200,
                                 bump_hunter = F)

# get table
raw_rand_200_table <- extractResults(raw_rand_200_models, 
                                     data_name = 'raw_rand_200')

# raw_500
raw_rand_500_models <- runModels(beta_raw, 
                                 random = T, 
                                 num_features = 500,
                                 bump_hunter = F)

# get table
raw_rand_500_table <- extractResults(raw_rand_500_models, 
                                     data_name = 'raw_rand_500')


# raw_1000
raw_rand_1000_models <- runModels(beta_raw, 
                                 random = T, 
                                 num_features = 1000,
                                 bump_hunter = F)

# get table
raw_rand_1000_table <- extractResults(raw_rand_1000_models, 
                                     data_name = 'raw_rand_1000')

# raw_2000
raw_rand_2000_models <- runModels(beta_raw, 
                                 random = T, 
                                 num_features = 2000,
                                 bump_hunter = F)

# get table
raw_rand_2000_table <- extractResults(raw_rand_2000_models, 
                                     data_name = 'raw_rand_2000')


# raw_10000
raw_rand_10000_models <- runModels(beta_raw, 
                                  random = T, 
                                  num_features = 10000,
                                  bump_hunter = F)

# get table
raw_rand_10000_table <- extractResults(raw_rand_10000_models, 
                                      data_name = 'raw_rand_10000')


###########
# rbind tables and save RDA file
###########
rand_table <- rbind(raw_rand_100_table, raw_rand_200_table, 
                    raw_rand_500_table, raw_rand_1000_table,
                    raw_rand_2000_table, raw_rand_10000_table)

# remove data 
rm(raw_rand_100_table, raw_rand_200_table, 
   raw_rand_500_table, raw_rand_1000_table,
   raw_rand_2000_table, raw_rand_10000_table)


#save table 
saveRDS(rand_table, 
        file = paste0(rand_folder, '/rand_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(raw_rand_100_models, 
        file = paste0(rand_folder, '/raw_rand_100_models.rda'))
saveRDS(raw_rand_200_models, 
        file = paste0(rand_folder, '/raw_rand_200_models.rda'))
saveRDS(raw_rand_500_models, 
        file = paste0(rand_folder, '/raw_rand_500_models.rda'))
saveRDS(raw_rand_1000_models, 
        file = paste0(rand_folder, '/raw_rand_1000_models.rda'))
saveRDS(raw_rand_2000_models, 
        file = paste0(rand_folder, '/raw_rand_2000_models.rda'))
saveRDS(raw_rand_10000_models, 
        file = paste0(rand_folder, '/raw_rand_10000_models.rda'))


###########
# swan
###########

# swan_100
swan_rand_100_models <- runModels(beta_swan, 
                                 random = T,
                                 num_features = 100,
                                 bump_hunter = F)

# get table
swan_rand_100_table <- extractResults(swan_rand_100_models, 
                                     data_name = 'swan_rand_100')


# swan_200
swan_rand_200_models <- runModels(beta_swan, 
                                 random = T, 
                                 num_features = 200,
                                 bump_hunter = F)

# get table
swan_rand_200_table <- extractResults(swan_rand_200_models, 
                                     data_name = 'swan_rand_200')

# swan_500
swan_rand_500_models <- runModels(beta_swan, 
                                 random = T, 
                                 num_features = 500,
                                 bump_hunter = F)

# get table
swan_rand_500_table <- extractResults(swan_rand_500_models, 
                                     data_name = 'swan_rand_500')


# swan_1000
swan_rand_1000_models <- runModels(beta_swan, 
                                  random = T, 
                                  num_features = 1000,
                                  bump_hunter = F)

# get table
swan_rand_1000_table <- extractResults(swan_rand_1000_models, 
                                      data_name = 'swan_rand_1000')

# swan_2000
swan_rand_2000_models <- runModels(beta_swan, 
                                  random = T, 
                                  num_features = 2000,
                                  bump_hunter = F)

# get table
swan_rand_2000_table <- extractResults(swan_rand_2000_models, 
                                      data_name = 'swan_rand_2000')


# swan_10000
swan_rand_10000_models <- runModels(beta_swan, 
                                   random = T, 
                                   num_features = 10000,
                                   bump_hunter = F)

# get table
swan_rand_10000_table <- extractResults(swan_rand_10000_models, 
                                       data_name = 'swan_rand_10000')


###########
# rbind tables and save RDA file
###########
rand_table <- rbind(swan_rand_100_table, swan_rand_200_table, 
                    swan_rand_500_table, swan_rand_1000_table,
                    swan_rand_2000_table, swan_rand_10000_table)

# remove data 
rm(swan_rand_100_table, swan_rand_200_table, 
   swan_rand_500_table, swan_rand_1000_table,
   swan_rand_2000_table, swan_rand_10000_table)


#save table 
saveRDS(rand_table, 
        file = paste0(rand_folder, '/rand_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(swan_rand_100_models, 
        file = paste0(rand_folder, '/swan_rand_100_models.rda'))
saveRDS(swan_rand_200_models, 
        file = paste0(rand_folder, '/swan_rand_200_models.rda'))
saveRDS(swan_rand_500_models, 
        file = paste0(rand_folder, '/swan_rand_500_models.rda'))
saveRDS(swan_rand_1000_models, 
        file = paste0(rand_folder, '/swan_rand_1000_models.rda'))
saveRDS(swan_rand_2000_models, 
        file = paste0(rand_folder, '/swan_rand_2000_models.rda'))
saveRDS(swan_rand_10000_models, 
        file = paste0(rand_folder, '/swan_rand_10000_models.rda'))


###########
# quan
###########

# quan_100
quan_rand_100_models <- runModels(beta_quan, 
                                 random = T,
                                 num_features = 100,
                                 bump_hunter = F)

# get table
quan_rand_100_table <- extractResults(quan_rand_100_models, 
                                     data_name = 'quan_rand_100')


# quan_200
quan_rand_200_models <- runModels(beta_quan, 
                                 random = T, 
                                 num_features = 200,
                                 bump_hunter = F)

# get table
quan_rand_200_table <- extractResults(quan_rand_200_models, 
                                     data_name = 'quan_rand_200')

# quan_500
quan_rand_500_models <- runModels(beta_quan, 
                                 random = T, 
                                 num_features = 500,
                                 bump_hunter = F)

# get table
quan_rand_500_table <- extractResults(quan_rand_500_models, 
                                     data_name = 'quan_rand_500')


# quan_1000
quan_rand_1000_models <- runModels(beta_quan, 
                                  random = T, 
                                  num_features = 1000,
                                  bump_hunter = F)

# get table
quan_rand_1000_table <- extractResults(quan_rand_1000_models, 
                                      data_name = 'quan_rand_1000')

# quan_2000
quan_rand_2000_models <- runModels(beta_quan, 
                                  random = T, 
                                  num_features = 2000,
                                  bump_hunter = F)

# get table
quan_rand_2000_table <- extractResults(quan_rand_2000_models, 
                                      data_name = 'quan_rand_2000')


# quan_10000
quan_rand_10000_models <- runModels(beta_quan, 
                                   random = T, 
                                   num_features = 10000,
                                   bump_hunter = F)

# get table
quan_rand_10000_table <- extractResults(quan_rand_10000_models, 
                                       data_name = 'quan_rand_10000')


###########
# rbind tables and save RDA file
###########
rand_table <- rbind(quan_rand_100_table, quan_rand_200_table, 
                    quan_rand_500_table, quan_rand_1000_table,
                    quan_rand_2000_table, quan_rand_10000_table)

# remove data 
rm(quan_rand_100_table, quan_rand_200_table, 
   quan_rand_500_table, quan_rand_1000_table,
   quan_rand_2000_table, quan_rand_10000_table)


#save table 
saveRDS(rand_table, 
        file = paste0(rand_folder, '/rand_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(quan_rand_100_models, 
        file = paste0(rand_folder, '/quan_rand_100_models.rda'))
saveRDS(quan_rand_200_models, 
        file = paste0(rand_folder, '/quan_rand_200_models.rda'))
saveRDS(quan_rand_500_models, 
        file = paste0(rand_folder, '/quan_rand_500_models.rda'))
saveRDS(quan_rand_1000_models, 
        file = paste0(rand_folder, '/quan_rand_1000_models.rda'))
saveRDS(quan_rand_2000_models, 
        file = paste0(rand_folder, '/quan_rand_2000_models.rda'))
saveRDS(quan_rand_10000_models, 
        file = paste0(rand_folder, '/quan_rand_10000_models.rda'))


###########
# funnorm
###########

# funnorm_100
funnorm_rand_100_models <- runModels(beta_funnorm, 
                                 random = T,
                                 num_features = 100,
                                 bump_hunter = F)

# get table
funnorm_rand_100_table <- extractResults(funnorm_rand_100_models, 
                                     data_name = 'funnorm_rand_100')


# funnorm_200
funnorm_rand_200_models <- runModels(beta_funnorm, 
                                 random = T, 
                                 num_features = 200,
                                 bump_hunter = F)

# get table
funnorm_rand_200_table <- extractResults(funnorm_rand_200_models, 
                                     data_name = 'funnorm_rand_200')

# funnorm_500
funnorm_rand_500_models <- runModels(beta_funnorm, 
                                 random = T, 
                                 num_features = 500,
                                 bump_hunter = F)

# get table
funnorm_rand_500_table <- extractResults(funnorm_rand_500_models, 
                                     data_name = 'funnorm_rand_500')


# funnorm_1000
funnorm_rand_1000_models <- runModels(beta_funnorm, 
                                  random = T, 
                                  num_features = 1000,
                                  bump_hunter = F)

# get table
funnorm_rand_1000_table <- extractResults(funnorm_rand_1000_models, 
                                      data_name = 'funnorm_rand_1000')

# funnorm_2000
funnorm_rand_2000_models <- runModels(beta_funnorm, 
                                  random = T, 
                                  num_features = 2000,
                                  bump_hunter = F)

# get table
funnorm_rand_2000_table <- extractResults(funnorm_rand_2000_models, 
                                      data_name = 'funnorm_rand_2000')


# funnorm_10000
funnorm_rand_10000_models <- runModels(beta_funnorm, 
                                   random = T, 
                                   num_features = 10000,
                                   bump_hunter = F)

# get table
funnorm_rand_10000_table <- extractResults(funnorm_rand_10000_models, 
                                       data_name = 'funnorm_rand_10000')


###########
# rbind tables and save RDA file
###########
rand_table <- rbind(funnorm_rand_100_table, funnorm_rand_200_table, 
                    funnorm_rand_500_table, funnorm_rand_1000_table,
                    funnorm_rand_2000_table, funnorm_rand_10000_table)

# remove data 
rm(funnorm_rand_100_table, funnorm_rand_200_table, 
   funnorm_rand_500_table, funnorm_rand_1000_table,
   funnorm_rand_2000_table, funnorm_rand_10000_table)


#save table 
saveRDS(rand_table, 
        file = paste0(rand_folder, '/rand_table.rda'))


###########
# save all model objects as RDA file
###########
saveRDS(funnorm_rand_100_models, 
        file = paste0(rand_folder, '/funnorm_rand_100_models.rda'))
saveRDS(funnorm_rand_200_models, 
        file = paste0(rand_folder, '/funnorm_rand_200_models.rda'))
saveRDS(funnorm_rand_500_models, 
        file = paste0(rand_folder, '/funnorm_rand_500_models.rda'))
saveRDS(funnorm_rand_1000_models, 
        file = paste0(rand_folder, '/funnorm_rand_1000_models.rda'))
saveRDS(funnorm_rand_2000_models, 
        file = paste0(rand_folder, '/funnorm_rand_2000_models.rda'))
saveRDS(funnorm_rand_10000_models, 
        file = paste0(rand_folder, '/funnorm_rand_10000_models.rda'))


