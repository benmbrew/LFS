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


###########
# rbind tables and save RDA file
###########
raw_table <- rbind(raw_bal_cancer_table, raw_bal_cancer_sig_table, raw_bal_counts_cancer_table, raw_bal_counts_cancer_sig_table,
                   raw_unbal_cancer_table, raw_unbal_cancer_sig_table, raw_bal_p53_table, raw_bal_p53_sig_table, raw_bal_counts_p53_table, 
                   raw_bal_counts_p53_sig_table, raw_unbal_p53_table, raw_unbal_p53_sig_table, raw_cancer_intersection_table, 
                   raw_cancer_intersection_sig_table, raw_cancer_union_table, raw_cancer_union_sig_table, raw_p53_intersection_table, 
                   raw_p53_intersection_sig_table, raw_p53_union_table, raw_p53_union_sig_table)

# remove data 
rm(raw_bal_cancer_table, raw_bal_cancer_sig_table, raw_bal_counts_cancer_table, raw_bal_counts_cancer_sig_table,
   raw_unbal_cancer_table, raw_unbal_cancer_sig_table, raw_bal_p53_table, raw_bal_p53_sig_table, raw_bal_counts_p53_table, 
   raw_bal_counts_p53_sig_table, raw_unbal_p53_table, raw_unbal_p53_sig_table, raw_cancer_intersection_table, 
   raw_cancer_intersection_sig_table, raw_cancer_union_table, raw_cancer_union_sig_table, raw_p53_intersection_table, 
   raw_p53_intersection_sig_table, raw_p53_union_table, raw_p53_union_sig_table)

# remove table
saveRDS(raw_table, 
        file = paste0(results_folder, '/raw_table.rda'))

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


# remove model objects
rm(raw_bal_cancer_models, raw_bal_cancer_sig_models, raw_bal_counts_cancer_models, raw_bal_counts_cancer_sig_models,
   raw_unbal_cancer_models, raw_unbal_cancer_sig_models, raw_bal_p53_models, raw_bal_p53_sig_models, raw_bal_counts_p53_models, 
   raw_bal_counts_p53_sig_models, raw_unbal_p53_models, raw_unbal_p53_sig_models, raw_cancer_intersection_models, 
   raw_cancer_intersection_sig_models, raw_cancer_union_models, raw_cancer_union_sig_models, raw_p53_intersection_models, 
   raw_p53_intersection_sig_models, raw_p53_union_models, raw_p53_union_sig_models)

###################################################################################################################################
## beta_swan #HERE
save.image('/home/benbrew/Desktop/temp_Results.RData')
load('/home/benbrew/Desktop/temp_Results.RData')

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



save.image('/home/benbrew/Desktop/temp_Results.RData')
load('/home/benbrew/Desktop/temp_Results.RData')

###########
# rbind tables and save RDA file
###########
swan_table <- rbind(swan_bal_cancer_table, swan_bal_cancer_sig_table, swan_bal_counts_cancer_table, swan_bal_counts_cancer_sig_table,
                   swan_unbal_cancer_table, swan_unbal_cancer_sig_table, swan_bal_p53_table, swan_bal_p53_sig_table, swan_bal_counts_p53_table, 
                   swan_bal_counts_p53_sig_table, swan_unbal_p53_table, swan_unbal_p53_sig_table, swan_cancer_intersection_table, 
                   swan_cancer_intersection_sig_table, swan_cancer_union_table, swan_cancer_union_sig_table, swan_p53_intersection_table, 
                   swan_p53_intersection_sig_table, swan_p53_union_table, swan_p53_union_sig_table)

# remove data 
rm(swan_bal_cancer_table, swan_bal_cancer_sig_table, swan_bal_counts_cancer_table, swan_bal_counts_cancer_sig_table,
   swan_unbal_cancer_table, swan_unbal_cancer_sig_table, swan_bal_p53_table, swan_bal_p53_sig_table, swan_bal_counts_p53_table, 
   swan_bal_counts_p53_sig_table, swan_unbal_p53_table, swan_unbal_p53_sig_table, swan_cancer_intersection_table, 
   swan_cancer_intersection_sig_table, swan_cancer_union_table, swan_cancer_union_sig_table, swan_p53_intersection_table, 
   swan_p53_intersection_sig_table, swan_p53_union_table, swan_p53_union_sig_table)

# remove table
saveRDS(swan_table, 
        file = paste0(results_folder, '/swan_table.rda'))

###########
# save all model objects as RDA file
###########
saveRDS(swan_bal_cancer_models, 
        file = paste0(results_folder, '/swan_bal_cancer_models.rda'))
saveRDS(swan_bal_cancer_sig_models, 
        file = paste0(results_folder, '/swan_bal_cancer_sig_models.rda'))
saveRDS(swan_bal_counts_cancer_models, 
        file = paste0(results_folder, '/swan_bal_counts_cancer_models.rda'))
saveRDS(swan_bal_counts_cancer_sig_models, 
        file = paste0(results_folder, '/swan_bal_counts_cancer_sig_models.rda'))
saveRDS(swan_unbal_cancer_models, 
        file = paste0(results_folder, '/swan_unbal_cancer_models.rda'))
saveRDS(swan_unbal_cancer_sig_models, 
        file = paste0(results_folder, '/swan_unbal_cancer_sig_models.rda'))
saveRDS(swan_bal_p53_models, 
        file = paste0(results_folder, '/swan_bal_p53_models.rda'))
saveRDS(swan_bal_p53_sig_models, 
        file = paste0(results_folder, '/swan_bal_p53_sig_models.rda'))
saveRDS(swan_bal_counts_p53_models, 
        file = paste0(results_folder, '/swan_bal_counts_p53_models.rda'))
saveRDS(swan_bal_counts_p53_sig_models, 
        file = paste0(results_folder, '/swan_bal_counts_p53_sig_models.rda'))
saveRDS(swan_unbal_p53_models, 
        file = paste0(results_folder, '/swan_unbal_p53_models.rda'))
saveRDS(swan_unbal_p53_sig_models, 
        file = paste0(results_folder, '/swan_unbal_p53_sig_models.rda'))
saveRDS(swan_cancer_intersection_models, 
        file = paste0(results_folder, '/swan_cancer_intersection_models.rda'))
saveRDS(swan_cancer_intersection_sig_models, 
        file = paste0(results_folder, '/swan_cancer_intersection_sig_models.rda'))
saveRDS(swan_cancer_union_models, 
        file = paste0(results_folder, '/swan_cancer_union_models.rda'))
saveRDS(swan_cancer_union_sig_models, 
        file = paste0(results_folder, '/swan_cancer_union_sig_models.rda'))
saveRDS(swan_p53_intersection_models, 
        file = paste0(results_folder, '/swan_p53_intersection_models.rda'))
saveRDS(swan_p53_intersection_sig_models, 
        file = paste0(results_folder, '/swan_p53_intersection_sig_models.rda'))
saveRDS(swan_p53_union_models, 
        file = paste0(results_folder, '/swan_p53_union_models.rda'))
saveRDS(swan_p53_union_sig_models, 
        file = paste0(results_folder, '/swan_p53_union_sig_models.rda'))


# remove model objects
rm(swan_bal_cancer_models, swan_bal_cancer_sig_models, swan_bal_counts_cancer_models, swan_bal_counts_cancer_sig_models,
   swan_unbal_cancer_models, swan_unbal_cancer_sig_models, swan_bal_p53_models, swan_bal_p53_sig_models, swan_bal_counts_p53_models, 
   swan_bal_counts_p53_sig_models, swan_unbal_p53_models, swan_unbal_p53_sig_models, swan_cancer_intersection_models, 
   swan_cancer_intersection_sig_models, swan_cancer_union_models, swan_cancer_union_sig_models, swan_p53_intersection_models, 
   swan_p53_intersection_sig_models, swan_p53_union_models, swan_p53_union_sig_models)

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



save.image('/home/benbrew/Desktop/temp_Results.RData')
load('/home/benbrew/Desktop/temp_Results.RData')

###########
# rbind tables and save RDA file
###########
quan_table <- rbind(quan_bal_cancer_table, quan_bal_cancer_sig_table, quan_bal_counts_cancer_table, quan_bal_counts_cancer_sig_table,
                   quan_unbal_cancer_table, quan_unbal_cancer_sig_table, quan_bal_p53_table, quan_bal_p53_sig_table, quan_bal_counts_p53_table, 
                   quan_bal_counts_p53_sig_table, quan_unbal_p53_table, quan_unbal_p53_sig_table, quan_cancer_intersection_table, 
                   quan_cancer_intersection_sig_table, quan_cancer_union_table, quan_cancer_union_sig_table, quan_p53_intersection_table, 
                   quan_p53_intersection_sig_table, quan_p53_union_table, quan_p53_union_sig_table)

# remove data 
rm(quan_bal_cancer_table, quan_bal_cancer_sig_table, quan_bal_counts_cancer_table, quan_bal_counts_cancer_sig_table,
   quan_unbal_cancer_table, quan_unbal_cancer_sig_table, quan_bal_p53_table, quan_bal_p53_sig_table, quan_bal_counts_p53_table, 
   quan_bal_counts_p53_sig_table, quan_unbal_p53_table, quan_unbal_p53_sig_table, quan_cancer_intersection_table, 
   quan_cancer_intersection_sig_table, quan_cancer_union_table, quan_cancer_union_sig_table, quan_p53_intersection_table, 
   quan_p53_intersection_sig_table, quan_p53_union_table, quan_p53_union_sig_table)

# remove table
saveRDS(quan_table, 
        file = paste0(results_folder, '/quan_table.rda'))

###########
# save all model objects as RDA file
###########
saveRDS(quan_bal_cancer_models, 
        file = paste0(results_folder, '/quan_bal_cancer_models.rda'))
saveRDS(quan_bal_cancer_sig_models, 
        file = paste0(results_folder, '/quan_bal_cancer_sig_models.rda'))
saveRDS(quan_bal_counts_cancer_models, 
        file = paste0(results_folder, '/quan_bal_counts_cancer_models.rda'))
saveRDS(quan_bal_counts_cancer_sig_models, 
        file = paste0(results_folder, '/quan_bal_counts_cancer_sig_models.rda'))
saveRDS(quan_unbal_cancer_models, 
        file = paste0(results_folder, '/quan_unbal_cancer_models.rda'))
saveRDS(quan_unbal_cancer_sig_models, 
        file = paste0(results_folder, '/quan_unbal_cancer_sig_models.rda'))
saveRDS(quan_bal_p53_models, 
        file = paste0(results_folder, '/quan_bal_p53_models.rda'))
saveRDS(quan_bal_p53_sig_models, 
        file = paste0(results_folder, '/quan_bal_p53_sig_models.rda'))
saveRDS(quan_bal_counts_p53_models, 
        file = paste0(results_folder, '/quan_bal_counts_p53_models.rda'))
saveRDS(quan_bal_counts_p53_sig_models, 
        file = paste0(results_folder, '/quan_bal_counts_p53_sig_models.rda'))
saveRDS(quan_unbal_p53_models, 
        file = paste0(results_folder, '/quan_unbal_p53_models.rda'))
saveRDS(quan_unbal_p53_sig_models, 
        file = paste0(results_folder, '/quan_unbal_p53_sig_models.rda'))
saveRDS(quan_cancer_intersection_models, 
        file = paste0(results_folder, '/quan_cancer_intersection_models.rda'))
saveRDS(quan_cancer_intersection_sig_models, 
        file = paste0(results_folder, '/quan_cancer_intersection_sig_models.rda'))
saveRDS(quan_cancer_union_models, 
        file = paste0(results_folder, '/quan_cancer_union_models.rda'))
saveRDS(quan_cancer_union_sig_models, 
        file = paste0(results_folder, '/quan_cancer_union_sig_models.rda'))
saveRDS(quan_p53_intersection_models, 
        file = paste0(results_folder, '/quan_p53_intersection_models.rda'))
saveRDS(quan_p53_intersection_sig_models, 
        file = paste0(results_folder, '/quan_p53_intersection_sig_models.rda'))
saveRDS(quan_p53_union_models, 
        file = paste0(results_folder, '/quan_p53_union_models.rda'))
saveRDS(quan_p53_union_sig_models, 
        file = paste0(results_folder, '/quan_p53_union_sig_models.rda'))


# remove model objects
rm(quan_bal_cancer_models, quan_bal_cancer_sig_models, quan_bal_counts_cancer_models, quan_bal_counts_cancer_sig_models,
   quan_unbal_cancer_models, quan_unbal_cancer_sig_models, quan_bal_p53_models, quan_bal_p53_sig_models, quan_bal_counts_p53_models, 
   quan_bal_counts_p53_sig_models, quan_unbal_p53_models, quan_unbal_p53_sig_models, quan_cancer_intersection_models, 
   quan_cancer_intersection_sig_models, quan_cancer_union_models, quan_cancer_union_sig_models, quan_p53_intersection_models, 
   quan_p53_intersection_sig_models, quan_p53_union_models, quan_p53_union_sig_models)

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



save.image('/home/benbrew/Desktop/temp_Results.RData')
load('/home/benbrew/Desktop/temp_Results.RData')

###########
# rbind tables and save RDA file
###########
funnorm_table <- rbind(funnorm_bal_cancer_table, funnorm_bal_cancer_sig_table, funnorm_bal_counts_cancer_table, funnorm_bal_counts_cancer_sig_table,
                   funnorm_unbal_cancer_table, funnorm_unbal_cancer_sig_table, funnorm_bal_p53_table, funnorm_bal_p53_sig_table, funnorm_bal_counts_p53_table, 
                   funnorm_bal_counts_p53_sig_table, funnorm_unbal_p53_table, funnorm_unbal_p53_sig_table, funnorm_cancer_intersection_table, 
                   funnorm_cancer_intersection_sig_table, funnorm_cancer_union_table, funnorm_cancer_union_sig_table, funnorm_p53_intersection_table, 
                   funnorm_p53_intersection_sig_table, funnorm_p53_union_table, funnorm_p53_union_sig_table)

# remove data 
rm(funnorm_bal_cancer_table, funnorm_bal_cancer_sig_table, funnorm_bal_counts_cancer_table, funnorm_bal_counts_cancer_sig_table,
   funnorm_unbal_cancer_table, funnorm_unbal_cancer_sig_table, funnorm_bal_p53_table, funnorm_bal_p53_sig_table, funnorm_bal_counts_p53_table, 
   funnorm_bal_counts_p53_sig_table, funnorm_unbal_p53_table, funnorm_unbal_p53_sig_table, funnorm_cancer_intersection_table, 
   funnorm_cancer_intersection_sig_table, funnorm_cancer_union_table, funnorm_cancer_union_sig_table, funnorm_p53_intersection_table, 
   funnorm_p53_intersection_sig_table, funnorm_p53_union_table, funnorm_p53_union_sig_table)

# remove table
saveRDS(funnorm_table, 
        file = paste0(results_folder, '/funnorm_table.rda'))

###########
# save all model objects as RDA file
###########
saveRDS(funnorm_bal_cancer_models, 
        file = paste0(results_folder, '/funnorm_bal_cancer_models.rda'))
saveRDS(funnorm_bal_cancer_sig_models, 
        file = paste0(results_folder, '/funnorm_bal_cancer_sig_models.rda'))
saveRDS(funnorm_bal_counts_cancer_models, 
        file = paste0(results_folder, '/funnorm_bal_counts_cancer_models.rda'))
saveRDS(funnorm_bal_counts_cancer_sig_models, 
        file = paste0(results_folder, '/funnorm_bal_counts_cancer_sig_models.rda'))
saveRDS(funnorm_unbal_cancer_models, 
        file = paste0(results_folder, '/funnorm_unbal_cancer_models.rda'))
saveRDS(funnorm_unbal_cancer_sig_models, 
        file = paste0(results_folder, '/funnorm_unbal_cancer_sig_models.rda'))
saveRDS(funnorm_bal_p53_models, 
        file = paste0(results_folder, '/funnorm_bal_p53_models.rda'))
saveRDS(funnorm_bal_p53_sig_models, 
        file = paste0(results_folder, '/funnorm_bal_p53_sig_models.rda'))
saveRDS(funnorm_bal_counts_p53_models, 
        file = paste0(results_folder, '/funnorm_bal_counts_p53_models.rda'))
saveRDS(funnorm_bal_counts_p53_sig_models, 
        file = paste0(results_folder, '/funnorm_bal_counts_p53_sig_models.rda'))
saveRDS(funnorm_unbal_p53_models, 
        file = paste0(results_folder, '/funnorm_unbal_p53_models.rda'))
saveRDS(funnorm_unbal_p53_sig_models, 
        file = paste0(results_folder, '/funnorm_unbal_p53_sig_models.rda'))
saveRDS(funnorm_cancer_intersection_models, 
        file = paste0(results_folder, '/funnorm_cancer_intersection_models.rda'))
saveRDS(funnorm_cancer_intersection_sig_models, 
        file = paste0(results_folder, '/funnorm_cancer_intersection_sig_models.rda'))
saveRDS(funnorm_cancer_union_models, 
        file = paste0(results_folder, '/funnorm_cancer_union_models.rda'))
saveRDS(funnorm_cancer_union_sig_models, 
        file = paste0(results_folder, '/funnorm_cancer_union_sig_models.rda'))
saveRDS(funnorm_p53_intersection_models, 
        file = paste0(results_folder, '/funnorm_p53_intersection_models.rda'))
saveRDS(funnorm_p53_intersection_sig_models, 
        file = paste0(results_folder, '/funnorm_p53_intersection_sig_models.rda'))
saveRDS(funnorm_p53_union_models, 
        file = paste0(results_folder, '/funnorm_p53_union_models.rda'))
saveRDS(funnorm_p53_union_sig_models, 
        file = paste0(results_folder, '/funnorm_p53_union_sig_models.rda'))


# remove model objects
rm(funnorm_bal_cancer_models, funnorm_bal_cancer_sig_models, funnorm_bal_counts_cancer_models, funnorm_bal_counts_cancer_sig_models,
   funnorm_unbal_cancer_models, funnorm_unbal_cancer_sig_models, funnorm_bal_p53_models, funnorm_bal_p53_sig_models, funnorm_bal_counts_p53_models, 
   funnorm_bal_counts_p53_sig_models, funnorm_unbal_p53_models, funnorm_unbal_p53_sig_models, funnorm_cancer_intersection_models, 
   funnorm_cancer_intersection_sig_models, funnorm_cancer_union_models, funnorm_cancer_union_sig_models, funnorm_p53_intersection_models, 
   funnorm_p53_intersection_sig_models, funnorm_p53_union_models, funnorm_p53_union_sig_models)






















##########
# most balanced
##########
bal_counts_cancer_intersection
bal_counts_cancer_intersection_sig
bal_counts_cancer_union
bal_counts_cancer_union_sig

bal_counts_p53_intersection
bal_counts_p53_intersection_sig
bal_counts_p53_union
bal_counts_p53_union_sig


#########
# cancer vs p53
#########
cancer_intersection
cancer_intersection_sig
cancer_union
cancer_union_sig

p53_intersection
p53_intersection_sig
p53_union
p53_union_sig


















# cancer balanced 
beta_swan_bal_cancer_models <- runModels(beta_swan, 
                                        random = F, 
                                        bump_hunter = T, 
                                        bump_hunter_data = beta_swan_bal_cancer_features)

beta_swan_bal_cancer_table <- extractResults(beta_swan_bal_cancer_models, 
                                            data_name = 'beta_swan_bal_cancer')


# cancer balanced counts
beta_swan_bal_counts_cancer_models <- runModels(beta_swan, 
                                               random = F, 
                                               bump_hunter = T, 
                                               bump_hunter_data = beta_swan_bal_counts_cancer_features)

beta_swan_bal_counts_cancer_table <- extractResults(beta_swan_bal_counts_cancer_models, 
                                                   data_name = 'beta_swan_bal_counts_cancer')

# cancer unbalanced 
beta_swan_unbal_cancer_models <- runModels(beta_swan, 
                                          random = F, 
                                          bump_hunter = T, 
                                          bump_hunter_data = beta_swan_unbal_cancer_features)

beta_swan_unbal_cancer_table <- extractResults(beta_swan_unbal_cancer_models, 
                                              data_name = 'beta_swan_unbal_cancer')

# p53 balanced 
beta_swan_bal_p53_models <- runModels(beta_swan, 
                                     random = F, 
                                     bump_hunter = T, 
                                     bump_hunter_data = beta_swan_bal_p53_features)

beta_swan_bal_p53_table <- extractResults(beta_swan_bal_p53_models, 
                                         data_name = 'beta_swan_bal_p53')


# p53 balanced counts
beta_swan_bal_counts_p53_models <- runModels(beta_swan, 
                                            random = F, 
                                            bump_hunter = T, 
                                            bump_hunter_data = beta_swan_bal_counts_p53_features)

beta_swan_bal_counts_p53_table <- extractResults(beta_swan_bal_counts_p53_models, 
                                                data_name = 'beta_swan_bal_counts_p53')

# p53 unbalanced 
beta_swan_unbal_p53_models <- runModels(beta_swan, 
                                       random = F, 
                                       bump_hunter = T, 
                                       bump_hunter_data = beta_swan_unbal_p53_features)

beta_swan_unbal_p53_table <- extractResults(beta_swan_unbal_p53_models, 
                                           data_name = 'beta_swan_unbal_p53')


##########
# rbind results
##########

beta_swan_results <- as.data.frame(rbind(beta_swan_bal_cancer_table, 
                                        beta_swan_bal_counts_cancer_table, 
                                        beta_swan_unbal_cancer_table,
                                        beta_swan_bal_p53_table, 
                                        beta_swan_bal_counts_p53_table, 
                                        beta_swan_unbal_p53_table))

# get rows and columns variable 
beta_swan_results <- getDims(beta_swan_results)

##########
# save results table and models for beta swan
##########
# first save results table or beta_swan
saveRDS(beta_swan_results, file = 
        paste0(results_folder, '/beta_swan_model_results.rda'))

# save models for beta_swan cancer
saveRDS(beta_swan_unbal_cancer_models, 
        file = paste0(results_folder, '/beta_swan_unbal_cancer_models.rda'))
saveRDS(beta_swan_bal_cancer_models, 
        file = paste0(results_folder, '/beta_swan_bal_cancer_models.rda'))
saveRDS(beta_swan_bal_counts_cancer_models,
        file = paste0(results_folder, '/beta_swan_bal_counts_cancer_models.rda'))

# save models for beta_swan p53
saveRDS(beta_swan_unbal_p53_models, 
        file = paste0(results_folder, '/beta_swan_unbal_p53_models.rda'))
saveRDS(beta_swan_bal_p53_models, 
        file = paste0(results_folder, '/beta_swan_bal_p53_models.rda'))
saveRDS(beta_swan_bal_counts_p53_models, 
        file = paste0(results_folder, '/beta_swan_bal_counts_p53_models.rda'))

##########
# remove unneeded objects
##########
rm(beta_swan_bal_cancer_features, 
   beta_swan_bal_counts_cancer_features, 
   beta_swan_unbal_cancer_features,
   beta_swan_bal_p53_features, 
   beta_swan_bal_counts_p53_features, 
   beta_swan_unbal_p53_features,
   beta_swan_bal_cancer_models, 
   beta_swan_bal_counts_cancer_models, 
   beta_swan_unbal_cancer_models,
   beta_swan_bal_p53_models, 
   beta_swan_bal_counts_p53_models, 
   beta_swan_unbal_p53_models,
   beta_swan_bal_cancer_table, 
   beta_swan_bal_counts_cancer_table, 
   beta_swan_unbal_cancer_table,
   beta_swan_bal_p53_table, 
   beta_swan_bal_counts_p53_table, 
   beta_swan_unbal_p53_table, 
   beta_swan_results)

######################################################################################################################

##########
# beta_swan
##########
# cancer balanced 
beta_swan_bal_cancer_models <- runModels(beta_swan, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_swan_bal_cancer_features)

beta_swan_bal_cancer_table <- extractResults(beta_swan_bal_cancer_models, 
                                             data_name = 'beta_swan_bal_cancer')


# cancer balanced counts
beta_swan_bal_counts_cancer_models <- runModels(beta_swan, 
                                                random = F, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_swan_bal_counts_cancer_features)

beta_swan_bal_counts_cancer_table <- extractResults(beta_swan_bal_counts_cancer_models, 
                                                    data_name = 'beta_swan_bal_counts_cancer')

# cancer unbalanced 
beta_swan_unbal_cancer_models <- runModels(beta_swan, 
                                           random = F, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_swan_unbal_cancer_features)

beta_swan_unbal_cancer_table <- extractResults(beta_swan_unbal_cancer_models, 
                                               data_name = 'beta_swan_unbal_cancer')

# p53 balanced 
beta_swan_bal_p53_models <- runModels(beta_swan, 
                                      random = F, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_swan_bal_p53_features)

beta_swan_bal_p53_table <- extractResults(beta_swan_bal_p53_models, 
                                          data_name = 'beta_swan_bal_p53')


# p53 balanced counts
beta_swan_bal_counts_p53_models <- runModels(beta_swan, 
                                             random = F, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_swan_bal_counts_p53_features)

beta_swan_bal_counts_p53_table <- extractResults(beta_swan_bal_counts_p53_models, 
                                                 data_name = 'beta_swan_bal_counts_p53')

# p53 unbalanced 
beta_swan_unbal_p53_models <- runModels(beta_swan, 
                                        random = F, 
                                        bump_hunter = T, 
                                        bump_hunter_data = beta_swan_unbal_p53_features)

beta_swan_unbal_p53_table <- extractResults(beta_swan_unbal_p53_models, 
                                            data_name = 'beta_swan_unbal_p53')

##########
# rbind results
##########
beta_swan_results <- as.data.frame(rbind(beta_swan_bal_cancer_table, 
                                         beta_swan_bal_counts_cancer_table, 
                                         beta_swan_unbal_cancer_table,
                                         beta_swan_bal_p53_table, 
                                         beta_swan_bal_counts_p53_table, 
                                         beta_swan_unbal_p53_table))


# get rows and columns variable 
beta_swan_results <- getDims(beta_swan_results)

##########
# save results table and models for beta swan
##########

# first save results table or beta_swan
saveRDS(beta_swan_results, file = paste0(results_folder, '/beta_swan_model_results.rda'))

# save models for beta_swan cancer
saveRDS(beta_swan_unbal_cancer_models, 
        file = paste0(results_folder, '/beta_swan_unbal_cancer_models.rda'))
saveRDS(beta_swan_bal_cancer_models, 
        file = paste0(results_folder, '/beta_swan_bal_cancer_models.rda'))
saveRDS(beta_swan_bal_counts_cancer_models, 
        file = paste0(results_folder, '/beta_swan_bal_counts_cancer_models.rda'))

# save models for beta_swan p53
saveRDS(beta_swan_unbal_p53_models, 
        file = paste0(results_folder, '/beta_swan_unbal_p53_models.rda'))
saveRDS(beta_swan_bal_p53_models, 
        file = paste0(results_folder, '/beta_swan_bal_p53_models.rda'))
saveRDS(beta_swan_bal_counts_p53_models, 
        file = paste0(results_folder, '/beta_swan_bal_counts_p53_models.rda'))

##########
# remove unneeded objects
##########

rm(beta_swan_bal_cancer_features, 
   beta_swan_bal_counts_cancer_features, 
   beta_swan_unbal_cancer_features,
   beta_swan_bal_p53_features, 
   beta_swan_bal_counts_p53_features, 
   beta_swan_unbal_p53_features,
   beta_swan_bal_cancer_models, 
   beta_swan_bal_counts_cancer_models, 
   beta_swan_unbal_cancer_models,
   beta_swan_bal_p53_models, 
   beta_swan_bal_counts_p53_models, 
   beta_swan_unbal_p53_models,
   beta_swan_bal_cancer_table, 
   beta_swan_bal_counts_cancer_table, 
   beta_swan_unbal_cancer_table,
   beta_swan_bal_p53_table, 
   beta_swan_bal_counts_p53_table, 
   beta_swan_unbal_p53_table, 
   beta_swan_results)


######################################################################################################################

##########
# beta_quan
##########

# cancer balanced 
beta_quan_bal_cancer_models <- runModels(beta_quan, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_quan_bal_cancer_features)

beta_quan_bal_cancer_table <- extractResults(beta_quan_bal_cancer_models, 
                                             data_name = 'beta_quan_bal_cancer')


# cancer balanced counts
beta_quan_bal_counts_cancer_models <- runModels(beta_quan, 
                                                random = F, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_quan_bal_counts_cancer_features)

beta_quan_bal_counts_cancer_table <- extractResults(beta_quan_bal_counts_cancer_models, 
                                                    data_name = 'beta_quan_bal_counts_cancer')

# cancer unbalanced 
beta_quan_unbal_cancer_models <- runModels(beta_quan, 
                                           random = F, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_quan_unbal_cancer_features)

beta_quan_unbal_cancer_table <- extractResults(beta_quan_unbal_cancer_models, 
                                               data_name = 'beta_quan_unbal_cancer')

# p53 balanced 
beta_quan_bal_p53_models <- runModels(beta_quan, 
                                      random = F, 
                                      bump_hunter = T, 
                                      bump_hunter_data = beta_quan_bal_p53_features)

beta_quan_bal_p53_table <- extractResults(beta_quan_bal_p53_models, 
                                          data_name = 'beta_quan_bal_p53')


# p53 balanced counts
beta_quan_bal_counts_p53_models <- runModels(beta_quan, 
                                             random = F, 
                                             bump_hunter = T, 
                                             bump_hunter_data = beta_quan_bal_counts_p53_features)

beta_quan_bal_counts_p53_table <- extractResults(beta_quan_bal_counts_p53_models, 
                                                 data_name = 'beta_quan_bal_counts_p53')

# p53 unbalanced 
beta_quan_unbal_p53_models <- runModels(beta_quan, 
                                        random = F, 
                                        bump_hunter = T, 
                                        bump_hunter_data = beta_quan_unbal_p53_features)

beta_quan_unbal_p53_table <- extractResults(beta_quan_unbal_p53_models, 
                                            data_name = 'beta_quan_unbal_p53')

##########
# rbind results
##########

beta_quan_results <- as.data.frame(rbind(beta_quan_bal_cancer_table, 
                                         beta_quan_bal_counts_cancer_table, 
                                         beta_quan_unbal_cancer_table,
                                         beta_quan_bal_p53_table, beta_quan_bal_counts_p53_table, 
                                         beta_quan_unbal_p53_table))

# get rows and columns variable 
beta_quan_results <- getDims(beta_quan_results)

##########
# save results table and models for beta quan
##########

# first save results table or beta_quan
saveRDS(beta_quan_results, file = paste0(results_folder, '/beta_quan_model_results.rda'))

# save models for beta_quan cancer
saveRDS(beta_quan_unbal_cancer_models, 
        file = paste0(results_folder, '/beta_quan_unbal_cancer_models.rda'))
saveRDS(beta_quan_bal_cancer_models, 
        file = paste0(results_folder, '/beta_quan_bal_cancer_models.rda'))
saveRDS(beta_quan_bal_counts_cancer_models, 
        file = paste0(results_folder, '/beta_quan_bal_counts_cancer_models.rda'))

# save models for beta_quan p53
saveRDS(beta_quan_unbal_p53_models, 
        file = paste0(results_folder, '/beta_quan_unbal_p53_models.rda'))
saveRDS(beta_quan_bal_p53_models, 
        file = paste0(results_folder, '/beta_quan_bal_p53_models.rda'))
saveRDS(beta_quan_bal_counts_p53_models, 
        file = paste0(results_folder, '/beta_quan_bal_counts_p53_models.rda'))

##########
# remove unneeded objects
##########

rm(beta_quan_bal_cancer_features, 
   beta_quan_bal_counts_cancer_features, 
   beta_quan_unbal_cancer_features,
   beta_quan_bal_p53_features, 
   beta_quan_bal_counts_p53_features, 
   beta_quan_unbal_p53_features,
   beta_quan_bal_cancer_models, 
   beta_quan_bal_counts_cancer_models, 
   beta_quan_unbal_cancer_models,
   beta_quan_bal_p53_models, 
   beta_quan_bal_counts_p53_models, 
   beta_quan_unbal_p53_models,
   beta_quan_bal_cancer_table, 
   beta_quan_bal_counts_cancer_table, 
   beta_quan_unbal_cancer_table,
   beta_quan_bal_p53_table, 
   beta_quan_bal_counts_p53_table, 
   beta_quan_unbal_p53_table, 
   beta_quan_results)

#########################################################################################################################

##########
# beta_funnorm
##########
# cancer balanced 
beta_funnorm_bal_cancer_models <- runModels(beta_funnorm, 
                                            random = F, 
                                            bump_hunter = T, 
                                            bump_hunter_data = beta_funnorm_bal_cancer_features)

beta_funnorm_bal_cancer_table <- extractResults(beta_funnorm_bal_cancer_models, 
                                                data_name = 'beta_funnorm_bal_cancer')


# cancer balanced counts
beta_funnorm_bal_counts_cancer_models <- runModels(beta_funnorm, 
                                                   random = F, 
                                                   bump_hunter = T, 
                                                   bump_hunter_data = beta_funnorm_bal_counts_cancer_features)

beta_funnorm_bal_counts_cancer_table <- extractResults(beta_funnorm_bal_counts_cancer_models, 
                                                       data_name = 'beta_funnorm_bal_counts_cancer')

# cancer unbalanced 
beta_funnorm_unbal_cancer_models <- runModels(beta_funnorm, 
                                              random = F, 
                                              bump_hunter = T, 
                                              bump_hunter_data = beta_funnorm_unbal_cancer_features)

beta_funnorm_unbal_cancer_table <- extractResults(beta_funnorm_unbal_cancer_models, 
                                                  data_name = 'beta_funnorm_unbal_cancer')

# p53 balanced 
beta_funnorm_bal_p53_models <- runModels(beta_funnorm, 
                                         random = F, 
                                         bump_hunter = T, 
                                         bump_hunter_data = beta_funnorm_bal_p53_features)

beta_funnorm_bal_p53_table <- extractResults(beta_funnorm_bal_p53_models, 
                                             data_name = 'beta_funnorm_bal_p53')


# p53 balanced counts
beta_funnorm_bal_counts_p53_models <- runModels(beta_funnorm, 
                                                random = F, 
                                                bump_hunter = T, 
                                                bump_hunter_data = beta_funnorm_bal_counts_p53_features)

beta_funnorm_bal_counts_p53_table <- extractResults(beta_funnorm_bal_counts_p53_models, data_name = 'beta_funnorm_bal_counts_p53')

# p53 unbalanced 
beta_funnorm_unbal_p53_models <- runModels(beta_funnorm, 
                                           random = F, 
                                           bump_hunter = T, 
                                           bump_hunter_data = beta_funnorm_unbal_p53_features)

beta_funnorm_unbal_p53_table <- extractResults(beta_funnorm_unbal_p53_models, 
                                               data_name = 'beta_funnorm_unbal_p53')

##########
# rbind results
##########
beta_funnorm_results <- as.data.frame(rbind(beta_funnorm_bal_cancer_table, 
                                            beta_funnorm_bal_counts_cancer_table, 
                                            beta_funnorm_unbal_cancer_table,
                                            beta_funnorm_bal_p53_table, 
                                            beta_funnorm_bal_counts_p53_table, 
                                            beta_funnorm_unbal_p53_table))


# get rows and columns variable 
beta_funnorm_results <- getDims(beta_funnorm_results)

##########
# save results table and models for beta funnorm
##########
# first save results table or beta_funnorm
saveRDS(beta_funnorm_results, file = paste0(results_folder, 
                                            '/beta_funnorm_model_results.rda'))

# save models for beta_funnorm cancer
saveRDS(beta_funnorm_unbal_cancer_models, 
        file = paste0(results_folder, '/beta_funnorm_unbal_cancer_models.rda'))
saveRDS(beta_funnorm_bal_cancer_models, 
        file = paste0(results_folder, '/beta_funnorm_bal_cancer_models.rda'))
saveRDS(beta_funnorm_bal_counts_cancer_models, 
        file = paste0(results_folder, '/beta_funnorm_bal_counts_cancer_models.rda'))

# save models for beta_funnorm p53
saveRDS(beta_funnorm_unbal_p53_models, 
        file = paste0(results_folder, '/beta_funnorm_unbal_p53_models.rda'))
saveRDS(beta_funnorm_bal_p53_models, 
        file = paste0(results_folder, '/beta_funnorm_bal_p53_models.rda'))
saveRDS(beta_funnorm_bal_counts_p53_models, 
        file = paste0(results_folder, '/beta_funnorm_bal_counts_p53_models.rda'))

##########
# remove unneeded objects
##########
rm(beta_funnorm_bal_cancer_features, 
   beta_funnorm_bal_counts_cancer_features, 
   beta_funnorm_unbal_cancer_features,
   beta_funnorm_bal_p53_features, 
   beta_funnorm_bal_counts_p53_features, 
   beta_funnorm_unbal_p53_features,
   beta_funnorm_bal_cancer_models, 
   beta_funnorm_bal_counts_cancer_models, 
   beta_funnorm_unbal_cancer_models,
   beta_funnorm_bal_p53_models, 
   beta_funnorm_bal_counts_p53_models, 
   beta_funnorm_unbal_p53_models,
   beta_funnorm_bal_cancer_table, 
   beta_funnorm_bal_counts_cancer_table, 
   beta_funnorm_unbal_cancer_table,
   beta_funnorm_bal_p53_table, 
   beta_funnorm_bal_counts_p53_table, 
   beta_funnorm_unbal_p53_table, 
   beta_funnorm_results)

######################################################################################################################

##########
# Now run each beta with random features - 100, 500, 1000, 2000, 10000
##########

##########
# 100 features 
##########
beta_swan_rand_models_100 <- runModels(beta_swan, 
                                  random = T, 
                                  bump_hunter = F,
                                  num_features = 100)

beta_swan_rand_table_100 <- extractResults(beta_swan_rand_models_100, 
                                          data_name = 'beta_swan_rand_100')

beta_swan_rand_models_100 <- runModels(beta_swan, 
                                   random = T, 
                                   bump_hunter = F,
                                   num_features = 100)

beta_swan_rand_table_100 <- extractResults(beta_swan_rand_models_100, 
                                       data_name = 'beta_swan_rand_100')

beta_quan_rand_models_100 <- runModels(beta_quan, 
                                   random = T, 
                                   bump_hunter = F,
                                   num_features = 100)

beta_quan_rand_table_100 <- extractResults(beta_quan_rand_models_100, 
                                           data_name = 'beta_quan_rand_100')

beta_funnorm_rand_models_100 <- runModels(beta_quan, 
                                          random = T, 
                                          bump_hunter = F,
                                          num_features = 100)

beta_funnorm_rand_table_100 <- extractResults(beta_funnorm_rand_models_100, 
                                              data_name = 'beta_funnorm_rand_100')


##########
# 500 features 
##########
beta_swan_rand_models_500 <- runModels(beta_swan, 
                                      random = T, 
                                      bump_hunter = F,
                                      num_features = 500)

beta_swan_rand_table_500 <- extractResults(beta_swan_rand_models_500, 
                                          data_name = 'beta_swan_rand_500')

beta_swan_rand_models_500 <- runModels(beta_swan, 
                                       random = T, 
                                       bump_hunter = F,
                                       num_features = 500)

beta_swan_rand_table_500 <- extractResults(beta_swan_rand_models_500, 
                                           data_name = 'beta_swan_rand_500')

beta_quan_rand_models_500 <- runModels(beta_quan, 
                                       random = T, 
                                       bump_hunter = F,
                                       num_features = 500)

beta_quan_rand_table_500 <- extractResults(beta_quan_rand_models_500, 
                                           data_name = 'beta_quan_rand_500')

beta_funnorm_rand_models_500 <- runModels(beta_quan, 
                                          random = T, 
                                          bump_hunter = F,
                                          num_features = 500)

beta_funnorm_rand_table_500 <- extractResults(beta_funnorm_rand_models_500, 
                                              data_name = 'beta_funnorm_rand_500')


##########
# 1000 features 
##########
beta_swan_rand_models_1000 <- runModels(beta_swan, 
                                      random = T, 
                                      bump_hunter = F,
                                      num_features = 1000)

beta_swan_rand_table_1000 <- extractResults(beta_swan_rand_models_1000, 
                                          data_name = 'beta_swan_rand_1000')

beta_swan_rand_models_1000 <- runModels(beta_swan, 
                                       random = T, 
                                       bump_hunter = F,
                                       num_features = 1000)

beta_swan_rand_table_1000 <- extractResults(beta_swan_rand_models_1000, 
                                           data_name = 'beta_swan_rand_1000')

beta_quan_rand_models_1000 <- runModels(beta_quan, 
                                       random = T, 
                                       bump_hunter = F,
                                       num_features = 1000)

beta_quan_rand_table_1000 <- extractResults(beta_quan_rand_models_1000, 
                                           data_name = 'beta_quan_rand_1000')

beta_funnorm_rand_models_1000 <- runModels(beta_quan, 
                                          random = T, 
                                          bump_hunter = F,
                                          num_features = 1000)

beta_funnorm_rand_table_1000 <- extractResults(beta_funnorm_rand_models_1000, 
                                              data_name = 'beta_funnorm_rand_1000')

##########
# 2000 features 
##########
beta_swan_rand_models_2000 <- runModels(beta_swan, 
                                      random = T, 
                                      bump_hunter = F,
                                      num_features = 2000)

beta_swan_rand_table_2000 <- extractResults(beta_swan_rand_models_2000, 
                                          data_name = 'beta_swan_rand_2000')

beta_swan_rand_models_2000 <- runModels(beta_swan, 
                                       random = T, 
                                       bump_hunter = F,
                                       num_features = 2000)

beta_swan_rand_table_2000 <- extractResults(beta_swan_rand_models_2000, 
                                           data_name = 'beta_swan_rand_2000')

beta_quan_rand_models_2000 <- runModels(beta_quan, 
                                       random = T, 
                                       bump_hunter = F,
                                       num_features = 2000)

beta_quan_rand_table_2000 <- extractResults(beta_quan_rand_models_2000, 
                                           data_name = 'beta_quan_rand_2000')

beta_funnorm_rand_models_2000 <- runModels(beta_quan, 
                                          random = T, 
                                          bump_hunter = F,
                                          num_features = 2000)

beta_funnorm_rand_table_2000 <- extractResults(beta_funnorm_rand_models_2000, 
                                              data_name = 'beta_funnorm_rand_2000')



##########
# 10000 features 
##########
beta_swan_rand_models_10000 <- runModels(beta_swan, 
                                      random = T, 
                                      bump_hunter = F,
                                      num_features = 10000)

beta_swan_rand_table_10000 <- extractResults(beta_swan_rand_models_10000, 
                                          data_name = 'beta_swan_rand_10000')

beta_swan_rand_models_10000 <- runModels(beta_swan, 
                                       random = T, 
                                       bump_hunter = F,
                                       num_features = 10000)

beta_swan_rand_table_10000 <- extractResults(beta_swan_rand_models_10000, 
                                           data_name = 'beta_swan_rand_10000')

beta_quan_rand_models_10000 <- runModels(beta_quan, 
                                       random = T, 
                                       bump_hunter = F,
                                       num_features = 10000)

beta_quan_rand_table_10000 <- extractResults(beta_quan_rand_models_10000, 
                                           data_name = 'beta_quan_rand_10000')

beta_funnorm_rand_models_10000 <- runModels(beta_quan, 
                                          random = T, 
                                          bump_hunter = F,
                                          num_features = 10000)

beta_funnorm_rand_table_10000 <- extractResults(beta_funnorm_rand_models_10000, 
                                              data_name = 'beta_funnorm_rand_10000')

#HERE
##########
# rbind results
##########
beta_rand_results <- as.data.frame(rbind(beta_swan_rand_table_100, 
                                         beta_swan_rand_table_100, 
                                         beta_quan_rand_table_100,
                                         beta_funnorm_rand_table_100,
                                         beta_swan_rand_table_500, 
                                         beta_swan_rand_table_500, 
                                         beta_quan_rand_table_500,
                                         beta_funnorm_rand_table_500,
                                         beta_swan_rand_table_1000, 
                                         beta_swan_rand_table_1000, 
                                         beta_quan_rand_table_1000,
                                         beta_funnorm_rand_table_1000,
                                         beta_swan_rand_table_2000, 
                                         beta_swan_rand_table_2000, 
                                         beta_quan_rand_table_2000,
                                         beta_funnorm_rand_table_2000,
                                         beta_swan_rand_table_10000, 
                                         beta_swan_rand_table_10000, 
                                         beta_quan_rand_table_10000,
                                         beta_funnorm_rand_table_10000))


# get rows and columns variable 
beta_rand_results <- getDims(beta_rand_results)

##########
# save results table and models for beta swan
##########
# first save results table or beta_swan
saveRDS(beta_rand_results, file = paste0(results_folder, 
                                         '/beta_rand_model_results.rda'))

# save models for 100
saveRDS(beta_swan_rand_models_100, 
        file = paste0(results_folder, '/beta_swan_rand_models_100.rda'))
saveRDS(beta_swan_rand_models_100, 
        file = paste0(results_folder, '/beta_swan_rand_models_100.rda'))
saveRDS(beta_quan_rand_models_100, 
        file = paste0(results_folder, '/beta_quan_rand_models_100.rda'))
saveRDS(beta_funnorm_rand_models_100, 
        file = paste0(results_folder, '/beta_funnorm_rand_models_100.rda'))

# save models for 500
saveRDS(beta_swan_rand_models_500, 
        file = paste0(results_folder, '/beta_swan_rand_models_500.rda'))
saveRDS(beta_swan_rand_models_500, 
        file = paste0(results_folder, '/beta_swan_rand_models_500.rda'))
saveRDS(beta_quan_rand_models_500, 
        file = paste0(results_folder, '/beta_quan_rand_models_500.rda'))
saveRDS(beta_funnorm_rand_models_500, 
        file = paste0(results_folder, '/beta_funnorm_rand_models_500.rda'))

# save models for 1000
saveRDS(beta_swan_rand_models_1000, 
        file = paste0(results_folder, '/beta_swan_rand_models_1000.rda'))
saveRDS(beta_swan_rand_models_1000, 
        file = paste0(results_folder, '/beta_swan_rand_models_1000.rda'))
saveRDS(beta_quan_rand_models_1000, 
        file = paste0(results_folder, '/beta_quan_rand_models_1000.rda'))
saveRDS(beta_funnorm_rand_models_1000, 
        file = paste0(results_folder, '/beta_funnorm_rand_models_1000.rda'))

# save models for 2000
saveRDS(beta_swan_rand_models_2000, 
        file = paste0(results_folder, '/beta_swan_rand_models_2000.rda'))
saveRDS(beta_swan_rand_models_2000, 
        file = paste0(results_folder, '/beta_swan_rand_models_2000.rda'))
saveRDS(beta_quan_rand_models_2000, 
        file = paste0(results_folder, '/beta_quan_rand_models_2000.rda'))
saveRDS(beta_funnorm_rand_models_2000, 
        file = paste0(results_folder, '/beta_funnorm_rand_models_2000.rda'))

# save models for 10000
saveRDS(beta_swan_rand_models_10000, 
        file = paste0(results_folder, '/beta_swan_rand_models_10000.rda'))
saveRDS(beta_swan_rand_models_10000, 
        file = paste0(results_folder, '/beta_swan_rand_models_10000.rda'))
saveRDS(beta_quan_rand_models_10000, 
        file = paste0(results_folder, '/beta_quan_rand_models_10000.rda'))
saveRDS(beta_funnorm_rand_models_10000, 
        file = paste0(results_folder, '/beta_funnorm_rand_models_10000.rda'))

