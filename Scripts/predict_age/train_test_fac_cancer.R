

##########
# load libraries
##########
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(biovizBase)m
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)
library(preprocessCore)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
# library(Metrics)
# library(ModelMetrics)
library(doParallel) 
library(nnet)
library(dplyr)
library(bumphunter)
library(sqldf)
library(e1071)
library(reshape2)


##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
idat_data_case <- paste0(methyl_data, '/raw_files')
idat_data_con <- paste0(methyl_data, '/controls')
idat_data_val <- paste0(data_folder, '/methyl_data/validation/idat_files')
model_data <- paste0(data_folder, '/model_data')
map_data <- paste0(data_folder, '/methyl_data/validation')


##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# read in clinical data
##########
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# clean clinical ids
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

##########
# cases 
##########

# cases batch1
id_map <- read.csv(paste0(methyl_data, '/ids_map.csv'), stringsAsFactors = F)

#cases batch2
id_map_other <- read.csv(paste0(methyl_data, '/batch_2014.csv'), stringsAsFactors = F)
id_map_other$Project <- NULL

# combine id_map and id_map_other
id_map <- rbind(id_map, id_map_other)
rm(id_map_other)

# clean id map
id_map <- cleanIdMap(id_map)


##########
# Controls batch1
##########
id_map_con <- read.csv(paste0(methyl_data, '/ids_map_controls.csv'), stringsAsFactors = F)

# clean idmap
id_map_con <- cleanIdMap(id_map_con)


##########
# valid
##########
id_map_val <- read.csv(paste0(map_data, '/id_map_validation.csv'), stringsAsFactors = F)

# homogenize valid map data with cases and controls
id_map_val <- id_map_val[, c('Sample.ID', 'Sample.Group', 'Sentrix.Barcode', 'Sample.Section',
                             'Project', 'Pool_ID', 'Sample_Well')]

# sub out '.' for '_'
colnames(id_map_val) <- gsub('.', '_', colnames(id_map_val), fixed = T)

# change 'Sample_ID' to 'Sample_Name' and 'Sentrix_Barcode' to 'Sentrix_ID'
colnames(id_map_val)[1] <- 'Sample_Name'
colnames(id_map_val)[3] <- 'Sentrix_ID'
colnames(id_map_val)[4] <- 'Sentrix_Position'
colnames(id_map_val)[5] <- 'Sample_Plate'

# clean idmap
id_map_val <- cleanIdMap(id_map_val)


##########
# read in meth array
##########
rgCases <- read.metharray.exp(idat_data_case)

rgControls <- read.metharray.exp(idat_data_con)

rgValid <- read.metharray.exp(idat_data_val)



##########
# remove outliers (previously determined) from rgset before normalization
##########
rgControls <- remove_outliers(rgSet = rgControls,
                              id_map = id_map_con,
                              method = 'doesnt_matter',
                              type = 'controls')

rgValid <- remove_outliers(rgSet = rgValid,
                           id_map = id_map_val,
                           method = 'doesnt_matter',
                           type = 'valid')

# save.image('~/Desktop/temp_rg_data.RData')
load('~/Desktop/temp_rg_data.RData')


##########
# subset data - remove controls probes on each data set only if raw preprocessing
##########

# cases 
rg_cases <- subset_rg_set(rg_set = rgCases, 
                          keep_gender = F,
                          keep_controls = F, 
                          keep_snps = F, 
                          get_island = "Island", 
                          get_chr = NULL, 
                          get_type = NULL)

# controls
rg_controls <- subset_rg_set(rg_set = rgControls, 
                             keep_gender = F,
                             keep_controls = F, 
                             keep_snps = F, 
                             get_island = "Island",
                             get_chr = NULL, 
                             get_type = NULL)

# valid 
rg_valid <- subset_rg_set(rg_set = rgValid, 
                          keep_gender = F,
                          keep_controls = F, 
                          keep_snps = F, 
                          get_island = "Island", 
                          get_chr = NULL, 
                          get_type = NULL)



# rg_cases = rg_cases
# rg_controls = rg_controls
# rg_valid = rg_valid
# k_folds = 5
# m_beta_thresh = 0.5
# method = 'raw'
# controls = 'normal'
# gender = T 
# p53 = F
# max_columns = 10000

full_pipeline_cancer <- function(rg_cases, 
                                 rg_controls, 
                                 rg_valid,
                                 gender, 
                                 p53,
                                 k_folds,
                                 method,
                                 max_columns) {
  
  # list to store cv results
  results_list <- list()
  clin_results <- list()

  # arguments or functions for subsetting data in different ways 
  # for example, remove chromosome, type I and type II probes, Chr17, remove Ch6
  
  # preprocess controls and valid
  m_cases <- preprocessMethod(rg_cases, preprocess = method, only_m_values = T)
  m_controls <- preprocessMethod(rg_controls, preprocess = method, only_m_values = T)
  m_valid <- preprocessMethod(rg_valid, preprocess = method, only_m_values = T)
  
  
  # do cases first (will return list of 2, second element is old controls)
  m_cases <- process_rg_set_single(beta_data = m_cases[1:max_columns,], 
                                         id_map = id_map, 
                                         clin = clin)
  # get controls
  m_controls <- process_rg_set_single(beta_data = m_controls[1:max_columns,], 
                                          id_map = id_map_con, 
                                          clin = clin)
  
  # get valid
  m_valid  <- process_rg_set_single(beta_data = m_valid[1:max_columns,], 
                                       id_map = id_map_val, 
                                       clin = clin)
  
  ##########
  # function to remove infinite values
  ##########
  
  # cases
  m_cases <- removeInf(m_cases, probe_start = 10)

  # controls 
  m_controls <- removeInf(m_controls, probe_start = 10)
  
  # valid 
  m_valid <- removeInf(m_valid, probe_start = 10)
  
  ##########
  # bumphunter 
  ##########
  intersect_names <- Reduce(intersect, list(colnames(m_cases)[10:ncol(m_cases)],
                                            colnames(m_controls)[10:ncol(m_controls)],
                                            colnames(m_valid)[10:ncol(m_valid)]))
  
  
  # train data 
  m_cases <- m_cases[, c('ids',
                         'p53_germline',
                         'cancer_diagnosis_diagnoses',
                         'age_diagnosis',
                         'age_sample_collection',
                         'gender',
                         'sentrix_id',
                         'family_name',
                         'tm_donor_',
                         intersect_names)]
  
  # controls data 
  m_controls <- m_controls[, c('ids',
                               'p53_germline',
                               'cancer_diagnosis_diagnoses',
                               'age_diagnosis',
                               'age_sample_collection',
                               'gender',
                               'sentrix_id',
                               'family_name',
                               'tm_donor_',
                               intersect_names)]
  
  # train data 
  m_valid <- m_valid[, c('ids',
                         'p53_germline',
                         'cancer_diagnosis_diagnoses',
                         'age_diagnosis',
                         'age_sample_collection',
                         'gender',
                         'sentrix_id',
                         'family_name',
                         'tm_donor_',
                         intersect_names)]
  
  # combined data sets
  m_full <- rbind(m_cases, m_controls, m_valid)
  
  # remove dups 
  m_full_mod <- m_full[!duplicated(m_full$ids),]
  m_full_mod <- m_full_mod[!duplicated(m_full_mod$tm_donor_),]
  
  rm(m_cases, m_controls, m_valid)
  rm(rg_cases, rg_controls, rg_valid, rgCases, rgControls, rgValid)
  
  # remove columns 
  m_full$ids <- m_full$sentrix_id <- m_full$family_name <- m_full$tm_donor_ <-
    m_full$age_diagnosis <- m_full$age_sample_collection <- NULL
  
  
  if(gender) {
    # get gender variable for each data set
    m_full <- cbind(as.data.frame(class.ind(m_full$gender)), m_full)
    
  }
  
  if (p53) {
    
    # get dummy variale for p53
    m_full <- cbind(as.data.frame(class.ind(m_full$p53_germline)), m_full)
    
  } else {
    
    # subset to LFS only
    m_full <- m_full[grepl('Mut', m_full$p53_germline),]
  }
  
  # remove NAs
  m_full <- m_full[complete.cases(m_full),]
  
  # get folds 
  fold_vec <- sample(1:k_folds, nrow(m_full), replace = T)
  
  
  # train and test
  for(i in 1:k_folds) {
    
    # get train and test index
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    # get training rgsets 
    m_train <- m_full[train_index,]
    m_test <- m_full[test_index,]
    
    # test clin
    test_clin <- as.data.frame(m_test[, c('p53_germline', 'cancer_diagnosis_diagnoses', 'gender')])
    
    # function to predict with all test, controls, controls old, and valid
    # HERE
    mod_result <- predCancer(training_dat = m_train, 
                             test_dat = m_test, 
                             clin_dat = test_clin,
                             bh_features = intersect_names,
                             gender = gender,
                             p53 = p53)
    
    
    # returns test_stats, test_stats_age, test_stats_controls, test_stats_valid, test_stats_age_valid
    temp.results <- get_class_results_cancer(mod_result, 
                                             dims_of_dat = length(intersect_names), 
                                             mod_name = 'enet', 
                                             gender = gender,
                                             p53 = p53)
    
    results_list[[i]] <- temp.results[[1]]
    clin_results[[i]] <- temp.results[[2]]
  
    print(i)
    
  }
  
  #combine all the folds 
  results_final <- do.call(rbind, results_list)
  results_final_clin <- do.call(rbind, clin_results)
  
  
  return(list(results_final, results_final_clin))
}


##########
# fixed variables
##########

method = 'raw'
k_folds = 5
gender = T
p53 = F
max_columns = 50000


#### TO FIX AND CHECK
# remove sex chromosomes and just control for gender
# check to see that you are removing cancer probes from bumphunter in the pipeline

# run full pipeline 
full_results <- full_pipeline_cancer(rg_cases = rgCases, 
                                     rg_controls = rgControls, 
                                     rg_valid = rgValid, 
                                     gender = T, 
                                     p53 = F, 
                                     k_folds = 5, 
                                     method = 'raw', 
                                     max_columns = 50000)

# save results 
saveRDS(full_results, paste0('~/Desktop/', method, '_', k_folds, '_', 
                             max_age, '_', max_columns, '_', 
                             m_beta_thresh, '_', combined, '_', 
                             combined_data_type, '_full_results_cancer.rda'))
