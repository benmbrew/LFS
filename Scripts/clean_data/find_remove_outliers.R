## this script will read in batch data from get_cases, get_controls, or get_valid and explore potential batches and outliers

##########
# initialize libraries
##########
library(tidyverse)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
clin_data <- paste0(data_folder, '/clin_data')

# get method 
method = 'funnorm'

##########
# load data
##########
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_batch_m.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_batch_m.rda')))
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_batch_m.rda')))

##########
# read in clinical data
##########
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# clean clinical ids
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

########## 
# remove NA
##########
betaCases <- removeNA(betaCases, probe_start = 8)
betaControls <- removeNA(betaControls, probe_start = 8)
betaValid <- removeNA(betaValid, probe_start = 8)


###########
# remove columns with Inf
###########
betaCases <- removeInf(betaCases, probe_start = 8)
betaControls <- removeInf(betaControls, probe_start = 8)
betaValid <- removeInf(betaValid, probe_start = 8)


##########
# if you use transform data, homogenize here
##########
intersect_names <- Reduce(intersect, list(colnames(betaCases)[8:ncol(betaCases)], 
                                          colnames(betaControls)[8:ncol(betaControls)], 
                                          colnames(betaValid)[8:ncol(betaValid)]))




# add column in
betaCases$batch <- 'cases'
betaControls$batch <- 'controls_mod'
betaValid$batch <- 'valid_mod'


# cases
betaCases <- betaCases[, c('ids',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender',
                           'batch',
                           'sentrix_id',
                           intersect_names)]

# controls
betaControls   <- betaControls[, c('ids',
                                       'age_diagnosis', 
                                       'age_sample_collection', 
                                       'cancer_diagnosis_diagnoses', 
                                       'gender', 
                                       'batch',
                                       'sentrix_id',
                                       intersect_names)]

#validation
betaValid <- betaValid [, c('ids',
                            'age_diagnosis', 
                            'age_sample_collection', 
                            'cancer_diagnosis_diagnoses', 
                            'gender',
                            'batch',
                            'sentrix_id',
                            intersect_names)]



# get features sites
pca_cases <- prcomp(betaCases[,8:ncol(betaCases)])
pca_controls <- prcomp(betaControls[,8:ncol(betaControls)])
pca_valid <- prcomp(betaValid[,8:ncol(betaValid)])


get_outlier <- function(pca_object, clin_vars){
  
  pca1 <- pca_object$x[,1]
  pca2 <- pca_object$x[,2]
  
  pca1_index <- which(abs(scale(pca1)) >= 2)
  pca2_index <- which(abs(scale(pca2)) >= 2)
  
  outlier_index <- append(pca1_index, pca2_index)
  
  # get coresponding ids
  outliers <- clin_vars[outlier_index,]
  
  return(outliers)
  
}

# get outliers
cases <- get_outlier(pca_object = pca_cases, clin_vars = betaCases[, 1:7])
controls <- get_outlier(pca_object = pca_controls, clin_vars = betaControls[, 1:7])
valid <- get_outlier(pca_object = pca_valid, clin_vars = betaValid[, 1:7])

# add together and get clinical info
cases_clin <- inner_join(cases, clin, by = 'ids')
controls_clin <- inner_join(controls, clin, by = 'ids')
valid_clin <- inner_join(valid, clin, by = 'ids')

full_data <- rbind(cases_clin,
                   controls_clin,
                   valid_clin)


# save identifier to suber rgset by good samples before normalization using funnorm
saveRDS(full_data, paste0(model_data, paste0('/', method, '_', 'outliers.rda')))

