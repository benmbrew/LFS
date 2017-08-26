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
method = 'raw'

##########
# load data
##########
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new.rda')))
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new.rda')))

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

#########
# binarize age variables 
#########
age_binary <- 
  function(dat, type, cutoff) {
    
    if (type == 'cases') {
      
      dat$age_diagnosis <- ifelse(dat$age_diagnosis > cutoff, 'yes', 'no')
      dat$age_sample_collection <- ifelse(dat$age_sample_collection > cutoff, 'yes', 'no')
      
      return(dat)
      
      
    } else {
      
      dat$age_sample_collection <- ifelse(dat$age_sample_collection > cutoff, 'yes', 'no')
      
      return(dat)
      
    }
    
}


##########
# if you use transform data, homogenize here
##########
intersect_names <- Reduce(intersect, list(colnames(betaCases)[8:ncol(betaCases)], 
                                          colnames(betaControls)[8:ncol(betaControls)], 
                                          colnames(betaValid)[8:ncol(betaValid)]))


# cases age of diagnosis 
betaCases <- age_binary(betaCases, 
                        'cases', 
                        cutoff = 48)

# controls age of diagnosis 
betaControls <- age_binary(betaControls, 
                           'controls', 
                            cutoff = 48)
# valid age of diagnosis 
betaValid <- age_binary(betaValid, 
                        'cases', 
                        cutoff = 48)

# get mod data
betaCasesMod <- getModData(betaCases)

# remove valid samples that are already present in betaCases
betaValidMod <- betaValid[!betaValid$ids %in% betaCases$ids,]

# homogenize data names by creating "Mod"
betaControlsMod <- betaControls
rm(betaControls, betaValid)

# add column in
betaCases$batch <- 'cases'
betaCasesMod$batch <- 'cases_mod'
betaControlsMod$batch <- 'controls_mod'
betaValidMod$batch <- 'valid_mod'


# cases
betaCases <- betaCases[, c('ids',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender',
                           'batch',
                           'sentrix_id',
                            intersect_names)]

# cases_mod
betaCasesMod <- betaCasesMod[, c('ids',
                                 'age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender',
                                 'batch',
                                 'sentrix_id',
                                  intersect_names)]
# controls
betaControlsMod <- betaControlsMod[, c('ids',
                                       'age_diagnosis', 
                                       'age_sample_collection', 
                                       'cancer_diagnosis_diagnoses', 
                                       'gender', 
                                       'batch',
                                       'sentrix_id',
                                        intersect_names)]

#validation
betaValidMod <- betaValidMod[, c('ids',
                                 'age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender',
                                 'batch',
                                 'sentrix_id',
                                 intersect_names)]


# combine for model data
full_data <- rbind(betaCases,
                   betaControlsMod,
                   betaValidMod)

# mod data
mod_data <- rbind(betaCasesMod,
                  betaControlsMod,
                  betaValidMod)


##########
# pca of age of onset, age of sample collection, gender. sentrix_id 
##########

# cases full

getPCA(pca_data = betaCases, 
       column_name = 'age_diagnosis', 
       gene_start = 7, 
       pca1 = 1, 
       pca2 = 2, 
       name = 'Full Cases, Age of Onset', 
       use_legend = T) # 3358

getPCA(pca_data = betaCases, 
       column_name = 'age_sample_collection', 
       gene_start = 7, 
       pca1 = 1, 
       pca2 = 2, 
       name = 'Full Cases, Age Collection', 
       use_legend = T) # 3358

getPCA(pca_data = betaCases, 
       column_name = 'cancer_diagnosis_diagnoses', 
       gene_start = 7, 
       pca1 = 1, 
       pca2 = 2, 
       name = 'Full Cases, Cancer', 
       use_legend = T) # 3358




