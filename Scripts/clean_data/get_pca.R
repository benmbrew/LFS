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
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_batch_m_sub.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_batch_m_sub.rda')))
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_batch_m_sub.rda')))

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


##########
# remove outliers
##########
betaCases <- removeOutlier(betaCases, 
                           cases = T, 
                           controls = F, 
                           val =F)

betaControls <- removeOutlier(betaControls, 
                              cases = F, 
                              controls = T, 
                              val =F)

betaValid <- removeOutlier(betaValid, 
                           cases = F, 
                           controls = F, 
                           val = T)


# ##########
# # convert to m values 
# ##########
# 
# betaCases <- get_m_values(betaCases, probe_start = 8)
# betaControls <- get_m_values(betaControls, probe_start = 8)
# betaValid <- get_m_values(betaValid, probe_start = 8)

# # read data for models
# betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m.rda')))
# betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new_m.rda')))
# betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new_m.rda')))

###########
# remove columns with Inf
###########
betaCases <- removeInf(betaCases, probe_start = 8)
betaControls <- removeInf(betaControls, probe_start = 8)
betaValid <- removeInf(betaValid, probe_start = 8)


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
setwd('~/Desktop/')
pdf(paste0(method, '_plots_m_values_funnorm_no_outliers.pdf'))
# cases full
data_cols <- c('cancer_diagnosis_diagnoses', 
               'gender')

# cases 

par(mfrow = c(1,1))

for (i in data_cols) {
  
  getPCA(pca_data = betaCases, 
         column_name = i, 
         gene_start = 8, 
         pca1 = 1, 
         pca2 = 2, 
         name = paste0('CasesFull_', i), 
         use_legend = F) # 3358
  
  print(i)
  
}

for (i in data_cols) {
  
  getPCA(pca_data = betaCasesMod, 
         column_name = i, 
         gene_start = 8, 
         pca1 = 1, 
         pca2 = 2, 
         name = paste0('CasesMod_', i), 
         use_legend = F) # 3358
  
  print(i)
  
}

for (i in data_cols[-1]) {
  
  getPCA(pca_data = betaControlsMod, 
         column_name = i, 
         gene_start = 8, 
         pca1 = 1, 
         pca2 = 2, 
         name = paste0('ControlsMod_', i), 
         use_legend = F) # 3358
  
  print(i)
  
}

for (i in data_cols) {
  
  getPCA(pca_data = betaValidMod, 
         column_name = i, 
         gene_start = 8, 
         pca1 = 1, 
         pca2 = 2, 
         name = paste0('ValidMod_', i), 
         use_legend = F) # 3358
  
  print(i)
  
}

# combined data cols
combined_cols <- c('cancer_diagnosis_diagnoses', 
                   'gender',
                   'batch')


# get full data combined
for (i in combined_cols) {
  
  getPCA(pca_data = full_data ,
         column_name = i, 
         gene_start = 8, 
         pca1 = 1, 
         pca2 = 2, 
         name = paste0('combined_data_full_', i, '_'), 
         use_legend = T) # 3358
  
  print(i)
  
}

# get mod data combined

for (i in combined_cols) {
  
  getPCA(pca_data = mod_data ,
         column_name = i, 
         gene_start = 8, 
         pca1 = 1, 
         pca2 = 2, 
         name = paste0('combined_data_mod_', i, '_'), 
         use_legend = F) # 3358
  
  print(i)
  
}

dev.off()


# determine the ranking of data non missingness for betaCasesFull (201) 

# age diagnosis
length(which(is.na(betaCases$age_diagnosis)))
length(which(is.na(betaCases$age_sample_collection)))
length(which(is.na(betaCases$cancer_diagnosis_diagnoses)))
length(which(is.na(betaCases$gender)))
length(which(is.na(betaCases$sentrix_id)))







