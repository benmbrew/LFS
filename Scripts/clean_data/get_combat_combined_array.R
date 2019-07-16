## this script will read in batch data from get_cases, get_controls, or get_valid and explore potential batches and outliers

##########
# initialize libraries
##########
library(BiocParallel)
library(doParallel)
library(sva)

registerDoParallel(1)

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
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))




##########
# load data
##########

# read in data
cases_con_beta_full <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_con_beta_m.rda')))
cases_val_beta_full <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_val_beta_m.rda')))

# seperate data
cases_con <- cases_con_beta_full[!grepl('^200', cases_con_beta_full$sentrix_id), ]
con_cases <- cases_con_beta_full[grepl('^200', cases_con_beta_full$sentrix_id), ]

cases_val <- cases_val_beta_full[!grepl('^20', cases_val_beta_full$sentrix_id), ]
val_cases <- cases_val_beta_full[grepl('^20', cases_val_beta_full$sentrix_id), ]


##########
# if you use transform data, homogenize here
##########
intersect_names <- colnames(cases_con)[9:ncol(cases_con)]


# add column in
cases_con$batch <- 'cases'
cases_val$batch <- 'cases'

con_cases$batch <- 'controls'
val_cases$batch <- 'valid'


# cases
cases_con <- cases_con[, c('ids',
                           'p53_germline',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender',
                           'batch',
                           'family_name',
                           intersect_names)]

# cases
cases_val<- cases_val[, c('ids',
                           'p53_germline',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender',
                           'batch',
                           'family_name',
                           intersect_names)]

# cases
con_cases<- con_cases[, c('ids',
                          'p53_germline',
                          'age_diagnosis', 
                          'age_sample_collection', 
                          'cancer_diagnosis_diagnoses', 
                          'gender',
                          'batch',
                          'family_name',
                          intersect_names)]

# cases
val_cases<- val_cases[, c('ids',
                          'p53_germline',
                          'age_diagnosis', 
                          'age_sample_collection', 
                          'cancer_diagnosis_diagnoses', 
                          'gender',
                          'batch',
                          'family_name',
                          intersect_names)]




# combine for model data
full_data_con<- rbind(cases_con,
                      con_cases)

full_data_val<- rbind(cases_val,
                      val_cases)






full_data_con_combat <- run_combat(full_data_con)
full_data_val_combat <- run_combat(full_data_val)



##########
# get mod data
##########

sep_combat_data <- function(combat_data, data_type){
  
  beta_cases <- combat_data[combat_data$batch == 'cases',]
  
  beta_cases <- getModData(beta_cases)
  
  if(data_type == 'controls') {
    beta_other <- combat_data[combat_data$batch == data_type,]
    
    
    beta_other <- subset(beta_other, p53_germline == 'Mut' & 
                              cancer_diagnosis_diagnoses == 'Unaffected')
    
  } else {
    beta_other <- combat_data[combat_data$batch == data_type,]
    beta_other <- beta_other[!beta_other$ids %in% beta_cases$ids,]
    
  }
 
  
  combined_data <- rbind(beta_cases,
                         beta_other)
  
  
  
  return(combined_data)
  
}

mod_data_con_combat <- sep_combat_data(full_data_con_combat, data_type = 'controls' )
mod_data_con <- sep_combat_data(full_data_con, data_type = 'controls')

mod_data_val_combat <- sep_combat_data(full_data_val_combat, data_type = 'valid' )
mod_data_val <- sep_combat_data(full_data_val, data_type = 'valid')


# save normal data and combat data
saveRDS(mod_data_con_combat, paste0(model_data, paste0('/', method, '_', 'mod_data_con_combat.rda')))
saveRDS(mod_data_con, paste0(model_data, paste0('/', method, '_', 'mod_data_con_m.rda')))


saveRDS(mod_data_val_combat, paste0(model_data, paste0('/', method, '_', 'mod_data_val_combat.rda')))
saveRDS(mod_data_val, paste0(model_data, paste0('/', method, '_', 'mod_data_val_m.rda')))





