###########################################
# this script will check the relationship between the three controls in 450 and 850k
library(dplyr)
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

##########################
# load 450k and 850k
##########################
load(paste0(idat_data, '/imputed_idat_betas_final_control.RData'))
beta_raw_control <- beta_raw
beta_quan_control <- beta_quan
beta_swan_control <- beta_swan
beta_funnorm_control <- beta_funnorm
rm(beta_raw, beta_quan, beta_swan, beta_funnorm)
load(paste0(idat_data, '/imputed_idat_betas_final.RData'))

##########################
# find the 3 controls that are both in 450 and 850
###########################

# see if the values are the same 

# find overalpping probes between 450k and 850k for those 3 individuals (?)

# for a few probes, plot three individuals

# 3 plots, 450k against 850k, ordering methylation values from smallest to largest

# check to see if the new idat data has enough mut unaffected to handle controls validation
