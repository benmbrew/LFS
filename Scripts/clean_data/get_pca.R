## this script will read in batch data from get_cases, get_controls, or get_valid and explore potential batches and outliers

############
# initialize libraries
##########
library(dplyr)
library(sva)
library(impute)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
clin_data <- paste0(data_folder, '/clin_data')

##########
# read in batch data
##########
betaCases <- readRDS(paste0(methyl_data, '/betaCasesBatch.rda'))
betaControls <- readRDS(paste0(methyl_data, '/betaControlsBatch.rda'))
betaValid <- readRDS(paste0(methyl_data, '/betaValidBatch.rda'))

##########
# run pca on data
##########
colnames(betaCases)[1:10]
