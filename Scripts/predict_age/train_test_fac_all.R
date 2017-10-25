########################################
# this script will read in and preprocess idat

##########
# load libraries
##########
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(biovizBase)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)
library(preprocessCore)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

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

##########
# source all_functions.R script
##########

##########
# fixed variables
##########
method = 'raw'
k_folds = 5

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


load('~/Desktop/temp_full_pipeline.RData')
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

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


##########
# read in idate for cases, controls, and validation set
########## 
rg_cases <- rgCases
rg_controls <- rgControls
rg_valid <- rgValid

full_pipeline <- function(rg_cases, rg_controls, rg_valid) {
  
  # get vector of random folds
  fold_vec <- sample(1:k_folds, dim(rg_cases)[2], replace = T)
  
  # get train and test index
  train_index <- 
  
}

##########
# get preprocedssing method
##########
betaControls <- preprocessMethod(rgControls, preprocess = method, only_m_values = T)
rm(rgControls)


##########
# get preprocedssing method
##########
# use noob on beta then conver to m
betaCases <- preprocessMethod(rgCases, preprocess = method, only_m_values = T)
# remove rgset
rm(rgCases)

##########
# get preprocedssing method
##########
betaValid <- preprocessMethod(rgValid, preprocess = method, only_m_values = T)
rm(rgValid)



