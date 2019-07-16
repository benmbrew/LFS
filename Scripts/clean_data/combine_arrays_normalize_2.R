
##########
# load libraries
##########
library(tidyverse)
library(minfi)
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
idat_data <- paste0(methyl_data, '/raw_files')
idat_data_con <- paste0(methyl_data, '/controls')
idat_data_valid <- paste0(methyl_data, '/validation/idat_files')
model_data <- paste0(data_folder, '/model_data')
map_data <- paste0(data_folder, '/methyl_data/validation')


##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'funnorm'

##########
# read in clinical data
##########
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# clean clinical ids
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

##########
# cases batch1
##########
id_map <- read.csv(paste0(methyl_data, '/ids_map.csv'), stringsAsFactors = F)

##########
#cases batch2
##########
id_map_other <- read.csv(paste0(methyl_data, '/batch_2014.csv'), stringsAsFactors = F)
id_map_other$Project <- NULL

##########
# combine id_map and id_map_other
##########
id_map <- rbind(id_map, id_map_other)
rm(id_map_other)

##########
# clean idmap
##########
id_map <- cleanIdMap(id_map)


##########
# valid
##########
id_map_valid <- read.csv(paste0(map_data, '/id_map_validation.csv'), stringsAsFactors = F)

# homogenize valid map data with cases and controls
id_map_valid <- id_map_valid[, c('Sample.ID', 'Sample.Group', 'Sentrix.Barcode', 'Sample.Section',
                                 'Project', 'Pool_ID', 'Sample_Well')]

# sub out '.' for '_'
colnames(id_map_valid) <- gsub('.', '_', colnames(id_map_valid), fixed = T)

# change 'Sample_ID' to 'Sample_Name' and 'Sentrix_Barcode' to 'Sentrix_ID'
colnames(id_map_valid)[1] <- 'Sample_Name'
colnames(id_map_valid)[3] <- 'Sentrix_ID'
colnames(id_map_valid)[4] <- 'Sentrix_Position'
colnames(id_map_valid)[5] <- 'Sample_Plate'

##########
# clean idmap
##########
id_map_valid <- cleanIdMap(id_map_valid)


##########
# read in idat 
##########t
rgSetCase <- read.metharray.exp(idat_data)
rgSetVal <- read.metharray.exp(idat_data_valid)

# remove outliers 
rgSetCase <- remove_outliers(rgSet = rgSetCase, method = 'noob', id_map = id_map, type = 'cases')
rgSetVal <- remove_outliers(rgSet = rgSetVal, method = 'noob', id_map = id_map_valid, type = 'valid')

case_val <- combineArrays(rgSetCase, 
                          rgSetVal, 
                          outType = "IlluminaHumanMethylation450k")

rm(rgSetCase, rgSetVal)
# 
# save.image('~/Desktop/temp_cases_valid.RData')
load('~/Desktop/temp_cases_valid.RData')


##########
# get preprocedssing method
##########
# use noob on beta then conver to matrix of
cases_val_beta <- preprocessMethod(case_val, preprocess = method, only_m_values = T)

# saveRDS(cases_val_beta, '~/Desktop/cases_val_m.rda')
# cases_val_beta <- readRDS('~/Desktop/cases_val_m.rda')


cases_val_beta_full <- process_rg_set(beta = cases_val_beta, 
                                     id_map_1 = id_map, 
                                     id_map_2 = id_map_valid, 
                                     clinical_dat = clin,
                                     controls = F)
##########
# save version of data to explore batches on pca
##########
saveRDS(cases_val_beta_full, paste0(model_data, paste0('/', method, '_', 'cases_val_beta_m.rda')))




