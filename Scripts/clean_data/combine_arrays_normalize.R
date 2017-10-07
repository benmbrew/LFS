
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
# Controls batch1
##########
id_map_con <- read.csv(paste0(methyl_data, '/ids_map_controls.csv'), stringsAsFactors = F)

##########
# clean idmap
##########
id_map_con <- cleanIdMap(id_map_con)

##########
# read in idat 
##########t
rgSetCase <- read.metharray.exp(idat_data)
rgSetCon <- read.metharray.exp(idat_data_con)

# remove outliers 
rgSetCase <- remove_outliers(rgSet = rgSetCase, method = 'noob', id_map = id_map, type = 'cases')
rgSetCon <- remove_outliers(rgSet = rgSetCon, method = 'noob', id_map = id_map_con, type = 'controls')


# combine array - (1) cases with con and (2) cases with valid
case_con <- combineArrays(rgSetCase, 
                          rgSetCon, 
                          outType = "IlluminaHumanMethylation450k")


rm(rgSetCase, rgSetCon)

##########
# get preprocedssing method
##########
# use noob on beta then conver to matrix of
cases_con_beta <- preprocessMethod(case_con, preprocess = method, only_m_values = T)

cases_con_beta <- readRDS('~/Desktop/cases_con.rda')
dim(cases_con_beta)

cases_con_beta_full <- process_rg_set(beta = cases_con_beta, 
                                     id_map_1 = id_map, 
                                     id_map_2 = id_map_con, 
                                     clinical_dat = clin,
                                     controls = T)

length(which(grepl('^200', cases_con_beta_full$sentrix_id)))
##########
# save version of data to explore batches on pca
##########
saveRDS(cases_con_beta_full, paste0(model_data, paste0('/', method, '_', 'cases_con_beta.rda')))



