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
idat_data <- paste0(methyl_data, '/raw_files')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'raw'

##########
# read in idate for cases, controls, and validation set
##########
rgCases <- read.metharray.exp(idat_data)

##########
# get preprocedssing method
##########
betaCases <- preprocessMethod(rgCases, preprocess = method)
rm(rgCases)

##########
# read in clinical data
##########
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# clean clinical idss
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

###########
# id functions
###########
# cases
betaCases <- findIds(betaCases, id_map = id_map)

# get id name (only cases)
betaCases <- getIdName(betaCases)

# clean ids
betaCases <- cleanIds(betaCases)

# remove 'ch' from column names
betaCases <- betaCases[, !grepl('ch', colnames(betaCases))]

##########
# join data
##########

# inner join
betaCases <- inner_join(clin, betaCases, by = 'ids')

# remove NAs from tm_donor 
betaCases <- betaCases[!is.na(betaCases$tm_donor_),]

# remove duplicates
betaCases <- betaCases[!duplicated(betaCases$tm_donor_),]

##########
# get data in format for saving
##########

# get cg_sites
cg_sites <- colnames(betaCases)[grepl('cg', colnames(betaCases))]

# subset data by colmns of interest and cg_sites
betaCases <- betaCases[, c('ids', 
                           'p53_germline', 
                           'cancer_diagnosis_diagnoses', 
                           'age_diagnosis',
                           'age_sample_collection',
                           'gender',
                           'sentrix_id',
                           cg_sites)]


##########
# save version of data to explore batches on pca
##########
# saveRDS(betaCases, paste0(methyl_data, '/betaCasesBatch.rda'))
# betaCases <- readRDS(paste0(methyl_data, '/betaCasesBatch.rda'))


##########
# remove NA
##########
betaCases <- removeNA(betaCases, probe_start = 8)

##########
# remove outliers
##########
betaCases <- removeOutlier(betaCases, 
                           cases = T, 
                           controls = F, 
                           val =F)


##########
# scale data
##########
betaCases <- scaleData(betaCases)

##########
# batch correction
########## 


