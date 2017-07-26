########################################
# this script will read in and preprocess idat for controls

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
idat_data <- paste0(methyl_data, '/controls')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'raw'

##########
# read in idate for Controls, controls, and validation set
##########
rgControls <- read.metharray.exp(idat_data)

##########
# get preprocedssing method
##########
betaControls <- preprocessMethod(rgControls, preprocess = method)
rm(rgControls)

##########
# read in clinical data
##########
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# clean clinical idss
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

##########
# Controls batch1
##########
id_map <- read.csv(paste0(methyl_data, '/ids_map_controls.csv'), stringsAsFactors = F)


##########
# clean idmap
##########
id_map <- cleanIdMap(id_map)

###########
# id functions
###########
# Controls
betaControls <- findIds(betaControls, id_map = id_map)

# get id name (only Controls)
betaControls <- getIdName(betaControls)

# clean ids
betaControls <- cleanIds(betaControls)

# remove 'ch' from column names
betaControls <- betaControls[, !grepl('ch', colnames(betaControls))]


##########
# read in cases and get features
##########
betaCases <- readRDS(paste0(methyl_data, '/betaCasesBatch.rda'))

# get intersecting features
intersect_feat <- Reduce(intersect, list(colnames(betaCases)[7:ncol(betaCases)],
                                         colnames(betaControls)[1:ncol(betaControls)]))

# subset betaCases and save
betaCases <- betaCases[, c('ids', 
                           'p53_germline', 
                           'cancer_diagnosis_diagnoses', 
                           'age_diagnosis',
                           'age_sample_collection',
                           'gender',
                           intersect_feat)]

saveRDS(betaCases, paste0(methyl_data, '/betaCasesBatch.rda'))
rm(betaCases)

# subset betaControls 
betaControls <- betaControls[, c('ids', 
                                 'sentrix_id', 
                                  intersect_feat)]

# save.image('/home/benbrew/Desktop/raw_Controls_temp.RData')
# load('/home/benbrew/Desktop/raw_Controls_temp.RData')
##########
# join data
##########

# inner join
betaControls <- inner_join(clin, betaControls, by = 'ids')

# remove NAs from tm_donor 
betaControls <- betaControls[!is.na(betaControls$tm_donor_),]

# remove duplicates
betaControls <- betaControls[!duplicated(betaControls$tm_donor_),]

##########
# get data in format for saving
##########

# get cg_sites
cg_sites <- colnames(betaControls)[grepl('cg', colnames(betaControls))]

# subset data by colmns of interest and cg_sites
betaControls <- betaControls[, c('ids', 
                                  'p53_germline', 
                                  'cancer_diagnosis_diagnoses', 
                                  'age_diagnosis',
                                  'age_sample_collection',
                                  'gender',
                                  cg_sites)]


##########
# save version of data to explore batches on pca
##########
saveRDS(betaControls, paste0(methyl_data, '/betaControlsBatch.rda'))

