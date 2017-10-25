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
model_data <- paste0(data_folder, '/model_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'noob'

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


##########
# read in idate for Controls, controls, and validation set
##########
rgControls <- read.metharray.exp(idat_data)

###########
# remove outliers (previously determined) from rgset before normalization
###########
# 
# rgControls <- remove_outliers(rgSet = rgControls, 
#                               id_map = id_map, 
#                               method = 'funnorm', 
#                               type = 'controls')
# ##########
# get preprocedssing method
##########
betaControls <- preprocessMethod(rgControls, preprocess = method, only_m_values = T)
rm(rgControls)

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
# join data
##########

cg_sites <-  readRDS(paste0(model_data, '/four_fifty_feats.rda'))

intersect_cg_cites <- intersect(cg_sites, colnames(betaControls))


# subset data by colmns of interest and cg_sites
betaControls <- betaControls[, c('ids', 
                                 'sentrix_id',
                                 intersect_cg_cites)]

# inner join
betaControls <- inner_join(clin, betaControls, by = 'ids')

# remove NAs from tm_donor 
betaControls <- betaControls[!is.na(betaControls$tm_donor_),]

# remove duplicates
betaControls <- betaControls[!duplicated(betaControls$tm_donor_),]

# get cg_sites
cg_sites <- colnames(betaControls)[grepl('cg', colnames(betaControls))]

# saveRDS(cg_sites, paste0(model_data, '/four_fifty_feats.rda'))


# subset data by colmns of interest and cg_sites
betaControls <- betaControls[, c('ids', 
                           'p53_germline', 
                           'cancer_diagnosis_diagnoses', 
                           'age_diagnosis',
                           'age_sample_collection',
                           'gender',
                           'sentrix_id',
                           'family_name',
                           cg_sites)]



##########
# remove outliers
##########
betaControls <- removeOutlier(betaControls, 
                              cases = F, 
                              controls = T, 
                              val = F)

##########
# saved unscaled data
##########

saveRDS(betaControls, paste0(model_data, paste0('/', method, '_', 'beta_controls_m.rda')))

