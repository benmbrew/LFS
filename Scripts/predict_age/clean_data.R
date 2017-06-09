######
# This script will act as the "model pipeline" grabbing bumphunter features on training set only

##########
# initialize libraries
##########
library(minfi)
library(bumphunter)
library(dplyr)
library(dgof)
library(graphics)
library(stringr)
library(impute)
library(data.table)
library(GenomicRanges)
library(biovizBase)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)
library(preprocessCore)


##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
map_data <- paste0(data_folder, '/methyl_data/validation')
idat_data <- paste0(methyl_data, '/raw_files')
idat_data_controls <- paste0(methyl_data, '/controls')
idat_data_valid <- paste0(methyl_data, '/validation/idat_files')
clin_data <- paste0(data_folder, '/clin_data')
model_data <- paste0(data_folder, '/model_data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'quan'

#################################################################################
# Load id_map, which has ids to merge with methylation - do for cases and controls

##########
# read in cg locations (which is generated later in this script)
##########

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
# controls
##########
id_map_control <- read.csv(paste0(methyl_data, '/ids_map_controls.csv'), stringsAsFactors = F)

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
# combine id_map and id_map_other
##########
id_map <- rbind(id_map, id_map_other)

##########
# clean idmap
##########
id_map <- cleanIdMap(id_map)

id_map_control <- cleanIdMap(id_map_control)
id_map_valid <- cleanIdMap(id_map_valid)

##########
# read in idate for cases, controls, and validation set
##########
rgCases <- read.metharray.exp(idat_data)
rgControls <- read.metharray.exp(idat_data_controls)
rgValid <- read.metharray.exp(idat_data_valid)

##########
# get preprocessing method
##########
betaCases <- preprocessMethod(rgCases, preprocess = method)
betaControls <- preprocessMethod(rgControls, preprocess = method)
betaValid <- preprocessMethod(rgValid, preprocess = method)

###########
# fix id functions
###########
# cases
betaCases <- findIds(betaCases, id_map = id_map)

# get id name (only cases)
betaCases <- getIdName(betaCases)

# controls and validation
betaControls <- findIds(betaControls, id_map_control)
betaValid <- findIds(betaValid, id_map_valid)

# remove identifier from betaValid and betaControls
# (cases has sentrix_id and ids) (controls sentrix_id, ids, identifier)
betaControls$identifier <- betaValid$identifier <- NULL

# clean ids
betaCases <- cleanIds(betaCases)
betaControls <- cleanIds(betaControls)
betaValid <- cleanIds(betaValid)

# # get cg locations
# cg_locations <- getIds()
# 
# # write.csv(cg_locations, paste0(model_data, '/cg_locations.csv'))
# 

###########
# join data 
###########

# join data to clinical 
betaCases <- joinData(betaCases, control = F)
betaControls <- joinData(betaControls, control = T)
betaValid <- joinData(betaValid, control = T)

##########
# get controls wild type- that is WT non cancer
##########
betaControlsWT <- betaCases[which(betaCases$cancer_diagnosis_diagnoses == 'Unaffected' & betaCases$p53_germline == 'WT'),]

##########
# remove cancers from controls data
##########
betaControls <- removeCancer(betaControls)

##########
# get model dat
##########
betaCases <- getModData(betaCases)



# # save each dataset
# saveRDS(betaCases, paste0(model_data, '/betaCases', method,'.rda'))
# saveRDS(betaControls, paste0(model_data, '/betaControls', method,'.rda'))
# saveRDS(betaValid, paste0(model_data, '/betaValid', method,'.rda'))
# 
# betaCases <- readRDS(paste0(model_data, '/betaCases', method,'.rda'))
# betaControls <- readRDS(paste0(model_data, '/betaControls', method,'.rda'))
# betaValid <- readRDS(paste0(model_data, '/betaValid', method,'.rda'))




