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
library(sva)
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
method = 'raw'

#################################################################################
# Load id_map, which has ids to merge with methylation - do for cases and controls

##########
# read in cg locations (which is generated later in this script)
##########
cg_locations <- read.csv(paste0(model_data, '/cg_locations.csv'))


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
# getwd()

# save.image('/home/benbrew/Desktop/temp_raw.RData')
# load('/home/benbrew/Desktop/temp_raw.RData')
###########
# scale and impute
###########
betaCases <- scaleImputeDat(dat = betaCases, scale = F)
betaControls <- scaleImputeDat(dat = betaControls, scale = F)
betaValid <- scaleImputeDat(dat = betaValid, scale = F)

###########
# id functions
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
betaValid <- joinData(betaValid, control = F)


##########
# get controls wild type- that is WT non cancer
##########
betaControlsWT <- betaCases[which(betaCases$cancer_diagnosis_diagnoses == 'Unaffected' & 
                                    betaCases$p53_germline == 'WT'),]

##########
# get old controls
##########
betaControlsOld <- betaCases[which(betaCases$cancer_diagnosis_diagnoses == 'Unaffected' &
                                     betaCases$p53_germline == 'Mut'),]


##########
# remove cancers from controls data
##########
betaControls <- removeCancer(betaControls)
betaControlsOld <- removeCancer(betaControlsOld)


##########
# get model dat
##########
betaCases <- getModData(betaCases)

# save.image('/home/benbrew/Desktop/temp_clean_funnorm.RData')
# load('/home/benbrew/Desktop/temp_clean_funnorm.RData')

betaCases <- betaCases[, !grepl('ch', colnames(betaCases))]
betaControls <- betaControls[, !grepl('ch', colnames(betaControls))]
betaControlsWT <- betaControlsWT[, !grepl('ch', colnames(betaControlsWT))]
betaControlsOld <- betaControlsOld[, !grepl('ch', colnames(betaControlsOld))]
betaValid <- betaValid[, !grepl('ch', colnames(betaValid))]

##########
# run combat
##########
# correct for the 9 samples done in montreal
betaCases <- runCombat(betaCases)

# # scale data
# # get row statistics
# colMean <- apply(dat, 2, mean, na.rm=TRUE)
# colSd <- apply(dat, 2, sd, na.rm=TRUE)
# # constantInd <- rowSd==0
# # rowSd[constantInd] <- 1
# rowStats <- list(mean=colMean, sd=colSd)
# 
# # apply normilization
# dat  <- (dat - rowStats$mean) / rowStats$sd
betaCases[, 8:ncol(betaCases)] <- scale(betaCases[, 8:ncol(betaCases)])

##########
# remove outliers
##########
betaControls <- removeOutlier(betaControls, wt = F, val = F)
betaControlsWT <- removeOutlier(betaControlsWT, wt = T, val = T)
betaValid <- removeOutlier(betaControls, wt = F, val = T)


# save each dataset
saveRDS(betaCases, paste0(model_data, '/betaCases', method,'.rda'))
saveRDS(betaControls, paste0(model_data, '/betaControls', method,'.rda'))
saveRDS(betaControlsWT, paste0(model_data, '/betaControlsWT', method,'.rda'))
saveRDS(betaControlsOld, paste0(model_data, '/betaControlsOld', method,'.rda'))
saveRDS(betaValid, paste0(model_data, '/betaValid', method,'.rda'))
# 
# betaCases <- readRDS(paste0(model_data, '/betaCases', method,'.rda'))
# betaControls <- readRDS(paste0(model_data, '/betaControls', method,'.rda'))
# betaValid <- readRDS(paste0(model_data, '/betaValid', method,'.rda'))


# ##########
# # plot pca for all three data types by gender, sentrix_id 
# ##########
# 
# cases, gender
getPCA(pca_data = betaCases,
       column_name = 'gender',
       name = 'betaCases gender',
       gene_start = 8,
       pca1 = 1,
       pca2 = 2)

# cases, sentrix_id
getPCA(pca_data = betaCases,
       column_name = 'sentrix_id',
       name = 'betaCases sentrix_id',
       gene_start = 8,
       pca1 = 1,
       pca2 = 2)

# Controls, gender
getPCA(pca_data = betaControls,
       column_name = 'gender',
       name = 'betaControls gender',
       gene_start = 8,
       pca1 = 1,
       pca2 = 2)

# Controls, sentrix_id
getPCA(pca_data = betaControls,
       column_name = 'sentrix_id',
       name = 'betaControls sentrix_id',
       gene_start = 8,
       pca1 = 1,
       pca2 = 2)


# Controls, gender
getPCA(pca_data = betaControlsWT,
       column_name = 'gender',
       name = 'betaControlsWT gender',
       gene_start = 8,
       pca1 = 1,
       pca2 = 2)

# Controls, sentrix_id
getPCA(pca_data = betaControlsWT,
       column_name = 'sentrix_id',
       name = 'betaControlsWt sentrix_id',
       gene_start = 8,
       pca1 = 1,
       pca2 = 2)

# Valid, gender
getPCA(pca_data = betaValid,
       column_name = 'gender',
       name = 'betaValid gender',
       gene_start = 8,
       pca1 = 1,
       pca2 = 2)

# Valid, sentrix_id
getPCA(pca_data = betaValid,
       column_name = 'sentrix_id',
       name = 'betaValid sentrix_id',
       gene_start = 8,
       pca1 = 1,
       pca2 = 2)


# 
# 
# 
