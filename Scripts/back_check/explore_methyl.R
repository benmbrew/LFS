# this script will do a full analysis of methylation data


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
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
# source("https://bioconductor.org/biocLite.R")
# biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
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
rm(rgCases)
betaControls <- preprocessMethod(rgControls, preprocess = method)
rm(rgControls)
betaValid <- preprocessMethod(rgValid, preprocess = method)
rm(rgValid)
# getwd()

###########
# scale and impute
###########
betaCases <- scaleImputeDat(dat = betaCases, scale = T)
betaControls <- scaleImputeDat(dat = betaControls, scale = T)
betaValid <- scaleImputeDat(dat = betaValid, scale = T)
