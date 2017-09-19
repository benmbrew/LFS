########################################
# this script will read in and preprocess idat for validation 

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
idat_data <- paste0(methyl_data, '/validation/idat_files')
map_data <- paste0(data_folder, '/methyl_data/validation')
model_data <- paste0(data_folder, '/model_data')


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

# clean clinical idss
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

##########
# valid
##########
id_map <- read.csv(paste0(map_data, '/id_map_validation.csv'), stringsAsFactors = F)

# homogenize valid map data with cases and controls
id_map <- id_map[, c('Sample.ID', 'Sample.Group', 'Sentrix.Barcode', 'Sample.Section',
                     'Project', 'Pool_ID', 'Sample_Well')]

# sub out '.' for '_'
colnames(id_map) <- gsub('.', '_', colnames(id_map), fixed = T)

# change 'Sample_ID' to 'Sample_Name' and 'Sentrix_Barcode' to 'Sentrix_ID'
colnames(id_map)[1] <- 'Sample_Name'
colnames(id_map)[3] <- 'Sentrix_ID'
colnames(id_map)[4] <- 'Sentrix_Position'
colnames(id_map)[5] <- 'Sample_Plate'

##########
# clean idmap
##########
id_map <- cleanIdMap(id_map)

##########
# read in idate for Valid, controls, and validation set
##########
rgValid <- read.metharray.exp(idat_data)

###########
# remove outliers (previously determined) from rgset before normalization
###########

rgValid <- remove_outliers(rgSet = rgValid, 
                           id_map = id_map, 
                           method = 'funnorm', 
                           type = 'valid')

##########
# get preprocedssing method
##########
betaValid <- preprocessMethod(rgValid, preprocess = method, only_m_values = T)
rm(rgValid)


###########
# id functions
###########
# Valid
betaValid <- findIds(betaValid, id_map = id_map)

# get id name (only Valid)
betaValid <- getIdName(betaValid)

# clean ids
betaValid <- cleanIds(betaValid)

# remove 'ch' from column names
betaValid <- betaValid[, !grepl('ch', colnames(betaValid))]

##########
# join data
##########
cg_sites <-  readRDS(paste0(model_data, '/four_fifty_feats.rda'))

intersect_cg_cites <- intersect(cg_sites, colnames(betaValid))


# subset data by colmns of interest and cg_sites
betaValid <- betaValid[, c('ids', 
                          'sentrix_id',
                          intersect_cg_cites)]


# inner join
betaValid <- inner_join(clin, betaValid, by = 'ids')

# order betaValid by tm_donor_ and age of sample collection
betaValid <- betaValid[order(betaValid$tm_donor_, betaValid$age_sample_collection),]

# remove NAs from tm_donor 
betaValid <- betaValid[!is.na(betaValid$tm_donor_),]

# remove duplicates
betaValid <- betaValid[!duplicated(betaValid$tm_donor_),]

##########
# get data in format for saving
##########

# get cg_sites
cg_sites <- colnames(betaValid)[grepl('cg', colnames(betaValid))]

# subset data by colmns of interest and cg_sites
betaValid <- betaValid[, c('ids', 
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
saveRDS(betaValid, paste0(model_data, paste0('/', method, '_', 'valid_batch_m_sub.rda')))

##########
# remove NA
##########
betaValid <- removeNA(betaValid, probe_start = 8) #450168

##########
# remove outliers
##########
betaValid <- removeOutlier(betaValid, 
                           cases = F, 
                           controls = F, 
                           val = T)

##########
# remove inf
##########
if(method == 'raw') {
  betaValid <- removeInf(betaValid, probe_start = 8)
  
}

##########
# saved unscaled data
##########

saveRDS(betaValid, paste0(model_data, paste0('/', method, '_', 'valid_new_m.rda')))

