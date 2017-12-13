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
# source all_functions.R script
##########
source('LFS/Scripts/predict_age/all_functions.R')

##########
# fixed variables
##########
method = 'raw'

##########
# read in clinical data
##########
clin <- read.csv('LFS/Data/clin_data/clinical_two.csv', stringsAsFactors = F)

# clean clinical ids
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

##########
# cases 
##########

# cases batch1
id_map_tor <- read.csv('LFS/Data/methyl_data/cases_toronto/SampleSheet.csv', stringsAsFactors = F)

#cases batch2
id_map_mon <- read.csv('LFS/Data/methyl_data/cases_montreal/SampleSheet.csv', stringsAsFactors = F)
id_map_mon$Project <- NULL

# combine id_map and id_map_other
id_map_cases <- rbind(id_map_tor, id_map_mon)
rm(id_map_tor, id_map_mon)

# clean id map
id_map_cases <- cleanIdMap(id_map_cases)
##########
# read in idate for cases, controls, and validation set
##########t
rgCases <- read.metharray.exp('LFS/Data/methyl_data/cases_montreal', recursive = T)


##########
# get preprocedssing method
##########
# use noob on beta then conver to m
betaCases <- preprocessMethod(rgCases, preprocess = method)
# remove rgset
rm(rgCases)

# save.image('~/Desktop/temp_cases_noob_m.RData')
# load('~/Desktop/temp_cases_noob_m.RData')
full_data <- readRDS('LFS/Data/model_data/raw_full_mod.rda')
full_data$family_name.1 <- NULL
full_data <- full_data[full_data$type != 'controls_wt_450k',]

cases <- subset(full_data, type == 'cases_450k' & 
                  p53_germline == 'Mut')
cases <- cases[!is.na(cases$age_sample_collection),]
controls <- subset(full_data, type == 'controls_850k' & 
                     p53_germline == 'Mut')
controls <- controls[!duplicated(controls$ids),]
controls_old <- subset(full_data, type == 'controls_450k' & 
                     p53_germline == 'Mut')
valid <- subset(full_data, type == 'valid_850k' & 
                     p53_germline == 'Mut')

full_data <- rbind(cases, controls, valid)

saveRDS(full_data, '~/Desktop/full_data_raw.rda')

length(which(duplicated(full_data$ids)))
###########
# id functions
###########
# cases
betaCases <- findIds(betaCases, id_map = id_map_cases)

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

# saveRDS(cg_sites, paste0(model_data, '/four_fifty_feats.rda'))


# subset data by colmns of interest and cg_sites
betaCases <- betaCases[, c('ids', 
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
betaCases <- removeOutlier(betaCases, 
                           cases = T, 
                           controls = F, 
                           val =F)


##########
# get old controls 
##########
betaControlsOld <- subset(betaCases, p53_germline == 'Mut' &
                            cancer_diagnosis_diagnoses == 'Unaffected')

##########
# saved unscaled data
##########

saveRDS(betaControlsOld, paste0(model_data, paste0('/', method, '_', 'beta_controls_old_m.rda')))

saveRDS(betaCases, paste0(model_data, paste0('/', method, '_', 'beta_cases_m.rda')))
















