########################################
# this script will read in and preprocess idat

##########
# load libraries
##########
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

# remove duplicates 

# store ids and sentrix id
ids <- cbind(malkin_ids = betaCases$ids, sentrix_id = betaCases$sentrix_id)

# remove sentrix_id
betaCases$sentrix_id <- NULL

# make betaCases a dat
keys <- colnames(betaCases)[!grepl('cg', colnames(betaCases))]

X <- as.data.table(betaCases)

betaCasesDup <- X[,lapply(.SD,mean),keys]


# remove duplicate ids 


# group by ids and get means
betaCases_dup <- betaCases %>%
  group_by(ids) %>%
  summarise_all(funs(mean))

##########
# join data
##########
temp <- inner_join(betaCases$ids)
length(unique(betaCases$ids))
length(unique(clin$ids))

###########
# scale and impute
###########

if (method == 'raw') {
  
  betaCases <- scaleImputeDat(dat = betaCases, scale = T)
  
}



