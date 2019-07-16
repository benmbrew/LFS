# this script will read in Valli's data and will normalize training and test set seperately. 

# load functions 
source('../predict_age/all_functions.R')

# read in data 
dat1 <- readRDS('~/Desktop/train_test.rds')
dat2 <- readRDS('~/Desktop/train_test2.rds')

# read in all methylation data
##########
# get base directory for 4 batch
##########
path_to_cases_tor <- '../../Data/methyl_data/cases_toronto'
path_to_cases_mon <- '../../Data/methyl_data/cases_montreal'
path_to_controls <- '../../Data/methyl_data/controls'
path_to_valid <- '../../Data/methyl_data/validation'

##########
# read in meth array - Data/methyl_data/cases_toronto, cases_montreal, controls, validation
##########

# cases 
rgCasesT <- read.metharray.exp(path_to_cases_tor, recursive = T)
rgCasesM <- read.metharray.exp(path_to_cases_mon, recursive = T)

# combine cases arrays 
rgCases <- combineArrays(rgCasesT, rgCasesM)
rm(rgCasesT, rgCasesM)

# controls
rgControls <- read.metharray.exp(path_to_controls, recursive = T)

rgValid <- read.metharray.exp(path_to_valid, recursive = T)

##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
##########
ratio_set <- readRDS('../../Data/model_data/raw_ratio_set.rda')

# get granges object
g_ranges <- as.data.frame(getLocations(ratio_set))

# get probes from rownames
g_ranges$probe <- rownames(g_ranges)

# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]

##########
# read in clinical data
##########
clin <- read.csv('../../Data/clin_data/new_clin.csv', stringsAsFactors = F)

# clean clinical ids
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab)

##########
# cases 
##########

# cases batch1
id_map_tor <- read.csv(paste0(path_to_cases_tor, '/SampleSheet.csv'), stringsAsFactors = F)

#cases batch2
id_map_mon <- read.csv(paste0(path_to_cases_mon, '/SampleSheet.csv'), stringsAsFactors = F)
id_map_mon$Project <- NULL

# combine id_map and id_map_other
id_map_cases <- rbind(id_map_tor, id_map_mon)
rm(id_map_tor, id_map_mon)

# clean id map
id_map_cases <- cleanIdMap(id_map_cases)


##########
# Controls batch1
##########
id_map_con <- read.csv(paste0(path_to_controls, '/SampleSheet.csv'), stringsAsFactors = F)

# clean idmap
id_map_con <- cleanIdMap(id_map_con)

##########
# valid
##########
id_map_val <- read.csv(paste0(path_to_valid, '/SampleSheet.csv'))

# homogenize valid map data with cases and controls
id_map_val <- id_map_val[, c('Sample.ID', 'Sample.Group', 'Sentrix.Barcode', 'Sample.Section',
                             'Project', 'Pool_ID', 'Sample_Well')]

# sub out '.' for '_'
colnames(id_map_val) <- gsub('.', '_', colnames(id_map_val), fixed = T)

# change 'Sample_ID' to 'Sample_Name' and 'Sentrix_Barcode' to 'Sentrix_ID'
colnames(id_map_val)[1] <- 'Sample_Name'
colnames(id_map_val)[3] <- 'Sentrix_ID'
colnames(id_map_val)[4] <- 'Sentrix_Position'
colnames(id_map_val)[5] <- 'Sample_Plate'

# clean idmap
id_map_val <- cleanIdMap(id_map_val)
