###############################################
# This script will run a PCA on the LFS methylation data and plot age of diagnosis.
# existing by gene methylation data we have. 
# this is the 6th step in the pipeline
library(dplyr)
library(stringr)
library(impute)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/Analyze')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

# load imputed methylation
setwd(data_folder)
load(paste0(data_folder, '/methyl_lsa.RData'))

# load functions
source(paste0(project_folder, '/Code/Functions/lsaImputation.R'))

################################################################
# Read in methyl and clinical data and join by ids
################################################################

# remove 'A' and '_' in methylation names
methyl_impute_raw$id <- gsub('_', '', methyl_impute_raw$id)
methyl_impute_raw$id <- gsub('A', '', methyl_impute_raw$id)

# Read in data (clinical or clinical_two)
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = TRUE)

# make clin id a factor so it joins with methylation data
clin$id <- as.factor(clin$blood_dna_malkin_lab_)

# inner_join clin
full_data <- inner_join(clin, methyl_impute_raw,
                        by = 'id')

# Save data to be used later
write.csv(full_data, paste0(data_folder, '/full_data.csv'))

