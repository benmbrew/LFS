####### This script will merge the existing clinical data we have with the 
# existing by gene methylation data we have. 

# Initialize folders
home_folder <- '/home/benbrew/Documents'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')

################################################################
# Read in methyl and clinical data and join by ids
################################################################

# Read in data 
methyl <- read.csv(paste0(methyl_data, '/methyl.csv'))
clin <- read.csv(paste0(clin_data, '/clinical.csv'))

# left_join 

