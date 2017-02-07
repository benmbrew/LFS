##########
# this script will look correlate random features to bh features and get union of all bh to remove from data for 
# radom models

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
model_data <- paste0(data_folder, '/model_data')
data_folder <- paste0(project_folder, '/Data')
scripts_folder <- paste0(project_folder, '/Scripts')

##########
# Read in cases and bumphunter features
##########
# bh_features
load(paste0(model_data, '/bh_features.RData'))

# cases - raw, swan, quan, funnorm, and clin
load(paste0(model_data, '/model_data_cases.RData'))

##########
# function taht gets data you want and removes other
##########
removeDat <- function(keep)
{
  data_set <- c('raw', 'quan', 'swan', 'funnorm')
  
  remove <- data_set[data_set != keep]
  
  # remove unwated
  rm(list=ls(pattern=remove[[1]]))
  rm(list=ls(pattern=remove[[2]]))
  rm(list=ls(pattern=remove[[3]]))
  
}

removeDat(keep = 'raw')


##########
# function that subsets data by bh and correlates to random
##########
beta_data <- beta_funnorm
bh_features <- 
bhRandFinder <- function(beta_data, bh_features)



##########
# get union of all bh features
##########

