# this script will get combat data on 2 tech types

# set preprocessing method
method <- 'noob'

# set type of data, beta or m
methyl_type <- 'm'

# source all_functions.R to load libraries and my functions
source('all_functions.R')

# set data directory
data_dir <- '../../Data/'

# get data
if(paste0(methyl_type, '_processed', '.RData') %in% dir(data_dir)) {
  load(paste0(data_dir,methyl_type, '_processed', '.RData'))
} else {
  source('get_data_new.R')
}