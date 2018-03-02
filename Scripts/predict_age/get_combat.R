# this script will get combat data on 2 tech types
library(sva)
# set preprocessing method
method <- 'noob'

# set type of data, beta or m
methyl_type <- 'm'

# get the type of data
data_used <- 'old'

# source all_functions.R to load libraries and my functions
source('all_functions.R')

# set data directory
data_dir <- '../../Data/'

# get data
if(paste0(data_used,'_',methyl_type, '_processed', '.RData') %in% dir(data_dir)) {
  load(paste0(data_dir,paste0(data_used,'_',methyl_type, '_processed', '.RData')))
} else {
  source(paste0('get_data_', data_used , '.R'))
}

if(data_used == 'new') {
  # combine data
  full_data <- rbind(data_cases_full,
                     data_controls_full)
  
  rm(data_cases_full, data_controls_full)
} else {
  # remove 'a' and 'b' from columns
  
  # combine data
  full_data <- rbind(data_cases,
                     data_controls_mod,
                     data_valid_mod)
  rm(data_cases, data_controls_mod, data_valid_mod)
}

# run combat on technology
full_data_combat <- run_combat(full_data)

# remove tech variable from full_data_combat
full_data_combat$tech <- NULL

save.image(paste0(data_dir,paste0(data_used,'_',methyl_type, '_final', '.RData')))




