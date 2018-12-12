# this scrip will predict on data that was linear transformed and age regressed out.

# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  
methyl_type = 'beta'
combat = TRUE
gender = TRUE
tech = FALSE 
how_many_seeds = 100
how_many_folds = 5
use_offset = TRUE
age_removal = TRUE
transform = TRUE


