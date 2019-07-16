
source('helper_functions.R')
# create fixed objects to model and pipeline inputs and saving  
data_type = 'beta'

if(data_type == 'beta'){
  beta_thresh = 0.05
} else {
  beta_thresh = 0.5
}


gender = FALSE
tech = FALSE
how_many_seeds = 20
how_many_folds = 5
age_cutoff = 72
model_type = 'enet'

if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}

num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds

# create range for random sampling for cross validation
seed_range <- c(1:how_many_seeds)


# # save data
final_dat <- readRDS(paste0('transform_data_cv/results/', data_type, '_',
                            num_seeds, '_', k_folds, '_', is_gen, '_',model_type,'_',age_cutoff,'.rda'))


