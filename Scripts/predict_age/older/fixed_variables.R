# set fixed variables
size = 'full'
model_type = 'rf'
null_450= TRUE
null_450_all = FALSE
use_p53 = FALSE
gender = TRUE
use_cancer = FALSE
method = 'noob'
include_under_6 = FALSE
combat = 'combat_1'
alpha_val = 0.1
train_alpha = FALSE
train_lambda = FALSE
train_cutoff = TRUE
mean_lambda = FALSE
which_methyl = 'beta'
beta_thresh = 0.1
age_cutoff =  72
tech = FALSE
standardize = FALSE
boot_num = 50


if(null_450 & !null_450_all){
  use_null_450 <- 'used_null_450_mut'
} else  if(!null_450 & null_450_all){
  use_null_450 <- 'used_null_450_all'
} else {
  use_null_450 <- 'no_null_450'
  
}

if(include_under_6){
  used_under_6 <- 'used_under_6'
} else {
  used_under_6 <- 'no_under_6'
}
if(use_cancer){
  removed_cancer <- 'removed_cancer'
} else {
  removed_cancer <- 'no_removed_cancer'
}

if(use_p53){
  used_p53 <- 'used_p53'
} else {
  used_p53 <- 'no_p53'
}

if(use_cancer){
  removed_cancer <- 'removed_cancer'
} else {
  removed_cancer <- 'no_removed_cancer'
}
if(standardize){
  standardize_data <- 'standardized'
} else {
  standardize_data <- 'not_standardized'
}

if(train_lambda){
  is_lambda <- 'lambda_test'
} else {
  is_lambda <- 'lambda_train'
}

if(train_alpha){
  is_alpha <- 'alpha_test'
} else {
  is_alpha <- 'alpha_train'
}


if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}
how_many_seeds = 50
how_many_folds = 5


num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds

# create range for random sampling for cross validation
seed_range <- c(1:how_many_seeds)
