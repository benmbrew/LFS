##########
# initialize libraries
##########
library(tidyverse)



registerDoParallel(1)

##########
# initialize folders
##########

home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
results_data <- paste0(data_folder, '/results_data')


##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'noob'
k = 5
max_age = 850



##########
# read in class results
##########

# random increase all ages
rand_increase <- readRDS(paste0(results_data, paste0('/', method, '_', max_age, '_','rand_increase.rda')))

# with bumphunter features all ages 
class_results <- readRDS(paste0(results_data, '/', method,'_','full_results_class_collapsed.rda'))

# sort each one by increasing number of features
rand_increase <- rand_increase[order(rand_increase$feature_num, decreasing = F),]
class_results <- class_results[order(class_results$feature_num, decreasing = F),]
