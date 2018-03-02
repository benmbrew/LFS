# this script will visualize batch effects with PCA 

library(RColorBrewer)
# set preprocessing method
method <- 'noob'

# old or new data
data_used <- 'old'

# set type of data, beta or m
methyl_type <- 'm'

# source all_functions.R to load libraries and my functions
source('all_functions.R')

# set data directory
data_dir <- '../../Data/'

# load data 
load(paste0(data_dir,paste0(data_used,'_',methyl_type, '_final', '.RData')))

##########
# run PCAs- function creates pca variables on the fly and plots those - age_fac, age_cat, cancer_fac, and tech
##########

# Age as a binary factor
# pca on normal data
get_pca(pca_data = full_data, 
        column_name = 'age_fac',
        pc_x = 1,
        pc_y = 2,
        main_title = 'No combat, Age')

# pca on combat data
get_pca(pca_data = full_data_combat, 
        column_name = 'age_fac',
        pc_x = 1,
        pc_y = 2,
        main_title = 'With combat, Age')

# Age as a multinomial factor
get_pca(pca_data = full_data, 
        column_name = 'age_cat',
        pc_x = 1,
        pc_y = 2,
        main_title = 'No combat, Age category')

# pca on combat data
get_pca(pca_data = full_data_combat, 
        column_name = 'age_cat',
        pc_x = 1,
        pc_y = 2,
        main_title = 'With combat, Age category')


# gender 
# pca on normal data
get_pca(pca_data = full_data, 
        column_name = 'gender',
        pc_x = 1,
        pc_y = 2,
        main_title = 'No combat, gender')

# pca on combat data
get_pca(pca_data = full_data_combat, 
        column_name = 'gender',
        pc_x = 1,
        pc_y = 2,
        main_title = 'With combat, gender')

# cancer as binary
# pca on normal data
get_pca(pca_data = full_data, 
        column_name = 'cancer_fac',
        pc_x = 1,
        pc_y = 2,
        main_title = 'No combat, cancer')

# pca on combat data
get_pca(pca_data = full_data_combat, 
        column_name = 'cancer_fac',
        pc_x = 1,
        pc_y = 2,
        main_title = 'With combat, cancer')


# just tech
# pca on normal data
get_pca(pca_data = full_data, 
        column_name = 'tech',
        pc_x = 1,
        pc_y = 2,
        main_title = 'No combat, tech')

# pca on combat data
get_pca(pca_data = full_data_combat, 
        column_name = 'tech',
        pc_x = 1,
        pc_y = 2,
        main_title = 'With combat, tech')
