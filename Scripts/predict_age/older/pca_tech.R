# this script will visualize batch effects with PCA 

library(RColorBrewer)
# set preprocessing method
method <- 'noob'

# old or new data
data_used <- 'new'

# set type of data, beta or m
methyl_type <- 'beta'

# source all_functions.R to load libraries and my functions
source('all_functions.R')

# set data directory
data_dir <- '../../Data/'

# load data 
load(paste0(data_dir,paste0(data_used,'_',methyl_type, '_final_beta_first_last_gen', '.RData')))

##########
# run PCAs- function creates pca variables on the fly and plots those - age_fac, age_cat, cancer_fac, and tech
##########
full_data_first <- full_data_first[full_data_first$tm_donor_ != '3955',]
full_data_first_combat <- full_data_first_combat[full_data_first_combat$tm_donor_ != '3955',]

full_data_last <- full_data_last[full_data_last$tm_donor_ != '3955',]
full_data_last_combat <- full_data_last_combat[full_data_last_combat$tm_donor_ != '3955',]

full_data_last <- full_data_last[!duplicated(full_data_last$tm_donor_, fromLast = TRUE),]
full_data_first <- full_data_first[!duplicated(full_data_first$tm_donor_, fromLast = FALSE),]
# Age as a multinomial factor

get_pca(pca_data = sample_tm, 
        age_cutoff = 72,
        column_name = 'age_fac',
        pc_x = 1,
        pc_y = 2,
        main_title = 'PC 1 and 2 for 72 months')


get_pca(pca_data = full_data_first,
        age_cutoff = 100,
        column_name = 'gender',
        pc_x = 1,
        pc_y = 2,
        main_title = 'Gender')


# gender 
# pca on normal data
get_pca(pca_data = full_data_first, 
        column_name = 'gender',
        age_cutoff = 72,
        pc_x = 1,
        pc_y = 2,
        main_title = 'PC 1 and PC 2 by gender')

# pca on combat data
get_pca(pca_data = full_data_last, 
        column_name = 'gender',
        age_cutoff = 72,
        pc_x = 1,
        pc_y = 2,
        main_title = 'PC 1 and 2 with combat')

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
