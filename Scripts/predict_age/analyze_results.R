####################################################################
# This script will analyze the results table from models on idat data.
# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
scripts_folder <- paste0(project_folder, '/Scripts')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
idat_data <- paste0(methyl_data, '/raw_files')
clin_data <- paste0(data_folder, '/clin_data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')
model_data <- paste0(data_folder, '/model_data')


# source model_functions to get functions to run models 
source(paste0(scripts_folder, '/predict_age/model_functions.R'))

# load in results table 
load(paste0(model_data, '/idat_beta_table_results.RData'))

# load in control table 

# this one filters out resid
temp <- beta_idat_results %>% 
  filter(p53_status == 'Mut') %>% 
  filter(type == 'normal') %>% 
  group_by(data, age) %>% 
  summarise(mean_score = mean(score))

temp <- temp[order(temp$mean_score, decreasing = T),]


# this one filters out normal
# this one filters out resid
temp <- beta_idat_results %>% 
  filter(p53_status == 'Mut') %>% 
  filter(type == 'resid') %>% 
  group_by(data, age) %>% 
  summarise(mean_score = mean(score))

temp <- temp[order(temp$mean_score, decreasing = T),]

