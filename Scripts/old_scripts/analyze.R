########################
# this script will analyze table results for rf and enet
library(dplyr)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
mod_data_folder <- paste0(data_folder, '/model_data')
results_folder <- paste0(project_folder, '/Results')
quan_folder <- paste0(results_folder, '/quan')
rand_folder <- paste0(results_folder, '/rand')
scripts_folder <- paste0(project_folder, '/Scripts')

##########
# read data
##########

# enet
enet_bat <- readRDS(paste0(quan_folder, '/enet_bat.rda'))
enet_cases <- readRDS(paste0(quan_folder, '/enet_cases.rda'))
enet_sam <- readRDS(paste0(quan_folder, '/enet_sam.rda'))
enet_sen <- readRDS(paste0(quan_folder, '/enet_sen.rda'))
enet_sam_bat <- readRDS(paste0(quan_folder, '/enet_sam_bat.rda'))
enet_sen_bat <- readRDS(paste0(quan_folder, '/enet_sen_bat.rda'))

# rf
rf_bat <- readRDS(paste0(quan_folder, '/rf_bat.rda'))
rf_cases <- readRDS(paste0(quan_folder, '/rf_cases.rda'))
rf_sam <- readRDS(paste0(quan_folder, '/rf_sam.rda'))
rf_sen <- readRDS(paste0(quan_folder, '/rf_sen.rda'))
rf_sam_bat <- readRDS(paste0(quan_folder, '/rf_sam_bat.rda'))
rf_sen_bat <- readRDS(paste0(quan_folder, '/rf_sen_bat.rda'))

##########
# combine data
##########

result_table <- rbind(enet_bat,
                      enet_cases,
                      enet_sam,
                      enet_sen,
                      enet_sam_bat,
                      enet_sen_bat,
                      rf_bat,
                      rf_cases,
                      rf_sam,
                      rf_sen,
                      rf_sam_bat,
                      rf_sen_bat)

rm(enet_bat,
   enet_cases,
   enet_sam,
   enet_sen,
   enet_sam_bat,
   enet_sen_bat,
   rf_bat,
   rf_cases,
   rf_sam,
   rf_sen,
   rf_sam_bat,
   rf_sen_bat)

##########
# order by score
##########
result_table <- result_table[order(result_table$score),]
