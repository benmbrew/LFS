####### Script will give summary stats of model data

##########
# initialize libraries
##########
library(dplyr)
library(sva)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
clin_data <- paste0(data_folder, '/clin_data')


##########
# load data batch data
##########
# read cases
quan_cases <- readRDS(paste0(model_data, '/quan_cases.rda'))

# read batch corrected data for gender
quan_cases_gen <- readRDS(paste0(model_data, '/quan_cases_gen.rda'))

# read batch corrected data for sentrix id and SAM
quan_cases_sen <- readRDS(paste0(model_data, '/quan_cases_sen.rda'))

quan_cases_sam <- readRDS(paste0(model_data, '/quan_cases_sam.rda'))

# read batch corrected data for sentrix id and SAM and gender!
quan_cases_sen_gen <- readRDS(paste0(model_data, '/quan_cases_sen_gen.rda'))

quan_cases_sam_gen <- readRDS(paste0(model_data, '/quan_cases_sam_gen.rda'))

# load features
load(paste0(model_data, '/bh_feat.RData'))


##########
# plot age of diagnosis vs age of sample collection
##########
x_axis <- quan_cases$age_diagnosis
y_axis <- quan_cases$age_sample_collection

plot_dat <- as.data.frame(cbind(x_axis, y_axis))

ggplot(plot_dat, aes(x_axis, y_axis)) + geom_point(alpha= 0.7) + xlab('Age of Diagnosis') + ylab('Age of Sample Collection') +
  ggtitle('Age Correlation') + theme_bw() + geom_abline(intercept = 0, slope = 1)  
  


##########
# look at model data (cases) summary stats 
##########

##########
# difference btw cases and controls
##########

