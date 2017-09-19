
# this script will force read in all methylation data and get column intersection and remove outliers before normalization

##########
# load libraries
##########
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(biovizBase)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)
library(preprocessCore)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
idat_data <- paste0(methyl_data, '/raw_files')
idat_data_con <- paste0(methyl_data, '/controls')
idat_data_valid <- paste0(methyl_data, '/validation/idat_files')
model_data <- paste0(data_folder, '/model_data')


##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'funnorm'


##########
# read in idat 
##########t
rgSet <- read.metharray.exp(idat_data)
rgSetCon <- read.metharray.exp(idat_data_con)
rgSetVal <- read.metharray.exp(idat_data_valid)


  
  
## S4 method for signature 'RGChannelSet,RGChannelSet'
temp_rg <- combineArrays(rgSetCon,
                         outType = "IlluminaHumanMethylation450k",
                         verbose = TRUE)
dim(temp_rg)

# now ass validation set
rgSetFull <- combineArrays(rgSetVal,
                           outType = "IlluminaHumanMethylation450k",
                           verbose = TRUE)

