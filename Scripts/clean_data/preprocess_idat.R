########## this script will read in raw idat and apply a preprocessing method, plot pca and remove outliers


##########
# load libraries
##########
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
idat_data <- paste0(methyl_data, '/raw_files')
idat_data_controls <- paste0(methyl_data, '/controls')
idat_data_valid <- paste0(methyl_data, '/validation/idat_files')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'raw'

##########
# read in idate for cases, controls, and validation set
##########
rgCases <- read.metharray.exp(idat_data)
rgControls <- read.metharray.exp(idat_data_controls)
rgValid <- read.metharray.exp(idat_data_valid)

##########
# get preprocessing method
##########
betaCases <- preprocessMethod(rgCases, preprocess = method)
betaControls <- preprocessMethod(rgControls, preprocess = method)
betaValid <- preprocessMethod(rgValid, preprocess = method)

###########
# scale and impute
###########

if (method == 'raw') {
  
  betaCases <- scaleImputeDat(dat = betaCases, scale = T)
  betaControls <- scaleImputeDat(dat = betaControls, scale = T)
  betaValid <- scaleImputeDat(dat = betaValid, scale = T)
  
}



