## Script for cleaning IDAT files
# this data could be used to replace methylation from nardine
# this script will be documented heavily as to provide a guide for IDAT files in minfi package to lab.

# initialize libraries
library(minfi)

# Use this https://www.bioconductor.org/help/course-materials/2014/BioC2014/minfi_BioC2014.pdf

# initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/classification_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
idat_data <- paste0(methyl_data, '/raw_files')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

# this function will read in each idat file using read.450k.exp from minifi package and store them in a 
# list

readIDAT <- function(path) {
  
  rgSet <- list()
  directory <- dir(path)
  
  for (folder in directory){
    # read into rgSet
    rgSet[[folder]] <- read.450k.exp(paste(path, folder, sep = '/'))
  }
  
  return(rgSet)
  
}

rgSetList <- readIDAT(idat_data)
# problem here at 5760666027
temp <- read.metharray.exp(paste(idat_data, "5760666027", sep = '/'))
rm(temp)
directory
