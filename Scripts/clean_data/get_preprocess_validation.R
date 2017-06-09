###### Script for cleaning IDAT files for validation set

##########
# initialize libraries
##########
library(minfi)
library(dplyr)


# Use this https://www.bioconductor.org/help/course-materials/2014/BioC2014/minfi_BioC2014.pdf

##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <-  paste0(data_folder, '/methyl_data')
map_data <- paste0(data_folder, '/methyl_data/validation')
idat_data <- paste0(methyl_data, '/validation/idat_files')

##########
# Load id_map
##########

# batch1 
id_map <- read.csv(paste0(map_data, '/id_map_validation.csv'), stringsAsFactors = F)

##########
# clean IdMap
##########
cleanIdMap <- function(data) {
  
  data <- as.data.frame(data)
  
  # new colnames, lowercase
  colnames(data) <- tolower(colnames(data))
  
  # replace '.' with '_'
  colnames(data) <- gsub('.', '_', colnames(data), fixed = T)
  
  # combine sentrix_id and sentrix_position 
  data$identifier <- paste(data$sentrix_barcode, data$sample_section, sep = '_')
  data$identifier <- as.factor(data$identifier)
  
  # keep only necessary columnes
  data <- data[, c('sample_id','identifier', 'sentrix_barcode', 'sample_section',
                   'sample_plate', 'sample_well')]
  
  return(data)
  
}

id_map <- cleanIdMap(id_map)

##########
# load idat files
##########

rgSetList <- list()

readIDAT <- function(path) {
  
  rgSet <- list()
  directory <- dir(path)
  
  for (folder in directory){
    # read into rgSet
    rgSet[[folder]] <- read.metharray.exp(paste(path, folder, sep = '/'))
  }
  
  return(rgSet)
  
}


rgSetList <- readIDAT(idat_data)

##########
# function that Loops through list, preprocesses, and convert to beta, m, and cn values 
##########
preprocessMethod <- function(data, preprocess) {
  
  Mset <- list()
  ratioSet <- list()
  gset <- list()
  beta <- list()
  
  for (dat in 1:length(data)) {
    
    if (preprocess == 'raw') {
      Mset[[dat]] <-preprocessRaw(data[[dat]])
      ratioSet[[dat]] <- ratioConvert(Mset[[dat]], what = 'both', keepCN = T)
    }
    
    if (preprocess == 'quan') {
      ratioSet[[dat]] <- preprocessQuantile(data[[dat]], fixOutliers = TRUE,
                                            removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                            quantileNormalize = TRUE, stratified = TRUE,
                                            mergeManifest = FALSE, sex = NULL)
    }
    
    
    gset[[dat]] <- mapToGenome(ratioSet[[dat]]) 
    beta[[dat]] <- getBeta(gset[[dat]])
    print(dat)
  }
  return(beta)
  
}

rgSet <- preprocessMethod(rgSetList, preprocess = 'raw')
rgSetQuan <- preprocessMethod(rgSetList, preprocess = 'quan')



##########
# Function that combines elements of a list into a data frame
##########
combineList <- function(list) {
  
  for (dat in 1:length(list)) {
    list[[dat]] <- t(list[[dat]])
  }
  data <- do.call(rbind, list)
  return(data)
}

raw_dat <- combineList(rgSet)
quan_dat <- combineList(rgSetQuan)


##########
# Function that combines methylation matrices with id_map, to get ids for methylation
##########
findIds <- function(data_methyl, id_map) {
  
  data_methyl <- as.data.frame(data_methyl)
  data_methyl$identifier <- rownames(data_methyl)
  data_methyl$identifier <- as.factor(data_methyl$identifier)
  # loop to combine identifiers, without merging large table
  data_methyl$ids <- NA
  data_methyl$sentrix_id <- NA
  
  for (i in data_methyl$identifier) {
    data_methyl$ids[data_methyl$identifier == i] <- id_map$sample_id[id_map$identifier == i]
    data_methyl$sentrix_id[data_methyl$identifier == i] <- id_map$sample_section[id_map$identifier == i]
    
    print(i)
  }
  
  return(data_methyl)
}

data_methyl <- findIds(raw_dat, id_map)
data_methyl_quan <- findIds(quan_dat, id_map)


# save data 
saveRDS(data_methyl, paste0(methyl_data, '/valid_raw.rda'))
saveRDS(data_methyl_quan, paste0(methyl_data, '/valid_quan.rda'))

