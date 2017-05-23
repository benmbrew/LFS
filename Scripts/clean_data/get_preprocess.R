###### Script for cleaning IDAT files
# this is 2nd step in pipeline

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
methyl_data <- paste0(data_folder, '/methyl_data')
idat_data <- paste0(methyl_data, '/raw_files')
idat_data_controls <- paste0(methyl_data, '/controls')
clin_data <- paste0(data_folder, '/clin_data')


##########
# Load id_map, which has ids to merge with methylation - do for cases and controls
##########
# read in id map for model data

# batch1 
id_map <- read.csv(paste0(methyl_data, '/ids_map.csv'), stringsAsFactors = F)

#batch2
id_map_other <- read.csv(paste0(methyl_data, '/batch_2014.csv'), stringsAsFactors = F)
id_map_other$Project <- NULL

# controls
id_map_control <- read.csv(paste0(methyl_data, '/ids_map_controls.csv'), stringsAsFactors = F)

##########
# combine id_map and id_map_other
##########
id_map <- rbind(id_map, id_map_other)


##########
# clean IdMap
##########

cleanIdMap <- function(data) {
  
  data <- as.data.frame(data)
  
  # new colnames, lowercase
  colnames(data) <- tolower(colnames(data))
  
  # combine sentrix_id and sentrix_position 
  data$identifier <- paste(data$sentrix_id, data$sentrix_position, sep = '_')
  data$identifier <- as.factor(data$identifier)
  
  return(data)
  
}

id_map <- cleanIdMap(id_map)
id_map_control <- cleanIdMap(id_map_control)


##########
# this function will read in each idat file using read.450k.exp from minifi package and store them in a list
##########

rgSetList <- list()
rgSetListControls <- list()

readIDAT <- function(path) {
  
  rgSet <- list()
  directory <- dir(path)
  
  for (folder in directory){
    # read into rgSet
    rgSet[[folder]] <- read.metharray.exp(paste(path, folder, sep = '/'))
  }
  
  return(rgSet)
  
}

rgSetListControls <- readIDAT(idat_data_controls)

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
    

    gset[[dat]] <- mapToGenome(ratioSet[[dat]]) 
    beta[[dat]] <- getBeta(gset[[dat]])
    print(dat)
  }
  return(beta)
   
}


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
    data_methyl$ids[data_methyl$identifier == i] <- id_map$sample_name[id_map$identifier == i]
    data_methyl$sentrix_id[data_methyl$identifier == i] <- id_map$sentrix_id[id_map$identifier == i]
    
    print(i)
  }
  
  return(data_methyl)
}

##########
# Function that gets malkin id from sample_name
##########

getIdName <- function(data) {
  
  column_split <- strsplit(as.character(data$ids), '#')
  last_digits <- lapply(column_split, function(x) x[length(x)])
  sub_ids <- unlist(last_digits)
  sub_ids <- gsub('RD-', '', sub_ids)
  data$ids <- sub_ids
  data$sentrix_position <- NULL
  data$pool_id <- NULL
  data$sample_group <- NULL
  data$sample_plate <- NULL
  data$sample_well <- NULL
  data$sample_name <- NULL
  data$identifier <- NULL
  return(data)
  
}
##########
# Main function that specifies a preprocessing method and get beta
##########
getMethyl <- function(data_list,control, method) {
  
  processed_list <-preprocessMethod(data_list, preprocess = method)
  
  # save.image('/home/benbrew/Desktop/temp_process.RData')
  # load('/home/benbrew/Desktop/temp_process.RData')
  # 
  # combinelist
  beta_methyl <- combineList(processed_list)
  # m_methyl <- combineList(processed_list[[2]])
  # cn_methyl <- combineList(processed_list[[3]])
  
  if (!control) {
    
    beta_methyl <- findIds(beta_methyl, id_map)
    # clean ids
    beta_methyl <- getIdName(beta_methyl)
    # m_methyl <- getIdName(m_methyl)
    # cn_methyl <- getIdName(cn_methyl)
  } else {
    
    # find ids
    beta_methyl <- findIds(beta_methyl, id_map_control)
    
    # m_methyl <- findIds(m_methyl, id_map_control)
    # cn_methyl <- findIds(cn_methyl, id_map_control)
  }
  
  
  return(beta_methyl)
  
}

##########
# apply functions to get both cases and controls data
##########
# raw
beta_raw <- getMethyl(rgSetList, control = F, method = 'raw')
beta_raw_controls <- getMethyl(rgSetListControls, control = T, method = 'raw')

# scale betas




##########
# new variable called sen_batch
#########

beta_raw$sen_batch <- ifelse(grepl('9721365183', rownames(beta_raw)), 'mon', 'tor_1')

# save data
saveRDS(beta_raw, paste0(methyl_data, '/beta_raw.rda'))
saveRDS(beta_raw_controls, paste0(methyl_data, '/beta_raw_controls.rda'))




