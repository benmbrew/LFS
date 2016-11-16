## Script for cleaning IDAT files
# this data could be used to replace methylation from nardine
# this script will be documented heavily as to provide a guide for IDAT files in minfi package to lab.

# initialize libraries
library(minfi)
library(dplyr)

# Use this https://www.bioconductor.org/help/course-materials/2014/BioC2014/minfi_BioC2014.pdf

# initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/classification_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
idat_data <- paste0(methyl_data, '/raw_files')
idat_data_controls <- paste0(methyl_data, '/controls')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


####################
# Load id_map to later merge with methylation values. clean it.
####################
id_map <- read.csv(paste0(methyl_data, '/ids_map.csv'), stringsAsFactors = F)
id_map_control <- read.csv(paste0(methyl_data, '/ids_map_controls.csv'), stringsAsFactors = F)


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




####################
# this function will read in each idat file using read.450k.exp from minifi package and store them in a list
####################

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

rgSetList <- readIDAT(idat_data_controls)

####################
# function that Loops through list, preprocesses, and convert to beta, m, and cn values 
####################

preprocessMethod <- function(data, preprocess) {
  
  Mset <- list()
  ratioSet <- list()
  gset <- list()
  beta <- list()
  m <- list()
  cn <- list()
  
  for (dat in 1:length(data)) {
    
    if (preprocess == 'raw') {
      Mset[[dat]] <-preprocessRaw(data[[dat]])
      ratioSet[[dat]] <- ratioConvert(Mset[[dat]], what = 'both', keepCN = T)
    } else if (preprocess == 'illumina') {
      Mset[[dat]] <-preprocessIllumina(data[[dat]])
      ratioSet[[dat]] <- ratioConvert(Mset[[dat]], what = 'both', keepCN = T)
    } else if (preprocess == 'swan') {
      Mset[[dat]] <-preprocessSWAN(data[[dat]])
      ratioSet[[dat]] <- ratioConvert(Mset[[dat]], what = 'both', keepCN = T)
    } else if (preprocess == 'funnorm') {
      ratioSet[[dat]] <-preprocessFunnorm(data[[dat]])
    } else if (preprocess == 'quantile') {
      ratioSet[[dat]] <- preprocessQuantile(data[[dat]], fixOutliers = TRUE,
                                            removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                            quantileNormalize = TRUE, stratified = TRUE,
                                            mergeManifest = FALSE, sex = NULL)
    }
    gset[[dat]] <- mapToGenome(ratioSet[[dat]]) 
    beta[[dat]] <- getBeta(gset[[dat]])
    m[[dat]] <- getM(gset[[dat]])
    cn[[dat]] <- getCN(gset[[dat]])
    print(dat)
  }
  return(list(beta, m, cn))
}


####################
# Function that combines elements of a list into a data frame
####################

combineList <- function(list) {
  
  for (dat in 1:length(list)) {
    list[[dat]] <- t(list[[dat]])
  }
  data <- do.call(rbind, list)
  return(data)
}

####################
# Function that combines methylation matrices with id_map, to get ids for methylation
####################

findIds <- function(data_methyl, id_map) {
  
  data_methyl <- as.data.frame(data_methyl)
  data_methyl$identifier <- rownames(data_methyl)
  data_methyl$identifier <- as.factor(data_methyl$identifier)
  # loop to combine identifiers, without merging large table
  data_methyl$ids <- NA
  for (i in data_methyl$identifier) {
    data_methyl$ids[data_methyl$identifier == i] <- id_map$sample_name[id_map$identifier == i]
    print(i)
  }
  
  return(data_methyl)
}

####################
# Function that gets malkin id from sample_name
####################
getIdName <- function(data) {
  
  column_split <- strsplit(as.character(data$ids), '#')
  last_digits <- lapply(column_split, function(x) x[length(x)])
  sub_ids <- unlist(last_digits)
  sub_ids <- gsub('RD-', '', sub_ids)
  data$ids <- sub_ids
  data$sentrix_position <- NULL
  data$sentrix_id <- NULL
  data$pool_id <- NULL
  data$sample_group <- NULL
  data$sample_plate <- NULL
  data$sample_well <- NULL
  data$sample_name <- NULL
  data$identifier <- NULL
  return(data)
  
}

###########################
# Main function that specifies a preprocessing method and get beta, m, cn_values
###########################

getMethyl <- function(data_list, method, control) {
 
  processed_list <- list()
  processed_list <-preprocessMethod(data_list, preprocess = method)
  
  # combinelist
  beta_methyl <- combineList(processed_list[[1]])
  m_methyl <- combineList(processed_list[[2]])
  cn_methyl <- combineList(processed_list[[3]])
  
  # find ids
  beta_methyl <- findIds(beta_methyl, id_map_control)
  m_methyl <- findIds(m_methyl, id_map_control)
  cn_methyl <- findIds(cn_methyl, id_map_control)
  
  if(!control) {
    
    # clean ids
    beta_methyl <- getIdName(beta_methyl)
    m_methyl <- getIdName(m_methyl)
    cn_methyl <- getIdName(cn_methyl)
  }
 
  
  return(list(beta_methyl, m_methyl, cn_methyl))
  
  
}

# run the 5 different preprocessing and normilization methods
raw <- getMethyl(rgSetList, method = 'raw', control = T)
beta_raw <- raw[[1]]
m_raw <- raw[[2]]
cn_raw <- raw[[3]]
rm(raw)
save.image(paste0(methyl_data, '/raw_control.RData'))
rm(beta_raw, m_raw, cn_raw)

# illum <- getMethyl(rgSetList, method = 'illumina', control = T)
# beta_illumina <- illum[[1]]
# m_illumina <- illum[[2]]
# cn_illumina <- illum[[3]]
# rm(illum)
# save.image(paste0(methyl_data, '/illum.RData'))
# rm(m_illumina, beta_illumina, cn_illumina)

swan <- getMethyl(rgSetList, method = 'swan', control = T)
beta_swan <- swan[[1]]
m_swan <- swan[[2]]
cn_swan <- swan[[3]]
rm(swan)
save.image(paste0(methyl_data, '/swan_control.RData'))
rm(beta_swan, m_swan, cn_swan)

funnorm <- getMethyl(rgSetList, method = 'funnorm', control = T)
beta_funnorm <- funnorm[[1]]
m_funnorm <- funnorm[[2]]
cn_funnorm <- funnorm[[3]]
rm(funnorm)
save.image(paste0(methyl_data, '/funnorm_control.RData'))
rm(beta_funnorm, m_funnorm, cn_funnorm)

quan <- getMethyl(rgSetList, method = 'quantile', control = T)
beta_quan <- quan[[1]]
m_quan <- quan[[2]]
cn_quan <- quan[[3]]
rm(quan)
save.image(paste0(methyl_data, '/quan_control.RData'))
rm(beta_quan, m_quan, cn_quan)












