####### Script will combine methylation and clinical data
# this is 4th step in pipeline

##########
# initialize libraries
##########
library(dplyr)
library(stringr)
library(impute)
library(data.table)
library(impute)
library(GenomicRanges)
library(biovizBase)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)
library(preprocessCore)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
imputed_data <- paste0(data_folder, '/imputed_data')
idat_data <- paste0(methyl_data, '/raw_files')
model_data <- paste0(data_folder, '/model_data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')
clin_data <- paste0(data_folder, '/clin_data')

##########
# Read in methylation probe and gene
##########
beta_raw <- readRDS(paste0(methyl_data, '/beta_raw.rda'))
beta_raw_controls <- readRDS(paste0(methyl_data, '/beta_raw_controls.rda'))

##########
# make data frames
##########
beta_raw <- as.data.frame(beta_raw)
beta_raw_controls <- as.data.frame(beta_raw_controls)


##########
# read in clinical data
##########
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# clean clinical ids
clin$id <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

##########
# function to clean id names and get methyl_indicator for clin
##########
getMethylVar <- function(dat_cases, dat_controls) {
  
  # get ids from cases and controls
  cases_names <- as.character(dat_cases$id)
  controls_names <- as.character(dat_controls$id)
  
  # combine ids
  methylation_names <- append(cases_names, controls_names)
  
  # remove 'A' and '_' in methylation names
  methylation_names <- gsub('A|B|_|-', '', methylation_names)
  
  # add '1' and make data frame 
  methylation_names <- as.data.frame(cbind(methylation_names, rep.int(1, length(methylation_names))))
  methylation_names$methylation_names <- as.character(methylation_names$methylation_names)
  methylation_names$V2 <- as.character(methylation_names$V2)
  names(methylation_names) <- c("id", "methyl_indicator")
  methylation_names <- methylation_names[!duplicated(methylation_names),]
  
  # keep only first 4 characters of methylation_names$id
  methylation_names$id <- substr(methylation_names$id, 1,4)
  
  # add methyl_indicator column to clin
  clin$methyl_indicator <- NA
  for (i in methylation_names$id) {
    clin$methyl_indicator[clin$id == i] <- methylation_names$methyl_indicator[methylation_names$id == i]
  }
  
  # recode 1 = TRUE, FALSE otherwise
  clin$methyl_indicator <- ifelse(clin$methyl_indicator == 1, TRUE, FALSE)
  clin$methyl_indicator[is.na(clin$methyl_indicator)] <- FALSE
  
  return(clin)
  
  # summary(clin$methyl_indicator) # 193
}


##########
# functions to clean and join data
##########
# clean ids in each data set 
cleanIDs <- function(data){
  
  data$id <- gsub('A|B|_|-', '', data$id)
  data$id <- substr(data$id, 1,4) 
  return(data)
}

# get probe locations 

getIDAT <- function(cg_locations) {
  
  #idat files
  idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
  sapply(idatFiles, gunzip, overwrite = TRUE)
  # read into rgSet
  rgSet <- read.450k.exp("GSE68777/idat")
  # preprocess quantil
  rgSet <- preprocessQuantile(rgSet)
  # get rangers 
  rgSet <- granges(rgSet)
  cg_locations <- as.data.frame(rgSet)
  # make rownames probe column
  cg_locations$probe <- rownames(cg_locations)
  rownames(cg_locations) <- NULL
  return(cg_locations)
}


# function that takes each methylation and merges with clinical - keep id, family, p53 status, age data
joinData <- function(data, control) {
  
  # get intersection of clin ids and data ids
  intersected_ids <- intersect(data$id, clin$id)
  features <- colnames(data)[2:(length(colnames(data)))]
  
  # loop to combine identifiers, without merging large table
  data$p53_germline <- NA
  data$age_diagnosis <- NA
  data$cancer_diagnosis_diagnoses <- NA
  data$age_sample_collection <- NA
  data$tm_donor_ <- NA
  data$gender <- NA
  
  if (!control) {
    
    for (i in intersected_ids) {
      
      data$p53_germline[data$id == i] <- clin$p53_germline[which(clin$id == i)]
      data$age_diagnosis[data$id == i] <- clin$age_diagnosis[which(clin$id == i)]
      data$cancer_diagnosis_diagnoses[data$id == i] <- clin$cancer_diagnosis_diagnoses[which(clin$id == i)]
      data$age_sample_collection[data$id == i] <- clin$age_sample_collection[which(clin$id == i)]
      data$tm_donor_[data$id == i] <- clin$tm_donor_[which(clin$id == i)]
      data$gender[data$id == i] <- clin$gender[which(clin$id == i)]
      
      
      print(i)
    } 
    data <- data[!is.na(data$p53_germline),]
    data <- data[!duplicated(data$id),]
    data <- data[!duplicated(data$tm_donor_),]
    data <- data[, c('id', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses',
                     'age_sample_collection', 'gender', features)]
  } else {
    
    for (i in intersected_ids) {
      
      data$p53_germline[data$id == i] <- clin$p53_germline[which(clin$id == i)]
      data$cancer_diagnosis_diagnoses[data$id == i] <- clin$cancer_diagnosis_diagnoses[which(clin$id == i)]
      data$age_sample_collection[data$id == i] <- clin$age_sample_collection[which(clin$id == i)]
      data$tm_donor_[data$id == i] <- clin$tm_donor_[which(clin$id == i)]
      data$gender[data$id == i] <- clin$gender[which(clin$id == i)]
      
      
      print(i)
    } 
    data <- data[!is.na(data$p53_germline),]
    data <- data[!duplicated(data$id),]
    # data <- data[!duplicated(data$tm_donor_),]
    data <- data[, c('id', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses',
                     'age_sample_collection', 'gender', features)]
  }
  
  return(data)
}

# take p53 germline column and relevel factors to get rid of NA level
relevelFactor <- function (data) {
  
  data$p53_germline <- factor(data$p53_germline, levels = c('Mut', 'WT'))
  return(data)
}

# # Function to convert all genes/probe columns to numeric
# makeNum <- function (model_data) {
#   
#   model_data[, 6:ncol(model_data)] <- apply (model_data[, 6:ncol(model_data)], 2, function(x) as.numeric(as.character(x)))
#   
#   return(model_data)
# }


##########
# apply functions to idat data - cases and controls and save to model_data folder
##########
# get clinical methylation indicator
clin <- getMethylVar(beta_raw, beta_raw_controls)

# get cg locations
cg_locations <- getIDAT()

# write.csv(cg_locations, paste0(model_data, '/cg_locations.csv'))
# write.csv(clin, paste0(clin_data, '/clinical_two.csv'))
##########
# First do cases
##########
# first clean ids
beta_raw <- cleanIDs(beta_raw)

# second join data
beta_raw <- joinData(beta_raw, control = F)

# thrid relevel factors
beta_raw <- relevelFactor(beta_raw)

##########
# 2nd do controls
##########
# first clean ids
beta_raw_controls <- cleanIDs(beta_raw_controls)

# second join data
beta_raw_controls <- joinData(beta_raw_controls, control = T)

#########
# save data
#########
saveRDS(beta_raw, paste0(methyl_data, '/beta_raw.rda'))

saveRDS(beta_raw_controls, paste0(methyl_data, '/beta_raw_controls.rda'))
