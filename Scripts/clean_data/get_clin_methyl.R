####### Script will combine methylation and clinical data
# this is 3th step in pipeline - check to see we are taking the first diagnosis

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
# quan
beta_quan <- readRDS(paste0(methyl_data, '/beta_quan.rda'))
beta_quan_controls <- readRDS(paste0(methyl_data, '/beta_quan_controls.rda'))

# funnorm
beta_funnorm <- readRDS(paste0(methyl_data, '/beta_funnorm.rda'))
beta_funnorm_controls <- readRDS(paste0(methyl_data, '/beta_funnorm_controls.rda'))

##########
# make data frames
##########
#quan
beta_quan <- as.data.frame(beta_quan, stringsAsFactors = F)
beta_quan_controls <- as.data.frame(beta_quan_controls, stringAsFactors = F)

#funnorm
beta_funnorm <- as.data.frame(beta_funnorm, stringsAsFactors = F)
beta_funnorm_controls <- as.data.frame(beta_funnorm_controls, stringAsFactors = F)


##########
# read in clinical data
##########
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# clean clinical idss
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

##########
# function to clean ids names and get methyl_indicator for clin
##########
getMethylVar <- function(dat_cases, dat_controls) {
  
  # get idss from cases and controls
  cases_names <- as.character(dat_cases$ids)
  controls_names <- as.character(dat_controls$ids)
  
  # combine idss
  methylation_names <- append(cases_names, controls_names)
  
  # remove 'A' and '_' in methylation names
  methylation_names <- gsub('A|B|_|-', '', methylation_names)
  
  # add '1' and make data frame 
  methylation_names <- as.data.frame(cbind(methylation_names, rep.int(1, length(methylation_names))))
  methylation_names$methylation_names <- as.character(methylation_names$methylation_names)
  methylation_names$V2 <- as.character(methylation_names$V2)
  names(methylation_names) <- c("ids", "methyl_indicator")
  methylation_names <- methylation_names[!duplicated(methylation_names),]
  
  # keep only first 4 characters of methylation_names$ids
  methylation_names$ids <- substr(methylation_names$ids, 1,4)
  
  # add methyl_indicator column to clin
  clin$methyl_indicator <- NA
  for (i in methylation_names$ids) {
    clin$methyl_indicator[clin$ids == i] <- methylation_names$methyl_indicator[methylation_names$ids == i]
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

# clean idss in each data set 
cleanids <- function(data){
  
  data$ids <- gsub('A|B|_|-', '', data$ids)
  data$ids <- substr(data$ids, 1,4) 
  return(data)
}
# get probe locations 

getIds <- function(cg_locations) {
  
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

# function that takes each methylation and merges with clinical - keep ids, family, p53 status, age data
joinData <- function(data, control) {
  
  # get intersection of clin idss and data idss
  intersected_ids <- intersect(data$ids, clin$ids)
  features <- colnames(data)[1:(length(colnames(data)) - 2)]
  
  # loop to combine idsentifiers, without merging large table
  data$p53_germline <- NA
  data$age_diagnosis <- NA
  data$cancer_diagnosis_diagnoses <- NA
  data$age_sample_collection <- NA
  data$tm_donor_ <- NA
  data$gender <- NA
  
  if (!control) {
    
    for (i in intersected_ids) {
      
      data$p53_germline[data$ids == i] <- clin$p53_germline[which(clin$ids == i)]
      data$age_diagnosis[data$ids == i] <- clin$age_diagnosis[which(clin$ids == i)]
      data$cancer_diagnosis_diagnoses[data$ids == i] <- clin$cancer_diagnosis_diagnoses[which(clin$ids == i)]
      data$age_sample_collection[data$ids == i] <- clin$age_sample_collection[which(clin$ids == i)]
      data$tm_donor_[data$ids == i] <- clin$tm_donor_[which(clin$ids == i)]
      data$gender[data$ids == i] <- clin$gender[which(clin$ids == i)]

      
      
      print(i)
    } 
    data <- data[!is.na(data$p53_germline),]
    data <- data[!duplicated(data$ids),]
    data <- data[!duplicated(data$tm_donor_),]
    # data <- data[!is.na(data$age_diagnosis),]
    # data <- data[!is.na(data$age_sample_collection), ]
    
    data <- data[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses',
                     'age_sample_collection', 'gender','sentrix_id','sen_batch', features)]
  } else {
    
    for (i in intersected_ids) {
      
      data$p53_germline[data$ids == i] <- clin$p53_germline[which(clin$ids == i)]
      data$cancer_diagnosis_diagnoses[data$ids == i] <- clin$cancer_diagnosis_diagnoses[which(clin$ids == i)]
      data$age_sample_collection[data$ids == i] <- clin$age_sample_collection[which(clin$ids == i)]
      data$tm_donor_[data$ids == i] <- clin$tm_donor_[which(clin$ids == i)]
      data$gender[data$ids == i] <- clin$gender[which(clin$ids == i)]
      
      
      print(i)
    } 
    data <- data[!is.na(data$p53_germline),]
    # data <- data[!duplicated(data$ids),]
    # data <- data[!duplicated(data$tm_donor_),]
    data <- data[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses',
                     'age_sample_collection', 'gender', 'sentrix_id', features)]
  }
  
  return(data)
}

# take p53 germline column and relevel factors to get rids of NA level
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
# ##########
# # get clinical methylation indicator
# clin <- getMethylVar(beta_quan, beta_quan_controls)
# 
# # get cg locations
# cg_locations <- getIds()
# 
# write.csv(cg_locations, paste0(model_data, '/cg_locations.csv'))
# write.csv(clin, paste0(clin_data, '/clinical_two.csv'))
##########
# First do cases
##########

# quan
# first clean idss
beta_quan <- cleanids(beta_quan)

options(warn=1)

# second join data
beta_quan <- joinData(beta_quan, control = F)

# thrids relevel factors
beta_quan <- relevelFactor(beta_quan)


# funnorm
# first clean idss
beta_funnorm <- cleanids(beta_funnorm)

options(warn=1)

# second join data
beta_funnorm <- joinData(beta_funnorm, control = F)

# thrids relevel factors
beta_funnorm <- relevelFactor(beta_funnorm)

##########
# 2nd do controls
##########

# quan
# first clean idss
beta_quan_controls <- cleanids(beta_quan_controls)

# second join data
beta_quan_controls <- joinData(beta_quan_controls, control = T)

# funnorm
# first clean idss
beta_funnorm_controls <- cleanids(beta_funnorm_controls)

# second join data
beta_funnorm_controls <- joinData(beta_funnorm_controls, control = T)


#########
# save data
#########
#quan
saveRDS(beta_quan, paste0(methyl_data, '/beta_quan.rda'))


saveRDS(beta_quan_controls, paste0(methyl_data, '/beta_quan_controls.rda'))

#funnorm
saveRDS(beta_funnorm, paste0(methyl_data, '/beta_funnorm.rda'))


saveRDS(beta_funnorm_controls, paste0(methyl_data, '/beta_funnorm_controls.rda'))

