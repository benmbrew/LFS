####### Script will prepare data for pca 
# This is 5th step (A)

##########
# initialize libraries
##########
library(dplyr)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')

##########
# load beta values
##########
load(paste0(model_data, '/model_data_cases.RData'))
load(paste0(model_data, '/model_data_controls.RData'))

##########
# read in clinical data
##########
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# clean clinical ids
clin$id <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

##########
# functions to clean and join data
##########
# clean ids in each data set 
cleanIDs <- function(data)
{
  
  data$id <- gsub('A|B|_|-', '', data$id)
  data$id <- substr(data$id, 1,4) 
  return(data)
}

# get probe locations 

getIDAT <- function(cg_locations) 
{
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
joinData <- function(data, control) 
{
  # get intersection of clin ids and data ids
  intersected_ids <- intersect(data$id, clin$id)
  features <- colnames(data)[2:(length(colnames(data)))]
  
  # loop to combine identifiers, without merging large table
  data$p53_germline <- NA
  data$age_diagnosis <- NA
  data$cancer_diagnosis_diagnoses <- NA
  data$age_sample_collection <- NA
  data$tm_donor_ <- NA
  data$gdna.exon.intron <- NA
  data$gdna.base.change <- NA
  data$gdna.codon <- NA
  data$protein.codon.change <- NA
  data$protein.codon.num <- NA
  data$splice.delins.snv <- NA
  data$codon72.npro <- NA
  data$mdm2.nG <- NA
  data$gender <- NA
  
  if (!control) {
    
    for (i in intersected_ids) {
      
      data$p53_germline[data$id == i] <- clin$p53_germline[which(clin$id == i)]
      data$age_diagnosis[data$id == i] <- clin$age_diagnosis[which(clin$id == i)]
      data$cancer_diagnosis_diagnoses[data$id == i] <- clin$cancer_diagnosis_diagnoses[which(clin$id == i)]
      data$age_sample_collection[data$id == i] <- clin$age_sample_collection[which(clin$id == i)]
      data$tm_donor_[data$id == i] <- clin$tm_donor_[which(clin$id == i)]
      data$gdna.exon.intron[data$id == i] <- clin$gdna.exon.intron[which(clin$id == i)]
      data$gdna.base.change[data$id == i] <- clin$gdna.base.change[which(clin$id == i)]
      data$protein.codon.change[data$id == i] <- clin$protein.codon.change[which(clin$id == i)]
      data$protein.codon.num[data$id == i] <- clin$protein.codon.num[which(clin$id == i)]
      data$splice.delins.snv[data$id == i] <- clin$splice.delins.snv[which(clin$id == i)]
      data$codon72.npro[data$id == i] <- clin$codon72.npro[which(clin$id == i)]
      data$mdm2.nG[data$id == i] <- clin$mdm2.nG[which(clin$id == i)]
      data$gender[data$id == i] <- clin$gender[which(clin$id == i)]
      
      print(i)
    } 
    data <- data[!is.na(data$p53_germline),]
    data <- data[!duplicated(data$id),]
    data <- data[!duplicated(data$tm_donor_),]
    data <- data[, c('id', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses','age_sample_collection',
                     'gdna.exon.intron', 'gdna.base.change', 'protein.codon.change', 'protein.codon.num',
                     'splice.delins.snv', 'codon72.npro', 'mdm2.nG', 'gender',features)]
  } else {
    
    for (i in intersected_ids) {
      
      data$p53_germline[data$id == i] <- clin$p53_germline[which(clin$id == i)]
      data$cancer_diagnosis_diagnoses[data$id == i] <- clin$cancer_diagnosis_diagnoses[which(clin$id == i)]
      data$age_sample_collection[data$id == i] <- clin$age_sample_collection[which(clin$id == i)]
      data$tm_donor_[data$id == i] <- clin$tm_donor_[which(clin$id == i)]
      data$gdna.exon.intron[data$id == i] <- clin$gdna.exon.intron[which(clin$id == i)]
      data$gdna.base.change[data$id == i] <- clin$gdna.base.change[which(clin$id == i)]
      data$protein.codon.change[data$id == i] <- clin$protein.codon.change[which(clin$id == i)]
      data$protein.codon.num[data$id == i] <- clin$protein.codon.num[which(clin$id == i)]
      data$splice.delins.snv[data$id == i] <- clin$splice.delins.snv[which(clin$id == i)]
      data$codon72.npro[data$id == i] <- clin$codon72.npro[which(clin$id == i)]
      data$mdm2.nG[data$id == i] <- clin$mdm2.nG[which(clin$id == i)]
      data$gender[data$id == i] <- clin$gender[which(clin$id == i)]
      
      print(i)
    } 
    data <- data[!is.na(data$p53_germline),]
    data <- data[!duplicated(data$id),]
    # data <- data[!duplicated(data$tm_donor_),]
    data <- data[, c('id', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses','age_sample_collection',
                     'gdna.exon.intron', 'gdna.base.change', 'protein.codon.change', 'protein.codon.num',
                     'splice.delins.snv', 'codon72.npro', 'mdm2.nG', 'gender',features)]  }
  
  return(data)
}

# take p53 germline column and relevel factors to get rid of NA level
relevelFactor <- function (data) 
{
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

# test this with smaller data
makeNumFast <- function(model_data) 
{
  temp_list <- list()
  for (i in 14:ncol(model_data)) {
    temp.row <- model_data[,i]
    temp.row <-   as.numeric(levels(temp.row))[temp.row]
    temp_list[[i]] <- temp.row
    print(i)
    
  }
  
  num_features <- as.data.frame(do.call(cbind, temp_list))
  feature_names <- names(model_data)[14:ncol(model_data)]
  model_data <- cbind(id = model_data$id, 
                      p53_germline = model_data$p53_germline, 
                      age_diagnosis = model_data$age_diagnosis,
                      cancer_diagnosis_diagnoses = model_data$cancer_diagnosis_diagnoses, 
                      age_sample_collection = model_data$age_sample_collection, 
                      gdna.exon.intron = model_data$gdna.exon.intron, 
                      gdna.base.change = model_data$gdna.base.change,    
                      protein.codon.change = model_data$protein.codon.change,
                      protein.codon.num = model_data$protein.codon.num,
                      splice.delins.snv = model_data$splice.delins.snv,
                      codon72.npro = model_data$codon72.npro,
                      mdm2.nG = model_data$mdm2.nG,
                      gender = model_data$gender,
                      num_features)
  names(model_data)[14:ncol(model_data)] <- feature_names
  
  return(model_data)
}

##########
# apply functions to idat data - cases and controls and save to model_data folder
##########

##########
# First do cases
##########
getCases <- function(cases, make_num) 
{
  # first clean ids
  cases <- cleanIDs(cases)
  
  # second join data
  cases <- joinData(cases, control = F)
  
  # thrid relevel factors
  cases <- relevelFactor(cases)
  
  if (make_num) {
    # only for raw preprocessing method.
    #  4th make numeric - only cases, since that was imputed on and made into a factor
    cases <- makeNumFast(cases)
  }
  
  
  return(cases)
}

beta_raw <- getCases(beta_raw, make_num = T)
beta_swan <- getCases(beta_swan, make_num = F)
beta_quan <- getCases(beta_quan, make_num = F)
beta_funnorm <- getCases(beta_funnorm, make_num = F)

##########
# 2nd do controls
##########
getControls <- function(controls, make_num)
{
  # first clean ids
  controls <- cleanIDs(controls)
  
  # second join data
  controls <- joinData(controls, control = T)
  
  if (make_num) {
    # only for raw preprocessing
    # # 4th make numeric - only model_data, since that was imputed on and made into a factor
    controls <- makeNumFast(controls)
  }
  
  return(controls)
}

beta_raw_controls <- getControls(beta_raw_controls, make_num = T)
beta_swan_controls <- getControls(beta_swan_controls, make_num = F)
beta_quan_controls <- getControls(beta_quan_controls, make_num = F)
beta_funnorm_controls <- getControls(beta_funnorm_controls, make_num = F)


 
# # # save.image('/home/benbrew/Desktop/temp.num.RData')
# # load('/home/benbrew/Desktop/temp.num.RData')
# load(paste0(methyl_data, '/pca_num.RData'))


