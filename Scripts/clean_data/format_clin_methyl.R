###############################################
# This script will subset features and merge with clin saving different versions of the data to be
# this is the 6th step in the pipeline
library(dplyr)
library(stringr)
library(impute)
library(mlbench)
library(caret)
library(FactoMineR)
library(genefilter)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(impute)
library(GenomicRanges)
library(biovizBase)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)


# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
imputed_data <- paste0(data_folder, '/imputed_data')
idat_data <- paste0(methyl_data, '/raw_files')
model_data <- paste0(data_folder, '/model_data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')
clin_data <- paste0(data_folder, '/clin_data')


#################################################################################################
# Read in methyl and clinical data and join by ids
#################################################################################################

# Read in data (clinical or clinical_two)
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# clin clinical ids
clin$id <-  gsub('A|B|_', '', clin$blood_dna_malkin_lab_)

# load methylation data - probe, gene, knn, lsa
load(paste0(imputed_data, '/imputed_gene_probe.RData'))

# load idat methylation imputerd- probe, knn - raw, swan, quan, funnorm
load(paste0(idat_data, '/imputed_idat_betas.RData'))
load(paste0(idat_data, '/imputed_idat_betas_control.RData'))

# reformat to look like other data so it fits in to function

# remove unneeded object
rm(methyl_gene, methyl_probe, data, gene, probe)

#############################################
# functions to clean and join data
#############################################

# clean ids in each data set 
cleanIDs <- function(data){
  
  data <- as.data.frame(data)
  data$id <- gsub('A|B|_', '', data$id)
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
  
 
  if(!control) {
    # loop to combine identifiers, without merging large table
    data$p53_germline <- NA
    data$age_diagnosis <- NA
    data$cancer_diagnosis_diagnoses <- NA
    data$age_sample_collection <- NA
    data$tm_donor_ <- NA
    
    for (i in intersected_ids) {
      
      data$p53_germline[data$id == i] <- clin$p53_germline[which(clin$id == i)]
      data$age_diagnosis[data$id == i] <- clin$age_diagnosis[which(clin$id == i)]
      data$cancer_diagnosis_diagnoses[data$id == i] <- clin$cancer_diagnosis_diagnoses[which(clin$id == i)]
      data$age_sample_collection[data$id == i] <- clin$age_sample_collection[which(clin$id == i)]
      data$tm_donor_[data$id == i] <- clin$tm_donor_[which(clin$id == i)]
      
      print(i)
    } 
    data <- data[!is.na(data$p53_germline),]
    data <- data[!duplicated(data$id),]
    data <- data[!duplicated(data$tm_donor_),]
    data <- data[, c('id', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses','age_sample_collection', features)]
  } else {
    # loop to combine identifiers, without merging large table
    data$p53_germline <- NA
    data$cancer_diagnosis_diagnoses <- NA
    data$age_sample_collection <- NA
    data$tm_donor_ <- NA
    
    for (i in intersected_ids) {
      
      data$p53_germline[data$id == i] <- clin$p53_germline[which(clin$id == i)]
      data$cancer_diagnosis_diagnoses[data$id == i] <- clin$cancer_diagnosis_diagnoses[which(clin$id == i)]
      data$age_sample_collection[data$id == i] <- clin$age_sample_collection[which(clin$id == i)]
      data$tm_donor_[data$id == i] <- clin$tm_donor_[which(clin$id == i)]
      
      print(i)
    } 
    data <- data[!is.na(data$p53_germline),]
    data <- data[!duplicated(data$id),]
    data <- data[!duplicated(data$tm_donor_),]
    data <- data[, c('id', 'p53_germline', 'cancer_diagnosis_diagnoses','age_sample_collection', features)]
  }
  
  return(data)
}

# take p53 germline column and relevel factors to get rid of NA level
relevelFactor <- function(data) {
  
  data$p53_germline <- factor(data$p53_germline, levels = c('Mut', 'WT'))
  return(data)
}

# Function to convert all genes/probe columns to numeric
makeNum <- function(model_data) {
  
  model_data[, 6:ncol(model_data)] <- apply(model_data[, 6:ncol(model_data)], 2, function(x) as.numeric(as.character(x)))
  
  return(model_data)
}

reformatData <- function(data) {
  data <- as.data.frame(data)
  names(data)[1] <- 'id'
  return(data)
}


###########################
# apply functions to gene and probe original dta
###########################

gene_knn <- cleanIDs(gene_knn)
gene_lsa <- cleanIDs(gene_lsa)
probe_knn <- cleanIDs(probe_knn)
probe_lsa <- cleanIDs(probe_lsa)

gene_knn <- joinData(gene_knn)
gene_lsa <- joinData(gene_lsa)
probe_lsa <- joinData(probe_lsa)
probe_knn <- joinData(probe_knn)

gene_knn <- relevelFactor(gene_knn)
gene_lsa <- relevelFactor(gene_lsa)
probe_knn <- relevelFactor(probe_knn)
probe_lsa <- relevelFactor(probe_lsa)

gene_knn <- makeNum(gene_knn)
gene_lsa <- makeNum(gene_lsa)
probe_knn <- makeNum(probe_knn)
probe_lsa <- makeNum(probe_lsa)

cg_locations <- getIDAT(cg_locations)

# save image file 
save.image(paste0(model_data, '/model_data.RData'))


###########################
# apply functions to idat data
###########################
beta_raw <- reformatData(beta_raw)
beta_swan <- reformatData(beta_swan)
beta_quan <- reformatData(beta_quan)
beta_funnorm <- reformatData(beta_funnorm)

beta_raw <- cleanIDs(beta_raw)
beta_quan <- cleanIDs(beta_quan)
beta_swan <- cleanIDs(beta_swan)
beta_funnorm <- cleanIDs(beta_funnorm)

beta_raw <- joinData(beta_raw, control = T)
beta_quan <- joinData(beta_quan, control = T)
beta_swan <- joinData(beta_swan)
beta_funnorm <- joinData(beta_funnorm)

beta_raw <- relevelFactor(beta_raw)
beta_quan <- relevelFactor(beta_quan)
beta_swan <- relevelFactor(beta_swan)
beta_funnorm <- relevelFactor(beta_funnorm)

# load(paste0(imputed_data, '/imputed_gene_probe.RData'))
beta_raw <- makeNum(beta_raw)
beta_quan <- makeNum(beta_quan)
beta_swan <- makeNum(beta_swan)
beta_funnorm <- makeNum(beta_funnorm)

# save image file 
# save.image(paste0(idat_data, '/imputed_idat_betas_final.RData'))
save.image(paste0(idat_data, '/imputed_idat_betas_final_control.RData'))


# ###################################################################################################
# # Using correlation
# # scale methyl
# #methyl <- scale(methyl[, -1])
# rm(agg1, d1,d2,d3,d4, )
# # make a correlation matrix 
# cor_mat <- cor(methyl[, -1])
# 
# # find attributes that are highly correlated
# # 0.6 for methyl_cor, 0.4, methyl_cor_small and full_data_cor_small
# highly_cor <- findCorrelation(cor_mat, cutoff = 0.5, names = TRUE)
# cor_index <- names(methyl) %in% highly_cor[2:length(highly_cor)]
# 
# # remove highly correlated attributes
# methyl_cor <- methyl[, !cor_index]
# names(methyl_cor)[1] <- 'id'
# 
# # write.csv(methyl_cor, paste0(data_folder, '/methyl_cor_small.csv'))
# # write.csv(methyl_cor, paste0(data_folder, '/methyl_cor.csv'))
# 
# 
# # inner_join clin
# full_data_cor <- inner_join(clin, methyl_cor,
#                         by = 'id')
# 
# # Save data to be used later
# # write.csv(full_data_cor, paste0(data_folder, '/full_data_cor_small.csv'))
# # write.csv(full_data_cor, paste0(data_folder, '/full_data_cor.csv'))
# 
# ######################################################################################################
# # Random forest recursive feature elimination
# # prepare training scheme
# control <- trainControl(method="repeatedcv", number=5, repeats=3)
# 
# # subset full data for the model 
# x_mat <- full_data[, c(6, 30:ncol(full_data))]
# x_mat <- x_mat[complete.cases(x_mat),]
# y <- as.numeric(x_mat$age_diagnosis)
# 
# # train the model
# model <- train(x = x_mat[,-1],
#                y = y,
#                preProcess = 'scale',
#                importance = TRUE,
#                trControl = control)
# 
# # estimate variable importance
# importance <- varImp(model)
# 
# # get vector of importance 
# importance <- importance$importance
# 
# # make rownames a column
# importance$gene <- rownames(importance)
# rownames(importance) <- NULL
# 
# # sort importance vector
# importance <- importance[order(importance$Overall, decreasing = T),]
# 
# # save importance vector
# # write.csv(importance, paste0(data_folder, '/importance.csv'))
# # plot importance, maybe cutoff around 20%
# hist(importance$Overall)
# 
# # subset data by top features 
# subset <- importance[importance$Overall > 26,]
# full_data_rf <- full_data[, c(names(full_data)[1:30], subset$gene)]
# 
# write.csv(full_data_rf, paste0(data_folder, '/full_data_rf.csv'))
# 
# # #########################################################################################
# # # Using low variance selection
# # subset full data for the model 
# x_mat <- full_data[, c(6, 30:ncol(full_data))]
# x_mat <- x_mat[complete.cases(x_mat),]
# y <- as.numeric(x_mat$age_diagnosis)
# 
# # #  First find the desired quantile breaks for the entire matrix
# qt <- quantile(data.matrix(x_mat[, -1]) , 0.1 )
# # # 20%  80%
# # #5.17 6.62
# # #  Next get a logical vector of the rows that have any values outside these breaks
# columns <- apply(x_mat[, -1] , 2, function(x) any( x < qt[1]))
# # #  Subset on this vector
# temp <- x_mat[ ,columns]
# #
# # # Use genefilter
# temp <- varFilter(t(temp[, -1]))
# 
# ########################################################################################################
# # # Using PCA
# # pca <- PCA(x_mat[,-1])
# #
# # #This line of code will sort the variables the most linked to each PC.
# # # It is very useful when you have many variables.
# # temp <- dimdesc(pca)
# # temp_new <- temp$Dim.1
# # temp1 <- temp_new$quanti
# 
# #########################################################################################################
# # # nearZeroVar function
# # temp <- nearZeroVar(x_mat[, -1], freqCut = 50/5, saveMetrics = TRUE)
# #
# # temp <- apply(methyl[, -1], 2, function(i) var(i))
# 
# 
# #########################################################################
# # subset to just three columns- id, methlation_indicator, age_sample_collection, p53_germline, cancer_diagnosis
# dat <- clin[, c('p53_germline', 'cancer_diagnosis_diagnoses', 'age_sample_collection', 'blood_dna_malkin_lab_',
#                 'methyl_indicator')]
# 
# range <- 12
# samples <- list()
# for ( i in 1:nrow(dat)) {
#   
#   if(dat$methyl_indicator[i] == 'Yes'){
#     
#     temp <- dat$age_sample_collection[i]
#     
# 
#     samples[[i]] <- dat[(dat$age_sample_collection > temp & dat$age_sample_collection < (temp + range)) |
#                   (dat$age_sample_collection < temp & dat$age_sample_collection > (temp - range)),] 
#   }
# 
# }
# 
# 
# samples <- do.call('rbind', samples)
# samples <- samples[!duplicated(samples$blood_dna_malkin_lab_),]
# 
# # subset to mutant, unaffected, no methyl 
# ids <- samples[samples$cancer_diagnosis_diagnoses == 'Unaffected' & samples$p53_germline == 'Mut' & 
#                  samples$methyl_indicator == 'No',]
# 
# ids <- ids[complete.cases(ids),]
# 
