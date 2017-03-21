####### Script will correct for batches by gender for cases and controls, and combine and correct for diff tech
# # This is 4th step

##########
# initialize libraries
##########
library(dplyr)
library(sva)
library(RColorBrewer)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
clin_data <- paste0(data_folder, '/clin_data')

##########
# load cases and controls data
##########
beta_quan <- readRDS(paste0(methyl_data, '/beta_quan.rda'))
beta_quan_controls <- readRDS(paste0(methyl_data, '/beta_quan_controls.rda'))

##########
# remove duplicated id columns
##########
beta_quan$ids.1 <- NULL
beta_quan_controls$ids.1 <- NULL
beta_quan$sentrix_id.1 <- NULL
beta_quan_controls$sentrix_id.1 <- NULL


##########
# get model data 
##########
getModData <- function(data) 
{
  # subset data by not na in age of diagnosis and mut
  data <- data[!is.na(data$age_diagnosis),]
  data <- data[data$p53_germline == 'Mut',]
  return(data)
}

beta_quan <- getModData(beta_quan)

##########
# get intersection of controls and normal
##########

getFeatInt <- function(data, data_controls)
{
  features <- names(data)[8:ncol(data)]
  features_controls <- names(data_controls)[8:ncol(data_controls)]
  
  # code for cases vs controls
  data$type <- 'cases'
  data_controls$type <- 'controls'
  
  # code fo sam vs not
  data$batch <- ifelse(grepl('9721365183', data$sentrix_id), 'yes', 'no')
  
  data_controls$batch <- 'no'
  
  data$type <- as.factor(data$type)
  data_controls$type <- as.factor(data_controls$type)
  
  
  intersect_features <- intersect(features, features_controls)
  
  data <- data[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                   'age_sample_collection' ,'gender', 'sentrix_id', 'batch', 'type', intersect_features)]
  
  data_controls <- data_controls[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                                     'age_sample_collection' ,'gender', 'sentrix_id', 'batch', 'type', intersect_features)]
  
  data <- rbind(data, data_controls)
  
  return(data)
  
}

quan <- getFeatInt(beta_quan, beta_quan_controls)

##########
# create a factor variable for age 
##########
features <- colnames(quan)[10:ncol(quan)]

quan$age_diagnosis_fac <- ifelse(quan$age_diagnosis > 0 & quan$age_diagnosis <= 6, 'min', 
                                 ifelse(quan$age_diagnosis > 6 & quan$age_diagnosis <= 27, '1st',
                                        ifelse(quan$age_diagnosis > 27 & quan$age_diagnosis <= 59, 'med',
                                               ifelse(quan$age_diagnosis > 59 & quan$age_diagnosis <= 175, 'mean',
                                                      ifelse(quan$age_diagnosis > 175 & quan$age_diagnosis <= 315, '3rd', 'max')))))


quan$age_sample_collection_fac <- ifelse(quan$age_sample_collection > 0 & quan$age_sample_collection <= 6, 'min', 
                                 ifelse(quan$age_sample_collection > 6 & quan$age_sample_collection <= 27, '1st',
                                        ifelse(quan$age_sample_collection > 27 & quan$age_sample_collection <= 59, 'med',
                                               ifelse(quan$age_sample_collection > 59 & quan$age_sample_collection <= 175, 'mean',
                                                      ifelse(quan$age_sample_collection > 175 & quan$age_sample_collection <= 315, '3rd', 'max')))))

##########
# create variable for ACC indicator
##########
quan$can_fac <- ifelse(grepl('ACC', quan$cancer_diagnosis_diagnoses), 'acc',
                       ifelse(grepl('breast', quan$cancer_diagnosis_diagnoses), 'breast_ca',
                              ifelse(grepl('ERMS', quan$cancer_diagnosis_diagnoses), 'erms', 'other')))

##########
# reorganize data
##########
quan <- quan[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                 'age_sample_collection', 'gender', 'sentrix_id', 'batch', 'type', 
                 'age_diagnosis_fac', 'age_sample_collection_fac', 'can_fac', features)]


########## 
# PCA of each data type and cases vs controls
##########
# functionp needs to take a clinical column, remove others, and plot pcas
getPCA <- function(pca_data, 
                   column_name, 
                   name, 
                   gene_start, 
                   controls, 
                   cases,
                   pca1,
                   pca2,
                   xlab,
                   ylab,
                   leg_title) 
{
  
  if (controls) {
    
    pca_data <- pca_data[pca_data$type == 'controls',]
    rownames(pca_data) <-pca_data$ids
    
  } else if (cases) {
    
    pca_data <- pca_data[pca_data$type == 'cases',]
    rownames(pca_data) <-pca_data$ids
    
    
  } 
  
  pca_data[, column_name] <- as.factor(pca_data[, column_name])
  
  # # subet data so only p53 mut
  # pca_data <- pca_data[pca_data$p53_germline == 'Mut',]
  # 
  # get features sites
  cg_sites <- colnames(pca_data)[gene_start:ncol(pca_data)]
  
  # subset by no NAs for column_name
  pca_data <- pca_data[!is.na(pca_data[, column_name]), ]
  
  stopifnot(!any(is.na(pca_data[, column_name])))
  
  # put column name with cg_sites 
  pca_data <- pca_data[ ,c(column_name, cg_sites)]
  
  # run pca
  data_length <- ncol(pca_data)
  pca <- prcomp(pca_data[,2:data_length])

  # plot data
  #fill in factors with colors 
  col_vec <- c('black','red' , 'darkblue', 'brown', 'darkgreen', 'yellow', 'orange', 
               'green', 'blueviolet', 'bisque', 'cyan', 'coral',
               'grey', 'bisque2', 'lightblue', 'purple','darkred', 
               'bisque1', 'darkorchid', 'gold', 'darkorange', 'deeppink',
               'greenyellow', 'beige')
  
  colors <- col_vec[pca_data[, column_name]]
  

  plot <- plot(pca$x[, pca1], 
               pca$x[, pca2],
               xlab = xlab,
               ylab = ylab,
               cex = 1,
               bty = 'n',
               main = name,
               pch = 16,
               col = adjustcolor(colors, alpha.f = 0.5)
  )
  
  # legend("bottomright", # position
  #        legend = c('Male', 'Female'),
  #        title = leg_title,
  #        fill = c('red', 'black'),
  #        cex = 0.56,
  #        bty = "n")
  abline(v = c(0,0),
         h = c(0,0))
  
  return(plot)
}


##########
# combine cases and controls and plot based on that
##########

##########
# cases
##########

##########
# cancer diagnosis
##########
# cancer
##########
# 1,2
getPCA(quan, 
       'can_fac', 
       'PCA',
       gene_start = 13,
       cases = T,
       controls = F,
       pca1 = 1,
       pca2 = 2,
       xlab = '1st PC',
       ylab = '2nd PC',
       leg_title = 'Cancer')

legend("bottomright", # position
       legend = unique(quan$can_fac),
       title = 'Gender',
       fill = c('red', 'black', 'blue'),
       cex = 0.7,
       bty = "n")

# 2,3
getPCA(quan, 
       'can_fac', 
       'PCA',
       gene_start = 13,
       cases = T,
       controls = F,
       pca1 = 2,
       pca2 = 3,
       xlab = '2nd PC',
       ylab = '3rd PC',
       leg_title = 'Cancer')

legend("bottomright", # position
       legend = unique(quan$can_fac),
       title = 'Gender',
       fill = c('red', 'black', 'blue'),
       cex = 0.7,
       bty = "n")



##########
# gender

#1,2
getPCA(quan, 
       'gender', 
       'PCA',
       gene_start = 12,
       cases = T,
       controls = F,
       pca1 = 1,
       pca2 = 2,
       xlab = '1st PC',
       ylab = '2nd PC')

legend("bottomright", # position
       legend = c('Male', 'Female'),
       title = 'Gender',
       fill = c('red', 'black'),
       cex = 0.7,
       bty = "n")

# 2,3
getPCA(quan, 
       'gender', 
       'PCA',
       gene_start = 12,
       cases = T,
       controls = F,
       pca1 = 2,
       pca2 = 3,
       xlab = '2nd PC',
       ylab = '3rd PC')

legend("bottomright", # position
       legend = c('Male', 'Female'),
       title = 'Gender',
       fill = c('red', 'black'),
       cex = 0.7,
       bty = "n")

##########
# batch
##########
# 1,2
getPCA(quan, 
       'batch', 
       'PCA',
       gene_start = 12,
       cases = T,
       controls = F,
       pca1 = 1,
       pca2 = 2,
       xlab = '1st PC',
       ylab = '2nd PC')

legend("bottomright", # position
       legend = c('Toronto', 'Quebec'),
       title = 'Batch Location',
       fill = c('red', 'black'),
       cex = 0.7,
       bty = "n")

# 2,3
getPCA(quan, 
       'batch', 
       'PCA',
       gene_start = 12,
       cases = T,
       controls = F,
       pca1 = 2,
       pca2 = 3,
       xlab = '2nd PC',
       ylab = '3rd PC',
       leg_title = 'Batch')

legend("bottomright", # position
       legend = c('Toronto', 'Quebec'),
       title = 'Batch Location',
       fill = c('red', 'black'),
       cex = 0.7,
       bty = "n")


##########
# age_diagnosis
##########
# 1,2
getPCA(quan, 
       'age_diagnosis_fac', 
       'PCA',
       gene_start = 12,
       cases = T,
       controls = F,
       pca1 = 1,
       pca2 = 2,
       xlab = '1st PC',
       ylab = '2nd PC',
       leg_title = 'age diagnosis')

legend("bottomright", # position
       legend = c('age1', 'age2', 'age3', 'age4', 'age5', 'age6'),
       title = 'Age',
       fill = c('black','red' , 'darkblue', 'brown', 'darkgreen', 'yellow'),
       cex = 0.5,
       pch = 16,
       bty = "n")


# 2,3
getPCA(quan, 
       'age_diagnosis_fac', 
       'PCA',
       gene_start = 12,
       cases = T,
       controls = F,
       pca1 = 2,
       pca2 = 3,
       xlab = '2nd PC',
       ylab = '3rd PC',
       leg_title = 'age diagnosis')

legend("bottomright", # position
       legend = c('age1', 'age2', 'age3', 'age4', 'age5', 'age6'),
       title = 'Age',
       fill = c('black','red' , 'darkblue', 'brown', 'darkgreen', 'yellow'),
       cex = 0.5,
       pch = 16,
       bty = "n")



##########
# age_sample_collection_fac
##########
# 1,2
getPCA(quan, 
       'age_sample_collection_fac', 
       'PCA',
       gene_start = 12,
       cases = T,
       controls = F,
       pca1 = 1,
       pca2 = 2,
       xlab = '1st PC',
       ylab = '2nd PC')
legend("bottomright", # position
       legend = c('age1', 'age2', 'age3', 'age4', 'age5', 'age6'),
       title = 'Age',
       fill = c('black','red' , 'darkblue', 'brown', 'darkgreen', 'yellow'),
       cex = 0.5,
       pch = 16,
       bty = "n")

# 2,3
getPCA(quan, 
       'age_sample_collection_fac', 
       'PCA',
       gene_start = 12,
       cases = T,
       controls = F,
       pca1 = 2,
       pca2 = 3,
       xlab = '2nd PC',
       ylab = '3rd PC')

legend("bottomright", # position
       legend = c('age1', 'age2', 'age3', 'age4', 'age5', 'age6'),
       title = 'Age',
       fill = c('black','red' , 'darkblue', 'brown', 'darkgreen', 'yellow'),
       cex = 0.5,
       pch = 16,
       bty = "n")







##########
# controls
##########

getPCA(quan, 
       'gender', 
       'PCA quan controls gender',
       gene_start = 10,
       cases = F,
       controls = T,
       pca1 = 1,
       pca2 = 2)


# cases sentrix
getPCA(quan, 
       'sentrix_id', 
       'PCA quan cases sentrix',
       gene_start = 10,
       cases = T,
       controls = F,
       pca1 = 1,
       pca2 = 2)

# cases sentrix
getPCA(quan, 
       'sentrix_id', 
       'PCA quan cases sentrix',
       gene_start = 10,
       cases = T,
       controls = F,
       pca1 = 2,
       pca2 = 3)


