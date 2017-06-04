####### Script will correct for batches by gender for cases and controls, and combine and correct for diff tech
# # This is 4th step

##########
# initialize libraries
##########
library(dplyr)
library(sva)
library(impute)

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
# Read in methylation probe and gene
##########
# raw
beta_raw <- readRDS(paste0(methyl_data, '/beta_raw.rda'))
beta_raw_controls <- readRDS(paste0(methyl_data, '/beta_raw_controls.rda'))
beta_raw_valid <- readRDS(paste0(methyl_data, '/beta_raw_valid.rda'))

##########
# remove id.1
##########
beta_raw_controls$ids.1 <- beta_raw$ids.1 <-NULL

##########
# remove identifier in controls
##########
beta_raw_controls$identifier <- NULL

##########t
# make data frames
##########
#raw
beta_raw <- as.data.frame(beta_raw, stringsAsFactors = F)
beta_raw_controls <- as.data.frame(beta_raw_controls, stringAsFactors = F)
beta_raw_valid <- as.data.frame(beta_raw_valid, stringAsFactors = F)

##########
# noramlize data after subsetting featres
##########

# data_cases <- beta_raw
# data_controls <- beta_raw_controls
normalizeDat <- function(data_cases, data_controls)
{
  features <- names(data_cases)[9:ncol(data_cases)]
  features_controls <- names(data_controls)[8:ncol(data_controls)]
  data_controls$sen_batch <- 'tor_2'
  
  
  intersect_features <- intersect(features, features_controls)
  data_cases <- data_cases[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                   'age_sample_collection' ,'gender', 'sentrix_id', 'sen_batch', intersect_features)]
  
  data_controls <- data_controls[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                                     'age_sample_collection' ,'gender', 'sentrix_id', 'sen_batch', intersect_features)]
  
  scaleDat <- function(dat) {
    
    clin_dat <- dat[, c(1:8)]
    dat <- t(dat[, c(9:ncol(dat))])
    # get row statistics
    rowMean <- apply(dat, 1, mean, na.rm=TRUE)
    rowSd <- apply(dat, 1, sd, na.rm=TRUE)
    # constantInd <- rowSd==0
    # rowSd[constantInd] <- 1
    rowStats <- list(mean=rowMean, sd=rowSd)
    
    # apply normilization
    dat  <- (dat - rowStats$mean) / rowStats$sd
    
    # make matrix
    dat <- as.matrix(dat)

    # impute with knn
    dat_knn <-  impute.knn(dat, k = 10)$data
    
    # transpose 
    dat_knn <- t(dat_knn)
    
    # get clin data back and return final data
    final_dat <- cbind(clin_dat, dat_knn)
    
    return(final_dat)
  }
  
  # apply function
  data_cases <- scaleDat(data_cases)
  data_controls <- scaleDat(data_controls)

  
  return(list(data_cases, data_controls))
}

norm_dat <- normalizeDat(beta_raw, beta_raw_controls)
cases <- norm_dat[[1]]
controls <- norm_dat[[2]]
rm(norm_dat, beta_raw, beta_raw_controls)

##########
# get WT no cancer
##########
controls_wt <- cases[which(cases$cancer_diagnosis_diagnoses == 'Unaffected' & 
                                cases$p53_germline == 'WT'),]

##########
# remove cancers from controls
##########
removeCancer <- function(data_controls) 
{
  data_controls <- data_controls[grepl('Unaffected', data_controls$cancer_diagnosis_diagnoses),]
  data_controls <- data_controls[!duplicated(data_controls$ids)]
  return(data_controls)
}

controls <- removeCancer(controls)

##########
# get extra controls
##########
getExtraCon <- function(data) 
{
  # subset data by not na in age of diagnosis and mut
  data <- data[data$p53_germline == 'Mut' & data$cancer_diagnosis_diagnoses == 'Unaffected',]
  data <- data[!is.na(data$age_sample_collection),]
  return(data)
}

# extra controls
controls_2 <- getExtraCon(cases)

# combine controls and controls_2 for full_controls
controls_full <- rbind(controls, controls_2)
rm(controls_2)

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

# funnorm
cases <- getModData(cases)

# get sub cases
cases_sub <- cases[!grepl('mon', cases$sen_batch),]


########## 
# PCA of each data type and cases vs controls
# ##########
# pca_data <- controls_wt
# column_name <- 'gender'
# name <- 'l'
# gene_start <- 9
# pca1<- 'x'
# pca2<- 'y'
# functionp needs to take a clinical column, remove others, and plot pcas
getPCA <- function(pca_data, 
                   column_name, 
                   name, 
                   gene_start, 
                   pca1,
                   pca2) 
{
  
  
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
  temp <- pca$x[,1]
  
  temp

  # plot data 4257,  94, 93
  #fill in factors with colors 
  col_vec <- c('black','red' , 'green', 'bisque', 'bisque1', 'bisque2', 'lightblue', 
               'blueviolet', 'brown', 'cyan', 'coral',
               'grey', 'orange', 'yellow', 'darkblue','darkred', 
               'darkgreen', 'darkorchid', 'gold', 'darkorange', 'deeppink',
               'greenyellow', 'purple')
  
  colors <- col_vec[pca_data[, column_name]]
  
  
  plot <- plot(pca$x[, pca1], 
               pca$x[, pca2],
               xlab = 'pca',
               ylab = 'pca',
               cex = 1,
               main = name,
               pch = 16,
               col = adjustcolor(colors, alpha.f = 0.5)
  )
  abline(v = c(0,0),
         h = c(0,0))
  
  return(plot)
}

#########
# cases
#########

# gender 
getPCA(cases, 
       'gender', 
       'PCA cases gender',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)

#cases
getPCA(cases, 
       'sen_batch', 
       'PCA cases sen_batch',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)


#cases
getPCA(cases, 
       'sentrix_id', 
       'PCA cases sentrix_id',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)


#########
# cases_full
#########

# gender 
getPCA(cases_sub, 
       'gender', 
       'PCA cases_sub gender',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)

#cases_full
getPCA(cases_sub, 
       'sen_batch', 
       'PCA cases_sub sen_batch',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)


#cases_full
getPCA(cases_sub, 
       'sentrix_id', 
       'PCA cases_sub sentrix_id',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)

#########
# controls
#########

# gender 
getPCA(controls, 
       'gender', 
       'PCA controls gender',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(controls, 
       'sen_batch', 
       'PCA controls sen_batch',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)



#controls
getPCA(controls, 
       'sentrix_id', 
       'PCA controls sentrix_id',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)

#########
# controls_full
#########

# gender 
getPCA(controls_full, 
       'gender', 
       'PCA controls_full gender',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(controls_full, 
       'sen_batch', 
       'PCA controls_full sen_batch',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)



#controls
getPCA(controls_full, 
       'sentrix_id', 
       'PCA controls_full sentrix_id',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)

#########
# controls_wt
#########

# gender 
getPCA(controls_wt, 
       'gender', 
       'PCA controls_wt gender',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(controls_wt, 
       'sen_batch', 
       'PCA controls_wt sen_batch',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)



#controls
getPCA(controls_wt, 
       'sentrix_id', 
       'PCA controls_wt sentrix_id',
       gene_start = 9,
       pca1 = 1,
       pca2 = 2)

##########
# remove outliers  4257 cases, 3391, 3392 controls
##########
removeOutlier <- function(data, wt) {
  

  #controls outlier
  data <- data[data$ids != '3391',]
  data <- data[data$ids != '3392',]
  
  if(wt) {
    data <- data[data$ids != '2564',]
    
  }
  
  # data <- data[data$ids != '4257',]
  # data <- data[data$ids != '94',]
  # data <- data[data$ids != '93',]
  # data <- data[data$ids != '2414',]
  
  
  return(data)
}

# raw
cases <- removeOutlier(cases, wt = F)
controls <- removeOutlier(controls, wt = F)
controls_full <- removeOutlier(controls_full, wt = F)
controls_wt <- removeOutlier(controls_wt, wt = T)

##########
# batch
# ##########
# 
# # # HERE
# getBatch <- function(data, cases, data_size)
# {
# 
#   if(cases) {
#     # subset to cases
#     data <- data[!grepl('Unaffected', data$cancer_diagnosis_diagnoses),]
# 
#     if(data_size == 'sub') {
# 
#       data <- data[!grepl('mon', data$sen_batch),]
# 
#     } else {
#       # get batch
#       batch_indicator <- as.character(data$sen_batch)
#       batch_indicator <- as.factor(batch_indicator)
# 
#       gender <- data$gender
#       sentrix_id <- data$sentrix_id
#       sen_batch <- data$sen_batch
#       ids <- data$ids
#       sample_collection <- data$age_sample_collection
#       diagnosis <- data$age_diagnosis
#       p53 <- data$p53_germline
#       cancer_diagnosis <- data$cancer_diagnosis_diagnoses
#       # put model ids in rownames and remove columns
#       mat_data <- data[, 9:ncol(data)]
#       # get features
#       features <- colnames(mat_data)
#       mat_data <- t(mat_data)
# 
#       # get intercept
#       modcombat <- model.matrix(~1, data = data)
#       combat <- ComBat(dat = mat_data, batch = batch_indicator, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
# 
#       # transpose and add back columns
#       final_dat <- as.data.frame(t(combat))
#       final_dat$gender <- gender
#       final_dat$sen_batch <- sen_batch
#       final_dat$ids <- ids
#       final_dat$age_sample_collection <- sample_collection
#       final_dat$age_diagnosis <- diagnosis
#       final_dat$p53_germline <- p53
#       final_dat$sentrix_id <- sentrix_id
#       final_dat$cancer_diagnosis_diagnoses <- cancer_diagnosis
#       final_dat <- final_dat[, c('ids', 'p53_germline', 'age_diagnosis','cancer_diagnosis_diagnoses' ,
#                                  'age_sample_collection', 'gender', 'sentrix_id', 'sen_batch', features)]
#       rownames(final_dat) <- NULL
# 
# 
#     }
# 
#   } else {
#     # subset to cases
#     data <- data[grepl('Unaffected', data$cancer_diagnosis_diagnoses),]
# 
#     if(data_size == 'sub') {
# 
#       data <- data[grepl('tor_2', data$sen_batch),]
# 
#     } else {
#       # get batch
#       batch_indicator <- as.character(data$sen_batch)
#       batch_indicator <- as.factor(batch_indicator)
# 
#       gender <- data$gender
#       sentrix_id <- data$sentrix_id
#       sen_batch <- data$sen_batch
#       ids <- data$ids
#       sample_collection <- data$age_sample_collection
#       diagnosis <- data$age_diagnosis
#       p53 <- data$p53_germline
#       cancer_diagnosis <- data$cancer_diagnosis_diagnoses
#       # put model ids in rownames and remove columns
#       mat_data <- data[, 10:ncol(data)]
#       # get features
#       features <- colnames(mat_data)
#       mat_data <- t(mat_data)
# 
#       # get intercept
#       modcombat <- model.matrix(~1, data = data)
#       combat <- ComBat(dat = mat_data, batch = batch_indicator, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
# 
#       # transpose and add back columns
#       final_dat <- as.data.frame(t(combat))
#       final_dat$gender <- gender
#       final_dat$sen_batch <- sen_batch
#       final_dat$ids <- ids
#       final_dat$age_sample_collection <- sample_collection
#       final_dat$age_diagnosis <- diagnosis
#       final_dat$p53_germline <- p53
#       final_dat$sentrix_id <- sentrix_id
#       final_dat$cancer_diagnosis_diagnoses <- cancer_diagnosis
#       final_dat <- final_dat[, c('ids', 'p53_germline', 'age_diagnosis','cancer_diagnosis_diagnoses' ,
#                                 'age_sample_collection', 'gender', 'sentrix_id', 'sen_batch', features)]
#       rownames(final_dat) <- NULL
# 
# 
#     }
# 
#   }
# 
#   if(data_size == 'sub') {
#     return(data)
#   } else {
#     return(final_dat)
# 
#   }
# 
# 
# }
# #
# 

# ##########
# # full data
# ##########
# 
# # cases
# quan_cases_full <- getBatch(quan, cases = T, data_size = 'full')
# funnorm_cases_full <- getBatch(funnorm, cases = T, data_size = 'full')
# raw_cases_full <- getBatch(raw, cases = T, data_size = 'full')
# 
# 
# # cases
# quan_controls_full <- getBatch(quan, cases = F, data_size = 'full')
# funnorm_controls_full <- getBatch(funnorm, cases = F, data_size = 'full')
# raw_controls_full <- getBatch(raw, cases = F, data_size = 'full')

##########
# # sub data
# ##########
# 
# # cases
# cases_sub <- getBatch(cases, cases = T, data_size = 'sub')
# funnorm_cases_sub <- getBatch(funnorm, cases = T, data_size = 'sub')
# raw_cases_sub <- getBatch(raw, cases = T, data_size = 'sub')
# 
# # cases
# quan_controls_sub <- getBatch(quan, cases = F, data_size = 'sub')
# funnorm_controls_sub <- getBatch(funnorm, cases = F, data_size = 'sub')
# raw_controls_sub <- getBatch(raw, cases = F, data_size = 'sub')
# 
# # remove unneeded data
# rm(quan, funnorm, raw)



##########
# full 
##########

# save cases 
saveRDS(cases, paste0(model_data, '/cases.rda'))

# save contorls
saveRDS(controls, paste0(model_data, '/controls.rda'))

# sub

# save cases 
saveRDS(cases_sub, paste0(model_data, '/cases_sub.rda'))

# save contorls
saveRDS(controls_full, paste0(model_data, '/controls_full.rda'))

# save controls wt
saveRDS(controls_wt, paste0(model_data, '/controls_wt.rda'))




