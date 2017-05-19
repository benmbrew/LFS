####### Script will correct for batches by gender for cases and controls, and combine and correct for diff tech
# # This is 4th step

##########
# initialize libraries
##########
library(dplyr)
library(sva)

##########
# test scale and normalization
##########
# temp <- as.data.frame(replicate(10,sample(0:10,1000,rep=TRUE)))
# 
# 
# 
# temp_scale <- scale(temp)
# data <- temp
# testNorm <- function(data) {
#   
#   data <- t(data)
#   # get col statistics
#   colMean <- apply(data, 1, mean, na.rm=TRUE)
#   colSd <- apply(data, 1, sd, na.rm=TRUE)
#   constantInd <- colSd==0
#   colSd[constantInd] <- 1
#   colStats <- list(mean=colMean, sd=colSd, ind=constantInd)
#   
#   # apply normilization
#   data  <- (data - colStats$mean) / colStats$sd
#   # data <- data[!colStats$ind, ]
#   
#   data <- t(data)
#   return(data)
#   
# }
# temp_row <- testNorm(temp)
# temp_col <- testNorm(temp)

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
# quan
beta_quan <- readRDS(paste0(methyl_data, '/beta_quan.rda'))
beta_quan_controls <- readRDS(paste0(methyl_data, '/beta_quan_controls.rda'))

# funnorm
beta_funnorm <- readRDS(paste0(methyl_data, '/beta_funnorm.rda'))
beta_funnorm_controls <- readRDS(paste0(methyl_data, '/beta_funnorm_controls.rda'))

# raw
beta_raw <- readRDS(paste0(methyl_data, '/beta_raw.rda'))
beta_raw_controls <- readRDS(paste0(methyl_data, '/beta_raw_controls.rda'))

##########
# remove id.1
##########
beta_quan_controls$ids.1 <- beta_quan$ids.1 <- beta_funnorm_controls$ids.1 <- beta_funnorm$ids.1 <- 
  beta_raw_controls$ids.1 <- beta_raw$ids.1 <- NULL

##########
# remove identifier in controls
##########
beta_funnorm_controls$identifier <- beta_quan_controls$identifier <- beta_raw_controls$identifier <- NULL

##########
# normalize raw 
##########
data <- beta_raw
controls <- T
normDat <- function(data, controls) 
{
  if(controls){
    # first get clin dat 
    clin_dat <- data[, 1:7]
    
    # now remove clin part of nov_dat and convert to matrix for scale
    data <- as.matrix(t(data[, -c(1:7)]))
    

    stopifnot(!is.na(data))
  } else {
    # first get clin dat 
    clin_dat <- data[, 1:8]
    
    # now remove clin part of nov_dat and convert to matrix for scale
    data <- as.matrix(t(data[, -c(1:8)]))
  }
  
  # data <- apply(data, 2, function(x) as.numeric(x))
  
  # get row statistics
  rowMean <- apply(data, 1, mean, na.rm=TRUE)
  rowSd <- apply(data, 1, sd, na.rm=TRUE)
  constantInd <- rowSd==0
  rowSd[constantInd] <- 1
  rowStats <- list(mean=rowMean, sd=rowSd, ind=constantInd)
  
  # apply normilization
  data  <- (data - rowStats$mean) / rowStats$sd
  # data <- data[!rowStats$ind, ]
  
  data <- t(data)
  
  data <- cbind(clin_dat, data)
  
  return(data)
}

beta_raw <- normDat(beta_raw, controls = F)
beta_raw_controls <- normDat(beta_raw_controls, controls = T)


##########
# make data frames
##########
#quan
beta_quan <- as.data.frame(beta_quan, stringsAsFactors = F)
beta_quan_controls <- as.data.frame(beta_quan_controls, stringAsFactors = F)

#funnorm
beta_funnorm <- as.data.frame(beta_funnorm, stringsAsFactors = F)
beta_funnorm_controls <- as.data.frame(beta_funnorm_controls, stringAsFactors = F)

#raw
beta_raw <- as.data.frame(beta_raw, stringsAsFactors = F)
beta_raw_controls <- as.data.frame(beta_raw_controls, stringAsFactors = F)


##########
# remove cancers from controls
##########
removeCancer <- function(data_controls) 
{
  data_controls <- data_controls[grepl('Unaffected', data_controls$cancer_diagnosis_diagnoses),]
  data_controls <- data_controls[!duplicated(data_controls$ids)]
  return(data_controls)
}

beta_quan_controls <- removeCancer(beta_quan_controls)
beta_funnorm_controls <- removeCancer(beta_funnorm_controls)
beta_raw_controls <- removeCancer(beta_raw_controls)

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
beta_quan_controls_2 <- getExtraCon(beta_quan)
beta_funnorm_controls_2 <- getExtraCon(beta_funnorm)
beta_raw_controls_2 <- getExtraCon(beta_raw)


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

# quan
beta_quan <- getModData(beta_quan)

# funnorm
beta_funnorm <- getModData(beta_funnorm)

# funnorm
beta_raw <- getModData(beta_raw)


##########
# get intersection of controls and normal
##########
# data <- beta_quan
# data_controls <- beta_quan_controls
# data_controls_2 <- beta_quan_controls_2

getFeatInt <- function(data, data_controls, data_controls_2)
{
  features <- names(data)[9:ncol(data)]
  features_controls <- names(data_controls)[8:ncol(data_controls)]
  features_controls_2 <- names(data_controls_2)[9:ncol(data_controls_2)]
  
  data_controls$sen_batch <- 'tor_2'
  
  
  intersect_features <- Reduce(intersect, list(features, features_controls, features_controls_2))
  
  data <- data[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                   'age_sample_collection' ,'gender', 'sentrix_id', 'sen_batch', intersect_features)]
  
  data_controls <- data_controls[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                                     'age_sample_collection' ,'gender', 'sentrix_id', 'sen_batch', intersect_features)]
  
  data_controls_2 <- data_controls_2[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                                     'age_sample_collection' ,'gender', 'sentrix_id', 'sen_batch', intersect_features)]
  
  data <- rbind(data, data_controls, data_controls_2)
  
  data <- data[!duplicated(data$ids),]

  return(data)
  
}

# quan
quan <- getFeatInt(beta_quan, beta_quan_controls, beta_quan_controls_2)

# funnorm
funnorm <- getFeatInt(beta_funnorm, beta_funnorm_controls, beta_funnorm_controls_2)

# raw
raw <- getFeatInt(beta_raw, beta_raw_controls, beta_raw_controls_2)

# remove uneeded objects
rm(list=ls(pattern="beta*"))

nrow(quan[quan$sen_batch == 'tor_1',])
nrow(quan[quan$sen_batch == 'tor_2',])

########## 
# PCA of each data type and cases vs controls
##########
# functionp needs to take a clinical column, remove others, and plot pcas
pca_data <- raw
getPCA <- function(pca_data, 
                   column_name, 
                   name, 
                   gene_start, 
                   case_control,
                   batch,
                   pca1,
                   pca2) 
{
  
  if (case_control == "cases") {
    
    pca_data <- pca_data[!grepl('Unaffected', pca_data$cancer_diagnosis_diagnoses),]
    rownames(pca_data) <-pca_data$ids
    
    
  } 
  if (case_control == "controls") {
    
    pca_data <- pca_data[grepl('Unaffected', pca_data$cancer_diagnosis_diagnoses),]
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
  temp <- pca$x[,1]

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


##########
# combine cases and controls and plot based on that
##########

#####################################################################
##########
# quan
##########
# cases
getPCA(quan, 
       'gender', 
       'PCA quan cases gender',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(quan, 
       'gender', 
       'PCA quan controls gender',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)


#cases
getPCA(quan, 
       'sen_batch', 
       'PCA quan cases sen_batch',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(quan, 
       'sen_batch', 
       'PCA quan controls sen_batch',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)

#cases
getPCA(quan, 
       'sentrix_id', 
       'PCA quan cases sentrix_id',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(quan, 
       'sentrix_id', 
       'PCA quan controls sentrix_id',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)


#####################################################################
##########
# raw
##########
# cases
getPCA(raw, 
       'gender', 
       'PCA raw cases gender',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(raw, 
       'gender', 
       'PCA raw controls gender',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)


#cases
getPCA(raw, 
       'sen_batch', 
       'PCA raw cases sen_batch',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(raw, 
       'sen_batch', 
       'PCA raw controls sen_batch',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)

#cases
getPCA(raw, 
       'sentrix_id', 
       'PCA raw cases sentrix_id',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(raw, 
       'sentrix_id', 
       'PCA raw controls sentrix_id',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)


##########
# remove raw outliers
##########



#############################
#cases
getPCA(funnorm, 
       'sen_batch', 
       'PCA funnorm cases sen_batch',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(funnorm, 
       'sen_batch', 
       'PCA funnorm controls sen_batch',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)

#cases
getPCA(funnorm, 
       'sentrix_id', 
       'PCA funnorm cases sentrix_id',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(funnorm, 
       'sentrix_id', 
       'PCA funnorm controls sentrix_id',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)




##########
# remove outliers  4257 cases, 3391, 3392 controls
##########
removeOutlier <- function(data, funnorm) {
  

  #controls outlier
  data <- data[data$ids != '3391',]
  data <- data[data$ids != '3392',]
  data <- data[data$ids != '4257',]
  data <- data[data$ids != '94',]
  data <- data[data$ids != '93',]
  data <- data[data$ids != '2414',]
  
  
  return(data)
}

# quan
quan <- removeOutlier(quan, funnorm = F)

# funnorm
funnorm <- removeOutlier(funnorm, funnorm = F)

# raw
raw <- removeOutlier(raw, funnorm  = F)

# ##########
# # create variable to loop through to correct for sentrix_id
# ##########
# rmSentrixBatch <- function(data) 
# {
#   data_types <- unique(data$sen_batch)
#   
#   result_list <- list()
#   
#   for(i in data_types) {
#     # get sub data
#     sub_data <- data[data$sen_batch == i,]
#     
#     batch_indicator <- as.factor(as.character(sub_data$sentrix_id))
#     
#     # get columns
#     gender <- sub_data$gender
#     sentrix_id <- sub_data$sentrix_id
#     sen_batch <- sub_data$sen_batch
#     ids <- sub_data$ids
#     sample_collection <- sub_data$age_sample_collection
#     diagnosis <- sub_data$age_diagnosis
#     p53 <- sub_data$p53_germline
#     cancer_diagnosis <- sub_data$cancer_diagnosis_diagnoses
#     # put model ids in rownames and remove columns
#     mat_sub_data <- sub_data[, 9:ncol(sub_data)]
#     # get features 
#     features <- colnames(mat_sub_data)
#     mat_sub_data <- t(mat_sub_data)
#     
#     # get intercept
#     modcombat <- model.matrix(~1, data = sub_data)
#     combat <- ComBat(dat = mat_sub_data, batch = batch_indicator, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
#     
#     # transpose and add back columns
#     final_dat <- as.data.frame(t(combat))
#     final_dat$gender <- gender
#     final_dat$sen_batch <- sen_batch
#     final_dat$ids <- ids
#     final_dat$age_sample_collection <- sample_collection
#     final_dat$age_diagnosis <- diagnosis
#     final_dat$p53_germline <- p53
#     final_dat$sentrix_id <- sentrix_id
#     final_dat$cancer_diagnosis_diagnoses <- cancer_diagnosis
#     final_dat <- final_dat[, c('ids', 'gender', 'age_sample_collection',
#                                'age_diagnosis','cancer_diagnosis_diagnoses' ,'p53_germline', 'sen_batch', 'sentrix_id',features)]
#     rownames(final_dat) <- NULL
#     
#     result_list[[data_types]] <- final_dat
#     
#     print(data_types)
# 
#     
#   }
#   
#   return(result_list)
#   
# }
# 
# quan <- rmSentrixBatch(quan)
# 
# funnorm <- rmSentrixBatch(funnorm)
# 

##########
# batch
##########
getBatch <- function(data, cases, data_size)
{
  
  if(cases) {
    # subset to cases
    data <- data[!grepl('Unaffected', data$cancer_diagnosis_diagnoses),]
    
    if(data_size == 'sub') {
      
      data <- data[!grepl('mon', data$sen_batch),]
      
    } else {
      # get batch
      batch_indicator <- as.character(data$sen_batch)
      batch_indicator <- as.factor(batch_indicator)
      
      gender <- data$gender
      sentrix_id <- data$sentrix_id
      sen_batch <- data$sen_batch
      ids <- data$ids
      sample_collection <- data$age_sample_collection
      diagnosis <- data$age_diagnosis
      p53 <- data$p53_germline
      cancer_diagnosis <- data$cancer_diagnosis_diagnoses
      # put model ids in rownames and remove columns
      mat_data <- data[, 9:ncol(data)]
      # get features 
      features <- colnames(mat_data)
      mat_data <- t(mat_data)
      
      # get intercept
      modcombat <- model.matrix(~1, data = data)
      combat <- ComBat(dat = mat_data, batch = batch_indicator, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
      
      # transpose and add back columns
      final_dat <- as.data.frame(t(combat))
      final_dat$gender <- gender
      final_dat$sen_batch <- sen_batch
      final_dat$ids <- ids
      final_dat$age_sample_collection <- sample_collection
      final_dat$age_diagnosis <- diagnosis
      final_dat$p53_germline <- p53
      final_dat$sentrix_id <- sentrix_id
      final_dat$cancer_diagnosis_diagnoses <- cancer_diagnosis
      final_dat <- final_dat[, c('ids', 'p53_germline', 'age_diagnosis','cancer_diagnosis_diagnoses' ,
                                 'age_sample_collection', 'gender', 'sentrix_id', 'sen_batch', features)]
      rownames(final_dat) <- NULL
      
      
    } 
    
  } else {
    # subset to cases
    data <- data[grepl('Unaffected', data$cancer_diagnosis_diagnoses),]
    
    if(data_size == 'sub') {
      
      data <- data[grepl('tor_2', data$sen_batch),]
      
    } else {
      # get batch
      batch_indicator <- as.character(data$sen_batch)
      batch_indicator <- as.factor(batch_indicator)
      
      gender <- data$gender
      sentrix_id <- data$sentrix_id
      sen_batch <- data$sen_batch
      ids <- data$ids
      sample_collection <- data$age_sample_collection
      diagnosis <- data$age_diagnosis
      p53 <- data$p53_germline
      cancer_diagnosis <- data$cancer_diagnosis_diagnoses
      # put model ids in rownames and remove columns
      mat_data <- data[, 10:ncol(data)]
      # get features 
      features <- colnames(mat_data)
      mat_data <- t(mat_data)
      
      # get intercept
      modcombat <- model.matrix(~1, data = data)
      combat <- ComBat(dat = mat_data, batch = batch_indicator, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
      
      # transpose and add back columns
      final_dat <- as.data.frame(t(combat))
      final_dat$gender <- gender
      final_dat$sen_batch <- sen_batch
      final_dat$ids <- ids
      final_dat$age_sample_collection <- sample_collection
      final_dat$age_diagnosis <- diagnosis
      final_dat$p53_germline <- p53
      final_dat$sentrix_id <- sentrix_id
      final_dat$cancer_diagnosis_diagnoses <- cancer_diagnosis
      final_dat <- final_dat[, c('ids', 'p53_germline', 'age_diagnosis','cancer_diagnosis_diagnoses' ,
                                'age_sample_collection', 'gender', 'sentrix_id', 'sen_batch', features)]
      rownames(final_dat) <- NULL
      
      
    } 
    
  }
  
  if(data_size == 'sub') {
    return(data)
  } else {
    return(final_dat)
    
  }
  
  
}



##########
# full data
##########

# cases
quan_cases_full <- getBatch(quan, cases = T, data_size = 'full')
funnorm_cases_full <- getBatch(funnorm, cases = T, data_size = 'full')
raw_cases_full <- getBatch(raw, cases = T, data_size = 'full')


# cases
quan_controls_full <- getBatch(quan, cases = F, data_size = 'full')
funnorm_controls_full <- getBatch(funnorm, cases = F, data_size = 'full')
raw_controls_full <- getBatch(raw, cases = F, data_size = 'full')

##########
# sub data
##########

# cases
quan_cases_sub <- getBatch(quan, cases = T, data_size = 'sub')
funnorm_cases_sub <- getBatch(funnorm, cases = T, data_size = 'sub')
raw_cases_sub <- getBatch(raw, cases = T, data_size = 'sub')

# cases
quan_controls_sub <- getBatch(quan, cases = F, data_size = 'sub')
funnorm_controls_sub <- getBatch(funnorm, cases = F, data_size = 'sub')
raw_controls_sub <- getBatch(raw, cases = F, data_size = 'sub')

# remove unneeded data
rm(quan, funnorm, raw)

# remove extra column in everything but controls full
quan_cases_full <- quan_cases_full[, colnames(quan_cases_full) != 'cg13869341']
funnorm_cases_full <- funnorm_cases_full[, colnames(funnorm_cases_full) != 'cg13869341']
raw_cases_full <- raw_cases_full[, colnames(raw_cases_full) != 'cg13869341']

# sub
quan_cases_sub <- quan_cases_sub[, colnames(quan_cases_sub) != 'cg13869341']
funnorm_cases_sub <- funnorm_cases_sub[, colnames(funnorm_cases_sub) != 'cg13869341']
raw_cases_sub <- raw_cases_sub[, colnames(raw_cases_sub) != 'cg13869341']

# controls sub
quan_controls_sub <- quan_controls_sub[, colnames(quan_controls_sub) != 'cg13869341']
funnorm_controls_sub <- funnorm_controls_sub[, colnames(funnorm_controls_sub) != 'cg13869341']
raw_controls_sub <- raw_controls_sub[, colnames(raw_controls_sub) != 'cg13869341']

save.image('/home/benbrew/Desktop/temp_quan_funnorm.RData')
# load('/home/benbrew/Desktop/temp_quan_funnorm.RData')

##########
# rerun pca and cases and controls
##########

##############################################################
##########
# quan
##########

# full
# cases
getPCA(quan_cases_full, 
       'gender', 
       'PCA quan cases gender',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(quan_controls_full, 
       'gender', 
       'PCA quan controls gender',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)


#cases
getPCA(quan_cases_full, 
       'sen_batch', 
       'PCA quan cases sen_batch',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(quan_controls_full, 
       'sen_batch', 
       'PCA quan controls sen_batch',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)

#cases
getPCA(quan_cases_full, 
       'sentrix_id', 
       'PCA quan cases sentrix_id',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(quan_controls_full, 
       'sentrix_id', 
       'PCA quan controls sentrix_id',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)




# sub
# cases
getPCA(quan_cases_sub, 
       'gender', 
       'PCA quan cases gender',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(quan_controls_sub, 
       'gender', 
       'PCA quan controls gender',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)


#cases
getPCA(quan_cases_sub, 
       'sen_batch', 
       'PCA quan cases sen_batch',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(quan_controls_sub, 
       'sen_batch', 
       'PCA quan controls sen_batch',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)

#cases
getPCA(quan_cases_sub, 
       'sentrix_id', 
       'PCA quan cases sentrix_id',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(quan_controls_sub, 
       'sentrix_id', 
       'PCA quan controls sentrix_id',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)



##############################################################
##########
# funnorm
##########

# full
# cases
getPCA(funnorm_cases_full, 
       'gender', 
       'PCA funnorm cases gender',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(funnorm_controls_full, 
       'gender', 
       'PCA funnorm controls gender',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)


#cases
getPCA(funnorm_cases_full, 
       'sen_batch', 
       'PCA funnorm cases sen_batch',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(funnorm_controls_full, 
       'sen_batch', 
       'PCA funnorm controls sen_batch',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)

#cases
getPCA(funnorm_cases_full, 
       'sentrix_id', 
       'PCA funnorm cases sentrix_id',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(funnorm_controls_full, 
       'sentrix_id', 
       'PCA funnorm controls sentrix_id',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)




# sub
# cases
getPCA(funnorm_cases_sub, 
       'gender', 
       'PCA funnorm cases gender',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(funnorm_controls_sub, 
       'gender', 
       'PCA funnorm controls gender',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)


#cases
getPCA(funnorm_cases_sub, 
       'sen_batch', 
       'PCA funnorm cases sen_batch',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(funnorm_controls_sub, 
       'sen_batch', 
       'PCA funnorm controls sen_batch',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)

#cases
getPCA(funnorm_cases_sub, 
       'sentrix_id', 
       'PCA funnorm cases sentrix_id',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(funnorm_controls_sub, 
       'sentrix_id', 
       'PCA funnorm controls sentrix_id',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)


##############################################################
##########
# raw
##########

# full
# cases
getPCA(raw_cases_full, 
       'gender', 
       'PCA raw cases gender',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(raw_controls_full, 
       'gender', 
       'PCA raw controls gender',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)


#cases
getPCA(raw_cases_full, 
       'sen_batch', 
       'PCA raw cases sen_batch',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(raw_controls_full, 
       'sen_batch', 
       'PCA raw controls sen_batch',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)

#cases
getPCA(raw_cases_full, 
       'sentrix_id', 
       'PCA raw cases sentrix_id',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(raw_controls_full, 
       'sentrix_id', 
       'PCA raw controls sentrix_id',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)




# sub
# cases
getPCA(raw_cases_sub, 
       'gender', 
       'PCA raw cases gender',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(raw_controls_sub, 
       'gender', 
       'PCA raw controls gender',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)


#cases
getPCA(raw_cases_sub, 
       'sen_batch', 
       'PCA raw cases sen_batch',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(raw_controls_sub, 
       'sen_batch', 
       'PCA raw controls sen_batch',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)

#cases
getPCA(raw_cases_sub, 
       'sentrix_id', 
       'PCA raw cases sentrix_id',
       gene_start = 9,
       case_control = 'cases',
       pca1 = 1,
       pca2 = 2)

#controls
getPCA(raw_controls_sub, 
       'sentrix_id', 
       'PCA raw controls sentrix_id',
       gene_start = 9,
       case_control = 'controls',
       pca1 = 1,
       pca2 = 2)




##########
# save non batch corrected
##########
##########
# quan
##########

# full 

# save cases 
saveRDS(quan_cases_full, paste0(model_data, '/quan_cases_full.rda'))

# save contorls
saveRDS(quan_controls_full, paste0(model_data, '/quan_controls_full.rda'))

# sub

# save cases 
saveRDS(quan_cases_sub, paste0(model_data, '/quan_cases_sub.rda'))

# save contorls
saveRDS(quan_controls_sub, paste0(model_data, '/quan_controls_sub.rda'))

#########
# funnorm
##########
# full 

# save cases 
saveRDS(funnorm_cases_full, paste0(model_data, '/funnorm_cases_full.rda'))

# save contorls
saveRDS(funnorm_controls_full, paste0(model_data, '/funnorm_controls_full.rda'))

# sub

# save cases 
saveRDS(funnorm_cases_sub, paste0(model_data, '/funnorm_cases_sub.rda'))

# save contorls
saveRDS(funnorm_controls_sub, paste0(model_data, '/funnorm_controls_sub.rda'))


#########
# raw
##########
# full 

# save cases 
saveRDS(raw_cases_full, paste0(model_data, '/raw_cases_full.rda'))

# save contorls
saveRDS(raw_controls_full, paste0(model_data, '/raw_controls_full.rda'))

# sub

# save cases 
saveRDS(raw_cases_sub, paste0(model_data, '/raw_cases_sub.rda'))

# save contorls
saveRDS(raw_controls_sub, paste0(model_data, '/raw_controls_sub.rda'))

