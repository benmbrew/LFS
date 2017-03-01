####### Script will correct for batches by gender for cases and controls, and combine and correct for diff tech
# # This is 4th step

##########
# initialize libraries
##########
library(dplyr)
library(sva)

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
  features <- names(data)[7:ncol(data)]
  features_controls <- names(data_controls)[7:ncol(data_controls)]
  
  data$type <- 'cases'
  data_controls$type <- 'controls'
  
  data$type <- as.factor(data$type)
  data_controls$type <- as.factor(data_controls$type)
  
  
  intersect_features <- intersect(features, features_controls)
  
  data <- data[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                   'age_sample_collection' ,'gender','type', intersect_features)]
  
  data_controls <- data_controls[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 
                                     'age_sample_collection' ,'gender','type', intersect_features)]
  
  data <- rbind(data, data_controls)
  
  return(data)
  
}

quan <- getFeatInt(beta_quan, beta_quan_controls)

########## 
# PCA of each data type and cases vs controls
##########
# functionp needs to take a clinical column, remove others, and plot pcas
getPCA <- function(pca_data, column_name, name, gene_start, controls, cases) 
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
  temp <- pca_data[,2:data_length]
  
  # plot data
  #fill in factors with colors 
  col_vec <- c('red', 'green', 'blue', 'orange', 'black', 'orange')
  colors <- col_vec[pca_data[, column_name]]
  
  
  plot <- plot(pca$x[,1], 
               pca$x[,2],
               xlab = 'pca 1',
               ylab = 'pca 2',
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
# cases
getPCA(quan, 
       'gender', 
       'PCA quan cases gender',
       gene_start = 8,
       cases = T,
       controls = F)


# controls
getPCA(quan, 
       'gender', 
       'PCA quan controls gender',
       gene_start = 8,
       cases = F,
       controls = T)


##########
# remove outliers  4257 cases, 3391, 3392 controls
##########
removeOutlier <- function(data, funnorm) {
  

  #controls outlier
  data <- data[data$ids != '3391',]
  data <- data[data$ids != '3392',]
  return(data)
}

quan <- removeOutlier(quan, funnorm = F)

##########
# rerun pca
##########
# cases
getPCA(quan, 
       'gender', 
       'PCA quan cases gender',
       gene_start = 8,
       cases = T,
       controls = F)


# controls
getPCA(quan, 
       'gender', 
       'PCA quan controls gender',
       gene_start = 8,
       cases = F,
       controls = T)


##########
# fix gender batch in cases and controls.
##########
getBatch <- function(data, cases)
{
  
  # make full just two batches
  if (cases) {
    data <- data[data$type == 'cases',]
    
    
  } else {
    data <- data[data$type == 'controls',]
  }
  
  batch <- as.factor(data$gender)
  type <- data$type
  ids <- data$ids
  sample_collection <- data$age_sample_collection
  diagnosis <- data$age_diagnosis
  p53 <- data$p53_germline
  cancer_diagnosis <- data$cancer_diagnosis_diagnoses
  # put model ids in rownames and remove columns
  mat_data <- data[, 8:ncol(data)]
  # get features 
  features <- colnames(mat_data)
  mat_data <- t(mat_data)
  
  # get intercept
  modcombat <- model.matrix(~1, data = data)
  combat <- ComBat(dat = mat_data, batch = batch, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  # transpose and add back columns
  final_dat <- as.data.frame(t(combat))
  final_dat$model_id <- rownames(final_dat)
  final_dat$gender <- batch
  final_dat$ids <- ids
  final_dat$type <- type
  final_dat$age_sample_collection <- sample_collection
  final_dat$age_diagnosis <- diagnosis
  final_dat$p53_germline <- p53
  final_dat$cancer_diagnosis_diagnoses <- cancer_diagnosis
  final_dat <- final_dat[, c('ids', 'type', 'gender', 'age_sample_collection',
                             'age_diagnosis','cancer_diagnosis_diagnoses' ,'p53_germline', features)]
  rownames(final_dat) <- NULL
  
  return(final_dat)
  
  
}

# get non batch corrected cases and controls
quan_cases <- quan[quan$type == 'cases',]
quan_controls <- quan[quan$type == 'controls',]


# get batch corrected cases and controls
quan_cases_batch <- getBatch(quan, cases = T)
quan_controls_batch <- getBatch(quan, cases = F)


##########
# rerun pca and cases and controls
##########

# cases
getPCA(quan_cases, 
       'gender', 
       'PCA quan cases gender no batch',
       gene_start = 8,
       cases = T,
       controls = F)


# controls
getPCA(quan_controls, 
       'gender', 
       'PCA quan controls gender no batch',
       gene_start = 8,
       cases = F,
       controls = T)

# cases
getPCA(quan_cases_batch, 
       'gender', 
       'PCA quan cases gender batch',
       gene_start = 8,
       cases = T,
       controls = F)


# controls
getPCA(quan_controls_batch, 
       'gender', 
       'PCA quan controls gender batch',
       gene_start = 8,
       cases = F,
       controls = T)



##########
# save batch corrected data
##########

# save cases 
saveRDS(quan_cases_batch, paste0(model_data, '/quan_cases_batch.rda'))

# save contorls
saveRDS(quan_controls_batch, paste0(model_data, '/quan_controls_batch.rda'))

##########
# save non batch corrected
##########

# save cases 
saveRDS(quan_cases, paste0(model_data, '/quan_cases.rda'))

# save contorls
saveRDS(quan_controls, paste0(model_data, '/quan_controls.rda'))
