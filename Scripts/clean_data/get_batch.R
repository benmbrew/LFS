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
# create a batch variable for containing 9721365183 or not (indicator for if it was run at sickkidsp)
##########

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
                   pca2) 
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
# cases
getPCA(quan, 
       'gender', 
       'PCA quan cases gender',
       gene_start = 10,
       cases = T,
       controls = F,
       pca1 = 1,
       pca2 = 2)

# cases
getPCA(quan, 
       'gender', 
       'PCA quan cases gender',
       gene_start = 10,
       cases = T,
       controls = F,
       pca1 = 2,
       pca2 = 3)

# cases
getPCA(quan, 
       'batch', 
       'PCA quan cases batch',
       gene_start = 10,
       cases = T,
       controls = F,
       pca1 = 1,
       pca2 = 2)

# cases
getPCA(quan, 
       'batch', 
       'PCA quan cases batch',
       gene_start = 10,
       cases = T,
       controls = F,
       pca1 = 2,
       pca2 = 3)




# controls
getPCA(quan, 
       'gender', 
       'PCA quan controls gender',
       gene_start = 10,
       cases = F,
       controls = T,
       pca1 = 1,
       pca2 = 2)


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
# fix gender batch in cases and controls.
##########
data<- quan

getBatch <- function(data, cases, batch)
{
  
  # make full just two batches
  if (cases) {
    data <- data[data$type == 'cases',]
    
    
  } else {
    data <- data[data$type == 'controls',]
  }
  
  if(batch == 'gender') {
    batch_indicator <- as.factor(data$gender)
  }
  
  if(batch == 'batch') {
    batch_indicator <- as.factor(data$batch)
    
  }
  
  if(batch == 'sentrix_id') {
    batch_indicator <- as.factor(data$sentrix_id)
    
  }
  
  gender <- data$gender
  sentrix_id <- data$sentrix_id
  batch <- data$batch
  type <- data$type
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
  final_dat$batch <- batch
  final_dat$ids <- ids
  final_dat$type <- type
  final_dat$age_sample_collection <- sample_collection
  final_dat$age_diagnosis <- diagnosis
  final_dat$p53_germline <- p53
  final_dat$sentrix_id <- sentrix_id
  final_dat$cancer_diagnosis_diagnoses <- cancer_diagnosis
  final_dat <- final_dat[, c('ids', 'type', 'gender', 'age_sample_collection',
                             'age_diagnosis','cancer_diagnosis_diagnoses' ,'p53_germline', 'batch', 'sentrix_id',features)]
  rownames(final_dat) <- NULL
  
  return(final_dat)
  
  
}

# get non batch corrected cases and controls
quan_cases <- quan[quan$type == 'cases',]
quan_controls <- quan[quan$type == 'controls',]


# get batch corrected cases and controls for gender
quan_cases_gen <- getBatch(quan, cases = T, batch = 'gender')
quan_controls_gen <- getBatch(quan, cases = F, batch = 'gender')

# get batch corrected data for SAM and not
quan_cases_sam <-  getBatch(quan, cases = T, batch = 'batch')
quan_cases_sen <-  getBatch(quan, cases = T, batch = 'sentrix_id')

# apply gender correction to sam and sentrix
quan_cases_sam_gen <-getBatch(quan_cases_sam, cases = T, batch = 'gender')
quan_cases_sen_gen <-getBatch(quan_cases_sen, cases = T, batch = 'gender')



##########
# rerun pca and cases and controls
##########

# cases
getPCA(quan_cases, 
       'gender', 
       'PCA quan cases gender no batch',
       gene_start = 10,
       cases = T,
       controls = F,
       pca1 = 1, 
       pca2 = 2)

# cases
getPCA(quan_cases_gen, 
       'gender', 
       'PCA quan cases gender gender batch',
       gene_start = 10,
       cases = T,
       controls = F,
       pca1 = 1, 
       pca2 = 2)

# cases
getPCA(quan_cases_sen_gen, 
       'gender', 
       'PCA quan cases gender and sentrix batch',
       gene_start = 10,
       cases = T,
       controls = F,
       pca1 = 1, 
       pca2 = 2)

# cases
getPCA(quan_cases_sam_gen, 
       'gender', 
       'PCA quan cases gender and sam batch',
       gene_start = 10,
       cases = T,
       controls = F,
       pca1 = 1, 
       pca2 = 2)



##########
# save non batch corrected
##########

# save cases 
saveRDS(quan_cases, paste0(model_data, '/quan_cases.rda'))

# save contorls
saveRDS(quan_controls, paste0(model_data, '/quan_controls.rda'))

##########
# save batch corrected data for gender
##########

# save cases 
saveRDS(quan_cases_gen, paste0(model_data, '/quan_cases_gen.rda'))

# save contorls
saveRDS(quan_controls_gen, paste0(model_data, '/quan_controls_gen.rda'))


##########
# save batch corrected data for sentrix id and SAM
##########

# save cases 
saveRDS(quan_cases_sen, paste0(model_data, '/quan_cases_sen.rda'))

# save contorls
saveRDS(quan_controls_sam, paste0(model_data, '/quan_controls_sam.rda'))

##########
# save batch corrected data for sentrix id and SAM and gender!
##########

# save cases 
saveRDS(quan_cases_sen_gen, paste0(model_data, '/quan_cases_sen_gen.rda'))

# save contorls
saveRDS(quan_controls_sam_gen, paste0(model_data, '/quan_controls_sam_gen.rda'))
