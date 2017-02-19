####### Script will correct for batches by gender for cases and controls, and combine and correct for diff tech
# # This is 5th step (B)

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
beta_raw <- readRDS(paste0(methyl_data, '/beta_raw.rda'))
beta_raw_controls <- readRDS(paste0(methyl_data, '/beta_raw_controls.rda'))

beta_swan <- readRDS(paste0(methyl_data, '/beta_swan.rda'))
beta_swan_controls <- readRDS(paste0(methyl_data, '/beta_swan_controls.rda'))

beta_quan <- readRDS(paste0(methyl_data, '/beta_quan.rda'))
beta_quan_controls <- readRDS(paste0(methyl_data, '/beta_quan_controls.rda'))

beta_funnorm <- readRDS(paste0(methyl_data, '/beta_funnorm.rda'))
beta_funnorm_controls <- readRDS(paste0(methyl_data, '/beta_funnorm_controls.rda'))

##########
# remove column at end of swan, quan, and funnorm
##########
beta_swan$id.1 <- NULL
beta_quan$id.1 <- NULL
beta_funnorm$id.1 <- NULL

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

beta_raw <- getModData(beta_raw)

beta_quan <- getModData(beta_quan)

beta_swan <- getModData(beta_swan)

beta_funnorm <- getModData(beta_funnorm)

##########
# get intersection of controls and normal
##########
# data <- beta_raw
# data_controls <- beta_raw_controls
getFeatInt <- function(data, data_controls)
{
  features <- names(data)[14:ncol(data)]
  features_controls <- names(data_controls)[14:ncol(data_controls)]
  
  data$type <- 'cases'
  data_controls$type <- 'controls'
  
  data$type <- as.factor(data$type)
  data_controls$type <- as.factor(data_controls$type)
  
  
  intersect_features <- intersect(features, features_controls)
  
  data <- data[, c('id', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 'age_sample_collection' ,
                   'gdna.exon.intron', 'gdna.base.change', 'protein.codon.change',  'protein.codon.num', 
                   'splice.delins.snv', 'codon72.npro' ,'mdm2.nG','gender','type', intersect_features)]
  
  data_controls <- data_controls[, c('id', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 'age_sample_collection' ,
                                     'gdna.exon.intron', 'gdna.base.change', 'protein.codon.change',  'protein.codon.num', 
                                     'splice.delins.snv', 'codon72.npro' ,'mdm2.nG','gender','type', intersect_features)]
  
  data <- rbind(data, data_controls)
  
  return(data)
  
}

raw <- getFeatInt(beta_raw, beta_raw_controls)

quan <- getFeatInt(beta_quan, beta_quan_controls)

swan <- getFeatInt(beta_swan, beta_swan_controls)

funnorm <- getFeatInt(beta_funnorm, beta_funnorm_controls)


########## 
# PCA of each data type and cases vs controls
##########
pca_data <- funnorm_cases_batch
column_name <- 'batch'
gene_start <- 4
# functionp needs to take a clinical column, remove others, and plot pcas
getPCA <- function(pca_data, column_name, name, gene_start, controls, cases) 
{
  
  if (controls) {
    
    pca_data <- pca_data[pca_data$type == 'controls',]
    rownames(pca_data) <-pca_data$id
    
  } else if (cases) {
    
    pca_data <- pca_data[pca_data$type == 'cases',]
    rownames(pca_data) <-pca_data$id
    
    
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
# use type
##########
# raw
getPCA(raw, 
      'type', 
      'PCA raw type', 
      gene_start = 15, 
      cases = F, 
      controls = F)

# swan
getPCA(swan, 
       'type', 
       'PCA swan type', 
       gene_start = 15,
       cases = F, 
       controls = F)

# quan
getPCA(quan, 
       'type', 
       'PCA quan type', 
       gene_start = 15,
       cases = F, 
       controls = F)

# funnorm
getPCA(funnorm, 
       'type', 
       'PCA funnorm type', 
       gene_start = 15,
       cases = F, 
       controls = F)


##########
# combine cases and controls and plot based on that
##########

# cases
getPCA(raw, 
       'gender', 
       'PCA raw cases gender',
        gene_start = 15,
        cases = T,
        controls = F)


getPCA(quan, 
       'gender', 
       'PCA quan cases gender', 
       gene_start = 15,
       cases = T,
       controls = F)

getPCA(swan, 
       'gender', 
       'PCA swan cases gender', 
       gene_start = 15,
       cases = T,
       controls = F)


getPCA(funnorm, 
       'gender', 
       'PCA funnorm cases gender', 
       gene_start = 15,
       cases = T,
       controls = F)

# controls
getPCA(raw, 
       'gender', 
       'PCA raw controls gender',
       gene_start = 15,
       cases = F,
       controls = T)


getPCA(quan, 
       'gender', 
       'PCA quan controls gender', 
       gene_start = 15,
       cases = F,
       controls = T)

getPCA(swan, 
       'gender', 
       'PCA swan controls gender', 
       gene_start = 15,
       cases = F,
       controls = T)


getPCA(funnorm, 
       'gender', 
       'PCA funnorm controls gender', 
       gene_start = 15,
       cases = F,
       controls = T)

##########
# remove outliers  4257 cases, 3391, 3392 controls
##########
removeOutlier <- function(data, funnorm) {
  
  if(funnorm) {
    data <- data[data$id != '3646',]
    
  }
  # cases outlier
  data <- data[data$id != '4257',]
  
  #controls outlier
  data <- data[data$id != '3391',]
  data <- data[data$id != '3392',]
  return(data)
}

raw <- removeOutlier(raw, funnorm = F)
swan <- removeOutlier(swan, funnorm = F)
quan <- removeOutlier(quan, funnorm = F)
funnorm <- removeOutlier(funnorm, funnorm = F)

##########
# rerun pca
##########
# cases
getPCA(raw, 
       'gender', 
       'PCA raw cases gender',
       gene_start = 15,
       cases = T,
       controls = F)


getPCA(quan, 
       'gender', 
       'PCA quan cases gender', 
       gene_start = 15,
       cases = T,
       controls = F)

getPCA(swan, 
       'gender', 
       'PCA swan cases gender', 
       gene_start = 15,
       cases = T,
       controls = F)


getPCA(funnorm, 
       'gender', 
       'PCA funnorm cases gender', 
       gene_start = 15,
       cases = T,
       controls = F)

# controls
getPCA(raw, 
       'gender', 
       'PCA raw controls gender',
       gene_start = 15,
       cases = F,
       controls = T)


getPCA(quan, 
       'gender', 
       'PCA quan controls gender', 
       gene_start = 15,
       cases = F,
       controls = T)

getPCA(swan, 
       'gender', 
       'PCA swan controls gender', 
       gene_start = 15,
       cases = F,
       controls = T)


getPCA(funnorm, 
       'gender', 
       'PCA funnorm controls gender', 
       gene_start = 15,
       cases = F,
       controls = T)


##########
# fix gender batch in cases and controls.
##########
data <- raw
getBatch <- function(data, cases)
{
  
  # make full just two batches
  if (cases) {
    data$gender <- ifelse(data$gender == 'M', 'a', 'b')
    data <- data[data$type == 'cases',]
    
    
  } else {
    data$gender <- ifelse(data$gender == 'M', 'a', 'b')
    data <- data[data$type == 'controls',]
  }
    
  batch <- as.factor(data$gender)
  type <- data$type
  id <- data$id
  sample_collection <- data$age_sample_collection
  diagnosis <- data$age_diagnosis
  p53 <- data$p53_germline
  cancer_diagnosis <- data$cancer_diagnosis_diagnoses
  # put model ids in rownames and remove columns
  mat_data <- data[, 15:ncol(data)]
  # get features 
  features <- colnames(mat_data)
  mat_data <- t(mat_data)
  
  # get intercept
  modcombat <- model.matrix(~1, data = data)
  combat <- ComBat(dat = mat_data, batch = batch, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  # transpose and add back columns
  final_dat <- as.data.frame(t(combat))
  final_dat$model_id <- rownames(final_dat)
  final_dat$batch <- batch
  final_dat$id <- id
  final_dat$type <- type
  final_dat$age_sample_collection <- sample_collection
  final_dat$age_diagnosis <- diagnosis
  final_dat$p53_germline <- p53
  final_dat$cancer_diagnosis_diagnoses <- cancer_diagnosis
  final_dat <- final_dat[, c('id', 'type', 'batch', 'age_sample_collection',
                             'age_diagnosis','cancer_diagnosis_diagnoses' ,'p53_germline',features)]
  rownames(final_dat) <- NULL
  
  return(final_dat)

  
}

# raw
raw_cases_batch <- getBatch(raw, cases = T)
raw_controls_batch <- getBatch(raw, cases = F)

# quan
quan_cases_batch <- getBatch(quan, cases = T)
quan_controls_batch <- getBatch(quan, cases = F)

# swan
swan_cases_batch <- getBatch(swan, cases = T)
swan_controls_batch <- getBatch(swan, cases = F)

# funnorm
funnorm_cases_batch <- getBatch(funnorm, cases = T)
funnorm_controls_batch <- getBatch(funnorm, cases = F)

save.image('/home/benbrew/Desktop/batch_temp.RData')
# load('/home/benbrew/Desktop/batch_temp.RData')

##########
# rerun pca and cases and controls
##########

# cases
getPCA(raw_cases_batch, 
       'batch', 
       'PCA raw cases gender batch',
       gene_start = 7,
       cases = T,
       controls = F)


getPCA(quan_cases_batch, 
       'batch', 
       'PCA quan cases gender batch', 
       gene_start = 7,
       cases = T,
       controls = F)

getPCA(swan_cases_batch, 
       'batch', 
       'PCA swan cases gender batch', 
       gene_start = 7,
       cases = T,
       controls = F)


getPCA(funnorm_cases_batch, 
       'batch', 
       'PCA funnorm cases gender batch', 
       gene_start = 7,
       cases = T,
       controls = F)

# controls
getPCA(raw_controls_batch, 
       'batch', 
       'PCA raw controls gender batch',
       gene_start = 7,
       cases = F,
       controls = T)


getPCA(quan_controls_batch, 
       'batch', 
       'PCA quan controls gender batch', 
       gene_start = 7,
       cases = F,
       controls = T)

getPCA(swan_controls_batch, 
       'batch', 
       'PCA swan controls gender batch', 
       gene_start = 7,
       cases = F,
       controls = T)


getPCA(funnorm_controls_batch, 
       'batch', 
       'PCA funnorm controls gender batch', 
       gene_start = 7,
       cases = F,
       controls = T)

##########
# remove outlier from funnorm and swan
##########
funnorm_cases_batch <- removeOutlier(funnorm_cases_batch, funnorm = T)
swan_cases_batch <- removeOutlier(swan_cases_batch, funnorm = T)


##########
# save data
##########

# save cases 
saveRDS(raw_cases_batch, paste0(model_data, '/raw_cases_batch.rda'))
saveRDS(swan_cases_batch, paste0(model_data, '/swan_cases_batch.rda'))
saveRDS(quan_cases_batch, paste0(model_data, '/quan_cases_batch.rda'))
saveRDS(funnorm_cases_batch, paste0(model_data, '/funnorm_cases_batch.rda'))

# save contorls
saveRDS(raw_controls_batch, paste0(model_data, '/raw_controls_batch.rda'))
saveRDS(swan_controls_batch, paste0(model_data, '/swan_controls_batch.rda'))
saveRDS(quan_controls_batch, paste0(model_data, '/quan_controls_batch.rda'))
saveRDS(funnorm_controls_batch, paste0(model_data, '/funnorm_controls_batch.rda'))


