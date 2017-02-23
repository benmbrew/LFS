####### Script will load 4 preprocessing methods for cases and controls and impute with knn
# this is 3rd step in pipeline

##########
# initialize libraries
##########
library(dplyr)
library(stringr)
library(impute)
library(data.table)


##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/Impute')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
idat_data <- paste0(methyl_data, '/raw_files')
imputed_data <- paste0(data_folder, '/imputed_data')
results_folder <- paste0(test, '/Results')

##########
# construct files for imputation
##########
results_file <- "file.txt"
constructPath <- function(intermediate_folders, parent=results_folder,
                          file=results_file) {
  paste(parent, intermediate_folders, file, sep="/")
}

incomplete_file <- constructPath("Incomplete")
imputed_file <- constructPath("Imputed")
jvmGBLimit <- 8

##########
# load knn and lsa (now just using knn)
##########
source(paste0(project_folder, '/Code/Functions/lsaImputation.R'))
source(paste0(project_folder, '/Code/Functions/knnImputation.R'))

##########
# Read in methylation probe and gene
##########
beta_raw <- readRDS(paste0(methyl_data, '/beta_raw.rda'))
beta_raw_controls <- readRDS(paste0(methyl_data, '/beta_raw_controls.rda'))

##########
# function to remove non cg sites 
##########
keepCG <- function(data) {
  
  features <- colnames(data)[1:(length(data) -1)]
  data <- data[, c('ids', features)]
  keep <- names(data)[-1][grepl('cg', names(data)[-1])]
  data <- data[, c('ids', keep)]
  return(data)
  
}

# apply to cases
beta_raw <- keepCG(beta_raw)

# apply to controls
beta_raw_controls <- keepCG(beta_raw_controls)

##########
# function that takes the avgs of duplicates, as they are controls
##########
avgBeta <- function(data){
  
  duplicate_ids <- data$ids[duplicated(data$ids)]
  mean_values <- list()
  for (ids in duplicate_ids) {
    sub_dat <- data[data$ids == ids,]
    mean_values[[ids]] <- apply(sub_dat[,2:ncol(sub_dat)], 2, mean)
    print(ids)
  }
  
  mean_val <- do.call(rbind, mean_values)
  mean_val <- as.data.frame(mean_val)
  mean_val$ids <- rownames(mean_val)
  
  duplicate_ids <- paste0('^', duplicate_ids, '$')
  duplicate_ids <- paste(duplicate_ids, collapse = '|')
  data <- data[!grepl(duplicate_ids, data$ids),]
  
  data <- rbind(data, mean_val)
  
  return(data)
  
}

# apply to cases
beta_raw <- avgBeta(beta_raw)

# apply to controls
beta_raw_controls <- avgBeta(beta_raw_controls)

##########
# function that imputes just knn on original idat data
##########

imputeIdat <- function (data) {
  
  # make data a matrix and put ids in rownames, remove id column
  rownames(data) <- data$ids
  data$ids <- NULL
  data <- as.matrix(data)
  
  # # run lsaImputaion of methylation data
  # data_lsa <- lsaImputation(incomplete_data = data, sample_rows = TRUE)
  # impute with knn as well
  data_knn <- knnImputation(data, sample_rows = TRUE)
  
  # # join rownames and methyl_impute and then erase rownames
  # data_lsa <- cbind(id = rownames(data_lsa), data_lsa)
  # rownames(data_lsa) <- NULL
  # 
  data_knn <- as.data.frame(data_knn)
  features <- colnames(data_knn)
  data_knn$id <- rownames(data_knn)
  data_knn <- data_knn[, c('id', features)]
  # join rownames and methyl_impute and then erase rownames
  # rownames(data_knn) <- NULL
  # join rownames and methyl_impute and then erase rownames

  return(data_knn)
  
}

##########
# apply to cases
##########
beta_raw <- imputeIdat(beta_raw)#303

##########
# apply to controls
##########
beta_raw_controls <- imputeIdat(beta_raw_controls)#552

#########
# save data
#########
saveRDS(beta_raw, paste0(methyl_data, '/beta_raw.rda'))

saveRDS(beta_raw_controls, paste0(methyl_data, '/beta_raw_controls.rda'))



