###############################################
# This script will run a PCA on the LFS methylation data and plot age of diagnosis.
# existing by gene methylation data we have. 
# This is the 4th step in the pipeline
library(dplyr)
library(stringr)
library(impute)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/Impute')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
imputed_data <- paste0(data_folder, '/imputed_data')
results_folder <- paste0(test, '/Results')


# construct files for imputation
results_file <- "file.txt"
constructPath <- function(intermediate_folders, parent=results_folder,
                          file=results_file) {
  paste(parent, intermediate_folders, file, sep="/")
}

incomplete_file <- constructPath("Incomplete")
imputed_file <- constructPath("Imputed")
jvmGBLimit <- 8

# load knn and lsa
source(paste0(project_folder, '/Code/Functions/lsaImputation.R'))
source(paste0(project_folder, '/Code/Functions/knnImputation.R'))

################################################################
# Read in methylation probe and gene
################################################################
# load idat files to be imputed with knn
load(paste0(methyl_data, '/raw.RData'))
load(paste0(methyl_data, '/illum.RData'))
load(paste0(methyl_data, '/swan.RData'))
load(paste0(methyl_data, '/funnorm.RData'))
load(paste0(methyl_data, '/quan.RData'))
rm(id_map)


# Read in methylation data
methyl_gene <- read.csv(paste0(methyl_data, '/methyl.csv'), stringsAsFactors = FALSE)

# read in probe level methylation 
methyl_probe <- read.csv(paste0(methyl_data, '/methylation.csv'), header = TRUE, check.names = FALSE)

#################
# function that imputes both knn and lsa on methylation data
#################
imputeMethyl <- function(data,
                         probe){
  
  
  if (probe) {
    
    # transpose probe methylation
    col_names <- data$Probe
    data <- as.data.frame(t(data))
    names(data)<- col_names
    data <- data[!duplicated(rownames(data)),]
    data <- cbind(x = rownames(data), data) 
    data <- data[2:nrow(data),]
    names(data)[1] <- 'id'
    rownames(data) <- NULL
    data[, 2:ncol(data)] <-
      apply(data[,2:ncol(data)], 2, function(x){as.numeric(as.character(x))})
    data <- as.data.frame(data)
    
    # drop duplicates from methylation so LSA work
    data <- data[!duplicated(data$id),]
    data <- data[!is.na(data$id),]
    
  }
  
  # put ids in rownames for imputation
  rownames(data) <- data[,1]
  data <- data[, -1]
  
  # remove duplicate rownames,so it can impute
  data <- data[!grepl('y', rownames(data)),]
  data <- as.matrix(data)
  
  # run lsaImputaion of methylation data
  data_lsa <- lsaImputation(incomplete_data = data, sample_rows = TRUE)
  # impute with knn as well
  data_knn <- knnImputation(data, sample_rows = TRUE)
  
  # join rownames and methyl_impute and then erase rownames
  data_lsa <- cbind(id = rownames(data_lsa), data_lsa)
  rownames(data_lsa) <- NULL
  
  # join rownames and methyl_impute and then erase rownames
  data_knn <- cbind(id = rownames(data_knn), data_knn)
  rownames(data_knn) <- NULL
  
  return(list(data_lsa, data_knn))
  
}

# run function on gene data
gene <- imputeMethyl(methyl_gene,
                     probe = F)
#extract data
gene_lsa <- gene[[1]]
gene_knn <- gene[[2]]

# write csv
write.csv(gene_lsa, paste0(imputed_data, '/methyl_impute_gene_lsa.csv'))
write.csv(gene_knn, paste0(imputed_data, '/methyl_impute_gene_knn.csv'))

# run function on probe data
probe <- imputeMethyl(methyl_probe,
                      probe = T)

#extract data
probe_lsa <- probe[[1]]
probe_knn <- probe[[2]]

# write csv
write.csv(probe_lsa, paste0(imputed_data, '/methyl_impute_probe_lsa.csv'))
write.csv(probe_knn, paste0(imputed_data, '/methyl_impute_probe_knn.csv'))

save.image(paste0(imputed_data, '/imputed_gene_probe.RData'))


#################
# function that imputes just knn on original idat data
#################
imputeIdat <- function(data) {
  
  # put ids in rownames for imputation
  rownames(data) <- data[,1]
  data <- data[, -1]
  
  # remove duplicate rownames,so it can impute
  data <- data[!grepl('y', rownames(data)),]
  data <- as.matrix(data)
  
  # run lsaImputaion of methylation data
  data_lsa <- lsaImputation(incomplete_data = data, sample_rows = TRUE)
  # impute with knn as well
  data_knn <- knnImputation(data, sample_rows = TRUE)
  
  # join rownames and methyl_impute and then erase rownames
  data_lsa <- cbind(id = rownames(data_lsa), data_lsa)
  rownames(data_lsa) <- NULL
  
  # join rownames and methyl_impute and then erase rownames
  data_knn <- cbind(id = rownames(data_knn), data_knn)
  rownames(data_knn) <- NULL
  
  return(list(data_lsa, data_knn))
  
}


