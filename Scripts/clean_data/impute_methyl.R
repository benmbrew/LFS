####### Script will load 4 preprocessing methods for cases and controls and impute with knn
# this is 3nd step in pipeline

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
# load idat files to be imputed with knn
load(paste0(methyl_data, '/raw.RData'))
load(paste0(methyl_data, '/swan.RData'))
load(paste0(methyl_data, '/quan.RData'))
load(paste0(methyl_data, '/funnorm.RData'))
rm(id_map, id_map_control, rgSetList, rgSetListControls)


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
beta_swan <- keepCG(beta_swan)
beta_quan <- keepCG(beta_quan)
beta_funnorm <- keepCG(beta_funnorm)

# apply to controls
beta_raw_controls <- keepCG(beta_raw_controls)
beta_swan_controls <- keepCG(beta_swan_controls)
beta_quan_controls <- keepCG(beta_quan_controls)
beta_funnorm_controls <- keepCG(beta_funnorm_controls)

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
beta_swan <- avgBeta(beta_swan)
beta_quan <- avgBeta(beta_quan)
beta_funnorm <- avgBeta(beta_funnorm)

# apply to controls
beta_raw_controls <- avgBeta(beta_raw_controls)
beta_swan_controls <- avgBeta(beta_swan_controls)
beta_quan_controls <- avgBeta(beta_quan_controls)
beta_funnorm_controls <- avgBeta(beta_funnorm_controls)


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
  
  # join rownames and methyl_impute and then erase rownames
  data_knn <- cbind(id = rownames(data_knn), data_knn)
  # rownames(data_knn) <- NULL
  # join rownames and methyl_impute and then erase rownames

  return(data_knn)
  
}
save.image(paste0(methyl_data, '/imputed_others.RData'))


# Remove swan, quan, and funnorm because they dont have missing values
rm(beta_swan, beta_quan, beta_funnorm, beta_swan_controls, beta_quan_controls, beta_funnorm_controls)
##########
# apply to cases
##########
beta_raw <- imputeIdat(beta_raw)#303
# beta_swan <- imputeIdat(beta_swan)#0
# beta_quan <- imputeIdat(beta_quan)#0
# beta_funnorm <- imputeIdat(beta_funnorm)#0 

##########
# apply to controls
##########
beta_raw_controls <- imputeIdat(beta_raw_controls)#552
# beta_swan_controls <- imputeIdat(beta_swan_controls)#0
# beta_quan_controls <- imputeIdat(beta_quan_controls)#0
# beta_funnorm_controls <- imputeIdat(beta_funnorm_controls)#0

# dont save illumina because it has too many missing values
save.image(paste0(methyl_data, '/imputed_raw.RData'))


# #########################################################################################
# # # Read in methylation data
# # methyl_gene <- read.csv(paste0(methyl_data, '/methyl.csv'), stringsAsFactors = FALSE)
# # 
# # # read in probe level methylation 
# # methyl_probe <- read.csv(paste0(methyl_data, '/methylation.csv'), header = TRUE, check.names = FALSE)
# 
# #################
# # function that imputes both knn and lsa on methylation data
# #################
# imputeMethyl <- function(data,
#                          probe){
#   
#   
#   if (probe) {
#     
#     # transpose probe methylation
#     col_names <- data$Probe
#     data <- as.data.frame(t(data))
#     names(data)<- col_names
#     data <- data[!duplicated(rownames(data)),]
#     data <- cbind(x = rownames(data), data) 
#     data <- data[2:nrow(data),]
#     names(data)[1] <- 'id'
#     rownames(data) <- NULL
#     data[, 2:ncol(data)] <-
#       apply(data[,2:ncol(data)], 2, function(x){as.numeric(as.character(x))})
#     data <- as.data.frame(data)
#     
#     # drop duplicates from methylation so LSA work
#     data <- data[!duplicated(data$id),]
#     data <- data[!is.na(data$id),]
#     
#   }
#   
#   # put ids in rownames for imputation
#   rownames(data) <- data[,1]
#   data <- data[, -1]
#   
#   # remove duplicate rownames,so it can impute
#   data <- data[!grepl('y', rownames(data)),]
#   data <- as.matrix(data)
#   
#   # run lsaImputaion of methylation data
#   data_lsa <- lsaImputation(incomplete_data = data, sample_rows = TRUE)
#   # impute with knn as well
#   data_knn <- knnImputation(data, sample_rows = TRUE)
#   
#   # join rownames and methyl_impute and then erase rownames
#   data_lsa <- cbind(id = rownames(data_lsa), data_lsa)
#   rownames(data_lsa) <- NULL
#   
#   # join rownames and methyl_impute and then erase rownames
#   data_knn <- cbind(id = rownames(data_knn), data_knn)
#   rownames(data_knn) <- NULL
#   
#   return(list(data_lsa, data_knn))
#   
# }
# 
# # run function on gene data
# gene <- imputeMethyl(methyl_gene,
#                      probe = F)
# #extract data
# gene_lsa <- gene[[1]]
# gene_knn <- gene[[2]]
# 
# # write csv
# write.csv(gene_lsa, paste0(imputed_data, '/methyl_impute_gene_lsa.csv'))
# write.csv(gene_knn, paste0(imputed_data, '/methyl_impute_gene_knn.csv'))
# 
# # run function on probe data
# probe <- imputeMethyl(methyl_probe,
#                       probe = T)
# 
# #extract data
# probe_lsa <- probe[[1]]
# probe_knn <- probe[[2]]
# 
# # write csv
# write.csv(probe_lsa, paste0(imputed_data, '/methyl_impute_probe_lsa.csv'))
# write.csv(probe_knn, paste0(imputed_data, '/methyl_impute_probe_knn.csv'))
# 
# save.image(paste0(imputed_data, '/imputed_gene_probe.RData'))


