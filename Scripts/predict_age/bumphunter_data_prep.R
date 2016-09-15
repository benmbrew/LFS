#######################################################################################
# This script will run bumphunter on methylation probes 
library(minfi)
library(bumphunter)
library(dplyr)
library(FDb.InfiniumMethylation.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(impute)
library(GenomicRanges)
library(biovizBase)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
clin_data <- paste0(data_folder, '/clin_data')

data_name <- '/Chr'

# results_file <- "file.txt"
# constructPath <- function(intermediate_folders, parent=results_folder,
#                           file=results_file) {
#   paste(parent, intermediate_folders, file, sep="/")
# }
# 
# incomplete_file <- constructPath("Incomplete")
# imputed_file <- constructPath("Imputed")
# jvmGBLimit <- 8

source(paste0(project_folder, '/Code/Functions/lsaImputation.R'))
source(paste0(project_folder, '/Code/Functions/knnImputation.R'))


################################################################
# read in and clean raw methylation. 
################################################################

# set empy vectors to fill with methylation
num_sets <- 23 
methylation_raw <-  vector('list', num_sets)

# read in raw methylation files. 
for(i in (1:num_sets)){
  methylation_raw[[i]] <- read.delim(paste(data_folder, data_name, i, '.txt', sep = ''), header = TRUE)
  sub_dat <- methylation_raw[[i]]
  column_split <- strsplit(as.character(sub_dat$X), '#')
  last_digits <- lapply(column_split, function(x) x[length(x)])
  sub_ids <- unlist(last_digits)
  sub_ids <- gsub('RD-', '', sub_ids)
  methylation_raw[[i]][,1] <- sub_ids
  if(i > 1){
    methylation_raw[[i]] <- methylation_raw[[i]][, -1]
  }
  print(i)
}

# concatenate methylation_raw
methylation <- do.call('cbind', methylation_raw)

# make id column for join later
names(methylation)[1] <- 'id'


# clean probe names 
cleanProbe <- function(data){
  
  col_names <- colnames(data)
  data <- data[,!grepl('ch|_|X', colnames(data))]
  data$id <- gsub('-', '_', data$id)
  rownames(data) <- NULL
  
  
  return(data)
  
}

# apply function and get new methylation. also return index for removed methylation
methylation_new <- cleanProbe(methylation)

###################################################################
# impute on methylation_new
###################################################################

# remove duplicated ids to put in rownames of methylation
methylation_new <- methylation_new[!duplicated(methylation_new$id),]
rownames(methylation_new) <- methylation_new$id
methylation_new$id <- NULL

# impute on methyl with knn 
methylation_new <- knnImputation(methylation_new, sample_rows = TRUE)

# join rownames and methyl_impute and then erase rownames
methylation_new <- cbind(id = rownames(methylation_new), methylation_new)
rownames(methylation_new) <- NULL

# make methylation a data frame 
methylation_new <- as.data.frame(methylation_new)

# write.csv(methylation_new, paste0(data_folder, '/methyl_knn.csv'))

#########################################################
# read in clinical data and merge with methylation_new
#########################################################

# read in clinical
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = TRUE)

# standardize ids
methylation_new <- methylation_new[!grepl('A|B|_', methylation_new$id),]
clin$id <- clin$blood_dna_malkin_lab_
clin <- clin[!grepl('A|B|_', clin$blood_dna_malkin_lab_),]

# inner join methylation and clinical
full_data <- inner_join(methylation_new, clin, by = 'id')

# get mut/wt vector
full_data$p53_germline <- factor(full_data$p53_germline, levels = c('WT', 'Mut'))

# remove unnecesaary objects 
rm(clin, methylation, methylation_raw, sub_dat, column_split, last_digits, num_sets)

# save full_data to be used in different bumphunter scripts 
save.image(paste0(data_folder, '/full_data_bh.RData'))

