## Script for reading, cleaning, and combining methylation data.
# this is the second step in the pipeline
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/classification_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

data_name <- '/Chr'

# #################################################################################
# # Take Nardine's methyalation data and clean the column names
# #################################################################################
# 
# # Read in methylation data 
# methylation24 <- read.table(paste(data_folder,'/methyl.txt', sep = ''), header = TRUE)
# methylation17 <- read.csv(paste(data_folder, '/methyl_17.csv', sep = ''), header = TRUE)
# 
# # This function takes the mehtylation data and removes the 'X' from the colnames 
# 
# cleanColNames <- function(data){
#   for(i in 1:ncol(data)){
#     sub_dat <- colnames(data)[i]
#     if(grepl('X', sub_dat)){
#       split <- strsplit(sub_dat, 'X')
#       last_split <- lapply(split, function(x) x[length(x)])
#       colnames(data)[i] <- unlist(last_split)
#     }else{
#       colnames(data)[i] <- sub_dat
#     }
#   }
#   return(data)
# }
# 
# methylation24 <- cleanColNames(methylation24)
# methylation17 <- cleanColNames(methylation17)

##############################################################################
# Read in raw methylation data from Nardine and 
# Combine the list and put in same format as previous 
# methylation with probes in first column and IDs in the following columns
##############################################################################

# Read in raw mehtylation data
num_sets <- 23 
methylation_raw <-  vector('list', num_sets)

# read in raw methylation files. 
for(i in (1:num_sets)){
  methylation_raw[[i]] <- read.delim(paste(data_folder, data_name, i, '.txt', sep = ''), header = TRUE)
}

# iterate through each data set and retrieve the IDs then combine them into one data table in same format as 
# Previous methylation data

for (i in 1:num_sets){
  sub_dat <- methylation_raw[[i]]
  column_split <- strsplit(as.character(sub_dat$X), '#')
  last_digits <- lapply(column_split, function(x) x[length(x)])
  sub_ids <- unlist(last_digits)
  sub_ids <- gsub('RD-', '', sub_ids)
  methylation_raw[[i]][,1] <- sub_ids
  if(i > 1){
    methylation_raw[[i]] <- methylation_raw[[i]][, -1]
  }
}
methylation <- do.call('cbind', methylation_raw)
methylation <- as.data.frame(t(methylation), stringsAsFactors = FALSE)
methylation <- cbind(x = rownames(methylation), methylation)
colnames(methylation) <- methylation[1,]
methylation <- methylation[-1,]
names(methylation)[1] <- 'Probe'
rownames(methylation) <- NULL  


cleanProbe <- function(data){
  for(i in 2:ncol(data)){
    data[,i] <- as.numeric(data[,i])
  }
  data <- data[!grepl('ch', data$Probe),]
  data <- data[!grepl('_', data$Probe),]
  colnames(data) <- gsub('-', '_', colnames(data))
  data <- data[data$Probe != 'X.1',]
  rownames(data) <- NULL
  

return(data)

}

methylation <- cleanProbe(methylation)
# write methylation_raw to methyl_data folder as csv
write.csv(methylation, paste(methyl_data, '/methylation.csv', sep = ''), row.names= FALSE)

# ###########################################################
# # Read in methylation data from tanya and clean it
# ###########################################################
# 
# methyl_tumor <- read.delim(paste0(methyl_data, '/methylation_tumor.txt'), check.names = FALSE)
# 
# # subset columns that just contain beta values 
# methyl_beta <- methyl_tumor[, grepl('TargetID|Beta',  names(methyl_tumor))]
# col_names <- names(methyl_beta)
# column_split <- strsplit(col_names, '.', fixed = TRUE)
# first_digits <- lapply(column_split, function(x) x[1])
# col_names <- unlist(first_digits)
# 
# # add new column names
# colnames(methyl_beta) <- col_names
# colnames(methyl_beta)[1] <- 'Probe'
# 
# write.csv(methyl_beta, paste(methyl_data, '/methylation_tumor.csv', sep = ''), row.names= FALSE)
# 
# 
# 
