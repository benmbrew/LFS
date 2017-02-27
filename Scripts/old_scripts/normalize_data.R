####### Script will normalize cases, controls, both for batch corrected and not
# # This is 6th step and final step of data cleaning and preparation

##########
# initialize libraries
##########
library(minfi)
library(preprocessCore)

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
# read batch corrected data
##########

# read cases 
raw_cases_batch <- readRDS(paste0(model_data, '/raw_cases_batch.rda'))

# read contorls
raw_controls_batch <- readRDS(paste0(model_data, '/raw_controls_batch.rda'))

##########
# read non batch corrected
##########

# read cases 
raw_cases <- readRDS(paste0(model_data, '/raw_cases.rda'))

# read contorls
raw_controls<- readRDS(paste0(model_data, '/raw_controls.rda'))


##########
# noramlize data 
##########

normalizeData <- function(data)
{
  features <- colnames(data)[8:ncol(data)]
  clin_data <- data[,1:7]
  methyl_data <- t(data[, 8:ncol(data)])
  methyl_data_norm  <- normalize.quantiles(methyl_data)
  methyl_data_norm <- as.data.frame(t(methyl_data_norm), stringsAsFactors = T)
  colnames(methyl_data_norm) <- features
  methyl_full <- cbind(clin_data, methyl_data_norm)
  return(methyl_full)
  
}

# get normalized cases and controls for batch
quan_cases_batch <- normalizeData(raw_cases_batch)
quan_controls_batch <- normalizeData(raw_controls_batch)


# get normalized cases and controls for non batch
quan_cases <- normalizeData(raw_cases)
quan_controls <- normalizeData(raw_controls)

# remove unneeded objects
rm(normalizeData)

# save R data image
save.image(paste0(model_data, '/model_data.RData'))




