## Script for reading, cleaning, and combining methylation data and clinical data.

# Initialize folders
home_folder <- '/home/benbrew/Documents'
project_folder <- paste(home_folder, 'LFS', sep = '/')
data_folder <- paste(project_folder, 'Data/', sep = '/')
clin_data <- paste(data_folder, 'clin_data', sep ='/')
methyl_data <- paste(data_folder, 'methyl_data/', sep = '/')
data_name <- 'Chr'

# Read in clinical data 
clin <- read.csv(paste(data_folder, 'clin_all.csv', sep = ''), header = FALSE)

# Identify rows with two ids and create a new row for each one. One row has 
# duplicates 1087 and 1087/1094

splitClinRows <- function(data, duplicates){
  
  for (i in 1:nrow(data)){
  sub_data <- data[i,]
  if(grepl('/', sub_data$V1)){
    split_id <- strsplit(as.character(sub_data$V1), '/')
    split_id <- cbind(unlist(split_id))
    duplicate <- cbind(split_id, sub_data)
    duplicates <- rbind(duplicates, duplicate)
  }
  
}
data <- data[!grepl('/', data$V1),]
duplicates$V1 <- NULL
colnames(duplicates)[1] <- 'V1'
data <- rbind(data, duplicates)
rownames(data) <- NULL
return(data)

}
duplicates <- data.frame(matrix(ncol = 12, nrow = 0))
clin <- splitClinRows(clin, duplicates)

# There is a duplicate id 1087. remove the first one because it has less data in row 62
clin <- clin[!duplicated(clin$V1, fromLast = TRUE),]
rownames(clin) <- NULL

# name columns according to the tables in google drive  
colnames(clin) <- c("id", "tp53", "cancer", "cancer_indicator", "age_of_onset",
                           "gdna", "protein", "codon_72", "pin_3", "mdm2","date", "gender")

# Recode variables.  For now, don't treat 'blank' fields as NA, but preserve them as random effects. 
clin$pin_3 <- as.character(clin$pin_3)
clin$pin_3 <- as.factor(ifelse(clin$pin_3 == 'A1/A1?', 'A1/A1',
                                    ifelse(clin$pin_3 == '', 'missing', clin$pin_3)))

clin$protein <- as.character(clin$protein)
clin$protein <- as.factor(ifelse(clin$protein == 'N/A', 'missing', 
                                 ifelse(clin$protein == '', 'missing',
                                        ifelse(clin$protein == 'n/a', 'missing', clin$protein))))

clin$codon_72 <- as.character(clin$codon_72)
clin$codon_72 <- as.factor(ifelse(clin$codon_72 == 'arg?', 'missing',
                                       ifelse(clin$codon_72 == 'arg homo', 'arg/arg',
                                              ifelse(clin$codon_72 == '', 'missing', clin$codon_72))))

clin$gender<- as.factor(ifelse(clin$gender == 1, 'female', 'male'))
clin$gender <- as.factor(clin$gender)

clin$mdm2 <- as.character(clin$mdm2)
clin$mdm2 <- as.factor(ifelse(clin$mdm2 == '', 'missing', clin$mdm2)) 

# write clin to data_folder so it can be loaded to database
write.csv(clin, paste(clin_data,'clinical.csv', sep ='/'), row.names = FALSE)

#####################################################################

# Read in methylation data 
methylation24 <- read.table(paste(data_folder,'methyl.txt', sep = ''), header = TRUE)
methylation17 <- read.csv(paste(data_folder, 'methyl_17.csv', sep = ''), header = TRUE)



# Change into format so that it can be combined with clin 
cleanMethyl <- function(data){
  
  data <- as.data.frame(t(data), stringsAsFactors = FALSE)
  names(data) <- data[1,]
  data <- data[-1,]
  data <- cbind(id = rownames(data), data)
  rownames(data) <- NULL
  split <- strsplit(as.character(data$id), 'X')
  last_split <- lapply(split, function(x) x[length(x)])
  data$id <- unlist(last_split)
  data <- t(data)
  return(data)
  
}

methylation24 <- cleanMethyl(methylation24)
methylation17 <- cleanMethyl(methylation17)

# Read in raw mehtylation data
num_sets <- 23 
methylation_raw <-  vector('list', num_sets)

# read in raw methylation files. 
for(i in (1:num_sets)){
  methylation_raw[[i]] <- read.delim(paste(data_folder, data_name, i, '.txt', sep = ''))
  
}

# Aggregate list 
# Attach methylation data to clin_full 
combineMethyl <- function(methylation_data, num_sets){
  
  for (i in 1:num_sets){
    sub_dat <- methylation_data[[i]]
    column_split <- strsplit(as.character(sub_dat$X), '#')
    last_digits <- lapply(column_split, function(x) x[length(x)])
    sub_ids <- unlist(last_digits)
    sub_ids <- gsub('RD-', '', sub_ids)
    methylation_data[[i]][,1] <- sub_ids
    if(i > 1){
      methylation_data[[i]] <- methylation_data[[i]][, -1]
    }
  }
  methylation <- do.call('cbind', methylation_data)
  methylation <- as.data.frame(t(methylation), stringsAsFactors = FALSE)
  data <- cbind(x = rownames(methylation), methylation)
  colnames(data) <- data[1,]
  data <- data[-1,]
  names(data)[1] <- 'probe'
  rownames(data) <- NULL  
  return(data)
}


methylation <- combineMethyl(methylation_raw, num_sets)

cleanProbe <- function(data){
  for(i in 2:ncol(data)){
  data[,i] <- as.numeric(data[,i])
  }
data <- data[!grepl('ch', data$probe),]
data <- data[!grepl('_', data$probe),]
colnames(data) <- gsub('-', '_', colnames(data))
rownames(data) <- NULL
data <- data[data$probe != 'X.1',]

return(data)
}

methylation <- cleanProbe(methylation)

# Clean raw methylation- getting rid of Ch information and blanks

# write methylation_raw to methyl_data folder as csv
write.csv(methylation, paste(methyl_data, 'methylation.csv', sep = ''), row.names= FALSE)

# save data image
rm(duplicates)
setwd('/home/benbrew/Documents/LFS/Data')
save.image('cleaned.RData')



