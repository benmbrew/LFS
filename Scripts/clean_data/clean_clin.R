### This Script will read in clinical data and clean it.

# Initialize folders
home_folder <- '/home/benbrew/Documents'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
clin_data <- paste0(data_folder, '/clin_data')

# Read in clinical data 
clin <- read.csv(paste(data_folder, '/clin_all.csv', sep = ''), header = FALSE)

# Identify rows with two ids and create a new row for each one. One row has 
# duplicates 1087 and 1087/1094

splitRows <- function(data, duplicate_table){
  
  for (i in 1:nrow(data)){
    sub_data <- data[i,]
    if(grepl('/', sub_data$V1)){
      split_id <- strsplit(as.character(sub_data$V1), '/')
      split_id <- cbind(unlist(split_id))
      duplicate <- cbind(split_id, sub_data)
      duplicate_table <- rbind(duplicate_table, duplicate)
    }
    
  }
  data <- data[!grepl('/', data$V1),]
  duplicate_table$V1 <- NULL
  colnames(duplicate_table)[1] <- 'V1'
  data <- rbind(data, duplicate_table)
  return(data)
  
}
empty_table <- data.frame(matrix(ncol = ncol(clin), nrow = 0))
clin <- splitRows(clin, empty_table)

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
                                  ifelse(clin$codon == 'arg pro', 'arg/pro',
                                         ifelse(clin$codon == 'arg/pro?', 'arg/pro',
                                                ifelse(clin$codon_72 == 'arg homo', 'arg/arg',
                                                       ifelse(clin$codon_72 == '', 'missing', clin$codon_72))))))

clin$gender<- as.factor(ifelse(clin$gender == 1, 'female', 'male'))
clin$gender <- as.factor(clin$gender)

clin$mdm2 <- as.character(clin$mdm2)
clin$mdm2 <- as.factor(ifelse(clin$mdm2 == '', 'missing', clin$mdm2)) 

# write clin to data_folder so it can be loaded to database
write.csv(clin, paste(clin_data,'clinical.csv', sep ='/'), row.names = FALSE)
