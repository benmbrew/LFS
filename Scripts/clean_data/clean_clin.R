### This Script will read in clinical data and clean it.

home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/classification_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


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

clin$cancer_indicator <- ifelse(clin$cancer_indicator == 1, TRUE, FALSE)

clin$age_fac <- ifelse(clin$age_of_onset > 5, 'older_5', 'younger_5')

# write clin to data_folder so it can be loaded to database
write.csv(clin, paste(clin_data,'clinical.csv', sep ='/'), row.names = FALSE)

#######################################################################
# clean clinical data sent from Ana on 1/29/2016 
library(xlsx)
library(gsubfn)
clin1 <- read.xlsx(paste0(clin_data, '/malkin_jan_29th.xlsx'), 1, stringsAsFactors = FALSE) # read the first sheet
clin2 <- read.xlsx(paste0(clin_data, '/malkin_jan_29th.xlsx'), 2, stringsAsFactors = FALSE) # read the second sheet

# For the time being drop family from clin1 
clin1$NA. <- NULL

# combine the two data sets 
clin <- rbind(clin1, clin2)

# remove all columns that are completely NA 
clin <- clin[, colSums(is.na(clin)) < nrow(clin)]

# remove all rows that are completely NA
clin <- clin[rowSums(is.na(clin)) < ncol(clin),]

#########################################
# make column names lower case and remove any '.' and replace with '_'
cleanColNames <- function(data) {
  col_names <- colnames(data)
  col_names <- tolower(col_names)
    if (grepl('.', col_names)) {
      col_names <- gsub("([.])\\1+", "_", col_names)
      col_names <- gsub('.', '_', col_names, fixed = TRUE)
  }
  colnames(data) <- col_names
  return(data)
}

clin <- cleanColNames(clin)


# clean NAs and N/As
clin <- as.data.frame(apply(clin, 2, function(x) gsub('N/A', 'NA', x)))

# clean white spaces 
clin <- as.data.frame(apply(clin, 2, function(x) gsub('\\s+', '', x)))

##############################################
# take any rows with double samples and split them into the number of rows equal to the number of samples.

splitRows <- function(data, duplicate_table){
  
  for (i in 1:nrow(data)) {
    sub_data <- data[i,]
    
    if (grepl('/', sub_data$blood_dna_malkin_lab_)) {
      split_malkin <- strsplit(as.character(sub_data$blood_dna_malkin_lab_), '/')
      split_age <- strsplit(as.character(sub_data$age_sample_collection), '/')
      split_malkin <- cbind(unlist(split_malkin))
      split_age <- cbind(unlist(split_age))
      duplicate <- cbind(split_malkin, split_age, sub_data)
      duplicate_table <- rbind(duplicate_table, duplicate)
    }
  }
  
  data <- data[!grepl('/', data$blood_dna_malkin_lab_),]
  duplicate_table$blood_dna_malkin_lab_ <- NULL
  duplicate_table$age_sample_collection <- NULL
  colnames(duplicate_table)[1:2] <- c('blood_dna_malkin_lab_', 'age_sample_collection')
  duplicate_table <- duplicate_table[, colnames(data)]
  data <- rbind(duplicate_table, data)
  return(data)
  
}

empty_table <- data.frame(matrix(ncol = ncol(clin), nrow = 0))
clin <- splitRows(clin, empty_table)

# Clean age of diagnosis column
clin$age_diagnosis <- as.character(clin$age_diagnosis)
clin$age_sample_collection <- as.character(clin$age_sample_collection)

##########################################################
# Clean the age variable so that is expressed in years 

cleanAge <- 
  
  function(data, column_name) {
  
  age_vector <- data[, column_name]
    
  for(i in 1:length(age_vector)){
      
    temp_age <- age_vector[i]
      
      if (grepl('~', temp_age, fixed = TRUE)){
        temp_age <- strsplit(temp_age, '~', fixed = TRUE)
        temp.2_age <- unlist(temp_age)
        temp.3_age <- temp.2_age[2]
        temp_age <- temp.3_age
      }
      if (grepl('(', temp_age, fixed = TRUE) || grepl(')', temp_age, fixed = TRUE)){
        temp_age <- strsplit(temp_age, '(', fixed = TRUE)
        temp.2_age <- unlist(temp_age)
        temp.3_age <- temp.2_age[1]
        temp_age <- temp.3_age
      }
     
      if (grepl('>', temp_age) || grepl('<', temp_age)) {
        temp_age <- gsub('<', '' , temp_age)
        temp_age <- gsub('>', '', temp_age)
      }
      
      if (grepl('k', temp_age)) {
        temp_age <- NA
      }
      
      if (grepl('NA', temp_age)) {
        temp_age <- NA
      }
      
      if (grepl('y', temp_age) && grepl('m', temp_age)) {
        temp_age <- gsubfn('([y,m])', list('y' = '.', 'm' = ''), as.character(temp_age))
        temp.2_age <- strsplit(temp_age, '.', fixed = TRUE)
        temp.3_age <- do.call('rbind', temp.2_age)
        temp.4_age <- gsub('.*\\.', paste0(temp.3_age[1], '.') , as.numeric(temp.3_age[2])/12)
        temp_age <- temp.4_age
   
        } else if (grepl('y', temp_age) && !grepl('m', temp_age)) {
        temp_age  <- gsubfn('y', '.0', as.character(temp_age))
    
        } else if (!grepl('y', temp_age) && grepl('m', temp_age) && !grepl('w', temp_age)) {
        temp_age <- gsubfn('m', '', as.character(temp_age))
        temp_age <- as.numeric(temp_age)/12
    
        } else if (grepl('w', temp_age) && grepl('m', temp_age)) {
        temp_age <- gsubfn('([m,w])', list('m' = '.', 'w' = ''), as.character(temp_age))
        temp_age <- ceiling(as.numeric(temp_age))
   
        } else if (grepl('w', temp_age) && !grepl('m', temp_age)) {
        temp_age <- gsub('w', '', temp_age)
    } 
     
    age_vector[i] <- temp_age
  
  }
  
  data[, column_name] <- age_vector
  return(data)

} 

clin <- cleanAge(clin, column_name = 'age_diagnosis')
clin <- cleanAge(clin, column_name = 'age_sample_collection')


############################################################
# clean date column and create variable for patient age 
clin$dob <- ifelse(grepl('unknown', as.character(clin$dob)), NA, 
                  ifelse(grepl('notknown', as.character(clin$dob)), NA, as.character(clin$dob)))
  
# extract only the last four characters to get year of birth

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

clin$yob <- as.numeric(substrRight(as.character(clin$dob), 4))
clin$age <- 2016 - clin$yob

############################################################

