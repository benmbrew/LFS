# initialize folders
library(gsubfn)
library(xlsx)
library(gsheet)

### This Script will read in clinical data and clean it.
# This is the first step in the pipeline
# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/classification_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

#######################################################################
# clean clinical data sent from Ana on 1/29/2016 
# clin1 <- as.data.frame(gsheet2tbl('https://docs.google.com/spreadsheets/d/1zUOYEXFh9RAQFpNPFjpPNBSLbFgnqQItucuNzuznivk/edit#gid=0'))
# write.csv(clin1, '/home/benbrew/Desktop/clin1.csv')
# clin2 <- gsheet2tbl('https://docs.google.com/spreadsheets/d/1zUOYEXFh9RAQFpNPFjpPNBSLbFgnqQItucuNzuznivk/edit#gid=621772204')
# write.csv(clin2, '/home/benbrew/Desktop/clin2.csv')

clin1 <- read.csv(paste0(clin_data, '/clin1.csv'), na.strings=c("","NA"), 
                 stringsAsFactors = FALSE, sep = ',') # read the first sheet
clin2 <- read.csv(paste0(clin_data, '/malkin_june_2.csv'), na.strings=c("","NA"),
                stringsAsFactors = FALSE) # read the second sheet

clin1$X <- NULL

# For the time being drop family from clin1
clin2$X <- NULL
# clin1$Family.Name <- NULL
clin2$Pin53 <- NULL

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

# clean white all leading and trailing white spaces 

convertColumn <- function(data) {
  
  for( i in names(data)) {
  
  data[,i] <- gsub("^\\s+|\\s+$", "", data[, i])
  
  }
  return(data)
}


clin <- convertColumn(clin)

##########################################################################
# clean malkin ids 
clin$blood_dna_malkin_lab_ <- as.character(clin$blood_dna_malkin_lab_)
NAs <- 'NA|nogermlineDNA|nosample|Nosample|samplefailed&discarded'
ABs <- 'A|B'

cleanIds <- function(data, column_name) { 
    
    id_vector <- data[, column_name]
    
    for (i in 1:length(id_vector)) {
      
      temp_id <- id_vector[i]
      
      if (grepl(NAs, temp_id)) {
        temp_id <- NA
      }
      
      if (grepl(ABs, temp_id)) {
        temp_id <- substring(temp_id, 1, 4)
      }
      
      if (grepl('receivedBMtransplant', temp_id)) {
        
        temp.2_id <- strsplit(temp_id, '(receivedBMtransplant)', , fixed = TRUE)
        temp.3_id <- strsplit(temp.2_id[[1]][2], '-')
        temp.4_id <- paste(unlist(temp.2_id[[1]][1]), temp.3_id[[1]][1], sep = '/')
        temp_id <- temp.4_id
        
      }

      if (grepl('receivedBMtransplantx2', temp_id)) {
        
        temp.2_id <- strsplit(temp_id, '(receivedBMtransplantx2)', , fixed = TRUE)
        temp.3_id <- strsplit(temp.2_id[[1]][2], '-')
        temp.4_id <- paste(unlist(temp.2_id[[1]][1]), temp.3_id[[1]][1], sep = '/')
        temp_id <- temp.4_id
        
      }
      
      id_vector[i] <- temp_id
    }
    
    data[, column_name] <- id_vector
    return(data)
    
}

clin <- cleanIds(clin, column_name = 'blood_dna_malkin_lab_')

##############################################
# take any rows with double samples and split them into the number of rows equal to the number of samples.

splitRows <- function(data, duplicate_table){
  
  for (i in 1:nrow(data)) {
    sub_data <- data[i,]
    
    if (grepl('/', sub_data$blood_dna_malkin_lab_) | 
        grepl('/', sub_data$age_sample_collection)) {
      
      split_malkin <- strsplit(as.character(sub_data$blood_dna_malkin_lab_), '/')
      split_age <- strsplit(as.character(sub_data$age_sample_collection), '/')
      split_malkin <- cbind(unlist(split_malkin))
      split_age <- cbind(unlist(split_age))
      duplicate <- cbind(split_malkin, split_age, sub_data)
      duplicate_table <- rbind(duplicate_table, duplicate)
      
    }
  }
  
  data <- data[!grepl('/', data$blood_dna_malkin_lab_),]
  data <- data[!grepl('/', data$age_sample_collection),]
  
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

cleanAge <-  function(data, column_name) {
  
  age_vector <- data[, column_name]
    
  for (i in 1:length(age_vector)) {
      
    temp_age <- age_vector[i]
    
    if(!grepl('y|m|d', temp_age)){
      temp_age <- NA
    
    } else {
        
      
      if (grepl('~', temp_age, fixed = TRUE)) {
        
        temp_age <- strsplit(temp_age, '~', fixed = TRUE)
        
        temp.2_age <- unlist(temp_age)
        
        temp.3_age <- temp.2_age[2]
        
        temp_age <- temp.3_age
        
      }
    
      if (grepl('(', temp_age, fixed = TRUE)) {
        
        temp_age <- strsplit(temp_age, '(', fixed = TRUE)
        
        temp.2_age <- unlist(temp_age)
       
        temp.3_age <- temp.2_age[1]
       
        temp_age <- temp.3_age
        
      }
     
      if (grepl('>', temp_age) || grepl('<', temp_age)) {
        
        temp_age <- gsub('<', '' , temp_age)
        
        temp_age <- gsub('>', '', temp_age)
        
      }
      
      # if (grepl('k|NA|noinfo|no info|Unaffected|unable to retest', temp_age)) {
      #   
      #   temp_age <- NA
      #   
      # }
    
      
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
    
        } else if (grepl('d', temp_age)) {
    
          temp_age <- gsub('2d', '0.005', temp_age)
        } 
    }
    age_vector[i] <- temp_age
    
  }
  
  data[, column_name] <- age_vector
  data[, column_name] <- as.numeric(as.character(data[, column_name]))
  return(data)

} 

clin <- cleanAge(clin, column_name = 'age_diagnosis')
clin<- cleanAge(clin, column_name = 'age_sample_collection')

# Convert age of diagnosis and sample collection to months 
clin$age_diagnosis <- clin$age_diagnosis*12
clin$age_sample_collection <- clin$age_sample_collection*12

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
# clean p53
clin$p53_germline <- as.character(clin$p53_germline)

cleanp53 <- 
  
  function(data, column_name) {
    
    p53_vector <- data[, column_name]
    
    for(i in 1:length(p53_vector)){
      
      temp_p53 <- p53_vector[i]
      
      if (grepl('noresult', temp_p53, fixed = TRUE) || grepl('nottested', temp_p53, fixed = TRUE) && 
            !grepl('WT', temp_p53) || 
            grepl('Nottested', temp_p53, fixed = TRUE)) {
        temp_p53 <- NA
        }
      
      if (grepl('paraffinblocktestingfailed|obligate|Declinedtesting|not tested|Not tested|paraffin|Declined testing', 
                temp_p53)) {
        temp_p53 <- NA
      }
      
      if (grepl('Mut', temp_p53, fixed = TRUE)) {
        temp_p53 <- 'Mut'
      }
      
       if (grepl('WT', temp_p53, fixed = TRUE) && grepl('nottested', temp_p53)) {
        temp_p53 <- 'WT'
      }
      
      if (grepl('WT', temp_p53, fixed = TRUE) && !grepl('nottested', temp_p53)) {
        temp_p53 <- 'WT'
      }
     
      p53_vector[i] <- temp_p53
    }
  data[, column_name] <- p53_vector
  return(data)
}

clin <- cleanp53(clin, column_name = 'p53_germline')

############################################################
# clean cancer_diagnoses
# clean p53
clin$cancer_diagnosis_diagnoses <- as.character(clin$cancer_diagnosis_diagnoses)

cleanCancer <- 
  
  function(data, column_name) {
    
    cancer_vector <- data[, column_name]
    
    for(i in 1:length(cancer_vector)){
      
      temp_cancer <- cancer_vector[i]
      
      if (grepl('P', temp_cancer) && grepl('-', temp_cancer) && grepl('OS', temp_cancer) && !grepl('/', temp_cancer)) {
        temp.2_cancer <- strsplit(temp_cancer, '-')
        temp.3_cancer <- unlist(temp.2_cancer[[1]][2])
        temp.4_cancer <- gsub('?Lowgrade', '', temp.3_cancer, fixed = TRUE)
        temp_cancer <- temp.4_cancer
      }
      
      if (grepl('P', temp_cancer) && grepl('-', temp_cancer) && !grepl('OS', temp_cancer) && !grepl('/', temp_cancer)) {
        temp.2_cancer <- strsplit(temp_cancer, '-')
        temp.3_cancer <- unlist(temp.2_cancer[[1]][2])
        temp.4_cancer <- gsub('?Lowgrade', '', temp.3_cancer, fixed = TRUE)
        temp_cancer <- temp.4_cancer
      }
      
      if (grepl('P', temp_cancer) && grepl('OS', temp_cancer) && grepl('-', temp_cancer) && grepl('left', temp_cancer)) {
        temp.2_cancer <- strsplit(temp_cancer, '-')
        temp.3_cancer <- strsplit(temp.2_cancer[[1]][2], ',')
        temp.4_cancer <- temp.3_cancer[[1]][1]
        temp_cancer <- temp.4_cancer
      }
      
      if (grepl('P', temp_cancer) && grepl('OS', temp_cancer) && grepl('-', temp_cancer) && grepl('right', temp_cancer)) {
        temp.2_cancer <- strsplit(temp_cancer, '-')
        temp.3_cancer <- strsplit(temp.2_cancer[[1]][2], ',')
        temp.4_cancer <- temp.3_cancer[[1]][1]
        temp_cancer <- temp.4_cancer
      }
      
      if (grepl('P', temp_cancer) && grepl('-', temp_cancer) && grepl('/', temp_cancer)) {
        temp.2_cancer <- strsplit(temp_cancer, '-')
        temp.3_cancer <- unlist(temp.2_cancer[[1]][2])
        temp.4_cancer <- strsplit(temp.3_cancer, '/')
        temp.5_cancer <- paste(unlist(temp.4_cancer[[1]][1]),unlist(temp.4_cancer[[1]][2]), sep = '_')
        temp_cancer <- temp.5_cancer
      }
      
      if (grepl('P', temp_cancer) && grepl('-', temp_cancer) && !grepl('/', temp_cancer) && !grepl('?', temp_cancer)) {
        temp.2_cancer <- strsplit(temp_cancer, '-')
        temp.3_cancer <- unlist(temp.2_cancer[[1]][2])
        temp_cancer <- temp.3_cancer
      }
      
      
      if (grepl(',', temp_cancer, fixed = TRUE)) {
        temp_cancer <- gsub(',', '_', temp_cancer)
      }
      
      if (grepl('?', temp_cancer, fixed = TRUE)) {
        temp_cancer <- gsub('?', '', temp_cancer)
      }
      
      if (grepl('&', temp_cancer, fixed = TRUE)) {
        temp_cancer <- gsub('&', '_', temp_cancer)
      }
      
      if (grepl('notdocumented', temp_cancer, fixed = TRUE)) {
        temp_cancer <- NA
      }
      
      if (grepl('unknown', temp_cancer, fixed = TRUE)) {
        temp_cancer <- NA
      }
      
      if (grepl('Unaffected', temp_cancer)) {
        temp_cancer <- 'Unaffected'
      }
    
      cancer_vector[i] <- temp_cancer
    }
    
    data[, column_name] <- cancer_vector
    return(data)
}


clin <- cleanCancer(clin, column_name = 'cancer_diagnosis_diagnoses')

# reclean cancer


# clean relationship column

############################################################
# clean p53
clin$gender <- as.character(clin$gender)
clin$gender <- ifelse(clin$gender == 1, 'M', 
                      ifelse(clin$gender == 0, 'F', 
                             ifelse(clin$gender == 'unknown', NA, clin$gender)))


###################################################################
    
# write clin to data_folder so it can be loaded to database
write.csv(clin, paste(clin_data,'clinical_two.csv', sep ='/'), row.names = FALSE)

