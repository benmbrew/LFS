###########################################################33
# This script will run the initial pass of the Lmm model, trained on the non methylation and 
# tested on methylation. 

# initialize folders
library(gsubfn)
library(dplyr)


### This Script will explore clinical data 
# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/classification_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

# read in clinical data 
# clin <- read.csv(paste0(clin_data, '/clinical_mut.csv'), stringsAsFactors = F)
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)


# read in methylation column names
methylation_names <- colnames(read.csv(paste0(methyl_data, '/methylation.csv'), 
                              stringsAsFactors = F, check.names = FALSE, nrows = 1))

# remove 'A' and '_' in methylation names
methylation_names <- gsub('_', '', methylation_names)
methylation_names <- gsub('A', '', methylation_names)

##############################################################
# Get clin in format to run a model

# remove leading and trailing white space 
clin$blood_dna_malkin_lab_ <- gsub("^\\s+|\\s+$", "", clin$blood_dna_malkin_lab_)


# differences in family in regards to clinical variables 

# clean family column 
for (i in 1:nrow(clin)) {
  
  sub_clin <- clin$family_name[i]
  
  if(is.na(sub_clin)) {
    temp <- NA
    
  } else {
    
    if(nchar(sub_clin) == 11) {
      temp <- substr(sub_clin, 1, 8)
    }
    
    if(nchar(sub_clin) == 12) {
      temp <- substr(sub_clin, 1, 9)
    }
    
    if(nchar(sub_clin) == 13) {
      temp <- substr(sub_clin, 1, 10)
    }
    
  }
  clin$family_name[i] <- temp
  
}

clin$family_name <- tolower(gsub(" ", "_", clin$family_name))

############################################################################
# add a column to clin to indicate if methylation data is available
# add '1' and make data frame 
methylation_names <- as.data.frame(cbind(methylation_names, rep.int(1, length(methylation_names))))
methylation_names$methylation_names <- as.character(methylation_names$methylation_names)
methylation_names$V2 <- as.character(methylation_names$V2)
names(methylation_names) <- c("blood_dna_malkin_lab_", "methyl_indicator")
methylation_names <- methylation_names[!duplicated(methylation_names),]

# join methylation names to clin
clin <- left_join(clin, methylation_names)

# recode true and false 
clin$methyl_indicator <- ifelse(is.na(clin$methyl_indicator), 'No', 'Yes')

# write to data as clin_only 
write.csv(clin, paste(clin_data,'clin_only.csv', sep ='/'), row.names = FALSE)
#############################################################################
# subset to those variables 
subset <- c("family_name", "relationship", "age_diagnosis", "p53_germline","gdna", 
            "protein", "codon72", "mdm2", "gender", 'methyl_indicator')

clin <- clin[, subset]

# Try the model with all different selection of features based on number of missinginess. 
temp <- clin[complete.cases(clin),]

# convert characters to factors 
for ( i in 1:ncol(clin)){
  
  if (class(clin[,i]) %in% c('character')) {
     clin[,i] <- as.factor(clin[,i])
     
  } 
}



