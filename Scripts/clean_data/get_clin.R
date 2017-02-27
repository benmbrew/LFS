##### This Script will read in clinical data and clean it.
# This is the first step in the pipeline

##########
# initialize folders
##########
library(gsubfn)
library(xlsx)
library(gsheet)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
idat_data <- paste0(methyl_data, '/raw_files')
clin_data <- paste0(data_folder, '/clin_data')

##########
# read in clinical data, sheets 1 and 2
##########
clin1 <- read.csv(paste0(clin_data, '/clin1.csv'), na.strings=c("","NA"), 
                  stringsAsFactors = FALSE) # read the first sheet
clin2 <- read.csv(paste0(clin_data, '/clin2.csv'), na.strings=c("","NA"),
                  stringsAsFactors = FALSE) # read the second sheet


##########
# remove unneeded fields and combine
##########
clin1$X <- NULL
clin2$X <- NULL
clin1$X.1 <- NULL
clin2$Pin53 <- NULL
# remove the unverified WT in clin2 (after row 136)
clin2 <- clin2[1:136,]

# combine the two data sets
clin <- rbind(clin1, clin2)

# remove all columns that are completely NA 
clin <- clin[, colSums(is.na(clin)) < nrow(clin)]

# remove all rows that are completely NA
clin <- clin[rowSums(is.na(clin)) < ncol(clin),]

#########
# clean NAs and N/As
#########
clin <- as.data.frame(apply(clin, 2, function(x) gsub('N/A', NA, x)))

##########
# function to make column names lower case and remove any '.' and replace with '_'
##########
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

##########
# function clean white all leading and trailing white spaces 
##########
convertColumn <- function(data) {
  
  for( i in names(data)) {
    
    data[,i] <- gsub("^\\s+|\\s+$", "", data[, i])
    
  }
  return(data)
}


clin <- convertColumn(clin)

#########
# function to split rows with multiple sample ids
#########
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

#########
# Clean age of diagnosis column
#########
clin$age_diagnosis <- as.character(clin$age_diagnosis)
clin$age_sample_collection <- as.character(clin$age_sample_collection)

##########
# Clean the age variable so that is expressed in years
##########
cleanAge <-  function(data, column_name) {
  
  age_vector <- data[, column_name]
  
  for (i in 1:length(age_vector)) {
    
    temp_age <- age_vector[i]
    
    if(!grepl('1|2|3|4|5|6|7|8|9|0', temp_age)) {
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

##########
# Convert age of diagnosis and sample collection to months 
##########
clin$age_diagnosis <- clin$age_diagnosis*12
clin$age_sample_collection <- clin$age_sample_collection*12

##########
# clean date column and create variable for patient age 
##########
clin$dob <- ifelse(grepl('unknown', as.character(clin$dob)), NA, 
                   ifelse(grepl('notknown', as.character(clin$dob)), NA, as.character(clin$dob)))

##########
# extract only the last four characters to get year of birth
##########
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

clin$yob <- as.numeric(substrRight(as.character(clin$dob), 4))
clin$age <- 2016 - clin$yob


##########
# clean p53
##########
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

##########
# function to clean cancer_diagnoses
##########
clin$cancer_diagnosis_diagnoses <- as.character(clin$cancer_diagnosis_diagnoses)

cleanCancer <- 
  
  function (data, column_name) {
    
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

##########
# clean p53
##########
clin$gender <- as.character(clin$gender)
clin$gender <- ifelse(clin$gender == 1, 'M', 
                      ifelse(clin$gender == 0, 'F', 
                             ifelse(clin$gender == 'unknown', NA, clin$gender)))


##########
# remove leading and trailing white space 
##########
clin$blood_dna_malkin_lab_ <- gsub("^\\s+|\\s+$", "", clin$blood_dna_malkin_lab_)

# replace LFS with family
clin$family_name <- gsub('LFS', 'Family', clin$family_name)


##########
# functon to clean family column 
##########
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

##########################################################################################
##########
# lauren's changes
##########
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

###
#     GDNA VARS
###
clin[,'gdna.exon.intron'] <- NA
clin[,'gdna.exon.intron'][is.na(clin$gdna) == FALSE] <- 'no clear intron or exon'
clin[,'gdna.exon.intron'] <- sub(":.*$", "", clin$gdna )
clin[,'gdna.exon.intron'][grep(pattern = "Deletion",x = clin[,'gdna.exon.intron'])] <- 'no clear intron or exon'
clin[,'gdna.exon.intron'][grep(pattern = "insertion",x = clin[,'gdna.exon.intron'])] <- 'no clear intron or exon'
clin[,'gdna.exon.intron'][grep(pattern = "c.",x = clin[,'gdna.exon.intron'])] <- 'no clear intron or exon'
clin[,'gdna.exon.intron'] <- trim(clin[,'gdna.exon.intron'])

head(clin$gdna.exon.intron)
unique(clin$gdna.exon.intron)
clin$gdna.exon.intron.fac <- factor(clin$gdna.exon.intron,levels=c("no clear intron or exon","Exon 2","Intron 3","Exon 4","Intron 4",
                                                                   "Exon 5","Intron 5","Exon 6","Exon 7","Exon 8","Intron 8","Exon 10"))
unique(clin$gdna.exon.intron.fac)

clin[,'gdna.base.change'] <- NA
clin[,'gdna.base.change'][is.na(clin$gdna) == FALSE] <- 'no bp change'
clin[,'gdna.base.change'][grep(pattern = '>',x = clin$gdna)] <- unlist(lapply(X = strsplit(x = clin$gdna,split = '>'),
                                                                              FUN = function(x){if(length(x) == 2){
                                                                                str1 <- substr(x[1],nchar(x[1]),nchar(x[1]))
                                                                                str2 <- substr(x[2],1,1)
                                                                                out.val <- paste0(str1,'>',str2)
                                                                                return(out.val)
                                                                              }
                                                                              }))
clin[,'gdna.base.change'] <- trim(clin[,'gdna.base.change'])

clin$gdna.base.change.fac <- factor(clin$gdna.base.change,levels=c("no bp change","T>A","A>G","A>C","G>A","C>A","C>G",
                                                                   "G>C","C>T","T>C","G>T","T>G"))
unique(clin$gdna.base.change.fac)
table(clin$gdna.base.change.fac)

clin[,'gdna.codon'] <- NA
clin[,'gdna.codon'][is.na(clin$gdna) == FALSE] <- 'no clear location'
clin[,'gdna.codon'][which(substr(clin$gdna,1,4) == 'Exon')] <- substr(x = sub(">.*","",gsub("^[^.]*.", "", clin$gdna[which(substr(clin$gdna,1,4) == 'Exon')])),start = 1,stop = nchar(sub(">.*","",gsub("^[^.]*.", "", clin$gdna[which(substr(clin$gdna,1,4) == 'Exon')])))-1)
clin[,'gdna.codon'][which(substr(clin$gdna,1,6) == 'Intron')] <- substr(x = sub(">.*","",gsub("^[^.]*.", "", clin$gdna[which(substr(clin$gdna,1,6) == 'Intron')])),start = 1,stop = nchar(sub(">.*","",gsub("^[^.]*.", "", clin$gdna[which(substr(clin$gdna,1,6) == 'Intron')])))-1)
clin[,'gdna.codon'][grep(pattern = "d",clin[,'gdna.codon'])] <- 'no clear location'
clin[,'gdna.codon'][grep(pattern = "_",clin[,'gdna.codon'])] <- 'no clear location'
clin[,'gdna.codon'][grep(pattern = "-",clin[,'gdna.codon'])] <- 'no clear location'
clin[,'gdna.codon'][grep(pattern = "/+",clin[,'gdna.codon'])] <- 'no clear location'
clin[,'gdna.codon'][clin[,'gdna.codon'] == ""] <- 'no clear location'
clin[,'gdna.codon'] <- trim(clin[,'gdna.codon'])

unique(clin$gdna.codon)

###
#     PROTEIN VARS
###

clin[,'protein.codon.change'] <- NA
clin[,'protein.codon.change'][is.na(clin$protein) == FALSE] <- trim(unlist(lapply(strsplit(x = sub(pattern = " / splice$",replacement = "",
                                                                                                   x = sub(pattern = 'p.',replacement = '',
                                                                                                           x = clin$protein[is.na(clin$protein) == FALSE])),split = '[0-9]+'),
                                                                                  function(x){if(length(x) == 1){
                                                                                    return('no_codon_change')
                                                                                  }
                                                                                    if(length(x) == 2){
                                                                                      return(paste0(x[1],'>',x[2]))
                                                                                    }                     
                                                                                    else{
                                                                                      return(NA)
                                                                                    }
                                                                                    
                                                                                  })))


unique(clin$protein.codon.change)
clin$protein.codon.change <- unlist(lapply(strsplit(clin$protein.codon.change,split = " "),FUN = function(x){x[1]}))
unique(unlist(lapply(strsplit(clin$protein.codon.change,split = " "),FUN = function(x){x[1]})))
clin$protein.codon.change.fac <- factor(clin$protein.codon.change,levels=c("no_codon_change","Cys>Ter","Phe>Tyr","Tyr>Cys","His>Pro","Thr>Thr","Arg>His",
                                                                           "Gly>Ser","Arg>Ser","Pro>fs","Glu>fs","Gln>Pro","Ser>Tyr","Arg>Gln","Trp>Ter",
                                                                           "Arg>Cys","Thr>fs","Trp>X","Pro>Leu","Pro>Ser","Gly>Arg","Gln>fs","Cys>Arg",        
                                                                           "Arg>Leu","Ile>Thr","Arg>Pro","Gly>Cys","Glu>Lys","Arg>X","Tyr>Stop","Tyr>X",          
                                                                           "Val>Met","Thr>Hisfs*","Arg>Trp","Leu>Pro","Ile>Ser","Val>Ile","His>Tyr","Tyr>*",          
                                                                           "Ser>Gly","Leu>Gln","Pro>Ala","Thr>Pro","Cys>Tyr","Phe>Ser","Ser>fs"))


clin[,'protein.codon.num'] <- NA
clin[,'protein.codon.num'][is.na(clin$protein) == FALSE] <- as.numeric(gsub("\\D", "", clin$protein[is.na(clin$protein) == FALSE]))
unique(clin$protein.codon.num)

clin[,'splice.delins.snv'] <- NA
clin[,'splice.delins.snv'][grep(pattern = '>',clin$gdna.base.change)] <- 'SNV'
clin[,'splice.delins.snv'][grep(pattern = 'deletion',clin$protein,ignore.case = T)] <- 'Deletion'
clin[,'splice.delins.snv'][grep(pattern = 'splice',clin$protein,ignore.case = T)] <- 'Splice'
clin[,'splice.delins.snv'][grep(pattern = 'dup',clin$gdna,ignore.case = T)] <- 'Duplication'
table(clin$splice.delins.snv)

clin$codon72.npro <- NA
clin$codon72.npro[clin$codon72 == 'arg/arg'] <- 0
clin$codon72.npro[clin$codon72 == 'arg/pro'] <- 1
clin$codon72.npro[clin$codon72 == 'pro/pro'] <- 2
table(clin$codon72)
table(clin$codon72.npro)

clin$mdm2.nG <- NA
clin$mdm2.nG[clin$mdm2 == 'T/T'] <- 0
clin$mdm2.nG[clin$mdm2 == 'T/G'] <- 1
clin$mdm2.nG[clin$mdm2 == 'G/G'] <- 2
table(clin$mdm2)
table(clin$mdm2.nG)

######
# # find errors
# temp_age <- clin$age_diagnosis
# temp_sample <- clin$age_sample_collection
# temp_cancer <- clin$cancer_diagnosis_diagnoses
# temp_p53  <- clin$p53_germline
# temp_id <- clin$blood_dna_malkin_lab_
# 
# temp_id[grepl('y|m|w', temp_id)]
# temp_id[grepl("[[:digit:]]", temp_id)]
# write clin to data_folder so it can be loaded to database
write.csv(clin, paste(clin_data,'clinical_two.csv', sep ='/'), row.names = FALSE)

