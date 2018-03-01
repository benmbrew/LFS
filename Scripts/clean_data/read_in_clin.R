
##### This Script will read in clinical data and clean it.
# This is the first step in the pipeline

##########
# initialize folders
##########
library(dplyr)
library(gsubfn)
library(gsheet)
library(readr)
library(Hmisc)
library(memisc)
library(readxl)
library(rowr)
library(stringr)


if('new_clin.RData' %in% dir()){
  
  load('new_clin.RData')
} else {
  # read in data - DWT - dead weight tonnage (tons) - how much the ship can carry
  clin_1  <- read_excel('../../Data/clin_data/Malkin_clinical.xlsx', sheet = 1)
  clin_2  <- read_excel('../../Data/clin_data/Malkin_clinical.xlsx', sheet = 2)
  
  # add MDM2 to LFS family
  clin_1$MDM2 <- NA
  
  # combine the two data sets
  clin <- rbind(clin_1, clin_2)
  
  # remove all columns that are completely NA 
  clin <- clin[, colSums(is.na(clin)) < nrow(clin)]
  
  # remove all rows that are completely NA
  clin <- clin[rowSums(is.na(clin)) < ncol(clin),]
  
  #########
  # clean NAs and N/As
  #########
  clin <- as.data.frame(apply(clin, 2, function(x) gsub('N/A', NA, x)))
  
  # make column names lower case and remove any '.' and replace with '_'
  names(clin) <- tolower(names(clin))
  names(clin) <- gsub('@', '', names(clin), fixed = TRUE)
  names(clin) <- gsub('#', '', names(clin), fixed = TRUE)
  names(clin) <- trimws(names(clin), 'both')
  names(clin) <- gsub('/', '_', names(clin), fixed = TRUE)
  names(clin) <- gsub(" ", '_', names(clin), fixed = TRUE)
  names(clin) <- gsub("__", '_', names(clin), fixed = TRUE)
  
  
  
 # remove leading and trailing white spaces to all columns
  clin <- as.data.frame(apply(clin, 2, function(x){
    trimws(x, 'both')
  }), stringsAsFactors = FALSE
  )
  
  ##########
  # check if properly read in
  ##########
  
  # replace no_sample in malkin id wiht MNA
  clin$blood_dna_malkin_lab <- ifelse(grepl('sample', clin$blood_dna_malkin_lab), 
                                      NA, clin$blood_dna_malkin_lab)
  
  # remove .0 at the end of ids 
  clin$tm_donor <- gsub('.0', '', clin$tm_donor, fixed = TRUE)
  clin$blood_dna_malkin_lab <- gsub('.0', '', clin$blood_dna_malkin_lab, fixed = TRUE)
  
  # remove if NA in tm_donor or 
  
  # clean unknown's from age_sample collection
  clin$age_sample_collection <- gsub('unknown|NA', NA, clin$age_sample_collection, ignore.case = TRUE)
  clin$age_diagnosis <- gsub('unknown|NA', NA, clin$age_diagnosis, ignore.case = TRUE)
  #########
  # function to split rows with multiple sample ids
  #########
  
  # HERE important note - write data and then read back in new version that has hand edits.
  # write_csv(clin, '../../Data/clin_data/temp_clin.csv')
  clin <- read.csv('../../Data/clin_data/temp_clin.csv')
  names(clin) <- gsub('.', '_', names(clin), fixed = TRUE)
  

  splitRows <- function(data, duplicate_table){
    
    for (ids in unique(clin$tm_donor)) {
      
      # get all diagnoses - if they have multiple first 4 columns are NA
      sub_data <- data %>% dplyr::filter(tm_donor %in% ids)
      
      if (any(grepl('/', sub_data$blood_dna_malkin_lab) | 
          grepl('/', sub_data$age_sample_collection))) {
        
        # see how many NAs relative to '/'
        dashes <- str_count(sub_data$blood_dna_malkin_lab, '/')[1]
        NAs <- (nrow(sub_data) -1)
        
        # split data
        split_malkin <- unlist(lapply(strsplit(as.character(sub_data$blood_dna_malkin_lab), '/'), 
                                      function(x) trimws(x, 'both')))
        split_age <- unlist(lapply(strsplit(as.character(sub_data$age_sample_collection), '/'), 
                                   function(x) trimws(x, 'both')))
        
        if(NAs == dashes){
          # remove NAs 
          split_malkin <- split_malkin[!is.na(split_malkin)]
          split_age <- split_age[!is.na(split_age)]
        }
        
        duplicate <- cbind.fill(split_malkin, split_age, sub_data)
        duplicate_table <- rbind(duplicate_table, duplicate)
      }
      
      print(ids)
    }
    data <- data[!grepl('/', data$blood_dna_malkin_lab, fixed = TRUE),]
    data <- data[!grepl('/', data$age_sample_collection, fixed = TRUE),]
    duplicate_table$blood_dna_malkin_lab <- NULL
    duplicate_table$age_sample_collection <- NULL
    colnames(duplicate_table)[1:2] <- c('blood_dna_malkin_lab', 'age_sample_collection')
    duplicate_table <- duplicate_table[, colnames(data)]
    data <- rbind(duplicate_table, data)
    return(data)
    
  }
  
  empty_table <- data.frame(matrix(ncol = ncol(clin), nrow = 0))
  clin <- splitRows(clin, empty_table)
  
  # Clean age of diagnosis column
  clin$age_diagnosis <- as.character(clin$age_diagnosis)
  clin$age_sample_collection <- as.character(clin$age_sample_collection)
  
  # vectorized cleaning 
  clin$age_diagnosis <- gsub('<', '', clin$age_diagnosis)

  # Clean the age variable so that is expressed in years
  
  cleanAge <-  function(data, column_name) {
    age_vector <- data[, column_name]
    for (i in 1:length(age_vector)) {
      temp_age <- age_vector[i]
      if(!grepl('[0-9]', temp_age)) {
        temp_age <- NA
      } else {
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
          temp_age <- gsub('12d', '0.032', temp_age)
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
  
  # clean date column and create variable for patient age 
  clin$dob <- ifelse(grepl('known', as.character(clin$dob)), NA,  as.character(clin$dob))
  
  # extract only the last four characters to get year of birth
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  clin$yob <- as.numeric(substrRight(as.character(clin$dob), 4))
  clin$age <- 2016 - clin$yob
  
  # clean p53
  clin$p53_germline <- as.character(clin$p53_germline)
  clin$p53_germline <- ifelse(is.na(clin$p53_germline), NA,
                              ifelse(grepl('WT', clin$p53_germline), 'WT',
                                     ifelse(grepl('Mut', clin$p53_germline), 'MUT', NA)))
  
  # create indicator for how many times an individual has had cancer 
  clin$cancer_diagnosis_diagnoses <- as.character(tolower(clin$cancer_diagnosis_diagnoses))
  
  clin$cancer_num <- ifelse(grepl('^p', clin$cancer_diagnosis_diagnoses), 
                            substr(clin$cancer_diagnosis_diagnoses, 2,2),
                            ifelse(grepl('unaffected', clin$cancer_diagnosis_diagnoses), '0',
                                   ifelse(is.na(clin$cancer_diagnosis_diagnoses), NA, '1')))
  
  # function to clean cancer_diagnoses
  clin$cancer_diagnosis_diagnoses <- gsub('?', '', clin$cancer_diagnosis_diagnoses, fixed = T)
  clin$cancer_diagnosis_diagnoses <- gsub(".*-","", clin$cancer_diagnosis_diagnoses)
  clin$cancer_diagnosis_diagnoses <- gsub(" ca","", clin$cancer_diagnosis_diagnoses, fixed = T)
  clin$cancer_diagnosis_diagnoses <- gsub("&",",", clin$cancer_diagnosis_diagnoses, fixed = T)
  clin$cancer_diagnosis_diagnoses <- gsub(";",",", clin$cancer_diagnosis_diagnoses, fixed = T)
  clin$cancer_diagnosis_diagnoses <- gsub("/",",", clin$cancer_diagnosis_diagnoses, fixed = T)
  clin$cancer_diagnosis_diagnoses <- gsub(" , ",",", clin$cancer_diagnosis_diagnoses, fixed = T)
  clin$cancer_diagnosis_diagnoses <- gsub("anaplastic","", clin$cancer_diagnosis_diagnoses, fixed = T)
  clin$cancer_diagnosis_diagnoses <- trimws(clin$cancer_diagnosis_diagnoses, which = 'both')
  clin$cancer_diagnosis_diagnoses <- ifelse(grepl('ffected', clin$cancer_diagnosis_diagnoses),
                                            'Unaffected',
                                            ifelse(grepl('adeno', clin$cancer_diagnosis_diagnoses), 
                                                   'adenocarcinoma', 
                                                   ifelse(grepl('os', clin$cancer_diagnosis_diagnoses), 
                                                          'os', 
                                                          ifelse(grepl('unknown|type not', clin$cancer_diagnosis_diagnoses),
                                                                 NA, 
                                                                 ifelse(grepl('hunting',clin$cancer_diagnosis_diagnoses),
                                                                        'colon',
                                                                        ifelse(grepl('paget', clin$cancer_diagnosis_diagnoses),
                                                                               'dcis;idc',clin$cancer_diagnosis_diagnoses))))))
  

  
  ##########
  # clean p53
  ##########
  clin$gender <- as.character(clin$gender)
  clin$gender <- ifelse(clin$gender == 'M', 'Male', 
                        ifelse(clin$gender == 'F', 'Female', 
                               ifelse(clin$gender == 'unknown', NA, clin$gender)))
  
  # replace LFS with family
  clin$family_name <- gsub('LFS', 'Family', clin$family_name)
  
  
  ##########
  # functon to clean family column 
  ##########
  clin$family_name <- ifelse(!grepl('Family', clin$family_name), 
                             NA, clin$family_name)
  clin$family_name<- unlist(lapply(strsplit(clin$family_name, '-'), function(x) x[1]))
  clin$family_name <- tolower(gsub(" ", "_", clin$family_name))
  
  # lauren's changes
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  
  # clean gdna VARS
  clin[,'gdna.exon.intron'] <- NA
  clin[,'gdna.exon.intron'][is.na(clin$gdna) == FALSE] <- 'no clear intron or exon'
  clin[,'gdna.exon.intron'] <- sub(":.*$", "", clin$gdna )
  clin[,'gdna.exon.intron'][grep(pattern = "Deletion",x = clin[,'gdna.exon.intron'])] <- 'no clear intron or exon'
  clin[,'gdna.exon.intron'][grep(pattern = "insertion",x = clin[,'gdna.exon.intron'])] <- 'no clear intron or exon'
  clin[,'gdna.exon.intron'][grep(pattern = "c.",x = clin[,'gdna.exon.intron'])] <- 'no clear intron or exon'
  clin[,'gdna.exon.intron'] <- trim(clin[,'gdna.exon.intron'])
  
  
  # clean base change
  clin[,'gdna.base.change'] <- NA
  clin[,'gdna.base.change'][is.na(clin$gdna) == FALSE] <- 'no bp change'
  clin[,'gdna.base.change'][grep(pattern = '>',x = clin$gdna)] <- unlist(lapply(X = strsplit(x = as.character(clin$gdna),split = '>'),
                                                                                FUN = function(x){if(length(x) == 2){
                                                                                  str1 <- substr(x[1],nchar(x[1]),nchar(x[1]))
                                                                                  str2 <- substr(x[2],1,1)
                                                                                  out.val <- paste0(str1,'>',str2)
                                                                                  return(out.val)
                                                                                }
                                                                                }))
  clin[,'gdna.base.change'] <- trim(clin[,'gdna.base.change'])
  
  # gdna.codon
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
  
  
  # clean protein vars
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
  
  
  clin$protein.codon.change <- unlist(lapply(strsplit(clin$protein.codon.change,split = " "),FUN = function(x){x[1]}))
  
  clin[,'protein.codon.num'] <- NA
  clin[,'protein.codon.num'][is.na(clin$protein) == FALSE] <- as.numeric(gsub("\\D", "", clin$protein[is.na(clin$protein) == FALSE]))
  unique(clin$protein.codon.num)
  
  clin[,'splice.delins.snv'] <- NA
  clin[,'splice.delins.snv'][grep(pattern = '>',clin$gdna.base.change)] <- 'SNV'
  clin[,'splice.delins.snv'][grep(pattern = 'deletion',clin$protein,ignore.case = T)] <- 'Deletion'
  clin[,'splice.delins.snv'][grep(pattern = 'splice',clin$protein,ignore.case = T)] <- 'Splice'
  clin[,'splice.delins.snv'][grep(pattern = 'dup',clin$gdna,ignore.case = T)] <- 'Duplication'
  clin$codon72.npro <- NA
  clin$codon72.npro[clin$codon72 == 'arg/arg'] <- 0
  clin$codon72.npro[clin$codon72 == 'arg/pro'] <- 1
  clin$codon72.npro[clin$codon72 == 'pro/pro'] <- 2
  
  # mdm2
  clin$mdm2.nG <- NA
  clin$mdm2.nG[clin$mdm2 == 'T/T'] <- 0
  clin$mdm2.nG[clin$mdm2 == 'T/G'] <- 1
  clin$mdm2.nG[clin$mdm2 == 'G/G'] <- 2
  
  
  write_csv(clin, '../../Data/clin_data/newest_clin.csv')
  
}


methy_dat <- load('../../Data/new_m_processed_temp.RData')
full_data <- rbind(data_cases_full,
                   data_controls_full)
clin_old <- read_csv('../../Data/clin_data/clinical_two.csv')
clin <- read_csv('../../Data/clin_data/newest_clin.csv')

# remove methylation 
full_data <- full_data[, c('tm_donor_', 'ids')]
names(full_data)[1] <- 'tm_donor'

# join data 
clin_full <- left_join(clin, full_data)

# change last variable to indicator method
names(clin_full)[ncol(clin_full)] <- 'have_methylation'
clin_full$have_methylation <- ifelse(is.na(clin_full$have_methylation), FALSE, TRUE)

# get new data identfiers 
new_data_ids <- as.data.frame(unique(clin$tm_donor)[!unique(clin$tm_donor) %in% clin_old$tm_donor_])
new_data_ids$new_var <- 'new_var'
colnames(new_data_ids)[1] <- 'tm_donor'


# get these from clinical and see if we have methylation for any
clin_full <- left_join(clin_full, new_data_ids)
names(clin_full)[ncol(clin_full)] <- 'recently_added_sample'
clin_full$recently_added_sample <- ifelse(is.na(clin_full$recently_added_sample), FALSE, TRUE)

write_csv(clin, '../../Data/clin_data/newest_clin.csv')








