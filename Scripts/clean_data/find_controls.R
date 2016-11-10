
###############################################
# This script will search through clinical data for samples for which we dont have methylation data for 
# that are close in age of sample collection to samples that we have methylation for. AND that do not have cancer yet,
# but are p53 Mutants.
# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/Analyze')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

library(sqldf)
library(dplyr)

#################################################################################################
# Read in clinical data and join by ids
#################################################################################################

# Read in data (clinical or clinical_two)
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = TRUE)

###################################################################
# subset to just three columns- id, methlation_indicator, age_sample_collection, p53_germline, cancer_diagnosis
dat <- clin[, c('p53_germline', 'cancer_diagnosis_diagnoses', 'age_sample_collection', 'blood_dna_malkin_lab_',
                'methyl_indicator')]

temp <- dat[dat$p53_germline == 'Mut' & dat$cancer_diagnosis_diagnoses == 'Unaffected' 
            & dat$methyl_indicator == 'No',]

temp1 <- temp[complete.cases(temp),]
# dat <- dat[!is.na(dat$cancer_diagnosis_diagnoses),]
# & dat$cancer_diagnosis_diagnoses[i] != 'Unaffected'
# 131 yes in methylation indicator
range <- 50
samples <- list()

for ( i in 1:nrow(dat) ) {
  
  if (dat$methyl_indicator[i] == 'Yes') {
    
      temp <- dat$age_sample_collection[i]
      
      
      samples[[i]] <- dat[(dat$age_sample_collection > temp & dat$age_sample_collection < (temp + range)) |
                            (dat$age_sample_collection < temp & dat$age_sample_collection > (temp - range)),] 
  } else {
    
    samples[[i]] <- NULL
  }
  
}


samples <- do.call('rbind', samples)
samples <- samples[!duplicated(samples$blood_dna_malkin_lab_),]

# subset to mutant, unaffected, no methyl 
ids_new <- samples[samples$cancer_diagnosis_diagnoses == 'Unaffected' & samples$p53_germline == 'Mut' & 
                 samples$methyl_indicator == 'No',]

ids_new <- ids_new[complete.cases(ids_new),]

# write.csv(ids_new, file = '/home/benbrew/Desktop/ids.csv')

############################################
# Now search through clinical data for samples that we have methylation for
# that are close in age of sample collection to samples we sent Tany AND that have cancer, and are p53 Mut
############################################
# criteria1: p53= Mut, methylation_indicator = yes, cancer_diagnosis != 'Unaffected'
# criteria2: As close as possible to samples we sent Tanya - age of sample collection - then any other clinical variables 

# first remove samples 3880, 3879, 2320, as well as 5064, 4176 from ids_new
ids_new <- ids_new[!grepl('3880|3879|2320|5064|4176', ids_new$blood_dna_malkin_lab_),]

# rejoin ids_new with clinical data to get full sample characteristics
profile_data <- inner_join(ids_new, clin, by = 'blood_dna_malkin_lab_')

# now subset clin by criteria1 : p53= Mut, methylation_indicator = yes, cancer_diagnosis = 'Unaffected'
crit1 <- clin[clin$p53_germline == 'Mut' & clin$methyl_indicator == 'Yes' & clin$cancer_diagnosis_diagnoses != 'Unaffected',]

# remove rows in crit1 that dont have age of sample collection
crit1 <- crit1[!is.na(crit1$age_sample_collection),]


# now search through crit1 and keep samples that are within 50 months of profile data in terms of age of sample collection
# use similar for loop from above
range <- 22
samples <- list()

for ( i in 1:nrow(profile_data) ) {
  

  temp <- profile_data$age_sample_collection.x[i]
  
  
  samples[[i]] <- crit1[(crit1$age_sample_collection > temp & crit1$age_sample_collection < (temp + range)) |
                        (crit1$age_sample_collection < temp & crit1$age_sample_collection > (temp - range)),] 

 print(i)
}


samples <- do.call('rbind', samples)
samples <- samples[!duplicated(samples$blood_dna_malkin_lab_),]

##################################
# Now match crit1 to tanya's data (profile_data) as well as possible 
##################################

# add an indicator in samples columns to show it is a part of samples 
colnames(samples) <- paste0(colnames(samples), '_', 'samples')

# first need to find which samples in profile_data crit1 was close to in age of sample collection
# create column in profile_data that gives begining and end of age diagnosis 
profile_data$beg <- profile_data$age_sample_collection.x - range
profile_data$end <- profile_data$age_sample_collection.x + range

# use sql data table to merger samples with profile_data based on age range
result = sqldf("select * from samples
                left join profile_data
                on samples.age_sample_collection_samples between profile_data.beg and profile_data.end")

####################################
# now find other clinical variables 
####################################

# first look at amount of NAs in each varible
length(which(is.na(result$gdna.exon.intron_samples))) #0 
length(which(is.na(result$gdna.base.change_samples)))# 0
length(which(is.na(result$gdna.codon_samples))) # 0
length(which(is.na(result$protein.codon.change_samples))) # 0
length(which(is.na(result$splice.delins.snv_samples))) # 10, 18
length(which(is.na(result$protein.codon.num_samples))) # 44, 31
length(which(is.na(result$codon72.npro_samples))) # 51, 183
length(which(is.na(result$mdm2.nG_samples))) # 129, 243

# leave out gnda.codon and protein.codon.num because they are exact locations and difficult to match

# temp <- result[result$gender_samples == result$gender & 
#               result$gdna.base.change_samples == result$gdna.base.change &
#               result$gdna.codon_samples == result$gdna.codon & 
#               result$protein.codon.change_samples == result$protein.codon.change &
#               result$splice.delins.snv_samples == result$splice.delins.snv &
#               result$protein.codon.num_samples == result$protein.codon.num &
#               result$codon72.npro_samples == result$codon72.npro & 
#               result$mdm2.nG_samples == result$mdm2.nG,]

# match on gender, gdna.base.change, protein.codon.change, splice.delins.sng
final <- result[result$gender_samples == result$gender & 
                 result$gdna.base.change_samples == result$gdna.base.change &
                 result$protein.codon.change_samples == result$protein.codon.change &
                 result$splice.delins.snv_samples == result$splice.delins.snv,]

# keep only matched variables and criteria variables to double check
final <- final[, c('blood_dna_malkin_lab__samples', 'blood_dna_malkin_lab_' ,'p53_germline_samples', 'p53_germline.x', 'gender_samples', 'gender', 'gdna.base.change_samples',
                   'gdna.base.change', 'protein.codon.change_samples', 'protein.codon.change',
                   'splice.delins.snv_samples', 'splice.delins.snv')]

"2263" "2318" "3432" "2257"
