
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

write.csv(ids_new, file = '/home/benbrew/Desktop/ids.csv')

