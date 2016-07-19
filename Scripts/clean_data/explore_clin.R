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

# read in data 
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# potential questions
# subsetting by p53 status Mut ?
# all Mut, unaffected WT are relatives
# fields : family_indicator, relationship to proband,
# tm_donor, p53_germline, cancer, age_diagnosis, malkin id, 
# age of sample collection, gdna, protein, codon_72, dob, gender

###########################################################
# group by p53 and cancer to see if any WT dont have cancer
temp <- clin %>%
  group_by(p53_germline, cancer_diagnosis_diagnoses) %>%
  summarise(counts = n())

# who are the 82 Mut who are unaffected 
temp <- clin[clin$p53_germline == 'Mut' & 
               clin$cancer_diagnosis_diagnoses == 'Unaffected',]

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

# get distribution of family size
temp <- clin %>%
  filter(!is.na(family_name)) %>%
  group_by(family_name) %>%
  summarise(counts = n())

hist(temp$counts)

mean(temp$counts)
median(temp$counts)

# difference between age of diagnosis and age of sample collection 
plot(clin$age_diagnosis, clin$age_sample_collection)

# avg difference bettween age diagnosis and age sample collection
diff <- clin$age_sample_collection - clin$age_diagnosis
mean(diff, na.rm = T)
hist(diff)
# 4 years

