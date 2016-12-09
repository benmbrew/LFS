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

##########
# read in data 
##########
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

##########
# trim columns
##########
clin$cancer_diagnosis_diagnoses <- trimws(clin$cancer_diagnosis_diagnoses, which = 'both')
clin$p53_germline <- trimws(clin$p53_germline, which = 'both')


##########
# find mean age of diagnosis among mut and WT
##########
temp <- clin %>%
  group_by(p53_germline, cancer_diagnosis_diagnoses) %>%
  summarise(mean_age = mean(age_diagnosis, na.rm = T),
            counts = n())

temp.1 <- clin %>%
  filter(!is.na(p53_germline)) %>%
  group_by(p53_germline) %>%
  summarise(counts = n())
