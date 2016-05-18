library(dplyr)
######################################### Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/classification_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

# Load in clinical data
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = FALSE)

# get counts for cancer
clin$age <- as.numeric(clin$age)
cancer <- clin %>%
  filter(!is.na(cancer_diagnosis_diagnoses)) %>%
  filter(!is.na(p53_germline)) %>%
  group_by(cancer_diagnosis_diagnoses, p53_germline) %>%
  summarise(mean_age = mean(age, na.rm = T),
            counts = n())

# order counts
cancer <- cancer[order(cancer$counts, decreasing = TRUE),]
  