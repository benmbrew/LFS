##################################################################################################
# This scripts will get ids to be sent for a validation test. That is, samples which are most 
# similar to the samples we already have, but don't have methylation. 
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
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = FALSE)

#################################
# get clinical data that is being used in model
#################################
# subset to mutant, with cancer, no methyl 
model_data <- clin[clin$cancer_diagnosis_diagnoses != 'Unaffected' & clin$p53_germline == 'Mut' & 
                     clin$methyl_indicator == 'Yes',]

# remove duplicates 
model_data <- model_data[!duplicated(model_data$blood_dna_malkin_lab_),]

model_data <- model_data[!duplicated(model_data$tm_donor_),]

# remove NA 
model_data <- model_data[!is.na(model_data$age_diagnosis),]

###################################
# set criteria for validators
####################################
temp <- clin[clin$p53_germline == 'Mut' & clin$cancer_diagnosis_diagnoses != 'Unaffected' 
            & clin$methyl_indicator == 'No',]

# remove duplicates 
temp <- temp[!duplicated(temp$blood_dna_malkin_lab_),]

temp <- temp[!duplicated(temp$tm_donor_),]

# remove NA 
temp <- temp[!is.na(temp$age_diagnosis),]

validators <- temp

###################################
# subset validators to 48 and fit the population characteristics 
# similar to model_data
###################################
range <- 18
# add an indicator in validators columns to show it is a part of validators 
colnames(validators) <- paste0(colnames(validators), '_', 'validators')

# first need to find which validators in model_data crit1 was close to in age of sample collection
# create column in model_data that gives begining and end of age diagnosis 
model_data$beg <- model_data$age_sample_collection - range
model_data$end <- model_data$age_sample_collection + range

# use sql data table to merger validators with model_data based on age range
result = sqldf("select * from validators
                left join model_data
                on validators.age_sample_collection_validators between model_data.beg and model_data.end")


result <- result[!duplicated(result$blood_dna_malkin_lab__validators),]
result <- result[!is.na(result$blood_dna_malkin_lab_),]


temp_r <- cbind(result$blood_dna_malkin_lab__validators, result$age_sample_collection_validators, 
               result$blood_dna_malkin_lab_, result$beg, result$end)

####################################
# now find other clinical variables 
####################################

# first look at amount of NAs in each varible
length(which(is.na(result$gdna.exon.intron_validators))) #0 
length(which(is.na(result$gdna.base.change_validators)))# 0
length(which(is.na(result$gdna.codon_validators))) # 0
length(which(is.na(result$protein.codon.change_validators))) # 0
length(which(is.na(result$splice.delins.snv_validators))) # 10, 18
length(which(is.na(result$protein.codon.num_validators))) # 44, 31
length(which(is.na(result$codon72.npro_validators))) # 51, 183
length(which(is.na(result$mdm2.nG_validators))) # 129, 243

# leave out gnda.codon and protein.codon.num because they are exact locations and difficult to match

# temp <- result[result$gender_validators == result$gender & 
#               result$gdna.base.change_validators == result$gdna.base.change &
#               result$gdna.codon_validators == result$gdna.codon & 
#               result$protein.codon.change_validators == result$protein.codon.change &
#               result$splice.delins.snv_validators == result$splice.delins.snv &
#               result$protein.codon.num_validators == result$protein.codon.num &
#               result$codon72.npro_validators == result$codon72.npro & 
#               result$mdm2.nG_validators == result$mdm2.nG,]

# match on gender, gdna.base.change, protein.codon.change, splice.delins.sng
final <- result[result$gender_validators == result$gender,]

final <- cbind(id_validators = final$blood_dna_malkin_lab__validators, 
               age_sample_validators = final$age_sample_collection_validators, 
                id_model_datat = final$blood_dna_malkin_lab_, range_beg = final$beg, range_end = final$end)

final_ids <- final[, 1]

write.csv(final_ids, '/home/benbrew/Desktop/validation_ids.csv')

