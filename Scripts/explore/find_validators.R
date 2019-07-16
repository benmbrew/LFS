##################################################################################################
# This scripts will get ids to be sent for a validation test. That is, samples which are most 
# similar to the samples we already have, but don't have methylation. 
# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
imputed_data <- paste0(data_folder, '/imputed_data')
idat_data <- paste0(methyl_data, '/raw_files')
model_data <- paste0(data_folder, '/model_data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')
clin_data <- paste0(data_folder, '/clin_data')

library(sqldf)
library(dplyr)

#################################################################################################
# Read in clinical data and join by ids
#################################################################################################
load(paste0(model_data, '/model_data_cases.RData')) #check if age is in id column - may need to fix clean clin
# # ids are ok
rm(beta_funnorm, beta_illumina, beta_quan, beta_raw, beta_swan, cg_locations, rgSetList)
# 
# load(paste0(idat_data, '/imputed_idat_betas.RData')) # use this to get ids we have methylation data for
# methyl <- as.data.frame(beta_raw[1:229, 1:10])
# rm(gene_knn, gene_lsa, probe_knn, probe_lsa)
# 
# save.image(paste0(idat_data, '/find_validators_data.RData')) 

#################################
# Take only necessary data and save image
#################################
# load(paste0(idat_data, '/find_validators_data.RData')) 
methyl <- beta_funnorm
# clean ids in each data set 
cleanIDs <- function(data){
  
  data <- as.data.frame(data)
  data$id <- gsub('A|B|_', '', data$id)
  data$id <- substr(data$id, 1,4) 
  return(data)
}

methyl <- cleanIDs(methyl)

#################################
# for loop to add a yes or no indicator in clin on whether or not mehtylation data exists for each patient
#################################
intersected_ids <- intersect(methyl$id, clin$id)
clin$methyl_fac <- NA
methyl$methyl_fac <- 'yes'

# loop 
for(i in intersected_ids) {
  clin$methyl_fac[which(clin$id == i)] <- methyl$methyl_fac[which(methyl$id == i)]
}

# change NA to no
clin$methyl_fac[is.na(clin$methyl_fac)] <- 'no'

#################################
# get clinical data that is being used in model
#################################
# subset to mutant, with cancer, no methyl 
model_data <- clin[clin$cancer_diagnosis_diagnoses != 'Unaffected' & clin$p53_germline == 'Mut' & 
                     clin$methyl_fac == 'yes',]

# remove duplicates 
model_data <- model_data[!duplicated(model_data$blood_dna_malkin_lab_),]

model_data <- model_data[!duplicated(model_data$tm_donor_),]

# remove NA 
model_data <- model_data[!is.na(model_data$age_diagnosis),]

###################################
# set criteria for validators
####################################
temp <- clin[clin$p53_germline == 'Mut' & clin$cancer_diagnosis_diagnoses != 'Unaffected' 
            & clin$methyl_fac == 'no',]

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
range <- 36
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

###################################
# check distribution of validat
hist(model_data$age_sample_collection)
hist(result$age_sample_collection_validators)

# check cancers
summary(as.factor(model_data$cancer_diagnosis_diagnoses))
summary(as.factor(result$cancer_diagnosis_diagnoses))

# check genders
summary(as.factor(model_data$gender))
summary(as.factor(result$gender))

write.csv(temp, '/home/benbrew/Desktop/validators.csv' )

# ####################################
# # now find other clinical variables 
# ####################################
# 
# # first look at amount of NAs in each varible
# length(which(is.na(result$gdna.exon.intron_validators))) #0 
# length(which(is.na(result$gdna.base.change_validators)))# 0
# length(which(is.na(result$gdna.codon_validators))) # 0
# length(which(is.na(result$protein.codon.change_validators))) # 0
# length(which(is.na(result$splice.delins.snv_validators))) # 10, 18
# length(which(is.na(result$protein.codon.num_validators))) # 44, 31
# length(which(is.na(result$codon72.npro_validators))) # 51, 183
# length(which(is.na(result$mdm2.nG_validators))) # 129, 243
# 
# 
# ##############
# # check age of sample collection distribution with model_data
# ##############
# 
# result_list <- as.data.frame(result$id_validators)
# 
# 
# # leave out gnda.codon and protein.codon.num because they are exact locations and difficult to match
# 
# # temp <- result[result$gender_validators == result$gender & 
# #               result$gdna.base.change_validators == result$gdna.base.change &
# #               result$gdna.codon_validators == result$gdna.codon & 
# #               result$protein.codon.change_validators == result$protein.codon.change &
# #               result$splice.delins.snv_validators == result$splice.delins.snv &
# #               result$protein.codon.num_validators == result$protein.codon.num &
# #               result$codon72.npro_validators == result$codon72.npro & 
# #               result$mdm2.nG_validators == result$mdm2.nG,]

# # match on gender, gdna.base.change, protein.codon.change, splice.delins.sng
# final <- result[result$gender_validators == result$gender,]
# 
# final <- cbind(id_validators = final$blood_dna_malkin_lab__validators, 
#                age_sample_validators = final$age_sample_collection_validators, 
#                 id_model_datat = final$blood_dna_malkin_lab_, range_beg = final$beg, range_end = final$end)
# 
# final_ids <- final[, 1]
# 
# write.csv(final_ids, '/home/benbrew/Desktop/validation_ids.csv')

