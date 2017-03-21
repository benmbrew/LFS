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
# temp.counts 
##########
temp.dup <- clin[!duplicated(clin$id),]
temp.counts <- temp.dup %>%
  filter(!is.na(id)) %>%
  group_by(p53_germline) %>%
  summarise(counts = n())

##########
# temp.age 
##########
temp.dup <- clin[!duplicated(clin$id),]
temp.age <- temp.dup %>%
  filter(!is.na(p53_germline)) %>%
  group_by(p53_germline, cancer_diagnosis_diagnoses) %>%
  summarise(mean_age = mean(age_diagnosis, na.rm = T))

##########
# temp.meth
##########
temp.dup <- clin[!duplicated(clin$id),]
temp.methyl <- temp.dup %>%
  filter(!is.na(p53_germline)) %>%
  group_by(p53_germline, methyl_indicator) %>%
  summarise(counts = n())

##########
# temp.unaff
##########
temp.dup <- clin[!duplicated(clin$id),]
temp.unaff <- temp.dup %>%
  filter(!is.na(p53_germline)) %>%
  group_by(p53_germline, cancer_diagnosis_diagnoses) %>%
  summarise(counts = n())



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


##########
# make cancer indicator
##########
clin$cancer <- ifelse(grepl('Unaffected', clin$cancer_diagnosis_diagnoses), FALSE, TRUE)

summary(as.factor(clin$p53_germline))

# subset 
temp <- clin[which(clin$p53_germline == 'Mut'),]

summary(as.factor(clin$methyl_indicator)[which(temp$p53_germline == 'Mut')])








##########
# Simple Pie Chart
##########
# get model data
dat_mod <- clin[clin$methyl_indicator == 'TRUE' & clin$p53_germline == 'Mut',]
dat_mod <- dat_mod[!is.na(dat_mod$age_diagnosis),]

# write.csv(dat_mod, '/home/benbrew/Desktop/model_temp.csv')
dat_mod <- read.csv('/home/benbrew/Desktop/model_temp.csv')

# make cancer indeicator
dat_mod$cancer <- ifelse(grepl('Unaffected', dat_mod$cancer_diagnosis_diagnoses), FALSE, TRUE)
# 
# # group by cancer_diagnosis
# dat_mod <- dat_mod %>%
#   group_by(cancer_diagnosis_diagnoses) %>%
#   summarise(counts = n())


# first recode for new cancer data

pie(dat_mod$counts, dat_mod$cancer_diagnosis_diagnoses, main="Distribution of cancers")

# get distribution of age 
dat_mod$age_diagnosis <- as.numeric(dat_mod$age_diagnosis)
hist(dat_mod$age_diagnosis, main = 'Age of onset', xlab = 'Age (in months)', col = 'lightblue')


# at 48
qplot(dat_mod$age_diagnosis, geom="histogram", bins = 25, fill=I("blue"), 
      col=I("black"), alpha=I(.2),) + xlab('Age of onset') +
  ylab('Frequency') + theme_bw() + theme(text = element_text(size = 15)) + geom_vline(xintercept = 72, 
                                                                                      colour = 'red')


# at 48
qplot(dat_mod$age_sample_collection, geom="histogram", bins = 25, fill=I("blue"), 
      col=I("black"), alpha=I(.2),) + xlab('Age of sample collection') +
  ylab('Frequency') + theme_bw() + theme(text = element_text(size = 15)) 
+ geom_vline(xintercept = 72, 
                                                                                      colour = 'red')




##########
# group by 
##########