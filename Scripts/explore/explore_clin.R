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
# get model data and correlate ages
##########
cases <- temp[!is.na(temp$age_diagnosis),]
cases <- cases[cases$methyl_indicator == TRUE,]

ggplot(cases, aes(age_diagnosis, age_sample_collection)) + geom_point() +
  xlab('') + ylab('') + xlim(c(0,800)) + ylim(c(0,800)) +  theme(text = element_text(size = 20),
                                                           axis.text = element_text(colour = 'black'),
                                                           panel.grid.major = element_line(colour = "grey", linetype = 3),
                                                           panel.background = element_rect(fill = "white"))
##########
  
##########
# Simple Pie Chart
##########
# remove nas from cancer_diagnosises

sub_clin <- clin[!is.na(clin$cancer_diagnosis_diagnoses),]
sub_clin <- sub_clin[!grepl('Unaffected', sub_clin$cancer_diagnosis_diagnoses),]

#recode sub_clin cancers
summary(as.factor(sub_clin$cancer_diagnosis_diagnoses))

# write.csv(sub_clin, '/home/benbrew/Desktop/sub_clin.csv')
sub_clin_new <- read.csv('/home/benbrew/Desktop/sub_clin_new.csv')

##########
# trim columns
##########
sub_clin$cancer_diagnosis_diagnoses <- trimws(sub_clin$cancer_diagnosis_diagnoses, which = 'both')
sub_clin$p53_germline <- trimws(sub_clin$p53_germline, which = 'both')

# group by cancer_diagnosis
dat_mod <- sub_clin_new %>%
  group_by(cancer_diagnosis_diagnoses) %>%
  summarise(counts = n())

# get color column in dat mod that is # of cancers (11)
dat_mod$col <- c('grey', 'red', 'blue', 'green', 'orange', 'yellow', 'purple',
                 'white', 'black', 'lightblue', 'darkgreen')


# first recode for new cancer data
pie(dat_mod$counts, dat_mod$cancer_diagnosis_diagnoses, main="Distribution of cancers",
    col = adjustcolor(dat_mod$col, alpha.f =0.8))


#########
# get distribution of age
#########
# get distribution of age 
sub_clin$age_diagnosis <- as.numeric(sub_clin$age_diagnosis)
hist(sub_clin$age_diagnosis, main = 'Age of onset', xlab = 'Age (in months)', 
     col = adjustcolor('blue'), alpha.f = 0.8)


# at 48
qplot(sub_clin$age_diagnosis, geom="histogram", bins = 25, fill=I("blue"), 
      col=I("black"), alpha=I(.8)) + xlab('Age of onset in months') +
  ylab('Frequency') + theme_bw() + theme(text = element_text(size = 15)) 

# at 48
qplot(sub_clin$age_sample_collection, geom="histogram", bins = 25, fill=I("blue"), 
      col=I("black"), alpha=I(.2),) + xlab('Age of sample collection') +
  ylab('Frequency') + theme_bw() + theme(text = element_text(size = 15)) 




##########
# generate random heatmap
##########

nba_heatmap <- heatmap(nba_matrix, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))