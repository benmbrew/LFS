# initialize folders
library(gsubfn)
librarty(tidyverse)



path_to_clin <- '../../Data/clin_data'
##########
# read in data 
##########
clin <- read_csv(paste0(path_to_clin, '/clinical_two.csv'))

##########
# trim columns
##########
clin$cancer_diagnosis_diagnoses <- trimws(clin$cancer_diagnosis_diagnoses, which = 'both')
clin$p53_germline <- trimws(clin$p53_germline, which = 'both')

##########
# get summary statisctics for clinical data
##########

# keep only necessary columns 
keep_cols <- c('tm_donor_', 'ids', 'family_name', 'relationship', 'gender', 'p53_germline', 'cancer_diagnosis_diagnoses', 'age_diagnosis', 'age_sample_collection',
               'gdna.exon.intron' ,'gdna.base.change', 'gdna.codon', 'protein.codon.change', 'protein.codon.num')
clin_sub <- clin[, keep_cols]

length(which(duplicated((clin$family_name))))

##########
# age of onset chart for p53 status 
##########
p53_age <- clin %>%
  filter(!is.na(p53_germline)) %>%
  group_by(p53_germline) %>%
  summarise(mean_age_onset = (mean(age_diagnosis, na.rm = T)))

temp <- clin[!is.na(clin$p53_germline),]

ggplot(temp, aes(age_diagnosis, fill = p53_germline)) + 
  geom_histogram(bins = 20, alpha=.5, position="identity") + 
  xlab('Age of diagnosis') + ylab('Frequency in clinical data') + 
  scale_fill_manual(name = 'TP53 Germline',
                    breaks = c('Mut', 'WT'),
                    labels = c('Mutation', 'Wild Type'),
                    values = c('Red', 'Blue'))  + 
  theme_bw() + theme(axis.text.y = element_text(size = 12),
                     axis.text.x = element_text(size = 12))

##########
# get top cancers age distribution 
##########
# group by cancer and get counts and mean age onset
cancer_age <- clin %>%
  filter(!is.na(cancer_diagnosis_diagnoses)) %>%
  filter(cancer_diagnosis_diagnoses != 'Unaffected') %>%
  group_by(cancer_diagnosis_diagnoses) %>%
  summarise(counts = n(),
            mean_age = mean(age_diagnosis, na.rm = T)) %>%
  arrange(-counts)

cancer_age_names <- paste(cancer_age$cancer_diagnosis_diagnoses[1:5], collapse = '|')
cancer_age_names <- "^ACC$|^CPC$|^OS$|^ERMS$|^Breast ca$"

temp_clin <- clin[grepl(cancer_age_names, clin$cancer_diagnosis_diagnoses),]

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")

ggplot(temp_clin, aes(age_diagnosis, fill = cancer_diagnosis_diagnoses)) + 
  geom_histogram(bins = 30, alpha=.5, position="identity") + 
  xlab('Age of diagnosis') + ylab('Frequency in clinical data') + 
  scale_fill_manual(name = 'Most prevelant cancers',
                    breaks = c('ACC', 'CPC', 'OS', 'ERMS', 'Breast ca'),
                    labels = c('ACC', 'CPC', 'OS', 'ERMS', 'Breast'),
                    values = cbPalette)  + 
  theme_bw() + theme(axis.text.y = element_text(size = 12),
                     axis.text.x = element_text(size = 12))


# subset to acc and do it on lfs
temp_acc <- temp_clin[grepl('ACC', temp_clin$cancer_diagnosis_diagnoses),]

ggplot(temp_acc, aes(age_diagnosis, fill = p53_germline)) + 
  geom_histogram(bins = 15, alpha=.5, position="identity") + 
  xlab('Age of diagnosis') + ylab('Frequency in clinical data') + 
  scale_fill_manual(name = 'TP53 Germline',
                    breaks = c('Mut', 'WT'),
                    labels = c('Mutation', 'Wild Type'),
                    values = c('Red', 'Blue')) + 
  theme_bw() + theme(axis.text.y = element_text(size = 12),
                     axis.text.x = element_text(size = 12))
##########
# get model data and check reoccuring cancers
##########

# first remove unaffected 
temp_mod <- clin[!grepl('Unaffected', clin$cancer_diagnosis_diagnoses),]

# ones with have methylation for
temp_mod <- temp_mod[temp_mod$methyl_indicator == TRUE,]

# keep tp53 mutant
temp_mod <- temp_mod[grepl('Mut', temp_mod$p53_germline),]

# figure out duplicates and reoccuring
length(which(duplicated(temp_mod$blood_dna_malkin_lab_)))

##########
# how many familes
##########
temp.fam <- clin[!is.na(clin$family_name),]
length(unique(temp.fam))

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
# 

# trim columns
##########
sub_clin$cancer_diagnosis_diagnoses <- trimws(sub_clin$cancer_diagnosis_diagnoses, which = 'both')
sub_clin$p53_germline <- trimws(sub_clin$p53_germline, which = 'both')

# group by cancer_diagnosis
dat_mod <- sub_clin %>%
  group_by(cancer_diagnosis_diagnoses) %>%
  summarise(counts = n())

dat_mod <- dat_mod[order(dat_mod$counts, decreasing = T),]
write.csv(dat_mod, '/home/benbrew/Desktop/dat_mod.csv')
# sub_clin_new <- read.csv('/home/benbrew/Desktop/sub_clin_new.csv')


##########
library(RColorBrewer)
dat_mod <- read.csv('/home/benbrew/Desktop/dat_mod.csv')
summary(dat_mod$name)

display.brewer.all()

# get color column in dat mod that is # of cancers (11)
cbPalette <- brewer.pal(7,"Set1")
dat_mod$col <- cbPalette


dat_mod$name <- sub(' ', '\n', dat_mod$name)
# first recode for new cancer data
pie(dat_mod$counts, dat_mod$name, main="Distribution of cancers",
    col = adjustcolor(dat_mod$col, alpha.f =0.6), cex = 1.2)

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