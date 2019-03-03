# initialize folders
library(gsubfn)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(plotly)
# validation or combined data 
data_used <- 'new'

# get cg regions
cg_gene_regions = 'Body'

# set preprocessing method
method <- 'noob'

# set type of data, beta or m
methyl_type <- 'beta'

# set data directory
data_dir <- '../../Data/'

load(paste0(data_dir,paste0(data_used,'_',methyl_type, '_final_beta_first_last', '.RData')))

# source all_functions.R to load libraries and my functions
source('all_functions.R')

clin <- read.csv('../../Data/clin_data/new_clin.csv')

# keep only Mut, Cancer, not NA in age, and have methylation
clin <- clin[clin$methyl_status == 'yes',]
clin <- clin[clin$p53 == 'MUT',]
clin <- clin[clin$cancer_diagnosis_diagnoses == 'Unaffected',]
# clin <- clin[!is.na(clin$age_diagnosis),]
clin <- clin[!duplicated(clin$tm_donor),]
clin <- clin[!duplicated(clin$blood_dna_malkin_lab),]
clin <- clin[!is.na(clin$age_sample_collection),]


summary(as.factor(clin$cancer_diagnosis_diagnoses))
clin$cancer_diagnosis_diagnoses <- as.character(clin$cancer_diagnosis_diagnoses)
clin$cancer_diagnosis_diagnoses <- ifelse(!grepl('ACC|Breast|Erms|OS|CPC', clin$cancer_diagnosis_diagnoses), 
                                          'Other', clin$cancer_diagnosis_diagnoses)


# read in wt data
data_wt <- readRDS(paste0(data_dir,paste0('new','_',methyl_type, '_wild_type', '.rda')))

full_data_first <- full_data_first[full_data_first$tm_donor_ != '3955',]
full_data_last <- full_data_last[full_data_last$tm_donor_ != '3955',]

# bar chart for p53_germline and gender, with mean age of onset, sample collection,
clin$age_diagnosis <- round((clin$age_diagnosis/12),2)
clin$age_sample_collection <- round((clin$age_sample_collection/12),2)
full_data_first$age_diagnosis <- round((full_data_first$age_diagnosis/12),2)
full_data_first$age_sample_collection <- round((full_data_first$age_sample_collection/12),2)



temp <- clin %>% filter(!is.na(p53_germline)) %>% 
  filter(p53_germline != 'no result') %>% group_by(p53_germline) %>% 
  summarise(mean_onset = mean(age_diagnosis, na.rm = TRUE),
            mean_sample = mean(age_sample_collection, na.rm = TRUE),
            counts = n(),
            per_cancer = round(sum(!grepl('Unaffected', cancer_diagnosis_diagnoses))/counts, 2),
            per_female = round(sum(!grepl('M', gender))/counts, 2))



clin$cancer_diagnosis_diagnoses <- gsub('Anaplastic ERMS', 'ERMS', clin$cancer_diagnosis_diagnoses)
clin$cancer_diagnosis_diagnoses <- gsub('Breast ca', 'Breast', clin$cancer_diagnosis_diagnoses)

temp <- clin %>% filter(!is.na(cancer_diagnosis_diagnoses)) %>%
  group_by(cancer_diagnosis_diagnoses) %>% 
  summarise(counts = n()) %>% arrange(-counts) %>% filter(!grepl('Unaffected', cancer_diagnosis_diagnoses))

temp$cancer_cat <- ifelse(temp$counts < 9, 'Other', temp$cancer_diagnosis_diagnoses)

temp_1<- temp %>% group_by(cancer_cat) %>%
  summarise(sum_counts = sum(counts, na.rm = T)) %>% arrange(-sum_counts)

library(RColorBrewer)
library(wordcloud)
temp_1 <- temp_1[1:10,]

ggplot(temp_1, aes(x= '', y = sum_counts, fill = cancer_cat)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +  scale_fill_brewer(name = '',palette="Blues")+
  theme_minimal(base_size = 20, base_family = 'Ubuntu') + 
  labs(x = '', y = '')

temp_p

pie(temp_1$sum_counts, labels = temp_1$cancer_cat, col = brewer.pal(8, 'Set1'))

temp_cancer <- clin[!grepl('Unaffected', clin$cancer_diagnosis_diagnoses),]
wordcloud(temp_cancer$cancer_diagnosis_diagnoses, random.order = FALSE,rot.per=0.35,
          colors = brewer.pal(8, "Dark2"), use.r.layout=TRUE)



# get pie chart for Mut and WT and unknown 
temp_p53 <- 
  clin %>% 
  group_by(p53_germline) %>%
  summarise(counts = n())

temp_p53$p53_germline <- gsub('no result', NA, temp_p53$p53_germline)

pie(temp_p53$counts, labels = temp_p53$p53_germline, col = brewer.pal(8, 'Set1'))
temp_p53 <- temp_p53[complete.cases(temp_p53),]
temp_p53$counts<- ifelse(temp_p53$counts == 280, 277, temp_p53$counts)



f <- list(
  family = "Ubuntu",
  size = 20,
  color = "white"
)

colors <- c('rgb(171,104,87)', 'rgb(114,147,203)')

p1 <-  plot_ly(temp_p53,labels = ~p53_germline, values = ~counts,
               type ='pie',
               textposition = 'inside',
               textinfo = 'label+percent+value',
               insidetextfont = f,
               hoverinfo = 'label+value',
               marker = list(colors = colors,
                             line = list(color = '#FFFFFF', width = 1)))  %>%
  
  config(displayModeBar = F) %>%
  
  layout(title = 'LFS', showlegend = F,
         annotations = list(
           showarrow = FALSE,
           text = '',
           font = list(color = '#1F2023',
                       family = 'sans serif')), 
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


ggplot(temp_p53, aes(x= '', y = counts, fill = p53_germline)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +  scale_fill_brewer(palette="Dark2")+
  theme_minimal(base_size = 12, base_family = 'Ubuntu') + 
  labs(x = '', y = '')


# plot age diagnosis against age of sample collection
temp_mod <- full_data_first[!duplicated(full_data_first$tm_donor_), 1:17]

ggplot(temp_mod, aes(age_diagnosis, age_sample_collection)) + 
  geom_point(size = 4, col = 'black', alpha = 0.6) +
  labs(x = 'Age of onset', y = 'Age of sample collection') +
  theme_minimal(base_size = 18, base_family = 'Ubuntu')  
  


# temp_dat <- clin
# temp_dat <- full_data_first[, 1:15]
# var1 = 'gender'
# var2 = 'p53_germline'
# var3=  'age_sample_collection'
# title = 'P53'
get_bar_plot <- function(temp_dat, var1, var2, var3, title) {
  
  temp <- temp_dat[, c(var1, var2, var3)]
  
  orig_labs <- names(temp)[1:3]
  orig_labs <- gsub('_', ' ', orig_labs)
  orig_labs <- Hmisc::capitalize(orig_labs)
  
  names(temp)[1:3] <- c('V1', 'V2', 'V3')
  
  temp <- 
    temp %>%
    group_by(V1, V2) %>%
    summarise(mean_value = round(mean(V3, na.rm = T),2),
                                counts = n())
  
  temp <- temp[complete.cases(temp),]
  
  
  cols <- brewer.pal(8, name = 'Set1')
  
  
  p <- ggplot(temp,
              aes(V1, mean_value, 
              fill = V2)) + geom_bar(stat = 'identity', 
                                    position = 'dodge',
                                    alpha = 0.8) +
    scale_fill_manual(name = '', values = cols) +
    labs(x = orig_labs[1], y = orig_labs[3]) +
    ggtitle(title) +
    geom_text(aes(label=mean_value), position=position_dodge(width=0.9), 
              vjust=-0.25, size = 5) +
    theme_minimal(base_size = 14, base_family = 'Ubuntu')
    
  return(p)

}

# get a temporary clinical data set used in model
temp_mod_dat <- full_data_first[, 1:17]
temp_mod_con <- temp_mod_dat[grepl('Unaffected', temp_mod_dat$cancer_diagnosis_diagnoses),]
temp_con_last <- temp_mod_con[!duplicated(temp_mod_con$tm_donor_, fromLast = T),]
temp_con_first <- temp_mod_con[!duplicated(temp_mod_con$tm_donor_, fromLast = F),]


# recode for just a cancer classifier in clin and modle data
clin$cancer_fac <- ifelse(grepl('Unaffected', clin$cancer_diagnosis_diagnoses), 
                          'No', 
                          ifelse(is.na(clin$cancer_diagnosis_diagnoses),
                                 NA, 'Yes'))

# recode for just a cancer classifier in temp_mod_dat and modle data
temp_mod_dat$cancer_fac <- ifelse(grepl('Unaffected', temp_mod_dat$cancer_diagnosis_diagnoses), 
                          'No', 
                          ifelse(is.na(temp_mod_dat$cancer_diagnosis_diagnoses),
                                 NA, 'Yes'))


# p53 and gender 
get_bar_plot(clin, var1 = 'gender', var2 = 'p53_germline', var3 = 'age_diagnosis', title = 'All clinical data')
get_bar_plot(temp_mod_dat, var1 = 'gender', var2 = 'p53_germline', var3 = 'age_diagnosis', title = 'Data used in model')

# gender and cases controls  
get_bar_plot(clin, var1 = 'gender', var2 = 'cancer_fac', var3 = 'age_diagnosis', title = 'All clinical data')
get_bar_plot(temp_mod_dat, var1 = 'gender', var2 = 'cancer_fac', var3 = 'age_diagnosis', title = 'All clinical data')

# p53 and gender 
get_bar_plot(clin, var1 = 'gender', var2 = 'p53_germline', var3 = 'age_sample_collection', title = 'All clinical data')
get_bar_plot(temp_mod_dat, var1 = 'gender', var2 = 'p53_germline', var3 = 'age_sample_collection', title = 'Data used in model')

# gender and cases controls  
get_bar_plot(clin, var1 = 'gender', var2 = 'cancer_fac', var3 = 'age_sample_collection', title = 'All clinical data')
get_bar_plot(temp_mod_dat, var1 = 'gender', var2 = 'cancer_fac', var3 = 'age_sample_collection', title = 'All clinical data')


# histogram for onset, and onset by cases and controls
# Histogram overlaid with kernel density curve
# ggplot(clin, aes(x=age_sample_collection)) + 
#   geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
#                  binwidth=.5,
#                  colour="black", fill="white") +
#   geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot

length(which(is.na(clin$cancer_diagnosis_diagnoses)))
length(which(is.na(clin$age_sample_collection)))
length(which(is.na(clin$age_diagnosis)))


temp_mod_dat <- temp_mod_dat[!duplicated(temp_mod_dat$tm_donor_),]

g1 <- ggplot(clin, aes(x=age_diagnosis)) + 
  geom_histogram(aes(y=..count..),      # Histogram with density instead of count on y-axis
                 binwidth=5,
                 colour="black", fill="grey") +
  geom_density(aes(y = ..density.. *(427*5)),
               alpha=.2, fill="blue") +
  labs(title = '', x = 'Age of diagnosis', y = 'Counts') +
  theme_minimal(base_size = 12, base_family = 'Ubuntu')

g2 <- ggplot(clin, aes(x=age_sample_collection)) + 
  geom_histogram(aes(y=..count..),      # Histogram with density instead of count on y-axis
                 binwidth=5,
                 colour="black", fill="grey") +
  geom_density(aes(y = ..density.. *(430*5)),
               alpha=.2, fill="blue") +
  labs(title = '', x = 'Age of sample collection', y = 'Counts') +
  theme_minimal(base_size = 12, base_family = 'Ubuntu')

library(gridExtra)

grid.arrange(g1, g2, nrow= 2)
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



temp_p53 <- clin %>%
  filter(!is.na(p53_germline)) %>%
  group_by(p53_germline) %>%
  summarise(counts = n())

pie(temp_p53$counts,temp_p53$p53_germline , main="",border = 'black', radius = 1,
    col = adjustcolor(c('grey', 'blue'), alpha.f =0.6), cex = 1.2)

temp_gen <- clin %>%
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  summarise(counts = n())

pie(temp_gen$counts,temp_gen$gender , main="",border = 'black',clockwise = TRUE,
    col = adjustcolor(c('darkgrey', 'red'), alpha.f =0.6), cex = 1.2)

##########
# generate random heatmap
##########

nba_heatmap <- heatmap(nba_matrix, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))