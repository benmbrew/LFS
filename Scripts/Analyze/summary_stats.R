########################################################
# This script will generate summary statistics on the clinical data.
# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/classification_template')
data_folder <- paste0(project_folder, '/Data')
clin_data <- paste0(data_folder, '/clin_data')

library(reshape2)
source('/home/benbrew/Desktop/budget/Lib/helpers.R')

########################### Load in clinical data
clin <- read.csv(paste0(data_folder, '/clin.csv'))

########################### give data the right structure

columnStructure <- function(data, col_names) {

  for (i in 1:length(col_names)) {
    
    data_column <- col_names[i]
    
    if (grepl('age', data_column)) {
      data[,data_column] <- as.numeric(as.character(data[,data_column]))
    } else {
      data[,data_column] <- as.factor(data[,data_column])
    }
    
  }
  return(data)
}

clin <- columnStructure(clin, col_names = names(clin))
##########################
# create age difference variable for time between diagnosis and sample
clin$age_diff <- clin$age_sample_collection - clin$age_diagnosis
hist(clin$age_diff, 
     xlab = 'Difference in Years',
     main = 'Distribution of Time Between 
     Sample Collection and Diagnosis',
     col = 'lightblue',
     breaks = 20)

plot(clin$age_sample_collection, 
     clin$age_diagnosis,
     col = adjustcolor('black', alpha.f = 0.6),
     pch = 16, 
     ylab = 'Age of Diagnosis',
     xlab = 'Age of Sample Collection',
     bty = 'n')
abline(0,1)



###########################
# Groupy by variables and get counts.
cancer <- clin %>%
  group_by(cancer_diagnosis_diagnoses, p53_germline) %>%
  summarise(counts = n(),
            mean_age = mean(age_diagnosis, na.rm = T))


###########################
# histogram on age variables 
hist(clin$age_diagnosis, 
     xlab = 'Age of Diagnosis',
     main = '',
     col = 'lightblue',
     breaks = 20)

hist(clin$age_sample_collection, 
     xlab = 'Age Sample of Collection',
     main = '',
     col = 'lightblue',
     breaks = 20)

hist(clin$difference, 
     xlab = 'Difference in Years',
     main = 'Distribution of Time Between 
     Sample Collection and Diagnosis',
     col = 'lightblue',
     breaks = 20)

###########################
# barplot of factor variables
# cancer or not and with mut and WT
clin$cancer <- ifelse(clin$cancer_diagnosis_diagnoses != 'Unaffected', 'Cancer', 'No_Cancer')

cancer_p53 <- clin %>%
  group_by(cancer) %>%
  summarise(Mut = sum(p53_germline == 'Mut', na.rm = T),
            WT = sum(p53_germline == 'WT', na.rm = T))


cancer_p53 <- melt(cancer_p53, id = 'cancer')

ggplot(data = cancer_p53, aes(cancer, value)) + 
  geom_bar(aes(fill = variable), position = 'dodge', stat = 'identity') +
  ylab('Count') + 
  theme(panel.background=element_rect(fill="#F0F0F0"), 
        plot.background=element_rect(fill="#F0F0F0"), 
        # panel.border=element_rect(colour="#F0F0F0"),
        panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
        legend.position="right",  plot.title=element_text(face="bold",hjust=-.08,vjust=2,size=20),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) 

