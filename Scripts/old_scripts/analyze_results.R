#######################################################################
# This scirpt will evaluate the results from the 3 basic test.

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
code_folder <- paste0(project_folder, '/Code')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


# Load Libraries 
library(dplyr)
source(paste0(code_folder, '/Lib/helpers.R'))



# read in data 
clin <- read.csv(paste0(data_folder, '/clin_model_results.csv'))
clin_train <- read.csv(paste0(data_folder, '/clin_train_model_results.csv'))
clin_methyl <- read.csv(paste0(data_folder, '/clin_methyl_model_results.csv'))

# remove unnecessarry columns and add columns to indicate data set and data
clin$X <- NULL
clin_train$X <- NULL
clin_methyl$X <- NULL

# add + signs in front of each cell in variable 
clin$variables <- paste0('+', clin$variables)
clin_train$variables <- paste0('+', clin_train$variables)
clin_methyl$variables <- paste0('+', clin_methyl$variables)

# change level of factors to match number of observations
clin$variables <- factor(clin$variables, levels = c("+gender_and_gdna.base.change", "+gdna.codon",
                                                    "+protein.codon.change", "+gdna.exon.intron", "+codon72.npro",
                                                    "+splice.delins.snv", "+protein.codon.num","+mdm2.nG"))
clin_train$variables <- factor(clin_train$variables, levels = c("+gender_and_gdna.base.change", "+gdna.codon",
                                                    "+protein.codon.change", "+gdna.exon.intron", "+codon72.npro",
                                                    "+splice.delins.snv", "+protein.codon.num","+mdm2.nG"))
clin_methyl$variables <- factor(clin_methyl$variables, levels = c("+just_methyl", "+gender_and_gdna.base.change", "+gdna.codon",
                                                    "+protein.codon.change", "+gdna.exon.intron", "+codon72.npro",
                                                    "+splice.delins.snv", "+protein.codon.num","+mdm2.nG"))

# First look at just clin. that is train and tested on clin 

# group by model and get mean rmse 
model <- clin %>%
  group_by(model) %>%
  summarise(rmse = mean(rmse))

# subset by random forest and plot
rf_clin <- clin[clin$model == 'rand_forest',]

ggplot(rf_clin, aes(variables, rmse)) + geom_bar(stat = 'identity') + 
  xlab('Additional Variables') + ylab('RMSE') + 
  geom_text(aes(variables, rmse, ymax=rmse, label= paste0('n =', observations), vjust = 0)) +
  theme(text = element_text(size=8),
  axis.text.x = element_text(angle=30)) 



# look at just clin_train. that is train on clin w/out methyl, test on clin with methyl

# group by model and get mean rmse 
model <- clin_train %>%
  group_by(model) %>%
  summarise(rmse = mean(rmse))

# subset by random forest and plot
rf_clin_train <- clin_train[clin_train$model == 'rand_forest',]

ggplot(rf_clin_train, aes(variables, rmse)) + geom_bar(stat = 'identity') + 
  xlab('Additional Variables') + ylab('RMSE') + 
  geom_text(aes(variables, rmse, ymax=rmse, label= paste0('n =', observations), vjust = 0)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=30)) 


# look at just clin_methyl. that is train and tested on methyl with clin variables

# group by model and get mean rmse 
model <- clin_methyl %>%
  group_by(model) %>%
  summarise(rmse = mean(rmse))

# subset by random forest and plot
rf_clin_methyl <- clin_methyl[clin_methyl$model == 'rand_forest',]

ggplot(rf_clin_methyl, aes(variables, rmse, group = data, fill = data)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  xlab('Additional Variables') + ylab('RMSE') + 
  geom_text(aes(variables, rmse, ymax=rmse, label= paste0('n =', observations), vjust = 0)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=30)) 


#####
# make plot comparing the three different RMSE for each test
clin$data <- 'just_clinical'
clin$dataset <- 'just_clinical'
clin_train$data <- 'train_clinical'
clin_train$dataset <- 'train_clinical'
clin_methyl$dataset <- 'methyl'

dat <- rbind(clin, clin_train, clin_methyl)


# group by dataset 
dataset <- dat %>%
  group_by(data) %>%
  summarise(rmse = mean(rmse))

ggplot(dataset, aes(data, rmse)) + 
  geom_bar(stat = 'identity') + 
  xlab('Data') + ylab('RMSE') 
