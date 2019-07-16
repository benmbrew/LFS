####### This script will load tables and models and analyze results
# this is 8th step in pipeline

##########
# initialize libraries
##########
library(dplyr)
library(ggplot2)
library(ggthemes)
library(plotly)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
results_folder <- paste0(project_folder, '/Results')
raw_folder <- paste0(results_folder, '/raw')
model_data <- paste0(data_folder, '/model_data')
scripts_folder <- paste0(project_folder, '/Scripts')


# source model_functions to get functions to run models 
source(paste0(scripts_folder, '/predict_age/functions.R'))

##########
# load table rda files
##########
result_table <- readRDS(paste0(model_data, '/result_table.rds'))
full_table <- readRDS(paste0(model_data, '/even_full.rds'))
rand_table <- readRDS(paste0(model_data, '/rand.rda'))


result_table <- rbind(result_table,
                      full_table)

#########
# random - group by data and features (ungen, gen, sam, etc)
#########

# first make features a factor to group by
rand_table$features <- as.character(rand_table$features)

# group by
rand_group <- rand_table %>%
  group_by(data, features) %>%
  summarise(mean_score = mean(score),
            sd_score = sd(score))

# order by features 
rand_group <- rand_group[order(rand_group$features),]

# ci

error <- qnorm(0.975)*rand_group$sd_score/sqrt(10)
rand_group$upper <- rand_group$mean_score + error
rand_group$lower <- rand_group$mean_score - error


##########
# not random - group by data and features (ungen, gen, sam, etc)
##########

# order by score
result_table <- result_table[order(result_table$score, decreasing = T),]

