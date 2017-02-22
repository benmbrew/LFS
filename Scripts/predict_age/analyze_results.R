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
swan_folder <- paste0(results_folder, '/swan')
quan_folder <- paste0(results_folder, '/quan')
funnorm_folder <- paste0(results_folder, '/funnorm')
rand_folder <- paste0(results_folder, '/rand')
scripts_folder <- paste0(project_folder, '/Scripts')


# source model_functions to get functions to run models 
source(paste0(scripts_folder, '/predict_age/model_functions.R'))

##########
# load table rda files
##########

# 4 methods
raw_table <- readRDS(paste0(raw_folder, '/raw_table.rda'))
swan_table <- readRDS(paste0(swan_folder, '/swan_table.rda'))
quan_table <- readRDS(paste0(quan_folder, '/quan_table.rda'))
funnorm_table <- readRDS(paste0(funnorm_folder, '/funnorm_table.rda'))

# random
raw_rand_table <- readRDS(paste0(rand_folder, '/raw_rand_table.rda'))
swan_rand_table <- readRDS(paste0(rand_folder, '/swan_rand_table.rda'))
quan_rand_table <- readRDS(paste0(rand_folder, '/quan_rand_table.rda'))
funnorm_rand_table <- readRDS(paste0(rand_folder, '/funnorm_rand_table.rda'))

##########
# combine results
##########

# 4 methods 
result_table <- rbind(raw_table, 
                      swan_table,
                      quan_table,
                      funnorm_table)

# rand
result_table_rand <- rbind(raw_rand_table,
                           swan_rand_table,
                           quan_rand_table,
                           funnorm_rand_table)
   
##########                        
# remove objects
##########

# 4 methods
rm(raw_table, 
   swan_table,
   quan_table,
   funnorm_table)

# rand
rm(raw_rand_table,
   swan_rand_table,
   quan_rand_table,
   funnorm_rand_table)

#####################################################################################################################

# ##########
# # recode features variable to have number of features (strsplit)
# ##########
# result_table$features_num <- gsub(')', '', unlist(lapply(strsplit(as.character(result_table$features), 
#                                      ' |:'), function(x) x[2])), fixed = T)
# 
# result_table_rand$features_num <- gsub(')', '', unlist(lapply(strsplit(as.character(result_table_rand$data), 
#                                                                   '_'), function(x) x[3])), fixed = T)

##########
# group by age, type, data, features to get rand like normal
##########
result_table_rand<- result_table_rand %>%
  group_by(age, type, data, features) %>%
  summarise(score = mean(score))

##########
# recode data type to grab normalization method
##########
result_table$norm <- unlist(lapply(strsplit(as.character(result_table$data), '_'), function(x) x[1]))
result_table_rand$norm <- unlist(lapply(strsplit(as.character(result_table_rand$data), '_'), function(x) x[1]))


##########
# order the score column and get top performers
##########
result_table <- result_table[order(result_table$score, decreasing = T),]
result_table_rand <- result_table_rand[order(result_table_rand$score, decreasing = T),]


##########
# seperate by regression or classification
##########
# normal
result_table <- result_table[result_table$age == 'regression',]

# random
result_table_rand <- result_table_rand[result_table_rand$age == 'regression',]

##########
# combine random and nromal
##########

# add random or not column to both data sets
result_table$feature_type <- 'bh'
result_table_rand$feature_type <- 'random'

# change into data frame
result_table <- as.data.frame(result_table)
result_table_rand <- as.data.frame(result_table_rand)

# get feature number
result_table$feature_num <- as.numeric(unlist(lapply(strsplit(as.character(result_table$features), '_'), 
                                                     function(x) x[2])))
result_table_rand$feature_num<- as.numeric(unlist(lapply(strsplit(as.character(result_table_rand$features), '_'), 
                                              function(x) x[2])))


# combine reg
result <- rbind(result_table, result_table_rand)

# drop union features 
result <- result[!grepl('union|07|08|09|10', result$data),]

##########
# histograms for regression
##########

# histrogram of regression scores for both bh and random
ggplot(result, aes(score)) + 
  geom_histogram(data = subset(result, type == 'normal' & feature_type == 'bh'), fill = 'green', alpha = 0.8) +
  geom_histogram(data = subset(result, type == 'normal' & feature_type == 'random'), fill = 'black', alpha = 0.8)

# histrogram of residual, regression scores for both bh and random
ggplot(result, aes(score)) + 
  geom_histogram(data = subset(result, type == 'resid' & feature_type == 'bh'), fill = 'green', alpha = 0.8) +
  geom_histogram(data = subset(result, type == 'resid' & feature_type == 'random'), fill = 'black', alpha = 0.8)

# historgram by preprocessing mehtod 
ggplot(result, aes(score)) + 
  geom_histogram(data = subset(result, type == 'normal' & norm == 'raw'), fill = 'green', alpha = 0.8) +
  geom_histogram(data = subset(result, type == 'normal' & norm == 'quan'), fill = 'black', alpha = 0.8) +
  geom_histogram(data = subset(result, type == 'normal' & norm == 'swan'), fill = 'red', alpha = 0.8) +
  geom_histogram(data = subset(result, type == 'normal' & norm == 'funnorm'), fill = 'lightblue', alpha = 0.8)
  
# histrogram of funnorm, residual, regression scores for both bh and random
ggplot(result, aes(score)) + 
  geom_histogram(data = subset(result, norm == 'swan' & type == 'normal' & feature_type == 'bh'), fill = 'green', alpha = 0.8) +
  geom_histogram(data = subset(result, norm == 'swan' & type == 'normal' & feature_type == 'random'), fill = 'black', alpha = 0.8)

# subset result so that its normal and not resid
result <- result[result$type == 'normal',]
result <- result[result$norm == 'swan',]

##########
# recode feature_num to cateogries 
##########
hist(result$feature_num, breaks= 20)
summary(result$feature_num)

# create splits at 500 - 1000, 1001 - 3000, 
result$feature_cat <- ifelse(result$feature_num > 0 & result$feature_num <= 1600, '0_1600',
                             ifelse(result$feature_num > 1600 & result$feature_num <= 4000, '1600_4k',
                                    ifelse(result$feature_num > 4000 & result$feature_num <= 7000, '4k-7k',
                                           ifelse(result$feature_num > 7000 & result$feature_num <= 15000, '7k-15k',
                                                  ifelse(result$feature_num > 15000 & result$feature_num <= 35000,'15k-35k', '30k-80k')))))

result_group <- result %>%
  group_by(feature_cat,feature_type) %>%
  summarise(mean_score = max(score))


##########
# barplot with feature_cat as x and score as y
##########
ggplot(result_group, aes(feature_cat, mean_score, group = feature_type, fill = feature_type)) + 
  geom_bar(stat = 'identity', position = 'dodge') 

+ ylim(0.5, 1) 



##########
# read in a random model and see if top features are in union_features
##########


##########
# barplot of for bh and random, where the x axis is ordered by number of features and y axis is score.
##########
bh <- result[result$feature_type == 'bh',]
rand <- result[result$feature_type == 'random',]

plot_ly(x = ~bh$feature_num, y = ~bh$score, type = 'scatter', mode = 'lines', name = 'BumpHunter Feat', fill = 'tozeroy') %>%
  add_trace(x = ~rand$feature_num, y = ~rand$score, name = 'Random Feat', fill = 'tozeroy') %>%
  layout(xaxis = list(title = 'Features'),
         yaxis = list(title = 'Score'))


# ##########
# # load in model data summaries
# ##########
# model_data_stats <- readRDS(paste0(results_folder, '/model_data_stats.rda'))
# 
# # 
# 

# 
# 
# ##########
# # Get plot for beta_funnorm_cancer_bal_models
# ##########
# plot_object <- plotObject(beta_funnorm_cancer_bal_models, 
#                           residual = T, 
#                           p53_mut = T)
# plotModel(plot_object, 
#           main1 = 'Age of Onset',
#           main2 = 'Age of Sample Collection',
#           xlim = c(0,1000),
#           ylim = c(0,1000))
# 
# 
# ##########
# # use ggplot for age of onset
# ##########
# 
# # 4 and 6 are predictions and real age 
# pred <- unlist(plot_object[[4]])
# real <- unlist(plot_object[[6]])
# plotData <- as.data.frame(cbind(predictions = pred, real = real))
# 
# ggplot(data = plotData, aes(pred, real)) + 
#   geom_point(colour = 'black', alpha = 0.6) +
#   geom_abline(intercept = 0, slope = 1, alpha = 0.3) +
#   xlim(0, 800) +
#   ylim(0, 800) +
#   xlab('Model predictions') +
#   ylab('Real age of onset (months)') +
#   ggtitle('Age of onset') + theme(panel.background=element_rect(fill="#F0F0F0"), 
#                                   plot.background=element_rect(fill="#F0F0F0"), 
#                                   # panel.border=element_rect(colour="#F0F0F0"),
#                                   panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
#                                   legend.position="right",  plot.title=element_text(face="bold",hjust=-.08,vjust=2,colour="#3C3C3C",size=20),
#                                   axis.text.x=element_text(size=11,colour="#535353",face="bold"),
#                                   axis.text.y=element_text(size=11,colour="#535353",face="bold"),
#                                   axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
#                                   axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5))
#   
# ##########
# # use ggplot for age of sample collection
# ##########
# # 4 and 6 are predictions and real age 
# pred <- unlist(plot_object[[4]])
# real_samp <- unlist(plot_object[[8]])
# plotData <- as.data.frame(cbind(predictions = pred, real_samp = real_samp))
# 
# ggplot(data = plotData, aes(pred, real_samp)) + 
#   geom_point(colour = 'black', alpha = 0.6) +
#   geom_abline(intercept = 0, slope = 1, alpha = 0.3) +
#   xlim(0, 800) +
#   ylim(0, 800) +
#   xlab('Model predictions') +
#   ylab('Real age of sample collection (months)') +
#   ggtitle('Age of sample collection') + theme(panel.background=element_rect(fill="#F0F0F0"), 
#                                   plot.background=element_rect(fill="#F0F0F0"), 
#                                   # panel.border=element_rect(colour="#F0F0F0"),
#                                   panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
#                                   legend.position="right",  plot.title=element_text(face="bold",hjust=-.08,vjust=2,colour="#3C3C3C",size=20),
#                                   axis.text.x=element_text(size=11,colour="#535353",face="bold"),
#                                   axis.text.y=element_text(size=11,colour="#535353",face="bold"),
#                                   axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
#                                   axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5))
# 
# ##########
# # Get confusion matrix with beta_funnorm_cancer_bal - #74 obs, 989 features 
# ##########
# # train.predictions, test.predictions, train.ground_truth, test.ground_truth, train.sample_collection,
# # test.sample_collection, test_acc, test_stats, test_acc_samp, test_stats_samp, model, importance, dims)
# # 22 people each time, 49 sample, 25 and 24
# 
# # get real value frequency 22 in test set
# real_freq <- unlist(mat_object[[4]])
# 
# # loop through 10 times (number of iterations) and get avg # of times mat_object == X1
# for (i in 1:iterations) {
#   
# }
# 
# mat_object <- matObject(beta_funnorm_cancer_bal_models, age = 60, residual = T, p53_mut = T)
# conMatrix(mat_object)
# 
# conMatrix <- function(results) {
#   
#   # test acc for age of diagnosis
#   acc_age <- mean(unlist(results[[7]]))
#   
#   # test acc for age of sample collection
#   acc_samp <-mean(unlist(results[[9]]))
#   
#   # confustion matrix age of diagnosis 10
#   iterations <- 10
#   temp <- list()
#   for (i in 1:10){
#     temp[[i]] <- results[[8]][[i]]$table
#   }
#   mat <- unlist(temp)
#   new_mat <- matrix(, 2, 2)
#   
#   mat_index <- seq(1, length(mat), 4)
#   
#   new_mat[1,1] <- sum(mat[mat_index])/iterations
#   new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
#   new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
#   new_mat[2,2] <- sum(mat[mat_index + 3])/iterations
#   
#   return(list(new_mat, acc_age, acc_samp))
#   
# }
# 
# 
# ##########
# # Get top probes - beta_funnorm_cancer_bal_models
# ##########
# # classification - 4(48,60, 72, 84), 2 (mut,wt), 13 (12), 10
# # regression - 2(mut, wt), 15 (14), 10, 1
# 
# # get
# # temp <- beta_funnorm_cancer_bal_models[[3]]
# # length(temp)
# # temp1 <- temp[[1]]
# # length(temp1)
# # temp2 <- temp1[[1]]
# # length(temp2)
# # temp3 <- temp2[[1]]
# # length(temp3)
# # length(beta_funnorm_cancer_bal_models)
# # 1 data_result, 2 data_resid_result, 3 data_fac_result, 4 data_resid_fac_result
# # 1 mut, 2 WT
# # 12 for fac, 14 for reg
# # function to grab mutant, 48 months. 
# topProbeObject <- function(results, 
#                            normal, 
#                            fac) {
#   
#   if (normal) {
#     normal_reg <- results[[1]]
#     normal_fac <- results[[3]]
#     
#     if (fac) {
#       data_48 <- normal_fac[[1]]
#       data_48_mut <- data_48[[1]]
#       data_48_mut_features <- data_48_mut[[12]]
#     } else {
#       data_48_mut <- nomal_reg[[1]]
#       data_48_mut_features<- data_48_mut[[14]]
#     }
#     
#   } else {
#     resid_reg <- results[[2]]
#     resid_fac <- results[[4]]
#     
#     if (fac) {
#       # slot 2 is for 68
#       data_48 <- resid_fac[[2]]
#       data_48_mut <- data_48[[1]]
#       data_48_mut_features <- data_48_mut[[12]]
#     } else {
#       data_48_mut <- resid_reg[[1]]
#       data_48_mut_features <- data_48_mut[[14]]
#     }
#   }
#   
#   return(data_48_mut_features)
#   
# }
# 
# ###########
# # get for residual fac and residual normal
# ###########
# 
# ###########
# # apply to our best model
# ###########
# 
# # classification
# beta_funnorm_cancer_bal_feature_importance_fac <- topProbeObject(beta_funnorm_cancer_bal_models,
#                                                              normal = F,
#                                                              fac = T)
# # regression
# beta_funnorm_cancer_bal_feature_importance_reg <- topProbeObject(beta_funnorm_cancer_bal_models,
#                                                                  normal = F,
#                                                                  fac = F)
# 
# # function to get top 50 non duplicated importance probes
# getImportance <- function(importance_data) {
#   
#   importance_data <- as.data.frame(do.call(rbind, importance_data))
#   importance_data$V2 <- as.numeric(as.character(importance_data$V2))
#   importance_data <- importance_data[order(importance_data$V2, decreasing = T),]
#   importance_data <- importance_data[1:50,]
#   importance_data <- importance_data %>% group_by(V1) %>% summarise_each(funs(mean))
#   importance_data <- importance_data[order(importance_data$V2, decreasing = T),]
#   return(importance_data)
# }
# 
# fac_importance <- getImportance(beta_funnorm_cancer_bal_feature_importance_fac)
# names(fac_importance) <- c('probe', 'accuracy')
# reg_importance <- getImportance(beta_funnorm_cancer_bal_feature_importance_reg)
# names(reg_importance) <- c('probe', 'correlation')
# 
# 
# 
# ##########
# # compare both regression and classification
# ##########
# 
# # find how many are in both
# length(which(fac_importance$probe %in% reg_importance$probe))
# length(which(reg_importance$probe %in% fac_importance$probe))
# 
# # only 12 in both
# 
# # get probes that are in both
# both_importance <- inner_join(fac_importance, reg_importance, by = 'probe')
# 
# # # save all 3 data sets temporarily as csv in model_Data
# # write.csv(fac_importance, paste0(model_data, '/fac_importance.csv'))
# # write.csv(reg_importance, paste0(model_data, '/reg_importance.csv'))
# # write.csv(both_importance, paste0(model_data, '/both_importance.csv'))
# 
# ################################################################################################################################
# # load in regions data and find regions, genes that match the top spaces
# 
# ##########
# # get probes for bumphunter results
# ##########
# library(minfi)
# setwd('/home/benbrew/')
# 
# getRgSet <- function() {
#   
#   # #idat files
#   # idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
#   # sapply(idatFiles, gunzip, overwrite = TRUE)
#   # 
#   # read into rgSet
#   rgSet <- read.450k.exp("GSE68777/idat")
#   
#   # preprocess quantil
#   rgSet <- preprocessQuantile(rgSet)
#   
#   # get rangers 
#   rgSet <- granges(rgSet)
#   rgSet <- as.data.frame(rgSet)
#   return(rgSet)
# }
# 
# rgSet <- getRgSet()
# 
# ##########
# # combine rgSet with importance probes
# ##########
# # create probe variables to match with importance data by probe
# rgSet$probe <- rownames(rgSet)
# 
# # function to join and order
# joinData <- function(data) {
#   data <- inner_join(rgSet, data, by = 'probe')
#   data <- data[order(data[, 7], decreasing = T),]
#   return(data)
# }
# 
# reg_importance <- joinData(reg_importance)
# fac_importance <- joinData(fac_importance)
# both_importance <- joinData(both_importance)
