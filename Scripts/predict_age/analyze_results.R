####### This script will load tables and models and analyze results
# this is 8th step in pipeline

##########
# initialize libraries
##########
library(dplyr)
library(ggplot2)
library(ggthemes)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
results_folder <- paste0(project_folder, '/Results')
scripts_folder <- paste0(project_folder, '/Scripts')


# source model_functions to get functions to run models 
source(paste0(scripts_folder, '/predict_age/model_functions.R'))

##########
# load table rda files
##########
# load table results 
beta_raw <- readRDS(paste0(results_folder, '/beta_raw_model_results.rda'))
beta_quan <- readRDS(paste0(results_folder, '/beta_quan_model_results.rda'))
beta_swan <- readRDS(paste0(results_folder, '/beta_swan_model_results.rda'))
beta_funnorm <- readRDS(paste0(results_folder, '/beta_funnorm_model_results.rda'))
beta_rand <- readRDS(paste0(results_folder, '/beta_rand_model_results.rda'))

# load table results for random

# 100
beta_raw_rand_100 <- readRDS <- readRDS(paste0(results_folder, '/beta_raw_rand_models_100.rda'))
beta_swan_rand_100 <- readRDS(paste0(results_folder, '/beta_swan_rand_models_100.rda'))
beta_quan_rand_100 <- readRDS(paste0(results_folder, '/beta_quan_rand_models_100.rda'))
beta_funnorm_rand_100 <- readRDS(paste0(results_folder, '/beta_funnorm_rand_models_100.rda'))

# 500
beta_raw_rand_500 <- readRDS <- readRDS(paste0(results_folder, '/beta_raw_rand_models_500.rda'))
beta_swan_rand_500 <- readRDS(paste0(results_folder, '/beta_swan_rand_models_500.rda'))
beta_quan_rand_500 <- readRDS(paste0(results_folder, '/beta_quan_rand_models_500.rda'))
beta_funnorm_rand_500 <- readRDS(paste0(results_folder, '/beta_funnorm_rand_models_500.rda'))

# 1000
beta_raw_rand_1000 <- readRDS <- readRDS(paste0(results_folder, '/beta_raw_rand_models_1000.rda'))
beta_swan_rand_1000 <- readRDS(paste0(results_folder, '/beta_swan_rand_models_1000.rda'))
beta_quan_rand_1000 <- readRDS(paste0(results_folder, '/beta_quan_rand_models_1000.rda'))
beta_funnorm_rand_1000 <- readRDS(paste0(results_folder, '/beta_funnorm_rand_models_1000.rda'))

# 2000
beta_raw_rand_2000 <- readRDS <- readRDS(paste0(results_folder, '/beta_raw_rand_models_2000.rda'))
beta_swan_rand_2000 <- readRDS(paste0(results_folder, '/beta_swan_rand_models_2000.rda'))
beta_quan_rand_2000 <- readRDS(paste0(results_folder, '/beta_quan_rand_models_2000.rda'))
beta_funnorm_rand_2000 <- readRDS(paste0(results_folder, '/beta_funnorm_rand_models_2000.rda'))

# 10000
beta_raw_rand_10000 <- readRDS <- readRDS(paste0(results_folder, '/beta_raw_rand_models_10000.rda'))
beta_swan_rand_10000 <- readRDS(paste0(results_folder, '/beta_swan_rand_models_10000.rda'))
beta_quan_rand_10000 <- readRDS(paste0(results_folder, '/beta_quan_rand_models_10000.rda'))
beta_funnorm_rand_10000 <- readRDS(paste0(results_folder, '/beta_funnorm_rand_models_10000.rda'))

##########
# combine results
##########
result_table <- rbind(beta_raw, 
                      beta_quan, 
                      beta_swan, 
                      beta_funnorm,
                      beta_raw_rand, 
                      beta_quan_rand, 
                      beta_swan_rand, 
                      beta_funnorm_rand)

# remove objects
rm(beta_raw,
   beta_quan,
   beta_swan,
   beta_funnorm,
   beta_raw_rand, 
   beta_quan_rand, 
   beta_swan_rand, 
   beta_funnorm_rand)

##########
# find the ones that score high on regression -> then link to categorical.
##########
result_mut <- result_table %>% 
  filter(p53_status == 'Mut') %>% 
  group_by(data, type, age) %>% 
  summarise(mean_score = mean(score))

result_mut <- result_mut[order(result_mut$mean_score, decreasing = T),]


##########
# load in model data summaries
##########
model_data_stats <- readRDS(paste0(results_folder, '/model_data_stats.rda'))

# 


##########
# Get plot for beta_funnorm_cancer_bal_models
##########
plot_object <- plotObject(beta_funnorm_cancer_bal_models, 
                          residual = T, 
                          p53_mut = T)
plotModel(plot_object, 
          main1 = 'Age of Onset',
          main2 = 'Age of Sample Collection',
          xlim = c(0,1000),
          ylim = c(0,1000))


##########
# use ggplot for age of onset
##########

# 4 and 6 are predictions and real age 
pred <- unlist(plot_object[[4]])
real <- unlist(plot_object[[6]])
plotData <- as.data.frame(cbind(predictions = pred, real = real))

ggplot(data = plotData, aes(pred, real)) + 
  geom_point(colour = 'black', alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.3) +
  xlim(0, 800) +
  ylim(0, 800) +
  xlab('Model predictions') +
  ylab('Real age of onset (months)') +
  ggtitle('Age of onset') + theme(panel.background=element_rect(fill="#F0F0F0"), 
                                  plot.background=element_rect(fill="#F0F0F0"), 
                                  # panel.border=element_rect(colour="#F0F0F0"),
                                  panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
                                  legend.position="right",  plot.title=element_text(face="bold",hjust=-.08,vjust=2,colour="#3C3C3C",size=20),
                                  axis.text.x=element_text(size=11,colour="#535353",face="bold"),
                                  axis.text.y=element_text(size=11,colour="#535353",face="bold"),
                                  axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
                                  axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5))
  
##########
# use ggplot for age of sample collection
##########
# 4 and 6 are predictions and real age 
pred <- unlist(plot_object[[4]])
real_samp <- unlist(plot_object[[8]])
plotData <- as.data.frame(cbind(predictions = pred, real_samp = real_samp))

ggplot(data = plotData, aes(pred, real_samp)) + 
  geom_point(colour = 'black', alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.3) +
  xlim(0, 800) +
  ylim(0, 800) +
  xlab('Model predictions') +
  ylab('Real age of sample collection (months)') +
  ggtitle('Age of sample collection') + theme(panel.background=element_rect(fill="#F0F0F0"), 
                                  plot.background=element_rect(fill="#F0F0F0"), 
                                  # panel.border=element_rect(colour="#F0F0F0"),
                                  panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
                                  legend.position="right",  plot.title=element_text(face="bold",hjust=-.08,vjust=2,colour="#3C3C3C",size=20),
                                  axis.text.x=element_text(size=11,colour="#535353",face="bold"),
                                  axis.text.y=element_text(size=11,colour="#535353",face="bold"),
                                  axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
                                  axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5))

##########
# Get confusion matrix with beta_funnorm_cancer_bal - #74 obs, 989 features 
##########
# train.predictions, test.predictions, train.ground_truth, test.ground_truth, train.sample_collection,
# test.sample_collection, test_acc, test_stats, test_acc_samp, test_stats_samp, model, importance, dims)
# 22 people each time, 49 sample, 25 and 24

# get real value frequency 22 in test set
real_freq <- unlist(mat_object[[4]])

# loop through 10 times (number of iterations) and get avg # of times mat_object == X1
for (i in 1:iterations) {
  
}

mat_object <- matObject(beta_funnorm_cancer_bal_models, age = 60, residual = T, p53_mut = T)
conMatrix(mat_object)

conMatrix <- function(results) {
  
  # test acc for age of diagnosis
  acc_age <- mean(unlist(results[[7]]))
  
  # test acc for age of sample collection
  acc_samp <-mean(unlist(results[[9]]))
  
  # confustion matrix age of diagnosis 10
  iterations <- 10
  temp <- list()
  for (i in 1:10){
    temp[[i]] <- results[[8]][[i]]$table
  }
  mat <- unlist(temp)
  new_mat <- matrix(, 2, 2)
  
  mat_index <- seq(1, length(mat), 4)
  
  new_mat[1,1] <- sum(mat[mat_index])/iterations
  new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
  new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
  new_mat[2,2] <- sum(mat[mat_index + 3])/iterations
  
  return(list(new_mat, acc_age, acc_samp))
  
}


##########
# Get top probes - beta_funnorm_cancer_bal_models
##########
# classification - 4(48,60, 72, 84), 2 (mut,wt), 13 (12), 10
# regression - 2(mut, wt), 15 (14), 10, 1

# get
# temp <- beta_funnorm_cancer_bal_models[[3]]
# length(temp)
# temp1 <- temp[[1]]
# length(temp1)
# temp2 <- temp1[[1]]
# length(temp2)
# temp3 <- temp2[[1]]
# length(temp3)
# length(beta_funnorm_cancer_bal_models)
# 1 data_result, 2 data_resid_result, 3 data_fac_result, 4 data_resid_fac_result
# 1 mut, 2 WT
# 12 for fac, 14 for reg
# function to grab mutant, 48 months. 
topProbeObject <- function(results, 
                           normal, 
                           fac) {
  
  if (normal) {
    normal_reg <- results[[1]]
    normal_fac <- results[[3]]
    
    if (fac) {
      data_48 <- normal_fac[[1]]
      data_48_mut <- data_48[[1]]
      data_48_mut_features <- data_48_mut[[12]]
    } else {
      data_48_mut <- nomal_reg[[1]]
      data_48_mut_features<- data_48_mut[[14]]
    }
    
  } else {
    resid_reg <- results[[2]]
    resid_fac <- results[[4]]
    
    if (fac) {
      # slot 2 is for 68
      data_48 <- resid_fac[[2]]
      data_48_mut <- data_48[[1]]
      data_48_mut_features <- data_48_mut[[12]]
    } else {
      data_48_mut <- resid_reg[[1]]
      data_48_mut_features <- data_48_mut[[14]]
    }
  }
  
  return(data_48_mut_features)
  
}

###########
# get for residual fac and residual normal
###########

###########
# apply to our best model
###########

# classification
beta_funnorm_cancer_bal_feature_importance_fac <- topProbeObject(beta_funnorm_cancer_bal_models,
                                                             normal = F,
                                                             fac = T)
# regression
beta_funnorm_cancer_bal_feature_importance_reg <- topProbeObject(beta_funnorm_cancer_bal_models,
                                                                 normal = F,
                                                                 fac = F)

# function to get top 50 non duplicated importance probes
getImportance <- function(importance_data) {
  
  importance_data <- as.data.frame(do.call(rbind, importance_data))
  importance_data$V2 <- as.numeric(as.character(importance_data$V2))
  importance_data <- importance_data[order(importance_data$V2, decreasing = T),]
  importance_data <- importance_data[1:50,]
  importance_data <- importance_data %>% group_by(V1) %>% summarise_each(funs(mean))
  importance_data <- importance_data[order(importance_data$V2, decreasing = T),]
  return(importance_data)
}

fac_importance <- getImportance(beta_funnorm_cancer_bal_feature_importance_fac)
names(fac_importance) <- c('probe', 'accuracy')
reg_importance <- getImportance(beta_funnorm_cancer_bal_feature_importance_reg)
names(reg_importance) <- c('probe', 'correlation')



##########
# compare both regression and classification
##########

# find how many are in both
length(which(fac_importance$probe %in% reg_importance$probe))
length(which(reg_importance$probe %in% fac_importance$probe))

# only 12 in both

# get probes that are in both
both_importance <- inner_join(fac_importance, reg_importance, by = 'probe')

# # save all 3 data sets temporarily as csv in model_Data
# write.csv(fac_importance, paste0(model_data, '/fac_importance.csv'))
# write.csv(reg_importance, paste0(model_data, '/reg_importance.csv'))
# write.csv(both_importance, paste0(model_data, '/both_importance.csv'))

################################################################################################################################
# load in regions data and find regions, genes that match the top spaces

##########
# get probes for bumphunter results
##########
library(minfi)
setwd('/home/benbrew/')

getRgSet <- function() {
  
  # #idat files
  # idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
  # sapply(idatFiles, gunzip, overwrite = TRUE)
  # 
  # read into rgSet
  rgSet <- read.450k.exp("GSE68777/idat")
  
  # preprocess quantil
  rgSet <- preprocessQuantile(rgSet)
  
  # get rangers 
  rgSet <- granges(rgSet)
  rgSet <- as.data.frame(rgSet)
  return(rgSet)
}

rgSet <- getRgSet()

##########
# combine rgSet with importance probes
##########
# create probe variables to match with importance data by probe
rgSet$probe <- rownames(rgSet)

# function to join and order
joinData <- function(data) {
  data <- inner_join(rgSet, data, by = 'probe')
  data <- data[order(data[, 7], decreasing = T),]
  return(data)
}

reg_importance <- joinData(reg_importance)
fac_importance <- joinData(fac_importance)
both_importance <- joinData(both_importance)
