####################################################################
# This script will analyze the results table from models on idat data.
library(dplyr)
library(ggplot2)
library(ggthemes)
# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
scripts_folder <- paste0(project_folder, '/Scripts')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
idat_data <- paste0(methyl_data, '/raw_files')
clin_data <- paste0(data_folder, '/clin_data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')
model_data <- paste0(data_folder, '/model_data')


# source model_functions to get functions to run models 
source(paste0(scripts_folder, '/predict_age/model_functions.R'))

# load in results table 
load(paste0(model_data, '/idat_beta_table_results.RData'))


rm(beta_funnorm_cancer_bal_table, beta_funnorm_cancer_unbal_table, beta_funnorm_global_bal_table, beta_funnorm_global_unbal_table,
   beta_raw_cancer_bal_table, beta_raw_cancer_unbal_table, beta_raw_global_bal_table, beta_raw_global_unbal_table,
   beta_quan_cancer_bal_table, beta_quan_cancer_unbal_table, beta_quan_global_bal_table, beta_quan_global_unbal_table,
   beta_swan_cancer_bal_table, beta_swan_cancer_unbal_table, beta_swan_global_bal_table, beta_swan_global_unbal_table,
   beta_raw_rand_table, beta_quan_rand_table, beta_swan_rand_table, beta_funnorm_rand_table, data_thresholds)

# order results 
beta_idat_results <- beta_idat_results[order(beta_idat_results$score, decreasing = T),]
# beta_funnorm_global_bal
# beta_funnorm_cancer_unbal 
# beta_funnorm_cancer_bal 89 correlation, 84 accuracy

#############
# find the ones that score high on regression -> then link to categorical.
#############
temp <- beta_idat_results %>% 
  filter(p53_status == 'Mut') %>% 
  group_by(data, type, age) %>% 
  summarise(mean_score = mean(score))

temp <- temp[order(temp$mean_score, decreasing = T),]

################
# Get plot for beta_funnorm_cancer_bal_models
################
plot_object <- plotObject(beta_funnorm_cancer_bal_models, 
                          residual = T, 
                          p53_mut = T)
plotModel(plot_object, 
          main1 = 'Age of Onset',
          main2 = 'Age of Sample Collection',
          xlim = c(0,1000),
          ylim = c(0,1000))


#####################
# use ggplot for age of onset
#####################
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
  
#####################
# use ggplot for age of sample collection
#####################
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

################
# Get confusion matrix with beta_funnorm_cancer_bal - #74 obs, 989 features 
# ################
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




