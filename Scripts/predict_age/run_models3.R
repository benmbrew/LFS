####### Script will run models

##########
# initialize libraries
##########
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(doParallel)
library(nnet)
library(dplyr)
library(ggthemes)

registerDoParallel(1)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')
results_folder <- paste0(project_folder, '/Results')
quan_folder <- paste0(results_folder, '/quan')
rand_folder <- paste0(results_folder, '/rand')
scripts_folder <- paste0(project_folder, '/Scripts')


##########
# load data batch data
##########
# read cases
quan_cases <- readRDS(paste0(model_data, '/quan_cases.rda'))

quan_cases_sam <- readRDS(paste0(model_data, '/quan_cases_sam.rda'))

quan_cases_sam_gen <- readRDS(paste0(model_data, '/quan_cases_sam_gen.rda'))

# read controls
quan_controls <- readRDS(paste0(model_data, '/quan_controls.rda'))

quan_controls_gen <- readRDS(paste0(model_data, '/quan_controls_gen.rda'))

# load features
load(paste0(model_data, '/bh_feat.RData'))

# load gene locations
cg_locations <- read.csv(paste0(model_data, '/cg_locations.csv'))


##########
# change column names of cg_locations for merging
##########
cg_locations$X <- NULL


##########
# source model_functions to get functions to run models 
##########
source(paste0(scripts_folder, '/predict_age/functions.R'))

# ##########
# # set parameters 
# ##########
# data_thresholds <- c(48, 60, 72)
# p53 <- c('Mut', 'WT')


##########
# function to run models - subset, get residuals, get categorical, predict with regression and fac. 
##########
data <- quan_cases
data_controls <- quan_controls
bump_hunter_data <- quan_even_10
 model <- 'rf'
control = T
# random = F
gender = F
cutoff <- .7
# i = 1
# gender <- F
runModels <- function(data,
                      data_controls,
                      model,
                      control,
                      random = F,
                      bump_hunter = F,
                      bump_hunter_data,
                      num_feat = NULL,
                      gender,
                      residual,
                      seed_num) 
{
  
  # get differenct variations of data
  if (random) {
    
    data <- subsetDat(data, 
                      random = T, 
                      num_feat,
                      seed_num)
  } else {
    
    data <- subsetDat(data, 
                      random = F, 
                      gender,
                      num_feat = NULL)
    
    data_controls <- subsetDat(data_controls, 
                               random = F, 
                               gender,
                               num_feat = NULL)
  }
  
  
  #########
  # get regression data and run regression
  #########
  # get data
  if (bump_hunter) {
    
    data <- bhSubset(data, bh_data = bump_hunter_data, gender)
    data_controls <- bhSubset(data_controls, bh_data = bump_hunter_data, gender)
    
    
  }
  
  if (residual) {
    data <- getResidual(data, gender)
  }
  

  
  if (model == 'rf') {
    
    data_result <- rfPredictReg(data, data_controls, cutoff = .7, iterations = 10, control = F)
    # data_resid_result <- rfPredictReg(data_resid, cutoff = .7, iterations = 5)
  }
  
  if (model == 'enet') {
    
    data_result <- enetPredReg(data, data_controls, N_CV_REPEATS = 3, nfolds = 3, cutoff = .7, iterations = 10)
    data_result_con <- enetPredCon(data, data_controls, N_CV_REPEATS = 3, nfolds = 3, cutoff = .7, iterations = 10)
    
    data_fac_48 <- enetPredFac(data, N_CV_REPEATS = 3, nfolds = 3, cutoff = .7, iterations = 10, thresh = 48)
    data_fac_60 <- enetPredFac(data, N_CV_REPEATS = 3, nfolds = 3, cutoff = .7, iterations = 10, thresh = 60)
    data_fac_72 <- enetPredFac(data, N_CV_REPEATS = 3, nfolds = 3, cutoff = .7, iterations = 10, thresh = 72)

  }
  
  # if (model == 'lasso') {
  #   
  #   data_result <- regPred(data, alpha = 1, nfolds = 5, cutoff = .7, iterations = 5)
  #   # data_resid_result <- regPred(data_resid, alpha = 1, nfolds = 5, cutoff = .7, iterations = 5)
  #   
  # }
  # 
  # #########
  # # get classification data and run classifier
  # #########
  # data_fac <- makeFac(data, threshold = 48)
  # data_resid_fac <- makeFac(data_resid, threshold = 48)
  # 
  # data_fac_result  <- rfPredictFac(data_fac, cutoff = .7, iterations = 5)
  # data_resid_fac_result <- rfPredictFac(data_resid_fac, cutoff = .7, iterations = 5)
  
  return (list(data_result, data_fac_48, data_fac_60, data_fac_72, data_result_con)) 
  # data_fac_result, 
  # data_resid_fac_result))
  
  
}

########################################################################

# 
# # full data sam
# even_full_sam_rf <- runModels(quan_cases_sam, 
#                               model = 'rf',
#                               rand = F,
#                               bump_hunter = F,
#                               bump_hunter_data)
# 
# even_full_sam_table_rf <- extractResults(even_full_sam_rf, data_name = 'even sam_full rf', regularize = F)
# 
# even_full_sam_enet <- runModels(quan_cases_sam, 
#                                 model = 'enet',
#                                 rand = F,
#                                 bump_hunter = F,
#                                 bump_hunter_data)
# 
# even_full_sam_table_enet <- extractResults(even_full_sam_enet, data_name = 'even sam_full enet', regularize = F)
# 
# #here
# # get table 
# even_full_sam <- rbind(even_full_sam_table_rf,
#                        even_full_sam_table_enet)
# 
# 
# rm(even_full_sam_table_rf,
#    even_full_sam_table_enet)
# 
# 
# 
# saveRDS(even_full_sam, paste0(model_data, '/even_full_2.rds'))

########################################################################3
##########
# run models
##########


#######
# WITH gender as covariate


#### 10
mod <- runModels(quan_cases_sam, 
                 quan_controls,
                 model = 'enet',
                 rand = F,
                 bump_hunter = T,
                 bump_hunter_data = quan_even_fwer_10,
                 gender= T,
                 residual = F)


# get reg, and class
mod_reg <- mod[[1]]
mod_48 <- mod[[2]]
mod_60 <- mod[[3]]
mod_72 <- mod[[4]]

# get table
mod_reg_table <- extractResults(mod_reg, 
                                data_name = 'sam_fwer_gender_10_enet', 
                                regularize = T,
                                bh_data = quan_even_fwer_10)

mod_48
# get conmatrix
conMatrix(mod_48)
conMatrix(mod_60)
conMatrix(mod_72)


# # save.image('/home/benbrew/Desktop/temp/main_model.RData')
load('/home/benbrew/Desktop/temp/main_model.RData')


##########
# get features 
##########
feature_list_reg <- list()
feature_list_48 <- list()
feature_list_60 <- list()
feature_list_72 <- list()


for (i in 1:10) {
  feature_list_reg[[i]] <- mod[[1]][[2]][[i]]
  feature_list_48[[i]] <- mod[[2]][[2]][[i]]
  feature_list_60[[i]] <- mod[[3]][[2]][[i]]
  feature_list_72[[i]] <- mod[[4]][[2]][[i]]
  
}

##########
# give a point to each probe everytime it shows up in each list
##########

getProbeScore <- function(data_list,name) {
  
  # loop throgh as assign run to a column
  for(i in 1:10){
    temp <- as.data.frame(data_list[[i]])
    temp$dat <- i
    names(temp) <- c('probe', 'dat')
    data_list[[i]] <- temp
  }
  
  # combine data 
  reg_full <- do.call(rbind, data_list)
  
  # create new column for score (how many times it shows up)
  # group by probe
  reg_group <- reg_full %>%
    group_by(probe) %>%
    summarise(counts = n())
  
  # order counts column to get top performers 
  reg_group <- reg_group[order(reg_group$counts, decreasing = T),]
  reg_group[,3] <- name
  
  return(reg_group)
  
}

full_reg <- getProbeScore(feature_list_reg, name = 'reg')
full_48 <- getProbeScore(feature_list_48, name = '48')
full_60 <- getProbeScore(feature_list_60, name = '60')
full_72 <- getProbeScore(feature_list_72, name = '72')

full_feat <- rbind(full_reg, full_48, full_60, full_72)

saveRDS(full_feat, '/home/benbrew/Desktop/full_feat.rda')

##########
# plot data
##########
# plot mod_reg: 1 is predictions, 3 is real 
x_axis <- unlist(mod_reg[[1]])
y_axis <- unlist(mod_reg[[3]])

cor(x_axis, y_axis)

plot_dat <- as.data.frame(cbind(x_axis, y_axis))

ggplot(plot_dat, aes(x_axis, y_axis)) + geom_point(alpha= 0.7) + xlab('Predicted age of onset') + ylab('Real age of onset') +
  ggtitle('Elastic Net Predictions') + theme_bw() + geom_abline(intercept = 0, slope = 1) + xlim(0,800) +
  ylim(0,800) + theme(text = element_text(size = 15))
  

##########
# plot data controls - pred vs sample
##########
x_con <- unlist(mod_reg[[6]])
y_con <- unlist(mod_reg[[7]])

plot_dat_con <- as.data.frame(cbind(x_con, y_con))

# get correlation
cor(x_con, y_con)

ggplot(plot_dat_con, aes(x_con, y_con)) + geom_point(alpha= 0.7) + xlab('Predictions') + ylab('Real age of sample collection') +
  ggtitle('Elastic Net Predictions') + theme_bw() + geom_abline(intercept = 0, slope = 1) + xlim(0,800) +
  ylim(0,800) + theme(text = element_text(size = 13))



#######
# WITHOUT gender as covariate

#### 10
mod_2 <- runModels(quan_cases_sam_gen, 
                 model = 'enet',
                 rand = F,
                 bump_hunter = T,
                 bump_hunter_data = quan_even_fwer_10,
                 gender= F,
                 residual = F)


mod_reg_2 <- mod_2[[1]]
mod_48_2 <- mod_2[[2]]
mod_60_2 <- mod_2[[3]]
mod_72_2 <- mod_2[[4]]

mod_reg_table_2 <- extractResults(mod_reg_2, 
                                data_name = 'sam_fwer_gender_10_enet', 
                                regularize = T,
                                bh_data = quan_even_fwer_10)

conMatrix(mod_48_2)
conMatrix(mod_60_2)
conMatrix(mod_72_2)


##########
# get features 
##########
feature_list_2 <- list()

for (i in 1:10) {
  feature_list_2[[i]] <- mod_2[[1]][[2]][[i]]
}

enet_features_2 <- Reduce(intersect, feature_list_2)


##########
# plot data
##########
# plot mod_reg: 1 is predictions, 3 is real 
x_axis_2 <- unlist(mod_reg_2[[1]])
y_axis_2 <- unlist(mod_reg_2[[3]])

plot_dat_2 <- as.data.frame(cbind(x_axis_2, y_axis_2))

ggplot(plot_dat_2, aes(x_axis_2, y_axis_2)) + geom_point(alpha= 0.7) + xlab('Model Predictions') + ylab('Real Age') +
  ggtitle('Elastic Net Predictions') + theme_bw() + geom_abline(intercept = 0, slope = 1) + xlim(0,600) 

# get regions
temp <- as.data.frame(temp)
colnames(temp) <- 'probe'
regions <- inner_join(temp, cg_locations, by = 'probe')

########################################################################################
# controls
# change main funcion to control =T
# and comment out class models

#### 10
mod <- runModels(quan_controls_gen, 
                 model = 'enet',
                 rand = F,
                 bump_hunter = T,
                 bump_hunter_data = quan_even_fwer_10,
                 gender= F,
                 residual = F)


mod_reg <- mod[[1]]
mod_48 <- mod[[2]]
mod_60 <- mod[[3]]
mod_72 <- mod[[4]]

mod_reg_table <- extractResults(mod, 
                                data_name = 'sam_fwer_gender_10_enet', 
                                regularize = T,
                                bh_data = quan_even_fwer_10)

conMatrix(mod_48)
conMatrix(mod_60)
conMatrix(mod_72)


##########
# get features 
##########
feature_list <- list()

for (i in 1:10) {
  feature_list[[i]] <- mod[[1]][[2]][[i]]
}

enet_features <- Reduce(intersect, feature_list)

##########
# plot data
##########
# plot mod_reg: 1 is predictions, 3 is real 
x_axis <- unlist(mod[[1]])
y_axis <- unlist(mod[[3]])

plot_dat <- as.data.frame(cbind(x_axis, y_axis))

ggplot(plot_dat, aes(x_axis, y_axis)) + geom_point(alpha= 0.7) + xlab('Model Predictions') + ylab('Real Age') +
  ggtitle('Elastic Net Predictions on controls') + theme_bw() + geom_abline(intercept = 0, slope = 1) + xlim(0,600)



