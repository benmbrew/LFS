
##########
# initialize libraries
##########
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(ModelMetrics)
library(doParallel)
library(nnet)
library(dplyr)
library(bumphunter)
library(sqldf)
library(e1071)
library(reshape2)

registerDoParallel(1)

##########
# initialize folders
##########

home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'raw'
k = 5

##########
# load data
##########
# read each data set in
betaFull <- readRDS(paste0(model_data, paste0('/', method, '_', 'full_mod.rda')))

betaFull$family_name.1 <- NULL

betaFull <- removeNA(betaFull, probe_start = 9) 

##########
# get column names
##########
intersect_names <- colnames(betaFull)[9:ncol(betaFull)]




mod_result <- mixed_lasso(training_dat = beta_cases[train_index,], 
                          controls_dat = beta_controls,
                          valid_dat = beta_valid,
                          test_dat = beta_cases[test_index,], 
                          age_cutoff = age_cutoff,
                          pred_cutoff = pred_cutoff,
                          bh_features = mod_feats,
                          gender = F)

# 
#  #Not run:
# N <- 20 # number of groups
# p <- 80 # number of covariates (including intercept)
# q <- 2 # number of random effect covariates
# ni <- rep(6,N) # observations per group
# n <- sum(ni) # total number of observations
# grp <- factor(rep(1:N,ni)) # grouping variable
# grp=rbind(grp,grp)
# beta <- c(1,2,4,3,rep(0,p-3)) # fixed-effects coefficients
# x <- cbind(1,matrix(rnorm(n*p),nrow=n)) # design matrix
# u1=rnorm(N,0,sd=sqrt(2))
# u2=rnorm(N,0,sd=sqrt(2))
# bi1 <- rep(u1,ni)
# bi2 <- rep(u2,ni)
# # lmme 5
# bi <- rbind(bi1,bi2)
# z=x[,1:2,drop=FALSE]
# epsilon=rnorm(120)
# y <- numeric(n)
# for (k in 1:n) y[k] <- x[k,]%*%beta + t(z[k,])%*%bi[,k] + epsilon[k]
# ########
# #independent random effects
# # x is a numeric matrix n*p (120, 81)
# # y is outcome, length 120
# # z is random effects matrix n*q (120, 2)
# # grp variable length n
# # rand where z is in x
# # fix variables (in front of data) not submitted for selection - use 1 or 2
# rand
# fit=lassop(x,y,z,grp,D=1,mu=0.2,fix=1,rand=c(1,2))
# # #dependent random effects
# fit=lassop(x,y,z,grp,mu=0.2,fix=1,rand=c(1,2))
# # ## End(Not run)







# get results 
temp_results[[i]] <- mod_result



# feature_length <- c(5000)
seeds <- c(1,2,3, 4,5)


