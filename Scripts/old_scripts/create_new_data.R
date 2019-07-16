# This script will: 
# 1) implement the horvath age calculator and then and save the age calculations as a feature in the new data set. 
# 2) regress each probe on the horvath caclulated age and extact the residuls
# 3) calculate the pc again and remove the first PC from the original data
library(WGCNA)
library(sqldf)
library(impute)

# source funcitons 
source('all_functions.R')

# read in all data
# cases
all_cases_beta <- readRDS('../../Data/all_cases_beta.rda')
all_cases_beta_combat < readRDS('../../Data/all_cases_beta_combat.rda')

# controls
all_con_beta <- readRDS('../../Data/all_con_beta.rda')
all_con_beta_combat <- readRDS('../../Data/all_con_beta_combat.rda')

# read m data
# cases
all_cases_m <- readRDS('../../Data/all_cases_m.rda')
all_cases_m_combat <- readRDS('../../Data/all_cases_m_combat.rda')

# controls
all_con_m <- readRDS('../../Data/all_con_m.rda')
all_con_m_combat <- readRDS('../../Data/all_con_m_combat.rda')

###----------------------------------------------------------------------------------- 
# apply horvaths method 


#Age transformation and probe annotation functions

trafo <-   function(x,adult.age=20) { 
  x=(x+1)/(1+adult.age)
  y=ifelse(x<=1, log(x), x-1)
  return(y)
  }
###------------------------------------------------------------------------------------
# regress probe on age calculatipons

# create function that takes a data set and reaturns the corresponding probes with residules from a 
# regression of each probe on age of sample collection, probes start at 12
temp_methyl <- all_cases_beta
get_residuals <- function(temp_methyl){
  # create list to store results 
  result_list <- list()
  
  # subset by complete outcome variable 
  temp_methyl
  
  # get outcome variable 
  outcome_variable <- temp_methyl$age_sample_collection
  
  # loop through each probe and regress it on age of sample collection to get residuls
  for(i in 12:ncol(temp_methyl)){
    probe <- temp_methyl[, i]
    lm(outcome_variable ~ probe)
  }
  
  
  
}




###----------------------------------------------------------------------------------
# remove first PC from analysis
# X = iris[,1:4]
# mu = colMeans(X)
# 
# Xpca = prcomp(X)
# 
# nComp = 2
# Xhat = Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
# Xhat = scale(Xhat, center = -mu, scale = FALSE)
# 
# Xhat[1,]
