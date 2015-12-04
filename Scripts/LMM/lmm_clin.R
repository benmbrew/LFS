###### Run Linear mixed model with random effects only on clinical data 
# Use this script as a way to test the differences between lme3 and glmmLasso. 

# load cleaned clinical and methylation data (will only use clinical in this scripts)
setwd('/home/benbrew/Documents/LFS/Data')
load('cleaned.RData')