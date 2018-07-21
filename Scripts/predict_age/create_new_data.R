# This script will: 
# 1) implement the horvath age calculator and then and save the age calculations as a feature in the new data set. 
# 2) regress each probe on the horvath caclulated age and extact the residuls
# 3) calculate the pc again and remove the first PC from the original data
library(WGCNA)

# source funcitons 
source('functions.R')

###----------------------------------------------------------------------------------- 
# apply horvaths method 



source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")



###------------------------------------------------------------------------------------
# regress probe on age calculatipons




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
