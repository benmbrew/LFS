##########
# initialize libraries
##########
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(biovizBase)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)
library(preprocessCore)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(methyAnalysis)
library(reshape2)

registerDoParallel(1)

##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
idat_data <- paste0(methyl_data, '/raw_files')
model_data <- paste0(data_folder, '/model_data')


##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'raw'

##########
# read in idate for cases, controls, and validation set
##########
rgCases <- read.metharray.exp(idat_data)

##########
# get preprocedssing method
##########
betaCases <- preprocessMethod(rgCases, preprocess = method)
rm(rgCases)


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
# read in full m value data 
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m_scaled.rda')))

# samples need to be in the columns
cluster_methyl <- function(data, 
                           k_means_num) {
  
  
  # get clin and probe data
  clin_data <- betaCases[, 1:6]
  probe_data <- betaCases[,7:ncol(betaCases)]
  
  probe_data <- as.matrix(probe_data)
  
  
  # transpose so that probes are in the rows
  kmeans_data <- t(probe_data)
  
  #accumulator for cost results
  cost_df <- data.frame()
  
  #run kmeans for all clusters up to 100
  for(i in 1:k_means_num){
    #Run kmeans for each level of i, allowing up to 100 iterations for convergence
    kmeans<- kmeans(kmeans_data, centers= k_means_num, iter.max=100)
    
    #Combine cluster number and cost together, write to df
    cost_df<- rbind(cost_df, cbind(i, kmeans$tot.withinss))
    
    print(i)
  }
  
  return(list(cost_df, kmeans))
}

# get probes and corresponding clustre
results <- cluster_methyl(betaCases, k_means_num = 10)
cost_df <- results[[1]]
kmeans <- results[[2]]

names(cost_df) <- c("cluster", "cost")
plot(cost_df$cluster, cost_df$cost,
     bty = 'n', xlab = 'Clusters', ylab = 'Cost')

kmeans_labels <- cbind(kmeans$cluster)
summary(as.factor(kmeans_labels))

# change colnames
kmeans_labels <- as.data.frame(kmeans_labels)
colnames(kmeans_labels) <- 'count'
kmeans_labels$probe <- rownames(kmeans_labels)


saveRDS(kmeans_labels, paste0(model_data, paste0('/', method, '_', 'cases_kmeans_labs_scaled.rda')))

