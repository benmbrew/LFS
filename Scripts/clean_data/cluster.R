### 
# this script will cluster the imputed data 
# Step 5 in pipeline
library(clValid)
library(fpc)
library(cluster) 
library(dynamicTreeCut)
library(fpc)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/regression_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


# read in methyl impute raw
methyl_impute <- read.csv(paste0(methyl_data, '/methyl_impute_raw.csv'), stringsAsFactors = F)
methyl_cor <- read.csv(paste0(data_folder, '/methyl_cor.csv'), stringsAsFactors = F)
methyl_cor_small <- read.csv(paste0(data_folder, '/methyl_cor_small.csv'), stringsAsFactors = F)


# remove clinical variables 

# remove X
methyl_impute$X <- NULL
methyl_cor$X <- NULL
methyl_cor_small$X <- NULL


# put ids in rownames for imputation
rownames(methyl_impute) <- methyl_impute[,1]
methyl_impute <- methyl_impute[, -1]

methyl_impute <- as.matrix(methyl_impute)


rownames(methyl_cor) <- methyl_cor[,1]
methyl_cor <- methyl_cor[, -1]

methyl_cor <- as.matrix(methyl_cor)

rownames(methyl_cor_small) <- methyl_cor_small[,1]
methyl_cor_small <- methyl_cor_small[, -1]

methyl_cor_small <- as.matrix(methyl_cor_small)


# hierarchical clustering function
# samples need to be in the columns
clusterMethyl <- function(data, cluster_size) {
  
  
  # ...done in one step by function plotcluster.
  # fpc method 
  ancoord(data, clvecd, clnum=1, nn=50, method="mcd", countmode=1000)
  
  # kmeans
  kmeans_data <- t(data)
  
  #accumulator for cost results
  cost_df <- data.frame()
  
  #run kmeans for all clusters up to 100
  for(i in 2:cluster_size){
    #Run kmeans for each level of i, allowing up to 100 iterations for convergence
    kmeans<- kmeans(kmeans_data, centers=cluster_size, iter.max=100)
    
    #Combine cluster number and cost together, write to df
    cost_df<- rbind(cost_df, cbind(i, kmeans$tot.withinss))
    
    print(i)
  }
  
  names(cost_df) <- c("cluster", "cost")
  plot(cost_df$cluster, cost_df$cost,
       bty = 'n', xlab = 'Clusters', ylab = 'Cost')
  
  kmeans_labels <- cbind(kmeans$cluster)
  
  # nclust <- 2:10
  # # Perform clustering
  # clustMethods <- clValid(as.matrix(distance), 4,  clMethods = c("hierarchical"), 
  #                         validation = "internal")
  
  # Calculate the correlation between genes
  # The cor function takes correlations between columns
  correlation <- cor(data)
  
  # Convert the correlation to a distance object
  distance <- as.dist(1 - correlation)
  
  # #use clValid or silhouetee to choose optimal clustering
  # The main function is clValid(), and the
  # available validation measures fall into the three general categories of “internal”, “stability”,
  hclustFit <- hclust(distance, method="average")
  dynamicCut <- cutreeDynamic(hclustFit, minClusterSize=40, method="hybrid", 
                          distM = as.matrix(distance))
  labels <- cutree(hclustFit, k=415)
  hier_labels <- cbind(labels)
  
  return(list(hier_labels, kmeans_labels))
}




full <- clusterMethyl(methyl_impute, 415)
kmeans_full <- full[[1]]
hier_full <- full[[2]]

cor <- clusterMethyl(methyl_cor, 155)
kmeans_cor <- cor[[1]]
hier_cor <- cor[[2]]

cor_small <- clusterMethyl(methyl_cor_small, 30)
kmeans_cor_small <- cor_small[[1]]
hier_cor_small <- cor_small[[2]]


# # save labels 
# write.csv(kmeans_full, paste(data_folder,'kmeans_full.csv', sep ='/'))
# write.csv(hier_full, paste(data_folder,'hier_full.csv', sep ='/'))
# 
# # save labels 
# write.csv(kmeans_cor, paste(data_folder,'kmeans_cor.csv', sep ='/'))
# write.csv(hier_cor, paste(data_folder,'hier_cor.csv', sep ='/'))
# 
# save labels
# write.csv(kmeans_cor_small, paste(data_folder,'kmeans_cor_small.csv', sep ='/'))
# write.csv(hier_cor_small, paste(data_folder,'hier_cor_small.csv', sep ='/'))



