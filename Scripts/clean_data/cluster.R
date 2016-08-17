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
full_data <- read.csv(paste0(data_folder, '/full_data.csv'), stringsAsFactors = F)
full_data_cor <- read.csv(paste0(data_folder, '/full_data_cor.csv'), stringsAsFactors = F)
full_data_cor_small <- read.csv(paste0(data_folder, '/full_data_cor_small.csv'), stringsAsFactors = F)
full_data_rf <- read.csv(paste0(data_folder, '/full_data_rf.csv'), stringsAsFactors = F)
full_data_rf_small <- read.csv(paste0(data_folder, '/full_data_rf_small.csv'), stringsAsFactors = F)



# remove X
full_data$X <- NULL
full_data_cor$X <- NULL
full_data_rf$X <- NULL
full_data_cor_small$X <- NULL
full_data_rf_small$X <- NULL


# hierarchical clustering function
# samples need to be in the columns
clusterMethyl <- function(data, cluster_size_kmeans, cluster_size_hier) {
  
  data <- data[!duplicated(data$blood_dna_malkin_lab_),]
  rownames(data) <- data[,7]
  data <- data[,-7]
  data_genes <- data[, 29:ncol(data)]
  
  data_genes <- as.matrix(data_genes)
  

  kmeans_data <- t(data_genes)
  
  #accumulator for cost results
  cost_df <- data.frame()
  
  #run kmeans for all clusters up to 100
  for(i in 1:cluster_size_kmeans){
    #Run kmeans for each level of i, allowing up to 100 iterations for convergence
    kmeans<- kmeans(kmeans_data, centers=cluster_size_kmeans, iter.max=100)
    
    #Combine cluster number and cost together, write to df
    cost_df<- rbind(cost_df, cbind(i, kmeans$tot.withinss))
    
    print(i)
  }
  
  names(cost_df) <- c("cluster", "cost")
  plot(cost_df$cluster, cost_df$cost,
       bty = 'n', xlab = 'Clusters', ylab = 'Cost')
  
  kmeans_labels <- cbind(kmeans$cluster)
  summary(as.factor(kmeans_labels))
  # nclust <- 2:10
  # # Perform clustering
  # clustMethods <- clValid(as.matrix(distance), 4,  clMethods = c("hierarchical"), 
  #                         validation = "internal")
  
  # Calculate the correlation between genes
  # The cor function takes correlations between columns
  correlation <- cor(data_genes)
  
  # Convert the correlation to a distance object
  distance <- as.dist(1 - correlation)
  
  # #use clValid or silhouetee to choose optimal clustering
  # The main function is clValid(), and the
  # available validation measures fall into the three general categories of “internal”, “stability”,
  hclustFit <- hclust(distance, method="average")
  labels <- cutree(hclustFit, k=cluster_size_hier)
  hier_labels <- cbind(labels)
  
  return(list(hier_labels, kmeans_labels))
}


full <- clusterMethyl(full_data, 6, 3)
kmeans_full <- full[[2]]
hier_full <- full[[1]]

rf <- clusterMethyl(full_data_rf, 6, 3)
kmeans_rf <- rf[[2]]
hier_rf <- rf[[1]]

rf_small <- clusterMethyl(full_data_rf_small, 6, 3)
kmeans_rf_small <- rf_small[[2]]
hier_rf_small <- rf_small[[1]]

cor <- clusterMethyl(full_data_cor, 6, 3)
kmeans_cor <- cor[[2]]
hier_cor <- cor[[1]]

cor_small <- clusterMethyl(full_data_cor_small, 6,3)
kmeans_cor_small <- cor_small[[2]]
hier_cor_small <- cor_small[[1]]


# # save labels
write.csv(kmeans_full, paste(data_folder,'kmeans_full.csv', sep ='/'))
write.csv(hier_full, paste(data_folder,'hier_full.csv', sep ='/'))

write.csv(kmeans_rf, paste(data_folder,'kmeans_rf.csv', sep ='/'))
write.csv(hier_rf, paste(data_folder,'hier_rf.csv', sep ='/'))

write.csv(kmeans_rf_small, paste(data_folder,'kmeans_rf_small.csv', sep ='/'))
write.csv(hier_rf_small, paste(data_folder,'hier_rf_small.csv', sep ='/'))

write.csv(kmeans_cor, paste(data_folder,'kmeans_cor.csv', sep ='/'))
write.csv(hier_cor, paste(data_folder,'hier_cor.csv', sep ='/'))

write.csv(kmeans_cor_small, paste(data_folder,'kmeans_cor_small.csv', sep ='/'))
write.csv(hier_cor_small, paste(data_folder,'hier_cor_small.csv', sep ='/'))

