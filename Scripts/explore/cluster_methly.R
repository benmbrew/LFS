# this script we will remove the cancer signature sample about 1500 probes
# and make a clustered methylation heatmap

# initialize libraries
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(biovizBase)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)
library(preprocessCore)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(reshape2)
library(tissuesGeneExpression)
library(rafalib)

registerDoParallel(1)

# load data 
setwd('../predict_age/')

# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs
methyl_type <- 'm'
combat <- FALSE
if(combat){
  tech <- FALSE
} else {
  tech <- TRUE
}
gender <- TRUE
k_folds <- 10
beta_thresh <- 0.1
##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
##########
ratio_set <- readRDS('../../Data/g_ranges.rda')

# get g_ranges
g_ranges <- ratio_set

# get probes from rownames
g_ranges$probe <- rownames(ratio_set)

# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]


# read in all data
if(methyl_type == 'beta'){
  if(combat){
    all_cases <- readRDS('../../Data/all_cases_beta_combat.rda')
    all_con <- readRDS('../../Data/all_con_beta_combat.rda')
    message('loaded beta values, with combat correction')
    
  } else {
    all_cases <- readRDS('../../Data/all_cases_beta.rda')
    all_con <- readRDS('../../Data/all_con_beta.rda')
    
    # recode to homogenize with combat data
    all_cases$tech <- ifelse(all_cases$tech == '450k', 'batch_1', 'batch_2')
    all_con$tech <- ifelse(all_con$tech == '450k', 'batch_1', 'batch_2')
    
    message('loaded beta values, with no combat')
  }
} else {
  if(combat){
    all_cases <- readRDS('../../Data/all_cases_m_combat.rda')
    all_con <- readRDS('../../Data/all_con_m_combat.rda')
    message('loaded m values, with combat correction')
    
  } else {
    all_cases <- readRDS('../../Data/all_cases_m.rda')
    all_con <- readRDS('../../Data/all_con_m.rda')
    
    # recode to homogenize with combat data
    all_cases$tech <- ifelse(all_cases$tech == '450k', 'batch_1', 'batch_2')
    all_con$tech <- ifelse(all_con$tech == '450k', 'batch_1', 'batch_2')
    
    message('loaded m values, with no combat')
  }
}

# load bh features 
bh_feats <- readRDS('../../Data/bh_feats_cancer.rda')


# tech
all_cases <- cbind(as.data.frame(class.ind(all_cases$tech)), 
                   all_cases)

all_con <- cbind(as.data.frame(class.ind(all_con$tech)), 
                 all_con)

# rempove old tech variable 
all_cases$tech <- all_con$tech <- NULL

# gender
all_cases <- cbind(as.data.frame(class.ind(all_cases$gender)), 
                   all_cases)

all_con <- cbind(as.data.frame(class.ind(all_con$gender)), 
                 all_con)

# rempove old tech variable 
all_cases$gender <- all_con$gender <- NULL

# make ages back to months (times by 12)

# cases
all_cases$age_diagnosis <- 
  round(all_cases$age_diagnosis*12, 2)
all_cases$age_sample_collection <- 
  round(all_cases$age_sample_collection*12, 2)

# controls
all_con$age_diagnosis <- 
  round(all_con$age_diagnosis*12, 2)
all_con$age_sample_collection <- 
  round(all_con$age_sample_collection*12, 2)

# Remove NA in age of dianogsis for cases
all_cases <- all_cases[!is.na(all_cases$age_sample_collection),]
all_con <- all_con[!is.na(all_con$age_sample_collection),]


# data_cases <- all_cases
# data_controls <- all_con
# bh_features = bh_feats
# k_means_num = 10
# i = 1
# samples need to be in the columns
cluster_methyl <- function(data_cases, 
                           data_controls,
                           bh_features,
                           g_ranges,
                           k_means_num) {
  
  # get intersect_names
  intersect_names <- names(data_cases)[14:ncol(data_cases)]
  
  # get feature list
  colnames(bh_features)[1] <- 'chr'
  remove_features <- inner_join(bh_features, g_ranges)$probe
  
  # take remove features out of colnames 
  bh_features <- intersect_names[!intersect_names %in% remove_features]
  
  # get clin and probe data
  clin_data <- data_cases[, 1:13]
  probe_data <- data_cases[,14:ncol(data_cases)]
  
  # subset probes by bh_feats
  probe_data <- probe_data[, names(probe_data) %in% bh_features]
  
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
results <- cluster_methyl(data_cases = all_cases, 
                          data_controls = all_con, 
                          bh_features = bh_feats, 
                          g_ranges = g_ranges, 
                          k_means_num = 5)
cost_df <- results[[1]]
kmeans <- results[[2]]

names(cost_df) <- c("cluster", "cost")
plot(cost_df$cluster, cost_df$cost,
     bty = 'n', xlab = 'Clusters', ylab = 'Cost')

kmeans_labels <- cbind(kmeans$cluster)
summary(as.factor(kmeans_labels))

# change colnames
kmeans_labels <- as.data.frame(kmeans_labels)
colnames(kmeans_labels) <- 'cluster'
kmeans_labels$probe <- rownames(kmeans_labels)

# join with g_ranges
final_dat <- inner_join(kmeans_labels, g_ranges, by = 'probe')

set.seed(5)
final_dat <- final_dat[sample(nrow(final_dat), 2000),]


final_dat <- final_dat %>% group_by(seqnames) %>%
  mutate(counts = n())
final_dat$cluster <- as.factor(final_dat$cluster)

ggplot(final_dat, aes(seqnames, counts)) +
  geom_bar(stat = 'identity', fill = 'black', alpha = 0.7) +
  xlab('Chromosome') +
  ylab('# of probes in model') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(final_dat, aes(seqnames, counts, fill = cluster)) +
  geom_bar(stat = 'identity') + 
  scale_fill_manual(name = 'Cluster group',
                    values = brewer.pal(n = 5, name = 'Accent')) +
  xlab('Chromosome') +
  ylab('# of probes in model') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# make manhattan plot
library(qqman)
str(final_dat)
final_dat <- final_dat[sample(nrow(final_dat), 1500),]
final_dat$chrom <- as.numeric(gsub('chr', '', final_dat$seqnames))
final_dat$col_vec <- ifelse(final_dat$chrom %% 2 == 0, 'black', 'darkgrey')
ggplot(final_dat, aes(chrom, count, color = col_vec)) +
  geom_jitter(size = 1) + 
  scale_color_manual(name = '',
                     values = c('black', 'darkgrey'))


# ##################################################################
# # get intersect_names
# intersect_names <- names(all_cases)[14:ncol(all_cases)]
# 
# # get feature list
# colnames(bh_feats)[1] <- 'chr'
# remove_features <- inner_join(bh_feats, g_ranges)$probe
# 
# # take remove features out of colnames 
# bh_feats <- intersect_names[!intersect_names %in% remove_features]
# beta_mat <- all_cases[, bh_feats]
# beta_mat <- beta_mat[, sample(names(beta_mat), 1500)]
# beta_mat <- dist(t(beta_mat), method = 'euclidean')
# hc <- hclust(beta_mat, method = 'complete')
# hclusters <- cutree(hc, k = 5)
# 
# plot(hc, labels = hc$labels, cex = 0.01)
# heatmap.2(as.matrix(beta_mat))


# hclust
clin_dat <- names(all_cases)[1:13]
temp_plot <- all_cases[, c(clin_dat, sample(names(all_cases[,14:ncol(all_cases)]), 2000))]


d_temp <- dist(t(temp_plot[,14:ncol(temp_plot)]), method = 'euclidean') # method="man" # is a bit better
hc_temp <- hclust(d_temp, method = "complete")

dend <- as.dendrogram(hc_temp)
# order it the closest we can to the order of the observations:
library(dendextend)
dend <- rotate(dend, 1:2000)

# Color the branches based on the clusters:
dend <- color_branches(dend, k=5) #, groupLabels=iris_species)

# # Manually match the labels, as much as possible, to the real classification of the flowers:
# library(colorspace)
# labels_colors(dend) <-
#   rainbow_hcl(5)[sort_levels_values(
#     as.numeric(hc_cancer)[order.dendrogram(dend)]
#   )]

# # We shall add the flower type to the labels:
# labels(dend) <- paste(as.character(cancer_type)[order.dendrogram(dend)],
#                       "(",labels(dend),")", 
#                       sep = "")
# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
# dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend <- set(dend, "labels_cex", 0.5)

# And plot:
# par(mar = c(3,3,3,7))
# plot(dend, 
#      main = "Clustered methylation
#          (the labels give the true cancer type)", 
#      horiz =  TRUE,  nodePar = list(cex = .007))
# legend("topleft", legend = unique(cancer_type), fill = rainbow_hcl(3))

# library(gplots)


gplots::heatmap.2(as.matrix(t(temp_plot[, 14:ncol(temp_plot)])),  
                  main = "Heatmap for LFS",
                  srtCol = 20,
                  labRow = FALSE,
                  labCol = '',
                  xlab = 'Patient samples',
                  ylab = 'Probes used in model',
                  dendrogram = "row",
                  Rowv = dend,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  denscol = "grey",
                  key = F,
                  density.info = "density") # to add nice colored strips        

