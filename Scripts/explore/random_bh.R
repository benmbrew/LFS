##########
# this script will look correlate random features to bh features and get union of all bh to remove from data for 
# radom models
library(SNFtool)
library(ggplot2)
library(reshape2)
##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')
scripts_folder <- paste0(project_folder, '/Scripts')

##########
# Read in cases and bumphunter features
##########
# bh_features
load(paste0(model_data, '/bh_features.RData'))

# cases - raw, swan, quan, funnorm, and clin
load(paste0(model_data, '/model_data_cases.RData'))

##########
# function taht gets data you want and removes other
##########
removeDat <- function(keep)
{
  data_set <- c('raw', 'quan', 'swan', 'funnorm')
  
  remove <- data_set[data_set != keep]
  
  # remove unwated
  rm(list=ls(pattern=remove[[1]], envir = .GlobalEnv), envir = .GlobalEnv)
  rm(list=ls(pattern=remove[[2]], envir = .GlobalEnv), envir = .GlobalEnv)
  rm(list=ls(pattern=remove[[3]], envir = .GlobalEnv), envir = .GlobalEnv)
  
}

removeDat(keep = 'raw')


##########
# function that subsets data by bh and correlates to random
##########

bhRandFinder <- function(beta_data, bh_data, title = NULL)
{
  
  set.seed(1000)
  # get beta feautes
  beta_features <- colnames(beta_data[, 6:ncol(beta_data)])
  
  # get features from bh
  bh_features <- bh_data$probe
  
  # get random subset of those features
  random_features <- sample(beta_features, length(bh_features))
  
  # subset data by bh_features
  bh_dat <- as.matrix(t(beta_data[, bh_features]))
  
  # subset data by random_features
  rand_dat <- as.matrix(t(beta_data[, random_features]))
  
  # get distance matrix
  cor_mat <- cor(bh_dat, rand_dat)
  
  # melt cor_mat 
  cor_mat_melt <- melt(cor_mat)
  
  # heatmap of distance matrix 
  ggplot(data = cor_mat_melt, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile(title) + xlab('Patients with BH Features') + ylab('Patients with Random Features')
    
}

# get bal_counts 
bhRandFinder(beta_funnorm, beta_bal_counts_cancer_intersection_features)
bhRandFinder(beta_funnorm, beta_bal_counts_cancer_intersection_sig_features)
bhRandFinder(beta_funnorm, beta_bal_counts_cancer_union_features)
bhRandFinder(beta_funnorm, beta_bal_counts_cancer_union_sig_features)
bhRandFinder(beta_funnorm, beta_bal_counts_p53_intersection_features)
bhRandFinder(beta_funnorm, beta_bal_counts_p53_intersection_sig_features)
bhRandFinder(beta_funnorm, beta_bal_counts_p53_union_features)
bhRandFinder(beta_funnorm, beta_bal_counts_p53_union_sig_features)

# cancer and p53 intersection and union

bhRandFinder(beta_funnorm, beta_cancer_intersection_features)
bhRandFinder(beta_funnorm, beta_cancer_intersection_sig_features)
bhRandFinder(beta_funnorm, beta_cancer_union_features)
bhRandFinder(beta_funnorm, beta_cancer_union_sig_features)
bhRandFinder(beta_funnorm, beta_p53_intersection_features)
bhRandFinder(beta_funnorm, beta_p53_intersection_sig_features)
bhRandFinder(beta_funnorm, beta_p53_union_features)
bhRandFinder(beta_funnorm, beta_p53_union_sig_features)

# data specific
                   
bhRandFinder(beta_funnorm, beta_funnorm_bal_cancer_features)
bhRandFinder(beta_funnorm, beta_funnorm_bal_cancer_sig_features)
bhRandFinder(beta_funnorm, beta_funnorm_bal_counts_cancer_features)
bhRandFinder(beta_funnorm, beta_funnorm_bal_counts_cancer_sig_features)
bhRandFinder(beta_funnorm, beta_funnorm_bal_counts_p53_features)
bhRandFinder(beta_funnorm, beta_funnorm_bal_counts_p53_sig_features)
bhRandFinder(beta_funnorm, beta_funnorm_bal_p53_features)
bhRandFinder(beta_funnorm, beta_funnorm_bal_p53_sig_features)

# more data specific

bhRandFinder(beta_funnorm, beta_funnorm_cancer_intersection_features)
bhRandFinder(beta_funnorm, beta_funnorm_cancer_intersection_sig_features)
bhRandFinder(beta_funnorm, beta_funnorm_cancer_union_features)
bhRandFinder(beta_funnorm, beta_funnorm_cancer_union_sig_features)
bhRandFinder(beta_funnorm, beta_funnorm_unbal_cancer_features)
bhRandFinder(beta_funnorm, beta_funnorm_unbal_cancer_sig_features)
bhRandFinder(beta_funnorm, beta_funnorm_unbal_p53_features)
bhRandFinder(beta_funnorm, beta_funnorm_unbal_p53_sig_features)


##########
# get union of all bh features
##########

beta_funnorm_union_features <- Reduce(union, list(beta_bal_counts_cancer_intersection_features$probe_rgSet,
                                              beta_bal_counts_cancer_intersection_sig_features$probe_rgSet,
                                              beta_bal_counts_cancer_union_features$probe_rgSet,
                                              beta_bal_counts_cancer_union_sig_features$probe_rgSet,
                                              beta_bal_counts_p53_intersection_features$probe_rgSet,
                                              beta_bal_counts_p53_intersection_sig_features$probe_rgSet,
                                              beta_bal_counts_p53_union_features$probe_rgSet,
                                              beta_bal_counts_p53_union_sig_features$probe_rgSet,
                                              beta_cancer_intersection_features$probe_rgSet,
                                              beta_cancer_intersection_sig_features$probe_rgSet, 
                                              beta_cancer_union_features$probe_rgSet,
                                              beta_cancer_union_sig_features$probe_rgSet,
                                              beta_p53_intersection_features$probe_rgSet,
                                              beta_p53_intersection_sig_features$probe_rgSet,
                                              beta_p53_union_features$probe_rgSet, 
                                              beta_p53_union_sig_features$probe_rgSet,
                                              beta_funnorm_bal_cancer_features$probe,
                                              beta_funnorm_bal_cancer_sig_features$probe, 
                                              beta_funnorm_bal_counts_cancer_features$probe,
                                              beta_funnorm_bal_counts_cancer_sig_features$probe,
                                              beta_funnorm_bal_counts_p53_features$probe,
                                              beta_funnorm_bal_counts_p53_sig_features$probe,
                                              beta_funnorm_bal_p53_features$probe,
                                              beta_funnorm_bal_p53_sig_features$probe,
                                              beta_funnorm_cancer_intersection_features, 
                                              beta_funnorm_cancer_intersection_sig_features$probe_rgSet,
                                              beta_funnorm_cancer_union_features$probe_rgSet,
                                              beta_funnorm_cancer_union_sig_features$probe_rgSet,
                                              beta_funnorm_unbal_cancer_features$probe_rgSet,
                                              beta_funnorm_unbal_cancer_sig_features$probe_rgSet,
                                              beta_funnorm_unbal_p53_features$probe_rgSet,
                                              beta_funnorm_unbal_p53_sig_features$probe_rgSet))

saveRDS(beta_funnorm_union_features, paste0(model_data, 'beta_funnorm_union_features.rda'))





