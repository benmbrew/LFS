# this script will cluster on the methylation data and see if it is clustering by genomic location

source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


# set fixed variables
method = 'swan'
combat = 'combat_1'

# condition on fixed objects to get saving identifiers
which_methyl = 'beta'
beta_thresh = 0.05

cases_450 <- readRDS(paste0('../../Data/', method,'/cases_450_', combat,'.rda'))
cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_', combat,'.rda'))
con_850 <- readRDS(paste0('../../Data/', method,'/con_850_', combat,'.rda'))
con_mut <- readRDS(paste0('../../Data/', method,'/con_450_', combat,'.rda'))
con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_', combat,'.rda'))

##########
# read in age probes
##########
age_probes <- readRDS('../../Data/age_probes.rda')

##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
##########
g_ranges <- readRDS('../../Data/g_ranges.rda')

# get probes from rownames
g_ranges$probe <- rownames(g_ranges)

# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]

names(g_ranges)[1] <- 'chr'

##########
# create variables
##########

# load cases
cases_450 <- cbind(as.data.frame(class.ind(cases_450$gender)), 
                   cases_450)

# rempove old tech variable 
cases_450$gender <- NULL

# gender
con_850 <- cbind(as.data.frame(class.ind(con_850$gender)), 
                 con_850)

# rempove old tech variable 
con_850$gender <- NULL

# gender
cases_850 <- cbind(as.data.frame(class.ind(cases_850$gender)), 
                   cases_850)

# rempove old tech variable 
cases_850$gender <- NULL



# ge tgender 
con_wt <- cbind(as.data.frame(class.ind(con_wt$gender)), 
                con_wt)
con_mut <- cbind(as.data.frame(class.ind(con_mut$gender)), 
                 con_mut)

# rempove old tech variable 
con_wt$gender <- NULL
con_mut$gender <- NULL

# subset to get controls lfs and wild type
names(con_wt)[3] <- 'ids'
names(con_mut)[3] <- 'ids'

names(con_850)[3] <- 'ids'
names(cases_850)[3] <- 'ids'

# remove age from literature
clin_names <- names(cases_450)[!grepl('^cg', names(cases_450))]
feats <- names(cases_450)[grepl('^cg', names(cases_450))]
feats <- feats[!feats %in% age_probes]
cases_450 <- cases_450[, c(clin_names, feats)]
con_850 <- con_850[, c(clin_names, feats)]
cases_850 <- cases_850[, c(clin_names, feats)]
con_wt <- con_wt[, c(clin_names, feats)]
con_mut <- con_mut[, c(clin_names, feats)]


# add dummy tech variable for data sets with only one, replace family_name
names(cases_450)[9] <- 'tech'
names(con_850)[9] <- 'tech'
names(cases_850)[9] <- 'tech'

# fill them with Zero
cases_450$tech <- '450k'
con_850$tech <- '850k'
cases_850$tech <- '850k'

# do the same to con_mut and con_wt
names(con_mut)[9] <- 'tech'
names(con_wt)[9] <- 'tech'

# fill new variable with right tech indication
con_mut$tech <- '450k'
con_wt$tech <- '450k'

# run bumphunter on LFS healthy patients (LFS no cancer) and LFS cancer patients (LFS cancer)
bh_feats <- bump_hunter(dat_1 = con_wt, 
                        dat_2 = con_mut, 
                        bump = 'lfs', 
                        boot_num = 5, 
                        beta_thresh = beta_thresh,
                        methyl_type = methyl_type,
                        g_ranges = g_ranges)

# cases
cases_450 <- join_new_features(cases_450, new_features = bh_feats)
con_850 <- join_new_features(con_850, new_features = bh_feats)
cases_850 <- join_new_features(cases_850, new_features = bh_feats)
con_mut <- join_new_features(con_mut, new_features = bh_feats)
con_wt <- join_new_features(con_wt, new_features = bh_feats)

# lfs probes 
lfs_bump_probes <- colnames(cases_450)[grepl('^cg', colnames(cases_450))]

# temp_dat = cases_450
# bh_features = bh_feats
# k_means_num = 5

cluster_methyl <- function(temp_dat, 
                           bh_features,
                           g_ranges,
                           k_means_num) {
  
  
  clin_data <- temp_dat[, !grepl('^cg', names(temp_dat))]
  probe_data <- temp_dat[, grepl('^cg', names(temp_dat))]
  
  # get intersect_names
  intersect_names <- names(temp_dat)[grepl('^cg', names(temp_dat))]
  
  # get feature list
  # colnames(bh_features)[1] <- 'chr'
  # remove_features <- inner_join(bh_features, g_ranges)$probe
  # 
  # # take remove features out of colnames 
  # bh_features <- intersect_names[!intersect_names %in% remove_features]

  # 
  # # subset probes by bh_feats
  # probe_data <- probe_data[, names(probe_data) %in% bh_features]
  # 
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
results <- cluster_methyl(temp_dat =cases_450, 
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
final_dat <- inner_join(temp_importance, g_ranges, by = 'probe')

final_dat <- final_dat %>% group_by(chr) %>%
  summarise(mean_score = mean(score),
            counts = n())

set.seed(5)
# final_dat <- final_dat[sample(nrow(final_dat), 2000),]


final_dat <- final_dat %>% group_by(chr,cluster) %>%
  summarise(counts = n())
final_dat$cluster <- as.factor(final_dat$cluster)

options(scipen = 999)
ggplot(final_dat, aes(chr, counts)) +
  geom_bar(stat = 'identity', fill = 'black', alpha = 0.7) +
  xlab('Chromosome') +
  ylab('# of probes in model') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(final_dat, aes(chr, counts, fill = cluster)) +
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
# final_dat <- final_dat[sample(nrow(final_dat), 1500),]
final_dat$chrom <- as.numeric(gsub('chr', '', final_dat$chr))
final_dat$col_vec <- ifelse(final_dat$chrom %% 2 == 0, 'black', 'darkgrey')
ggplot(final_dat, aes(chrom, counts, color = col_vec)) +
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
k_num = 5

clin_dat <- cases_450[, !grepl('^cg', names(cases_450))]

temp_plot <-cases_450[, grepl('^cg', names(cases_450))]


d_temp <- dist(t(temp_plot), method = 'euclidean') # method="man" # is a bit better
hc_temp <- hclust(d_temp, method = "complete")

dend <- as.dendrogram(hc_temp)
# order it the closest we can to the order of the observations:
library(dendextend)
dend <- rotate(dend, order = 1:ncol(temp_plot))

# Color the branches based on the clusters:
dend <- color_branches(dend, k=k_num) #, groupLabels=iris_species)

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


gplots::heatmap.2(as.matrix(t(temp_plot)),  
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

