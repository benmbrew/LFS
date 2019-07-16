
######
# This point of this script is to generate a "pretty" LFS plot for Anna
# directions: The whole point for a visual is to be a) pretty; b) related to the project. I’ll defer to your judgement regarding what to plot :)

# load library
library(tidyverse)
library(RColorBrewer)

# load wild type 
dat <- readRDS('../../Data/cases_all.rda')

# remove duplicates
dat <- dat[!duplicated(dat$tm_donor_),]

# get only cancer patients
dat <- dat[dat$cancer_diagnosis_diagnoses != 'Unaffected',]

# remove NAs
dat <- dat[!is.na(dat$cancer_diagnosis_diagnoses),]
dat <- dat[!is.na(dat$p53_germline),]
dat <- dat[!is.na(dat$age_sample_collection),]



##########
# get pca function
##########
pca_data <- dat
column_name <- 'p53_germline'

pca_data <- as.data.frame(pca_data)
pca_data[, column_name] <- as.factor(pca_data[, column_name])

# get cg cits
cg_sites <- names(pca_data)[grepl('^cg', names(pca_data))]
pca_data <- pca_data[ ,c(column_name, cg_sites)]

# recode type
pca_data$p53_germline <- ifelse(pca_data$p53_germline == 'Mut', 'LFS', 'Wild type')

# run pca
pca <- prcomp(pca_data[, 2:ncol(pca_data)])

# get pca dataframe with results and factor to color
pca_results <- data.frame(pca$x[, 1:2], column_name = pca_data[, column_name])


# get color
cols <- colorRampPalette(brewer.pal(n = 9, 'Set1'))(length(unique(pca_results$column_name)))

 
ggplot(pca_results, 
       aes(PC1, PC2, 
           col = column_name)) +
geom_point(size = 3, 
           alpha = 0.7) +
xlab('PC1') + ylab('PC2') +
scale_color_manual(name = '',
                   values = cols) + 
geom_hline(yintercept= 0, linetype="dashed", 
           color = "grey", size=1) +
geom_vline(xintercept=0, linetype="dashed", 
           color = "grey", size=1) +
theme_minimal(base_size = 18, base_family = 'ubuntu')


