###############################################
# This script will subset features and merge with clin saving different versions of the data to be
# this is the 6th step in the pipeline
library(dplyr)
library(stringr)
library(impute)
library(mlbench)
library(caret)
library(FactoMineR)
library(genefilter)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')


#################################################################################################
# Read in methyl and clinical data and join by ids
#################################################################################################

# Read in data (clinical or clinical_two)
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = TRUE)

# load methylation data imputed with lsa on genes
methyl <- read.csv(paste0(data_folder, '/methyl_impute_raw.csv'))
methyl$X <- NULL

# load methylation data imputed with knn on probes
load(paste0(data_folder, '/methyl_knn.RData'))

# remove 'A, B and _ in methylation names and then drop the 5th character 
methyl$id <- gsub('A|B|_', '', methyl$id)
methyl$id <- substr(methyl$id, 1,4) 

methylation$id <- gsub('A|B|_', '', methylation$id)
methylation$id <- substr(methylation$id, 1,4) 

clin$id <- gsub('A|B_', '', clin$blood_dna_malkin_lab_)

# inner_join clin
full_data <- inner_join(clin, methyl,
                        by = 'id')

full_data <- full_data[!is.na(full_data$age_diagnosis),]
full_data <- full_data[!duplicated(full_data$blood_dna_malkin_lab_),]

full_data_mut <- full_data[full_data$p53_germline == 'Mut',]
full_data_mut <- full_data_mut[!duplicated(full_data_mut$tm_donor_),]
full_data_mut <- full_data_mut[, c(6,8,29:ncol(full_data_mut))]


# inner_join clin
full_data_probe <- inner_join(clin, methylation,
                        by = 'id')


# subset by complete age of diagnosisremove duplicates 
full_data <- full_data[!is.na(full_data$age_diagnosis),]
full_data <- full_data[!duplicated(full_data$id),]

full_data_probe <- full_data_probe[!is.na(full_data_probe$age_diagnosis),]
full_data_probe <- full_data_probe[!duplicated(full_data_probe$id),]

full_data_probe_mut <- full_data_probe[full_data_probe$p53_germline == 'Mut',]
full_data_probe_mut <- full_data_probe_mut[!duplicated(full_data_probe_mut$tm_donor_),]
full_data_probe_mut <- full_data_probe_mut[, c(6,8,29:ncol(full_data_probe_mut))]


# subset full_data to just have age data and methylation data 
full_data <- full_data[, c(6,8,29:ncol(full_data))]
full_data_probe <- full_data_probe[, c(6,8,29:ncol(full_data_probe))]


# Save data to be used later
write.csv(full_data, paste0(data_folder, '/full_data.csv'))
write.csv(full_data_probe, paste0(data_folder, '/full_data_probe.csv'))
write.csv(full_data, paste0(data_folder, '/full_data_mut.csv'))
write.csv(full_data_probe, paste0(data_folder, '/full_data_probe_mut.csv'))

save.image('/home/benbrew/Desktop/model_data2.RData')
###################################################################################################
# Using correlation
# scale methyl
#methyl <- scale(methyl[, -1])

# make a correlation matrix 
cor_mat <- cor(methyl[, -1])

# find attributes that are highly correlated
# 0.6 for methyl_cor, 0.4, methyl_cor_small and full_data_cor_small
highly_cor <- findCorrelation(cor_mat, cutoff = 0.5, names = TRUE)
cor_index <- names(methyl) %in% highly_cor[2:length(highly_cor)]

# remove highly correlated attributes
methyl_cor <- methyl[, !cor_index]
names(methyl_cor)[1] <- 'id'

# write.csv(methyl_cor, paste0(data_folder, '/methyl_cor_small.csv'))
# write.csv(methyl_cor, paste0(data_folder, '/methyl_cor.csv'))


# inner_join clin
full_data_cor <- inner_join(clin, methyl_cor,
                        by = 'id')

# Save data to be used later
# write.csv(full_data_cor, paste0(data_folder, '/full_data_cor_small.csv'))
# write.csv(full_data_cor, paste0(data_folder, '/full_data_cor.csv'))

######################################################################################################
# Random forest recursive feature elimination
# prepare training scheme
control <- trainControl(method="repeatedcv", number=5, repeats=3)

# subset full data for the model 
x_mat <- full_data[, c(6, 30:ncol(full_data))]
x_mat <- x_mat[complete.cases(x_mat),]
y <- as.numeric(x_mat$age_diagnosis)

# train the model
model <- train(x = x_mat[,-1],
               y = y,
               preProcess = 'scale',
               importance = TRUE,
               trControl = control)

# estimate variable importance
importance <- varImp(model)

# get vector of importance 
importance <- importance$importance

# make rownames a column
importance$gene <- rownames(importance)
rownames(importance) <- NULL

# sort importance vector
importance <- importance[order(importance$Overall, decreasing = T),]

# save importance vector
# write.csv(importance, paste0(data_folder, '/importance.csv'))
# plot importance, maybe cutoff around 20%
hist(importance$Overall)

# subset data by top features 
subset <- importance[importance$Overall > 26,]
full_data_rf <- full_data[, c(names(full_data)[1:30], subset$gene)]

write.csv(full_data_rf, paste0(data_folder, '/full_data_rf.csv'))

# #########################################################################################
# # Using low variance selection
# subset full data for the model 
x_mat <- full_data[, c(6, 30:ncol(full_data))]
x_mat <- x_mat[complete.cases(x_mat),]
y <- as.numeric(x_mat$age_diagnosis)

# #  First find the desired quantile breaks for the entire matrix
qt <- quantile(data.matrix(x_mat[, -1]) , 0.1 )
# # 20%  80%
# #5.17 6.62
# #  Next get a logical vector of the rows that have any values outside these breaks
columns <- apply(x_mat[, -1] , 2, function(x) any( x < qt[1]))
# #  Subset on this vector
temp <- x_mat[ ,columns]
#
# # Use genefilter
temp <- varFilter(t(temp[, -1]))

########################################################################################################
# # Using PCA
# pca <- PCA(x_mat[,-1])
#
# #This line of code will sort the variables the most linked to each PC.
# # It is very useful when you have many variables.
# temp <- dimdesc(pca)
# temp_new <- temp$Dim.1
# temp1 <- temp_new$quanti

#########################################################################################################
# # nearZeroVar function
# temp <- nearZeroVar(x_mat[, -1], freqCut = 50/5, saveMetrics = TRUE)
#
# temp <- apply(methyl[, -1], 2, function(i) var(i))


#########################################################################
# subset to just three columns- id, methlation_indicator, age_sample_collection, p53_germline, cancer_diagnosis
dat <- clin[, c('p53_germline', 'cancer_diagnosis_diagnoses', 'age_sample_collection', 'blood_dna_malkin_lab_',
                'methyl_indicator')]

range <- 12
samples <- list()
for ( i in 1:nrow(dat)) {
  
  if(dat$methyl_indicator[i] == 'Yes'){
    
    temp <- dat$age_sample_collection[i]
    

    samples[[i]] <- dat[(dat$age_sample_collection > temp & dat$age_sample_collection < (temp + range)) |
                  (dat$age_sample_collection < temp & dat$age_sample_collection > (temp - range)),] 
  }

}


samples <- do.call('rbind', samples)
samples <- samples[!duplicated(samples$blood_dna_malkin_lab_),]

# subset to mutant, unaffected, no methyl 
ids <- samples[samples$cancer_diagnosis_diagnoses == 'Unaffected' & samples$p53_germline == 'Mut' & 
                 samples$methyl_indicator == 'No',]

ids <- ids[complete.cases(ids),]

