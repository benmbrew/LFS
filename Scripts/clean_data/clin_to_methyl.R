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
test <- paste0(project_folder, '/Scripts/Analyze')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


#################################################################################################
# Read in methyl and clinical data and join by ids
#################################################################################################

# Read in data (clinical or clinical_two)
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = TRUE)
methyl <- read.csv(paste0(data_folder, '/methyl_impute_raw.csv'))
methyl$X <- NULL

# remove 'A' and '_' in methylation names
methyl$id <- gsub('_', '', methyl$id)
methyl$id <- gsub('A', '', methyl$id)

# make clin id a factor so it joins with methylation data
clin$id <- as.factor(clin$blood_dna_malkin_lab_)

# inner_join clin
full_data <- inner_join(clin, methyl,
                        by = 'id')

# Save data to be used later
# write.csv(full_data, paste0(data_folder, '/full_data.csv'))

###################################################################################################
# Using correlation
# scale methyl
#methyl <- scale(methyl[, -1])

# make a correlation matrix 
cor_mat <- cor(methyl[, -1])

# find attributes that are highly correlated
# 0.6 for methyl_cor, 0.4, methyl_cor_small and full_data_cor_small
highly_cor <- findCorrelation(cor_mat, cutoff = 0.4, names = TRUE)
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




