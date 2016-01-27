###########################################
# this script is used to subset the number of features in the methylation data so it's easier 
# to run the model.
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

##################################
# Using correlation

# Read in methylation data that was imputed on and clin to get label
methyl <- read.csv(paste0(data_folder, '/methyl_impute_raw.csv'))
clin <- read.csv(paste0(data_folder, '/clin.csv'))

# scale methyl
#methyl <- scale(methyl[, -1])

# make a correlation matrix 
cor_mat <- cor(methyl)

# find attributes that are highly correlated
highly_cor <- findCorrelation(cor_mat, cutoff = 0.8)

# remove highly correlated attributes
methyl_cor <- methyl[,-highly_cor]
names(methyl_cor)[1] <- 'id'
write.csv(methyl_cor, paste0(data_folder, '/methyl_cor.csv'), row.names = FALSE)

###########################################
# Recursive Feature elimination

# define the control using the random forest selection function
control <- rfeControl(functions = rfFuncs, method = 'cv', number = 10)

# make label
names(methyl)[1] <- 'id'
methyl$id <- as.factor(methyl$id)
# inner_join clin
model_data <- inner_join(clin, methyl,
                         by = 'id')
model_data$age_six <- ifelse(model_data$age_of_onset > 6, TRUE, FALSE)
model_data$age_six[is.na(model_data$age_six)] <- FALSE
label<- model_data$age_six
ground_truth <- as.factor(label)
table(ground_truth)

results <- rfe(model_data[,16:19806], ground_truth, size = c(1, 19791),
               rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))

##############################################
# Using low variance selection
#  First find the desired quantile breaks for the entire matrix
qt <- quantile(data.matrix(methyl[, -1]) , 0.1 )
# 20%  80% 
#5.17 6.62 
#  Next get a logical vector of the rows that have any values outside these breaks
columns <- apply( methyl[, -1] , 2, function(x) any( x < qt[1]))
#  Subset on this vector
methyl[ ,columns ]

# Use genefilter
temp <- varFilter(t(methyl[, -1]))

###############################################
# Using PCA
pca <- PCA(methyl[,-1])

#This line of code will sort the variables the most linked to each PC. 
# It is very useful when you have many variables.
temp <- dimdesc(pca)

###############################################
# nearZeroVar function
temp <- nearZeroVar(methyl[, -1], freqCut = 50/5, saveMetrics = TRUE)

temp <- apply(methyl[, -1], 2, function(i) var(i))


