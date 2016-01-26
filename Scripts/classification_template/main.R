
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
classification_folder <- paste0(project_folder, '/Scripts/classification_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(project_folder, '/Results')
# Load Libraries and Dependencies ---------------------------------
library("igraph")
library("foreach")
library("ROCR")
library("doParallel")
library("caret")
library("glmnet")
library("randomForest")
library("kernlab")
library("pROC")
library(preprocessCore)

#### Set the "hyper" parameters 
registerDoParallel(1)
NUM_OF_PARTITION <- 2
## TO MODIFY:
# read the data and lables here 
# Load in data 
setwd(data_folder)
load(paste0(data_folder, '/methyl_lsa.RData'))

# use methyl, because more observations. make fake label 

# separate label and methylation data
label <- full_data$cancer_indicator
x_matrix <- methyl_impute[1:126, 2:50]

## Using CpG methylation:
x.methyl <- scale(x_matrix)
dim(x.methyl)

## labels
temp <- matrix(rexp(126, rate=1), ncol=1)
label<- ifelse(temp > 1, 'YES', 'NO')
ground_truth <- as.factor(label)
table(ground_truth)

#### Generate random partitions ---------------------------------
temp.data_ind <- 1:dim(x.methyl)[1]#x.expr

# 1:NUM_OF_PARTITION have random partition of 70% data for training
partition <- foreach (temp.run_ind = 1:NUM_OF_PARTITION, .errorhandling="stop") %dopar% {
  temp.good_split <- FALSE
  while (! temp.good_split) {
    temp.train_index <- sample(temp.data_ind, length(temp.data_ind) * 0.70)
    temp.complement <- setdiff(temp.data_ind, temp.train_index)
    temp.l.train <- length(unique(ground_truth[temp.train_index]))
    temp.l.test <- length(unique(ground_truth[temp.complement]))
    temp.good_split <- ( temp.l.train == temp.l.test && temp.l.test == length(levels(ground_truth)) )
  }

  ## Make test set without validation set
  temp.test_index <- temp.complement
  ## Split to validation and test sets
  # temp.valid_index <- sample(temp.complement, length(temp.data_ind) * 0.5)
  # temp.test_index <- setdiff(temp.complement, temp.valid_index)
  
  stopifnot(length(temp.train_index) > 0)
  stopifnot(length(temp.test_index) > 0)
  stopifnot(length(intersect(temp.test_index, temp.train_index)) == 0)
  stopifnot(temp.train_index %in% temp.data_ind)
  stopifnot(temp.test_index %in% temp.data_ind)

  # stopifnot(length(temp.valid_index) > 0)
  # stopifnot(length(intersect(temp.test_index, temp.valid_index)) == 0)
  # stopifnot(length(intersect(temp.train_index, temp.valid_index)) == 0)
  # stopifnot(temp.valid_index %in% temp.data_ind)

  list(
    train_index = temp.train_index,
    test_index = temp.test_index
    # ,valid_index = temp.valid_index
  )
}
rm(list = ls(pattern="temp*"))

#### Standard models make predictions ---------------------------------
source(paste0(classification_folder,'/standard_model_predict.R'))
source(paste0(classification_folder,"/plots.R"))

# for (temp.run_ind in 1:NUM_OF_PARTITION) { #to run in series 
models.methyl <- foreach (temp.run_ind = 1:NUM_OF_PARTITION, .errorhandling="stop") %dopar% { # run in parallel
  print(Sys.time())
  print(temp.run_ind)
  standard_model_predict(data = x.methyl, 
                        ground_truth = ground_truth, 
                        partition = partition,
                        selected_features = NULL, 
                        NFOLDS = 5, 
                        N_CV_REPEATS = 2, 
                        run_ind = temp.run_ind)
}
print("completed!")
# save(models.methyl, file="trained_modelsNresults.RData")
plot_models_performance(models.methyl,
                        NUM_OF_PARTITION,
                        "",
                        paste0(results_folder,"/Classifiers_test_results.pdf"),
                        "Classification using methyl. markers"
                        )

