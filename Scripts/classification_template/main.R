### 
# this script will run the models

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/classification_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


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
library(dplyr)
library(preprocessCore)

#### Set the "hyper" parameters 
registerDoParallel(1)
NUM_OF_PARTITION <- 2
## TO MODIFY:

# Load in clinical data
clin <- read.csv(paste0(data_folder, '/clin.csv'))

##########################
# Load in methylation data. methyl_cor is a subset of features, excluding features with correlation 
# over .8

# methyl_impute <- read.csv(paste0(data_folder, '/methyl_impute.csv'))
# methyl_impute_raw <- read.csv(paste0(data_folder, '/methyl_impute_raw.csv'))
methyl_cor <- read.csv(paste0(data_folder, '/methyl_cor.csv'))
methyl_cor$id <- as.factor(methyl_cor$id)

#########################
## Create different variabes in clin to be used as label later
# acc or not 
# 
clin$acc_status <- ifelse(clin$cancer == 'ACC', TRUE, FALSE)
clin$age_six <- ifelse(clin$age_of_onset > 6, TRUE, FALSE)
clin$age_six[is.na(clin$age_six)] <- FALSE

# inner_join clin and methylation, retrieving only the ids in both
model_data <- inner_join(clin, methyl_cor,
                         by = 'id')

# Create label for the model using one of the variables created in clin

label<- model_data$age_six
ground_truth <- as.factor(label)
table(ground_truth)

# Select only the methylaion variables in model_data 
x_matrix <- model_data[1:71, 17:50]

# Scale data 
x.methyl <- scale(x_matrix)
x.methyl <- cbind(model_data$acc_status, x.methyl)
dim(x.methyl)

#### Generate random partitions ---------------------------------
temp.data_ind <- 1:dim(x.methyl)[1]

# 1:NUM_OF_PARTITION have random partition of 70% data for training
# foreach, run parallel 1:number of partitions, do the train and test split.
partition <- foreach (temp.run_ind = 1:NUM_OF_PARTITION, .errorhandling="stop") %dopar% {
  temp.good_split <- FALSE # set equal to false so while loop runs once
  while (! temp.good_split) {
    temp.train_index <- sample(temp.data_ind, length(temp.data_ind) * 0.70)# samples length of data 2/3 of data times
    temp.complement <- setdiff(temp.data_ind, temp.train_index) # gets index of obsverations not in train. 
    temp.l.train <- length(unique(ground_truth[temp.train_index])) # 2 labels
    temp.l.test <- length(unique(ground_truth[temp.complement])) # 2 labels
    temp.good_split <- ( temp.l.train == temp.l.test && temp.l.test == length(levels(ground_truth)) )
    # checks to see if label lengths are equal and that they are equal to ground_truth
    # while loops forces it to run until temp.good_split is true
  }

  ## Make test set without validation set
  temp.test_index <- temp.complement # this is the set to be tested on, the other 30%
  ## Split to validation and test sets
  # temp.valid_index <- sample(temp.complement, length(temp.data_ind) * 0.5)
  # temp.test_index <- setdiff(temp.complement, temp.valid_index)
  
  # stopifnot - if any of the expressions in ... are not all true, then stop
  stopifnot(length(temp.train_index) > 0) # train index must be greater than zero
  stopifnot(length(temp.test_index) > 0) # test index must be greater than zero
  stopifnot(length(intersect(temp.test_index, temp.train_index)) == 0) # there should be no overlap between
  # train and test inddex
  stopifnot(temp.train_index %in% temp.data_ind) # indexes must be in original index
  stopifnot(temp.test_index %in% temp.data_ind)

  # stopifnot(length(temp.valid_index) > 0)
  # stopifnot(length(intersect(temp.test_index, temp.valid_index)) == 0)
  # stopifnot(length(intersect(temp.train_index, temp.valid_index)) == 0)
  # stopifnot(temp.valid_index %in% temp.data_ind)

  list(
    train_index = temp.train_index, # concatenate train and test index into list
    test_index = temp.test_index
    # ,valid_index = temp.valid_index
  )
} # dopar agrument ends here and parition now has all the different train and test splits (equal to number of 
# partitions set.)
rm(list = ls(pattern="temp*")) # remove all temps

#### Standard models make predictions ---------------------------------
source(paste0(test,'/standard_model_predict.R'))
source(paste0(test,"/plots.R"))

# for (temp.run_ind in 1:NUM_OF_PARTITION) { #to run in series 
# another dopar that takes parition as an arugment. Refer to other script
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

