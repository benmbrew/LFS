####### This script will merge the existing clinical data we have with the 
# existing by gene methylation data we have. 
library(stringr)
library(dplyr)
library(glmnet)
library(lme4)
library(impute)


# Initialize folders
home_folder <- '/home/benbrew/Documents'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')

################################################################
# Read in methyl and clinical data and join by ids
################################################################

# Read in data 
methyl <- read.csv(paste0(methyl_data, '/methyl.csv'), stringsAsFactors = FALSE)
clin <- read.csv(paste0(clin_data, '/clinical.csv'), stringsAsFactors = FALSE)

# remove everthing by numbers in the ids for both data sets 
removePun <- function(data){
  replace <- str_replace_all(data$id, "[[:punct:]]", " ")
  split <- strsplit(replace, ' ')
  combine_split <- lapply(split, function(x) x[1])
  new_ids <- unlist(combine_split)
  return(data)
}

clin <- removePun(clin)
methyl <- removePun(methyl)

# left_join clin
full_data <- inner_join(clin, methyl,
                       by = 'id')

###############################################################
# Run preliminary models using regularization
###############################################################

# Turn data into list of 3 sets
folds <- 3
full_data$folds <- sample(folds, nrow(full_data), replace = TRUE)
data_list <- list()
for(i in 1:folds){
  data_list[[i]] <- full_data[full_data$folds == i,]
}

# Make classifier for cancer or not 
full_data$cancer_indicator <- ifelse(full_data$cancer_indicator == 1, TRUE, FALSE)

# include all variable except fold.
# n <- names(full_data)[1:length(names(full_data))-1]
# formula <- as.formula(paste("cancer_indicator ~", paste(n[!n %in% "cancer_indicator"], collapse = " + ")))
# formula
# mod <- glm(cancer_indicator ~ ZWINT , data = full_data, lambda = 10, famliy = gaussian)
# impute NAS
# first remove rows that have NAs in excess of 30
# run model 
# fit_lasso <- cv.glmnet(model.matrix(cancer_indicator ~ ZWINT, data = x_train), 
#                        y_train, 
#                        family = 'binomial',
#                        nfolds = 5, 
#                        alpha = 1,
#                        standardize = F)

# radomly remove some columns before running. Currently cannot process 19804 columns in model.
full_data <- full_data[!is.na(full_data$cancer_indicator),]
x_train <- full_data[, 1:1000]
y_train <- as.factor(x_train[,4])
x_train$cancer_indicator <- as.factor(x_train$cancer_indicator)
x_train[is.na(x_train)] <- 0
length(y_train)

n <- names(x_train)[1:length(names(x_train))]
formula <- as.formula(paste("cancer_indicator ~", paste(n[!n %in% c("cancer_indicator", 
                                                                    "date", 
                                                                    "cancer", 
                                                                    "id", 
                                                                    "folds")], 
                                                        collapse = " + ")))


fit_ridge <- cv.glmnet(model.matrix(formula, data = x_train), 
                       y_train, 
                       family = 'binomial',
                       nfolds = 3,
                       alpha = 0)
plot(fit_ridge)
plot(fit_ridge)
print(fit_ridge)
fit_ridge$lambda.min

# kFolds <- length(folds) # list of data split into 5 
# 
# for (test.i in 1:kFolds) {
#   test <- folds[[test.i]]
#   train.folds <- folds[-test.i]
#   for (valid.i in 1:(kFolds - 1)) {
#     valid <- train.folds[[valid.i]]
#     train <- do.call("rbind", train.folds[-valid.i])
#     for (parameter.i in 1:length(parameters)) {
#       # train model on "train" with the this parameter
#       # predict outcomes on "valid" and save the score
#       #mod <- glmmLasso(fix = formula , rnd = NULL, data = clinical, lambda = 10)
#       
#     }
#   }
#   train.total <- do.call("rbind", train.folds)
#   # train model on "train.total" using the parameter value with the best average score on the validation sets
#   # predict outcomes on "test" and save the score
# }
# # average the scores on the test set to get an estimate of your procedure's ability to generalize to new data
# 
# # Note: you might not want to follow this procedure if you have a lot of data to start with



