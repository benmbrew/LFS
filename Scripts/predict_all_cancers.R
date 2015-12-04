# Script for task one in LFS project

# 1. Will the tumor start before 5 years old?
# All cancers that start before 5 and all cancers that start after 5   
# (ACCearly+CPC+RMS vs the rest of the mutants) (UseLMM)
# this script is for reading in raw methylation files. 
setwd('/home/benbrew/Documents/LFS/Data/')
library(dplyr)
load('clinical.RData')

# clean clinical data
# plot cancer age variation 
# subset for tp53 mutants 
# age is response variable, do one with binary > 5 years and one as continous.
# protein and gdna are one to one mapping, use one, probably protein. 
# cancer is random effects, "whatever other clinical data that you don’t have for every 
# individual and potentially gender" summa
# turn data into list of k folds 
# female is 1, male is zero
# protein and rdna missing can be "not found"
# Different penalties- ridge, lasso, elastic net
# Treat missing as random effects in LMM
# start with just clinical and successivley add methylation
# compare lme4 and glmmLasso on basic model. 
# duplicate in clin because 1087 and 1087/1094.
# clean read in code and get rid of uncessary objects
# read paper in multiple random effects

duplicates <- data.frame(matrix(ncol = 12, nrow = 0))
# take clin ids that are separated by '/' and create a new row for them 
for (i in 1:nrow(clin_full)){
  sub_clin <- clin_full[i,]
  if(grepl('/', sub_clin$id)){
    split_id <- strsplit(as.character(sub_clin$id), '/')
    split_id <- cbind(unlist(split_id))
    duplicate <- cbind(split_id, sub_clin)
    duplicates <- rbind(duplicates, duplicate)
  }
}

# remove observations from clinical data that have '/' in id field and then rbind duplicates to 
# replace those observations 
clin_full <- clin_full[!grepl('/', clin_full$id),]
duplicates$id <- NULL
names(duplicates)[1] <- 'id'
clin_full <- rbind(clin_full, duplicates)

# Recode variables.  For now, don't treat 'blank' fields as NA. 
clin_full$pin_3 <- as.character(clin_full$pin_3)
clin_full$pin_3 <- as.factor(ifelse(clin_full$pin_3 == 'A1/A1?', 'A1/A1',
                          ifelse(clin_full$pin_3 == '', 'missing', clin_full$pin_3)))

clin_full$protein <- as.character(clin_full$protein)
clin_full$protein <- as.factor(ifelse(clin_full$protein == 'N/A', 'missing', clin_full$protein))

clin_full$codon_72 <- as.character(clin_full$codon_72)
clin_full$codon_72 <- as.factor(ifelse(clin_full$codon_72 == 'arg?', 'missing',
                                  ifelse(clin_full$codon_72 == 'arg homo', 'arg/arg',
                                    ifelse(clin_full$codon_72 == '', 'missing', clin_full$codon_72))))

clin_full$gender<- as.factor(ifelse(clin_full$gender == 1, 'female', 'male'))
clin_full$gender <- as.factor(clin_full$gender)

clin_full$mdm2 <- as.character(clin_full$mdm2)
clin_full$mdm2 <- as.factor(ifelse(clin_full$mdm2 == '', 'missing', clin_full$mdm2)) 

clin_full$age_five <- as.factor(ifelse(clin_full$age_of_onset > 5, TRUE, FALSE))

# Subset by tp53 mutants 
clin_mut <- clin_full[clin_full$tp53 == 'Mut',]


# Attach methylation data to clin_full 
combineMethyl <- function(methylation_data, num_sets){
for (i in 1:num_sets){
    sub_dat <- methylation_data[[i]]
    column_split <- strsplit(as.character(sub_dat$X), '#')
    last_digits <- lapply(column_split, function(x) x[length(x)])
    sub_ids <- unlist(last_digits)
    sub_ids <- gsub('RD-', '', sub_ids)
    methylation_data[[i]][,1] <- sub_ids
    if(i > 1){
      methylation_data[[i]] <- methylation_data[[i]][, -1]
    }
  }
methylation <- do.call('cbind', methylation_data)
return(methylation)
}

methylation <- combineMethyl(dat, num_sets = 23)

# can't left_join methylation to clin_full because it is too large
# subset methylation by clin_full ids, order, and cbind 
methylation_sub <- methylation[methylation[,1] %in% clin_mut$id,] # 60, but 55 because duplicates from 
# repeated samples on same id
clin_mut_sub <- clin_mut[clin_mut$id %in% methylation[,1],] # 56, but 55 unique because repeating from
# duplicated rows from '/'

#clinical has 78 observations where as previous cross with raw methylation showed 73. this 
# is because there are duplicates in clinical. The reasons there are duplicates is because 
# the raw methylation data has duplicates - same patient id, but different sample. 

# So the columns with multiple observations are id, cancer, gender and gender. Use
# missing values as multiple observation for codon_72, protein (gdna), pin_3, mdm2. 
# break into 5 groups (list) called folds 
num_folds <- 3
folds <- list()
clinical$num_folds <- as.factor(sample(1:num_folds, nrow(clinical), replace = TRUE))
for(i in 1:num_folds){
  folds[[i]] <- clinical[clinical$num_folds == i,]
}

kFolds <- length(folds)
library(glmnet)
library(lme4)
library(glmmLasso)

n <- names(clinical)
n <- n[2:10]
formula <- as.formula(paste("age_five ~", paste(n[!n %in% "age_five"], collapse = " + ")))
formula

mod <- glmmLasso(fix = formula , rnd = NULL, data = clinical, lambda = 10)

for (test.i in 1:kFolds) { # loops 3 times 
  test <- folds[[test.i]] # for test.i =1, this makes the first item in list the test set
  train.folds <- folds[-test.i] # for test.i = 1, this makes the other two the training set
  for (valid.i in 1:(kFolds - 1)) { # loops twice
    valid <- train.folds[[valid.i]] # for valid.i = 1, validation set is equl to one of the training sets
    train <- do.call("rbind", train.folds[-valid.i]) # or valid.i =1, training set equal to 
    # train.folds set not taken by validation. 
    for (parameter.i in 1:length(parameters)) {
      # train model on "train" with the this parameter
      # predict outcomes on "valid" and save the score
    }
  }
  train.total <- do.call("rbind", train.folds)
  # train model on "train.total" using the parameter value with the best average score on the validation sets
  # predict outcomes on "test" and save the score
}
# average the scores on the test set to get an estimate of your procedure's ability to generalize to new data

# Note: you might not want to follow this procedure if you have a lot of data to start with


