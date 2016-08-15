### 
# this script will run grouped lasso
library(gglasso)
library(caret)
library(glmnet)
library(Metrics)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/regression_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


# read in methyl impute raw and full_data_cor
full_data <- read.csv(paste0(data_folder, '/full_data.csv'), stringsAsFactors = F)
full_data_cor <- read.csv(paste0(data_folder, '/full_data_cor.csv'), stringsAsFactors = F)
full_data_cor_small <- read.csv(paste0(data_folder, '/full_data_cor_small.csv'), stringsAsFactors = F)


# read in cluster labels
kmeans_full <- read.csv(paste0(data_folder, '/kmeans_full.csv'), stringsAsFactors = F)
hier_full <- read.csv(paste0(data_folder, '/hier_full.csv'), stringsAsFactors = F)
kmeans_cor <- read.csv(paste0(data_folder, '/kmeans_cor.csv'), stringsAsFactors = F)
hier_cor <- read.csv(paste0(data_folder, '/hier_cor.csv'), stringsAsFactors = F)
kmeans_cor_small <- read.csv(paste0(data_folder, '/kmeans_cor_small.csv'), stringsAsFactors = F)
hier_cor_small <- read.csv(paste0(data_folder, '/hier_cor_small.csv'), stringsAsFactors = F)


# add 1 to each label because grplasso cant deal with label 1. 
kmeans_full <- kmeans_full$labels 
hier_full <- hier_full$V1 
kmeans_cor <- kmeans_cor$labels
hier_cor <- hier_cor$V1 
kmeans_cor_small <- kmeans_cor_small$labels
hier_cor_small <- hier_cor_small$V1

# remove unnecessary columns 
full_data$X <- NULL
full_data_cor$X <- NULL
full_data_cor_small$X <- NULL


##############################################################################################
# function that takes data and index and returns predictions using group lasso
# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       group,
                       subset, 
                       selected_features,
                       iterations,
                       max) {
  model <- list()
  predictions <- list()
  mse <- list()
  test.ground_truth <- list()
  
  vars <- ncol(data[c(6,27:ncol(data))])

  genes <- colnames(data)[27:ncol(data)]
  
  data <- data[, c(subset, genes)]
  
  # Try the model with all different selection of features based on number of missinginess. 
  data <- data[complete.cases(data),]
  
  obs <- nrow(data)
  
  # convert characters to factors 
  for ( i in 1:ncol(data)){
    
    if(typeof(data[,i]) == 'character' || grepl('num', names(data[i]))) {
      data[,i] <- as.factor(data[,i])
      data[,i] <- as.numeric(data[,i])
      print(i)
    } 
  }
  
  
  num_fixed <- abs(vars - ncol(data))
  group <- append(rep(13, num_fixed), group)
  
  # set group 
  group <- as.numeric(group)
  
  for (i in 5:iterations){
    
    set.seed(i)
    train_index <- sample(nrow(data), nrow(data) *.7, replace = F)
    #run model
    y <- data$age_diagnosis[train_index]
    
    print(dim(data[train_index, c(selected_features, genes)]))
    print(length(group))
    print(length(y))
    
    model[[i]] = cv.gglasso(as.matrix(data[train_index, c(selected_features, genes)]), 
                            y, group=group, loss="ls",
                            pred.loss="L2", lambda.factor=0.05, nfolds=5, maxit = max)
    
    predictions[[i]] <- predict(model[[i]], 
                                newx = data.matrix(data[-train_index, c(selected_features, genes)]),
                                s = model[[i]]$lambda.min, type = 'link')
    
    test.ground_truth[[i]] <- data$age_diagnosis[-train_index]
    mse[[i]] <- rmse(unlist(predictions[[i]]), unlist(test.ground_truth[[i]]))
    
    print(i)
    
  }
  
  return(list(mse, predictions, model, test.ground_truth, obs))
  
}


########################################################################################################gdna.exon.intron

# variables missing
# gender 0
# gdna.base.change 164
# gdna.codon 164
# protein.codon.change 177
# gdna.exon.intron 492
# codon72.npro 517
# splice.delins.snv 519
# protein.codon.num 549
# mdm2.nG 652
##################################################################################################################3
# group lass
#############################
# cor small
# just methylation kmeans
grp_methyl_kmeans <- predictAll(data = full_data_cor_small,
                                group = kmeans_cor_small,
                                subset <- c("age_diagnosis"), 
                                selected_features = NULL,
                                iterations = 7,
                                max = 80000)

plot(unlist(grp_methyl_kmeans[[2]]), unlist(grp_methyl_kmeans[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)
# just methylation hier
grp_methyl_hier <- predictAll(data = full_data_cor_small,
                              group = hier_cor_small,
                              subset <- c("age_diagnosis"), 
                              selected_features = NULL,
                              iterations = 7,
                              max = 80000)

plot(unlist(grp_methyl_hier[[2]]), unlist(grp_methyl_hier[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)
# gender and gdna.base.change
grp_mut_kmeans <- predictAll(data = full_data_cor_small,
                             group = kmeans_cor_small,
                             subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                             selected_features = c("gender", "gdna.base.change"), 
                             iterations = 7,
                             max = 80000)

plot(unlist(grp_mut_kmeans[[2]]), unlist(grp_mut_kmeans[[4]]), xlab = 'Predictions', ylab = 'Actual Age')

# gender and gdna.base.change
grp_mut_hier <- predictAll(data = full_data_cor_small,
                           group = hier_cor_small,
                           subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                           selected_features = c("gender", "gdna.base.change"), 
                           iterations = 7,
                           max = 80000)

plot(unlist(grp_mut_hier[[2]]), unlist(grp_mut_hier[[4]]), xlab = 'Predictions', ylab = 'Actual Age')


# add gdna.codon
grp_mut1_kmeans <- predictAll(data = full_data_cor_small,
                              group = kmeans_cor_small,
                              subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                              selected_features = c("gender", "gdna.base.change", "gdna.codon"), 
                              iterations = 7,
                              max = 80000)

plot(unlist(grp_mut1_kmeans[[2]]), unlist(grp_mut1_kmeans[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add gdna.codon
grp_mut1_hier <- predictAll(data = full_data_cor_small,
                            group = hier_cor_small,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon"), 
                            iterations = 7,
                            max = 80000)

plot(unlist(grp_mut1_hier[[2]]), unlist(grp_mut1_hier[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

#HERE
# add protein.codon.change
grp_mut2_kmeans <- predictAll(data = full_data_cor_small,
                              group = kmeans_cor_small,
                              subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                              selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                              iterations = 7,
                              max = 80000)

plot(unlist(grp_mut2_kmeans[[2]]), unlist(grp_mut2_kmeans[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.change
grp_mut2_hier <- predictAll(data = full_data_cor_small,
                            group = hier_cor_small,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                            iterations = 7,
                            max = 80000)

plot(unlist(grp_mut2_hier[[2]]), unlist(grp_mut2_hier[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add gdna.exon.intron
grp_mut3_kmeans <- predictAll(data = full_data_cor_small,
                              group = kmeans_cor_small,
                              subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                   "gdna.exon.intron"), 
                              selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                             "gdna.exon.intron"), 
                              iterations = 7,
                              max = 80000)

plot(unlist(grp_mut3_kmeans[[2]]), unlist(grp_mut3_kmeans[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)



# add gdna.exon.intron
grp_mut3_hier <- predictAll(data = full_data_cor_small,
                            group = hier_cor_small,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron"), 
                            iterations = 7,
                            max = 80000)

plot(unlist(grp_mut3_hier[[2]]), unlist(grp_mut3_hier[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
grp_mut4_kmeans <- predictAll(data = full_data_cor_small,
                              group = kmeans_cor_small,
                              subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                          "gdna.exon.intron", "codon72.npro"), 
                              selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                    "gdna.exon.intron", "codon72.npro"), 
                              iterations = 7,
                              max = 80000)

plot(unlist(grp_mut4_kmeans[[2]]), unlist(grp_mut4_kmeans[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add codon72.npro
grp_mut4_hier <- predictAll(data = full_data_cor_small,
                            group = hier_cor_small,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro"), 
                            iterations = 7,
                            max = 80000)

plot(unlist(grp_mut4_hier[[2]]), unlist(grp_mut4_hier[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add splice.delins.snv
grp_mut5_kmeans <- predictAll(data = full_data_cor_small,
                              group = kmeans_cor_small,
                       subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                   "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                       selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                             "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                       iterations = 7,
                       max = 80000)

plot(unlist(grp_mut5_kmeans[[2]]), unlist(grp_mut5_kmeans[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add splice.delins.snv
grp_mut5_hier <- predictAll(data = full_data_cor_small,
                            group = hier_cor_small,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                            iterations = 7,
                            max = 80000)

plot(unlist(grp_mut5_hier[[2]]), unlist(grp_mut5_hier[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add protein.codon.num
grp_mut6_kmeans <- predictAll(data = full_data_cor_small,
                              group = kmeans_cor_small,
                              subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                          "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                              selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                    "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                              iterations = 7,
                              max = 80000)

plot(unlist(grp_mut6_kmeans[[2]]), unlist(grp_mut6_kmeans[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add protein.codon.num
grp_mut6_hier <- predictAll(data = full_data_cor_small,
                            group = hier_cor_small,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                            iterations = 7,
                            max = 80000)

plot(unlist(grp_mut6_hier[[2]]), unlist(grp_mut6_hier[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)


# add mdm2.nG
grp_mut7_kmeans <- predictAll(data = full_data_cor_small,
                              group = kmeans_cor_small,
                              subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                          "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                              selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                    "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                              iterations = 7,
                              max = 80000)

plot(unlist(grp_mut7_kmeans[[2]]), unlist(grp_mut7_kmeans[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

# add mdm2.nG
grp_mut7_hier <- predictAll(data = full_data_cor_small,
                            group = hier_cor_small,
                                    subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
                            iterations = 7,
                            max = 80000)

plot(unlist(grp_mut7_hier[[2]]), unlist(grp_mut7_hier[[4]]), xlab = 'Predictions', ylab = 'Actual Age')
abline(0,1)

#############################
# Corelation data
# just methylation kmeans
grp_methyl_kmeans_cor <- predictAll(data = full_data_cor,
                                group = kmeans_cor,
                                subset <- c("age_diagnosis"), 
                                selected_features = NULL,
                                iterations = 6,
                                max = 100000)


# just methylation hier
grp_methyl_hier_cor <- predictAll(data = full_data_cor,
                              group = hier_cor,
                              subset <- c("age_diagnosis"), 
                              selected_features = NULL,
                              iterations = 7,
                              max = 150000)


# gender and gdna.base.change
grp_mut_kmeans_cor <- predictAll(data = full_data_cor,
                             group = kmeans_cor,
                             subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                             selected_features = c("gender", "gdna.base.change"), 
                             iterations = 6,
                             max = 100000)


# gender and gdna.base.change
grp_mut_hier_cor <- predictAll(data = full_data_cor,
                           group = hier_cor,
                           subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                           selected_features = c("gender", "gdna.base.change"), 
                           iterations = 6,
                           max = 100000)

# add gdna.codon
grp_mut1_kmeans_cor <- predictAll(data = full_data_cor,
                              group = kmeans_cor,
                              subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                              selected_features = c("gender", "gdna.base.change", "gdna.codon"), 
                              iterations = 6,
                              max = 100000)


# add gdna.codon
grp_mut1_hier_cor <- predictAll(data = full_data_cor,
                            group = hier_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon"), 
                            iterations = 6,
                            max = 100000)

#HERE
# add protein.codon.change
grp_mut2_kmeans_cor <- predictAll(data = full_data_cor,
                              group = kmeans_cor,
                              subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                              selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                              iterations = 6,
                              max = 100000)

# add protein.codon.change
grp_mut2_hier_cor <- predictAll(data = full_data_cor,
                            group = hier_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                            iterations = 6,
                            max = 100000)


# add gdna.exon.intron
grp_mut3_kmeans_cor <- predictAll(data = full_data_cor,
                              group = kmeans_cor,
                              subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                          "gdna.exon.intron"), 
                              selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                    "gdna.exon.intron"), 
                              iterations = 6,
                              max = 100000)

# add gdna.exon.intron
grp_mut3_hier_cor <- predictAll(data = full_data_cor,
                            group = hier_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron"), 
                            iterations = 6,
                            max = 200000)

# add codon72.npro
grp_mut4_kmeans_cor <- predictAll(data = full_data_cor,
                              group = kmeans_cor,
                              subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                          "gdna.exon.intron", "codon72.npro"), 
                              selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                    "gdna.exon.intron", "codon72.npro"), 
                              iterations = 6,
                              max = 100000)


# add codon72.npro
grp_mut4_hier_cor <- predictAll(data = full_data_cor,
                            group = hier_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro"), 
                            iterations = 6,
                            max = 100000)

# add splice.delins.snv
grp_mut5_kmeans_cor <- predictAll(data = full_data_cor,
                              group = kmeans_cor,
                              subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                          "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                              selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                    "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                              iterations = 6,
                              max = 100000)

# add splice.delins.snv
grp_mut5_hier_cor <- predictAll(data = full_data_cor,
                            group = hier_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                            iterations = 6,
                            max = 100000)

# add protein.codon.num
grp_mut6_kmeans_cor <- predictAll(data = full_data_cor,
                              group = kmeans_cor,
                              subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                          "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                              selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                    "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                              iterations = 6,
                              max = 100000)


# add protein.codon.num
grp_mut6_hier_cor <- predictAll(data = full_data_cor,
                            group = hier_cor,
                            subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                            selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                  "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                            iterations = 6,
                            max = 100000)

# # add mdm2.nG
# grp_mut7_kmeans_cor <- predictAll(data = full_data_cor,
#                               group = kmeans_cor,
#                               subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                           "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                               selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                                     "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                               iterations = 6,
#                               max = 100000)
# 
# # add mdm2.nG
# grp_mut7_hier_cor <- predictAll(data = full_data_cor,
#                             group = hier_cor,
#                             subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                         "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                             selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                                   "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                             iterations = 6,
#                             max = 100000)

# combine all tests 
save.image(paste0(data_folder, '/grp_lasso.RData'))
load(paste0(data_folder, '/grp_lasso.RData'))

grp_all <- rbind (
  append('just_methyl', c(mean(unlist(grp_methyl_kmeans[[1]]), na.rm = T), grp_methyl_kmeans[[5]], 'small_kmeans')),
  append('just_methyl', c(mean(unlist(grp_methyl_hier[[1]]), na.rm = T), grp_methyl_hier[[5]], 'small_hier')),
  append('gender_and_gdna.base.change', c(mean(unlist(grp_mut_kmeans[[1]]), na.rm = T), grp_mut_kmeans[[5]], 'small_kmeans')),
  append('gender_and_gdna.base.change', c(mean(unlist(grp_mut_hier[[1]]), na.rm = T), grp_mut_hier[[5]], 'small_hier')),
  append('gdna.codon', c(mean(unlist(grp_mut1_kmeans[[1]]), na.rm = T), grp_mut1_kmeans[[5]], 'small_kmeans')),
  append('gdna.codon', c(mean(unlist(grp_mut1_hier[[1]]), na.rm = T), grp_mut1_hier[[5]], 'small_hier')),
  append('protein.codon.change', c(mean(unlist(grp_mut2_kmeans[[1]]), na.rm = T), grp_mut2_kmeans[[5]], 'small_kmeans')),
  append('protein.codon.change', c(mean(unlist(grp_mut2_hier[[1]]), na.rm = T), grp_mut2_hier[[5]], 'small_hier')),
  append('gdna.exon.intron', c(mean(unlist(grp_mut3_kmeans[[1]]), na.rm = T), grp_mut3_kmeans[[5]], 'small_kmeans')),
  append('gdna.exon.intron', c(mean(unlist(grp_mut3_hier[[1]]), na.rm = T), grp_mut3_hier[[5]], 'small_hier')),
  append('codon72.npro', c(mean(unlist(grp_mut4_kmeans[[1]]), na.rm = T), grp_mut4_kmeans[[5]], 'small_kmeans')),
  append('codon72.npro', c(mean(unlist(grp_mut4_hier[[1]]), na.rm = T), grp_mut4_hier[[5]], 'small_hier')),
  append('splice.delins.snv', c(mean(unlist(grp_mut5_kmeans[[1]]), na.rm = T), grp_mut5_kmeans[[5]], 'small_kmeans')),
  append('splice.delins.snv', c(mean(unlist(grp_mut5_hier[[1]]), na.rm = T), grp_mut5_hier[[5]], 'small_hier')),
  append('protein.codon.num', c(mean(unlist(grp_mut5_kmeans[[1]]), na.rm = T), grp_mut5_kmeans[[5]], 'small_kmeans')),
  append('protein.codon.num', c(mean(unlist(grp_mut5_hier[[1]]), na.rm = T), grp_mut5_hier[[5]], 'small_hier')),
  append('protein.codon.num', c(mean(unlist(grp_mut6_kmeans[[1]]), na.rm = T), grp_mut6_kmeans[[5]], 'small_kmeans')),
  append('protein.codon.num', c(mean(unlist(grp_mut6_hier[[1]]), na.rm = T), grp_mut6_hier[[5]], 'small_hier')),
  append('mdm2.nG', c(mean(unlist(grp_mut7_kmeans[[1]]), na.rm = T), grp_mut7_kmeans[[5]], 'small_kmeans')),
  append('mdm2.nG', c(mean(unlist(grp_mut7_hier[[1]]), na.rm = T), grp_mut7_hier[[5]], 'small_hier')),
  
  append('just_methyl', c(mean(unlist(grp_methyl_kmeans_cor[[1]]), na.rm = T), grp_methyl_kmeans_cor[[5]], 'small_kmeans_cor')),
  append('just_methyl', c(mean(unlist(grp_methyl_hier_cor[[1]]), na.rm = T), grp_methyl_hier_cor[[5]], 'small_hier_cor')),
  append('gender_and_gdna.base.change', c(mean(unlist(grp_mut_kmeans_cor[[1]]), na.rm = T), grp_mut_kmeans_cor[[5]], 'small_kmeans_cor')),
  append('gender_and_gdna.base.change', c(mean(unlist(grp_mut_hier_cor[[1]]), na.rm = T), grp_mut_hier_cor[[5]], 'small_hier_cor')),
  append('gdna.codon', c(mean(unlist(grp_mut1_kmeans_cor[[1]]), na.rm = T), grp_mut1_kmeans_cor[[5]], 'small_kmeans_cor')),
  append('gdna.codon', c(mean(unlist(grp_mut1_hier_cor[[1]]), na.rm = T), grp_mut1_hier_cor[[5]], 'small_hier_cor')),
  append('protein.codon.change', c(mean(unlist(grp_mut2_kmeans_cor[[1]]), na.rm = T), grp_mut2_kmeans_cor[[5]], 'small_kmeans_cor')),
  append('protein.codon.change', c(mean(unlist(grp_mut2_hier_cor[[1]]), na.rm = T), grp_mut2_hier_cor[[5]], 'small_hier_cor')),
  append('gdna.exon.intron', c(mean(unlist(grp_mut3_kmeans_cor[[1]]), na.rm = T), grp_mut3_kmeans_cor[[5]], 'small_kmeans_cor')),
  append('gdna.exon.intron', c(mean(unlist(grp_mut3_hier_cor[[1]]), na.rm = T), grp_mut3_hier_cor[[5]], 'small_hier_cor')),
  append('codon72.npro', c(mean(unlist(grp_mut4_kmeans_cor[[1]]), na.rm = T), grp_mut4_kmeans_cor[[5]], 'small_kmeans_cor')),
  append('codon72.npro', c(mean(unlist(grp_mut4_hier_cor[[1]]), na.rm = T), grp_mut4_hier_cor[[5]], 'small_hier_cor')),
  append('splice.delins.snv', c(mean(unlist(grp_mut5_kmeans_cor[[1]]), na.rm = T), grp_mut5_kmeans_cor[[5]], 'small_kmeans_cor')),
  append('splice.delins.snv', c(mean(unlist(grp_mut5_hier_cor[[1]]), na.rm = T), grp_mut5_hier_cor[[5]], 'small_hier_cor')),
  append('protein.codon.num', c(mean(unlist(grp_mut5_kmeans_cor[[1]]), na.rm = T), grp_mut5_kmeans_cor[[5]], 'small_kmeans_cor')),
  append('protein.codon.num', c(mean(unlist(grp_mut5_hier_cor[[1]]), na.rm = T), grp_mut5_hier_cor[[5]], 'small_hier_cor')),
  append('protein.codon.num', c(mean(unlist(grp_mut6_kmeans_cor[[1]]), na.rm = T), grp_mut6_kmeans_cor[[5]], 'small_kmeans_cor')),
  append('protein.codon.num', c(mean(unlist(grp_mut6_hier_cor[[1]]), na.rm = T), grp_mut6_hier_cor[[5]], 'small_hier_cor'))
  
  
)

#############################
# Corelation data
# just methylation kmeans
grp_methyl_kmeans_full <- predictAll(data = full_data,
                                    group = kmeans_full,
                                    subset <- c("age_diagnosis"), 
                                    selected_features = NULL,
                                    iterations = 6,
                                    max = 100000)


# just methylation hier
grp_methyl_hier_full <- predictAll(data = full_data,
                                  group = hier_full,
                                  subset <- c("age_diagnosis"), 
                                  selected_features = NULL,
                                  iterations = 6,
                                  max = 100000)


# gender and gdna.base.change
grp_mut_kmeans_full <- predictAll(data = full_data,
                                 group = kmeans_full,
                                 subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                                 selected_features = c("gender", "gdna.base.change"), 
                                 iterations = 6,
                                 max = 100000)


# gender and gdna.base.change
grp_mut_hier_full <- predictAll(data = full_data,
                               group = hier_full,
                               subset <- c("age_diagnosis", "gender", "gdna.base.change"), 
                               selected_features = c("gender", "gdna.base.change"), 
                               iterations = 6,
                               max = 100000)

# add gdna.codon
grp_mut1_kmeans_full <- predictAll(data = full_data,
                                  group = kmeans_full,
                                  subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                                  selected_features = c("gender", "gdna.base.change", "gdna.codon"), 
                                  iterations = 6,
                                  max = 100000)


# add gdna.codon
grp_mut1_hier_full <- predictAll(data = full_data,
                                group = hier_full,
                                subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon"), 
                                selected_features = c("gender", "gdna.base.change", "gdna.codon"), 
                                iterations = 6,
                                max = 100000)


# add protein.codon.change
grp_mut2_kmeans_full <- predictAll(data = full_data,
                                  group = kmeans_full,
                                  subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                                  selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                                  iterations = 6,
                                  max = 100000)

# add protein.codon.change
grp_mut2_hier_full <- predictAll(data = full_data,
                                group = hier_full,
                                subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                                selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change"), 
                                iterations = 6,
                                max = 100000)


# add gdna.exon.intron
grp_mut3_kmeans_full <- predictAll(data = full_data,
                                  group = kmeans_full,
                                  subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron"), 
                                  selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                        "gdna.exon.intron"), 
                                  iterations = 6,
                                  max = 100000)

# add gdna.exon.intron
grp_mut3_hier_full <- predictAll(data = full_data,
                                group = hier_full,
                                subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron"), 
                                selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                      "gdna.exon.intron"), 
                                iterations = 6,
                                max = 200000)

# add codon72.npro
grp_mut4_kmeans_full <- predictAll(data = full_data,
                                  group = kmeans_full,
                                  subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron", "codon72.npro"), 
                                  selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                        "gdna.exon.intron", "codon72.npro"), 
                                  iterations = 6,
                                  max = 100000)


# add codon72.npro
grp_mut4_hier_full <- predictAll(data = full_data,
                                group = hier_full,
                                subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro"), 
                                selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                      "gdna.exon.intron", "codon72.npro"), 
                                iterations = 6,
                                max = 100000)

# add splice.delins.snv
grp_mut5_kmeans_full <- predictAll(data = full_data,
                                  group = kmeans_full,
                                  subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                                  selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                                  iterations = 6,
                                  max = 100000)

# add splice.delins.snv
grp_mut5_hier_full <- predictAll(data = full_data,
                                group = hier_full,
                                subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                                selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv"), 
                                iterations = 6,
                                max = 100000)

# add protein.codon.num
grp_mut6_kmeans_full <- predictAll(data = full_data,
                                  group = kmeans_full,
                                  subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                              "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                                  selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                        "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                                  iterations = 6,
                                  max = 100000)


# add protein.codon.num
grp_mut6_hier_full <- predictAll(data = full_data,
                                group = hier_full,
                                subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                            "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                                selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
                                                      "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num"), 
                                iterations = 6,
                                max = 100000)

# # add mdm2.nG
# grp_mut7_kmeans_full <- predictAll(data = full_data,
#                               group = kmeans_full,
#                               subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                           "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                               selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                                     "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                               iterations = 6,
#                               max = 100000)
# 
# # add mdm2.nG
# grp_mut7_hier_full <- predictAll(data = full_data,
#                             group = hier_full,
#                             subset <- c("age_diagnosis", "gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                         "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                             selected_features = c("gender", "gdna.base.change", "gdna.codon", "protein.codon.change",
#                                                   "gdna.exon.intron", "codon72.npro", "splice.delins.snv", "protein.codon.num", "mdm2.nG"), 
#                             iterations = 6,
#                             max = 100000)


# grp_all <- rbind (
#   append('just_methyl', c(mean(unlist(grp_methyl_kmeans[[1]])), grp_methyl_kmeans[[5]], 'small_kmeans')),
#   append('just_methyl', c(mean(unlist(grp_methyl_hier[[1]])), grp_methyl_hier[[5]], 'small_hier')),
#   append('gender_and_gdna.base.change', c(mean(unlist(grp_mut_kmeans[[1]])), grp_mut_kmeans[[5]], 'small_kmeans')),
#   append('gender_and_gdna.base.change', c(mean(unlist(grp_mut_hier[[1]])), grp_mut_hier[[5]], 'small_hier')),
#   append('gdna.codon', c(mean(unlist(grp_mut1_kmeans[[1]])), grp_mut1_kmeans[[5]], 'small_kmeans')),
#   append('gdna.codon', c(mean(unlist(grp_mut1_hier[[1]])), grp_mut1_hier[[5]], 'small_hier')),
#   append('protein.codon.change', c(mean(unlist(grp_mut2_kmeans[[1]])), grp_mut2_kmeans[[5]], 'small_kmeans')),
#   append('protein.codon.change', c(mean(unlist(grp_mut2_hier[[1]])), grp_mut2_hier[[5]], 'small_hier')),
#   append('gdna.exon.intron', c(mean(unlist(grp_mut3_kmeans[[1]])), grp_mut3_kmeans[[5]], 'small_kmeans')),
#   append('gdna.exon.intron', c(mean(unlist(grp_mut3_hier[[1]])), grp_mut3_hier[[5]], 'small_hier')),
#   append('codon72.npro', c(mean(unlist(grp_mut4_kmeans[[1]])), grp_mut4_kmeans[[5]], 'small_kmeans')),
#   append('codon72.npro', c(mean(unlist(grp_mut4_hier[[1]])), grp_mut4_hier[[5]], 'small_hier')),
#   append('splice.delins.snv', c(mean(unlist(grp_mut5_kmeans[[1]])), grp_mut5_kmeans[[5]], 'small_kmeans')),
#   append('splice.delins.snv', c(mean(unlist(grp_mut5_hier[[1]])), grp_mut5_hier[[5]], 'small_hier')),
#   append('protein.codon.num', c(mean(unlist(grp_mut5_kmeans[[1]])), grp_mut5_kmeans[[5]], 'small_kmeans')),
#   append('protein.codon.num', c(mean(unlist(grp_mut5_hier[[1]])), grp_mut5_hier[[5]], 'small_hier')),
#   append('protein.codon.num', c(mean(unlist(grp_mut6_kmeans[[1]])), grp_mut6_kmeans[[5]], 'small_kmeans')),
#   append('protein.codon.num', c(mean(unlist(grp_mut6_hier[[1]])), grp_mut6_hier[[5]], 'small_hier')),
#   append('mdm2.nG', c(mean(unlist(grp_mut7_kmeans[[1]])), grp_mut7_kmeans[[5]], 'small_kmeans')),
#   append('mdm2.nG', c(mean(unlist(grp_mut7_hier[[1]])), grp_mut7_hier[[5]], 'small_hier')),
#   
#   append('just_methyl', c(mean(unlist(grp_methyl_kmeans_cor[[1]])), grp_methyl_kmeans_cor[[5]], 'small_kmeans_cor')),
#   append('just_methyl', c(mean(unlist(grp_methyl_hier_cor[[1]])), grp_methyl_hier_cor[[5]], 'small_hier_cor')),
#   append('gender_and_gdna.base.change', c(mean(unlist(grp_mut_kmeans_cor[[1]])), grp_mut_kmeans_cor[[5]], 'small_kmeans_cor')),
#   append('gender_and_gdna.base.change', c(mean(unlist(grp_mut_hier_cor[[1]])), grp_mut_hier_cor[[5]], 'small_hier_cor')),
#   append('gdna.codon', c(mean(unlist(grp_mut1_kmeans_cor[[1]])), grp_mut1_kmeans_cor[[5]], 'small_kmeans_cor')),
#   append('gdna.codon', c(mean(unlist(grp_mut1_hier_cor[[1]])), grp_mut1_hier_cor[[5]], 'small_hier_cor')),
#   append('protein.codon.change', c(mean(unlist(grp_mut2_kmeans_cor[[1]])), grp_mut2_kmeans_cor[[5]], 'small_kmeans_cor')),
#   append('protein.codon.change', c(mean(unlist(grp_mut2_hier_cor[[1]])), grp_mut2_hier_cor[[5]], 'small_hier_cor')),
#   append('gdna.exon.intron', c(mean(unlist(grp_mut3_kmeans_cor[[1]])), grp_mut3_kmeans_cor[[5]], 'small_kmeans_cor')),
#   append('gdna.exon.intron', c(mean(unlist(grp_mut3_hier_cor[[1]])), grp_mut3_hier_cor[[5]], 'small_hier_cor')),
#   append('codon72.npro', c(mean(unlist(grp_mut4_kmeans_cor[[1]])), grp_mut4_kmeans_cor[[5]], 'small_kmeans_cor')),
#   append('codon72.npro', c(mean(unlist(grp_mut4_hier_cor[[1]])), grp_mut4_hier_cor[[5]], 'small_hier_cor')),
#   append('splice.delins.snv', c(mean(unlist(grp_mut5_kmeans_cor[[1]])), grp_mut5_kmeans_cor[[5]], 'small_kmeans_cor')),
#   append('splice.delins.snv', c(mean(unlist(grp_mut5_hier_cor[[1]])), grp_mut5_hier_cor[[5]], 'small_hier_cor')),
#   append('protein.codon.num', c(mean(unlist(grp_mut5_kmeans_cor[[1]])), grp_mut5_kmeans_cor[[5]], 'small_kmeans_cor')),
#   append('protein.codon.num', c(mean(unlist(grp_mut5_hier_cor[[1]])), grp_mut5_hier_cor[[5]], 'small_hier_cor')),
#   append('protein.codon.num', c(mean(unlist(grp_mut6_kmeans_cor[[1]])), grp_mut6_kmeans_cor[[5]], 'small_kmeans_cor')),
#   append('protein.codon.num', c(mean(unlist(grp_mut6_hier_cor[[1]])), grp_mut6_hier_cor[[5]], 'small_hier_cor')),
#   append('mdm2.nG', c(mean(unlist(grp_mut7_kmeans_cor[[1]])), grp_mut7_kmeans_cor[[5]], 'small_kmeans_cor')),
#   append('mdm2.nG', c(mean(unlist(grp_mut7_hier_cor[[1]])), grp_mut7_hier_cor[[5]], 'small_hier_cor')),
#   
#   append('just_methyl', c(mean(unlist(grp_methyl_kmeans_full[[1]])), grp_methyl_kmeans_full[[5]], 'small_kmeans_full')),
#   append('just_methyl', c(mean(unlist(grp_methyl_hier_full[[1]])), grp_methyl_hier_full[[5]], 'small_hier_full')),
#   append('gender_and_gdna.base.change', c(mean(unlist(grp_mut_kmeans_full[[1]])), grp_mut_kmeans_full[[5]], 'small_kmeans_full')),
#   append('gender_and_gdna.base.change', c(mean(unlist(grp_mut_hier_full[[1]])), grp_mut_hier_full[[5]], 'small_hier_full')),
#   append('gdna.codon', c(mean(unlist(grp_mut1_kmeans_full[[1]])), grp_mut1_kmeans_full[[5]], 'small_kmeans_full')),
#   append('gdna.codon', c(mean(unlist(grp_mut1_hier_full[[1]])), grp_mut1_hier_full[[5]], 'small_hier_full')),
#   append('protein.codon.change', c(mean(unlist(grp_mut2_kmeans_full[[1]])), grp_mut2_kmeans_full[[5]], 'small_kmeans_full')),
#   append('protein.codon.change', c(mean(unlist(grp_mut2_hier_full[[1]])), grp_mut2_hier_full[[5]], 'small_hier_full')),
#   append('gdna.exon.intron', c(mean(unlist(grp_mut3_kmeans_full[[1]])), grp_mut3_kmeans_full[[5]], 'small_kmeans_full')),
#   append('gdna.exon.intron', c(mean(unlist(grp_mut3_hier_full[[1]])), grp_mut3_hier_full[[5]], 'small_hier_full')),
#   append('codon72.npro', c(mean(unlist(grp_mut4_kmeans_full[[1]])), grp_mut4_kmeans_full[[5]], 'small_kmeans_full')),
#   append('codon72.npro', c(mean(unlist(grp_mut4_hier_full[[1]])), grp_mut4_hier_full[[5]], 'small_hier_full')),
#   append('splice.delins.snv', c(mean(unlist(grp_mut5_kmeans_full[[1]])), grp_mut5_kmeans_full[[5]], 'small_kmeans_full')),
#   append('splice.delins.snv', c(mean(unlist(grp_mut5_hier_full[[1]])), grp_mut5_hier_full[[5]], 'small_hier_full')),
#   append('protein.codon.num', c(mean(unlist(grp_mut5_kmeans_full[[1]])), grp_mut5_kmeans_full[[5]], 'small_kmeans_full')),
#   append('protein.codon.num', c(mean(unlist(grp_mut5_hier_full[[1]])), grp_mut5_hier_full[[5]], 'small_hier_full')),
#   append('protein.codon.num', c(mean(unlist(grp_mut6_kmeans_full[[1]])), grp_mut6_kmeans_full[[5]], 'small_kmeans_full')),
#   append('protein.codon.num', c(mean(unlist(grp_mut6_hier_full[[1]])), grp_mut6_hier_full[[5]], 'small_hier_full')),
#   append('mdm2.nG', c(mean(unlist(grp_mut7_kmeans_full[[1]])), grp_mut7_kmeans_full[[5]], 'small_kmeans_full')),
#   append('mdm2.nG', c(mean(unlist(grp_mut7_hier_full[[1]])), grp_mut7_hier_full[[5]], 'small_hier_full'))
#   
# )
# 
