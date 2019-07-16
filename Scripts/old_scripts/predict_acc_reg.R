
#!/hpf/tools/centos6/R/3.2.3/bin/Rscript

# script will predict age of onset within ACC only

argv <- as.numeric(commandArgs(T))

##########
# initialize libraries
##########
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(doParallel)
library(nnet)
library(dplyr)
library(bumphunter)
library(sqldf)
library(e1071)

registerDoParallel(1)
##########
# initialize folders
##########
home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'raw'
k = 10
##########
# load data
##########
# read in full m value data 
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m.rda')))
# betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new_m.rda'))) #34 449936
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new_m.rda')))

###########
# make id into ids
###########
colnames(betaCases)[1] <- 'ids'
# colnames(betaControls)[1] <- 'ids'
colnames(betaValid)[1] <- 'ids'

##########
# get model data
##########
betaCases <- getModData(betaCases)

##########
# remove inf
##########
betaCases <- removeInf(betaCases, probe_start = 8)
# betaControls <- removeInf(betaControls, probe_start = 8)
betaValid<- removeInf(betaValid, probe_start = 8)

#subset valid
betaValid <- betaValid[!betaValid$ids %in% betaCases$ids,]

##########
# get intersecting colnames and prepare data for modeling
##########
intersect_names <- Reduce(intersect, list(colnames(betaCases)[8:ncol(betaCases)],
                                          colnames(betaValid)[8:ncol(betaValid)]))

# assign dataframe identifier
betaCases$type <- '450k'

betaValid$type <- '850k'

# cases
betaCases <- betaCases[, c('ids',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender',
                           'type',
                           intersect_names)]


#validation
betaValid <- betaValid[, c('ids', 
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           'type',
                           intersect_names)]



###########
# get full data
###########
beta_full <- rbind(betaCases,
                   betaValid)

rm(betaCases,
   betaValid)


# summary of data types 
summary(as.factor(beta_full$type))

# get gender 
beta_full <- cbind(as.data.frame(class.ind(beta_full$gender)), beta_full)

##########

# betaCases <- betaCases[, c(1:2000, ncol(betaCases))]
trainTest <- function(cases, 
                      mod_feats,
                      diff_thresh,
                      feature_set,
                      seed_num, 
                      subset, 
                      k) 
{
  
  if(subset == '450k'){
    cases <- cases[grepl(subset, cases$type),]
  } 
  if(subset == '850k'){
    cases <- cases[grepl(subset, cases$type),]
  } 
  
  set.seed(seed_num)
  # get a column for each dataset indicating the fold
  cases <- getFolds(cases, seed_number = seed_num, k = k)
  
  # create variable to
  # summary(as.factor(cases$cancer_diagnosis_diagnoses))
  cases <- cases[grepl('ACC', cases$cancer_diagnosis_diagnoses),]
  
  # remove samples that dont have an age of sample collection
  cases <- cases[complete.cases(cases),]
  
  # list to store results
  model_results <- list()
  y_results <- list()
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, cases$folds)
    test_index <- !train_index
    
    
    # # get residuals
    # cases_resid <- getResidual(data = cases, 
    #                            bh_features = bh_feat_sig)
    
    mod_feats <- sample(mod_feats, feature_set, replace = T)
    
    mod_result <- runEnetDiff(training_dat = cases[train_index,], 
                              test_dat = cases[test_index,], 
                              bh_features = mod_feats,
                              gender = F)
    
    # predictions, age of onset, age of sample collection
    temp_results <- as.data.frame(t(do.call(rbind, mod_result)))
    temp_results$num_feats <- length(mod_feats)
    temp_results$seed_number <- j
    colnames(temp_results) <- c('preds', 'onset', 'age', 'num_feats', 'seed_num')
    y_results[[i]] <- temp_results
    
    # mod_result_resid <- runEnet(training_dat = cases_resid[train_index,], 
    #                             test_dat = cases_resid[test_index,], 
    #                             bh_features = bh_feat_sig,
    #                             gender = T)
    
    
    
  }
  
  return(y_results)
  
}

seed_results <- list()
full_results <- list()


feature_length <- c(100, 500, 1000, 5000, 10000, 20000)
seeds <- c(1, 2, 3)

for(l in 1:length(feature_length)) {
  
  for (j in 1:length(seeds)) {
    
    mod_results <- trainTest(cases = beta_full,
                             mod_feats = intersect_names,
                             diff_thresh = 12, 
                             feature_set = feature_length[l],
                             seed_num = seeds[j], 
                             subset = 'both',
                             k = k)
    
    seed_results[[j]] <- do.call(rbind, mod_results)
    
    print(paste0('done with ', j, ' seed'))
  }
  full_results[[l]] <- do.call(rbind, seed_results)
  print(paste0('done with ', l, ' random feature set'))
  
}


final_results <- do.call(rbind, full_results)

# saveRDS(final_results, paste0(model_data, '/predict_acc_results.rda'))

# change variable types 
names(final_results)
final_results$num_feats <- as.character(final_results$num_feats)
final_results$seed_num <- as.character(final_results$seed_num)

# remove duplicated onset age combo
final_results$dup <- paste0(final_results$num_feats, '_', final_results$onset, '_',
                            final_results$age, '_', final_results$seed_num)

# remove duplicate
final_results <- final_results[!duplicated(final_results$dup),]


results <- 
  final_results %>%
  group_by(num_feats, onset, age) %>%
  summarise(mean_pred = mean(preds),
            counts = n())

stopifnot(all(results$counts == length(seeds)))

##########
# plots 
##########

# abs
##########
# plot difference vs features

# get absolute value of difference
results$onset_diff <- abs(results$onset - results$mean_pred)
results$age_diff <- abs(results$age - results$mean_pred)

# # get rmse
# results$onset_diff <- sqrt((results$onset - results$mean_pred)^2)
# abs(results$onset - results$mean_pred)
# results$age_diff <- (results$age - results$mean_pred)^2

temp_result <- results[, c('num_feats', 'onset_diff', 'age_diff')]

# ave 
temp_result_group <- temp_result %>%
  group_by(num_feats) %>%
  summarise(mean_onset_diff = mean(onset_diff),
            mean_age_diff = mean(age_diff))

temp_result_group_melt <- melt(temp_result_group, id.vars = 'num_feats')
temp_result_group_melt$num_feats <- as.numeric(temp_result_group_melt$num_feats)

# plot difference
ggplot(temp_result_group_melt, aes(num_feats, value, group = variable, colour = variable)) +
  geom_point(alpha = 0.7, size = 3) + geom_line(size = 2, alpha=0.7) +
  xlab('Number of features') + ylab('Age errors (AV) (in months)') +
  scale_colour_manual(name = '',
                      breaks = c('mean_onset_diff', 'mean_age_diff'),
                      labels = c('Onset diff', 'Age diff'),
                      values = c('darkred', 'darkblue')) +
  theme_bw() + theme(axis.text.y = element_text(size = 12),
                     axis.text.x = element_text(size = 12)) 



# raw
##########
# plot difference vs features

# get absolute value of difference
results$onset_diff <- results$onset - results$mean_pred
results$age_diff <- results$age - results$mean_pred

temp_result <- results[, c('num_feats', 'onset_diff', 'age_diff')]

# ave 
temp_result_group <- temp_result %>%
  group_by(num_feats) %>%
  summarise(mean_onset_diff = mean(onset_diff),
            mean_age_diff = mean(age_diff))

temp_result_group_melt <- melt(temp_result_group, id.vars = 'num_feats')
temp_result_group_melt$num_feats <- as.numeric(temp_result_group_melt$num_feats)

# plot difference
ggplot(temp_result_group_melt, aes(num_feats, value, group = variable, colour = variable)) +
  geom_point(alpha = 0.7, size = 3) + geom_line(size = 2, alpha=0.7) +
  xlab('Number of features') + ylab('Age errors (AV) (in months)') +
  scale_colour_manual(name = '',
                      breaks = c('mean_onset_diff', 'mean_age_diff'),
                      labels = c('Onset diff', 'Age diff'),
                      values = c('darkred', 'darkblue')) +
  theme_bw() + theme(axis.text.y = element_text(size = 12),
                     axis.text.x = element_text(size = 12)) 

###########
# plot raw predicted and onset and age against number of features 
temp_results <- results[, 1:4]

# group by features and get means
temp_mean <- 
  temp_results %>%
  group_by(num_feats) %>%
  summarise(mean_onset = mean(onset),
            mean_age = mean(age),
            mean_pred = mean(mean_pred))


# melt it!
temp_melt <- melt(temp_mean, id.vars = 'num_feats')
temp_melt$num_feats <- as.numeric(temp_mean$num_feats)

# plot difference
ggplot(temp_melt, aes(num_feats, value, group = variable, colour = variable)) +
  geom_point(alpha = 0.7, size = 3) + geom_line(size = 2, alpha=0.7) +
  xlab('Number of features') + ylab('Age errors (AV) (in months)') +
  scale_colour_manual(name = '',
                      breaks = c('mean_onset', 'mean_age', 'mean_pred'),
                      labels = c('Avg age of diagnosis', 'Avg age', 'Avg prediction'),
                      values = c('darkblue', 'darkred', 'grey')) +
  theme_bw() + theme(axis.text.y = element_text(size = 12),
                     axis.text.x = element_text(size = 12)) 

############
# histograms
for (i in 1:length(feature_length)) {
  
  feature_sub <- feature_length[i]
  sub_dat <- temp_results[temp_results$num_feats == feature_sub, ]
  
  sub_dat$counts <-
    sub_dat$num_feats <- NULL
  
  sub_melt <- melt(sub_dat)
  title_name <- paste0(feature_sub, ' features')
  
  p <- ggplot(sub_melt, aes(value, fill=variable)) + 
    geom_histogram(bins = 25, alpha=0.6, position="identity", colour = 'grey') +
    scale_fill_manual(name = '',
                      breaks = c('onset', 'age', 'mean_pred'),
                      labels = c('age of diagnosis', 'age', 'Avg prediction'),
                      values = c('darkblue', 'darkred', 'grey')) +
    theme_bw() + theme(axis.text.y = element_text(size = 12),
                       axis.text.x = element_text(size = 12)) + ggtitle(title_name)
  
  print(p)
}






