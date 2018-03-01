
##########
# load libraries
##########
library(tidyverse)
library(IlluminaHumanMethylation450kmanifest)
library(preprocessCore)
library(caret)
library(pROC)
library(doParallel) 
library(bumphunter)
library(e1071)
library(nnet)
library(glmnet)
library(PRROC)
library(survival)
library(broom)
library(randomForestSRC)
library(data.table)
library(c060)
library(glmmLasso)


##########
# function that Loops through list, preprocesses, and convert to beta, m, and cn values 
##########

preprocessMethod <- function(data, preprocess, methyl_type) {
  if (preprocess == 'raw') {
    Mset <- preprocessRaw(data)
  }
  if (preprocess == 'quan') {
    Mset <- preprocessQuantile(data, fixOutliers = TRUE,
                               removeBadSamples = TRUE, badSampleCutoff = 10.5,
                               quantileNormalize = TRUE, stratified = TRUE,
                               mergeManifest = FALSE, sex = NULL)
  }
  if (preprocess == 'illumina') {
    Mset  <- preprocessIllumina(data)
  } 
  if (preprocess == 'swan') {
    Mset  <-preprocessSWAN(data)
  }
  if (preprocess == 'funnorm') {
    Mset <-preprocessFunnorm(data)
  }
  if (preprocess == 'noob') {
    Mset <-preprocessNoob(data)
  }
  # map methyl set to genome (funnorm already does this)
  Gset <- mapToGenome(Mset)
  # # get m values
  if(methyl_type == 'beta') {
    dat <- getBeta(Gset)
    
  } else {
    
    dat <- getM(Mset)
    
  }
  return(dat)
}

# functions to be used in model_pipeline script
# data <- id_map
cleanIdMap <- function(data) {
  data <- as.data.frame(data)
  colnames(data) <- tolower(colnames(data))
  # combine sentrix_id and sentrix_position 
  data$identifier <- paste(data$sentrix_id, data$sentrix_position, sep = '_')
  data$identifier <- as.factor(data$identifier)
  return(data)
}


process_rg_set <-
  function(beta, id_map_1, id_map_2, clinical_dat, controls) {
    
    # get ids
    beta <- findIdsCombined(beta, id_map_1 = id_map_1, id_map_2 = id_map_2, controls = controls)
    
    if(controls) {
      # seperate beta 
      beta_cases <- beta[!grepl('^200', rownames(beta)),]
      beta_controls <- beta[grepl('^200', rownames(beta)),]
      
    } else {
      # seperate beta 
      beta_cases <- beta[!grepl('^20', rownames(beta)),]
      beta_controls <- beta[grepl('^20', rownames(beta)),]
      
    }
    
    
    # get id name (only cases)
    beta_cases <- getIdName(beta_cases)
    
    # clean ids
    beta_cases <- cleanIds(beta_cases)
    beta_controls <- cleanIds(beta_controls)
    
    
    # remove 'ch' from column names
    beta_cases <- beta_cases[, !grepl('ch', colnames(beta_cases))]
    beta_controls <- beta_controls[, !grepl('ch', colnames(beta_controls))]
    
    ##########
    # join data
    ##########
    
    # inner join
    beta_cases <- dplyr::inner_join(clinical_dat, beta_cases, by = 'ids')
    beta_controls <- dplyr::inner_join(clinical_dat, beta_controls, by = 'ids')
    
    
    # remove NAs from tm_donor 
    beta_cases <- beta_cases[!is.na(beta_cases$tm_donor_),]
    beta_controls <- beta_controls[!is.na(beta_controls$tm_donor_),]
    
    
    # remove duplicates
    beta_cases <- beta_cases[!duplicated(beta_cases$tm_donor_),]
    beta_controls <- beta_controls[!duplicated(beta_controls$tm_donor_),]
    
    
    ##########
    # get data in format for saving
    ##########
    
    # get cg_sites
    cg_sites <- colnames(beta)[grepl('cg', colnames(beta))]
    
    
    # subset data by colmns of interest and cg_sites
    beta_cases <- beta_cases[, c('ids', 
                                 'p53_germline', 
                                 'cancer_diagnosis_diagnoses', 
                                 'age_diagnosis',
                                 'age_sample_collection',
                                 'gender',
                                 'sentrix_id',
                                 'family_name',
                                 'tm_donor_',
                                 cg_sites)]
    
    # subset data by colmns of interest and cg_sites
    beta_controls <- beta_controls[, c('ids', 
                                       'p53_germline', 
                                       'cancer_diagnosis_diagnoses', 
                                       'age_diagnosis',
                                       'age_sample_collection',
                                       'gender',
                                       'sentrix_id',
                                       'family_name',
                                       'tm_donor_',
                                       cg_sites)]
    
    # combine data
    beta <- rbind(beta_cases,
                  beta_controls)
    
    # remove duplcated tm donor
    beta <- beta[!duplicated(beta$tm_donor_),]
    
    return(beta)
  }


# # id function
# beta_data <- beta_valid[1:10000,]
# id_map <- id_map_val
process_rg_set_single <- function(beta_data, id_map, clin) {
  beta_data <- findIds(beta_data, id_map)
  beta_data <- getIdName(beta_data)
  beta_data <- cleanIds(beta_data)
  beta_data <- beta_data[, !grepl('ch', colnames(beta_data))]
  beta_data <- dplyr::inner_join(clin, beta_data, by = 'ids')
  beta_data <- beta_data[!is.na(beta_data$tm_donor_),]
  beta_data <- beta_data[!duplicated(beta_data$tm_donor_),]
  cg_sites <- colnames(beta_data)[grepl('cg', colnames(beta_data))]
  beta_data <- beta_data[, c('ids', 
                             'p53_germline', 
                             'cancer_diagnosis_diagnoses', 
                             'age_diagnosis',
                             'age_sample_collection',
                             'gender',
                             'sentrix_id',
                             'family_name',
                             'tm_donor_',
                             'gdna.exon.intron',
                             'gdna.base.change',
                             cg_sites)]
  return(beta_data)
}



findIds <- function(data_methyl, id_map) {
  data_methyl <- as.data.frame(t(data_methyl))
  data_methyl$identifier <- rownames(data_methyl)
  data_methyl$identifier <- as.factor(data_methyl$identifier)
  # loop to combine identifiers, without merging large table
  data_methyl$ids <- NA
  data_methyl$sentrix_id <- NA
  for (i in data_methyl$identifier) {
    data_methyl$ids[data_methyl$identifier == i] <- as.character(id_map$sample_name[id_map$identifier == i])
    data_methyl$sentrix_id[data_methyl$identifier == i] <- as.character(id_map$sentrix_id[id_map$identifier == i])
    print(i)
  }
  return(data_methyl)
}


##########
# Function that gets malkin id from sample_name
##########
getIdName <- function(data) {
  column_split <- strsplit(as.character(data$ids), '#')
  last_digits <- lapply(column_split, function(x) x[length(x)])
  sub_ids <- unlist(last_digits)
  sub_ids <- gsub('RD-', '', sub_ids)
  data$ids <- sub_ids
  data$identifier <- NULL
  return(data)
}

# clean ids in each data set 
cleanIds <- function(data){
  data$ids <- gsub('A|B|_|-', '', data$ids)
  data$ids <- substr(data$ids, 1,4) 
  return(data)
}




# functions to be used in model_pipeline script
# data <- id_map
cleanIdMap <- function(data) {
  data <- as.data.frame(data)
  # new colnames, lowercase
  colnames(data) <- tolower(colnames(data))
  # combine sentrix_id and sentrix_position 
  data$identifier <- paste(data$sentrix_id, data$sentrix_position, sep = '_')
  data$identifier <- as.factor(data$identifier)
  return(data)
}


##########
# function for removing outlier from rgset
##########
# 
# rgSet <- rgControls
# id_map_dat <- id_map_con
# type = 'controls'
remove_outliers <- function(rgSet, id_map_dat, method, type) {
  # get outlier ids
  outliers <- data.frame(ids =c('3010','3391','3392','3540'),
                         batch = c('cases', 'controls', 'controls', 'valid')
  )
  # clean sample name
  column_split <- strsplit(as.character(id_map_dat$sample_name), '#')
  last_digits <- lapply(column_split, function(x) x[length(x)])
  sub_ids <- unlist(last_digits)
  sub_ids <- gsub('RD-', '', sub_ids)
  id_map_dat$ids <- sub_ids
  id_map_dat$ids <- gsub('A|B|_|-', '', id_map_dat$ids)
  id_map_dat$ids <- substr(id_map_dat$ids, 1,4) 
  # combine outliers and id map by 
  temp <- dplyr::inner_join(id_map_dat, outliers, by = 'ids')
  temp <- temp[grepl(type, temp$batch),]
  # get identifier 
  temp_id <- as.character(temp$identifier)
  # keep only
  rg_names <- colnames(rgSet)
  print(paste0(length(which(rg_names %in% temp_id)), ' found'))
  # intersecting index
  int_index <- rg_names %in% temp_id
  # subset rgSet by temp_id
  rgSet <- rgSet[, !int_index]
  return(rgSet)
}


##########
# function for subsetting rgsetts
##########

subset_rg_set <- function(rg_set, keep_gender, keep_controls, keep_snps, get_island, get_type, get_chr, gene_probes) {
  # get annotation
  temp_rg_set <- getAnnotation(rg_set)
  # syubet set to get relevant columns
  rg_dat <- temp_rg_set[, c('chr', 'Name', 'Type', 'NextBase', 'Color', 'Relation_to_Island')]
  if(!keep_gender) {
    # remove x and y chromosomes
    rg_dat_sub <- rg_dat[!grepl('X|Y', rg_dat$chr),]
  } else {
    rg_dat_sub <- rg_dat
  }
  # subset by conditions
  if (!is.null(get_island)) {
    rg_dat_sub <- rg_dat_sub[grepl(get_island, rg_dat_sub$Relation_to_Island),]
  } else if (!is.null(get_type)) {
    rg_dat_sub <- rg_dat_sub[grepl(get_type, rg_dat_sub$Type),]
  } else if (!is.null(get_chr)) {
    rg_dat_sub <- rg_dat_sub[grepl(get_chr, rg_dat_sub$chr),]
  }
  # now get character vectors for what probes to include and not include 
  keep_probes <- as.character(rg_dat_sub$Name)
  keep_int <- intersect(keep_probes, gene_probes)
  remove_probes <- rg_dat$Name[!as.character(rg_dat$Name) %in% as.character(rg_dat_sub$Name)]
  # use subset function from minfi
  rg_set_new <- subsetByLoci(rg_set, 
                             includeLoci = keep_int, 
                             excludeLoci = remove_probes, 
                             keepControls = keep_controls, 
                             keepSnps = keep_snps)
  return(rg_set_new)
}



# data_combat <- full_data_con

run_combat <- function(data_combat) {
  data_combat <- as.data.frame(data_combat)
  # get batch
  batch_indicator <- as.character(data_combat$batch)
  batch_indicator <- as.factor(batch_indicator)
  # put model ids in rownames and remove columns
  mat_combat <- as.matrix(data_combat[, 9:ncol(data_combat)])
  rownames(mat_combat) <- NULL
  clin_combat <- as.data.frame(data_combat[, 1:8])
  # get features
  features <- colnames(mat_combat)
  mat_combat <- t(mat_combat)
  # get intercept
  modcombat <- model.matrix(~1, data = data_combat)
  combat <- ComBat(dat = mat_combat, batch = batch_indicator, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
  any(is.na(mat_combat))
  any(is.na(batch_indicator))
  any(is.na(modcombat))
  # transpose and add back columns
  final_dat <- as.data.frame(t(combat))
  final_dat <- as.data.frame(cbind(clin_combat, final_dat))
  rownames(final_dat) <- NULL
  return(final_dat)
}

##########
# get family cancer status and ratio
##########
cases_full <- data_cases_full
controls_full <- data_controls_full
get_family_cancer <- function(cases_full, controls_full){
  #combine data subseted data
  temp_full <- bind_rows(cases_full[, c('tm_donor_', 'cancer_diagnosis_diagnoses', 'family_name')],
                         controls_full[, c('tm_donor_', 'cancer_diagnosis_diagnoses', 'family_name')])
  # create cancer indicator
  temp_full$cancer_fac <- ifelse(grepl('Unaffected', temp_full$cancer_diagnosis_diagnoses), 'no_cancer', 'cancer')
  temp_full$cancer_diagnosis_diagnoses <- NULL
  
  temp <- temp_full %>%
    group_by(family_name) %>%
    tally() %>% 
    left_join(temp_full)
  
  temp$num_with_cancer <- NA
  temp$cancer_ratio <- NA
 for(fam_name in unique(temp$family_name)){
   sub_fam <- temp[temp$family_name == fam_name,]
   if(nrow(sub_fam) > 1) {
     num_cancer <- length(which(sub_fam$cancer_fac == 'cancer'))
     num_no <- length(which(sub_fam$cancer_fac == 'no_cancer'))
     sub_fam$num_with_cancer <- ifelse(sub_fam$cancer_fac == 'cancer',
                                       num_cancer - 1, num_cancer)
     # HERE sub_fam$cancer_ratio <- 
   } else {
     sub_fam$num_with_cancer <- 0
     sub_fam$cancer_ratio <- 0
   }
   temp[temp$family_name == fam_name,] <- sub_fam
   
 }

}

##########
# Function that combines methylation matrices with id_map, to get ids for methylation
##########


findIdsCombined <- function(data_methyl, id_map_1, id_map_2, controls) {
  
  
  data_methyl <- as.data.frame(t(data_methyl))
  data_methyl$identifier <- rownames(data_methyl)
  data_methyl$identifier <- as.factor(data_methyl$identifier)
  # loop to combine identifiers, without merging large table
  data_methyl$ids <- NA
  data_methyl$sentrix_id <- NA
  
  for (i in data_methyl$identifier) {
    
    if(controls) {
      if (grepl('^200', i)) {
        data_methyl$ids[data_methyl$identifier == i] <- id_map_2$sample_name[id_map_2$identifier == i]
        data_methyl$sentrix_id[data_methyl$identifier == i] <- id_map_2$sentrix_id[id_map_2$identifier == i]
      } else {
        data_methyl$ids[data_methyl$identifier == i] <- id_map_1$sample_name[id_map_1$identifier == i]
        data_methyl$sentrix_id[data_methyl$identifier == i] <- id_map_1$sentrix_id[id_map_1$identifier == i]
        
      }
    } else {
      if (grepl('^20', i)) {
        data_methyl$ids[data_methyl$identifier == i] <- id_map_2$sample_name[id_map_2$identifier == i]
        data_methyl$sentrix_id[data_methyl$identifier == i] <- id_map_2$sentrix_id[id_map_2$identifier == i]
      } else {
        data_methyl$ids[data_methyl$identifier == i] <- id_map_1$sample_name[id_map_1$identifier == i]
        data_methyl$sentrix_id[data_methyl$identifier == i] <- id_map_1$sentrix_id[id_map_1$identifier == i]
        
      }
    }
    
    print(i)
    
  }
  
  return(data_methyl)
}


##########
# function to remove columns that have any NAs
##########
removeNA <- function(data_frame, probe_start) {
  
  # get full probe (all non missing columns) data set
  temp_data <- 
    data_frame[, probe_start:ncol(data_frame)][sapply(data_frame[, probe_start:ncol(data_frame)], 
                                                      function(x) all(!is.na(x)))]
  
  # combine probes with clin
  full_data <- as.data.frame(cbind(data_frame[, 1:(probe_start-1)], temp_data))
  
  # check that it worked
  stopifnot(all(!is.na(full_data[, probe_start:ncol(full_data)])))
  
  return(full_data)
  
}

##########
# function to remove columns that have any NAs
##########

removeInf <- function(data_frame, probe_start) {
  
  # get full probe (all non missing columns) data set
  temp_data <- 
    data_frame[, probe_start:ncol(data_frame)][, sapply(data_frame[, probe_start:ncol(data_frame)], 
                                                        function(x) all(!is.infinite((x))))]
  
  # combine probes with clin
  full_data <- as.data.frame(cbind(data_frame[, 1:(probe_start -1)], temp_data))
  
  # check that it worked
  # stopifnot(all(!is.na(full_data[, probe_start:ncol(full_data)])))
  
  return(full_data)
  
}


findIds <- function(data_methyl, id_map) {
  
  data_methyl <- as.data.frame(t(data_methyl))
  data_methyl$identifier <- rownames(data_methyl)
  data_methyl$identifier <- as.factor(data_methyl$identifier)
  # loop to combine identifiers, without merging large table
  data_methyl$ids <- NA
  data_methyl$sentrix_id <- NA
  
  for (i in data_methyl$identifier) {
    
    
    data_methyl$ids[data_methyl$identifier == i] <- as.character(id_map$sample_name[id_map$identifier == i])
    data_methyl$sentrix_id[data_methyl$identifier == i] <- as.character(id_map$sentrix_id[id_map$identifier == i])
    
    
    print(i)
    
  }
  
  return(data_methyl)
}




##########
# get mutant
##########
getModData <- function(data) 
{
  # subset data by not na in age of diagnosis and mut
  data <- data[!is.na(data$age_diagnosis),]
  data <- data[data$p53_germline == 'Mut',]
  return(data)
}

testKS <- function(x, y)
{
  y <- y[!is.na(y)]
  x <- x[!is.na(x)]
  
  # Do x and y come from the same distribution?
  ks.test(jitter(x), jitter(y), alternative = 'two.sided')
  
}



##########
# get pca function
##########

getPCA <- function(pca_data, 
                   column_name, 
                   name, 
                   gene_start, 
                   pca1,
                   pca2,
                   use_legend) 
{
  
  
  pca_data[, column_name] <- as.factor(pca_data[, column_name])
  
  
  # get features sites
  cg_sites <- colnames(pca_data)[gene_start:ncol(pca_data)]
  
  # subset by no NAs for column_name
  pca_data <- pca_data[!is.na(pca_data[, column_name]), ]
  
  stopifnot(!any(is.na(pca_data[, column_name])))
  
  # put column name with cg_sites 
  pca_data <- pca_data[ ,c(column_name, cg_sites)]
  
  # run pca
  data_length <- ncol(pca_data)
  pca <- prcomp(pca_data[,2:data_length])
  
  #fill in factors with colors 
  col_vec <- c('blue','red' , 'green', 'brown', 'orange', 'purple', 'lightblue', 
               'blueviolet', 'bisque', 'cyan', 'deeppink',
               'grey', 'yellow', 'bisque1', 'darkblue','darkred', 
               'darkgreen', 'darkorchid', 'gold', 'darkorange', 'coral',
               'greenyellow', 'bisque2')
  
  colors <- col_vec[pca_data[, column_name]]
  
  
  plot <- plot(pca$x[, pca1], 
               pca$x[, pca2],
               xlab = 'pca',
               ylab = 'pca',
               bty = 'n',
               cex = 1.3,
               main = name,
               pch = 16,
               col = adjustcolor(colors, alpha.f = 0.7)
  )
  abline(v = c(0,0),
         h = c(0,0))
  
  if(use_legend) {
    legend('topright',  
           legend = unique(pca_data[, column_name]), 
           col=unique(colors), 
           pch=16,  
           cex = 0.7)
  }
  
  return(plot)
}


##########
# remove outliers  4257 cases, 3391, 3392 controls
##########

removeOutlier <- function(data, 
                          cases, 
                          controls, 
                          val) {
  
  if (cases) {
    data <- data[data$ids != '3010',]
    
  }
  #controls outlier
  if (controls) {
    
    data <- data[data$ids != '3391',]
    data <- data[data$ids != '3392',]
  }
  
  if(val){
    data <- data[data$ids != '3540',]
    
  }
  
  
  return(data)
}





##########
# predict cancer
##########
predCancer <- function(training_dat, 
                       test_dat,
                       bh_features,
                       gender,
                       tech,
                       intersect_names) 
{
  
  # get intersection of bh features and real data
  intersected_feats <- intersect(intersect_names, colnames(training_dat))
  
  
  if(gender) {
    
    intersected_feats <- append('M', intersected_feats)
    intersected_feats <- append('F', intersected_feats)
  }
  
  if (tech) {
    
    intersected_feats <- append('a', intersected_feats)
    intersected_feats <- append('b', intersected_feats)
  }
  
  # # get 
  train_y <- ifelse(!grepl('Unaffected',training_dat$cancer_diagnosis_diagnoses), 1, 0)
  test_y <- ifelse(!grepl('Unaffected',test_dat$cancer_diagnosis_diagnoses), 1, 0)
  
  
  # get clinical data
  cg_start <- which(grepl('cg', colnames(test_dat)))[1]
  test_clin <- test_dat[, 1:(cg_start - 1)]
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
  
  N_CV_REPEATS = 2
  nfolds = 5
  
  ###### ENET
  # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
  # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
  elastic_net.cv_error = vector()
  elastic_net.cv_model = list()
  elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  
  # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
  # or if you have a high number fo N_CV_REPEATS
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
                                                , y =  train_y
                                                , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                , type.measure = type_measure
                                                , family = type_family
                                                , standardize = FALSE 
                                                , nfolds = nfolds 
                                                , nlambda = 10
                                                , parallel = TRUE
      )
      elastic_net.cv_error[alpha] = min(elastic_net.cv_model[[alpha]]$cvm)
    }
    elastic_net.cv_error # stores 9 errors    
  }
  
  if (N_CV_REPEATS == 1) {
    temp.cv_error_mean = temp.cv_error_matrix
  } else {
    temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean) # take the mean of the 5 iterations  
    # as your value for alpha
  }
  
  # stop if you did not recover error for any models 
  stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
  
  # get index of best alpha (lowest error) - alpha is values 0.1-0.9
  temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
  print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
  best_alpha <- elastic_net.ALPHA[temp.best_alpha_index]
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
                                     , y =  train_y
                                     , alpha = elastic_net.ALPHA[temp.best_alpha_index]
                                     , type.measure = type_measure
                                     , family = type_family
                                     , standardize=FALSE
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    # get optimal lambda - the tuning parameter for ridge and lasso
    # THIS IS IMPORTANT BECAUSE WHEN YOU TRAIN THE MODEL ON 100 SEPERATE VALUES OF LAMBDA
    # AND WHEN YOU TEST THE MODEL IT WILL RETURN PREDCITION FOR ALL THOSE VALUES (1-100). YOU NEED TO 
    # GRAB THE PREDICTION WITH SAME LAMBDA THAT YOU TRAINED ON. ITS ALL IN THE CODE, BUT JUST WANTED TO 
    # GIVE YOU REASONS
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    temp.loop_count = temp.loop_count + 1
    
    # set seed for next loop iteration
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    if (temp.loop_count > 10) {
      print("diverged")
      temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
      break
    }
  }# while loop ends 
  print(temp.non_zero_coeff)  
  
  model  = glmnet(x = as.matrix(training_dat)
                  , y =  train_y
                  ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  
  # get predictions with corresponding lambda.
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  
  return(list(importance, test_stats, clin_data))
  
  
}





# run bumphunter on two populations
# dat_1 <- m_train_cases
# dat_2 <- m_controls
# bump<- 'cancer'
# boot_num = 3
# m_beta_thresh = 0.5
bump_hunter <- function(dat_1,
                        dat_2,
                        bump,
                        boot_num,
                        thresh,
                        g_ranges) {
  
  # combine data
  dat <- rbind(dat_1, dat_2)
  
  if(bump == 'cancer') {
    
    dat$type <- ifelse(grepl('Unaffected', dat$cancer_diagnosis_diagnoses), 'controls', 'cases')
    ##########
    # get indicator and put into design matrix with intercept 1
    #########
    indicator_vector <- as.factor(dat$type)
    designMatrix <- cbind(rep(1, nrow(dat)), indicator_vector)
    designMatrix <- as.matrix(designMatrix)
  }
  # 
  # if(bump == 'age') {
  #   
  #   dat$type <- ifelse(dat$age_sample_collection > , 'controls', 'cases')
  #   ##########
  #   # get indicator and put into design matrix with intercept 1
  #   #########
  #   indicator_vector <- as.factor(dat$type)
  #   designMatrix <- cbind(rep(1, nrow(dat)), indicator_vector)
  #   designMatrix <- as.matrix(designMatrix)
  # }
  # 
  ##########
  # Get genetic locations
  ##########
  cg_start <- which(grepl('cg', colnames(dat)))[1]
  dat <- dat[, cg_start:(ncol(dat) -1) ]
  
  
  # transpose methylation to join with cg_locations to get genetic location vector.
  dat <- as.data.frame(t(dat), stringsAsFactors = F)
  
  # make probe a column in methyl
  dat$probe <- rownames(dat)
  rownames(dat) <- NULL
  
  # get probe column in granges 
  g_ranges$probe <- rownames(g_ranges)
  
  # inner join methyl and cg_locations by probe
  methyl_cg <- dplyr::inner_join(dat, g_ranges, by = 'probe')
  
  # get chr and pos vector 
  chr <- methyl_cg$seqnames
  pos <- methyl_cg$start
  
  # create beta matrix
  beta <- methyl_cg[, 1:(ncol(methyl_cg) - 6)]
  
  # make beta numeric 
  for (i in 1:ncol(beta)) {
    beta[,i] <- as.numeric(beta[,i])
    print(i)
  } 
  
  beta <- as.matrix(beta)
  
  ##########
  # Run bumphunter
  ##########
  
  # check dimensions 
  stopifnot(dim(beta)[2] == dim(designMatrix)[1])
  stopifnot(dim(beta)[1] == length(chr))
  stopifnot(dim(beta)[1] == length(pos))
  
  # set paramenters 
  DELTA_BETA_THRESH = thresh # DNAm difference threshold
  NUM_BOOTSTRAPS = boot_num  # number of randomizations
  
  # create tab list
  tab <- list()
  bump_hunter_results <- list()
  for (i in 1:length(DELTA_BETA_THRESH)) {
    tab[[i]] <- bumphunter(beta, 
                           designMatrix, 
                           chr = chr, 
                           pos = pos,
                           nullMethod = "bootstrap",
                           cutoff = DELTA_BETA_THRESH,
                           B = NUM_BOOTSTRAPS,
                           type = "M")
    
    bump_hunter_results[[i]] <- tab[[i]]$table
    bump_hunter_results[[i]]$run <- DELTA_BETA_THRESH[i]
  }
  
  bh_results <- do.call(rbind, bump_hunter_results)
  
  return(bh_results)
  
}


# training_dat <- beta_train_cases
# controls_dat <- beta_controls_full
# test_dat <-  beta_test_cases
# age_cutoff <- age_cutoff
# bh_features <- bh_features
# gender = T
# tech = T

run_enet_450_850 <- function(training_dat,
                             controls_dat,
                             test_dat,
                             age_cutoff,
                             gender, 
                             tech,
                             fam_cancer,
                             bh_features) {
  
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    
    intersected_feats <- append('M', intersected_feats)
    intersected_feats <- append('F', intersected_feats)
  }
  
  if (tech) {
    
    intersected_feats <- append('a', intersected_feats)
    intersected_feats <- append('b', intersected_feats)
  }
  
  if (fam_cancer){
    
    
    intersected_feats <- append('no', intersected_feats)
    intersected_feats <- append('yes', intersected_feats)
    
  }
  
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- ifelse(training_dat$age_diagnosis < age_cutoff, 1, 0)
  test_y <-  ifelse(test_dat$age_diagnosis < age_cutoff, 1, 0)
  
  # controls
  controls_y <-  ifelse(controls_dat$age_sample_collection < age_cutoff, 1, 0)
  
  # get clinical data
  cg_start <- which(grepl('cg', colnames(test_dat)))[1]
  test_clin <- test_dat[, 1:(cg_start - 1)]
  test_controls <- controls_dat[, 1:(cg_start - 1)]
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  
  # start elastic net tuning
  N_CV_REPEATS = 2
  nfolds = 3
  
  ###### ENET
  # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
  # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
  elastic_net.cv_error = vector()
  elastic_net.cv_model = list()
  elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  
  # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
  # or if you have a high number fo N_CV_REPEATS
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
                                                , y =  train_y
                                                , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                , type.measure = type_measure
                                                , family = type_family
                                                , standardize = FALSE 
                                                , nfolds = nfolds 
                                                , nlambda = 10
                                                , parallel = TRUE
      )
      elastic_net.cv_error[alpha] = min(elastic_net.cv_model[[alpha]]$cvm)
    }
    elastic_net.cv_error # stores 9 errors    
  }
  
  if (N_CV_REPEATS == 1) {
    temp.cv_error_mean = temp.cv_error_matrix
  } else {
    temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean) # take the mean of the 5 iterations  
    # as your value for alpha
  }
  
  # stop if you did not recover error for any models 
  stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
  
  # get index of best alpha (lowest error) - alpha is values 0.1-0.9
  temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
  # print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
  best_alpha <- elastic_net.ALPHA[temp.best_alpha_index]
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
                                     , y =  train_y
                                     , alpha = elastic_net.ALPHA[temp.best_alpha_index]
                                     , type.measure = type_measure
                                     , family = type_family
                                     , standardize=FALSE
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    # get optimal lambda - the tuning parameter for ridge and lasso
    # THIS IS IMPORTANT BECAUSE WHEN YOU TRAIN THE MODEL ON 100 SEPERATE VALUES OF LAMBDA
    # AND WHEN YOU TEST THE MODEL IT WILL RETURN PREDCITION FOR ALL THOSE VALUES (1-100). YOU NEED TO 
    # GRAB THE PREDICTION WITH SAME LAMBDA THAT YOU TRAINED ON. ITS ALL IN THE CODE, BUT JUST WANTED TO 
    # GIVE YOU REASONS
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    temp.loop_count = temp.loop_count + 1
    
    # set seed for next loop iteration
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    if (temp.loop_count > 10) {
      print("diverged")
      temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
      break
    }
  }# while loop ends 
  # print(temp.non_zero_coeff)  
  
  model  = glmnet(x = as.matrix(training_dat)
                  , y =  train_y
                  ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  # Predictions on test data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  
  # get predictions with corresponding lambda.
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  
  
  # Predictions on controls data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_con <- predict(model, 
                                       data.matrix(controls_dat),
                                       type = 'response')
  
  # get predictions with corresponding lambda.
  test.predictions_con <- temp_test.predictions_con[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  temp_dat_con <- as.data.frame(cbind(controls_age_pred = test.predictions_con, controls_age_label = controls_y, test_controls))
  
  ###########################################################################################
  return(list(temp_dat, temp_dat_con, model, elastic_net.cv_model$lambda.min, best_alpha))
  
}



run_enet_450_850_test <- function(cases_dat,
                                  controls_dat,
                                  age_cutoff,
                                  gender, 
                                  tech,
                                  base_change,
                                  exon_intron,
                                  bh_features) {
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(cases_dat))
  
  if(gender) {
    
    intersected_feats <- append('M', intersected_feats)
    intersected_feats <- append('F', intersected_feats)
  }
  
  if (tech) {
    
    intersected_feats <- append('a', intersected_feats)
    intersected_feats <- append('b', intersected_feats)
  }
  
  # if (base_change){
  #   
  #   
  #   intersected_feats <- append('none', intersected_feats)
  #   intersected_feats <- append('A', intersected_feats)
  #   intersected_feats <- append('C', intersected_feats)
  #   intersected_feats <- append('G', intersected_feats)
  #   intersected_feats <- append('T', intersected_feats)
  # }
  # 
  # if(exon_intron) {
  #   
  #   intersected_feats <- append('exon', intersected_feats)
  #   intersected_feats <- append('intron', intersected_feats)
  #   intersected_feats <- append('not_clear', intersected_feats)
  #   
  #   
  # }
  # # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  cases_y <- ifelse(cases_dat$age_diagnosis < age_cutoff, 1, 0)
  # controls
  controls_y <-  ifelse(controls_dat$age_sample_collection < age_cutoff, 1, 0)
  
  # get clinical data
  cg_start <- which(grepl('cg', colnames(cases_dat)))[1]
  cases_clin <- cases_dat[, 1:(cg_start - 1)]
  controls_clin <- controls_dat[, 1:(cg_start - 1)]
  
  # get bumphunter features
  cases_dat <- cases_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  
  # start elastic net tuning
  N_CV_REPEATS = 2
  nfolds = 3
  
  ###### ENET
  # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
  # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
  elastic_net.cv_error = vector()
  elastic_net.cv_model = list()
  elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  
  # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
  # or if you have a high number fo N_CV_REPEATS
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(cases_dat)
                                                , y =  cases_y
                                                , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                , type.measure = type_measure
                                                , family = type_family
                                                , standardize = FALSE 
                                                , nfolds = nfolds 
                                                , nlambda = 10
                                                , parallel = TRUE
      )
      elastic_net.cv_error[alpha] = min(elastic_net.cv_model[[alpha]]$cvm)
    }
    elastic_net.cv_error # stores 9 errors    
  }
  
  if (N_CV_REPEATS == 1) {
    temp.cv_error_mean = temp.cv_error_matrix
  } else {
    temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean) # take the mean of the 5 iterations  
    # as your value for alpha
  }
  
  # stop if you did not recover error for any models 
  stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
  
  # get index of best alpha (lowest error) - alpha is values 0.1-0.9
  temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
  # print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
  best_alpha <- elastic_net.ALPHA[temp.best_alpha_index]
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(cases_dat)
                                     , y =  cases_y
                                     , alpha = elastic_net.ALPHA[temp.best_alpha_index]
                                     , type.measure = type_measure
                                     , family = type_family
                                     , standardize=FALSE
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    # get optimal lambda - the tuning parameter for ridge and lasso
    # THIS IS IMPORTANT BECAUSE WHEN YOU TRAIN THE MODEL ON 100 SEPERATE VALUES OF LAMBDA
    # AND WHEN YOU TEST THE MODEL IT WILL RETURN PREDCITION FOR ALL THOSE VALUES (1-100). YOU NEED TO 
    # GRAB THE PREDICTION WITH SAME LAMBDA THAT YOU TRAINED ON. ITS ALL IN THE CODE, BUT JUST WANTED TO 
    # GIVE YOU REASONS
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    temp.loop_count = temp.loop_count + 1
    
    # set seed for next loop iteration
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    if (temp.loop_count > 10) {
      print("diverged")
      temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
      break
    }
  }# while loop ends 
  # print(temp.non_zero_coeff)  
  
  model  = glmnet(x = as.matrix(cases_dat)
                  , y =  cases_y
                  ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_con <- predict(model, 
                                       data.matrix(controls_dat),
                                       type = 'response')
  
  # get predictions with corresponding lambda.
  test.predictions_con <- temp_test.predictions_con[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  temp_dat_con <- as.data.frame(cbind(controls_age_pred = test.predictions_con, controls_age_label = controls_y, controls_clin))
  
  ###########################################################################################
  return(list(temp_dat_con, cases_clin))
  
}

# 
# training_dat = full_train_cases
# test_dat = full_test_cases
# age_cutoff = age_cutoff
# gender = gender
# tech = tech
# base_change = base_change
# exon_intron = exon_intron

run_coxreg <- function(training_dat, 
                       controls_dat,
                       test_dat,
                       random_forest,
                       rf_surv_fac,
                       rf_surv_con,
                       age_cutoff,
                       gender,
                       tech,
                       base_change,
                       exon_intron,
                       intersect_names) {
  
  
  if(random_forest) {
    
    
    test_clin <- test_dat[, 14:23]
    
    intersect_names <- intersect(intersect_names, colnames(training_dat))
    
    if(gender) {
      
      intersect_names <- append('M', intersect_names)
      intersect_names <- append('F', intersect_names)
    }
    
    if (tech) {
      
      intersect_names <- append('a', intersect_names)
      intersect_names <- append('b', intersect_names)
    }
    
    if (base_change){
      
      
      intersect_names <- append('none', intersect_names)
      intersect_names <- append('A', intersect_names)
      intersect_names <- append('C', intersect_names)
      intersect_names <- append('G', intersect_names)
      intersect_names <- append('T', intersect_names)
    }
    
    if(exon_intron) {
      
      intersect_names <- append('exon', intersect_names)
      intersect_names <- append('intron', intersect_names)
      intersect_names <- append('not_clear', intersect_names)
      
    }
    
    
    # # get training dat surv
    training_dat$time_to_event <- training_dat$age_diagnosis
    missing_ind <- is.na(training_dat$time_to_event)
    training_dat$time_to_event[missing_ind] <- training_dat$age_sample_collection[missing_ind]
    
    # # get test dat surv
    test_dat$time_to_event <- test_dat$age_diagnosis
    missing_ind <- is.na(test_dat$time_to_event)
    test_dat$time_to_event[missing_ind] <- test_dat$age_sample_collection[missing_ind]
    
    if(rf_surv_fac) {
      
      training_dat$cancer_status <- as.factor(ifelse(training_dat$cancer_diagnosis_diagnoses != 'Unaffected', 'a', 'b'))
      test_dat$cancer_status <- as.factor(ifelse(test_dat$cancer_diagnosis_diagnoses != 'Unaffected', 'a', 'b'))
      
      training_dat <- training_dat[, c('cancer_status', intersect_names)]
      test_dat <- test_dat[, c('cancer_status', intersect_names)]
      # get survival object
      surv_mod <- rfsrc(cancer_status ~., data = training_dat, ntree = 1000, tree.err = TRUE)
      
      # surv_mod <- rfsrc(Surv(time_to_event, cancer_status, type = 'right') ~., data = training_dat, ntree = 100, tree.err = TRUE)
      surv_pred <- predict(surv_mod, test_dat,importance = FALSE)
      temp_dat <- cbind(pred_y = surv_pred$predicted[,1], test_clin)
      
    }
    
    if(rf_surv_con){
      
      training_dat$cancer_status <- ifelse(training_dat$cancer_diagnosis_diagnoses != 'Unaffected', 1, 0)
      test_dat$cancer_status <- ifelse(test_dat$cancer_diagnosis_diagnoses != 'Unaffected', 1, 0)
      
      training_dat <- training_dat[, c('cancer_status', 'time_to_event',intersect_names)]
      test_dat <- test_dat[, c('cancer_status', 'time_to_event',intersect_names)]
      surv_mod <- rfsrc(Surv(time_to_event, cancer_status)~., data = training_dat, ntree = 1000, tree.err = TRUE)
      surv_pred <- predict(surv_mod, test_dat, importance = FALSE)
      temp_dat <- cbind(pred_y = surv_pred$predicted, test_clin)
      
    }
    
    return(temp_dat)
    
    
    
    
  } else {
    # get survival time as days to onset and days to sample collection in one column for training dat
    time_to_event <- training_dat$age_diagnosis
    missing_ind <- is.na(time_to_event)
    time_to_collection <- training_dat$age_sample_collection
    
    time_to_event[missing_ind] <- time_to_collection[missing_ind]
    
    
    # get cancer status 
    cancer_status <- ifelse(training_dat$cancer_diagnosis_diagnoses != 'Unaffected', 1, 0)
    
    test_clin <- test_dat[, 14:23]
    
    intersected_feats <- intersect(intersect_names, colnames(training_dat))
    
    if(gender) {
      
      intersected_feats <- append('M', intersected_feats)
      intersected_feats <- append('F', intersected_feats)
    }
    
    if (tech) {
      
      intersected_feats <- append('a', intersected_feats)
      intersected_feats <- append('b', intersected_feats)
    }
    
    if (base_change){
      
      
      intersected_feats <- append('none', intersected_feats)
      intersected_feats <- append('A', intersected_feats)
      intersected_feats <- append('C', intersected_feats)
      intersected_feats <- append('G', intersected_feats)
      intersected_feats <- append('T', intersected_feats)
    }
    
    if(exon_intron) {
      
      intersected_feats <- append('exon', intersected_feats)
      intersected_feats <- append('intron', intersected_feats)
      intersected_feats <- append('not_clear', intersected_feats)
      
    }
    
    # get bumphunter features
    training_dat <- training_dat[, intersected_feats]
    test_dat <- test_dat[, intersected_feats]
    
    # get survival object
    surv_outcome <- Surv(time_to_event,cancer_status, type = 'right')
    
    # fit coxph
    
    # start elastic net tuning
    N_CV_REPEATS = 2
    nfolds = 3
    
    ###### ENET
    # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
    # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
    elastic_net.cv_error = vector()
    elastic_net.cv_model = list()
    elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
    
    # set parameters for training model
    type_family <- 'cox'
    type_measure <- 'deviance'
    
    # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
    # or if you have a high number fo N_CV_REPEATS
    temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
      for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
      {      
        elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
                                                  , surv_outcome
                                                  , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                  , type.measure = type_measure
                                                  , family = type_family
                                                  , standardize = FALSE 
                                                  , nfolds = nfolds 
                                                  , nlambda = 10
                                                  , parallel = TRUE
        )
        elastic_net.cv_error[alpha] = min(elastic_net.cv_model[[alpha]]$cvm)
      }
      elastic_net.cv_error # stores 9 errors    
    }
    
    if (N_CV_REPEATS == 1) {
      temp.cv_error_mean = temp.cv_error_matrix
    } else {
      temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean) # take the mean of the 5 iterations  
      # as your value for alpha
    }
    
    # stop if you did not recover error for any models 
    stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
    
    # get index of best alpha (lowest error) - alpha is values 0.1-0.9
    temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
    # print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
    best_alpha <- elastic_net.ALPHA[temp.best_alpha_index]
    temp.non_zero_coeff = 0
    temp.loop_count = 0
    # loop runs initially because temp.non_zero coefficient <3 and then stops 
    # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
    # it they are never greater than 1, then the model does not converge. 
    while (temp.non_zero_coeff < 1) { 
      elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
                                       , surv_outcome
                                       , alpha = elastic_net.ALPHA[temp.best_alpha_index]
                                       , type.measure = type_measure
                                       , family = type_family
                                       , standardize=FALSE
                                       , nlambda = 100
                                       , nfolds = nfolds
                                       , parallel = TRUE
      )
      
      # get optimal lambda - the tuning parameter for ridge and lasso
      # THIS IS IMPORTANT BECAUSE WHEN YOU TRAIN THE MODEL ON 100 SEPERATE VALUES OF LAMBDA
      # AND WHEN YOU TEST THE MODEL IT WILL RETURN PREDCITION FOR ALL THOSE VALUES (1-100). YOU NEED TO 
      # GRAB THE PREDICTION WITH SAME LAMBDA THAT YOU TRAINED ON. ITS ALL IN THE CODE, BUT JUST WANTED TO 
      # GIVE YOU REASONS
      temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
      
      # # number of non zero coefficients at that lambda    
      temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
      temp.loop_count = temp.loop_count + 1
      
      # set seed for next loop iteration
      as.numeric(Sys.time())-> t 
      set.seed((t - floor(t)) * 1e8 -> seed) 
      if (temp.loop_count > 10) {
        print("diverged")
        temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
        break
      }
    }# while loop ends 
    # print(temp.non_zero_coeff)  
    
    model  = glmnet(x = as.matrix(training_dat)
                    , surv_outcome
                    ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                    ,standardize=FALSE
                    ,nlambda = 100
                    ,family = type_family)
    
    
    # This returns 100 prediction with 1-100 lambdas
    temp_test.predictions <- predict.coxnet(model, 
                                            data.matrix(test_dat),
                                            type = 'response')
    
    # get predictions with corresponding lambda.
    test.predictions <- temp_test.predictions[, temp.min_lambda_index]
    
    # combine predictions and real labels 
    temp_dat <- as.data.frame(cbind(test_pred = test.predictions,  test_clin))
    
    
    ###########################################################################################
    return(temp_dat)
  }
  
}


run_enet_all <- function(training_dat = train_cases,
                         controls_dat = beta_controls_mod,
                         valid_dat = beta_valid_mod,
                         test_dat = test_cases,
                         age_cutoff = age_cutoff,
                         gender = gender,
                         tech = tech,
                         base_change = base_change,
                         exon_intron = exon_intron,
                         bh_features = bh_features) {
  
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    
    intersected_feats <- append('M', intersected_feats)
    intersected_feats <- append('F', intersected_feats)
  }
  
  # if (tech) {
  #   
  #   intersected_feats <- append('a', intersected_feats)
  #   intersected_feats <- append('b', intersected_feats)
  # }
  
  # if (base_change){
  #   
  #   
  #   intersected_feats <- append('none', intersected_feats)
  #   intersected_feats <- append('A', intersected_feats)
  #   intersected_feats <- append('C', intersected_feats)
  #   intersected_feats <- append('G', intersected_feats)
  #   intersected_feats <- append('T', intersected_feats)
  # }
  # 
  # if(exon_intron) {
  #   
  #   intersected_feats <- append('exon', intersected_feats)
  #   intersected_feats <- append('intron', intersected_feats)
  #   intersected_feats <- append('not_clear', intersected_feats)
  #   
  #   
  # }
  # 
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- ifelse(training_dat$age_diagnosis < age_cutoff, 1, 0)
  test_y <-  ifelse(test_dat$age_diagnosis < age_cutoff, 1, 0)
  # controls
  controls_y <-  ifelse(controls_dat$age_sample_collection < age_cutoff, 1, 0)
  valid_y <-  ifelse(valid_dat$age_sample_collection < age_cutoff, 1, 0)
  
  
  # get clinical data
  cg_start <- which(grepl('cg', colnames(test_dat)))[1]
  test_clin <- test_dat[, 1:(cg_start - 1)]
  controls_clin <- controls_dat[, 1:(cg_start - 1)]
  valid_clin <- valid_dat[, 1:(cg_start - 1)]
  
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
  
  # start elastic net tuning
  N_CV_REPEATS = 2
  nfolds = 3
  
  ###### ENET
  # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
  # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
  elastic_net.cv_error = vector()
  elastic_net.cv_model = list()
  elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  
  # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
  # or if you have a high number fo N_CV_REPEATS
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
                                                , y =  train_y
                                                , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                , type.measure = type_measure
                                                , family = type_family
                                                , standardize = FALSE 
                                                , nfolds = nfolds 
                                                , nlambda = 10
                                                , parallel = TRUE
      )
      elastic_net.cv_error[alpha] = min(elastic_net.cv_model[[alpha]]$cvm)
    }
    elastic_net.cv_error # stores 9 errors    
  }
  
  if (N_CV_REPEATS == 1) {
    temp.cv_error_mean = temp.cv_error_matrix
  } else {
    temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean) # take the mean of the 5 iterations  
    # as your value for alpha
  }
  
  # stop if you did not recover error for any models 
  stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
  
  # get index of best alpha (lowest error) - alpha is values 0.1-0.9
  temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
  # print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
  best_alpha <- elastic_net.ALPHA[temp.best_alpha_index]
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
                                     , y =  train_y
                                     , alpha = elastic_net.ALPHA[temp.best_alpha_index]
                                     , type.measure = type_measure
                                     , family = type_family
                                     , standardize=FALSE
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    # get optimal lambda - the tuning parameter for ridge and lasso
    # THIS IS IMPORTANT BECAUSE WHEN YOU TRAIN THE MODEL ON 100 SEPERATE VALUES OF LAMBDA
    # AND WHEN YOU TEST THE MODEL IT WILL RETURN PREDCITION FOR ALL THOSE VALUES (1-100). YOU NEED TO 
    # GRAB THE PREDICTION WITH SAME LAMBDA THAT YOU TRAINED ON. ITS ALL IN THE CODE, BUT JUST WANTED TO 
    # GIVE YOU REASONS
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    temp.loop_count = temp.loop_count + 1
    
    # set seed for next loop iteration
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    if (temp.loop_count > 10) {
      print("diverged")
      temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
      break
    }
  }# while loop ends 
  # print(temp.non_zero_coeff)  
  
  model  = glmnet(x = as.matrix(training_dat)
                  , y =  train_y
                  ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  # Predictions on test data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  
  # get predictions with corresponding lambda.
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  
  
  # Predictions on controls data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_con <- predict(model, 
                                       data.matrix(controls_dat),
                                       type = 'response')
  
  # get predictions with corresponding lambda.
  test.predictions_con <- temp_test.predictions_con[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  temp_dat_con <- as.data.frame(cbind(controls_age_pred = test.predictions_con, controls_age_label = controls_y, controls_clin))
  
  # Predictions on controls data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_valid <- predict(model, 
                                         data.matrix(valid_dat),
                                         type = 'response')
  
  # get predictions with corresponding lambda.
  test.predictions_valid <- temp_test.predictions_valid[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  temp_dat_valid <- as.data.frame(cbind(valid_age_pred = test.predictions_valid, valid_age_label = valid_y, valid_clin))
  
  ###########################################################################################
  return(list(temp_dat, temp_dat_con, temp_dat_valid, model, elastic_net.cv_model$lambda.min, best_alpha))
  
}

### for one shot train and test on family overlaps
run_enet_family <- function(training_dat,
                            test_dat,
                            age_cutoff,
                            gender, 
                            tech,
                            bh_features) {
  
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    
    intersected_feats <- append('M', intersected_feats)
    intersected_feats <- append('F', intersected_feats)
  }
  
  if (tech) {
    
    intersected_feats <- append('a', intersected_feats)
    intersected_feats <- append('b', intersected_feats)
  }
  
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- ifelse(training_dat$age_diagnosis < age_cutoff, 1, 0)
  test_y <-  ifelse(test_dat$age_diagnosis < age_cutoff, 1, 0)
  
  # get clinical data
  cg_start <- which(grepl('cg', colnames(test_dat)))[1]
  test_clin <- test_dat[, 1:(cg_start - 1)]
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
  # start elastic net tuning
  N_CV_REPEATS = 2
  nfolds = 3
  
  ###### ENET
  # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
  # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
  elastic_net.cv_error = vector()
  elastic_net.cv_model = list()
  elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  
  # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
  # or if you have a high number fo N_CV_REPEATS
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
                                                , y =  train_y
                                                , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                , type.measure = type_measure
                                                , family = type_family
                                                , standardize = FALSE 
                                                , nfolds = nfolds 
                                                , nlambda = 10
                                                , parallel = TRUE
      )
      elastic_net.cv_error[alpha] = min(elastic_net.cv_model[[alpha]]$cvm)
    }
    elastic_net.cv_error # stores 9 errors    
  }
  
  if (N_CV_REPEATS == 1) {
    temp.cv_error_mean = temp.cv_error_matrix
  } else {
    temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean) # take the mean of the 5 iterations  
    # as your value for alpha
  }
  
  # stop if you did not recover error for any models 
  stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
  
  # get index of best alpha (lowest error) - alpha is values 0.1-0.9
  temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
  # print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
  best_alpha <- elastic_net.ALPHA[temp.best_alpha_index]
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
                                     , y =  train_y
                                     , alpha = elastic_net.ALPHA[temp.best_alpha_index]
                                     , type.measure = type_measure
                                     , family = type_family
                                     , standardize=FALSE
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    # get optimal lambda - the tuning parameter for ridge and lasso
    # THIS IS IMPORTANT BECAUSE WHEN YOU TRAIN THE MODEL ON 100 SEPERATE VALUES OF LAMBDA
    # AND WHEN YOU TEST THE MODEL IT WILL RETURN PREDCITION FOR ALL THOSE VALUES (1-100). YOU NEED TO 
    # GRAB THE PREDICTION WITH SAME LAMBDA THAT YOU TRAINED ON. ITS ALL IN THE CODE, BUT JUST WANTED TO 
    # GIVE YOU REASONS
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    temp.loop_count = temp.loop_count + 1
    
    # set seed for next loop iteration
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    if (temp.loop_count > 10) {
      print("diverged")
      temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
      break
    }
  }# while loop ends 
  # print(temp.non_zero_coeff)  
  
  model  = glmnet(x = as.matrix(training_dat)
                  , y =  train_y
                  ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  # Predictions on test data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  
  # get predictions with corresponding lambda.
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  
  
  
  ###########################################################################################
  return(temp_dat)
  
}



run_glmm_lasso <- function(training_dat,
                           controls_dat,
                           test_dat,
                           age_cutoff,
                           gender, 
                           tech,
                           bh_features) {
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    
    intersected_feats <- append('M', intersected_feats)
    intersected_feats <- append('F', intersected_feats)
  }
  
  if (tech) {
    
    intersected_feats <- append('a', intersected_feats)
    intersected_feats <- append('b', intersected_feats)
  }
  
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- ifelse(training_dat$age_diagnosis < age_cutoff, 1, 0)
  test_y <-  ifelse(test_dat$age_diagnosis < age_cutoff, 1, 0)
  
  # controls
  controls_y <-  ifelse(controls_dat$age_sample_collection < age_cutoff, 1, 0)
  
  # get clinical data
  cg_start <- which(grepl('cg', colnames(test_dat)))[1]
  test_clin <- test_dat[, 1:(cg_start - 1)]
  test_controls <- controls_dat[, 1:(cg_start - 1)]
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  
  # start elastic net tuning
  N_CV_REPEATS = 2
  nfolds = 3
  
  ###### ENET
  # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
  # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
  elastic_net.cv_error = vector()
  elastic_net.cv_model = list()
  elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  
  # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
  # or if you have a high number fo N_CV_REPEATS
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
                                                , y =  train_y
                                                , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                , type.measure = type_measure
                                                , family = type_family
                                                , standardize = FALSE 
                                                , nfolds = nfolds 
                                                , nlambda = 10
                                                , parallel = TRUE
      )
      elastic_net.cv_error[alpha] = min(elastic_net.cv_model[[alpha]]$cvm)
    }
    elastic_net.cv_error # stores 9 errors    
  }
  
  if (N_CV_REPEATS == 1) {
    temp.cv_error_mean = temp.cv_error_matrix
  } else {
    temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean) # take the mean of the 5 iterations  
    # as your value for alpha
  }
  
  # stop if you did not recover error for any models 
  stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
  
  # get index of best alpha (lowest error) - alpha is values 0.1-0.9
  temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
  # print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
  best_alpha <- elastic_net.ALPHA[temp.best_alpha_index]
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
                                     , y =  train_y
                                     , alpha = elastic_net.ALPHA[temp.best_alpha_index]
                                     , type.measure = type_measure
                                     , family = type_family
                                     , standardize=FALSE
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    temp.loop_count = temp.loop_count + 1
    
    # set seed for next loop iteration
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    if (temp.loop_count > 10) {
      print("diverged")
      temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
      break
    }
  }# while loop ends 
  # print(temp.non_zero_coeff)  
  
  model  = glmnet(x = as.matrix(training_dat)
                  , y =  train_y
                  ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  # Predictions on test data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  
  # get predictions with corresponding lambda.
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  
  
  # Predictions on controls data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_con <- predict(model, 
                                       data.matrix(controls_dat),
                                       type = 'response')
  
  # get predictions with corresponding lambda.
  test.predictions_con <- temp_test.predictions_con[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  temp_dat_con <- as.data.frame(cbind(controls_age_pred = test.predictions_con, controls_age_label = controls_y, test_controls))
  
  ###########################################################################################
  return(list(temp_dat, temp_dat_con, model, elastic_net.cv_model$lambda.min, best_alpha))
  
}



run_rf <- function(training_dat,
                   controls_dat,
                   test_dat,
                   age_cutoff,
                   gender, 
                   tech,
                   fam_cancer,
                   bh_features) {
  
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    
    intersected_feats <- append('M', intersected_feats)
    intersected_feats <- append('F', intersected_feats)
  }
  
  if (tech) {
    
    intersected_feats <- append('a', intersected_feats)
    intersected_feats <- append('b', intersected_feats)
  }
  
  if (fam_cancer){
    
    
    intersected_feats <- append('no', intersected_feats)
    intersected_feats <- append('yes', intersected_feats)
    
  }
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- ifelse(training_dat$age_diagnosis < age_cutoff, 'yes', 'no')
  test_y <-  ifelse(test_dat$age_diagnosis < age_cutoff, 'yes', 'no')
  
  # controls
  controls_y <-  ifelse(controls_dat$age_sample_collection < age_cutoff, 'yes', 'no')
  
  # get clinical data
  cg_start <- which(grepl('cg', colnames(test_dat)))[1]
  test_clin <- test_dat[, 1:(cg_start - 1)]
  test_controls <- controls_dat[, 1:(cg_start - 1)]
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  
  # determines how you train the model.
  NFOLDS <- 2
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = min(10, NFOLDS),
    classProbs = TRUE,
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = twoClassSummary
    
  )
  
  # mtry: Number of variables randomly sampled as candidates at each split.
  # ntree: Number of trees to grow.
  
  mtry <- sqrt(ncol(train_cases[,colnames(training_dat)]))
  tunegrid <- expand.grid(.mtry=mtry)
  
  model <- train(x = training_dat
                 , y = train_y
                 , method = "rf"
                 , trControl = fitControl
                 , tuneGrid = tunegrid
                 , importance = T
                 , verbose = FALSE)
  
  temp <- varImp(model)[[1]]
  importance <- cbind(rownames(temp), temp$X1)
  
  test.predictions <- predict(model
                              , newdata = model_data[-train_index, selected_features]
                              , type = "prob")
  
  ###########################################################################################
  return(list(temp_dat, temp_dat_con, model, elastic_net.cv_model$lambda.min, best_alpha))
  
}

