
##########
# load libraries
##########
library(IlluminaHumanMethylation450kmanifest)
library(tidyverse)
library(outliers)
library(preprocessCore)
library(bumphunter)
library(caret)
library(pROC)
library(doParallel)
library(e1071)
library(nnet)
library(glmnet)
library(PRROC)
library(ROCR)
library(survival)
library(wateRmelon)
library(sva)
library(randomForest)
library(RPMM)
library(c060)
library(RColorBrewer)
##################

registerDoParallel(6)




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
    Mset  <- preprocessSWAN(data)
  }
  if (preprocess == 'funnorm') {
    Mset <- preprocessFunnorm(data)
  }
  if (preprocess == 'noob') {
    Mset <- preprocessNoob(data, dyeMethod = 'single')
    dat <- wateRmelon::BMIQ(Mset)
    
   }
  # # map methyl set to genome (funnorm already does this)
  # Gset <- mapToGenome(Mset)
  # # get m values
  if(methyl_type == 'beta') {
    # dat_beta <- getBeta(Mset)
    dat <- getBeta(Mset)
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

# get log of data

get_log <- function(temp_dat, base_num){
  
  # store clinical data
  clin_dat <- temp_dat[, !grepl('^cg', names(temp_dat))]
  beta_dat <- temp_dat[, grepl('^cg', names(temp_dat))]
  negatives <- any(beta_dat < 0 )
  rm(temp_dat)
  
  if(negatives){
    get_positives <- apply(beta_dat, 2, function(x) all(x > 0))
    beta_dat <- beta_dat[, get_positives]
  }
  
  beta_log <- log(beta_dat, base = base_num)
  final_dat <- as.data.frame(cbind(clin_dat, beta_log))
  
  return(final_dat)
}


# # id function
# beta_data <- data_cases[1:10000,]
# id_map <- id_map_cases
process_rg_set_single <- function(beta_data, id_map, clin) {
  colnames(clin)[7] <- 'ids'
  beta_data <- findIds(beta_data, id_map)
  beta_data <- getIdName(beta_data)
  beta_data <- cleanIds(beta_data)
  beta_data <- beta_data[, !grepl('ch', colnames(beta_data))]
  # there are duplicate ids - loop through ids and avg duplicates 
  # dup_ids <- beta_data$ids[duplicated(beta_data$ids)]
  # col_list <- list()
  # sample_list <- list()
  # for(i in 1:length(dup_ids)){
  #   id_name <- dup_ids[i]
  #   sub_dat <- beta_data[beta_data$ids == id_name, ]
  #   for(j in 1:ncol(sub_dat)){
  #     temp_col <- sub_dat[, j]
  #     col_list[[j]] <- mean(temp_col)
  #     print(j)
  #   }
  #   sample_list[[i]] <- append(id_name, unlist(col_list))
  #  print(i)
  # }
  # 
  # remove duplicated ids
  # beta_data <- beta_data[!duplicated(beta_data$ids), ]
  # clin <- clin[!duplicated(clin$ids), ]
  beta_data <- dplyr::inner_join(clin, beta_data, by = 'ids')
  beta_data <- beta_data[!is.na(beta_data$tm_donor),]
  beta_data <- beta_data[!duplicated(beta_data$tm_donor),]
  cg_sites <- colnames(beta_data)[grepl('cg', colnames(beta_data))]
  beta_data <- beta_data[, c('ids', 
                             'p53_germline', 
                             'cancer_diagnosis_diagnoses', 
                             'age_diagnosis',
                             'age_sample_collection',
                             'gender',
                             'sentrix_id',
                             'family_name',
                             'tm_donor',
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


# remove WT
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}

clean_dat <- function(temp_dat, 
                      tech,
                      cases_or_controls, 
                      mut_or_wt){
  
  
  if(mut_or_wt == 'mut'){
    temp_dat <- remove_wild_type(temp_dat)
    
  }
  temp_dat <- temp_dat[!is.na(temp_dat$age_sample_collection),]
  
  if(cases_or_controls == 'cases') {
    temp_dat <- temp_dat[!is.na(temp_dat$age_diagnosis),]
    temp_dat <- temp_dat[!grepl('Unaffected', temp_dat$cancer_diagnosis_diagnoses),]
    if(tech == '850k'){
      temp_dat <- temp_dat[!temp_dat$tm_donor %in% tm_donor_450,]
    }
  } else {
    if(tech == '850k'){
      temp_dat <- temp_dat[grepl('Unaffected', temp_dat$cancer_diagnosis_diagnoses),]
    }
  }
  
  return(temp_dat)
}

# get shared controls and valid
get_shared_data <- function(temp_450, temp_850, cases_or_controls){
  
  if(cases_or_controls == 'cases'){
    temp_450 <- temp_450[!is.na(temp_450$age_diagnosis),]
    temp_450 <- temp_450[!is.na(temp_450$age_sample_collection),]
    temp_450 <- temp_450[!is.na(temp_450$gender),]
    temp_450 <- remove_wild_type(temp_450)
    temp_450 <- temp_450[!duplicated(temp_450$tm_donor),]
    
    temp_850 <- temp_850[!is.na(temp_850$age_diagnosis),]
    temp_850 <- temp_850[!is.na(temp_850$age_sample_collection),]
    temp_850 <- temp_850[!is.na(temp_850$gender),]
    temp_850 <- remove_wild_type(temp_850)
    temp_850 <- temp_850[!duplicated(temp_850$tm_donor),]
    
  } else {
    # 450
    temp_450 <- temp_450[!is.na(temp_450$age_sample_collection),]
    temp_450 <- temp_450[!is.na(temp_450$gender),]
    temp_450 <- remove_wild_type(temp_450)
    temp_450 <- temp_450[!duplicated(temp_450$tm_donor),]
    
    # 850
    temp_850 <- temp_850[!is.na(temp_850$age_sample_collection),]
    temp_850 <- temp_850[!is.na(temp_850$gender),]
    temp_850 <- remove_wild_type(temp_850)
    temp_850 <- temp_850[!duplicated(temp_850$tm_donor),]
    
  }
  # get unique ids
  tm_donor_450 <- unique(temp_450$tm_donor)
  tm_donor_850 <- unique(temp_850$tm_donor)
  
  # get interseciton
  intersect_tm_donors <- intersect(tm_donor_450, tm_donor_850)
  
  # subset both data sets by intersected tm 
  temp_450 <- temp_450[temp_450$tm_donor %in% intersect_tm_donors,]
  temp_850 <- temp_850[temp_850$tm_donor %in% intersect_tm_donors,]
  
  return(temp_850)
  
}

##########
# function for removing outlier from rgset
##########
# 
# rgSet <- rgControls
# id_map_dat <- id_map_con
# type = 'controls'
remove_outliers <- function(rgSet, id_map_dat, method, type) {
  # get outlier ids'3740|4122',
  outliers <- data.frame(ids =c('3391','3646','3301', '3701'),
                         batch = c('cases','cases', 'valid', 'valid')
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

subset_rg_set <- function(rg_set, 
                          keep_gender, 
                          keep_controls, 
                          keep_snps, 
                          get_island, 
                          get_type, 
                          get_chr, 
                          gene_probes) {
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



# combat function
run_combat <- function(temp_data, type) {
  temp_data <- as.data.frame(temp_data)
  
  if(type == 'age'){
    # create age fac for plotting
    temp_data$age_fac <- ifelse(temp_data$age_sample_collection > 0 &
                                 temp_data$age_sample_collection <= 36, 'first_age', 
                                ifelse(temp_data$age_sample_collection > 36 &
                                         temp_data$age_sample_collection <= 348, 'second_age','third_age'))
    # get batch
    batch_indicator <- as.character(temp_data$age_fac)
    batch_indicator <- as.factor(batch_indicator)
  }
  
  if(type == 'age2'){
    # create age fac for plotting
    temp_data$age_fac <- ifelse(temp_data$age_sample_collection > 0 &
                                  temp_data$age_sample_collection <= 12, 'first_age', 
                                ifelse(temp_data$age_sample_collection > 12 &
                                         temp_data$age_sample_collection <= 36, 'second_age',
                                       ifelse(temp_data$age_sample_collection > 36 &
                                                temp_data$age_sample_collection <= 140, 'fourth_age',
                                              ifelse(temp_data$age_sample_collection > 140 &
                                                       temp_data$age_sample_collection <= 211, 'fifth_age',
                                                     ifelse(temp_data$age_sample_collection > 211 &
                                                              temp_data$age_sample_collection <= 348, 'sixth_age', 'seventh_age')))))
    # get batch
    batch_indicator <- as.character(temp_data$age_fac)
    batch_indicator <- as.factor(batch_indicator)
  }
  
  if(type == 'tech'){
    # get tech variable variable back to categories 
    temp_data$tech <- ifelse(temp_data$tech == '450k', 'batch_1', 'batch_2')
    
    # get batch
    batch_indicator <- as.character(temp_data$tech)
    batch_indicator <- as.factor(batch_indicator)
  } 
  
  if(type == 'sentrix_id'){
    # get tech variable variable back to categories 

    # get batch
    batch_indicator <- as.character(temp_data$sentrix_id)
    batch_indicator <- as.factor(batch_indicator)
  }
  
  if(type == 'gender'){
    # get batch
    temp_data <- temp_data[!is.na(temp_data$gender),]
    batch_indicator <- as.character(temp_data$gender)
    batch_indicator <- as.factor(batch_indicator)
  }
  
  # put model ids in rownames and remove columns
  mat_combat <- as.matrix(temp_data[, grepl('^cg', names(temp_data))])
  rownames(mat_combat) <- NULL
  clin_combat <- as.data.frame(temp_data[, !grepl('^cg', names(temp_data))])
  # get features
  features <- colnames(mat_combat)
  mat_combat <- t(mat_combat)
  # get intercept
  modcombat <- model.matrix(~1, data = temp_data)
  combat <- ComBat(dat = mat_combat, batch = batch_indicator, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  # transpose and add back columns
  final_dat <- as.data.frame(t(combat))
  final_dat <- as.data.frame(cbind(clin_combat, final_dat))
  rownames(final_dat) <- NULL
  return(final_dat)
}


##########
# get pca function
##########
# pca_data = con_450
# column_name = 'tech'
# show_variance = FALSE
# pc_x = 1
# pc_y = 2
# main_title = 'PC:cases bet'

get_pca <- function(pca_data, 
                    column_name,
                    show_variance,
                    pc_x,
                    pc_y,
                    main_title) {
  
  
  pca_data <- as.data.frame(pca_data)
  pca_data[, column_name] <- as.factor(pca_data[, column_name])
  
  # get other clinical data
  column_names <- names(pca_data)[!grepl('^cg', names(pca_data))]
    
  # remove column_name
  column_names <- column_names[!column_names %in% column_name]
  
  # get features sites
  cg_sites <- names(pca_data)[grepl('^cg', names(pca_data))]
  # subset by no NAs for column_name
  pca_data <- pca_data[!is.na(pca_data[, c(column_name)]), ]
  
  stopifnot(!any(is.na(pca_data[, column_name])))
  
  # put column name with cg_sites 
  pca_data <- pca_data[ ,c(column_name,column_names, cg_sites)]
  
  # run pca
  data_length <- ncol(pca_data)
  cg_start <- which(grepl('^cg', names(pca_data)))[1]
  pca <- prcomp(pca_data[,cg_start:data_length])
  
  if(show_variance){
    # plot lambda
    return(plot(pca, xlim = c(0,10), type = 'l', main = main_title))
  } else {
    # get pca dataframe with results and factor to color
    pca_results <- data.frame(pca$x[, 1:10], column_name = pca_data[, column_name],
                              pca_data[, column_names])
    
    # get actual PC
    pca_results <- pca_results[ c(pc_x,pc_y, 11:ncol(pca_results))]
    real_x_axis <- names(pca_results)[1]
    real_y_axis <- names(pca_results)[2]
    
    # now rename first two to var1, var2 so it can still plot
    names(pca_results)[1:2] <- c('V1', 'V2')
    
    # get color
    cols <- colorRampPalette(brewer.pal(n = 9, 'Set1'))(length(unique(pca_results$column_name)))
    
    plot <- 
      ggplot(pca_results, 
             aes(V1, V2, 
                 col = column_name)) +
      geom_point(size = 3, 
                 alpha = 0.7) +
      xlab(real_x_axis) + 
      ylab(real_y_axis) +
      scale_color_manual(name = '',
                         values = cols) + 
      ggtitle(main_title) +
      geom_hline(yintercept= 0, linetype="dashed", 
                 color = "grey", size=1) +
      geom_vline(xintercept=0, linetype="dashed", 
                 color = "grey", size=1) +
      theme_minimal(base_size = 18)
    
    # plot <- ggplotly(plot)
    
    return(plot)
  }
  
}



##########b
# predict cancer
##########
# 
# train_dat = train_full
# test_dat = test_full
# age_cutoff = 72
# gender = gender
# tech = tech
# fam_num = fam_num
# fam_ratio = fam_ratio
# bh_features = remaining_features

pred_cancer_enet <- function(train_dat, 
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
                             fam_num,
                             fam_ratio,
                             bh_features) {
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(train_dat))
  
  if(gender) {
    intersected_feats <- c('M', intersected_feats)
    intersected_feats <- c('F', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('a', intersected_feats)
    intersected_feats <- c('b', intersected_feats)
  }
  if (fam_num){
    intersected_feats <- c('fam_num_cancer', intersected_feats)
  }
  if (fam_ratio){
    intersected_feats <- c('fam_cancer_ratio', intersected_feats)
  }
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- ifelse(!grepl('Unaffected', train_dat$cancer_diagnosis_diagnoses), 1, 0)
  test_y <- ifelse(!grepl('Unaffected', test_dat$cancer_diagnosis_diagnoses), 1, 0)
  

  # get clinical data
  train_clin <- train_dat[, !grepl('^cg', colnames(train_dat))]
  test_clin <- test_dat[, !grepl('^cg', colnames(test_dat))]
  
  # get model data
  train_dat <- train_dat[, intersected_feats]
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
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(train_dat)
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
    elastic_net.cv_model = cv.glmnet(x = as.matrix(train_dat)
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
  
  model  = glmnet(x = as.matrix(train_dat)
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
  
  return(list(model, temp_dat, best_alpha))
  
  
}



##########
# predict cancer
##########

# train_dat = train_full
# test_dat = test_full
# age_cutoff = 72
# gender = gender
# tech = tech
# fam_num = fam_num
# fam_ratio = fam_ratio
# bh_features = remaining_features

pred_cancer_rf <- function(train_dat, 
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
                             fam_num,
                             fam_ratio,
                             bh_features) {
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(train_dat))
  
  if(gender) {
    intersected_feats <- c('M', intersected_feats)
    intersected_feats <- c('F', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('a', intersected_feats)
    intersected_feats <- c('b', intersected_feats)
  }
  if (fam_num){
    intersected_feats <- c('fam_num_cancer', intersected_feats)
  }
  if (fam_ratio){
    intersected_feats <- c('fam_cancer_ratio', intersected_feats)
  }
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- ifelse(!grepl('Unaffected', train_dat$cancer_diagnosis_diagnoses), 'yes', 'no')
  test_y <- ifelse(!grepl('Unaffected', test_dat$cancer_diagnosis_diagnoses), 'yes', 'no')
  
  
  # get clinical data
  train_clin <- train_dat[, !grepl('^cg', colnames(train_dat))]
  test_clin <- test_dat[, !grepl('^cg', colnames(test_dat))]
  
  # get model data
  train_dat <- train_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
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
  
  mtry <- sqrt(ncol(train_dat[,colnames(train_dat)]))
  tunegrid <- expand.grid(.mtry=mtry)
  
  model <- train(x = train_dat
                 , y = train_y
                 , metric = 'ROC'
                 , method = "rf"
                 , trControl = fitControl
                 , tuneGrid = tunegrid
                 , importance = T
                 , verbose = FALSE)
  
  temp <- varImp(model)[[1]]
  importance <- cbind(rownames(temp), temp$X1)
  
  # Predictions on test data
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions <- predict(model, 
                              data.matrix(test_dat),
                              type = 'prob')
  
  
  # combine predictions and real labels 
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  
  return(list(model, temp_dat, importance))
  
  
}

get_age_removal <- function(temp_dat){
  clin_names <- names(temp_dat)[1:11]
  cpg_names <- names(temp_dat)[12:ncol(temp_dat)]
  # remove age_probes 
  intersect_probes <- intersect(cpg_names, age_probes)
  new_probes <- cpg_names[!cpg_names %in% intersect_probes]
  temp_dat <- temp_dat[, c(clin_names, new_probes)]
  
  return(temp_dat)
}

# subset all dady by bh_feats
join_new_features <- function(temp_dat, new_features){
  # new_features$start <- as.character(new_features$start)
  # new_features$end <- as.character(new_features$end)
  clin_name <- colnames(temp_dat)[!grepl('^cg', colnames(temp_dat))]
  colnames(new_features)[1] <- 'chr'
  keep_features <- inner_join(new_features, g_ranges)$probe
  final_dat <- temp_dat[c(clin_name, keep_features)]
  return(final_dat)
  
}

remove_cancer_feats <- function(temp_dat, bh_feats){
  clin_names <- names(temp_dat)[!grepl('^cg', names(temp_dat))]
  final_dat <- temp_dat[, c(clin_names, bh_feats)]
  return(final_dat)
}
# dat_1 = con_wt
# dat_2 = con_mut
# bump = 'lfs'
# boot_num = 5 
# beta_thresh = beta_thresh
# methyl_type = methyl_type
# g_ranges = g_ranges

bump_hunter <- function(dat_1,
                        dat_2,
                        wild_type,
                        bump,
                        boot_num,
                        methyl_type,
                        beta_thresh,
                        g_ranges) {
  
  
  
  if(bump == 'lfs'){
    dat_wt <- dat_1
    # get mutual cgs
    wt_names <- colnames(dat_wt)[grepl('^cg', colnames(dat_wt))]
    dat_names <- colnames(dat_2)[grepl('^cg', colnames(dat_2))]
    wt_intersect <- intersect(wt_names, dat_names)
    
    # stor clinical data
    clin_wt_names <- colnames(dat_wt)[!grepl('^cg', colnames(dat_wt))]
    clin_con_names <- colnames(dat_2)[!grepl('^cg', colnames(dat_2))]
    
    # get data
    dat_wt <- dat_wt[, c(clin_wt_names, wt_intersect)]
    dat_2 <- dat_2[, c(clin_con_names, wt_intersect)]
    
    # combine data
    dat <- rbind(dat_wt, 
                 dat_2)
    
    # remove NAs 
    dat <- dat[!is.na(dat$p53_germline),]
    
    # get indicator (bump) vector
    dat$type <- dat$p53_germline
    indicator_vector <- as.factor(dat$type)
    designMatrix <- cbind(rep(1, nrow(dat)), indicator_vector)
    designMatrix <- as.matrix(designMatrix)
  } else if(bump == 'cancer') {
    # combine data
    dat <- rbind(dat_1, dat_2)
    dat$type <- ifelse(grepl('Unaffected', dat$cancer_diagnosis_diagnoses), 'controls', 'cases')
    indicator_vector <- as.factor(dat$type)
    designMatrix <- cbind(rep(1, nrow(dat)), indicator_vector)
    designMatrix <- as.matrix(designMatrix)
  }
  
  # get probe data 
  dat <- dat[, grepl('^cg', colnames(dat))]
  # features <- names(dat)[grepl('^cg', colnames(dat))]

  # transpose methylation to join with cg_locations to get genetic location vector.
  dat <- as.data.frame(t(dat), stringsAsFactors = F)
  
  # make probe a column in methyl
  dat$probe <- rownames(dat)


  # get probe column in granges 
  g_ranges$probe <- rownames(g_ranges)
  
  # inner join methyl and cg_locations by probe
  methyl_cg <- dplyr::inner_join(dat, g_ranges, by = 'probe')
  
  # get chr and pos vector 
  chr <- methyl_cg$chr
  pos <- as.numeric(methyl_cg$start)
  
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
  DELTA_BETA_THRESH = beta_thresh # DNAm difference threshold
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
                           type = methyl_type)
    
    bump_hunter_results[[i]] <- tab[[i]]$table
    bump_hunter_results[[i]]$run <- DELTA_BETA_THRESH[i]
  }
  
  bh_results <- do.call(rbind, bump_hunter_results)
  
  return(bh_results)
  
}



run_enet_all_test <- function(training_dat,
                              test_dat,
                              controls_dat,
                              valid_dat,
                              age_cutoff,
                              gender,
                              tech,
                              bh_features,
                              standardize) {
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    intersected_feats <- c('Female', 'Male', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('batch_1', 'batch_2', intersected_feats)
  }
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- as.factor(ifelse(training_dat$age_diagnosis < age_cutoff, 'positive', 'negative'))
  test_y <-  as.factor(ifelse(test_dat$age_diagnosis < age_cutoff, 'positive', 'negative'))
  
  
  # get clinical data
  test_clin <- test_dat[, !grepl('^cg', colnames(test_dat))]
  
  # if(use_offset){
  #   offsetted_train <- as.numeric(training_dat$age_sample_collection)
  #   offsetted_test <- as.numeric(test_dat$age_sample_collection)
  #   
  # } else {
  #   offsetted <- NULL
  # }
  
  # get model data
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
  
  # start elastic net tuning
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
                                                , standardize = standardize
                                                , nfolds = nfolds 
                                                , nlambda = 10
                                                , parallel = TRUE
                                                , offset = NULL
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
                                     , standardize=standardize
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    lambda_value <- elastic_net.cv_model$lambda.min
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
                  ,standardize= standardize
                  ,nlambda = 100
                  ,family = type_family)
  
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  library(DescTools)
  new_lambda <- Closest(model$lambda, lambda_value)
  lambda_value_2 <- model$lambda[model$lambda == new_lambda]
  
  # get predictions with corresponding lambda.
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  test_results <- as.data.frame(cbind(preds = test.predictions, real = test_y, test_clin))
  test_results$pred_class <- as.factor(ifelse(test_results$preds > .5, 'positive', 'negative'))
  
  # relevel both factore
  test_results$pred_class <- factor(test_results$pred_class, c('positive', 'negative'))
  test_results$real <- factor(test_results$real, c('positive', 'negative'))
  
  test_results$accuracy <- caret::confusionMatrix(table(test_results$pred_class, test_results$real))$overall[1]
  test_results$alpha <- best_alpha
  test_results$lambda <- elastic_net.cv_model$lambda.min
  test_results$non_zero <- temp.non_zero_coeff
  test_results$lambda_value <- lambda_value
  test_results$lambda_value_model <- lambda_value_2
  
  test_results$tot_probes <- ncol(training_dat)
  
  return(test_results)
  
}

run_enet_test <- function(training_dat,
                              test_dat,
                              controls_dat,
                              valid_dat,
                              age_cutoff,
                              alpha_num,
                              gender,
                              tech,
                              bh_features,
                              standardize) {
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    intersected_feats <- c('Female', 'Male', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('450k', '850k', intersected_feats)
  }
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- as.factor(ifelse(training_dat$age_diagnosis < age_cutoff, 'positive', 'negative'))
  test_y <-  as.factor(ifelse(test_dat$age_diagnosis < age_cutoff, 'positive', 'negative'))
  train_y <- factor(train_y, levels = c('positive', 'negative'))
  test_y <- factor(test_y, levels = c('positive', 'negative'))
  
  sub_con <- controls_dat[controls_dat$tech == '450k',]
  con_y <- as.factor(ifelse(sub_con$age_sample_collection < age_cutoff, 'positive', 'negative'))
  con_y <- factor(con_y, levels = c('positive', 'negative'))
  
  # get clinical data
  test_clin <- test_dat[, !grepl('^cg', colnames(test_dat))]
  con_clin <-sub_con[, !grepl('^cg', colnames(sub_con))]
  
  
  # if(use_offset){
  #   offsetted_train <- as.numeric(training_dat$age_sample_collection)
  #   offsetted_test <- as.numeric(test_dat$age_sample_collection)
  #   
  # } else {
  #   offsetted <- NULL
  # }
  
  # get model data
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  sub_con <- sub_con[, intersected_feats]
  
  
  # start elastic net tuning
  N_CV_REPEATS = 2
  nfolds = 5
  
  
  # print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
  best_alpha <- alpha_num
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  type_measure = 'auc'
  type_family = 'binomial'
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
                                     , y =  train_y
                                     , alpha = best_alpha
                                     , type.measure = type_measure
                                     , family = type_family
                                     , standardize=standardize
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    lambda_value <- elastic_net.cv_model$lambda.min
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
                  ,alpha = alpha_num
                  ,standardize= standardize
                  ,nlambda = 100
                  ,family = type_family)
  
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  # get predictions with corresponding lambda.
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  test_results <- as.data.frame(cbind(preds = test.predictions, real = test_y, test_clin))
  
  
  
  # controls
  # This returns 100 prediction with 1-100 lambdas
  temp_con.predictions <- predict(model, 
                                   data.matrix(sub_con),
                                   type = 'response')
  
  # get predictions with corresponding lambda.
  con.predictions <- temp_con.predictions[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  test_results_con <- as.data.frame(cbind(preds = con.predictions, real = con_y, con_clin))
  
  # get optimal threshold
  # creat real label
  
  
  
  test_results$alpha <- best_alpha
  test_results$lambda <- elastic_net.cv_model$lambda.min
  test_results$non_zero <- temp.non_zero_coeff
  test_results$lambda_value <- lambda_value

  test_results$tot_probes <- ncol(training_dat)
  
  test_results_con$alpha <- best_alpha
  test_results_con$lambda <- elastic_net.cv_model$lambda.min
  test_results_con$non_zero <- temp.non_zero_coeff
  test_results_con$lambda_value <- lambda_value
  
  test_results_con$tot_probes <- ncol(training_dat)
  
  return(list(test_results, test_results_con))
  
}


run_lasso_all_test <- function(training_dat,
                              test_dat,
                              controls_dat,
                              valid_dat,
                              age_cutoff,
                              gender,
                              tech,
                              bh_features,
                              standardize) {
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    intersected_feats <- c('Female', 'Male', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('batch_1', 'batch_2', intersected_feats)
  }
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- as.factor(ifelse(training_dat$age_diagnosis < age_cutoff, 'positive', 'negative'))
  test_y <-  as.factor(ifelse(test_dat$age_diagnosis < age_cutoff, 'positive', 'negative'))
  
  nfolds = 5
  # get clinical data
  test_clin <- test_dat[, !grepl('^cg', colnames(test_dat))]
  
  # if(use_offset){
  #   offsetted_train <- as.numeric(training_dat$age_sample_collection)
  #   offsetted_test <- as.numeric(test_dat$age_sample_collection)
  #   
  # } else {
  #   offsetted <- NULL
  # }
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  
  # get model data
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
  
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
                                     , y =  train_y
                                     , alpha = 1
                                     , type.measure = type_measure
                                     , family = type_family
                                     , standardize=standardize
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    lambda_value <- elastic_net.cv_model$lambda.min
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
                  ,alpha = 1
                  ,standardize= standardize
                  ,nlambda = 100
                  ,family = type_family)
  
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  library(DescTools)
  new_lambda <- Closest(model$lambda, lambda_value)
  lambda_value_2 <- model$lambda[model$lambda == new_lambda]
  
  # get predictions with corresponding lambda.
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  test_results <- as.data.frame(cbind(preds = test.predictions, real = test_y, test_clin))
  test_results$pred_class <- as.factor(ifelse(test_results$preds > .5, 'positive', 'negative'))
  
  # relevel both factore
  test_results$pred_class <- factor(test_results$pred_class, c('positive', 'negative'))
  test_results$real <- factor(test_results$real, c('positive', 'negative'))
  
  test_results$accuracy <- caret::confusionMatrix(table(test_results$pred_class, test_results$real))$overall[1]
  test_results$lambda <- elastic_net.cv_model$lambda.min
  test_results$non_zero <- temp.non_zero_coeff
  test_results$lambda_value <- lambda_value
  test_results$lambda_value_model <- lambda_value_2
  
  test_results$tot_probes <- ncol(training_dat)
  
  return(test_results)
  
}

# temp_450 <- shared_cases_450
# temp_850 <- shared_cases_850
# full_data <- cases_valid
linearTransform <- function (temp_450, 
                             temp_850, 
                             full_data,
                             pred_direction) {
  
  probe_model <- list()
  probe_850_result <- list()
  
  for (i in 12:ncol(full_data)) {
    
    probe_450 <- as.data.frame(temp_450[, i])
    probe_850<- as.data.frame(temp_850[, i])
    model_data <- data.frame(probe_850= probe_850, probe_450 = probe_450)
    names(model_data) <- c('probe_850', 'probe_450')
    model_data$probe_850 <- as.numeric(as.character(model_data$probe_850))
    model_data$probe_450 <- as.numeric(as.character(model_data$probe_450))
    
    probe_model[[i]] <- lm(probe_450 ~ probe_850, data = model_data)
    probe_full <- as.numeric(full_data[, i])
    model_data_new <- data.frame(probe_full = probe_full)
    names(model_data_new) <- 'probe_850'
    probe_850_result[[i]] <- predict(probe_model[[i]], 
                                     newdata = model_data_new, 
                                     type = 'response')
    
    print(i) 
  }
  
  # transpose results
  temp <- do.call(rbind, probe_850_result)
  transform_850 <- as.data.frame(t(temp))
  
  # add cg sites
  colnames(transform_850) <- colnames(full_data)[12:ncol(full_data)]
  transform_850 <- as.data.frame(transform_850)
  
  # add clinical variables
  transform_850 <- as.data.frame(cbind(id = full_data$id, 
                                       p53_germline = full_data$p53_germline, 
                                       cancer_diagnosis_diagnoses = full_data$cancer_diagnosis_diagnoses, 
                                       age_diagnosis = full_data$age_diagnosis, 
                                       age_sample_collection = full_data$age_sample_collection,
                                       gender = full_data$gender, 
                                       sentrix_id = full_data$sentrix_id, 
                                       family_name = full_data$family_name,
                                       tm_donor = full_data$tm_donor,
                                       transform_850))
  
  # make numeric
  transform_850[, 12:ncol(transform_850)] <- apply(transform_850[, 12:ncol(transform_850)], 
                                                   2, 
                                                   function(x) as.numeric(x))
  
  
  return(transform_850)
  
}


# training_dat = train_cases
# test_dat = test_cases
# controls_dat = controls_full
# valid_dat = cases_850
# age_cutoff = age_cutoff
# gender = gender
# tech = tech
# bh_features = bh_features

run_rf_all_test <- function(training_dat,
                            test_dat,
                            controls_dat,
                            valid_dat,
                            age_cutoff,
                            gender,
                            tech,
                            bh_features) {
  
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    intersected_feats <- c('Female', 'Male', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('batch_1', 'batch_2', intersected_feats)
  }
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- factor(as.factor(ifelse(training_dat$age_diagnosis < age_cutoff, 'positive', 'negative')), levels = c('positive', 'negative'))
  test_y <- factor( as.factor(ifelse(test_dat$age_diagnosis < age_cutoff, 'positive', 'negative')), levels = c('positive', 'negative'))
  con_y <-  factor(as.factor(ifelse(controls_dat$age_sample_collection < age_cutoff, 'positive', 'negative')), levels = c('positive', 'negative'))
  valid_y <- factor( as.factor(ifelse(valid_dat$age_diagnosis < age_cutoff, 'positive', 'negative')), levels = c('positive', 'negative'))
  
  
  # get clinical data
  test_clin <- test_dat[, !grepl('^cg', colnames(test_dat))]
  con_clin <- controls_dat[, !grepl('^cg', colnames(controls_dat))]
  valid_clin <- valid_dat[, !grepl('^cg', colnames(valid_dat))]
  
  # if(use_offset){
  #   offsetted_train <- as.numeric(training_dat$age_sample_collection)
  #   offsetted_test <- as.numeric(test_dat$age_sample_collection)
  #   
  # } else {
  #   offsetted <- NULL
  # }
  
  # get model data
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  controls_dat<- controls_dat[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
  
  # determines how you train the model.
  NFOLDS <- 5
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
  
  mtry <- sqrt(ncol(training_dat[,colnames(training_dat)]))
  tunegrid <- expand.grid(.mtry=mtry)
  
  model <- train(x = training_dat
                 , y = train_y
                 , metric = 'ROC'
                 , method = "rf"
                 , trControl = fitControl
                 , tuneGrid = tunegrid
                 , importance = T
                 , verbose = FALSE)
  
  temp <- varImp(model)[[1]]
  importance <- cbind(rownames(temp), temp$positive)
  importance <- as.data.frame(importance)
  importance$V2 <- round(as.numeric(as.character(importance$V2)), 2)
  names(importance) <- c('probe', 'score')
  
  # Predictions on test data
  # This returns 100 prediction with 1-100 lambdas
  test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'prob')
  con.predictions <- predict(model, 
                             data.matrix(controls_dat),
                             type = 'prob')
  
  valid.predictions <- predict(model, 
                             data.matrix(valid_dat),
                             type = 'prob')
  


  
  # combine predictions and real labels 
  test_results <- as.data.frame(cbind(test.predictions,real = test_y, test_clin))
  
  test_results$pred_class <- as.factor(ifelse(test_results$positive > .5, 'positive', 'negative'))
  test_results$pred_class_con <- as.factor(ifelse(test_results$positive > .5, 'positive', 'negative'))
  test_results$pred_class_val<- as.factor(ifelse(test_results$positive > .5, 'positive', 'negative'))
  
  # relevel both factore
  test_results$pred_class <- factor(test_results$pred_class, c('positive', 'negative'))
  test_results$real <- factor(test_results$real, c('positive', 'negative'))
  
  test_results$tot_probes <- ncol(training_dat)
  
  # controls
  test_results_con <- as.data.frame(cbind(con.predictions, real_con_y = con_y, con_clin))
  test_results_con$pred_class <- as.factor(ifelse(test_results_con$positive > .5, 'positive', 'negative'))

  # relevel both factore
  test_results_con$pred_class <- factor(test_results_con$pred_class, c('positive', 'negative'))
  test_results_con$real <- factor(test_results_con$real, c('positive', 'negative'))
  
  test_results_con$tot_probes <- ncol(controls_dat)
  
  
  # valid
  test_results_val <- as.data.frame(cbind(valid.predictions, real_valid_y = valid_y, valid_clin))
  
  test_results_val$pred_class <- as.factor(ifelse(test_results_val$positive > .5, 'positive', 'negative'))
  
  # relevel both factore
  test_results_val$pred_class <- factor(test_results_val$pred_class, c('positive', 'negative'))
  test_results_val$real <- factor(test_results_val$real, c('positive', 'negative'))
  
  test_results_val$tot_probes <- ncol(valid_dat)
  
  
  
  return(list(test_results, importance, test_results_con, test_results_val))
  

}



run_rf_test <- function(training_dat,
                            test_dat,
                            controls_dat,
                            valid_dat,
                            age_cutoff,
                            gender,
                            tech,
                            bh_features) {
  
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    intersected_feats <- c('Female', 'Male', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('450k', '850k', intersected_feats)
  }
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- as.factor(ifelse(training_dat$age < age_cutoff, 'positive', 'negative'))
  test_y <-  as.factor(ifelse(test_dat$age < age_cutoff, 'positive', 'negative'))
  train_y <- factor(train_y, levels = c('positive', 'negative'))
  test_y <- factor(test_y, levels = c('positive', 'negative'))
  
  sub_con <- controls_dat[controls_dat$tech == '850k',]
  con_y <- as.factor(ifelse(sub_con$age_sample_collection < age_cutoff, 'positive', 'negative'))
  con_y <- factor(con_y, levels = c('positive', 'negative'))
  
  # get clinical data
  test_clin <- test_dat[, !grepl('^cg', colnames(test_dat))]
  con_clin <-sub_con[, !grepl('^cg', colnames(sub_con))]
  
  
  # if(use_offset){
  #   offsetted_train <- as.numeric(training_dat$age_sample_collection)
  #   offsetted_test <- as.numeric(test_dat$age_sample_collection)
  #   
  # } else {
  #   offsetted <- NULL
  # }
  
  # get model data
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  sub_con <- sub_con[, intersected_feats]
  
  # determines how you train the model.
  NFOLDS <- 5
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
  
  mtry <- sqrt(ncol(training_dat[,colnames(training_dat)]))
  tunegrid <- expand.grid(.mtry=mtry)
  
  model <- train(x = training_dat
                 , y = train_y
                 , metric = 'ROC'
                 , method = "rf"
                 , trControl = fitControl
                 , tuneGrid = tunegrid
                 , importance = T
                 , verbose = FALSE)
  
  temp <- varImp(model)[[1]]
  importance <- cbind(rownames(temp), temp$positive)
  importance <- as.data.frame(importance)
  importance$V2 <- round(as.numeric(as.character(importance$V2)), 2)
  names(importance) <- c('probe', 'score')
  
  # Predictions on test data
  # This returns 100 prediction with 1-100 lambdas
  test.predictions <- predict(model, 
                              data.matrix(test_dat),
                              type = 'prob')
  
  
  
  
  # combine predictions and real labels 
  test_results <- as.data.frame(cbind(test.predictions, real = test_y, test_clin))
  
  # controls
  con.predictions <- predict(model, 
                              data.matrix(sub_con),
                              type = 'prob')
  
  
  
  
  # combine predictions and real labels 
  con_results <- as.data.frame(cbind(con.predictions, real = con_y, con_clin))
  
  
  
  return(list(test_results, importance, con_results))
  
  
}



# training_dat = all_450
# test_dat = all_850
# age_cutoff = age_cutoff
# gender = gender
# tech = tech
# control_age = control_age
# alpha_value = alpha_num
# bh_features = bh_features
run_enet_surv <- function(training_dat,
                         test_dat,
                         fit_type,
                         outcome_type,
                         age_cutoff,
                         gender, 
                         tech,
                         alpha_value,
                         control_age,
                         bh_features) {
  
  
  library(peperr)
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    intersected_feats <- c('Female', 'Male', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('batch_1', 'batch_2', intersected_feats)
  }
  
  if(control_age){
    intersected_feats <- c('age_1', 'age_2', intersected_feats)
  }
  
  # get survival time as days to onset and days to sample collection in one column for training dat
  time_to_event <- training_dat$age_diagnosis
  missing_ind <- is.na(time_to_event)
  time_to_collection <- training_dat$age_sample_collection
  
  time_to_event[missing_ind] <- time_to_collection[missing_ind]
  
  # get cancer status 
  if(outcome_type == 'cancer'){
    cancer_status <- ifelse(training_dat$cancer_diagnosis_diagnoses != 'Unaffected', 1, 0)
    
  } else {
    cancer_status <- training_dat$age_var
  }
  
  # for test dat
  time_to_event_test <- test_dat$age_diagnosis
  missing_ind_test <- is.na(time_to_event_test)
  time_to_collection_test <- test_dat$age_sample_collection
  
  time_to_event_test[missing_ind_test] <- time_to_collection_test[missing_ind_test]
  
  if(outcome_type == 'cancer'){
    # get cancer status 
    cancer_status_test <- ifelse(test_dat$cancer_diagnosis_diagnoses != 'Unaffected', 1, 0)
    
  } else {
    cancer_status_test <- test_dat$age_var
  }
  
  
  # # refactor
  
  # cancer_status <- factor(cancer_status, levels = c('positive', 'negative'))
  # cancer_status_test <- factor(cancer_status_test, levels = c('positive', 'negative'))
  # 
  
  # get clinical data
  test_clin <- test_dat[, !grepl('^cg', colnames(test_dat))]
  
  

  # get model data
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]

  # get survival object
  surv_outcome <- Surv(time_to_event,cancer_status, type = 'right')
  # get survival object
  surv_outcome_test <- Surv(time_to_event_test,cancer_status_test, type = 'right')
  
  
  
  if(fit_type == 'coxph'){
    training_dat <- as.data.frame(cbind(surv_outcome = surv_outcome, training_dat))
    # fit coxph
    library(survival)
    model <- coxph(surv_outcome ~., data = training_dat)
    test.predictions <- predict(model, newdata = test_dat, type = 'lp')
    
  }
  
  if(fit_type == 'rf'){
    library(rpart)
    
    training_dat <- as.data.frame(cbind(surv_outcome = surv_outcome, training_dat))
    # fit coxph
    model <- rpart(surv_outcome ~., data = training_dat)
    test.predictions <- predict(model, newdata = test_dat)
    
  }
  
  if(fit_type == 'enet'){
    # fit coxph
    # start elastic net tuning
    N_CV_REPEATS = 2
    nfolds = 5
    
    
    # set parameters for training model
    type_family <- 'cox'
    type_measure <- 'deviance'
    best_alpha <- alpha_value
    temp.non_zero_coeff = 0
    temp.loop_count = 0
    # loop runs initially because temp.non_zero coefficient <3 and then stops 
    # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
    # it they are never greater than 1, then the model does not converge. 
    while (temp.non_zero_coeff < 1) { 
      elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
                                       , surv_outcome
                                       , alpha = best_alpha
                                       , type.measure = type_measure
                                       , family = type_family
                                       , standardize=FALSE
                                       , nlambda = 100
                                       , nfolds = nfolds
                                       , parallel = TRUE
      )
      
      # lambda 
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
                    ,alpha = best_alpha
                    ,standardize=FALSE
                    ,nlambda = 100
                    ,family = type_family)
    
    
    # This returns 100 prediction with 1-100 lambdas
    temp_test.predictions <- predict.coxnet(model, 
                                            data.matrix(test_dat))
    
    
    # get predictions with corresponding lambda.
    test.predictions <- temp_test.predictions[, temp.min_lambda_index]
    
    
  }
  
 
    
    
    # combine predictions and real labels 
    temp_dat <- as.data.frame(cbind(test_pred = test.predictions,  test_clin))
    temp_dat$alpha = alpha_value
  
  
  
  ###########################################################################################
  return(list(temp_dat, model))

}


# cases = cases_450
# controls = con_all
# valid = cases_850
# null_450 = use_null_450
# age_cutoff = age_cutoff
# gender = gender
# tech = tech
# test_lambda = train_lambda
# alpha_value = alpha_num
# lambda_value = s_num
# control_age = FALSE
# bh_features = bh_features
test_model_enet <- function(cases, 
                            controls, 
                            valid, 
                            null_450,
                            use_p53,
                            use_6,
                            trained_lambda,
                            gender,
                            tech,
                            control_age,
                            age_cutoff,
                            test_lambda,
                            alpha_value, 
                            lambda_value,
                            bh_features) {
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(cases))
  
  if(gender) {
    intersected_feats <- c('Female', 'Male', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('batch_1', 'batch_2', intersected_feats)
  }
  if(control_age){
    intersected_feats <- c('age_1', 'age_2', intersected_feats)
  }
  
  if(use_p53){
    intersected_feats <- c('MUT', 'WT', intersected_feats)
    
    
  }
  
  if(null_450 == 'used_null_450_all'){
    null_450_all <- controls[controls$tech == '450k',]
    null_450_all$age <- null_450_all$age_sample_collection
    if(!use_6){
      null_450_all <- null_450_all[null_450_all$age > 72, ]
      
    }
    cases$age <- cases$age_diagnosis
    controls <- controls[controls$tech == '850k',] 
    cases <- rbind(cases,
                   null_450_all)
    cases_y <- as.factor(ifelse(cases$age < age_cutoff, 'positive', 'negative'))
    
    cases <- cases[!duplicated(cases$tm_donor),]
    con_y <- as.factor(ifelse(controls$age_sample_collection< age_cutoff, 'positive', 'negative'))
    valid_y <-as.factor(ifelse(valid$age_diagnosis < age_cutoff, 'positive', 'negative'))
    cases_y <- factor(cases_y, levels = c('positive', 'negative'))
    con_y <- factor(con_y, levels = c('positive', 'negative'))
    valid_y <- factor(valid_y, levels = c('positive', 'negative'))
    
    
    
  } else  if(null_450 == 'used_null_450_mut'){
    controls_wt <- controls[controls$p53_germline == 'WT',]
    null_450_all <- controls[controls$tech == '450k',]
    null_450_all <- null_450_all[null_450_all$p53_germline == 'MUT',]
    null_450_all$age <- null_450_all$age_sample_collection
    if(!use_6){
      null_450_all <- null_450_all[null_450_all$age > 72, ]
      
    }    
    cases$age <- cases$age_diagnosis
    controls <- controls[controls$tech == '850k',]
    controls <- rbind(controls, controls_wt)
    cases <- rbind(cases,
                   null_450_all)
    cases_y <- as.factor(ifelse(cases$age < age_cutoff, 'positive', 'negative'))
    cases <- cases[!duplicated(cases$tm_donor),]
    con_y <- as.factor(ifelse(controls$age_sample_collection< age_cutoff, 'positive', 'negative'))
    valid_y <-as.factor(ifelse(valid$age_diagnosis < age_cutoff, 'positive', 'negative'))
    cases_y <- factor(cases_y, levels = c('positive', 'negative'))
    con_y <- factor(con_y, levels = c('positive', 'negative'))
    valid_y <- factor(valid_y, levels = c('positive', 'negative'))
    
    
    
  } else {
    # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
    # # get y
    cases_y <- as.factor(ifelse(cases$age_diagnosis < age_cutoff, 'positive', 'negative'))
    con_y <- as.factor(ifelse(controls$age_sample_collection< age_cutoff, 'positive', 'negative'))
    valid_y <-as.factor(ifelse(valid$age_diagnosis < age_cutoff, 'positive', 'negative'))
    
    cases_y <- factor(cases_y, levels = c('positive', 'negative'))
    con_y <- factor(con_y, levels = c('positive', 'negative'))
    valid_y <- factor(valid_y, levels = c('positive', 'negative'))
    
    
  }
  
  
  # get clinical data
  con_clin <- controls[, !grepl('^cg', colnames(controls))]
  valid_clin <- valid[, !grepl('^cg', colnames(valid))]
  
  
  # if(use_offset){
  #   offsetted_train <- as.numeric(training_dat$age_sample_collection)
  #   offsetted_test <- as.numeric(test_dat$age_sample_collection)
  #   
  # } else {
  #   offsetted <- NULL
  # }
  # get model data
  cases <- cases[, intersected_feats]
  controls <- controls[, intersected_feats]
  valid <- valid[, intersected_feats]
  
  
  
  # store fixed values
  best_alpha <- alpha_value
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  nfolds = 5
  
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  elastic_net.cv_model = cv.glmnet(x = as.matrix(cases)
                                   , y =  cases_y
                                   , alpha = alpha_value
                                   , type.measure = type_measure
                                   , family = type_family
                                   , standardize=FALSE
                                   , nlambda = 100
                                   , nfolds = nfolds
                                   , parallel = TRUE)
  
  
  # get outcome variables and clin variables
  lambda_s <- elastic_net.cv_model$lambda.min
  lambda_s_train <- lambda_value
  temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
  trained_index <- which(abs(elastic_net.cv_model$lambda-lambda_s_train)==min(abs(elastic_net.cv_model$lambda-lambda_s_train)))
  # # number of non zero coefficients at that lambda    
  temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
  
  # print(temp.non_zero_coeff)  
  
  model  = glmnet(x = as.matrix(cases)
                  , y =  cases_y
                  ,alpha = best_alpha
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  
  # Predictions on controls data
  
  if(!test_lambda){
    # This returns 100 prediction with 1-100 lambdas
    # This returns 100 prediction with 1-100 lambdas
    temp_test.predictions_con <- predict(model, 
                                         data.matrix(controls),
                                         type = 'response')
    # get predictions with corresponding lambda.
    test.predictions_con <- temp_test.predictions_con[, trained_index]
    
    message('used cv lambda')
  } else {
    # This returns 100 prediction with 1-100 lambdas
    temp_test.predictions_con <- predict(model, 
                                         data.matrix(controls),
                                         type = 'response')
    # get predictions with corresponding lambda.
    test.predictions_con <- temp_test.predictions_con[, temp.min_lambda_index]
    
  }
  
  # combine predictions and real labels 
  temp_dat_con <- as.data.frame(cbind(controls_age_pred = test.predictions_con, controls_age_label = con_y, con_clin))
  temp_dat_con$alpha <- best_alpha
  temp_dat_con$non_zero <- temp.non_zero_coeff
  
  if(!test_lambda){
    # This returns 100 prediction with 1-100 lambdas
    # This returns 100 prediction with 1-100 lambdas
    temp_test.predictions_valid <- predict(model, 
                                           data.matrix(valid),
                                           type = 'response')
    # get predictions with corresponding lambda.
    test.predictions_valid <- temp_test.predictions_valid[, trained_index]
    
    
  } else {
    # This returns 100 prediction with 1-100 lambdas
    temp_test.predictions_valid <- predict(model, 
                                           data.matrix(valid),
                                           type = 'response')
    # get predictions with corresponding lambda.
    test.predictions_valid <- temp_test.predictions_valid[, temp.min_lambda_index]
    
  }
  
  # combine predictions and real labels 
  temp_dat_valid <- as.data.frame(cbind(valid_age_pred = test.predictions_valid, valid_age_label = valid_y, valid_clin))
  temp_dat_valid$alpha <- best_alpha
  temp_dat_valid$non_zero <- temp.non_zero_coeff
  
  
  
  return(list(temp_dat_con, temp_dat_valid, model))
  
  
}



# testing using random forest
test_model_rf <- function(cases, 
                          controls, 
                          valid, 
                          null_450,
                          use_6,
                          use_p53,
                          gender,
                          tech,
                          control_age,
                          optimal_cutoff,
                          age_cutoff,
                          bh_features) {
  
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(cases))
  
  if(gender) {
    intersected_feats <- c('Female', 'Male', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('batch_1', 'batch_2', intersected_feats)
  }
  if(control_age){
    intersected_feats <- c('age_1', 'age_2', intersected_feats)
  }
  
  if(use_p53){
    intersected_feats <- c('MUT', 'WT', intersected_feats)
   
    
  }
  
  if(null_450 == 'used_null_450_all'){
    null_450_all <- controls[controls$tech == '450k',]
    null_450_all$age <- 'negative'
    if(!use_6){
      null_450_all <- null_450_all[null_450_all$age > 72, ]
      cases$age <- cases$age_diagnosis
      controls <- controls[controls$tech == '850k',] 
      cases$age <- as.factor(ifelse(cases$age < age_cutoff, 'positive', 'negative'))
      cases <- rbind(cases,
                     null_450_all)
      
      cases <- cases[!duplicated(cases$tm_donor),]
      cases_y <- factor(cases$age, levels = c('positive', 'negative'))
      con_y <- as.factor(ifelse(controls$age_sample_collection< age_cutoff, 'positive', 'negative'))
      valid_y <-as.factor(ifelse(valid$age_diagnosis < age_cutoff, 'positive', 'negative'))
      cases_y <- factor(cases_y, levels = c('positive', 'negative'))
      con_y <- factor(con_y, levels = c('positive', 'negative'))
      valid_y <- factor(valid_y, levels = c('positive', 'negative'))
       
    } else {
      cases$age <-  as.factor(ifelse(cases$age < age_cutoff, 'positive', 'negative'))
      controls <- controls[controls$tech == '850k',] 
      cases <- rbind(cases,
                     null_450_all)
      
      cases <- cases[!duplicated(cases$tm_donor),]
      con_y <- as.factor(ifelse(controls$age_sample_collection< age_cutoff, 'positive', 'negative'))
      valid_y <-as.factor(ifelse(valid$age_diagnosis < age_cutoff, 'positive', 'negative'))
      cases_y <- factor(cases_y, levels = c('positive', 'negative'))
      con_y <- factor(con_y, levels = c('positive', 'negative'))
      valid_y <- factor(valid_y, levels = c('positive', 'negative'))
    }
    
    
    
    
  } else  if(null_450 == 'used_null_450_mut'){
    controls_wt <- controls[controls$p53_germline == 'WT',]
    null_450_all <- controls[controls$tech == '450k',]
    null_450_all <- null_450_all[null_450_all$p53_germline == 'MUT',]
    null_450_all$age <- 'negative'
    
    if(!use_6){
      null_450_all <- null_450_all[null_450_all$age > 72, ]
      cases$age <- cases$age_diagnosis
      controls <- controls[controls$tech == '850k',] 
      cases$age <- as.factor(ifelse(cases$age < age_cutoff, 'positive', 'negative'))
      cases <- rbind(cases,
                     null_450_all)
      
      cases <- cases[!duplicated(cases$tm_donor),]
      cases_y <- factor(cases$age, levels = c('positive', 'negative'))
      con_y <- as.factor(ifelse(controls$age_sample_collection< age_cutoff, 'positive', 'negative'))
      valid_y <-as.factor(ifelse(valid$age_diagnosis < age_cutoff, 'positive', 'negative'))
      cases_y <- factor(cases_y, levels = c('positive', 'negative'))
      con_y <- factor(con_y, levels = c('positive', 'negative'))
      valid_y <- factor(valid_y, levels = c('positive', 'negative'))
      
    } else {
      cases$age <- cases$age_diagnosis
      controls <- controls[controls$tech == '850k',] 
      cases <- rbind(cases,
                     null_450_all)
      cases_y <- as.factor(ifelse(cases$age < age_cutoff, 'positive', 'negative'))
      
      cases <- cases[!duplicated(cases$tm_donor),]
      con_y <- as.factor(ifelse(controls$age_sample_collection< age_cutoff, 'positive', 'negative'))
      valid_y <-as.factor(ifelse(valid$age_diagnosis < age_cutoff, 'positive', 'negative'))
      cases_y <- factor(cases_y, levels = c('positive', 'negative'))
      con_y <- factor(con_y, levels = c('positive', 'negative'))
      valid_y <- factor(valid_y, levels = c('positive', 'negative'))
    }
    
    
  } else {
    # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
    # # get y
    cases_y <- as.factor(ifelse(cases$age_diagnosis < age_cutoff, 'positive', 'negative'))
    con_y <- as.factor(ifelse(controls$age_sample_collection< age_cutoff, 'positive', 'negative'))
    valid_y <-as.factor(ifelse(valid$age_diagnosis < age_cutoff, 'positive', 'negative'))
    
    cases_y <- factor(cases_y, levels = c('positive', 'negative'))
    con_y <- factor(con_y, levels = c('positive', 'negative'))
    valid_y <- factor(valid_y, levels = c('positive', 'negative'))
    
    
  }
  
 
  # get clinical data
  con_clin <- controls[, !grepl('^cg', colnames(controls))]
  valid_clin <- valid[, !grepl('^cg', colnames(valid))]
  
  
  # if(use_offset){
  #   offsetted_train <- as.numeric(training_dat$age_sample_collection)
  #   offsetted_test <- as.numeric(test_dat$age_sample_collection)
  #   
  # } else {
  #   offsetted <- NULL
  # }
  # get model data
  cases <- cases[, intersected_feats]
  controls <- controls[, intersected_feats]
  valid <- valid[, intersected_feats]
  
  
  set.seed(1)
  # determines how you train the model.
  NFOLDS <- 5
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = min(10, NFOLDS),
    classProbs = TRUE,
    repeats = 0,
    allowParallel = TRUE,
    summaryFunction = twoClassSummary
    
  )
  
  # mtry: Number of variables randomly sampled as candidates at each split.
  # ntree: Number of trees to grow.
  
  mtry <- sqrt(ncol(cases[,colnames(cases)]))
  tunegrid <- expand.grid(.mtry=mtry)
  
  model <- train(x = cases
                 , y = cases_y
                 , metric = 'ROC'
                 , method = "rf"
                 , trControl = fitControl
                 , tuneGrid = tunegrid
                 , importance = T
                 , verbose = FALSE)
  
  temp <- varImp(model)[[1]]
  importance <- cbind(rownames(temp), temp$positive)
  importance <- as.data.frame(importance)
  importance$V2 <- round(as.numeric(as.character(importance$V2)), 2)
  names(importance) <- c('probe', 'score')
  
  # Predictions on controls
  test.predictions_con <- predict(model, 
                              data.matrix(controls),
                              type = 'prob')
  
  
  # combine predictions and real labels 
  test_results_con <- as.data.frame(cbind(test.predictions_con, real = con_y, con_clin))
  test_results_con$pred_class <- as.factor(ifelse(test_results_con$positive > optimal_cutoff, 'positive', 'negative'))
  test_results_con$optimal_cutoff <- optimal_cutoff
  # relevel both factore
  test_results_con$pred_class <- factor(test_results_con$pred_class, c('positive', 'negative'))
  test_results_con$real <- factor(test_results_con$real, c('positive', 'negative'))
  
  test_results_con$accuracy <- caret::confusionMatrix(table(test_results_con$pred_class, test_results_con$real))$overall[1]
  test_results_con$tot_probes <- ncol(cases)
  
  
  # Predictions on validation
  test.predictions_valid <- predict(model, 
                                  data.matrix(valid),
                                  type = 'prob')
  
  
  # combine predictions and real labels 
  test_results_valid <- as.data.frame(cbind(test.predictions_valid, real = valid_y, valid_clin))
  test_results_valid$pred_class <- as.factor(ifelse(test_results_valid$positive > optimal_cutoff, 'positive', 'negative'))
  test_results_valid$optimal_cutoff <- optimal_cutoff
  # relevel both factore
  test_results_valid$pred_class <- factor(test_results_valid$pred_class, c('positive', 'negative'))
  test_results_valid$real <- factor(test_results_valid$real, c('positive', 'negative'))
  
  test_results_valid$accuracy <- caret::confusionMatrix(table(test_results_valid$pred_class, test_results_valid$real))$overall[1]
  test_results_valid$tot_probes <- ncol(cases)
  
  return(list(test_results_valid, test_results_con,importance, model))
  
}




plot_pred <- function(dat, type, plot_type,strategy, log, other_title){
  if(type == 'val'){
    pred_short <- prediction(dat$valid_age_pred, dat$pred_label)
    # optimal cutoff
    cost.perf = performance(pred_short, "cost", cost.fp = 1, cost.fn = 1)
    optimal_cutoff <- pred_short@cutoffs[[1]][which.min(cost.perf@y.values[[1]])]
    
  }
  if(strategy == 'strategy_1'){
    if(log){
      g_title = 'log_'
      
    } else {
      g_title = 'no log_'
    }
  } else if(strategy == 'strategy_2'){
    g_title = 'transform_'
  } else if(strategy == 'strategy_3'){
    g_title = 'pc_'
  } else {
    g_title = 'transform_pc'
  }
  if(plot_type == 'age_pred'){
    if(type == 'con'){
      g_title = paste0(g_title,'null_set')
      g <- ggplot(dat, aes(age_sample_collection, controls_age_pred)) +
        geom_point(size = 1, col = 'black') + 
        labs(x = 'Age in months',
             y = 'Predictions',
             title = g_title)
      
    } else {
      g_title = paste0(g_title,'validation_set')
      g <- ggplot(dat, aes(age_sample_collection, valid_age_pred)) +
        geom_point(size = 1, col = 'black') + 
        labs(x = 'Age in months',
             y = 'Predictions',
             title = g_title)
    }
  } else if(plot_type == 'ROC'){
    
    if(type == 'con'){
      other_title = paste0(g_title, 'null')
      # get dataset of predictions and labels for both small and large data
      pred_short <- prediction(dat$controls_age_pred, dat$controls_age_label)
      # pred_long <- prediction($preds, final_dat$real)
      
      # get performace objects
      perf_s <- performance(prediction.obj = pred_short, measure = 'tpr', x.measure = 'fpr')
      # perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')
      
      # plot mean preds
      plot(perf_s)
      abline(a = 0, b =1)
      
    } else {
      other_title = paste0(g_title, 'validation_set')
      
      # get dataset of predictions and labels for both small and large data
      pred_short <- prediction(dat$valid_age_pred, dat$valid_age_label)
      # pred_long <- prediction($preds, final_dat$real)
      
      # get performace objects
      perf_s <- performance(prediction.obj = pred_short, measure = 'tpr', x.measure = 'fpr')
      # perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')
      
      # plot mean preds
      plot(perf_s)
      abline(a = 0, b =1)
      
    }
  } else {
    if(type == 'con'){
      other_title = paste0(g_title, 'null')
      
      pred_short <- prediction(dat$controls_age_pred, dat$controls_age_label)
      cost.perf = performance(pred_short, "cost", cost.fp = 1, cost.fn = 1)
      optimal_cutoff <- pred_short@cutoffs[[1]][which.min(cost.perf@y.values[[1]])]
      
      # get confusion matrix function for plotting 
      g <- ConfusionMatrixInfo(data = dat, 
                               predict = 'controls_age_pred', 
                               actual = 'controls_age_label', 
                               cutoff = .57,
                               get_plot = TRUE,
                               other_title = other_title)
    } else {
      other_title = paste0(g_title, 'validation_set')
      
      
      # get confusion matrix function for plotting 
      g <-ConfusionMatrixInfo(data = dat, 
                              predict = 'valid_age_pred', 
                              actual = 'pred_label', 
                              cutoff = .57,
                              get_plot = TRUE,
                              other_title = other_title)
      
    }
  }
  return(g)
}

get_acc_val <- function(temp_dat, thresh){
  all_alphas <- (1:10)/10
  result_list <- list()
  for(i in 1:length(all_alphas)){
    this_alpha <- all_alphas[i]
    sub_dat <- temp_dat[temp_dat$alpha == this_alpha,]
    sub_dat$pred_label <-as.factor(ifelse(sub_dat$valid_age_pred > thresh, 'positive', 'negative'))
    sub_dat$pred_label <- factor(sub_dat$pred_label, levels = c('positive', 'negative'))
    sub_dat$valid_age_label <- factor(sub_dat$valid_age_label, levels = c('positive', 'negative'))
    sub_dat$acc <- caret::confusionMatrix(sub_dat$pred_label, sub_dat$valid_age_label)$overall[1]
    result_list[[i]] <- sub_dat
    print(i)
  }
  temp <- do.call('rbind', result_list)
  return(temp)
}


get_acc_val <- function(temp_dat, thresh){
  all_alphas <- unique(temp_dat$alpha)
  result_list <- list()
  for(i in 1:length(all_alphas)){
    this_alpha <- all_alphas[i]
    sub_dat <- temp_dat[temp_dat$alpha == this_alpha,]
    sub_dat$pred_label <-as.factor(ifelse(sub_dat$preds > thresh, 'positive', 'negative'))
    sub_dat$pred_label <- factor(sub_dat$pred_label, levels = c('positive', 'negative'))
    sub_dat$real <- factor(sub_dat$real, levels = c('positive', 'negative'))
    sub_dat$acc <- caret::confusionMatrix(sub_dat$pred_label, sub_dat$real)$overall[1]
    result_list[[i]] <- sub_dat
    print(i)
  }
  temp <- do.call('rbind', result_list)
  return(temp)
}

get_young_labels <- function(temp_dat, thresh, age){
  all_alphas <- (1:10)/10
  result_list <- list()
  for(i in 1:length(all_alphas)){
    this_alpha <- all_alphas[i]
    sub_dat <- temp_dat[temp_dat$alpha == this_alpha,]
    sub_dat <- sub_dat[sub_dat$age_sample_collection < age,]
    sub_dat$pred_label <-as.factor(ifelse(sub_dat$controls_age_pred > thresh, 'positive', 'negative'))
    sub_dat$pred_label <- factor(sub_dat$pred_label, levels = c('positive', 'negative'))
    sub_dat$controls_age_label <- factor(sub_dat$controls_age_label, levels = c('positive', 'negative'))
    sub_dat <- sub_dat[, c('tm_donor','alpha','age_sample_collection','controls_age_pred', 'controls_age_label', 'pred_label')]
    result_list[[i]] <- sub_dat
    print(i)
  }
  temp <- do.call('rbind', result_list)
  return(temp)
}

plot_acc <- function(temp_dat, acc_column, column, bar) {
  
  column_name <- column
  
  
  if(bar){
    temp_dat <- temp_dat[, c(column, acc_column)]
    names(temp_dat) <- c('V1', 'Avg Accuracy')
    
    g1 <- ggplot(temp_dat, 
                 aes(reorder(V1, -`Avg Accuracy`),
                     `Avg Accuracy`)) +
      geom_bar(alpha = 0.6,
               color = 'black',
               fill = 'grey',
               stat = 'identity') +
      geom_smooth(method = 'lm',
                  color = 'red',
                  linetype = 1) +
      labs(title = paste0(column_name, ' and Accuracy'),
           x = column_name,
           y = 'Accuracy') +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.5, size = 5))
  } else {
    temp_dat <- temp_dat[, c(column, acc_column)]
    names(temp_dat) <- c('V1', 'Accuracy')
    
    g1 <- ggplot(temp_dat, 
                 aes(V1, Accuracy)) +
      geom_point(pch = 21,
                 size = 2,
                 alpha = 0.6,
                 color = 'black',
                 fill = 'grey') +
      geom_smooth(method = 'lm',
                  color = 'red',
                  linetype = 1) +
      labs(title = paste0(column_name, ' and Accuracy'),
           x = column_name,
           y = 'Accuracy')
    
  }
  
  return(g1)
  
}


plot_3d_model_means <- function(temp_dat){
  with(temp_dat, {
    s3d <- scatterplot3d(mean_lambda, mean_alpha, mean_acc,        # x y and z axis
                         color="darkgrey", 
                         pch=1,        # filled blue circles
                         type="h",
                         main="Alpha and Lambda choice",
                         xlab="Enet lambda",
                         ylab="Enet alpha",
                         zlab="Model accuracy")
    s3d.coords <- s3d$xyz.convert(mean_lambda, mean_alpha, mean_acc) # convert 3D coords to 2D projection
    my.lm <- lm(temp_dat$mean_acc ~ temp_dat$mean_lambda + temp_dat$mean_alpha )
    s3d$plane3d(my.lm)
    s3d$points3d(mean_lambda, mean_alpha, mean_acc,
                 col = adjustcolor("black", alpha.f = 0.8), type = 'p', pch = 16)
  })
}


plot_3d_model <- function(temp_dat, type){
  with(temp_dat, {
    s3d <- scatterplot3d(lambda, alpha, accuracy,        # x y and z axis
                         color=adjustcolor('grey', alpha.f = 0.2), 
                         pch=1,        # filled blue circles
                         type=type,
                         main="Alpha and Lambda choice",
                         xlab="Enet lambda",
                         ylab="Enet alpha",
                         zlab="Model accuracy")
    s3d.coords <- s3d$xyz.convert(lambda, alpha, accuracy) # convert 3D coords to 2D projection
    my.lm <- lm(temp_dat$accuracy ~ temp_dat$lambda + temp_dat$alpha )
    s3d$plane3d(my.lm)
    s3d$points3d(lambda, alpha, accuracy,
                 col = adjustcolor("black", alpha.f = 0.1), type = 'p', pch = 16)
  })
}


AccuracyCutoffInfo <- function( test, predict, actual )
{
  # change the cutoff value's range as you please 
  cutoff <- seq( .4, .8, by = .05 )
  
  accuracy <- lapply( cutoff, function(c)
  {
    # use the confusionMatrix from the caret package
    cm_test  <- confusionMatrix( as.numeric( test[[predict]]  > c ), test[[actual]]  )
    
    dt <- data.table( cutoff = c,
                      test   = cm_test$overall[["Accuracy"]] )
    return(dt)
  }) %>% rbindlist()
  
  # visualize the accuracy of the train and test set for different cutoff value 
  # accuracy in percentage.
  accuracy_long <- gather( accuracy, "data", "accuracy", -1 )
  
  plot <- ggplot( accuracy_long, aes( cutoff, accuracy) ) + 
    geom_line( size = 1 ) + geom_point( size = 3 ) +
    scale_y_continuous( label = percent ) +
    ggtitle( "Cutoff" )
  
  return( list( data = accuracy, plot = plot ) )
}


# ------------------------------------------------------------------------------------------
# [ConfusionMatrixInfo] : 
# Obtain the confusion matrix plot and data.table for a given
# dataset that already consists the predicted score and actual outcome.
# @data    : your data.table or data.frame type data that consists the column
#            of the predicted score and actual outcome 
# @predict : predicted score's column name
# @actual  : actual results' column name
# @cutoff  : cutoff value for the prediction score 
# return   : 1. data : a data.table consisting of three column
#            		   the first two stores the original value of the prediction and actual outcome from
#			 		   the passed in data frame, the third indicates the type, which is after choosing the 
#			 		   cutoff value, will this row be a true/false positive/ negative 
#            2. plot : plot that visualizes the data.table 

# 
# data <- temp_valid
# predict <- 'positive'
# actual = 'real'
# cutoff = 0.5
ConfusionMatrixInfo <- function( data, predict, actual, cutoff, get_plot, other_title, data_type)
{	
  # extract the column ;
  # relevel making 1 appears on the more commonly seen position in 
  # a two by two confusion matrix	
  predict <- data[[predict]]
  actual  <- relevel( as.factor( data[[actual]] ), "positive" )
  if(data_type == 'null') {
    age <- data$age_sample_collection
    
  } else {
    age <- data$age_diagnosis
    
  }
  result <- data.table( actual = actual, predict = predict, age = age)
  
  # caculating each pred falls into which category for the confusion matrix
  result[ , type := ifelse( predict >= cutoff & actual == 'positive', "TP",
                            ifelse( predict >= cutoff & actual == 'negative', "FP", 
                                    ifelse( predict <  cutoff & actual == 'positive', "FN", "TN" ) ) ) %>% as.factor() ]
  
  result <- as.data.frame(result)
  if(data_type == 'null'){
    result$actual <- as.factor(ifelse(result$actual == 'positive', 'Under age of 6', 'Over age of 6'))
    result$acual <- factor(result$actual, levels = c('Under age of 6', 'Over age of 6'))
    x_lab = 'Real age of sample collection'
  } else {
    result$actual <- as.factor(ifelse(result$actual == 'positive', 'Onset before 6', 'Onset after 6'))
    result$acual <- factor(result$actual, levels = c('Onset before 6', 'Onset after 6'))
    x_lab = 'Real age of onset'
    
  }
  result$actual<- with(result, reorder(actual, acual, function(x) length(x)))  
  library(ggthemes)
  library(ggrepel)
  # jittering : can spread the points along the x axis 
  result$age <- round(result$age/12)
  plot <- ggplot( result, aes( actual, predict, color = type ) ) + 
    geom_point(size = 0, show.legend = FALSE) +
    scale_color_manual(name = 'Result',
                       values = c('red', 'orange','green', 'blue'),
                       breaks = c( "TP", "FN", "FP", "TN" ))+
    geom_violin( fill = "white", color = NA ) +
    geom_hline( yintercept = cutoff, color = 'black', alpha = 0.6, linetype = 2 ) + 
    geom_text(aes(label = age),alpha = 0.7, cex = 5,position=position_jitter(width = 0.4, height = 0), show.legend = FALSE)+
    geom_vline(xintercept = 1.5, linetype = 2) +
    scale_y_continuous( limits = c( 0, 1 ) ) + 
    guides( col = guide_legend( nrow = 2 ) ) + # adjust the legend to have two rows  
    ggtitle( sprintf( other_title,"_Cutoff at %.2f", cutoff ) ) +
    theme(text = element_text(size=14)) + 
    xlab(x_lab) + ylab('Predicted probability of onset') + theme_base(base_size = 18)
  
  
  
  if(get_plot) {
    return(plot)
  } else {
    return(as.data.frame(result))
  }
}


# [ROCInfo] : 
# Pass in the data that already consists the predicted score and actual outcome.
# to obtain the ROC curve 
# @data    : your data.table or data.frame type data that consists the column
#            of the predicted score and actual outcome
# @predict : predicted score's column name
# @actual  : actual results' column name
# @cost.fp : associated cost for a false positive 
# @cost.fn : associated cost for a false negative 
# return   : a list containing  
#			 1. plot        : a side by side roc and cost plot, title showing optimal cutoff value
# 				 	   		  title showing optimal cutoff, total cost, and area under the curve (auc)
# 		     2. cutoff      : optimal cutoff value according to the specified fp/fn cost 
#		     3. totalcost   : total cost according to the specified fp/fn cost
#			 4. auc 		: area under the curve
#		     5. sensitivity : TP / (TP + FN)
#		     6. specificity : TN / (FP + TN)


# data = dat_sample 
# predict = 'mean_preds' 
# actual = 'real_label' 
# cost.fp = cost_fp
# cost.fn = cost_fn
ROCInfo <- function( data, predict, actual, cost.fp, cost.fn, other_title )
{
  # calculate the values using the ROCR library
  # true positive, false postive 
  pred <- prediction( data[[predict]], data[[actual]] )
  perf <- performance( pred, "tpr", "fpr" )
  roc_dt <- data.frame( fpr = perf@x.values[[1]], tpr = perf@y.values[[1]] )
  
  # cost with the specified false positive and false negative cost 
  # false postive rate * number of negative instances * false positive cost + 
  # false negative rate * number of positive instances * false negative cost
  cost <- perf@x.values[[1]] * cost.fp * sum( data[[actual]] == 'negative' ) + 
    ( 1 - perf@y.values[[1]] ) * cost.fn * sum( data[[actual]] == 'positive' )
  
  cost_dt <- data.frame( cutoff = pred@cutoffs[[1]], cost = cost )
  
  # optimal cutoff value, and the corresponding true positive and false positive rate
  best_index  <- which.min(cost)
  best_cost   <- cost_dt[ best_index, "cost" ]
  best_tpr    <- roc_dt[ best_index, "tpr" ]
  best_fpr    <- roc_dt[ best_index, "fpr" ]
  best_cutoff <- pred@cutoffs[[1]][ best_index ]
  
  # area under the curve
  auc <- performance( pred, "auc" )@y.values[[1]]
  
  # normalize the cost to assign colors to 1
  normalize <- function(v) ( v - min(v) ) / diff( range(v) )
  
  # create color from a palette to assign to the 100 generated threshold between 0 ~ 1
  # then normalize each cost and assign colors to it, the higher the blacker
  # don't times it by 100, there will be 0 in the vector
  col_ramp <- colorRampPalette( c( "green", "orange", "red", "black" ) )(100)   
  col_by_cost <- col_ramp[ ceiling( normalize(cost) * 99 ) + 1 ]
  
  roc_plot <- ggplot( roc_dt, aes( fpr, tpr ) ) + 
    geom_line( color = rgb( 0, 0, 1, alpha = 0.3 ) ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.2 ) + 
    geom_segment( aes( x = 0, y = 0, xend = 1, yend = 1 ), alpha = 0.8, color = "royalblue" ) + 
    labs( title = "ROC", x = "False Postive Rate", y = "True Positive Rate" ) +
    geom_hline( yintercept = best_tpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    geom_vline( xintercept = best_fpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" )				
  
  cost_plot <- ggplot( cost_dt, aes( cutoff, cost ) ) +
    geom_line( color = "blue", alpha = 0.5 ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.5 ) +
    ggtitle( "Cost" ) +
    scale_y_continuous( labels = comma ) +
    geom_vline( xintercept = best_cutoff, alpha = 0.8, linetype = "dashed", color = "steelblue4" )	
  
  options(scipen = '999')
  # the main title for the two arranged plot
  sub_title <- sprintf(other_title,  "Cutoff at %.2f - Total Cost = %a, AUC = %.3f", 
                       best_cutoff, best_cost, auc )
  
  # arranged into a side by side plot
  plot <- arrangeGrob( roc_plot, cost_plot, ncol = 2, 
                       top = textGrob( sub_title, gp = gpar( fontsize = 16, fontface = "bold" ) ) )
  
  return( list( plot 		  = plot, 
                cutoff 	  = best_cutoff, 
                totalcost   = best_cost, 
                auc         = auc,
                sensitivity = best_tpr, 
                specificity = 1 - best_fpr ) )
}

