
##########
# load libraries
##########
library(tidyverse)
library(ROCR)
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
library(GWAStools)



##########
# function that Loops through list, preprocesses, and convert to beta, m, and cn values 
##########

preprocessMethod <- function(data, preprocess) {
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
  m <- getM(Mset)
  # get beta values
  beta <- getBeta(Gset)
  return(beta)
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

# # id function
# beta_data <- beta_valid[1:10000,]
# id_map <- id_map_val
process_rg_set_single <- function(beta_data, id_map, clin) {
  beta_data <- findIds(beta_data, id_map)
  beta_data <- getIdName(beta_data)
  beta_data <- cleanIds(beta_data)
  beta_data <- beta_data[, !grepl('ch', colnames(beta_data))]
  beta_data <- inner_join(clin, beta_data, by = 'ids')
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

# remove ages 
remove_ages <- function(mod_dat, max_age){
    mod_dat <- subset(mod_dat, age_sample_collection <= max_age)
  return(mod_dat)
}

# id function
process_rg_set_single <- function(beta_data, id_map, clin) {
  # cases
  beta_data <- findIds(beta_data, id_map)
  # get id name (only cases)
  beta_data <- getIdName(beta_data)
  # clean ids
  beta_data <- cleanIds(beta_data)
  # remove 'ch' from column names
  beta_data <- beta_data[, !grepl('ch', colnames(beta_data))]
  # inner join
  beta_data <- inner_join(clin, beta_data, by = 'ids')
  # remove NAs from tm_donor 
  beta_data <- beta_data[!is.na(beta_data$tm_donor_),]
  # remove duplicates
  beta_data <- beta_data[!duplicated(beta_data$tm_donor_),]
  # get cg_sites
  cg_sites <- colnames(beta_data)[grepl('cg', colnames(beta_data))]
  # saveRDS(cg_sites, paste0(model_data, '/four_fifty_feats.rda'))
  # subset data by colmns of interest and cg_sites
  beta_data <- beta_data[, c('ids', 
                               'p53_germline', 
                               'cancer_diagnosis_diagnoses', 
                               'age_diagnosis',
                               'age_sample_collection',
                               'gender',
                               'sentrix_id',
                               'family_name',
                               'tm_donor_',
                               cg_sites)]
    return(beta_data)
}


########
# binarize age variables find
#########
age_binary <- function(dat, type, cutoff) {
    if (type == 'cases') {
      dat$age_diagnosis <- ifelse(dat$age_diagnosis > cutoff, 'yes', 'no')
      dat$age_sample_collection <- ifelse(dat$age_sample_collection > cutoff, 'yes', 'no')
    } else {
      dat$age_sample_collection <- ifelse(dat$age_sample_collection > cutoff, 'yes', 'no')
    }
  return(dat)
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
    temp <- inner_join(id_map_dat, outliers, by = 'ids')
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

subset_rg_set <- function(rg_set, keep_gender, keep_controls, keep_snps, get_island, get_type, get_chr) {
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
  remove_probes <- rg_dat$Name[!as.character(rg_dat$Name) %in% as.character(rg_dat_sub$Name)]
  # use subset function from minfi
  rg_set_new <- subsetByLoci(rg_set, 
                             includeLoci = keep_probes, 
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

# ##########
# # scale data
# ##########
# dat <- betaCases 
# probe_start <- 8
# scaleData <- function(dat, probe_start)
# {
#   
#   # get row statistics
#   colMean <- apply(dat[, probe_start:ncol(dat)], 2, mean)
#   colSd <- apply(dat[, probe_start:ncol(dat)], 2, sd)
#   # constantInd <- rowSd==0
#   # rowSd[constantInd] <- 1
#   colStats <- list(mean=colMean, sd=colSd)
#   
#   # apply normilization
#   dat[, probe_start:ncol(dat)]  <- (dat[, probe_start:ncol(dat)] - colStats$mean) / colStats$sd
#   
#   return(dat)
# }

##########
# impute and scale for raw data
##########
scaleImputeDat <- function(dat, scale) {
  if (scale) {
    # get row statistics
    rowMean <- apply(dat, 1, mean, na.rm=TRUE)
    rowSd <- apply(dat, 1, sd, na.rm=TRUE)
    # constantInd <- rowSd==0
    # rowSd[constantInd] <- 1
    rowStats <- list(mean=rowMean, sd=rowSd)
    # apply normilization
    dat  <- (dat - rowStats$mean) / rowStats$sd
    # make matrix
    dat <- as.matrix(dat)
    # impute with knn
    dat_knn <-  impute.knn(dat, k = 10)$data
    # get clin data back and return final data
    final_dat <- dat_knn
  } else {
    # make matrix
    dat <- as.matrix(dat)
    # impute with knn
    dat_knn <-  impute.knn(dat, k = 10)$data
    # get clin data back and return final data
    final_dat <- dat_knn
  }
  return(final_dat)
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

##########
# function for get m values 
##########

get_m_values <- 
  function(data_set, probe_start) {
    
    data_set[,probe_start:ncol(data_set)] <- apply(data_set[, probe_start:ncol(data_set)], 2, 
                                                   function(x) log(x/(1-x)))
    # log(data_set[, probe_start:ncol(data_set)]/(1- data_set[, probe_start:ncol(data_set)]))
    
    return(data_set)
    
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
##########
# Main function that specifies a preprocessing method and get beta
##########
getMethyl <- function(data_list, cases , method) {
  
  processed_list <-preprocessMethod(data_list, preprocess = method)
  
  # save.image('/home/benbrew/Desktop/temp_process.RData')
  # load('/home/benbrew/Desktop/temp_process.RData')
  # 
  # combinelist
  beta_methyl <- combineList(processed_list)
  # m_methyl <- combineList(processed_list[[2]])
  # cn_methyl <- combineList(processed_list[[3]])
  
  if (cases) {
    
    beta_methyl <- findIds(beta_methyl, id_map)
    # clean ids
    beta_methyl <- getIdName(beta_methyl)
    # m_methyl <- getIdName(m_methyl)
    # cn_methyl <- getIdName(cn_methyl)
    
    # make data frame
    beta_methyl <- as.data.frame(beta_methyl, stringsAsFactors = F)
    
  } else {
    
    # find ids
    beta_methyl <- findIds(beta_methyl, id_map_control)
    
    # make data frame
    beta_methyl <- as.data.frame(beta_methyl, stringsAsFactors = F)
    
    # m_methyl <- findIds(m_methyl, id_map_control)
    # cn_methyl <- findIds(cn_methyl, id_map_control)
  }
  
  
  return(beta_methyl)
  
}

# clean ids in each data set 
cleanIds <- function(data){
  
  data$ids <- gsub('A|B|_|-', '', data$ids)
  data$ids <- substr(data$ids, 1,4) 
  return(data)
}
# get probe locations 

getIds <- function(cg_locations) {
  
  #idat files
  idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
  sapply(idatFiles, gunzip, overwrite = TRUE)
  # read into rgSet
  rgSet <- read.450k.exp("GSE68777/idat")
  # preprocess quantil
  rgSet <- preprocessQuantile(rgSet)
  # get rangers 
  rgSet <- granges(rgSet)
  cg_locations <- as.data.frame(rgSet)
  # make rownames probe column
  cg_locations$probe <- rownames(cg_locations)
  rownames(cg_locations) <- NULL
  return(cg_locations)
}


# # get model data 
# combined_data <- m_test_con_full
# train_data <- m_train_mod

combine_clean_split <- function(combined_data, train_data, controls){
  
  # give temp indicator for training and testing 
  combined_data$train_test <- 'test'
  train_data$train_test <- 'train'
  
  # intersected feats 
  intersect_names <- Reduce(intersect, list(colnames(combined_data)[10:ncol(combined_data)],
                                            colnames(train_data)[10:ncol(train_data)]))
  
  # subset by intersected feats
  combined_data <- combined_data[, c('ids',
                                     'p53_germline',
                                     'cancer_diagnosis_diagnoses',
                                     'age_diagnosis',
                                     'age_sample_collection',
                                     'gender',
                                     'sentrix_id',
                                     'family_name',
                                     'tm_donor_',
                                     intersect_names)]
  
  # subset by intersected feats
  train_data <- train_data[, c('ids',
                               'p53_germline',
                               'cancer_diagnosis_diagnoses',
                               'age_diagnosis',
                               'age_sample_collection',
                               'gender',
                               'sentrix_id',
                               'family_name',
                               'tm_donor_',
                               intersect_names)]
  
  # combined data
  all_data <- rbind(combined_data, train_data)
  
  # remove inf 
  all_data <- removeInf(all_data, probe_start = 10)
  
  # remove duplicates 
  all_data <- all_data[!duplicated(all_data$ids),]
  
  # remove WT 
  all_data <- subset(all_data, p53_germline == 'Mut')
  
  # split data
  train_data <- subset(all_data, train_test == 'train')
  test_data <- subset(all_data, train_test == 'test')
  
  # remove column
  train_data$train_test <- test_data$train_test <- NULL
  
  if(controls) {
    
    # split test data 
    test_cases <- subset(test_data, cancer_diagnosis_diagnoses != 'Unaffected')
    test_other <- subset(test_data, cancer_diagnosis_diagnoses == 'Unaffected')
    
    # remove NA from sample controls
    test_other <- test_other[!is.na(test_other$age_sample_collection),]
    
  } else {
    
    # split test data 
    test_cases <- test_data[!grepl('^20', test_data$sentrix_id),]
    test_other <- test_data[grepl('^20', test_data$sentrix_id),]
    
    # remove NAs from validation
    test_other <- test_other[!is.na(test_other$age_diagnosis),]
    test_other <- test_other[!is.na(test_other$age_sample_collection),]
    
    # remove samples from validation that are in cases
    test_other <- test_other[!test_other$ids %in% test_cases$ids,]
    test_other <- test_other[!test_other$ids %in% train_data$ids,]
    
  }
  
  # remove na on training
  train_data <- train_data[!is.na(train_data$age_sample_collection),]
  train_data <- train_data[!is.na(train_data$age_diagnosis),]
  
  # remove na on test
  test_data <- test_data[!is.na(test_data$age_sample_collection),]
  test_data <- test_data[!is.na(test_data$age_diagnosis),]
  
  return(list(train_data, test_cases, test_other))
    
}

# function for processing rg set
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
    beta_cases <- inner_join(clinical_dat, beta_cases, by = 'ids')
    beta_controls <- inner_join(clinical_dat, beta_controls, by = 'ids')
    
    
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

# data <- betaControls
# function that takes each methylation and merges with clinical - keep ids, family, p53 status, age data
joinData <- function(data, control) {
  
  # get intersection of clin idss and data idss
  intersected_ids <- intersect(data$ids, clin$ids)
  features <- colnames(data)[1:(length(colnames(data)) - 3)]
  
  # loop to combine idsentifiers, without merging large table
  data$p53_germline <- NA
  data$age_diagnosis <- NA
  data$cancer_diagnosis_diagnoses <- NA
  data$age_sample_collection <- NA
  data$tm_donor_ <- NA
  data$gender <- NA
  data$family_name <- NA
  
  
  if (!control) {
    
    for (i in intersected_ids) {
      
      data$p53_germline[data$ids == i] <- clin$p53_germline[which(clin$ids == i)]
      data$age_diagnosis[data$ids == i] <- clin$age_diagnosis[which(clin$ids == i)]
      data$cancer_diagnosis_diagnoses[data$ids == i] <- clin$cancer_diagnosis_diagnoses[which(clin$ids == i)]
      data$age_sample_collection[data$ids == i] <- clin$age_sample_collection[which(clin$ids == i)]
      data$tm_donor_[data$ids == i] <- clin$tm_donor_[which(clin$ids == i)]
      data$gender[data$ids == i] <- clin$gender[which(clin$ids == i)]
      data$family_name[data$ids == i] <- clin$family_name[which(clin$ids == i)]
      
      
      
      print(i)
    } 
    data <- data[!is.na(data$p53_germline),]
    data <- data[!duplicated(data$ids),]
    data <- data[!duplicated(data$tm_donor_),]
    # data <- data[!is.na(data$age_diagnosis),]
    data <- data[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses',
                     'age_sample_collection', 'gender','sentrix_id', 'family_name', features)]
    
  } else {
    
    for (i in intersected_ids) {
      
      data$p53_germline[data$ids == i] <- clin$p53_germline[which(clin$ids == i)]
      data$cancer_diagnosis_diagnoses[data$ids == i] <- clin$cancer_diagnosis_diagnoses[which(clin$ids == i)]
      data$age_sample_collection[data$ids == i] <- clin$age_sample_collection[which(clin$ids == i)]
      data$tm_donor_[data$ids == i] <- clin$tm_donor_[which(clin$ids == i)]
      data$gender[data$ids == i] <- clin$gender[which(clin$ids == i)]
      data$family_name[data$ids == i] <- clin$family_name[which(clin$ids == i)]
      
      
      print(i)
    } 
    data <- data[!is.na(data$p53_germline),]
    data <- data[!duplicated(data$ids),]
    data <- data[!duplicated(data$tm_donor_),]
    data <- data[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses',
                     'age_sample_collection', 'gender', 'sentrix_id', 'family_name',features)]
  }
  
  return(data)
}

##########
# remove cancer from controls 
##########
removeCancer <- function(data_controls) 
{
  data_controls <- data_controls[grepl('Unaffected', data_controls$cancer_diagnosis_diagnoses),]
  data_controls <- data_controls[!duplicated(data_controls$ids)]
  return(data_controls)
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

##########
# get old controls 
##########
getControls <- function(data, mut) 
{
  data <- data[grepl('Unaffected', data$cancer_diagnosis_diagnoses),]
  
  if (mut) {
    # subset data by not na in age of diagnosis and mut
    data <- data[grepl('Mut', data$p53_germline),]
  } else {
    data <- data[grepl('WT', data$p53_germline),]
    
  }
  
  return(data)
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
# batch correction
##########




# data_controls <- controls_wt
getBalAge <- function(data_controls, full)
{
  
  # # # # balance age
  # hist(cases$age_sample_collection)
  # hist(data_controls$age_sample_collection)
  
  # remove a few from ranges 100-200, 300-400
  # randomly remove controls that have a less than 50 month age of diganosis to have balanced classes
  remove_index <- which((data_controls$age_sample_collection >= 100 & data_controls$age_sample_collection <= 250) |
                          (data_controls$age_sample_collection >= 300 & data_controls$age_sample_collection <= 400))
  
  if(full) {
    remove_index <- sample(remove_index, 8, replace = F)
    
  } else {
    remove_index <- sample(remove_index, 2, replace = F)
    
  }
  
  data_controls <- data_controls[-remove_index,]
  
  return(data_controls)
  
}


bumpHunterSurv <- function(dat_cases,
                           dat_controls) 
{
  
  
  # combine data
  dat <- rbind(dat_cases, dat_controls)
  
  
  ##########
  # get clinical dat 
  ##########
  bump_clin <- dat[,1:4]
  
  # recode type
  dat$type <- ifelse(grepl('Unaffected', dat$cancer_diagnosis_diagnoses), 'controls', 'cases')
  
  ##########
  # get indicator and put into design matrix with intercept 1
  #########
  indicator_vector <- as.factor(dat$type)
  designMatrix <- cbind(rep(1, nrow(dat)), indicator_vector)
  designMatrix <- as.matrix(designMatrix)
  
  ##########
  # Get genetic locations
  ##########
  dat$p53_germline <- dat$age_diagnosis <- dat$cancer_diagnosis_diagnoses <- dat$ids <- dat$batch <- 
    dat$age_sample_collection <- dat$id <- dat$type <- dat$gender <-  dat$sentrix_id <-  NULL
  # transpose methylation to join with cg_locations to get genetic location vector.
  dat <- as.data.frame(t(dat), stringsAsFactors = F)
  
  # make probe a column in methyl
  dat$probe <- rownames(dat)
  rownames(dat) <- NULL
  
  # inner join methyl and cg_locations by probe
  methyl_cg <- inner_join(dat, cg_locations, by = 'probe')
  
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
  DELTA_BETA_THRESH = 0.5 # DNAm difference threshold
  NUM_BOOTSTRAPS = 3  # number of randomizations
  
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
                           type = "Beta")
    
    bump_hunter_results[[i]] <- tab[[i]]$table
    bump_hunter_results[[i]]$run <- DELTA_BETA_THRESH[i]
  }
  
  bh_results <- do.call(rbind, bump_hunter_results)
  
  return(bh_results)
  
}

###########
# bumphunter for predictions - WT controls, p53 contorls
###########
bumpHunterPred <- function(dat_controls_wt,
                           dat_controls_mut) 
{
  
  # add columns indicating p53 status (in place of age of diagnosis
  dat_controls_wt$age_diagnosis <- 'WT'
  dat_controls_mut$age_diagnosis <- 'MUT'
  
  
  # combine data
  dat <- rbind(dat_controls_wt, dat_controls_mut)
  
  dat$ids <- NULL
  
  # change variable name for age of diagnosis
  colnames(dat)[1] <- 'status'
  
  ##########
  # get clinical dat 
  ##########
  bump_clin <- dat[,1:4]
  
  
  ##########
  # get indicator and put into design matrix with intercept 1
  #########
  indicator_vector <- as.factor(dat$status)
  designMatrix <- cbind(rep(1, nrow(dat)), indicator_vector)
  designMatrix <- as.matrix(designMatrix)
  
  ##########
  # Get genetic locations
  ##########
  dat$p53_germline <- dat$age_diagnosis <- dat$cancer_diagnosis_diagnoses <- dat$ids <- dat$batch <- 
    dat$age_sample_collection <- dat$id <- dat$status <- dat$gender <-  dat$sentrix_id <-  NULL
  # transpose methylation to join with cg_locations to get genetic location vector.
  dat <- as.data.frame(t(dat), stringsAsFactors = F)
  
  # make probe a column in methyl
  dat$probe <- rownames(dat)
  rownames(dat) <- NULL
  
  # inner join methyl and cg_locations by probe
  methyl_cg <- inner_join(dat, cg_locations, by = 'probe')
  
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
  DELTA_BETA_THRESH = 0.5 # DNAm difference threshold
  NUM_BOOTSTRAPS = 3   # number of randomizations
  
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
                           type = "Beta")
    
    bump_hunter_results[[i]] <- tab[[i]]$table
    bump_hunter_results[[i]]$run <- DELTA_BETA_THRESH[i]
  }
  
  bh_results <- do.call(rbind, bump_hunter_results)
  
  return(bh_results)
  
}


##########
# generate folds and set the seed for random samplings - referenced in arg[]
##########
getFolds <- function(model_dat, seed_number, k){
  # assign folds
  set.seed(seed_number)
  model_dat$folds <- sample(1:k, nrow(model_dat), replace = T)
  
  return(model_dat)
}



getProbe <- function(data) {
  
  results <- list()
  results_data <- list()
  colnames(cg_locations) <- paste0(colnames(cg_locations), '_', 'rgSet')
  
  
  # first make seqnames in cg_locations and chr in tab character vectors for identification
  cg_locations$seqnames_rgSet <- as.character(cg_locations$seqnames_rgSet)
  data$chr <- as.character(data$chr)
  
  # use sql data table to merger validators with model_data based on age range
  result = sqldf("select * from cg_locations
                 inner join data
                 on cg_locations.start_rgSet between data.start and data.end")
  
  # keep only necessary columns
  result <- result[, c('chr' , 'start_rgSet','end_rgSet', 'probe_rgSet', 'p.value', 'fwer', 'run')]
  
  # rename cols
  colnames(result) <- c('chr', 'start', 'end', 'probe', 'p.value', 'fwer', 'run')
  
  # get sig results
  result_sig <- result[result$p.value < 0.05,]
  
  # get fwer resultswf
  result_fwer <- result[result$fwer == 0,]
  
  return(list(result, result_sig, result_fwer))
  
}

getRun <- function(data, run_num)
{
  data <- data[data$run == run_num,]
  data <- data[!duplicated(data$probe),]
  data_feat <- as.character(data$probe)
  return(data_feat)
}


testKS <- function(x, y)
{
  y <- y[!is.na(y)]
  x <- x[!is.na(x)]
  
  # Do x and y come from the same distribution?
  ks.test(jitter(x), jitter(y), alternative = 'two.sided')
  
}

# # # #
# training_dat = cases[train_index,]
# test_dat = cases[test_index,]
# controls_dat = betaControls
# valid_dat = betaValid
# bh_features = bh_feat_sig
# gender = T

runEnet <- function(training_dat, 
                    test_dat,
                    controls_dat,
                    controls_dat_old,
                    controls_dat_full,
                    valid_dat,
                    bh_features,
                    gender) 
{
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  bh_features <- append('M', bh_features)
  bh_features <- append('F', bh_features)
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  # # get y
  train_y <- as.numeric(training_dat$age_diagnosis)
  test_y <- as.numeric(test_dat$age_diagnosis)
  test_y_controls <- as.numeric(controls_dat$age_sample_collection)
  test_y_controls_old <- as.numeric(controls_dat_old$age_sample_collection)
  test_y_controls_full <- as.numeric(controls_dat_full$age_sample_collection)
  
  test_y_valid <- as.numeric(valid_dat$age_diagnosis)
  patient_age <- as.numeric(test_dat$age_sample_collection)
  missing_ind <- !is.na(patient_age)
  patient_age <- patient_age[missing_ind]
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  controls_dat_old <- controls_dat_old[, intersected_feats]
  controls_dat_full <- controls_dat_full[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
  N_CV_REPEATS = 2
  nfolds = 3
  
  ###### ENET
  # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
  # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
  elastic_net.cv_error = vector()
  elastic_net.cv_model = list()
  elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
  
  # set parameters for training model
  type_family <- 'gaussian'
  
  # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
  # or if you have a high number fo N_CV_REPEATS
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
                                                , y =  train_y
                                                , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                , type.measure = 'deviance'
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
                                     , type.measure = 'deviance'
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
                                   type = 'response', )
  
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  # get controls
  temp_test.predictions_controls <- predict(model,
                                            data.matrix(controls_dat),
                                            type = 'response')
  
  
  # look at s = 'lambda.min'
  test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
  
  # get controls old
  temp_test.predictions_controls_old <- predict(model,
                                                data.matrix(controls_dat_old),
                                                type = 'response')
  
  
  
  test.predictions_controls_old <- temp_test.predictions_controls_old[, temp.min_lambda_index]
  
  # get controls full
  temp_test.predictions_controls_full <- predict(model,
                                                 data.matrix(controls_dat_full),
                                                 type = 'response')
  
  
  
  test.predictions_controls_full <- temp_test.predictions_controls_full[, temp.min_lambda_index]
  
  # valid
  temp_test.predictions_valid <- predict(model,
                                         data.matrix(valid_dat),
                                         type = 'response')
  
  
  
  test.predictions_valid  <- temp_test.predictions_valid[, temp.min_lambda_index]
  
  importance <- coef(model)
  
  lambda_value <- elastic_net.cv_model$lambda.min
  
  cases_cor  <- cor(test_y, test.predictions)
  
  age_cor <- cor(test.predictions[missing_ind], patient_age)
  
  # # for each iteration, this should always be the same.
  controls_cor <- cor(test_y_controls, test.predictions_controls)
  controls_cor_old <- cor(test_y_controls_old, test.predictions_controls_old)
  controls_cor_full <- cor(test_y_controls_full, test.predictions_controls_full)
  
  valid_cor  <- cor(test_y_valid, test.predictions_valid)
  
  alpha  <- best_alpha
  
  return(list(alpha, lambda_value, importance, 
              cases_cor, 
              age_cor, 
              controls_cor, 
              controls_cor_full, 
              controls_cor_old, 
              valid_cor, 
              temp.non_zero_coeff))
  
  
}



##########
# enet diff
###########
# # # #
# training_dat = cases[train_index,]
# test_dat = cases[test_index,]
# controls_dat = betaControls
# valid_dat = betaValid
# bh_features = bh_feat_sig
# gender = T

runEnetDiff <- function(training_dat, 
                        test_dat,
                        bh_features,
                        gender) 
{
  
  if(gender) {
    # get intersection of bh features and real data
    bh_features <- as.character(unlist(bh_features))
    bh_features <- append('M', bh_features)
    bh_features <- append('F', bh_features)
  }
  
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  # # get y
  train_y <- as.numeric(training_dat$age_diagnosis)
  test_y <- as.numeric(test_dat$age_diagnosis)
  # get test age of sample collection
  patient_age <- as.numeric(test_dat$age_sample_collection)
  
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
  type_family <- 'gaussian'
  
  # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
  # or if you have a high number fo N_CV_REPEATS
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
                                                , y =  train_y
                                                , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                , type.measure = 'deviance'
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
                                     , type.measure = 'deviance'
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
  
  test.predictions <- temp_test.predictions[, temp.min_lambda_index] 
  
  
  return(list(test.predictions, test_y, patient_age))
  
  
}


##########
# enet diff
###########

runEnetRand <- function(training_dat,
                        controls_dat,
                        valid_dat,
                        test_dat,
                        bh_features,
                        age_cutoff,
                        gender) 
{
  
  if(gender) {
    # get intersection of bh features and real data
    bh_features <- as.character(unlist(bh_features))
    bh_features <- append('M', bh_features)
    bh_features <- append('F', bh_features)
  }
  
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  # # get y
  train_y <- as.numeric(training_dat$age_diagnosis)
  
  # other y and age
  valid_y <- as.numeric(valid_dat$age_diagnosis)
  test_y <- as.numeric(test_dat$age_diagnosis)
  # get test age of sample collection
  patient_age <- as.numeric(test_dat$age_sample_collection)
  patient_age_controls <- as.numeric(controls_dat$age_sample_collection)
  patient_age_valid <- as.numeric(valid_dat$age_sample_collection)
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
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
  type_family <- 'gaussian'
  
  # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
  # or if you have a high number fo N_CV_REPEATS
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
                                                , y =  train_y
                                                , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                , type.measure = 'deviance'
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
                                     , type.measure = 'deviance'
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
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  test.predictions <- temp_test.predictions[, temp.min_lambda_index] 
  
  # get controls
  temp_test.predictions_controls <- predict(model,
                                            data.matrix(controls_dat),
                                            type = 'response')
  
  
  
  test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
  
  
  # get validation
  temp_test.predictions_valid <- predict(model,
                                         data.matrix(valid_dat),
                                         type = 'response')
  
  
  test.predictions_valid  <- temp_test.predictions_valid[, temp.min_lambda_index]
  
  
  
  return(list(test.predictions, test_y, patient_age, 
              test.predictions_valid, valid_y, patient_age_valid,
              test.predictions_controls, patient_age_controls))
  
  
}

##########
# lasso rand
###########

runLassoRand <- function(training_dat,
                         controls_dat,
                         controls_dat_old,
                         controls_dat_full,
                         valid_dat,
                         test_dat,
                         bh_features,
                         gender) 
{
  
  if(gender) {
    # get intersection of bh features and real data
    bh_features <- as.character(unlist(bh_features))
    bh_features <- append('M', bh_features)
    bh_features <- append('F', bh_features)
  }
  
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  # # get y
  train_y <- as.numeric(training_dat$age_diagnosis)
  
  # other y and age
  valid_y <- as.numeric(valid_dat$age_diagnosis)
  test_y <- as.numeric(test_dat$age_diagnosis)
  # get test age of sample collection
  patient_age <- as.numeric(test_dat$age_sample_collection)
  patient_age_controls <- as.numeric(controls_dat$age_sample_collection)
  patient_age_controls_old <- as.numeric(controls_dat_old$age_sample_collection)
  patient_age_controls_full <- as.numeric(controls_dat_full$age_sample_collection)
  patient_age_valid <- as.numeric(valid_dat$age_sample_collection)
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  controls_dat_old <- controls_dat_old[, intersected_feats]
  controls_dat_full <- controls_dat_full[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
  test_dat <- test_dat[, intersected_feats]
  
  
  nfolds <- 5
  type_family <- 'gaussian'
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
                                     , y =  train_y
                                     , alpha = 1
                                     , type.measure = 'deviance'
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
                  ,alpha = 1
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  test.predictions <- temp_test.predictions[, temp.min_lambda_index] 
  
  # get controls
  temp_test.predictions_controls <- predict(model,
                                            data.matrix(controls_dat),
                                            type = 'response')
  
  
  
  test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
  
  # get controls
  temp_test.predictions_controls_old <- predict(model,
                                                data.matrix(controls_dat_old),
                                                type = 'response')
  
  
  
  test.predictions_controls_old <- temp_test.predictions_controls_old[, temp.min_lambda_index]
  
  # get controls full
  temp_test.predictions_controls_full <- predict(model,
                                                 data.matrix(controls_dat_full),
                                                 type = 'response')
  
  
  
  test.predictions_controls_full <- temp_test.predictions_controls_full[, temp.min_lambda_index]
  
  
  
  # get validation
  temp_test.predictions_valid <- predict(model,
                                         data.matrix(valid_dat),
                                         type = 'response')
  
  
  test.predictions_valid  <- temp_test.predictions_valid[, temp.min_lambda_index]
  
  
  
  return(list(test.predictions, test_y, patient_age, 
              test.predictions_valid, valid_y, patient_age_valid,
              test.predictions_controls, patient_age_controls,
              test.predictions_controls_full, patient_age_controls_full,
              test.predictions_controls_old, patient_age_controls_old))
  
  
}

##########
# lasso rand
###########

runRidgeRand <- function(training_dat,
                         controls_dat,
                         controls_dat_old,
                         controls_dat_full,
                         valid_dat,
                         test_dat,
                         bh_features,
                         gender) 
{
  
  if(gender) {
    # get intersection of bh features and real data
    bh_features <- as.character(unlist(bh_features))
    bh_features <- append('M', bh_features)
    bh_features <- append('F', bh_features)
  }
  
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  # # get y
  train_y <- as.numeric(training_dat$age_diagnosis)
  
  # other y and age
  valid_y <- as.numeric(valid_dat$age_diagnosis)
  test_y <- as.numeric(test_dat$age_diagnosis)
  # get test age of sample collection
  patient_age <- as.numeric(test_dat$age_sample_collection)
  patient_age_controls <- as.numeric(controls_dat$age_sample_collection)
  patient_age_controls_old <- as.numeric(controls_dat_old$age_sample_collection)
  patient_age_controls_full <- as.numeric(controls_dat_full$age_sample_collection)
  patient_age_valid <- as.numeric(valid_dat$age_sample_collection)
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  controls_dat_old <- controls_dat_old[, intersected_feats]
  controls_dat_full <- controls_dat_full[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
  test_dat <- test_dat[, intersected_feats]
  
  
  nfolds <- 5
  type_family <- 'gaussian'
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
                                     , y =  train_y
                                     , alpha = 0
                                     , type.measure = 'deviance'
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
                  ,alpha = 0
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  test.predictions <- temp_test.predictions[, temp.min_lambda_index] 
  
  # get controls
  temp_test.predictions_controls <- predict(model,
                                            data.matrix(controls_dat),
                                            type = 'response')
  
  
  
  test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
  
  # get controls
  temp_test.predictions_controls_old <- predict(model,
                                                data.matrix(controls_dat_old),
                                                type = 'response')
  
  
  
  test.predictions_controls_old <- temp_test.predictions_controls_old[, temp.min_lambda_index]
  
  # get controls full
  temp_test.predictions_controls_full <- predict(model,
                                                 data.matrix(controls_dat_full),
                                                 type = 'response')
  
  
  
  test.predictions_controls_full <- temp_test.predictions_controls_full[, temp.min_lambda_index]
  
  
  
  # get validation
  temp_test.predictions_valid <- predict(model,
                                         data.matrix(valid_dat),
                                         type = 'response')
  
  
  test.predictions_valid  <- temp_test.predictions_valid[, temp.min_lambda_index]
  
  
  
  return(list(test.predictions, test_y, patient_age, 
              test.predictions_valid, valid_y, patient_age_valid,
              test.predictions_controls, patient_age_controls,
              test.predictions_controls_full, patient_age_controls_full,
              test.predictions_controls_old, patient_age_controls_old))
  
  
}


##########
# enet diff
###########

runRfRand <- function(training_dat,
                      controls_dat,
                      controls_dat_old,
                      controls_dat_full,
                      valid_dat,
                      test_dat,
                      bh_features,
                      gender) 
{
  
  if(gender) {
    # get intersection of bh features and real data
    bh_features <- as.character(unlist(bh_features))
    bh_features <- append('M', bh_features)
    bh_features <- append('F', bh_features)
  }
  
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  # # get y
  train_y <- as.numeric(training_dat$age_diagnosis)
  
  # other y and age
  valid_y <- as.numeric(valid_dat$age_diagnosis)
  test_y <- as.numeric(test_dat$age_diagnosis)
  # get test age of sample collection
  patient_age <- as.numeric(test_dat$age_sample_collection)
  patient_age_controls <- as.numeric(controls_dat$age_sample_collection)
  patient_age_controls_old <- as.numeric(controls_dat_old$age_sample_collection)
  patient_age_controls_full <- as.numeric(controls_dat_full$age_sample_collection)
  patient_age_valid <- as.numeric(valid_dat$age_sample_collection)
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  controls_dat_old <- controls_dat_old[, intersected_feats]
  controls_dat_full <- controls_dat_full[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
  test_dat <- test_dat[, intersected_feats]
  
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = 4,      
    repeats = 1,
    allowParallel = TRUE)
  
  # mtry: Number of variables randomly sampled as candidates at each split.
  # ntree: Number of trees to grow.
  mtry <- sqrt(ncol(training_dat))
  tunegrid <- expand.grid(.mtry=mtry)
  
  model  <- train(x = training_dat
                  , y =train_y
                  , method = "rf"
                  , trControl = fitControl
                  , tuneGrid = tunegrid
                  , importance = T
                  , verbose = FALSE)
  
  temp <- varImp(model)[[1]]
  importance <- cbind(rownames(temp), temp$Overall)
  
  # predict on test data
  test.predictions <- predict(model,
                              newdata = test_dat)
  
  
  
  # get controls
  test.predictions_controls <- predict(model,
                                       controls_dat)
  
  
  
  # get controls old
  test.predictions_controls_old <- predict(model,
                                           controls_dat_old)
  
  
  # get controls
  test.predictions_controls_full <- predict(model,
                                            controls_dat_full)
  
  
  
  # get controls
  test.predictions_valid <- predict(model,
                                    valid_dat)
  
  
  
  
  return(list(test.predictions, test_y, patient_age, 
              test.predictions_valid, valid_y, patient_age_valid,
              test.predictions_controls, patient_age_controls,
              test.predictions_controls_full, patient_age_controls_full,
              test.predictions_controls_old, patient_age_controls_old))
  
  
}



##########
# enet diff
###########

runSvmRand <- function(training_dat,
                       controls_dat,
                       controls_dat_old,
                       controls_dat_full,
                       valid_dat,
                       test_dat,
                       bh_features,
                       gender) 
{
  
  if(gender) {
    # get intersection of bh features and real data
    bh_features <- as.character(unlist(bh_features))
    bh_features <- append('M', bh_features)
    bh_features <- append('F', bh_features)
  }
  
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  # # get y
  train_y <- as.numeric(training_dat$age_diagnosis)
  
  # other y and age
  valid_y <- as.numeric(valid_dat$age_diagnosis)
  test_y <- as.numeric(test_dat$age_diagnosis)
  # get test age of sample collection
  patient_age <- as.numeric(test_dat$age_sample_collection)
  patient_age_controls <- as.numeric(controls_dat$age_sample_collection)
  patient_age_controls_old <- as.numeric(controls_dat_old$age_sample_collection)
  patient_age_controls_full <- as.numeric(controls_dat_full$age_sample_collection)
  patient_age_valid <- as.numeric(valid_dat$age_sample_collection)
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  controls_dat_old <- controls_dat_old[, intersected_feats]
  controls_dat_full <- controls_dat_full[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
  test_dat <- test_dat[, intersected_feats]
  # 6) SVM with radial kernel
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = 4,      
    repeats = 1,
    allowParallel = TRUE)
  
  model <- train(x = training_dat,
                 , y = train_y
                 , method = "svmRadial"
                 , trControl = fitControl
                 , verbose = FALSE
  )
  # predict on test data
  test.predictions <- predict(model,
                              newdata = test_dat)
  
  
  
  # get controls
  test.predictions_controls <- predict(model,
                                       controls_dat)
  
  
  
  # get controls old
  test.predictions_controls_old <- predict(model,
                                           controls_dat_old)
  
  
  # get controls
  test.predictions_controls_full <- predict(model,
                                            controls_dat_full)
  
  
  
  # get controls
  test.predictions_valid <- predict(model,
                                    valid_dat)
  
  
  
  
  return(list(test.predictions, test_y, patient_age, 
              test.predictions_valid, valid_y, patient_age_valid,
              test.predictions_controls, patient_age_controls,
              test.predictions_controls_full, patient_age_controls_full,
              test.predictions_controls_old, patient_age_controls_old))
  
  
}



##########
# enet case
###########

runEnetCase <- function(cases_data, cases_y, alpha_number) 
{
  
  set.seed(alpha_number)
  # get a column for each dataset indicating the fold
  fold_vec <- sample(1:5, nrow(cases_data), replace = T)
  
  cases_cor <- list()
  alpha_score_list <- list()
  for (i in 1:5) {
    
    # get x 
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    # get traiind and test data
    training_dat <- cases_data[train_index,]
    training_y <- cases_y[train_index]
    
    testing_dat <- cases_data[!train_index,]
    testing_y <- cases_y[!train_index]
    
    # get training and test data
    N_CV_REPEATS = 2
    nfolds = 3
    
    ###### ENET
    # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
    # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
    elastic_net.cv_error = vector()
    elastic_net.cv_model = list()
    elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
    
    # set parameters for training model
    type_family <- 'gaussian'
    
    # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
    # or if you have a high number fo N_CV_REPEATS
    temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
      for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
      {      
        elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
                                                  , y =  training_y
                                                  , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                  , type.measure = 'deviance'
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
                                       , y =  training_y
                                       , alpha = elastic_net.ALPHA[temp.best_alpha_index]
                                       , type.measure = 'deviance'
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
                    , y =  training_y
                    ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                    ,standardize=FALSE
                    ,nlambda = 100
                    ,family = type_family)
    
    # This returns 100 prediction with 1-100 lambdas
    temp_test.predictions <- predict(model, 
                                     data.matrix(testing_dat),
                                     type = 'response')
    
    test.predictions <- temp_test.predictions[, temp.min_lambda_index]
    
    
    importance <- coef(model)
    
    lambda_value <- elastic_net.cv_model$lambda.min
    
    cases_cor[[i]]  <- cor(testing_y, test.predictions)
    
    alpha_score_list[[i]]  <- best_alpha
  }
  
  mean_cor <- mean(unlist(cases_cor))
  mean_alpha <- mean(unlist(alpha_score_list))
  
  return(list(mean_alpha, 
              mean_cor))
  
  
}




# training_dat = cases[train_index,]
# test_dat = cases[test_index,]
# controls_dat = betaControls
# valid_dat = betaValid
# bh_features = bh_feat_sig
# gender = T
# cutoff = 48

runEnetFac <- function(training_dat, 
                       test_dat,
                       controls_dat,
                       valid_dat,
                       bh_features,
                       gender,
                       cutoff) 
{
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  bh_features <- append('M', bh_features)
  bh_features <- append('F', bh_features)
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  
  # # get y
  train_y <- factor(ifelse(training_dat$age_diagnosis <= cutoff, 'yes', 'no'), levels = c('yes', 'no'))
  test_y <- factor(ifelse(test_dat$age_diagnosis <= cutoff, 'yes', 'no'), levels = c('yes', 'no'))
  
  
  # test_y_controls <- as.numeric(controls_dat$age_sample_collection)
  # test_y_valid <- as.numeric(valid_dat$age_diagnosis)
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  # controls_dat <- controls_dat[, intersected_feats]
  # valid_dat <- valid_dat[, intersected_feats]
  
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
                                   type = 'class')
  
  
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  test.predictions <- factor(test.predictions, levels = c('yes', 'no'))
  
  test_stats <- confusionMatrix(test_y, test.predictions)
  importance <- coef(model)
  
  lambda_value <- elastic_net.cv_model$lambda.min
  
  alpha  <- best_alpha
  
  return(list(alpha, lambda_value, importance, test_stats, model, temp.non_zero_coeff))
  
  
}


##########
# get the variation of the probe that is orthogonal to age of sample collection
##########


getResidual <- function(data,
                        bh_features) 
  
{
  
  # get feature intersection
  intersected_feats <- intersect(bh_features, colnames(data))
  
  # get genes 
  data <- data[, c("age_diagnosis", 
                   "age_sample_collection",
                   "cancer_diagnosis_diagnoses", 
                   "M",
                   "F",
                   intersected_feats)]
  
  probes <- colnames(data)[6:ncol(data)]
  
  resid <- list()
  
  for (i in 6:ncol(data)){
    
    temp <- data[, i]
    temp1 <- data$age_sample_collection
    
    resid[[i]] <- lm(temp ~ temp1)$residuals
    
    print(i)
    
  }
  
  resid <- do.call('cbind', resid)
  resid <- apply(resid, 2, function(x) as.numeric(x))
  resid <- as.data.frame(resid)
  resid <- cbind(data$age_diagnosis, 
                 data$age_sample_collection, 
                 data$cancer_diagnosis_diagnoses,
                 data$M, 
                 data$F,
                 resid)
  
  # change colnames
  colnames(resid) <- c('age_diagnosis', 
                       'age_sample_collection',
                       'cancer_diagnosis_diagnoses',
                       'gender',
                       probes)
  
  
  
  return(resid)
  
}



##########
# get results
##########

getResults <- function(result_list)
{
  
  
  results_final <- list(result_list[[1]],
                        result_list[[2]],
                        result_list[[3]],
                        result_list[[4]],
                        result_list[[5]],
                        result_list[[6]],
                        result_list[[7]],
                        result_list[[8]],
                        result_list[[9]],
                        result_list[[10]])
  
  
  return(results_final)
  
}


##########
# get results cancer
##########

getResultsCancer <- function(result_list)
{
  
  
  results_final <- list(result_list[[1]],
                        result_list[[2]],
                        result_list[[3]],
                        result_list[[4]],
                        result_list[[5]],
                        result_list[[6]])
  
  return(results_final)
  
}


##########
# function for finding probes with threshold differenceW
##########

get_diff_probes <-
  function(cases, 
           other_850,
           thresh) {
    
    temp_names <- list()
    
    for (i in 1:nrow(cases)) {
      
      temp.sample_cases <- as.data.frame(t(cases[i, 8:ncol(cases)]))
      temp.sample_other_850 <- as.data.frame(t(other_850[i, 8:ncol(other_850)]))
      
      temp_diff <- abs(temp.sample_cases - temp.sample_other_850)
      
      temp_names[[i]] <- rownames(temp_diff)[temp_diff > thresh]
      
      print(i)
      
    }
    
    probe_names <- unlist(temp_names)
    # probe_names_dup <- probe_names[!duplicated(probe_names)]
    
    return(probe_names)
    
  }


##########
# function for finding probes with threshold differenceW
##########

get_diff_probes_keep <-
  function(cases, 
           other_850,
           thresh) {
    
    temp_names <- list()
    
    for (i in 1:nrow(cases)) {
      
      temp.sample_cases <- as.data.frame(t(cases[i, 8:ncol(cases)]))
      temp.sample_other_850 <- as.data.frame(t(other_850[i, 8:ncol(other_850)]))
      
      temp_diff <- as.data.frame(abs(temp.sample_cases - temp.sample_other_850))
      
      colnames(temp_diff)[1] <- 'diff'
      temp_diff$sample <- i
      temp_diff$names <- rownames(temp_diff)
      
      temp_names[[i]] <- temp_diff[temp_diff$diff > thresh,]
      
      print(i)
      
    }
    
    probe_data <- do.call(rbind, temp_names)
    
    
    return(probe_data)
    
  }

##########
# estimate linear model and transform
##########

linearTransform <- function (cases_12, 
                             controls_12, 
                             controls_full) {
  
  probe_model <- list()
  probe_control_result <- list()
  
  for (i in 8:ncol(controls_full)) {
    
    control <- as.data.frame(controls_12[, i])
    cases <- as.data.frame(cases_12[, i])
    model_data <- data.frame(control = control, cases = cases)
    names(model_data) <- c('control', 'cases')
    probe_model[[i]] <- lm(cases ~ control, data = model_data)
    control_full <- as.numeric(controls_full[, i])
    model_data_new <- data.frame(control = control_full)
    names(model_data_new) <- 'control'
    probe_control_result[[i]] <- predict(probe_model[[i]], newdata = model_data_new, type = 'response')
    
    print(i) 
  }
  
  # transpose results
  temp <- do.call(rbind, probe_control_result)
  transform_controls <- t(temp)
  
  # add cg sites
  colnames(transform_controls) <- colnames(controls_full[8:ncol(controls_full)])
  
  # add clinical variables
  transform_controls <- as.data.frame(cbind(id = controls_full$id, 
                                            p53_germline = controls_full$p53_germline, 
                                            cancer_diagnosis_diagnoses = controls_full$cancer_diagnosis_diagnoses, 
                                            age_diagnosis = controls_full$age_diagnosis, 
                                            age_sample_collection = controls_full$age_sample_collection,
                                            gender = controls_full$gender, 
                                            sentrix_id = controls_full$sentrix_id, 
                                            transform_controls))
  
  # make numeric
  transform_controls[, 8:ncol(transform_controls)] <- apply(transform_controls[, 8:ncol(transform_controls)], 
                                                            2, 
                                                            function(x) as.numeric(x))
  
  
  return(transform_controls)
  
}

##########
# function to plot each id against the other
##########

plotCaseCon <- 
  function (cases, controls, row_index) {
    
    cases_cg <- as.numeric(cases[row_index, 8:ncol(cases)])
    controls_cg <- as.numeric(controls[row_index, 8:ncol(controls)])
    
    smoothScatter(cases_cg, 
                  controls_cg, 
                  main = paste0(row_index, '_', 'sample'),
                  xlab = 'cases', 
                  ylab = 'controls',
                  xlim = c(0, 10),
                  ylim = c(0,10))
    
  }

# subset to get methylation data 
get_model_dat <- function(full_data, probe_start, seed_num, k) {
  
  beta_methyl <- as.matrix(t(full_data[ ,probe_start:ncol(full_data)]))
  rownames(beta_methyl) <- NULL
  attributes(beta_methyl)[2] <- NULL
  
  clin_data <- full_data[1:(probe_start-1)]
  
  fold_vec <- getFolds(full_data, seed_number = seed_num, k = k)
  
  feat_names <- colnames(full_data)[probe_start:ncol(full_data)]
  
  
  return(list(beta_methyl, clin_data, fold_vec$folds, feat_names))
} 



#########
# predict
#########

superpc.predict <- function (object, data, newdata, threshold, n.components = 3, 
                             prediction.type = c("continuous", "discrete", "nonzero"), 
                             n.class = 2) 
{
  this.call <- match.call()
  prediction.type <- match.arg(prediction.type)
  if (n.class > 3) {
    stop("Maximum number of survival classes is 3")
  }
  
  # get features that have feature.scores larger than the threshold 
  which.features <- (abs(object$feature.scores) >= threshold)
  
  # get x training data with which.featurs (on row index because matrix is p x n)
  x.sml <- data$x[which.features, ]
  
  # get number of PCAs
  n.pc <- n.components
  
  # calculate pca
  x.sml.svd <- mysvd(x.sml, n.components = n.components)
  if (prediction.type == "nonzero") {
    if (!is.null(data$featurenames)) {
      out <- data$featurenames[which.features]
    }
    else {
      
      # get features to be used in PCA
      out <- (1:nrow(data$x))[which.features]
    }
  }
  
  #### HERE
  if (prediction.type == "continuous" | prediction.type == 
      "discrete") {
    
    # test x, with chosen features
    xtemp = newdata$x[which.features, ]
    
    # this is all to get your pca from the test data 
    # x temp 9 test data centered by the means from 
    # xtemp1 = t(scale(t(xtemp), center = x.sml.svd$feature.means, 
    #                 scale = F))
    scal = apply(scale(abs(x.sml.svd$u), center = F, scale = x.sml.svd$d), 
                 2, sum)
    
    # pca for test data
    cur.v <- scale(t(xtemp) %*% x.sml.svd$u, center = FALSE, 
                   scale = scal * x.sml.svd$d)
    
    # x train with chosen features
    xtemp0 = data$x[which.features, ]
    
    # center x train with mean features
    xtemp0 = t(scale(t(xtemp0), center = x.sml.svd$feature.means, 
                     scale = F))
    # get PCA for training data
    cur.v0 <- scale(t(xtemp0) %*% x.sml.svd$u, center = FALSE, 
                    scale = scal * x.sml.svd$d)
  }
  
  # training pca and training data
  result <- superpc.fit.to.outcome(object, data, cur.v0, print = FALSE)$results
  
  if (object$type == "survival") {
    coef = result$coef
  }
  if (object$type == "regression") {
    
    # training coefficient
    coef = result$coef[-1]
  }
  if (prediction.type == "continuous") {
    
    # test data - flip sign on negative PCA coefficients 
    out <- scale(cur.v, center = FALSE, scale = sign(coef))
    
    # test data - flip sign on negative PCA coefficients 
    out_train <- scale(cur.v0, center = FALSE, scale = sign(coef))
    
    # each row, multiply coefficient by pca
    v.pred.1df = apply(scale(out, center = FALSE, scale = 1/abs(coef)), 
                       1, sum)
  }
  else if (prediction.type == "discrete") {
    out0 <- scale(cur.v0, center = FALSE, scale = sign(coef))
    v.pred0.1df = apply(scale(out0, center = FALSE, scale = 1/abs(coef)), 
                        1, sum)
    out <- scale(cur.v, center = FALSE, scale = sign(coef))
    v.pred.1df = apply(scale(out, center = FALSE, scale = 1/abs(coef)), 
                       1, sum)
    for (j in 1:ncol(out)) {
      # br = quantile(cur.v0[, j], (0:n.class)/n.class)
      br = quantile(out0[, j], (0:n.class)/n.class) ## yp
      # out[, j] <- cut(out[, j], breaks = br, n.class, labels = FALSE)
      out[,j] = ifelse(out[,j] <= br[2], 1, 2) ## yp
      #  out[is.na(out[, j]), j] <- 1
    }
    br = quantile(v.pred0.1df, (0:n.class)/n.class)
    # v.pred.1df <- cut(v.pred.1df, breaks = br, labels = FALSE)
    # v.pred.1df[is.na(v.pred.1df)] <- 1
    v.pred.1df = ifelse(v.pred.1df <= br[2], 1, 2)  ## yp
  }
  if (is.matrix(out)) {
    dimnames(out) = list(NULL, rep(prediction.type, ncol(out)))
  }
  junk <- list(v.pred = out, v.train = out_train, u = x.sml.svd$u, d = x.sml.svd$d, 
               which.features = which.features, v.pred.1df = v.pred.1df, 
               n.components = n.pc, coef = result$coef, call = this.call, 
               prediction.type = prediction.type)
  return(junk)
}

##########
# fit to outcome function
##########
# superpc.fit.to.outcome(train.obj, data.test, fit.cts$v.pred)
# fit <- fit.cts
# data.test <- data.test
# score <- fit.cts$v.pred
superpc.fit.to.outcome<- function(fit, data.test,score, competing.predictors=NULL,  print=TRUE, iter.max=5){
  
  
  type=fit$type
  
  if(type=="survival"){temp.list=makelist(data.test$y, data.test$censoring.status, score)}
  if(type=="regression"){temp.list=makelist(data.test$y,NULL, score)}
  
  if(!is.null(competing.predictors)){
    temp.list=c(temp.list,competing.predictors)
  }
  
  
  if(type=="survival"){
    require(survival)
    results<-coxph(Surv(y, censoring.status)~., data=temp.list, control=coxph.control(iter.max=iter.max))
  }
  
  else{
    # fit test y against 3 pcas
    results<-lm(data.test$y~.,  data=temp.list)
  }
  
  
  if(print){print(summary(results))}
  
  
  ss=summary(results)
  if(type=="survival"){ test.stat=ss$logtest[1]
  df=ss$logtest[2]
  pvalue=ss$logtest[3]
  }
  if(type=="regression"){ test.stat=ss$fstat[1]
  df=ss$fstat[2:3]
  pvalue=1-pf(test.stat,df[1],df[2])
  }
  
  teststat.table=matrix(c(test.stat, df, pvalue), nrow=1)
  if(length(df)==1){dflabel="df"}
  if(length(df)==2){dflabel=c("df1", "df2")}
  
  dimnames(teststat.table)=list(NULL,c("test statistic",dflabel,"p-value"))
  
  
  return(list(results=results, teststat.table=teststat.table,  coeftable=ss$coef))
}

##########
# fit to outcome function cross validation
##########

superpc_fit_lm_cv <- 
  function(y, 
           score, 
           cv){
    
    # combine y and score
    mod_data <-as.data.frame(cbind(y, score))
    
    # fix colnames
    colnames(mod_data) <- c('y', 'pca1', 'pca2', 'pca3')
    
    # feat_names 
    feat_names <- colnames(mod_data)[-1]
    
    # list to strore predictions,ground truth and model
    test_y <- 
      y_pred <- 
      trained_model <- list()
    
    
    if (cv == 'loocv') {
      # perform leave one out cross validation 
      for(i in 1:nrow(mod_data)){
        
        # get loocv training and test data
        train_x_y <- mod_data[-i,]
        test_x <- mod_data[i, feat_names]
        test_y[[i]] <- mod_data$y[i] 
        
        # fit test y against 3 pcas
        trained_model[[i]] <- lm(y ~.,  train_x_y)
        
        # predict 
        y_pred[[i]] <- predict(trained_model[[i]], test_x)
        
      }
      onset_cor <- cor(unlist(y_pred), unlist(test_y))
      
      return(list(onset_cor, trained_model))
      
    }
    
    if (cv == 'k_fold') {
      
      # get fold vector
      mod_data$folds <- sample(c(1,2,3), nrow(mod_data), replace = T)
      
      for(i in unique(sort(mod_data$folds))) {
        
        # train_set and test_set index
        train_set <- !grepl(i, mod_data$folds)
        test_set <- !train_set
        
        # get training and testing data
        train_x_y <- mod_data[train_set, c('y', feat_names)]
        test_x <- mod_data[test_set, feat_names]
        test_y[[i]] <- mod_data$y[test_set] 
        
        # fit test y against 3 pcas
        trained_model[[i]] <- lm(y ~.,  train_x_y)
        
        # predict 
        y_pred[[i]] <- predict(trained_model[[i]], test_x)
        
      }
      
      onset_cor <- cor(unlist(y_pred), unlist(test_y))
      
      return(list(onset_cor, trained_model))
    }
    
    if (cv =='auto') {
      auto_cv <- cv.lm(mod_data, formula(y ~.), m = 2)
      
      temp <- do.call(rbind, auto_cv)
      onset_cor <- cor(temp$cvpred, temp$y)
      
      return(list(onset_cor, auto_cv))
      
    }
    
    
  }


##########
# function for listing objects
##########
makelist=function (y, censoring.status, predictors)
{
  val = list(y = y)
  if (!is.null(censoring.status)) {
    val$censoring.status = censoring.status
  }
  if (!is.matrix(predictors)) {
    val$score.1 = predictors
  }
  
  if (is.matrix(predictors)) {
    if (ncol(predictors) > 3) {
      stop("Can't have > 3 principal components")
    }
    predictor.type=dimnames(predictors)[[2]]
    
    if(is.null(dimnames(predictors)[[2]])){
      predictor.type=rep("continuous",ncol(predictors))
    }
    score1 = predictors[, 1]
    if(predictor.type[1]=="factor") {
      score1 = as.factor(score1)
    }
    val$score.1 = score1
    if (ncol(predictors) > 1) {
      score2 = predictors[, 2]
      if(predictor.type[2]=="factor") {
        score2 = as.factor(score2)
      }
      val$score.2 = score2
    }
    if (ncol(predictors) > 2) {
      score3 = predictors[, 3]
      if(predictor.type[3]=="factor") {
        score3 = as.factor(score3)
      }
      val$score.3 = score3
    }
  }
  return(val)
}


mysvd <- function (x, n.components = NULL) {
  p <- nrow(x)
  n <- ncol(x)
  feature.means <- rowMeans(x)
  x <- t(scale(t(x), center = feature.means, scale = F))
  if (is.null(n.components)) {
    n.components = min(n, p)
  }
  if (p > n) {
    a <- eigen(t(x) %*% x)
    v <- a$vec[, 1:n.components, drop = FALSE]
    d <- sqrt(a$val[1:n.components, drop = FALSE])
    u <- scale(x %*% v, center = FALSE, scale = d)
    return(list(u = u, d = d, v = v, feature.means = feature.means))
  }
  else {
    junk <- svd(x, LINPACK = TRUE)
    nc = min(ncol(junk$u), n.components)
    return(list(u = junk$u[, 1:nc], d = junk$d[1:nc], v = junk$v[, 
                                                                 1:nc], feature.means = feature.means))
  }
}


##########
# enet diff
###########

runEnetRandResid <- function(training_dat,
                             test_dat,
                             bh_features,
                             gender) 
{
  
  if(gender) {
    # get intersection of bh features and real data
    bh_features <- as.character(unlist(bh_features))
    bh_features <- append('M', bh_features)
    bh_features <- append('F', bh_features)
  }
  
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  # # get y
  train_y <- as.numeric(training_dat$age_diagnosis)
  
  test_y <- as.numeric(test_dat$age_diagnosis)
  # get test age of sample collection
  patient_age <- as.numeric(test_dat$age_sample_collection)
  
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
  type_family <- 'gaussian'
  
  # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
  # or if you have a high number fo N_CV_REPEATS
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
                                                , y =  train_y
                                                , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                , type.measure = 'deviance'
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
                                     , type.measure = 'deviance'
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
  
  test.predictions <- temp_test.predictions[, temp.min_lambda_index] 
  
  
  
  return(list(test.predictions, test_y, patient_age))
  
  
}

##########
# predict cancer
##########
predCancer <- function(training_dat, 
                       test_dat,
                       clin_dat,
                       bh_features,
                       gender,
                       p53) 
{
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  # # get y
  train_y <- as.factor(ifelse(!grepl('Unaffected',training_dat$cancer_diagnosis_diagnoses), 'a', 'b'))
  test_y <- as.factor(ifelse(!grepl('Unaffected',test_dat$cancer_diagnosis_diagnoses), 'a', 'b'))
  
  train_y <- factor(train_y, levels = c('a', 'b'))
  test_y <- factor(test_y, levels = c('a', 'b'))
  
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
                                   type = 'class')
  
  
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  test.predictions <- factor(test.predictions, levels = c('a', 'b'))
  
  test_stats <- caret::confusionMatrix(test_y, test.predictions)
  importance <- coef(model)
  
  # get prediction in clinical data
  clin_dat$prediction <- test.predictions
  
  return(list(importance, test_stats, clin_data))
  
  
}




##########
# enet diff
###########
# # 
# training_dat = m_train_cases
# controls_dat = m_controls_mod
# valid_dat = m_valid_mod
# test_dat = m_test_cases
# age_cutoff = 72
# bh_features = bh_features
# rand_feats = rand_feats
# gender = T

runEnetRandFac <- function(training_dat,
                           controls_dat,
                           valid_dat,
                           test_dat,
                           age_cutoff,
                           bh_features,
                           gender) {
  
  
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    
    intersected_feats <- append('M', intersected_feats)
    intersected_feats <- append('F', intersected_feats)
  }
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- ifelse(training_dat$age_diagnosis < age_cutoff, 1, 0)
  test_y <-  ifelse(test_dat$age_diagnosis < age_cutoff, 1, 0)
  valid_y <-  ifelse(valid_dat$age_diagnosis < age_cutoff, 1, 0)

  patient_age <-  ifelse(test_dat$age_sample_collection < age_cutoff, 1, 0)
  patient_age_controls <-  ifelse(controls_dat$age_sample_collection < age_cutoff, 1, 0)
  patient_age_valid <-  ifelse(valid_dat$age_sample_collection < age_cutoff, 1, 0)
  
  
  # get train and test clinical data
  training_clin <- training_dat[,3:11]
  test_clin <- test_dat[,3:11]
  controls_clin <- controls_dat[, 3:11]
  valid_clin <- valid_dat[, 3:11]

  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
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
  
  
  
  ##########
  # Predictions on test data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  
  # original should be fine, something wrong with caret package
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  # original should be fine, something wrong with caret package
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  temp_cases_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y))
  
  temp_cases_test <- cbind(test_clin, temp_cases_dat)
  

  # get controls
  temp_test.predictions_controls <- predict(model,
                                            data.matrix(controls_dat),
                                            type = 'class')
  
  
  test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
  
  temp_controls_dat <- as.data.frame(cbind(test_pred = test.predictions_controls, 
                                           test_label = patient_age_controls))
  
  temp_controls <- cbind(controls_clin, temp_controls_dat)
  
 
  ##########
  # Predictions on valid dat
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_valid <- predict(model, 
                                         data.matrix(valid_dat),
                                         type = 'class')
  
  
  # original should be fine, something wrong with caret package
  test.predictions_valid <- temp_test.predictions_valid[, temp.min_lambda_index]
 
  temp_valid_dat <- as.data.frame(cbind(test_pred = test.predictions_valid, 
                                           test_label = valid_y))
  
  temp_valid <- cbind(valid_clin, temp_valid_dat)
 
  ###########################################################################################
  
  return(list(temp_cases_test, temp_controls, temp_valid))
  
  
}




##########
# enet diff
###########

runLassoL1RandFac <- function(training_dat,
                              controls_dat,
                              valid_dat,
                              test_dat,
                              age_cutoff,
                              bh_features,
                              rand_feats,
                              gender) {
  
  
  
  if(gender) {
    # get intersection of bh features and real data
    bh_features <- as.character(unlist(bh_features))
    bh_features <- append('M', bh_features)
    bh_features <- append('F', bh_features)
  }
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- as.factor(ifelse(training_dat$age_diagnosis < age_cutoff, 'yes', 'no'))
  test_y <- as.factor(ifelse(test_dat$age_diagnosis < age_cutoff, 'yes', 'no'))
  valid_y <- as.factor(ifelse(valid_dat$age_diagnosis < age_cutoff, 'yes', 'no'))
  
  # get test age of sample collection
  patient_age <-  as.factor(ifelse(test_dat$age_sample_collection < age_cutoff, 'yes', 'no'))
  patient_age_controls <- as.factor(ifelse(controls_dat$age_sample_collection < age_cutoff, 'yes', 'no'))
  patient_age_valid <- as.factor(ifelse(valid_dat$age_sample_collection < age_cutoff, 'yes', 'no'))
  
  
  # get random features
  training_dat_rand <- training_dat[, intersected_feats_rand]
  controls_dat_rand <- controls_dat[, intersected_feats_rand]
  valid_dat_rand <- valid_dat[, intersected_feats_rand]
  test_dat_rand <- test_dat[,intersected_feats_rand]
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
  
  # start elastic net tuning
  N_CV_REPEATS = 2
  nfolds = 3
  
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  
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
                  ,alpha = 1
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  
  
  ##########
  # Predictions on test data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'class')
  
  
  # original should be fine, something wrong with caret package
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  # test.predictions <- ifelse(test.predictions >= pred_cutoff, 'yes', 'no')
  test.predictions <- factor(test.predictions, levels = c('yes', 'no'))
  test_y <- factor(test_y, levels = c('yes','no'))
  patient_age <- factor(patient_age, levels = c('yes', 'no'))
  test_stats <- caret::confusionMatrix(test.predictions, test_y)
  test_stats_age <- caret::confusionMatrix(test.predictions, patient_age)
  
  ##########
  # Predictions against controls 
  ##########
  
  # get controls
  temp_test.predictions_controls <- predict(model,
                                            data.matrix(controls_dat),
                                            type = 'class')
  
  
  test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
  # test.predictions_controls <- ifelse(test.predictions_controls >= pred_cutoff, 'yes', 'no')
  test.predictions_controls <- factor(test.predictions_controls, levels = c('yes', 'no'))
  patient_age_controls <- factor(patient_age_controls, levels = c('yes', 'no'))
  test_stats_controls <- caret::confusionMatrix(test.predictions_controls, patient_age_controls)
  
  
  ##########
  # Predictions on valid dat
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_valid <- predict(model, 
                                         data.matrix(valid_dat),
                                         type = 'class')
  
  
  # original should be fine, something wrong with caret package
  test.predictions_valid <- temp_test.predictions_valid[, temp.min_lambda_index]
  # test.predictions_valid <- ifelse(test.predictions_valid >=  pred_cutoff, 'yes', 'no')
  test.predictions_valid <- factor(test.predictions_valid, levels = c('yes', 'no'))
  valid_y <- factor(valid_y, levels = c('yes', 'no'))
  patient_age_valid <- factor(patient_age_valid, levels = c('yes', 'no'))
  test_stats_valid <- caret::confusionMatrix(test.predictions_valid, valid_y)
  test_stats_age_valid <- caret::confusionMatrix(test.predictions_valid, patient_age_valid)
  
  
  test_stats_feats <- test_stats
  test_stats_age_feats <- test_stats_age
  test_stats_controls_feats <- test_stats_controls
  test_stats_valid_feats <- test_stats_valid
  test_stats_age_valid_feats <- test_stats_age_valid
  
  ###########################################################################################
  # RAND
  
  
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat_rand)
                                     , y =  train_y
                                     , alpha = 1
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
  
  model  = glmnet(x = as.matrix(training_dat_rand)
                  , y =  train_y
                  ,alpha = 1
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  
  
  ##########
  # Predictions on test data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat_rand),
                                   type = 'class')
  
  
  # original should be fine, something wrong with caret package
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  # test.predictions <- ifelse(test.predictions >= pred_cutoff, 'yes', 'no')
  test.predictions <- factor(test.predictions, levels = c('yes', 'no'))
  test_y <- factor(test_y, levels = c('yes','no'))
  patient_age <- factor(patient_age, levels = c('yes', 'no'))
  test_stats <- caret::confusionMatrix(test.predictions, test_y)
  test_stats_age <- caret::confusionMatrix(test.predictions, patient_age)
  
  ##########
  # Predictions against controls 
  ##########
  
  # get controls
  temp_test.predictions_controls <- predict(model,
                                            data.matrix(controls_dat_rand),
                                            type = 'class')
  
  
  test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
  # test.predictions_controls <- ifelse(test.predictions_controls >= pred_cutoff, 'yes', 'no')
  test.predictions_controls <- factor(test.predictions_controls, levels = c('yes', 'no'))
  patient_age_controls <- factor(patient_age_controls, levels = c('yes', 'no'))
  test_stats_controls <- caret::confusionMatrix(test.predictions_controls, patient_age_controls)
  
  
  ##########
  # Predictions on valid dat
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_valid <- predict(model, 
                                         data.matrix(valid_dat_rand),
                                         type = 'class')
  
  
  # original should be fine, something wrong with caret package
  test.predictions_valid <- temp_test.predictions_valid[, temp.min_lambda_index]
  # test.predictions_valid <- ifelse(test.predictions_valid >=  pred_cutoff, 'yes', 'no')
  test.predictions_valid <- factor(test.predictions_valid, levels = c('yes', 'no'))
  valid_y <- factor(valid_y, levels = c('yes', 'no'))
  patient_age_valid <- factor(patient_age_valid, levels = c('yes', 'no'))
  test_stats_valid <- caret::confusionMatrix(test.predictions_valid, valid_y)
  test_stats_age_valid <- caret::confusionMatrix(test.predictions_valid, patient_age_valid)
  
  ###########################################################################################
  
  return(list(test_stats_feats, test_stats_age_feats,
              test_stats_controls_feats, test_stats_valid_feats,
              test_stats_age_valid_feats, 
              test_stats, test_stats_age,
              test_stats_controls, test_stats_valid,
              test_stats_age_valid))
  
  
}

##########
# enet diff
###########

run_enet_rand <- function(training_dat,
                          controls_dat,
                          valid_dat,
                          test_dat,
                          rand_feats,
                          age_cutoff) {
  
  
  
  
  
  intersected_feats <- intersect(rand_feats, colnames(training_dat))
  
  # # get y
  train_y <- as.factor(ifelse(training_dat$age_diagnosis < age_cutoff, 'yes', 'no'))
  test_y <- as.factor(ifelse(test_dat$age_diagnosis < age_cutoff, 'yes', 'no'))
  valid_y <- as.factor(ifelse(valid_dat$age_diagnosis < age_cutoff, 'yes', 'no'))
  
  # get test age of sample collection
  patient_age <-  as.factor(ifelse(test_dat$age_sample_collection < age_cutoff, 'yes', 'no'))
  patient_age_controls <- as.factor(ifelse(controls_dat$age_sample_collection < age_cutoff, 'yes', 'no'))
  patient_age_valid <- as.factor(ifelse(valid_dat$age_sample_collection < age_cutoff, 'yes', 'no'))
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
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
  
  
  
  ##########
  # Predictions on test data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'class')
  
  
  # original should be fine, something wrong with caret package
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  # test.predictions <- ifelse(test.predictions >= pred_cutoff, 'yes', 'no')
  test.predictions <- factor(test.predictions, levels = c('yes', 'no'))
  test_y <- factor(test_y, levels = c('yes','no'))
  patient_age <- factor(patient_age, levels = c('yes', 'no'))
  test_stats <- caret::confusionMatrix(test.predictions, test_y)
  test_stats_age <- caret::confusionMatrix(test.predictions, patient_age)
  
  ##########
  # Predictions against controls 
  ##########
  
  # get controls
  temp_test.predictions_controls <- predict(model,
                                            data.matrix(controls_dat),
                                            type = 'class')
  
  
  test.predictions_controls <- temp_test.predictions_controls[, temp.min_lambda_index]
  # test.predictions_controls <- ifelse(test.predictions_controls >= pred_cutoff, 'yes', 'no')
  test.predictions_controls <- factor(test.predictions_controls, levels = c('yes', 'no'))
  patient_age_controls <- factor(patient_age_controls, levels = c('yes', 'no'))
  test_stats_controls <- caret::confusionMatrix(test.predictions_controls, patient_age_controls)
  
  
  ##########
  # Predictions on valid dat
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_valid <- predict(model, 
                                         data.matrix(valid_dat),
                                         type = 'class')
  
  
  # original should be fine, something wrong with caret package
  test.predictions_valid <- temp_test.predictions_valid[, temp.min_lambda_index]
  # test.predictions_valid <- ifelse(test.predictions_valid >=  pred_cutoff, 'yes', 'no')
  test.predictions_valid <- factor(test.predictions_valid, levels = c('yes', 'no'))
  valid_y <- factor(valid_y, levels = c('yes', 'no'))
  patient_age_valid <- factor(patient_age_valid, levels = c('yes', 'no'))
  test_stats_valid <- caret::confusionMatrix(test.predictions_valid, valid_y)
  test_stats_age_valid <- caret::confusionMatrix(test.predictions_valid, patient_age_valid)
  
  
  return(list(test_stats, test_stats_age,
              test_stats_controls, test_stats_valid,
              test_stats_age_valid))
  
  
}




##########
# Lasso
###########
# training_dat = beta_cases[train_index,]
# controls_dat = beta_controls
# valid_dat = beta_valid
# test_dat = beta_cases[test_index,]
# bh_features = mod_feats

runGlmLassoRandFac <- function(training_dat,
                               controls_dat,
                               valid_dat,
                               test_dat,
                               age_cutoff,
                               pred_cutoff,
                               bh_features,
                               gender) {
  
  
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- bh_features
  
  # # get y
  training_dat$age_diagnosis <- ifelse(training_dat$age_diagnosis > age_cutoff, 0,1)
  
  
  test_y <-ifelse(test_dat$age_diagnosis > age_cutoff, 'no', 'yes')
  valid_y <- ifelse(valid_dat$age_diagnosis > age_cutoff, 'no', 'yes')
  
  test_y <- factor(test_y, levels = c('yes', 'no'))
  valid_y <-factor(valid_y, levels = c('yes', 'no'))
  
  # get test age of sample collection
  patient_age  <-  ifelse(test_dat$age_sample_collection > age_cutoff, 'no', 'yes')
  patient_age_controls<-  ifelse(controls_dat$age_sample_collection > age_cutoff,'no', 'yes')
  patient_age_valid <- ifelse(valid_dat$age_sample_collection > age_cutoff, 'no', 'yes')
  
  patient_age <- factor(patient_age, levels = c('yes', 'no'))
  patient_age_controls <- factor(patient_age_controls, levels = c('yes', 'no'))
  patient_age_valid <- factor(patient_age_valid, levels = c('yes', 'no'))
  
  # get bumphunter features
  training_dat <- training_dat[, c('age_diagnosis', 'family_name', intersected_feats)]
  controls_dat <- controls_dat[, c('age_sample_collection', 'family_name', intersected_feats)]
  valid_dat <- valid_dat[, c('age_diagnosis', 'family_name', intersected_feats)]
  test_dat <- test_dat[, c('age_diagnosis', 'family_name', intersected_feats)]
  
  # rename sample to diagnosis
  colnames(controls_dat)[1] <- 'age_diagnosis'
  
  # make factor variabls 
  training_dat$family_name <- as.factor(training_dat$family_name)  
  
  fmla <- as.formula(paste("age_diagnosis ~ ", paste(colnames(training_dat[,-c(1:2)]), collapse= "+")))
  
  ##########
  # linear mixed model
  ##########
  model <- glmmLasso(fmla,  
                     rnd = list(family_name = ~1), 
                     lambda=10, 
                     data = training_dat, family=binomial(link = "logit"))
  
  
  ##########
  # test dat
  ##########
  temp.pred <- predict(model, 
                       test_dat, 
                       type = 'response')
  
  test.pred <- ifelse(temp.pred > pred_cutoff, 'no', 'yes') 
  test.predictions <- factor(test.pred, levels = c('yes', 'no'))
  
  test_stats <- confusionMatrix(test_y, test.predictions)
  test_stats_age <- confusionMatrix(patient_age, test.predictions)
  
  
  ##########
  # controls dat
  ##########
  temp.pred <- predict(model, 
                       controls_dat, 
                       type = 'response')
  
  test.pred <- ifelse(temp.pred > pred_cutoff, 'no', 'yes') 
  test.predictions <- factor(test.pred, levels = c('yes', 'no'))
  
  test_stats <- confusionMatrix(test.predictions, patient_age_controls)
  
  
  
  return(list(test_stats, test_stats_age,
              test_stats_controls, test_stats_valid,
              test_stats_age_valid))
  
  
}

##########
# RF fac
###########


runRfRandFac <- function(training_dat,
                         controls_dat,
                         valid_dat,
                         test_dat,
                         age_cutoff,
                         bh_features,
                         rand_feats,
                         pred_cutoff,
                         gender) {
  if(gender) {
    # get intersection of bh features and real data
    bh_features <- as.character(unlist(bh_features))
    bh_features <- append('M', bh_features)
    bh_features <- append('F', bh_features)
  }
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- as.factor(ifelse(training_dat$age_diagnosis < age_cutoff, 'yes', 'no'))
  test_y <- as.factor(ifelse(test_dat$age_diagnosis < age_cutoff, 'yes', 'no'))
  valid_y <- as.factor(ifelse(valid_dat$age_diagnosis < age_cutoff, 'yes', 'no'))
  
  # get test age of sample collection
  patient_age <-  as.factor(ifelse(test_dat$age_sample_collection < age_cutoff, 'yes', 'no'))
  patient_age_controls <- as.factor(ifelse(controls_dat$age_sample_collection < age_cutoff, 'yes', 'no'))
  patient_age_valid <- as.factor(ifelse(valid_dat$age_sample_collection < age_cutoff, 'yes', 'no'))
  
  # get bumphunter features
  training_dat_rand <- training_dat[, intersected_feats_rand]
  controls_dat_rand <- controls_dat[, intersected_feats_rand]
  valid_dat_rand <- valid_dat[, intersected_feats_rand]
  test_dat_rand <- test_dat[, intersected_feats_rand]
  
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
  summaryFunc <- twoClassSummary
  
  NFOLDS = 4
  # determines how you train the model.
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = min(10, NFOLDS), 
    classProbs = TRUE,     
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = summaryFunc)
  
  
  # mtry: Number of variables randomly sampled as candidates at each split.
  # ntree: Number of trees to grow.
  mtry <- sqrt(ncol(training_dat))
  tunegrid <- expand.grid(.mtry=mtry)
  
  model  <- train(x = training_dat
                  , y =train_y
                  , method = "rf"
                  , trControl = fitControl
                  , tuneGrid = tunegrid
                  , importance = T
                  , metric = "ROC"
                  , verbose = FALSE)
  
  temp <- varImp(model)[[1]]
  importance <- cbind(rownames(temp), temp$Overall)
  
  ##########
  # Predictions on test data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions <- predict(model, 
                              data.matrix(test_dat),
                              type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions <- test.predictions$yes
  test.predictions <- as.factor(ifelse(test.predictions >= pred_cutoff, 'yes', 'no'))
  test.predictions <- factor(test.predictions, levels = c('yes', 'no'))
  test_y <- factor(test_y, levels = c('yes', 'no'))
  patient_age <- factor(patient_age, levels = c('yes', 'no'))
  
  test_stats <- caret::confusionMatrix(test.predictions, test_y)
  test_stats_age <- caret::confusionMatrix(test.predictions, patient_age)
  
  ##########
  # Predictions on controls
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions_controls <- predict(model, 
                                       data.matrix(controls_dat),
                                       type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions_controls <- test.predictions_controls$yes
  test.predictions_controls <- as.factor(ifelse(test.predictions_controls >= pred_cutoff, 'yes', 'no'))
  test.predictions_controls <- factor(test.predictions_controls, levels = c('yes', 'no'))
  patient_age_controls <- factor(patient_age_controls, levels = c('yes', 'no'))
  test_stats_controls <-caret:: confusionMatrix(test.predictions_controls, patient_age_controls)
  
  ##########
  # Predictions on valid data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions_valid <- predict(model, 
                                    data.matrix(valid_dat),
                                    type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions_valid <- test.predictions_valid$yes
  test.predictions_valid <- as.factor(ifelse(test.predictions_valid >= pred_cutoff, 'yes', 'no'))
  test.predictions_valid <- factor(test.predictions_valid, levels = c('yes', 'no'))
  valid_y <- factor(valid_y, levels = c('yes', 'no'))
  patient_age_valid <- factor(patient_age_valid, levels = c('yes', 'no'))
  
  
  test_stats_valid <- caret::confusionMatrix(test.predictions_valid, valid_y)
  test_stats_age_valid <- caret::confusionMatrix(test.predictions_valid, patient_age_valid)
  
  test_stats_norm <- test_stats
  test_stats_age_norm <- test_stats_age
  test_stats_controls_norm <- test_stats_controls
  test_stats_valid_norm <- test_stats_valid
  test_stats_age_valid_norm <- test_stats_age_valid
  
  ####################################################################################################
  # random 
  
  # mtry: Number of variables randomly sampled as candidates at each split.
  # ntree: Number of trees to grow.
  mtry <- sqrt(ncol(training_dat))
  tunegrid <- expand.grid(.mtry=mtry)
  
  model  <- train(x = training_dat_rand
                  , y =train_y
                  , method = "rf"
                  , trControl = fitControl
                  , tuneGrid = tunegrid
                  , importance = T
                  , metric = "ROC"
                  , verbose = FALSE)
  
  temp <- varImp(model)[[1]]
  importance <- cbind(rownames(temp), temp$Overall)
  
  ##########
  # Predictions on test data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions <- predict(model, 
                              data.matrix(test_dat_rand),
                              type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions <- test.predictions$yes
  test.predictions <- as.factor(ifelse(test.predictions >= pred_cutoff, 'yes', 'no'))
  test.predictions <- factor(test.predictions, levels = c('yes', 'no'))
  test_y <- factor(test_y, levels = c('yes', 'no'))
  patient_age <- factor(patient_age, levels = c('yes', 'no'))
  
  test_stats <- caret::confusionMatrix(test.predictions, test_y)
  test_stats_age <- caret::confusionMatrix(test.predictions, patient_age)
  
  ##########
  # Predictions on controls
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions_controls <- predict(model, 
                                       data.matrix(controls_dat_rand),
                                       type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions_controls <- test.predictions_controls$yes
  test.predictions_controls <- as.factor(ifelse(test.predictions_controls >= pred_cutoff, 'yes', 'no'))
  test.predictions_controls <- factor(test.predictions_controls, levels = c('yes', 'no'))
  patient_age_controls <- factor(patient_age_controls, levels = c('yes', 'no'))
  test_stats_controls <-caret:: confusionMatrix(test.predictions_controls, patient_age_controls)
  
  ##########
  # Predictions on valid data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions_valid <- predict(model, 
                                    data.matrix(valid_dat_rand),
                                    type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions_valid <- test.predictions_valid$yes
  test.predictions_valid <- as.factor(ifelse(test.predictions_valid >= pred_cutoff, 'yes', 'no'))
  test.predictions_valid <- factor(test.predictions_valid, levels = c('yes', 'no'))
  valid_y <- factor(valid_y, levels = c('yes', 'no'))
  patient_age_valid <- factor(patient_age_valid, levels = c('yes', 'no'))
  
  
  test_stats_valid <- caret::confusionMatrix(test.predictions_valid, valid_y)
  test_stats_age_valid <- caret::confusionMatrix(test.predictions_valid, patient_age_valid)
  
  return(list(test_stats_norm, test_stats_age_norm,
              test_stats_controls_norm, test_stats_valid_norm,
              test_stats_age_valid_norm,
              test_stats, test_stats_age,
              test_stats_controls, test_stats_valid,
              test_stats_age_valid))
}






##########
# RF fac
###########

run_rf_rand <- function(training_dat,
                        controls_dat,
                        valid_dat,
                        test_dat,
                        rand_feats,
                        age_cutoff) {
  
  
  intersected_feats <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- as.factor(ifelse(training_dat$age_diagnosis < age_cutoff, 'yes', 'no'))
  test_y <- as.factor(ifelse(test_dat$age_diagnosis < age_cutoff, 'yes', 'no'))
  valid_y <- as.factor(ifelse(valid_dat$age_diagnosis < age_cutoff, 'yes', 'no'))
  
  # get test age of sample collection
  patient_age <-  as.factor(ifelse(test_dat$age_sample_collection < age_cutoff, 'yes', 'no'))
  patient_age_controls <- as.factor(ifelse(controls_dat$age_sample_collection < age_cutoff, 'yes', 'no'))
  patient_age_valid <- as.factor(ifelse(valid_dat$age_sample_collection < age_cutoff, 'yes', 'no'))
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
  summaryFunc <- twoClassSummary
  
  NFOLDS = 4
  # determines how you train the model.
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = min(10, NFOLDS), 
    classProbs = TRUE,     
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = summaryFunc)
  
  
  # mtry: Number of variables randomly sampled as candidates at each split.
  # ntree: Number of trees to grow.
  mtry <- sqrt(ncol(training_dat))
  tunegrid <- expand.grid(.mtry=mtry)
  
  model  <- train(x = training_dat
                  , y =train_y
                  , method = "rf"
                  , trControl = fitControl
                  , tuneGrid = tunegrid
                  , importance = T
                  , metric = "ROC"
                  , verbose = FALSE)
  
  temp <- varImp(model)[[1]]
  importance <- cbind(rownames(temp), temp$Overall)
  
  
  ##########
  # Predictions on test data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions <- predict(model, 
                              data.matrix(test_dat),
                              type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions <- test.predictions$yes
  test.predictions <- as.factor(ifelse(test.predictions >= 0.5, 'yes', 'no'))
  test.predictions <- factor(test.predictions, levels = c('yes', 'no'))
  test_y <- factor(test_y, levels = c('yes', 'no'))
  patient_age <- factor(patient_age, levels = c('yes', 'no'))
  
  test_stats <- caret::confusionMatrix(test.predictions, test_y)
  test_stats_age <- caret::confusionMatrix(test.predictions, patient_age)
  
  ##########
  # Predictions on controls
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions_controls <- predict(model, 
                                       data.matrix(controls_dat),
                                       type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions_controls <- test.predictions_controls$yes
  test.predictions_controls <- as.factor(ifelse(test.predictions_controls >= 0.5, 'yes', 'no'))
  test.predictions_controls <- factor(test.predictions_controls, levels = c('yes', 'no'))
  patient_age_controls <- factor(patient_age_controls, levels = c('yes', 'no'))
  test_stats_controls <-caret:: confusionMatrix(test.predictions_controls, patient_age_controls)
  
  ##########
  # Predictions on valid data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions_valid <- predict(model, 
                                    data.matrix(valid_dat),
                                    type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions_valid <- test.predictions_valid$yes
  test.predictions_valid <- as.factor(ifelse(test.predictions_valid >= 0.5, 'yes', 'no'))
  test.predictions_valid <- factor(test.predictions_valid, levels = c('yes', 'no'))
  valid_y <- factor(valid_y, levels = c('yes', 'no'))
  patient_age_valid <- factor(patient_age_valid, levels = c('yes', 'no'))
  
  
  test_stats_valid <- caret::confusionMatrix(test.predictions_valid, valid_y)
  test_stats_age_valid <- caret::confusionMatrix(test.predictions_valid, patient_age_valid)
  
  
  
  return(list(test_stats, test_stats_age,
              test_stats_controls, test_stats_valid,
              test_stats_age_valid))
}



##########
# enet diff
###########

runSvmRandFac <- function(training_dat,
                          controls_dat,
                          valid_dat,
                          test_dat,
                          age_cutoff,
                          pred_cutoff,
                          bh_features,
                          gender) 
{
  
  
  if(gender) {
    # get intersection of bh features and real data
    bh_features <- as.character(unlist(bh_features))
    bh_features <- append('M', bh_features)
    bh_features <- append('F', bh_features)
  }
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  # # get y
  train_y <- ifelse(training_dat$age_diagnosis > age_cutoff, 'no', 'yes')
  test_y <-ifelse(test_dat$age_diagnosis > age_cutoff, 'no', 'yes')
  valid_y <- ifelse(valid_dat$age_diagnosis > age_cutoff, 'no', 'yes')
  
  train_y <- factor(train_y, levels = c('yes', 'no'))
  test_y <- factor(test_y, levels = c('yes', 'no'))
  valid_y <-factor(valid_y, levels = c('yes', 'no'))
  
  # get test age of sample collection
  patient_age <-  ifelse(test_dat$age_sample_collection > age_cutoff, 'no', 'yes')
  patient_age_controls <-  ifelse(controls_dat$age_sample_collection > age_cutoff, 'no', 'yes')
  patient_age_valid <- ifelse(valid_dat$age_sample_collection > age_cutoff, 'no', 'yes')
  
  patient_age <- factor(patient_age, levels = c('yes', 'no'))
  patient_age_controls <- factor(patient_age_controls, levels = c('yes', 'no'))
  patient_age_valid <- factor(patient_age_valid, levels = c('yes', 'no'))
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
  summaryFunc <- twoClassSummary
  
  NFOLDS = 4
  # determines how you train the model.
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = min(10, NFOLDS), 
    classProbs = TRUE,     
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = summaryFunc)
  
  
  
  model <- train(x = training_dat
                 , y = train_y
                 , method = "svmRadial"
                 , trControl = fitControl
                 , verbose = FALSE)
  
  ##########
  # Predictions on test data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions <- predict(model, 
                              data.matrix(test_dat),
                              type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions <- test.predictions$yes
  test.predictions <- as.factor(ifelse(test.predictions >= pred_cutoff, 'yes', 'no'))
  test.predictions <- factor(test.predictions, levels = c('yes', 'no'))
  
  test_stats <- confusionMatrix(test_y, test.predictions)
  test_stats_age <- confusionMatrix(patient_age, test.predictions)
  
  
  ##########
  # Predictions on controls
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions_controls <- predict(model, 
                                       data.matrix(controls_dat),
                                       type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions_controls <- test.predictions_controls$yes
  test.predictions_controls <- as.factor(ifelse(test.predictions_controls >= pred_cutoff, 'yes', 'no'))
  test.predictions_controls <- factor(test.predictions_controls, levels = c('yes', 'no'))
  
  test_stats_controls <- confusionMatrix(patient_age_controls, test.predictions_controls)
  
  ##########
  # Predictions on valid data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions_valid <- predict(model, 
                                    data.matrix(valid_dat),
                                    type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions_valid <- test.predictions_valid$yes
  test.predictions_valid <- as.factor(ifelse(test.predictions_valid >= pred_cutoff, 'yes', 'no'))
  test.predictions_valid <- factor(test.predictions_valid, levels = c('yes', 'no'))
  
  test_stats_valid <- confusionMatrix(valid_y, test.predictions_valid)
  test_stats_age_valid <- confusionMatrix(patient_age_valid, test.predictions_valid)
  
  
  
  return(list(test_stats, test_stats_age,
              test_stats_controls, test_stats_valid,
              test_stats_age_valid))
  
}



##########
# mixed lasso family as random effects
###########

run_mixed_lasso <- function(training_dat,
                            controls_dat,
                            valid_dat,
                            test_dat,
                            age_cutoff,
                            pred_cutoff,
                            bh_features,
                            gender) 
{
  
  
  if(gender) {
    # get intersection of bh features and real data
    bh_features <- as.character(unlist(bh_features))
    bh_features <- append('M', bh_features)
    bh_features <- append('F', bh_features)
  }
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  # # get y
  train_y <- ifelse(training_dat$age_diagnosis > age_cutoff, 0, 1)
  test_y <-ifelse(test_dat$age_diagnosis > age_cutoff, 0, 1)
  valid_y <- ifelse(valid_dat$age_diagnosis > age_cutoff, 0, 1)
  
  # train_y <- factor(train_y, levels = c('yes', 'no'))
  # test_y <- factor(test_y, levels = c('yes', 'no'))
  # valid_y <-factor(valid_y, levels = c('yes', 'no'))
  
  # get test age of sample collection
  patient_age <-  ifelse(test_dat$age_sample_collection > age_cutoff, 0, 1)
  patient_age_controls <-  ifelse(controls_dat$age_sample_collection > age_cutoff, 0, 1)
  patient_age_valid <- ifelse(valid_dat$age_sample_collection > age_cutoff, 0, 1)
  # 
  # patient_age <- factor(patient_age, levels = c('yes', 'no'))
  # patient_age_controls <- factor(patient_age_controls, levels = c('yes', 'no'))
  # patient_age_valid <- factor(patient_age_valid, levels = c('yes', 'no'))
  
  # get bumphunter features
  training_dat <- training_dat[, c('family_name', intersected_feats)]
  controls_dat <- controls_dat[, c('family_name', intersected_feats)]
  valid_dat <- valid_dat[, c('family_name', intersected_feats)]
  test_dat <- test_dat[, c('family_name', intersected_feats)]
  
  # recode family name (our groups) as numeric
  training_dat$family_name <- as.numeric(as.factor(training_dat$family_name))
  controls_dat$family_name <- as.numeric(as.factor(controls_dat$family_name))
  valid_dat$family_name <- as.numeric(as.factor(valid_dat$family_name))
  test_dat$family_name <- as.numeric(as.factor(test_dat$family_name))
  
  # add an intercept and make it a matrix
  training_dat <- as.matrix(cbind(1, training_dat))
  test_dat <- as.matrix(cbind(1, test_dat))
  controls_dat <- as.matrix(cbind(1, controls_dat))
  valid_dat <- as.matrix(cbind(1, valid_dat))
  
  # get random effects matrix which is just first two columns of data
  z_train <- as.matrix(training_dat[, 1])
  
  # get grp variable length of data
  g_train <- t(factor(training_dat[,2]))
  
  fit <-lassop(training_dat, 
               train_y, 
               z_train, 
               g_train, 
               mu=1.5, 
               fix=1, 
               rand=1)
  
  predict(fit)
  
  
  #HERE
  ########
  #independent random effects
  # x is a numeric matrix n*p (120, 81)  = training_dat (67*p)
  # y is outcome, length 120 = train_y (factor variable)
  # z is random effects matrix n*q (120, 2) 
  # grp variable length n
  # rand where z is in x
  # fix variables (in front of data) not submitted for selection - use 1 or 2
  rand
  fit=lassop(x,y,z,grp,D=1,mu=0.2,fix=1,rand=c(1,2))
  
  ##########
  # Predictions on test data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions <- predict(model, 
                              data.matrix(test_dat),
                              type = 'prob')
  
  
  # original should be fine, something wrong with caret package
  test.predictions <- test.predictions$yes
  test.predictions <- as.factor(ifelse(test.predictions >= pred_cutoff, 'yes', 'no'))
  test.predictions <- factor(test.predictions, levels = c('yes', 'no'))
  
  test_stats <- caret::confusionMatrix(test_y, test.predictions)
  test_stats_age <- caret::confusionMatrix(patient_age, test.predictions)
  
  
  
  
  
  return(list(test_stats, test_stats_age,
              test_stats_controls, test_stats_valid,
              test_stats_age_valid))
  
}

# mod_results_list <- mod_result
# dims_of_dat <- 22
# feat_name <- 'bh'
# mod_name <- 'rf'
# seed_number <- 1
# mod_results_list <- mod_result
get_class_results <- function(mod_results_list, dims_of_dat, mod_name, feat_name, seed_number) {
  
  class_onset_norm <- as.data.frame(t(mod_results_list[[1]]$byClass))
  class_onset_norm$age_type <- 'cases_onset'
  class_onset_norm$seed_number <- seed_number
  class_onset_norm$feature_num <- dims_of_dat
  class_onset_norm$model_method <-  mod_name
  class_onset_norm$feat_name <-  feat_name
  
  
  # cases age
  class_age_norm <-as.data.frame(t(mod_results_list[[2]]$byClass))
  class_age_norm$age_type <- 'cases_age'
  class_age_norm$seed_number <- seed_number
  class_age_norm$feature_num <- dims_of_dat
  class_age_norm$model_method <-  mod_name
  class_age_norm$feat_name <-  feat_name
  
  # valid age
  class_controls_age_norm <- as.data.frame(t(mod_results_list[[3]]$byClass))
  class_controls_age_norm$age_type <- 'controls_age'
  class_controls_age_norm$seed_number <- seed_number
  class_controls_age_norm$feature_num <- dims_of_dat
  class_controls_age_norm$model_method <- mod_name
  class_controls_age_norm$feat_name <-  feat_name
  
  # valid onset
  class_valid_onset_norm <- as.data.frame(t(mod_results_list[[4]]$byClass))
  class_valid_onset_norm$age_type <- 'valid_onset'
  class_valid_onset_norm$seed_number <- seed_number
  class_valid_onset_norm$feature_num <- dims_of_dat
  class_valid_onset_norm$model_method <-  mod_name
  class_valid_onset_norm$feat_name <-  feat_name
  
  
  # valid age
  class_valid_age_norm <- as.data.frame(t(mod_results_list[[5]]$byClass))
  class_valid_age_norm$age_type <- 'valid_age'
  class_valid_age_norm$seed_number <- seed_number
  class_valid_age_norm$feature_num <- dims_of_dat
  class_valid_age_norm$model_method <- mod_name
  class_valid_age_norm$feat_name <-  feat_name
  
  # get clinical data
  clin_data <- mod_results_list[[6]]
  
  # collapse
  clin_data_full <- do.call(rbind, clin_data)
  
  # combine class
  class_results_norm <- rbind(class_onset_norm,
                              class_age_norm,
                              class_valid_age_norm,
                              class_valid_onset_norm, 
                              class_controls_age_norm)
  
  

  return(list(class_results_norm, clin_data_full))
}

get_class_results_cancer <- function(mod_results_list, dims_of_dat, mod_name, gender, p53) {
  
  
  # save importance each time 
  
  
  # get stat results
  class_cancer <- as.data.frame(t(mod_results_list[[2]]$byClass))
  class_cancer$p53 <- p53
  class_cancer$gender <- gender
  class_cancer$feature_num <- dims_of_dat
  class_cancer$model_name <-  mod_name


  # get clinical data
  clin_data <- mod_results_list[[3]]

  
  return(list(class_cancer, clin_data))
}


get_class_results_test <- function(mod_results_list, dims_of_dat, mod_name, feat_name) {
  
  mat_onset_norm <-  mod_results_list[[1]]$table
  class_onset_norm <- as.data.frame(t(mod_results_list[[1]]$byClass))
  class_onset_norm$age_type <- 'cases_onset'
  class_onset_norm$feature_num <- dims_of_dat
  class_onset_norm$model_method <-  mod_name
  class_onset_norm$feat_name <-  feat_name
  
  
  # cases age
  mat_age_norm <-  mod_results_list[[2]]$table
  class_age_norm <- as.data.frame(t(mod_results_list[[2]]$byClass))
  class_age_norm$age_type <- 'cases_age'
  class_age_norm$feature_num <- dims_of_dat
  class_age_norm$model_method <-  mod_name
  class_age_norm$feat_name <-  feat_name
  
  # valid age
  mat_controls_age_norm <-mod_results_list[[3]]$table
  class_controls_age_norm <- as.data.frame(t(mod_results_list[[3]]$byClass))
  class_controls_age_norm$age_type <- 'controls_age'
  class_controls_age_norm$feature_num <- dims_of_dat
  class_controls_age_norm$model_method <- mod_name
  class_controls_age_norm$feat_name <-  feat_name
  
  # valid onset
  mat_valid_onset_norm <- mod_results_list[[4]]$table
  class_valid_onset_norm <- as.data.frame(t(mod_results_list[[4]]$byClass))
  class_valid_onset_norm$age_type <- 'valid_onset'
  class_valid_onset_norm$feature_num <- dims_of_dat
  class_valid_onset_norm$model_method <-  mod_name
  class_valid_onset_norm$feat_name <-  feat_name
  
  
  # valid age
  mat_valid_age_norm <-  mod_results_list[[5]]$table
  class_valid_age_norm <- as.data.frame(t(mod_results_list[[5]]$byClass))
  class_valid_age_norm$age_type <- 'valid_age'
  class_valid_age_norm$feature_num <- dims_of_dat
  class_valid_age_norm$model_method <- mod_name
  class_valid_age_norm$feat_name <-  feat_name
  
  
  # combine class
  class_results_norm <- rbind(class_onset_norm,
                              class_age_norm,
                              class_valid_age_norm,
                              class_valid_onset_norm, 
                              class_controls_age_norm)
  
  
  ##################################################################################################
  # RAND
  
  mat_onset_rand <-  mod_results_list[[6]]$table
  class_onset_rand <- as.data.frame(t(mod_results_list[[6]]$byClass))
  class_onset_rand$age_type <- 'cases_onset'
  class_onset_rand$feature_num <- dims_of_dat
  class_onset_rand$model_method <-  mod_name
  class_onset_rand$feat_name <-  feat_name
  
  
  # cases age
  mat_age_rand <-  mod_results_list[[7]]$table
  class_age_rand <- as.data.frame(t(mod_results_list[[7]]$byClass))
  class_age_rand$age_type <- 'cases_age_rand'
  class_age_rand$feature_num <- dims_of_dat
  class_age_rand$model_method <-  mod_name
  class_age_rand$feat_name <-  feat_name
  
  # valid age
  mat_controls_age_rand <-mod_results_list[[8]]$table
  class_controls_age_rand <- as.data.frame(t(mod_results_list[[8]]$byClass))
  class_controls_age_rand$age_type <- 'controls_age_rand'
  class_controls_age_rand$feature_num <- dims_of_dat
  class_controls_age_rand$model_method <- mod_name
  class_controls_age_rand$feat_name <-  feat_name
  
  # valid onset
  mat_valid_onset_rand <- mod_results_list[[9]]$table
  class_valid_onset_rand <- as.data.frame(t(mod_results_list[[9]]$byClass))
  class_valid_onset_rand$age_type <- 'valid_onset_rand'
  class_valid_onset_rand$feature_num <- dims_of_dat
  class_valid_onset_rand$model_method <-  mod_name
  class_valid_onset_rand$feat_name <-  feat_name
  
  
  # valid age
  mat_valid_age_rand <-  mod_results_list[[10]]$table
  class_valid_age_rand <- as.data.frame(t(mod_results_list[[10]]$byClass))
  class_valid_age_rand$age_type <- 'valid_age_rand'
  class_valid_age_rand$feature_num <- dims_of_dat
  class_valid_age_rand$model_method <- mod_name
  class_valid_age_rand$feat_name <-  feat_name
  
  
  # combine class
  class_results_rand <- rbind(class_onset_rand,
                              class_age_rand,
                              class_valid_age_rand,
                              class_valid_onset_rand, 
                              class_controls_age_rand)
  
  
  
  class_results <- rbind(class_results_norm,
                         class_results_rand)
  
  return(list(class_results, list(mat_onset_norm, mat_age_norm,
                                  mat_valid_onset_norm,  mat_valid_age_norm, 
                                  mat_controls_age_norm)))
}






# mod_results_list <- mod_result
get_class_results_rand <- function(mod_results_list, mod_name, dims_of_dat) {
  
  # NORM
  
  mat_onset_norm <-  Reduce('+',  lapply(lapply(mod_results_list, function(x) x[[1]]$table), function(x) rbind(x)))/5
  class_onset_norm <- as.data.frame(t(apply(do.call(rbind, lapply(mod_results_list, function(x) x[[1]]$byClass)), 2, function(x) mean(x, na.rm =T))))
  class_onset_norm$age_type <- 'cases_onset'
  class_onset_norm$feature_num <- dims_of_dat
  class_onset_norm$model_method <-  mod_name
  
  # cases age
  mat_age_norm <-  Reduce('+',  lapply(lapply(mod_results_list, function(x) x[[2]]$table), function(x) rbind(x)))/5
  class_age_norm <- as.data.frame(t(apply(do.call(rbind, lapply(mod_results_list, function(x) x[[2]]$byClass)), 2, function(x) mean(x, na.rm =T))))
  class_age_norm$age_type <- 'cases_age'
  class_age_norm$feature_num <- dims_of_dat
  class_age_norm$model_method <-  mod_name
  
  
  
  # valid age
  mat_controls_age_norm <-  Reduce('+',  lapply(lapply(mod_results_list, function(x) x[[3]]$table), function(x) rbind(x)))/5
  class_controls_age_norm <- as.data.frame(t(apply(do.call(rbind, lapply(mod_results_list, function(x) x[[3]]$byClass)), 2, function(x) mean(x, na.rm =T))))
  class_controls_age_norm$age_type <- 'controls_age'
  class_controls_age_norm$feature_num <- dims_of_dat
  class_controls_age_norm$model_method <- mod_name
  
  # valid onset
  mat_valid_onset_norm <-  Reduce('+',  lapply(lapply(mod_results_list, function(x) x[[4]]$table), function(x) rbind(x)))/5
  class_valid_onse_norm <- as.data.frame(t(apply(do.call(rbind, lapply(mod_results_list, function(x) x[[4]]$byClass)), 2, function(x) mean(x, na.rm =T))))
  class_valid_onset_norm$age_type <- 'valid_onset'
  class_valid_onset_norm$feature_num <- dims_of_dat
  class_valid_onset_norm$model_method <-  mod_name
  
  # valid age
  mat_valid_age_norm <-  Reduce('+',  lapply(lapply(mod_results_list, function(x) x[[5]]$table), function(x) rbind(x)))/5
  class_valid_age_norm <-as.data.frame(t(apply(do.call(rbind, lapply(mod_results_list, function(x) x[[5]]$byClass)), 2, function(x) mean(x, na.rm =T))))
  class_valid_age_norm$age_type <- 'valid_age'
  class_valid_age_norm$feature_num <- dims_of_dat
  class_valid_age_norm$model_method <- mod_name
  
  
  # combine class
  class_results <- rbind(class_onset_norm,
                         class_age_norm,
                         class_valid_age_norm,
                         class_valid_onse_norm, 
                         class_controls_age_norm)
  
  return(list(class_results, list(mat_onset, mat_age_norm,
                                  mat_valid_onset_norm,  mat_valid_age_norm, 
                                  mat_controls_age_norm,
                                  mat_age_rand,
                                  mat_valid_onset_rand,  mat_valid_age_rand, 
                                  mat_controls_age_rand)))
}


# train_data <- m_train_cases
# test_data <- m_test_cases
# remove duplicates
remove_dups <- function(train_data, test_data) {
  
  # give data indicator columns
  train_data$type <- 'train'
  test_data$type <- 'test'
  
  # combine
  full_data <- rbind(train_data,
                     test_data)
  
  # remove duplicates
  full_data <- full_data[!duplicated(full_data$ids),]
  
  train_data <- full_data[grepl('train', full_data$type),]
  test_data <- full_data[grepl('test', full_data$type),]
  
  train_data$type <- NULL
  test_data$type <- NULL
  
  
  return(list(train_data, test_data))
  
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
  
  if(bump == 'age') {
    
    dat$type <- ifelse(dat$age_sample_collection > , 'controls', 'cases')
    ##########
    # get indicator and put into design matrix with intercept 1
    #########
    indicator_vector <- as.factor(dat$type)
    designMatrix <- cbind(rep(1, nrow(dat)), indicator_vector)
    designMatrix <- as.matrix(designMatrix)
  }
  
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
  methyl_cg <- inner_join(dat, g_ranges, by = 'probe')
  
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
                           type = "Beta")
    
    bump_hunter_results[[i]] <- tab[[i]]$table
    bump_hunter_results[[i]]$run <- DELTA_BETA_THRESH[i]
  }
  
  bh_results <- do.call(rbind, bump_hunter_results)
  
  return(bh_results)
  
}




runEnetRoc <- function(training_dat,
                       test_dat,
                       age_cutoff,
                       bh_features,
                       gender) {
  
  

  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  if(gender) {
    
    intersected_feats <- append('M', intersected_feats)
    intersected_feats <- append('F', intersected_feats)
  }
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- ifelse(training_dat$age_diagnosis < age_cutoff, 1, 0)
  test_y <-  ifelse(test_dat$age_diagnosis < age_cutoff, 1, 0)
  
  # get clinical data
  test_clin <- test_dat[, 1:9]

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
  
  
  
  ##########
  # Predictions on test data
  ##########
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  
  # original should be fine, something wrong with caret package
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y))
  
  ###########################################################################################
  return(temp_dat)
  
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
                             base_change,
                             exon_intron,
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
  return(list(temp_dat, temp_dat_con))
  
}
