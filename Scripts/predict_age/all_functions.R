
########
# load libraries
##########
library(IlluminaHumanMethylation450kmanifest)
library(tidyverse)
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
library(randomForest)
library(RPMM)
library(RColorBrewer)
##################

registerDoParallel(2)


get_diff_dups <- function(temp_data) {
  # get controls - HERE
  # get controls that are deduped from both sides
  
  data_cases <- temp_data[!grepl('Unaffected', temp_data$cancer_diagnosis_diagnoses),]
  data_controls <- temp_data[grepl('Unaffected', temp_data$cancer_diagnosis_diagnoses),]
  # remove duplicates from controls, from both sides 
  data_controls_from_last <- data_controls[!duplicated(data_controls$tm_donor_, fromLast = TRUE),]
  data_controls_from_first <- data_controls[!duplicated(data_controls$tm_donor_, fromLast = FALSE),]
  
  # combine data
  temp_data_from_first <- rbind(data_cases,
                                data_controls_from_first)
  
  # combine data
  temp_data_from_last <- rbind(data_cases,
                               data_controls_from_last)
  
  return(list(temp_data_from_first, temp_data_from_last))
  
}

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
   }
  # # map methyl set to genome (funnorm already does this)
  # Gset <- mapToGenome(Mset)
  # # get m values
  if(methyl_type == 'beta') {
    # dat_beta <- getBeta(Mset)
    dat <- wateRmelon::BMIQ(Mset)
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
    beta_cases <- beta_cases[!is.na(beta_cases$tm_donor),]
    beta_controls <- beta_controls[!is.na(beta_controls$tm_donor),]
    
    
    # remove duplicates
    beta_cases <- beta_cases[!duplicated(beta_cases$tm_donor),]
    beta_controls <- beta_controls[!duplicated(beta_controls$tm_donor),]
    
    
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
                                 'tm_donor',
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
                                       'tm_donor',
                                       cg_sites)]
    
    # combine data
    beta <- rbind(beta_cases,
                  beta_controls)
    
    # remove duplcated tm donor
    beta <- beta[!duplicated(beta$tm_donor),]
    
    return(beta)
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
run_combat <- function(temp_data) {
  temp_data <- as.data.frame(temp_data)
  # get tech variable variable back to categories 
  temp_data$tech <- ifelse(temp_data$tech == '450k', 'batch_1', 'batch_2')
  
  # get batch
  batch_indicator <- as.character(temp_data$tech)
  batch_indicator <- as.factor(batch_indicator)
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
# get family cancer status and ratio
##########

get_family_cancer <- function(cases_full, controls_full){
  #combine data subseted data
  temp_full <- bind_rows(cases_full[, c('tm_donor', 'cancer_diagnosis_diagnoses', 'family_name')],
                         controls_full[, c('tm_donor', 'cancer_diagnosis_diagnoses', 'family_name')])
  # create cancer indicator
  temp_full$cancer_fac <- ifelse(grepl('Unaffected', temp_full$cancer_diagnosis_diagnoses), 'no_cancer', 'cancer')
  temp_full$cancer_diagnosis_diagnoses <- NULL
  
  temp <- temp_full %>%
    group_by(family_name) %>%
    tally() %>% 
    left_join(temp_full)
  
  temp$fam_num_cancer <- NA
  temp$fam_cancer_ratio <- NA
  for(fam_name in unique(temp$family_name)){
    sub_fam <- temp[temp$family_name == fam_name,]
    if(nrow(sub_fam) > 1) {
      num_cancer <- length(which(sub_fam$cancer_fac == 'cancer'))
      num_no <- length(which(sub_fam$cancer_fac == 'no_cancer'))
      tot_fam_num <- nrow(sub_fam)
      stopifnot((num_cancer + num_no) == tot_fam_num)
      
      # get number of family members with cancer 
      sub_fam$fam_num_cancer <- ifelse(sub_fam$cancer_fac == 'cancer',
                                       num_cancer - 1, num_cancer)
      
      # condition for if the denominator is zero - that is no non cancers in family - just put number of family members with cancer 
      # this applies to the situation where there might be a non cancer present, 
      if(num_no == 0) {
        sub_fam$fam_cancer_ratio <- num_cancer - 1
        
      } else if(num_no == 1){
        # get ratio of family members that have cancer to those that don't have cancer
        sub_fam$fam_cancer_ratio <- ifelse(sub_fam$cancer_fac == 'cancer',
                                           (num_cancer - 1)/num_no, num_cancer)
      } else {
        
        # conidtion if num_cancer is zero
        if(num_cancer == 0) {
          # get ratio of family members that have cancer to those that don't have cancer
          sub_fam$fam_cancer_ratio <- 0
        } else {
          # get ratio of family members that have cancer to those that don't have cancer
          sub_fam$fam_cancer_ratio <- ifelse(sub_fam$cancer_fac == 'cancer',
                                             (num_cancer - 1)/num_no, num_cancer/(num_no - 1))
        }
        
      }
      
    } else {
      sub_fam$fam_num_cancer <- 0
      sub_fam$fam_cancer_ratio <- 0
    }
    temp[temp$family_name == fam_name,] <- sub_fam
    
  }
  # remove columns not needed in join
  temp <- temp[, c('tm_donor', 'family_name', 'fam_num_cancer', 'fam_cancer_ratio')]
  
  # join temp back with cases and controls 
  cases_full <- inner_join(temp, cases_full)
  controls_full <- inner_join(temp, controls_full)
  
  return(list(cases_full, controls_full))
  
}


get_family_cancer_old <- function(cases, controls, valid){
  #combine data subseted data
  temp_full <- bind_rows(cases[, c('tm_donor', 'cancer_diagnosis_diagnoses', 'family_name')],
                         controls[, c('tm_donor', 'cancer_diagnosis_diagnoses', 'family_name')],
                         valid[,c('tm_donor', 'cancer_diagnosis_diagnoses', 'family_name') ])
  # create cancer indicator
  temp_full$cancer_fac <- ifelse(grepl('Unaffected', temp_full$cancer_diagnosis_diagnoses), 'no_cancer', 'cancer')
  temp_full$cancer_diagnosis_diagnoses <- NULL
  
  temp <- temp_full %>%
    group_by(family_name) %>%
    tally() %>% 
    left_join(temp_full)
  
  temp$fam_num_cancer <- NA
  temp$fam_cancer_ratio <- NA
  for(fam_name in unique(temp$family_name)){
    sub_fam <- temp[temp$family_name == fam_name,]
    if(nrow(sub_fam) > 1) {
      num_cancer <- length(which(sub_fam$cancer_fac == 'cancer'))
      num_no <- length(which(sub_fam$cancer_fac == 'no_cancer'))
      tot_fam_num <- nrow(sub_fam)
      stopifnot((num_cancer + num_no) == tot_fam_num)
      
      # get number of family members with cancer 
      sub_fam$fam_num_cancer <- ifelse(sub_fam$cancer_fac == 'cancer',
                                       num_cancer - 1, num_cancer)
      
      # condition for if the denominator is zero - that is no non cancers in family - just put number of family members with cancer 
      # this applies to the situation where there might be a non cancer present, 
      if(num_no == 0) {
        sub_fam$fam_cancer_ratio <- num_cancer - 1
        
      } else if(num_no == 1){
        # get ratio of family members that have cancer to those that don't have cancer
        sub_fam$fam_cancer_ratio <- ifelse(sub_fam$cancer_fac == 'cancer',
                                           (num_cancer - 1)/num_no, num_cancer)
      } else {
        
        # conidtion if num_cancer is zero
        if(num_cancer == 0) {
          # get ratio of family members that have cancer to those that don't have cancer
          sub_fam$fam_cancer_ratio <- 0
        } else {
          # get ratio of family members that have cancer to those that don't have cancer
          sub_fam$fam_cancer_ratio <- ifelse(sub_fam$cancer_fac == 'cancer',
                                             (num_cancer - 1)/num_no, num_cancer/(num_no - 1))
        }
        
      }
      
    } else {
      sub_fam$fam_num_cancer <- 0
      sub_fam$fam_cancer_ratio <- 0
    }
    temp[temp$family_name == fam_name,] <- sub_fam
    
  }
  # remove columns not needed in join
  temp <- temp[, c('tm_donor', 'family_name', 'fam_num_cancer', 'fam_cancer_ratio')]
  
  # join temp back with cases and controls 
  cases <- inner_join(temp, cases)
  controls <- inner_join(temp, controls)
  valid <- inner_join(temp, valid)
  
  return(list(cases, controls, valid))
  
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
  data <- data[data$p53_germline == 'MUT',]
  return(data)
}

testKS <- function(x, y)
{
  y <- y[!is.na(y)]
  x <- x[!is.na(x)]
  
  # Do x and y come from the same distribution?
  ks.test(jitter(x), jitter(y), alternative = 'two.sided')
  
}


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
  transform_850 <- t(temp)
  
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
                                       cancer_status = full_data$cancer_status,
                                       tech = full_data$tech,
                                       transform_850))
  
  # make numeric
  transform_850[, 12:ncol(transform_850)] <- apply(transform_850[, 12:ncol(transform_850)], 
                                                   2, 
                                                   function(x) as.numeric(x))
  
  
  return(transform_850)
  
}



##########
# get pca function
##########
# pca_data = all_cases_beta_combat
# age_cutoff = 72
# column_name = 'tech'
# show_variance = FALSE
# pc_x = 1
# pc_y = 2
# main_title = 'PC:cases bet'

get_pca <- function(pca_data, 
                    column_name,
                    age_cutoff,
                    show_variance,
                    pc_x,
                    pc_y,
                    main_title) {
  # # make function to get pca variables 
  # pca_data$tech <- ifelse(pca_data$a == 1, '405k', '850k')
  # 
  # # get cancer fac
  # pca_data$cancer_fac <- ifelse(!grepl('Unaffected', pca_data$cancer_diagnosis_diagnoses), 'cancer', 'no_cancer')
  # 
  # # get an age category and fac
  # pca_data$age_cat <- ifelse(pca_data$age_sample_collection >= 0 & pca_data$age_sample_collection <= 50, '0_50', 
  #                            ifelse(pca_data$age_sample_collection > 50 & pca_data$age_sample_collection <= 100, '51_100',
  #                                   ifelse(pca_data$age_sample_collection > 100 & pca_data$age_sample_collection <= 150, '101_150',
  #                                          ifelse(pca_data$age_sample_collection > 150 & pca_data$age_sample_collection <= 200, '151_200', 
  #                                                 ifelse(pca_data$age_sample_collection > 200 & pca_data$age_sample_collection <= 250,'201_250',
  #                                                        ifelse(pca_data$age_sample_collection > 250 & pca_data$age_sample_collection <= 300, '251_300', '300+'))))))
  # 
  # # relevel this data
  # pca_data$age_cat <- factor(pca_data$age_cat, levels = c('0_50', '51_100', '101_150', '151_200', '201_250', '251_300', '300+'))
  # 
  # 
  # pca_data$age_fac <- ifelse(pca_data$age_sample_collection >= age_cutoff, 'adult', 'not_adult')
  # 
  
  pca_data <- as.data.frame(pca_data)
  pca_data[, column_name] <- as.factor(pca_data[, column_name])
  
  # get other clinical data
  column_names <- c('tm_donor','ids' ,'cancer_diagnosis_diagnoses', 
                    'age_diagnosis', 'age_sample_collection', 'gender')
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
  pca <- prcomp(pca_data[,8:data_length])
  
  if(show_variance){
    # plot lambda
    return(plot(pca, xlim = c(0,10), type = 'l', main = 'PCs and variance'))
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
    cols <- colorRampPalette(brewer.pal(n = 9, 'Spectral'))(length(unique(pca_results$column_name)))
    
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
      theme_minimal(base_size = 18, base_family = 'ubuntu')
    
    # plot <- ggplotly(plot)
    
    return(plot)
  }
  
}

# add age as dummy variable 
get_age_dummy <- function(temp_dat, the_age){
  
  temp_dat$age_dum_ <- ifelse(temp_dat$age_sample_collection > the_age, 1, 0)
  temp_dat$age_dum_young <- ifelse(temp_dat$age_sample_collection <= the_age, 1, 0)
  
  return(temp_dat)
  
}

get_age_cat <- function(temp_dat){
  temp_dat$age_var <- ntile(temp_dat$age_sample_collection, 5)
  return(temp_dat)
}


##########
# remove outliers  4257 cases, 3391, 3392 controls
##########

removeOutlier <- function(data, 
                          cases, 
                          controls, 
                          val) {
  
  if (cases) {
    data <- data[data$tm_donor != '3955',]
    
  }
  #controls outlier
  if (controls) {
    
    data <- data[data$tm_donor != '3847',]
    data <- data[data$ids != '3484',]
  }
  
  if(val){
    data <- data[data$ids != '3540',]
    
  }
  
  
  return(data)
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



# dat = train_cases, 
# wild_type = wt_data,
# boot_num = 5, 
# thresh = beta_thresh,
# g_ranges = g_ranges
bump_hunter_lfs <- function(dat,
                        wild_type,
                        boot_num,
                        thresh,
                        g_ranges) {
  

  
    
  # get mutual cgs
  dat_names <- colnames(dat)[grepl('^cg', colnames(dat))]
  wt_names <- colnames(wild_type)[grepl('^cg', colnames(wild_type))]
  wt_intersect <- intersect(wt_names, dat_names)
  
  # stor clinical data
  clin_dat_names <- colnames(dat)[!grepl('^cg', colnames(dat))]
  clin_wt_names <- colnames(wild_type)[!grepl('^cg', colnames(wild_type))]
  
  # get data
  dat <- dat[, c(clin_dat_names, wt_intersect)]
  wild_type <- wild_type[, c(clin_wt_names, wt_intersect)]
  
  # drop columns so they can match
  dat$fam_cancer_ratio <- dat$fam_num_cancer <- NULL
  
  # combine data
  dat <- rbind(dat,
               wild_type)
  
  # remove NAs 
  dat <- dat[!is.na(dat$p53_germline),]
  
  # get indicator (bump) vector
  dat$type <- dat$p53_germline
  indicator_vector <- as.factor(dat$type)
  designMatrix <- cbind(rep(1, nrow(dat)), indicator_vector)
  designMatrix <- as.matrix(designMatrix)
  
  
  # get probe data 
  dat <- dat[, grepl('^cg', colnames(dat))]
  
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



# dat = train_cases
# wild_type = wt_data
# boot_num = 5
# thresh = beta_thresh
# g_ranges = g_ranges

bump_hunter_surv <- function(dat,
                            wild_type,
                            boot_num,
                            thresh,
                            g_ranges) {
  
  # get mutual cgs
  dat_names <- colnames(dat)[grepl('^cg', colnames(dat))]
  wt_names <- colnames(wild_type)[grepl('^cg', colnames(wild_type))]
  wt_intersect <- intersect(wt_names, dat_names)
  
  # stor clinical data
  clin_dat_names <- colnames(dat)[!grepl('^cg', colnames(dat))]
  clin_wt_names <- colnames(wild_type)[!grepl('^cg', colnames(wild_type))]
  
  # get data
  dat <- dat[, c(clin_dat_names, wt_intersect)]
  wild_type <- wild_type[, c(clin_wt_names, wt_intersect)]
  
  # drop columns so they can match
  dat$fam_cancer_ratio <- dat$fam_num_cancer <- NULL
  
  # combine data
  dat <- rbind(dat,
               wild_type)
  
  # remove NAs 
  dat <- dat[!is.na(dat$p53_germline),]
  
  # get indicator (bump) vector
  dat$type <- dat$p53_germline
  indicator_vector <- as.factor(dat$type)
  designMatrix <- cbind(rep(1, nrow(dat)), indicator_vector)
  designMatrix <- as.matrix(designMatrix)
  
  
  # get probe data 
  dat <- dat[, grepl('^cg', colnames(dat))]
  
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

# training_dat = train_cases
# controls_dat = controls_full
# valid = valid_full
# test_dat = test_cases
# age_cutoff = 72
# gender = T
# tech = T
# fam_num = T
# fam_ratio = T
# bh_features = remaining_features

run_rf <- function(training_dat,
                   controls_dat,
                   test_dat,
                   age_cutoff,
                   age_dum,
                   gender, 
                   tech,
                   fam_num,
                   fam_ratio,
                   bh_features) {
  
  
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
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
  
  
  if (age_dum){
    intersected_feats <- c('first', 'second', 'third',intersected_feats)
  }
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- ifelse(training_dat$age_diagnosis < age_cutoff, 'yes', 'no')
  test_y <-  ifelse(test_dat$age_diagnosis < age_cutoff, 'yes', 'no')
  controls_y <-  ifelse(controls_dat$age_sample_collection < age_cutoff, 'yes', 'no')
  
  # get clinical data
  test_clin <- test_dat[, !grepl('^cg', colnames(test_dat))]
  controls_clin <- controls_dat[, !grepl('^cg', colnames(controls_dat))]
  
  # get model data
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  
  # determines how you train the model.
  NFOLDS <- 5
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = min(10, NFOLDS),
    classProbs = TRUE,
    repeats = 5,
    allowParallel = TRUE,
    summaryFunction = twoClassSummary
    
  )
  
  # mtry: Number of variables randomly sampled as candidates at each split.
  # ntree: Number of trees to grow.
  
  mtry <- sqrt(ncol(training_dat[,colnames(training_dat)]))
  tunegrid <- expand.grid(.mtry=mtry)
  
  data_size <- ncol(training_dat)
  model <- train(x = training_dat
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
  
  
  # Predictions on controls data
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions_con <- predict(model, 
                                       data.matrix(controls_dat),
                                       type = 'prob')

  # combine predictions and real labels 
  temp_dat_con <- as.data.frame(cbind(controls_age_pred = test.predictions_con, controls_age_label = controls_y, controls_clin))
  
  
  return(list(temp_dat, temp_dat_con,  model, importance, data_size))
  
  
  
}


run_enet_all_test <- function(training_dat,
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
                                                , standardize = FALSE 
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
                                     , standardize=FALSE
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
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  lambda_value_2 <- model$lambda[model$lambda == lambda_value]
  
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

# training_dat = train_cases
# test_dat = test_cases
# controls_dat = con_transform
# valid_dat = valid_transform
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
  test_results$pred_class <- as.factor(ifelse(test_results$positive > .5, 'positive', 'negative'))
  
  # relevel both factore
  test_results$pred_class <- factor(test_results$pred_class, c('positive', 'negative'))
  test_results$real <- factor(test_results$real, c('positive', 'negative'))
  
  test_results$accuracy <- caret::confusionMatrix(table(test_results$pred_class, test_results$real))$overall[1]
  test_results$tot_probes <- ncol(training_dat)
  
  
  return(list(test_results, importance))
  

}



# get age dummy category 
get_age_cat_dummy <- function(temp_dat) {
  
  temp_dat$temp_age_var <- ntile(temp_dat$age_sample_collection, 5)
  temp_dat <- cbind(as.data.frame(class.ind(temp_dat$temp_age_var)), temp_dat)
  colnames(temp_dat)[1:5] <- c('first', 'second', 'third', 'fourth', 'fifth')
  temp_dat$temp_age_var <- NULL
  return(temp_dat)
}
#
# training_dat = train_cases
# controls_dat = controls_full
# test_dat = test_cases
# age_cutoff = 72
# gender = FALSE
# tech = TRUE
# bh_features
#

run_enet_450_850 <- function(training_dat,
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
    intersected_feats <- c('Female', 'Male', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('batch_1', 'batch_2', intersected_feats)
  }
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- as.factor(ifelse(training_dat$age_diagnosis < age_cutoff, 'positive', 'negative'))
  test_y <-  as.factor(ifelse(test_dat$age_diagnosis < age_cutoff, 'positive', 'negative'))
  controls_y <-  as.factor(ifelse(controls_dat$age_sample_collection < age_cutoff, 'positive', 'negative'))
  
  # get clinical data
  test_clin <- test_dat[, !grepl('^cg', colnames(test_dat))]
  controls_clin <- controls_dat[, !grepl('^cg', colnames(controls_dat))]
  
  # get model data
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
  temp.non_zero_coeff_min = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff_min < 1) { 
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
    temp.1se_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.1se) 
    
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff_min = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    temp.non_zero_coeff_1se = elastic_net.cv_model$nzero[temp.1se_lambda_index] 
    
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
  
  
  return(list(temp_dat, temp_dat_con, model, temp.non_zero_coeff_min,
              temp.non_zero_coeff_1se, best_alpha))
  
}

# run_enet_all <- function(training_dat,
#                         controls_dat,
#                         valid_dat,
#                         age_cutoff,
#                         gender,
#                         tech,
#                         bh_features) {
#   
#   
#   # get intersection of bh features and real data
#   bh_features <- as.character(unlist(bh_features))
#   
#   intersected_feats <- intersect(bh_features, colnames(training_dat))
#   
#   if(gender) {
#     intersected_feats <- c('Female', 'Male', intersected_feats)
#   }
#   if (tech) {
#     intersected_feats <- c('batch_1', 'batch_2', intersected_feats)
#   }
#   
#   # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
#   # # get y
#   train_y <- as.factor(ifelse(training_dat$age_diagnosis < age_cutoff, 'positive', 'negative'))
#   valid_y <-  as.factor(ifelse(valid_dat$age_diagnosis < age_cutoff, 'positive', 'negative'))
#   controls_y <-  as.factor(ifelse(controls_dat$age_sample_collection < age_cutoff, 'positive', 'negative'))
#   
#   # get clinical data
#   valid_clin <- test_dat[, !grepl('^cg', colnames(test_dat))]
#   controls_clin <- controls_dat[, !grepl('^cg', colnames(controls_dat))]
#   
#   # get model data
#   training_dat <- training_dat[, intersected_feats]
#   controls_dat <- controls_dat[, intersected_feats]
#   valid_dat <- valid_dat[, intersected_feats]
#   
# 
#   # get bumphunter features
#   training_dat <- training_dat[, intersected_feats]
#   controls_dat <- controls_dat[, intersected_feats]
#   valid_dat <- valid_dat[, intersected_feats]
#   
#   
#   # start elastic net tuning
#   N_CV_REPEATS = 2
#   nfolds = 3
#   
#   ###### ENET
#   # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
#   # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
#   elastic_net.cv_error = vector()
#   elastic_net.cv_model = list()
#   elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
#   
#   # set parameters for training model
#   type_family <- 'binomial'
#   type_measure <- 'auc'
#   
#   # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
#   # or if you have a high number fo N_CV_REPEATS
#   temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
#     for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
#     {      
#       elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
#                                                 , y =  train_y
#                                                 , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
#                                                 , type.measure = type_measure
#                                                 , family = type_family
#                                                 , standardize = FALSE 
#                                                 , nfolds = nfolds 
#                                                 , nlambda = 10
#                                                 , parallel = TRUE
#       )
#       elastic_net.cv_error[alpha] = min(elastic_net.cv_model[[alpha]]$cvm)
#     }
#     elastic_net.cv_error # stores 9 errors    
#   }
#   
#   if (N_CV_REPEATS == 1) {
#     temp.cv_error_mean = temp.cv_error_matrix
#   } else {
#     temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean) # take the mean of the 5 iterations  
#     # as your value for alpha
#   }
#   
#   # stop if you did not recover error for any models 
#   stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
#   
#   # get index of best alpha (lowest error) - alpha is values 0.1-0.9
#   temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
#   # print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
#   best_alpha <- elastic_net.ALPHA[temp.best_alpha_index]
#   temp.non_zero_coeff = 0
#   temp.loop_count = 0
#   # loop runs initially because temp.non_zero coefficient <3 and then stops 
#   # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
#   # it they are never greater than 1, then the model does not converge. 
#   while (temp.non_zero_coeff < 1) { 
#     elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
#                                      , y =  train_y
#                                      , alpha = elastic_net.ALPHA[temp.best_alpha_index]
#                                      , type.measure = type_measure
#                                      , family = type_family
#                                      , standardize=FALSE
#                                      , nlambda = 100
#                                      , nfolds = nfolds
#                                      , parallel = TRUE
#     )
#     
#     # get optimal lambda - the tuning parameter for ridge and lasso
#     # THIS IS IMPORTANT BECAUSE WHEN YOU TRAIN THE MODEL ON 100 SEPERATE VALUES OF LAMBDA
#     # AND WHEN YOU TEST THE MODEL IT WILL RETURN PREDCITION FOR ALL THOSE VALUES (1-100). YOU NEED TO 
#     # GRAB THE PREDICTION WITH SAME LAMBDA THAT YOU TRAINED ON. ITS ALL IN THE CODE, BUT JUST WANTED TO 
#     # GIVE YOU REASONS
#     temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
#     
#     # # number of non zero coefficients at that lambda    
#     temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
#     temp.loop_count = temp.loop_count + 1
#     
#     # set seed for next loop iteration
#     as.numeric(Sys.time())-> t 
#     set.seed((t - floor(t)) * 1e8 -> seed) 
#     if (temp.loop_count > 10) {
#       print("diverged")
#       temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
#       break
#     }
#   }# while loop ends 
#   # print(temp.non_zero_coeff)  
#   
#   model  = glmnet(x = as.matrix(training_dat)
#                   , y =  train_y
#                   ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
#                   ,standardize=FALSE
#                   ,nlambda = 100
#                   ,family = type_family)
#   
#   # Predictions on test data
# 
#   # Predictions on controls data
#   
#   # This returns 100 prediction with 1-100 lambdas
#   temp_test.predictions_con <- predict(model, 
#                                        data.matrix(controls_dat),
#                                        type = 'response')
#   
#   # get predictions with corresponding lambda.
#   test.predictions_con <- temp_test.predictions_con[, temp.min_lambda_index]
#   
#   # combine predictions and real labels 
#   temp_dat_con <- as.data.frame(cbind(controls_age_pred = test.predictions_con, controls_age_label = controls_y, controls_clin))
#   
#   # Predictions on controls data
#   
#   # This returns 100 prediction with 1-100 lambdas
#   temp_test.predictions_valid <- predict(model, 
#                                          data.matrix(valid_dat),
#                                          type = 'response')
#   
#   # get predictions with corresponding lambda.
#   test.predictions_valid <- temp_test.predictions_valid[, temp.min_lambda_index]
#   
#   # combine predictions and real labels 
#   temp_dat_valid <- as.data.frame(cbind(valid_age_pred = test.predictions_valid, valid_age_label = valid_y, valid_clin))
#   
#   ###########################################################################################
#   return(list(temp_dat, temp_dat_con, temp_dat_valid, model, elastic_net.cv_model$lambda.min, best_alpha))
#   
# }
# 

# cases_dat = cases_full
# controls_dat = controls_full
# age_cutoff = 72
# gender = gender
# tech = tech
# bh_features = remaining_features

run_enet_test <- function(cases_dat,
                          controls_dat,
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
  cases_y <- ifelse(cases_dat$age_diagnosis < age_cutoff, 'positive', 'negative')
  # controls
  controls_y <-  ifelse(controls_dat$age_sample_collection < 'positive', 'negative')

  
  # get clinical data
  cg_start <- which(grepl('^cg', colnames(cases_dat)))[1]
  cases_clin <- cases_dat[, 1:(cg_start - 1)]
  controls_clin <- controls_dat[, 1:(cg_start - 1)]

  
  # get bumphunter features
  cases_dat <- cases_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]

  # start elastic net tuning
  N_CV_REPEATS = 5
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
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(cases_dat)
                                                , y =  cases_y
                                                , alpha = elastic_net.ALPHA[alpha] 
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
  temp.non_zero_coeff_min = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff_min < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(cases_dat)
                                     , y =  cases_y
                                     , alpha =0.1 #elastic_net.ALPHA[temp.best_alpha_index]
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
    temp.1se_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.1se) 
    
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff_min = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    temp.non_zero_coeff_1se = elastic_net.cv_model$nzero[temp.1se_lambda_index] 
    
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
  return(list(model, temp_dat_con, temp.non_zero_coeff_min, temp.non_zero_coeff_1se, best_alpha))
  
}
# training_dat = train_cases[, 1:100]
# test_dat = test_cases[, 1:100]
# controls_dat = controls_full[, 1:100]
# valid_dat = valid_full[, 1:100]
# age_cutoff = age_cutoff
# gender = gender
# tech = tech
# offset = use_offset
# bh_features = bh_features


# 
# training_dat = all_train
# test_dat = all_test
# age_cutoff = age_cutoff
# gender = gender
# tech = tech
# alpha_value = alpha_num
# bh_features = bh_features
run_enet_surv <- function(training_dat,
                         test_dat,
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
  
  
  # get survival time as days to onset and days to sample collection in one column for training dat
  time_to_event <- training_dat$age_diagnosis
  missing_ind <- is.na(time_to_event)
  time_to_collection <- training_dat$age_sample_collection
  
  time_to_event[missing_ind] <- time_to_collection[missing_ind]
  
  # get cancer status 
  cancer_status <- ifelse(training_dat$cancer_diagnosis_diagnoses != 'Unaffected', TRUE, FALSE)
  
  # for test dat
  time_to_event_test <- test_dat$age_diagnosis
  missing_ind_test <- is.na(time_to_event_test)
  time_to_collection_test <- test_dat$age_sample_collection
  
  time_to_event_test[missing_ind_test] <- time_to_collection_test[missing_ind_test]
  
  # get cancer status 
  cancer_status_test <- ifelse(test_dat$cancer_diagnosis_diagnoses != 'Unaffected', TRUE, FALSE)
  
  # get clinical data
  test_clin <- test_dat[, !grepl('^cg', colnames(test_dat))]

  # get model data
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]

  # get survival object
  surv_outcome <- Surv(time_to_event,cancer_status, type = 'right')
  # get survival object
  surv_outcome_test <- Surv(time_to_event_test,cancer_status_test, type = 'right')
  
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
  
  # combine predictions and real labels 
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions,  test_clin))
  
  ###########################################################################################
  return(temp_dat)

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



### for one shot train and test on family overlaps
run_enet_family <- function(training_dat,
                            test_dat,
                            age_cutoff,
                            age_dum,
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
  
  
  if (age_dum){
    intersected_feats <- c('age_dum_young', 'age_dum_old' ,intersected_feats)
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
                           age_dum,
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
  
  
  if (age_dum){
    intersected_feats <- c('age_dum_young', 'age_dum_old' ,intersected_feats)
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

