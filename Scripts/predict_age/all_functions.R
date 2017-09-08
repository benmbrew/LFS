# functions to be used in model_pipeline script

cleanIdMap <- function(data, valid) {
  
  data <- as.data.frame(data)
  
  # new colnames, lowercase
  colnames(data) <- tolower(colnames(data))
  
  # combine sentrix_id and sentrix_position 
  data$identifier <- paste(data$sentrix_id, data$sentrix_position, sep = '_')
  data$identifier <- as.factor(data$identifier)
  
  return(data)
  
}

##########
# function that Loops through list, preprocesses, and convert to beta, m, and cn values 
##########
preprocessMethod <- function(data, preprocess) {
  
  if (preprocess == 'raw') {
    
    Mset <- preprocessRaw(data)
  }
  
  if (preprocess == 'quan') {
    Mset   <-preprocessQuantile(data, fixOutliers = TRUE,
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
  
  # gset <- mapToGenome(ratioSet) 
  m <- getM(Mset)
  # beta <- getBeta(gset)
  return(m)
  
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
  full_data <- as.data.frame(cbind(data_frame[, 1:7], temp_data))
  
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
  full_data <- as.data.frame(cbind(data_frame[, 1:7], temp_data))
  
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
# run combat on cases, controls, and valid
##########

runCombat <- function(data,data_controls, data_valid , batch_var)
{
  
  if(batch_var == 'mon'){
    data$batch <- ifelse(grepl('9721365183', data$sentrix_id), 'tor', 'mon')
    
  } else {
    
    # make type column
    data$batch <- 'cases'
    data_controls$batch <- 'controls'
    data_valid$batch <- 'valid'
    
    # get common features
    intersected_feats <- Reduce(intersect, list(colnames(data), colnames(data_controls), colnames(data_valid)))
    
    # subset data by common featrues 
    data <- data[, intersected_feats]
    data_controls <- data_controls[,  intersected_feats]
    data_valid <- data_valid[, intersected_feats]
    
    #combine data
    combined_dat <- rbind(data, data_controls, data_valid)
    
    # rename as data
    data <- combined_dat
    
  }
  # get batch
  batch_indicator <- as.character(data$batch)
  batch_indicator <- as.factor(batch_indicator)
  
  gender <- data$gender
  sentrix_id <- data$sentrix_id
  sen_batch <- data$sen_batch
  ids <- data$ids
  sample_collection <- data$age_sample_collection
  diagnosis <- data$age_diagnosis
  p53 <- data$p53_germline
  cancer_diagnosis <- data$cancer_diagnosis_diagnoses
  # put model ids in rownames and remove columns
  mat_data <- data[, (8:ncol(data) - 1)]
  # get features
  features <- colnames(mat_data)
  mat_data <- t(mat_data)
  
  # get intercept
  modcombat <- model.matrix(~1, data = data)
  combat <- ComBat(dat = mat_data, batch = batch_indicator, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  # transpose and add back columns
  final_dat <- as.data.frame(t(combat))
  final_dat$gender <- gender
  final_dat$ids <- ids
  final_dat$age_sample_collection <- sample_collection
  final_dat$age_diagnosis <- diagnosis
  final_dat$p53_germline <- p53
  final_dat$sentrix_id <- sentrix_id
  final_dat$cancer_diagnosis_diagnoses <- cancer_diagnosis
  final_dat <- final_dat[, c('ids', 'p53_germline', 'age_diagnosis','cancer_diagnosis_diagnoses' ,
                             'age_sample_collection', 'gender', 'sentrix_id', features)]
  final_dat$sentrix_id.1 <- NULL
  rownames(final_dat) <- NULL
  
  return(final_dat)
  
}


##########
# Function that combines methylation matrices with id_map, to get ids for methylation
##########
# data_methyl <- betaCases

findIds <- function(data_methyl, id_map) {
  
  
  data_methyl <- as.data.frame(t(data_methyl))
  data_methyl$identifier <- rownames(data_methyl)
  data_methyl$identifier <- as.factor(data_methyl$identifier)
  # loop to combine identifiers, without merging large table
  data_methyl$ids <- NA
  data_methyl$sentrix_id <- NA
  
  for (i in data_methyl$identifier) {
    data_methyl$ids[data_methyl$identifier == i] <- id_map$sample_name[id_map$identifier == i]
    data_methyl$sentrix_id[data_methyl$identifier == i] <- id_map$sentrix_id[id_map$identifier == i]
    
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
  
  if (!control) {
    
    for (i in intersected_ids) {
      
      data$p53_germline[data$ids == i] <- clin$p53_germline[which(clin$ids == i)]
      data$age_diagnosis[data$ids == i] <- clin$age_diagnosis[which(clin$ids == i)]
      data$cancer_diagnosis_diagnoses[data$ids == i] <- clin$cancer_diagnosis_diagnoses[which(clin$ids == i)]
      data$age_sample_collection[data$ids == i] <- clin$age_sample_collection[which(clin$ids == i)]
      data$tm_donor_[data$ids == i] <- clin$tm_donor_[which(clin$ids == i)]
      data$gender[data$ids == i] <- clin$gender[which(clin$ids == i)]
      
      
      
      print(i)
    } 
    data <- data[!is.na(data$p53_germline),]
    data <- data[!duplicated(data$ids),]
    data <- data[!duplicated(data$tm_donor_),]
    # data <- data[!is.na(data$age_diagnosis),]
    data <- data[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses',
                     'age_sample_collection', 'gender','sentrix_id', features)]
    
  } else {
    
    for (i in intersected_ids) {
      
      data$p53_germline[data$ids == i] <- clin$p53_germline[which(clin$ids == i)]
      data$cancer_diagnosis_diagnoses[data$ids == i] <- clin$cancer_diagnosis_diagnoses[which(clin$ids == i)]
      data$age_sample_collection[data$ids == i] <- clin$age_sample_collection[which(clin$ids == i)]
      data$tm_donor_[data$ids == i] <- clin$tm_donor_[which(clin$ids == i)]
      data$gender[data$ids == i] <- clin$gender[which(clin$ids == i)]
      
      
      print(i)
    } 
    data <- data[!is.na(data$p53_germline),]
    data <- data[!duplicated(data$ids),]
    data <- data[!duplicated(data$tm_donor_),]
    data <- data[, c('ids', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses',
                     'age_sample_collection', 'gender', 'sentrix_id', features)]
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
  
  # get fwer results
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
# predict cancer
##########
predCancer <- function(training_dat, 
                       test_dat,
                       bh_features,
                       gender,
                       cutoff) 
{
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(training_dat))
  
  # # get y
  train_y <- as.factor(ifelse(grepl('Unaffected',training_dat$cancer_diagnosis_diagnoses), 'no', 'yes'))
  test_y <- as.factor(ifelse(grepl('Unaffected',test_dat$cancer_diagnosis_diagnoses), 'no', 'yes'))
  
  train_y <- factor(train_y, levels = c('yes', 'no'))
  test_y <- factor(test_y, levels = c('yes', 'no'))
  
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
  training_dat <- training_dat[, 1:100]
  test_dat <- test_dat[, 1:100]
  
  
  # controls_dat <- controls_dat[, intersected_feats]
  # valid_dat <- valid_dat[, intersected_feats]
  
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
  test.predictions <- factor(test.predictions, levels = c('yes', 'no'))
  
  test_stats <- confusionMatrix(test_y, test.predictions)
  importance <- coef(model)
  
  lambda_value <- elastic_net.cv_model$lambda.min
  
  alpha  <- best_alpha
  
  return(list(alpha, lambda_value, importance, test_stats, model, temp.non_zero_coeff))
  
  
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
# test model
# ##########

testModelThresh <- function(cases_dat,
                            controls_dat,
                            controls_dat_old,
                            controls_dat_full,
                            valid_dat,
                            features,
                            alpha)
{
  
  
  intersected_feats <- features
  # get y
  train_y <- as.numeric(cases_dat$age_diagnosis)
  test_y_controls <- as.numeric(controls_dat$age_sample_collection)
  test_y_controls_old <- as.numeric(controls_dat_old$age_sample_collection)
  test_y_controls_full <- as.numeric(controls_dat_full$age_sample_collection)
  test_y_valid <- as.numeric(valid_dat$age_diagnosis)
  
  # get bumphunter features
  cases_dat <- cases_dat[, intersected_feats]
  controls_dat <- controls_dat[, intersected_feats]
  controls_dat_old <- controls_dat_old[, intersected_feats]
  controls_dat_full <- controls_dat_full[, intersected_feats]
  valid_dat <- valid_dat[, intersected_feats]
  
  N_CV_REPEATS = 2
  nfolds = 3
  
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  elastic_net.cv_model = cv.glmnet(x = as.matrix(cases_dat)
                                   , y =  train_y
                                   , alpha = alpha
                                   , type.measure = 'deviance'
                                   , family = 'gaussian'
                                   , standardize=FALSE
                                   , nlambda = 100
                                   , nfolds = nfolds
                                   , parallel = TRUE
  )
  
  
  temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
  
  # # number of non zero coefficients at that lambda    
  temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
  print(temp.non_zero_coeff)  
  
  
  model  = glmnet(x = as.matrix(cases_dat)
                  , y =  train_y
                  ,alpha = alpha
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = 'gaussian')
  
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
  
  importance <- coef(model)
  
  lambda_value <- elastic_net.cv_model$lambda.min
  
  # # for each iteration, this should always be the same.
  controls_cor <- cor(test_y_controls, test.predictions_controls)
  
  controls_cor_old <- cor(test_y_controls_old, test.predictions_controls_old)
  
  
  controls_cor_full <- cor(test_y_controls_full, test.predictions_controls_full)
  
  valid_cor  <- cor(test_y_valid, test.predictions_valid)
  
  return(list(controls_cor, controls_cor_old, controls_cor_full, valid_cor))
  
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
# lasso rand
###########

runLassoRandResid <- function(training_dat,
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
  
  
  
  return(list(test.predictions, test_y, patient_age))
  
  
}

##########
# lasso rand
###########

runRidgeRandResid <- function(training_dat,
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
  
  
  
  return(list(test.predictions, test_y, patient_age))
  
  
}


##########
# enet diff
###########

runRfRandResid <- function(training_dat,
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
  
  
  
  return(list(test.predictions, test_y, patient_age))
  
  
}



##########
# enet diff
###########

runSvmRandResid <- function(training_dat,
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
  
  
  
  return(list(test.predictions, test_y, patient_age))
  
  
}


