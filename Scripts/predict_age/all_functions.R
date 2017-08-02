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
    ratioSet <- ratioConvert(Mset, what = 'both', keepCN = T)
  }
  
  if (preprocess == 'quan') {
    ratioSet   <- preprocessQuantile(data, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, stratified = TRUE,
                                     mergeManifest = FALSE, sex = NULL)
  }
  
  if (preprocess == 'swan') {
    Mset  <-preprocessSWAN(data)
    ratioSet <- ratioConvert(Mset, what = 'both', keepCN = T)
  } 
  
  if (preprocess == 'funnorm') {
    ratioSet <-preprocessFunnorm(data)
  }
  
  gset <- mapToGenome(ratioSet) 
  beta <- getBeta(gset)
  return(beta)
  
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
  
  if(scale){
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
    data_frame[, probe_start:ncol(data_frame)][, sapply(data_frame[, probe_start:ncol(data_frame),], 
                                                        function(x) all(!is.na(x)))]
  
  # combine probes with clin
  full_data <- as.data.frame(cbind(data_frame[, 1:7], temp_data))
  
  # check that it worked
  stopifnot(all(!is.na(full_data[, probe_start:ncol(full_data)])))
  
  return(full_data)
  
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
# get model data 
##########
getModData <- function(data) 
{
  # subset data by not na in age of diagnosis and mut
  data <- data[!is.na(data$age_diagnosis),]
  data <- data[data$p53_germline == 'Mut',]
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
    legend('bottomright',  
           legend = unique(pca_data$batch), 
           col=1:length(unique(pca_data$batch)), 
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
  DELTA_BETA_THRESH = .30 # DNAm difference threshold
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
  DELTA_BETA_THRESH = .20 # DNAm difference threshold
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
getFolds <- function(model_dat, seed_number, k_num){
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

# # #

runEnet <- function(training_dat, 
                    test_dat,
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
  patient_age <- as.numeric(test_dat$age_sample_collection)
  missing_ind <- !is.na(patient_age)
  patient_age <-patient_age[missing_ind]
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  

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
                                   type = 'response')
  
  
  
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  importance <- coef(model)
  
  lambda_value <- elastic_net.cv_model$lambda.min
  
  cases_cor  <- cor(test_y, test.predictions)
  
  age_cor <- cor(test.predictions[missing_ind], patient_age)
  
  alpha  <- best_alpha
  
  return(list(alpha, lambda_value, importance, cases_cor, age_cor, temp.non_zero_coeff))


}
# training_dat = cases[train_index,]
# test_dat = cases[test_index,]
# bh_features = bh_feat_sig
# gender = T
# cutoff = 48

runEnetFac <- function(training_dat, 
                       test_dat,
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
  
  if(gender) {
    
    intersected_feats <- append('gender', intersected_feats)
    training_dat$gender <- as.numeric(as.factor(training_dat$gender))
    test_dat$gender <- as.numeric(as.factor(test_dat$gender))
    # controls_dat$gender <- as.numeric(as.factor(controls_dat$gender))
    # valid_dat$gender <- as.numeric(as.factor(valid_dat$gender))
    
  }
  
  
  # # get y
  train_y <- factor(training_dat$cancer_diagnosis_diagnoses, levels = c('yes', 'no'))
  test_y <- factor(test_dat$cancer_diagnosis_diagnoses, levels = c('yes', 'no'))
  
  
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
# get results
##########

getResults <- function(result_list, result_list_resid)
{
  
  
  results_final <- list(result_list[[1]],
                        result_list[[2]],
                        result_list[[3]],
                        result_list[[4]],
                        result_list[[5]],
                        result_list[[6]])
  
  results_final_resid <- list(result_list_resid[[1]],
                              result_list_resid[[2]],
                              result_list_resid[[3]],
                              result_list_resid[[4]],
                              result_list_resid[[5]],
                              result_list_resid[[6]])
  
  return(list(results_final, results_final_resid))
  
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


