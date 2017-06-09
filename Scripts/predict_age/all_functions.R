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
                   pca2) 
{
  
  
  pca_data[, column_name] <- as.factor(pca_data[, column_name])
  
  # # subet data so only p53 mut
  # pca_data <- pca_data[pca_data$p53_germline == 'Mut',]
  # 
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
  temp <- pca$x[,1]
  
  temp
  
  # plot data 4257,  94, 93
  #fill in factors with colors 
  col_vec <- c('black','red' , 'green', 'bisque', 'bisque1', 'bisque2', 'lightblue', 
               'blueviolet', 'brown', 'cyan', 'coral',
               'grey', 'orange', 'yellow', 'darkblue','darkred', 
               'darkgreen', 'darkorchid', 'gold', 'darkorange', 'deeppink',
               'greenyellow', 'purple')
  
  colors <- col_vec[pca_data[, column_name]]
  
  
  plot <- plot(pca$x[, pca1], 
               pca$x[, pca2],
               xlab = 'pca',
               ylab = 'pca',
               cex = 1,
               main = name,
               pch = 16,
               col = adjustcolor(colors, alpha.f = 0.5)
  )
  abline(v = c(0,0),
         h = c(0,0))
  
  return(plot)
}
