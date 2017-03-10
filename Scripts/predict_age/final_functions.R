
getSubset <- function(model_data) {
  # homogenize model_data to model_data_resdi
  feat_names <- colnames(model_data)[10:ncol(model_data)]
  
  model_data <- model_data[,c('age_diagnosis', 'age_sample_collection', 'gender', 'type', feat_names)] 
  
  return(model_data)
}


getResidual <- function(model_data) 
{
  
  # subset by mut, and complete cases for age diagnosis and age sample collection
  model_data <- model_data[complete.cases(model_data),]

  feature_names <- colnames(model_data)[10:ncol(model_data)]
  
  resid <- list()
  
  for (i in 10:ncol(model_data)) {
    
    temp_response <- model_data[, i]
    temp_var <- model_data$age_sample_collection
    
    resid[[i]] <- lm(temp_response ~ temp_var)$residuals
    
    print(i)
    
  }
  
    resid_data <- as.data.frame(do.call('cbind', resid))
    model_data <- cbind(model_data$age_diagnosis, model_data$age_sample_collection, model_data$gender, 
                        model_data$type, resid_data)
    colnames(model_data) <- c('age_diagnosis', 'age_sample_collection', 'gender', 'type', feature_names)

  return(model_data)
}


getBal <- function(dat_cases, dat_controls) {
  # combine data
  dat <- rbind(dat_cases, dat_controls)
  
  # # # balance age 
  # hist(dat$age_sample_collection[dat$type == 'cases'])
  # hist(dat$age_sample_collection[dat$type == 'controls'])
  
  # remove a few from ranges 100-200, 300-400
  # randomly remove controls that have a less than 50 month age of diganosis to have balanced classes
  remove_index <- which(dat$type == 'controls' & ((dat$age_sample_collection >= 100 & dat$age_sample_collection <= 200) |
                                                    (dat$age_sample_collection >= 300 & dat$age_sample_collection <= 400)))
  
  remove_index <- sample(remove_index, 10, replace = F )
  dat <- dat[-remove_index,]
  
}


bumpHunter<- function(model_data) {
  
  
 
  ##########
  # get clinical model_data 
  ##########
  bump_clin <- model_data[,1:4]
  
  ##########
  # get indicator and put into design matrix with intercept 1
  ##########
  indicator_vector <- as.factor(model_data$type)
  designMatrix <- cbind(rep(1, nrow(model_data)), indicator_vector)
  designMatrix <- as.matrix(designMatrix)
  
  ##########
  # Get genetic locations
  ##########
  model_data$age_diagnosis <- model_data$age_sample_collection <- 
    model_data$gender <- model_data$type <-  NULL
  # transpose methylation to join with cg_locations to get genetic location vector.
  model_data <- as.data.frame(t(model_data), stringsAsFactors = F)
  
  # make probe a column in methyl
  model_data$probe <- rownames(model_data)
  rownames(model_data) <- NULL
  
  # inner join methyl and cg_locations by probe
  methyl_cg <- inner_join(model_data, cg_locations, by = 'probe')
  
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
  DELTA_BETA_THRESH = .10 # DNAm difference threshold
  NUM_BOOTSTRAPS = 3   # number of randomizations
  
  # create tab list
  tab <- bumphunter(beta, 
                    designMatrix, 
                    chr = chr, 
                    pos = pos,
                    nullMethod = "bootstrap",
                    cutoff = DELTA_BETA_THRESH,
                    B = NUM_BOOTSTRAPS,
                    type = "Beta")
    
    bump_hunter_results <- tab$table
    bump_hunter_results$run <- DELTA_BETA_THRESH
  
  

  return(bump_hunter_results)
  
}


##########
# create function that grabs probe site and gene name for results from bumphunter
##########
# all bh features have some data with region length > 0. 
# cg_locations have no retions with length > 0.

getProbe <- function(data, cgs) {
  
  colnames(cgs) <- paste0(colnames(cgs), '_', 'rgSet')
  
  results <- list()
  results_data <- list()
  
  # first make seqnames in cgs and chr in tab character vectors for identification
  cgs$seqnames <- as.character(cgs$seqnames)
  data$chr <- as.character(data$chr)
  
  # use sql data table to merger validators with model_data based on age range
  result = sqldf("select * from cgs
                 inner join data
                 on cgs.start_rgSet between data.start and data.end")
  
  # keep only necessary columns
  result <- result[, c('chr' , 'start_rgSet','end_rgSet', 'probe_rgSet', 'p.value', 'fwer', 'run')]
  
  # rename cols
  colnames(result) <- c('chr', 'start', 'end', 'probe', 'p.value', 'fwer', 'run')
  
  # get sig results
  # result_sig <- result[result$p.value < 0.05,]
  
  # get fwer results
  result_fwer <- result[result$fwer == 0,]
  
  return(result_fwer)
  
}
