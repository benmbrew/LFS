###############################################
# This script will run a PCA on the LFS methylation data and plot age of diagnosis.
# existing by gene methylation data we have. 
library(dplyr)
library(stringr)
library(impute)

# Initialize folders
home_folder <- '/home/benbrew/Documents'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(project_folder, '/Results')
setwd(data_folder)

if('methyl_lsa.RData' %in% dir()){
  
  load(paste0(data_folder, '/methyl_lsa.RData'))

  }else{

  results_file <- "file.txt"
  constructPath <- function(intermediate_folders, parent=results_folder,
                            file=results_file) {
    paste(parent, intermediate_folders, file, sep="/")
  }
  
  incomplete_file <- constructPath("Incomplete")
  imputed_file <- constructPath("Imputed")
  jvmGBLimit <- 8
  
  source(paste0(project_folder, '/Code/Functions/lsaImputation.R'))
  
  ################################################################
  # Read in methyl and clinical data and join by ids
  ################################################################
  
  # Read in data 
  methyl <- data.matrix(read.csv(paste0(methyl_data, '/methyl.csv'), stringsAsFactors = FALSE))
  clin <- read.csv(paste0(clin_data, '/clinical.csv'), stringsAsFactors = TRUE)
  
  # drop duplicates from methylation so LSA work
  methyl <- as.data.frame(methyl)
  methyl <- methyl[!duplicated(methyl$id),]
  methyl <- data.matrix(methyl)
  
  # Scale and impute NA with 0s
  scaleData <- function(matrix_data){ 
    id <- matrix_data[,1]
    scaled_matrix <- scale(matrix_data[,-1])  
    scaled_data <- cbind(id, scaled_matrix)
    scaled_data[,1] <- as.character(scaled_data[,1])
    return(scaled_data)
  }
  
  methyl <- scaleData(methyl)
  
  # remove na from ids
  methyl <- methyl[!is.na(methyl[,1]),]
  
  # make id rownames
  methyl <- apply(methyl, 2, function(x) as.numeric(x))
  rownames(methyl) <- methyl[,1]
  methyl <- methyl[, -1]
  
  #sub_methyl <- methyl[1:126, 1:10000]# 110 is the cutoff. Something about column 111.
  
  # run lsaImputaion of methylation data
  methyl_impute <- lsaImputation(incomplete_data = methyl, sample_rows = TRUE)
  
  
  # join rownames and methyl_impute and then erase rownames
  methyl_impute <- cbind(id = rownames(methyl_impute), methyl_impute)
  rownames(methyl_impute) <- NULL
  
  # remove everthing but numbers in the ids for both data sets 
  removePun <- function(data){
    replace <- str_replace_all(data$id, "[[:punct:]]", " ")
    split <- strsplit(replace, ' ')
    combine_split <- lapply(split, function(x) x[1])
    new_ids <- unlist(combine_split)
    return(data)
  }
  
  clin <- removePun(clin)
  methyl_impute <- removePun(methyl_impute)
  
  # inner_join clin
  full_data <- inner_join(clin, methyl_impute,
                          by = 'id')
  
  # PCA function that takes only methylation data
  exclude <- ncol(clin)
  pca <- function(data, exclude){
    pca <- prcomp(data[,-c(1:exclude)])
    return(pca)
  }
  
  # Run PCA
  pca_methyl <- pca(full_data, exclude)

  save.image(paste0(data_folder, '/methyl_lsa.RData'))

}

# make character vectors factors true and false
pcaPlot <- function(pca, 
                    data, 
                    clin_var, 
                    name,
                    numeric = TRUE,
                    PCA1 = 1,
                    PCA2 = 2){
  data <- data[!is.na(data[, clin_var]),]
  data <- data[!is.na(data$cancer_indicator),]
  data$cancer_indicator <- factor(data$cancer_indicator, levels = unique(data$cancer_indicator))
  if(numeric){
    
    data[, clin_var] <- (data[, clin_var])/mean(data[, clin_var])
    
    plot(pca$x[,PCA1], 
         pca$x[,PCA2],
         xlab = paste0('PCA', PCA1),
         ylab = paste0('PCA', PCA2),
         cex = ((data[, clin_var])/1.1),
         main = name,
         pch = 16,
         col = as.factor(data$cancer_indicator)
    )
    
    abline(v = c(0,0),
           h = c(0,0))
    legend('bottomleft',
           legend = unique(data$cancer_indicator),
           col=1:length(data$cancer_indicator),
           pch=16,
           cex = 0.7)


  }else{
    
    if(clin_var != 'cancer_indicator' && clin_var != 'cancer'){
      
      data[, clin_var] <- interaction(data[, clin_var], data$cancer_indicator)  
      data[, clin_var] <- as.character(data[, clin_var])
      data[, clin_var] <- ifelse(grepl('FALSE', data[, clin_var]), 'None', data[, clin_var])
      data[, clin_var] <- factor(data[ ,clin_var])
      data[, clin_var] <- factor(data[, clin_var], levels = unique(data[, clin_var]))
    }
    plot(pca$x[, PCA1], 
         pca$x[, PCA2],
         xlab = paste0('PCA', PCA1),
         ylab = paste0('PCA', PCA2),
         cex = 1,
         main = name,
         pch = 16,
         col = as.factor(data[, clin_var])
    )
    
    abline(v = c(0,0),
           h = c(0,0))
    
    legend('bottomleft',
           legend = unique(data[, clin_var]),
           col=1:length(data[, clin_var]),
           pch=16,
           cex = 0.7)
  
  }

}

plot(pca_methyl, type = 'l')

# Plot PCA 1 and 2
pcaPlot(pca_methyl, 
        data = full_data, 
        clin_var = 'age_of_onset', 
        name = 'age_of_onset')

pcaPlot(pca_methyl, 
        data = full_data, 
        clin_var = 'age_fac', 
        name = 'age_fac',
        numeric = FALSE)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'cancer', 
        name = 'cancer', 
        numeric = FALSE)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'cancer_indicator', 
        name = 'cancer_indicator', 
        numeric = FALSE)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'tp53', 
        name = 'tp53', 
        numeric = FALSE)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'gender', 
        name = 'gender', 
        numeric = FALSE)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'mdm2', 
        name = 'mdm2', 
        numeric = FALSE)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'pin_3', 
        name = 'pin_3', 
        numeric = FALSE)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'codon_72', 
        name = 'codon_72', 
        numeric = FALSE)

# Plot PCA 1 and 3
pcaPlot(pca_methyl, 
        data = full_data, 
        clin_var = 'age_of_onset', 
        name = 'age_of_onset',
        PCA1 = 1, 
        PCA2 = 3)

pcaPlot(pca_methyl, 
        data = full_data, 
        clin_var = 'age_fac', 
        name = 'age_fac',
        numeric = FALSE,
        PCA1 = 1,
        PCA2 = 3)


pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'cancer', 
        name = 'cancer', 
        numeric = FALSE,
        PCA1 = 1, 
        PCA2 = 3)


pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'cancer_indicator', 
        name = 'cancer_indicator', 
        numeric = FALSE,
        PCA1 = 1, 
        PCA2 = 3)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'tp53', 
        name = 'tp53', 
        numeric = FALSE,
        PCA1 = 1, 
        PCA2 = 3)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'gender', 
        name = 'gender', 
        numeric = FALSE,
        PCA1 = 1, 
        PCA2 = 3)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'mdm2', 
        name = 'mdm2', 
        numeric = FALSE,
        PCA1 = 1, 
        PCA2 = 3)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'pin_3', 
        name = 'pin_3', 
        numeric = FALSE, 
        PCA1 = 1, 
        PCA2 = 3)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'codon_72', 
        name = 'codon_72', 
        numeric = FALSE, 
        PCA1 = 1, 
        PCA2 = 3)

# Plot PCA 2 and 3
pcaPlot(pca_methyl, 
        data = full_data, 
        clin_var = 'age_of_onset', 
        name = 'age_of_onset',
        PCA1 = 2, 
        PCA2 = 3)


pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'cancer', 
        name = 'cancer', 
        numeric = FALSE,
        PCA1 = 2, 
        PCA2 = 3)


pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'cancer_indicator', 
        name = 'cancer_indicator', 
        numeric = FALSE,
        PCA1 = 2, 
        PCA2 = 3)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'tp53', 
        name = 'tp53', 
        numeric = FALSE,
        PCA1 = 2, 
        PCA2 = 3)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'gender', 
        name = 'gender', 
        numeric = FALSE,
        PCA1 = 2, 
        PCA2 = 3)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'mdm2', 
        name = 'mdm2', 
        numeric = FALSE,
        PCA1 = 2, 
        PCA2 = 3)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'pin_3', 
        name = 'pin_3', 
        numeric = FALSE, 
        PCA1 = 2, 
        PCA2 = 3)

pcaPlot(pca_methyl, data = full_data, 
        clin_var = 'codon_72', 
        name = 'codon_72', 
        numeric = FALSE, 
        PCA1 = 2, 
        PCA2 = 3)

