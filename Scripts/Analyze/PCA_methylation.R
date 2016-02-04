###############################################
# This script will run a PCA on the LFS methylation data and plot age of diagnosis.
# existing by gene methylation data we have. 
library(dplyr)
library(stringr)
library(impute)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/Analyze')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

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
  
  # Read in data (clinical or clinical_two)
  methyl <- data.matrix(read.csv(paste0(methyl_data, '/methyl.csv'), stringsAsFactors = FALSE))
  clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = TRUE)
  
  # put ids in rownames for imputation
  rownames(methyl) <- methyl[,1]
  methyl <- methyl[, -1]
  
  # run lsaImputaion of methylation data
  methyl_impute_raw <- lsaImputation(incomplete_data = methyl, sample_rows = TRUE)

  # join rownames and methyl_impute and then erase rownames
  methyl_impute_raw <- cbind(id = rownames(methyl_impute_raw), methyl_impute_raw)
  rownames(methyl_impute_raw) <- NULL
  
  # make clin id a factor so it joins with methylation data
  clin$id <- as.factor(clin$blood_dna_malkin_lab_)
  
  # inner_join clin
  full_data <- inner_join(clin, methyl_impute,
                          by = 'id')
  
  # Save data to be used later
  write.csv(full_data, paste0(data_folder, '/full_data.csv'))
  write.csv(methyl_impute_raw, paste0(data_folder, '/methyl_impute_raw.csv'))
  write.csv(clin, paste0(data_folder, '/clin.csv'), row.names = FALSE)
  
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

pcaPlot <- function(pca, 
                    data, 
                    clin_var, 
                    name,
                    PCA1 = 1,
                    PCA2 = 2){
  data <- data[!is.na(data[, clin_var]),]
  data[, clin_var] <- factor(data[, clin_var], levels = unique(data[, clin_var]))
  
          
    plot(pca$x[,PCA1], 
         pca$x[,PCA2],
         xlab = paste0('PCA', PCA1),
         ylab = paste0('PCA', PCA2),
         cex = 1, #((data[, clin_var])/1.1),
         main = name,
         pch = 16,
         col = as.factor(data[, clin_var])
    )
    
    abline(v = c(0,0),
           h = c(0,0))
    legend('bottomleft',
           legend = unique(data[,clin_var]),
           col=1:length(data[,clin_var]),
           pch=16,
           cex = 0.7)
  }
plot(pca_methyl, type = 'l')

pcaPlot(pca_methyl,
        data = full_data,
        clin_var = 'p53_germline',
        name = 'P53',
        PCA1 = 1,
        PCA2 = 3)


full_data$cancer <- ifelse(full_data$cancer_diagnosis_diagnoses != 'Unaffected', TRUE, FALSE)
pcaPlot(pca_methyl, 
        data = full_data, 
        clin_var = 'cancer', 
        name = 'cancer',
        PCA1 = 1, 
        PCA2 = 3)

pcaPlot(pca_methyl, 
        data = full_data, 
        clin_var = 'cancer_diagnosis_diagnoses', 
        name = 'cancer_diagnosis_diagnoses',
        PCA1 = 1, 
        PCA = 3)

full_data$acc <- ifelse(full_data$cancer_diagnosis_diagnoses == 'ACC', TRUE, FALSE)
pcaPlot(pca_methyl, 
        data = full_data, 
        clin_var = 'acc', 
        name = 'acc',
        PCA1 = 1,
        PCA2 =3)

