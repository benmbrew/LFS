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

load(paste0(data_folder, '/methyl_lsa.RData'))


# PCA function that takes only methylation data
exclude <- ncol(clin)
pca <- function(data, exclude){
  pca <- prcomp(data[,-c(1:exclude)])
  return(pca)
}

# Run PCA
pca_methyl <- pca(full_data, exclude)


pcaPlotFac <- function(pca, 
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
  legend('bottomright',
         legend = unique(data[,clin_var]),
         col=1:length(data[,clin_var]),
         pch=16,
         cex = 0.7)
}
plot(pca_methyl, type = 'l')

pcaPlotFac(pca_methyl,
      data = full_data,
      clin_var = 'p53_germline',
      name = 'P53',
      PCA1 = 1,
      PCA2 = 3)


full_data$cancer <- ifelse(full_data$cancer_diagnosis_diagnoses != 'Unaffected', TRUE, FALSE)
pcaPlotFac(pca_methyl, 
      data = full_data, 
      clin_var = 'cancer', 
      name = 'cancer',
      PCA1 = 1, 
      PCA2 = 3)

pcaPlotFac(pca_methyl, 
      data = full_data, 
      clin_var = 'cancer_diagnosis_diagnoses', 
      name = 'cancer_diagnosis_diagnoses',
      PCA1 = 1, 
      PCA = 3)

full_data$acc <- ifelse(full_data$cancer_diagnosis_diagnoses == 'ACC', TRUE, FALSE)
pcaPlotFac(pca_methyl, 
      data = full_data, 
      clin_var = 'acc', 
      name = 'acc',
      PCA1 = 1,
      PCA2 =3)

full_data$age_diagnosis <- as.numeric(as.character(full_data$age_diagnosis))
full_data$age_six <- ifelse(full_data$age_diagnosis > 6, TRUE, FALSE)
pcaPlotFac(pca_methyl, 
         data = full_data, 
         clin_var = 'age_six', 
         name = 'age_6',
         PCA1 = 1,
         PCA2 =3)


pcaPlotNum <- function(pca, 
                     data, 
                     clin_var, 
                     name,
                     PCA1 = 1,
                     PCA2 = 2){
data <- data[!is.na(data[, clin_var]),]
data[, clin_var] <- (data[, clin_var] - mean(data[, clin_var]))/10


plot(pca$x[,PCA1], 
     pca$x[,PCA2],
     xlab = paste0('PCA', PCA1),
     ylab = paste0('PCA', PCA2),
     cex = ((data[, clin_var])/1.1),
     main = name,
     pch = 16,
     col = adjustcolor('blue', alpha.f = 0.6)
)

abline(v = c(0,0),
       h = c(0,0))
#   legend('bottomright',
#          legend = unique(data[,clin_var]),
#          col=1:length(data[,clin_var]),
#          pch=16,
#          cex = 0.7)
}
plot(pca_methyl, type = 'l')


pcaPlotNum(pca_methyl, 
         data = full_data, 
         clin_var = 'age_diagnosis', 
         name = 'age_diagnosis',
         PCA1 = 1,
         PCA2 =3)
