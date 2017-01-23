####### Script will load prepared data and run pca and visualize
# original idat, and controls

##########
# initialize libraries
##########
library(dplyr)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')

##########
# load cases and controls data
##########
beta_raw <- readRDS(paste0(methyl_data, '/beta_raw.rda'))
beta_raw_controls <- readRDS(paste0(methyl_data, '/beta_raw_controls.rda'))

beta_swan <- readRDS(paste0(methyl_data, '/beta_swan.rda'))
beta_swan_controls <- readRDS(paste0(methyl_data, '/beta_swan_controls.rda'))

beta_quan <- readRDS(paste0(methyl_data, '/beta_quan.rda'))
beta_quan_controls <- readRDS(paste0(methyl_data, '/beta_quan_controls.rda'))

beta_funnorm <- readRDS(paste0(methyl_data, '/beta_funnorm.rda'))
beta_funnorm_controls <- readRDS(paste0(methyl_data, '/beta_funnorm_controls.rda'))

########## 
# PCA of each data type and cases vs controls
##########

# subset control data by features in normal data.

# function needs to take a clinical column, remove others, and plot pcas
getPCA <- function(pca_data, column_name, name) 
{
  # subet data so only p53 mut
  pca_data <- pca_data[pca_data$p53_germline == 'Mut',]
  
  
  # get features sites
  cg_sites <- colnames(pca_data)[14:ncol(pca_data)]
  
  # subset by no NAs for column_name
  pca_data <- pca_data[!is.na(pca_data[, column_name]), ]
  
  stopifnot(!any(is.na(pca_data[, column_name])))
  
  # put column name with cg_sites 
  pca_data <- pca_data[ ,c(column_name, cg_sites)]
  
  # run pca
  data_length <- ncol(pca_data)
  pca <- prcomp(pca_data[,2:data_length])
  
  # plot data
  #fill in factors with colors 
  col_vec <- c('red', 'green', 'blue', 'orange', 'black', 'orange')
  colors <- col_vec[pca_data[, column_name]]
  min_x <- min(pca$x[,1])
  max_x <- max(pca$x[,1])
  min_y <- min(pca$x[,2])
  max_y <- max(pca$x[,2])
  
  max <- max(pca$x[, 2])
  plot <- plot(pca$x[,1], 
               pca$x[,2],
               xlab = 'pca 1',
               ylab = 'pca 2',
               cex = 1,
               main = name,
               pch = 16,
               xlim= c(min_x, max_x),
               ylim = c(min_y, max_y),
               col = adjustcolor(colors, alpha.f = 0.5)
  )
  abline(v = c(0,0),
         h = c(0,0))
  
  return(plot)
}

##########
# use gender
##########
# raw
plot_raw_gender <- getPCA(beta_raw, 'gender', 'PCA raw Gender')

# swan
plot_swan_gender <- getPCA(beta_swan, 'gender', 'PCA swan Gender')

# quan
plot_quan_gender <- getPCA(beta_quan, 'gender', 'PCA quan Gender')

# funnorm
plot_funnorm_gender <-getPCA(beta_funnorm, 'gender', 'PCA funnorm Gender')


##########
# use cancer diagnosis
##########
# raw
plot_raw_canc <- getPCA(beta_raw, 'cancer_diagnosis_diagnoses', 'PCA raw cancer_diagnosis_diagnoses')

# swan
plot_swan_canc <- getPCA(beta_swan, 'cancer_diagnosis_diagnoses', 'PCA swan cancer_diagnosis_diagnoses')

# quan
plot_quan_canc <- getPCA(beta_quan, 'cancer_diagnosis_diagnoses', 'PCA quan cancer_diagnosis_diagnoses')

# funnorm
plot_funnorm_canc <-getPCA(beta_funnorm, 'cancer_diagnosis_diagnoses', 'PCA funnorm cancer_diagnosis_diagnoses')


##########
# use gdna.exon.intron
##########
# raw
getPCA(beta_raw, 'gdna.exon.intron', 'PCA raw gdna.exon.intron')

# swan
getPCA(beta_swan, 'gdna.exon.intron', 'PCA swan gdna.exon.intron')

# quan
getPCA(beta_quan, 'gdna.exon.intron', 'PCA quan gdna.exon.intron')

# funnorm
getPCA(beta_funnorm, 'gdna.exon.intron', 'PCA funnorm gdna.exon.intron')


##########
# use gdna.base.change
##########
# raw
getPCA(beta_raw, 'gdna.base.change', 'PCA raw gdna.base.change')

# swan
getPCA(beta_swan, 'gdna.base.change', 'PCA swan gdna.base.change')

# quan
getPCA(beta_quan, 'gdna.base.change', 'PCA quan gdna.base.change')

# funnorm
getPCA(beta_funnorm, 'gdna.base.change', 'PCA funnorm gdna.base.change')


##########
# use protein.codon.change
##########
# raw
getPCA(beta_raw, 'protein.codon.change', 'PCA raw protein.codon.change')

# swan
getPCA(beta_swan, 'protein.codon.change', 'PCA swan protein.codon.change')

# quan
getPCA(beta_quan, 'protein.codon.change', 'PCA quan protein.codon.change')

# funnorm
getPCA(beta_funnorm, 'protein.codon.change', 'PCA funnorm protein.codon.change')

##########
# use protein.codon.num
##########
# raw
getPCA(beta_raw, 'protein.codon.num', 'PCA raw protein.codon.num')

# swan
getPCA(beta_swan, 'protein.codon.num', 'PCA swan protein.codon.num')

# quan
getPCA(beta_quan, 'protein.codon.num', 'PCA quan protein.codon.num')

# funnorm
getPCA(beta_funnorm, 'protein.codon.num', 'PCA funnorm protein.codon.num')


##########
# use splice.delins.snv
##########
# raw
getPCA(beta_raw, 'splice.delins.snv', 'PCA raw splice.delins.snv')

# swan
getPCA(beta_swan, 'splice.delins.snv', 'PCA swan splice.delins.snv')

# quan
getPCA(beta_quan, 'splice.delins.snv', 'PCA quan splice.delins.snv')

# funnorm
getPCA(beta_funnorm, 'splice.delins.snv', 'PCA funnorm splice.delins.snv')


##########
# use codon72.npro
##########
# raw
getPCA(beta_raw, 'codon72.npro', 'PCA raw codon72.npro')

# swan
getPCA(beta_swan, 'codon72.npro', 'PCA swan codon72.npro')

# quan
getPCA(beta_quan, 'codon72.npro', 'PCA quan codon72.npro')

# funnorm
getPCA(beta_funnorm, 'codon72.npro', 'PCA funnorm codon72.npro')


##########
# use mdm2.nG
##########
# raw
getPCA(beta_raw, 'mdm2.nG', 'PCA raw mdm2.nG')

# swan
getPCA(beta_swan, 'mdm2.nG', 'PCA swan mdm2.nG')

# quan
getPCA(beta_quan, 'mdm2.nG', 'PCA quan mdm2.nG')

# funnorm
getPCA(beta_funnorm, 'mdm2.nG', 'PCA funnorm mdm2.nG')


##########
# combine cases and controls and plot based on that
##########

