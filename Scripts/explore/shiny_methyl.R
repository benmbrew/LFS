##########
# initialize libraries
##########
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(biovizBase)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)
library(preprocessCore)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(methyAnalysis)
library(reshape2)
library(shinyMethyl)
library(biovizBase)
library(GenomicFeatures)
library(devtools)
library(minfiData)

registerDoParallel(1)

# LOOK AT CONTROL TYPES, PROBE TYPES, AND ALL FUNCTIONS IN MINFI
# https://bioconductor.org/packages/devel/bioc/manuals/minfi/man/minfi.pdf
# http://127.0.0.1:17958/library/shinyMethyl/doc/shinyMethyl.pdf
##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
idat_data <- paste0(methyl_data, '/raw_files')
idat_data_controls <- paste0(methyl_data, '/controls')
idat_data_valid <- paste0(methyl_data, '/validation/idat_files')
model_data <- paste0(data_folder, '/model_data')

###############

# Array: One sample
# Slide: Physical slide containing 12 arrays (6 × 2 grid)
# Plate: Physical plate containing at most 8 slides (96 arrays). For this tutorial, we use batch
# and plate interchangeably

#################

##########
# read in idate for cases, controls, and validation set
##########
rgCases <- read.metharray.exp(idat_data)
rgControls <- read.metharray.exp(idat_data_controls)
rgValid <- read.metharray.exp(idat_data_valid)

temp <- getProbeInfo(rgCases, type = 'I-Green')
getManifest(rgCases)

# large RGChannelSet


# what functions can be used on RGset

summary <- shinySummarize(rgCases)
summary_controls <- shinySummarize(rgControls)
summary_valid <- shinySummarize(rgValid)

# save.image('~/Desktop/shiny_methyl.RData')
# load('~/Desktop/shiny_methyl.RData')


runShinyMethyl(summary_controls)


##########
# get preprocedssing method
##########
betaCases <- preprocessMethod(rgCases, preprocess = method)
rm(rgCases)

browseVignettes("shinyMethyl")

# 
# Usage
# plotCpg(dat, cpg, pheno, type = c("categorical", "continuous"),
#         measure = c("beta", "M"), ylim = NULL, ylab = NULL, xlab = "",
#         fitLine = TRUE, mainPrefix = NULL, mainSuffix = NULL)
# Arguments
# dat An RGChannelSet, a MethylSet or a matrix. We either use the getBeta (or
# getM for measure="M") function to get Beta values (or M-values) (for the first
#                                                                                                                                           two) or we assume the matrix contains Beta values (or M-values).
# cpg A character vector of the CpG position identifiers to be plotted.
# pheno A vector of phenotype values.
# type Is the phenotype categorical or continuous?
# measure Should Beta values or log-ratios (M) be plotted?
# ylim y-axis limits.
# ylab y-axis label.
# xlab x-axis label.
# fitLine Fit a least-squares best fit line when using a continuous phenotype.
# mainPrefix Text to prepend to the CpG name in the plot main title.
# mainSuffix Text to append to the CpG name in the plot main title.


##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'raw'
k = 5

##########
# load data
##########
# read in full m value data 
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new_m.rda'))) #34 449936
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new_m.rda')))
#35 449783

##########
# read in cluster labels
##########
kmeans_lab <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_kmeans_labs.rda')))

###########
# make id into ids
###########
colnames(betaCases)[1] <- 'ids'
colnames(betaControls)[1] <- 'ids'
colnames(betaValid)[1] <- 'ids'

##########
# remove inf
##########
betaCases <- removeInf(betaCases, probe_start = 8)
betaControls <- removeInf(betaControls, probe_start = 8)
betaValid<- removeInf(betaValid, probe_start = 8)


# get old controls - Mut and 'Unaffected'
betaControlsOld <- subset(betaCases, p53_germline == 'Mut' & 
                            cancer_diagnosis_diagnoses == 'Unaffected')

# get p53, not 'Unaffected'
betaCases <- getModData(betaCases)

# get rid of cancer samples in controls 
betaControls <- betaControls[grepl('Unaffected', betaControls$cancer_diagnosis_diagnoses),]

#subset valid
betaValid <- betaValid[!betaValid$ids %in% betaCases$ids,]

##########
# get intersecting colnames and prepare data for modeling
##########

intersect_names <- Reduce(intersect, list(colnames(betaCases)[8:ncol(betaCases)], 
                                          colnames(betaControls)[8:ncol(betaControls)], 
                                          colnames(betaValid)[8:ncol(betaValid)]))
# assign dataframe identifier
betaCases$type <- '450k'
betaControls$type <- '850k'
betaControlsOld$type <- '450k'
betaValid$type <- '850k'


# cases
betaCases <- betaCases[, c('ids',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender',
                           'type',
                           intersect_names)]
# controls
betaControls <- betaControls[, c('ids',
                                 'age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender', 
                                 'type',
                                 intersect_names)]

# controls
betaControlsOld <- betaControlsOld[, c('ids',
                                       'age_diagnosis', 
                                       'age_sample_collection', 
                                       'cancer_diagnosis_diagnoses', 
                                       'gender', 
                                       'type',
                                       intersect_names)]

#validation
betaValid <- betaValid[, c('ids', 
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           'type',
                           intersect_names)]



# get controls full
betaControlsFull <- rbind(betaControls,
                          betaControlsOld)


# remove duplicates from betaControlsFull
length(which(duplicated(betaControlsFull$ids)))
betaControlsFull <- betaControlsFull[!duplicated(betaControlsFull$ids),]
# 
# function for plotcpg

plot_cpg <- function(beta, 
                     num_feats,
                     type, 
                     column_name, 
                     probe_start){
  
  
  mat <- as.matrix(t(beta[, probe_start:ncol(beta)]))
  cpg_list <- colnames(betaCases)[probe_start:ncol(betaCases)]
  cpg_list <- sample(cpg_list, num_feats, replace = T)
  pheno_dat <- betaCases[, 1:(probe_start -1)]
  
  plotCpg(mat, 
          cpg_list, 
          pheno_dat[, column_name], 
          type = type) 
  
}

plot_cpg(beta = betaCases, 
         num_feats = 100,
         type = 'continuous', 
         column_name = 'age_diagnosis', 
         probe_start = 7)

