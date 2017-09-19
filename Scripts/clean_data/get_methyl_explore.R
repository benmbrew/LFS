########################################
# this script will examine Mset 

# The first step is usually to preprocess the data, using a number of functions including
# 
# preprocessRaw() : do nothing.
# preprocessIllumina() : use Illumina’s standard processing choices.
# preprocessQuantile() : use a version of quantile normalization adapted to methylation arrays.
# preprocessNoob() : use the NOOB background correction method.
# preprocessSWAN() : use the SWAN method.
# preprocessFunnorm() : use functional normalization.
# These functions output different types of objects.
# 
# The class hierarchy in minfi is as follows: data can be stored in an Methylation and 
# Unmethylation channel or in a percent methylation (called Beta) channel. For the first case we have 
# the class MethylSet, for the second case we have the class RatioSet. When you have methylation / unmethylation 
# values you can still compute Beta values on the fly. You convert from a MethylSet to a RatioSet with  ratioConvert().
# 
# In addition to these two classes, we have GenomicMethylSet and GenomicRatioSet. 
# The Genomic indicates that the data has been associated with genomic coordinates using the mapToGenome() function.
# 
# The starting point for most analyses ought to be a GenomicRatioSet class. 
# If your preprocessing method of choice does not get you there, use ratioConvert() and mapToGenome() to go the last steps.
# 

# two minfi objects - 
# 1) genomic ratio set - data organized by the CpG locus level, but not mapped to a genome. The data has at least one of two
# channels: Beta and/or M (logratio of Beta). It may optionally include a CN channel (copy number), mapped to genome.
# 2) genomic methyl set - data organized by the CpG locus level, but not mapped to a genome. This data has two channels:
# Meth (methylated) and Unmeth (unmethylated), mapped to genome.


##########
# load libraries
##########
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(biovizBase)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)
library(preprocessCore)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
idat_data <- paste0(methyl_data, '/raw_files')
idat_data_con <- paste0(methyl_data, '/controls')
idat_data_valid <- paste0(methyl_data, '/validation/idat_files')
model_data <- paste0(data_folder, '/model_data')


##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
dat = 'valid'
method = 'raw'


##########
# read in idat 
##########t
if(dat == 'cases') {
  rgSet <- read.metharray.exp(idat_data)
  
} else if(dat == 'controls') {
  rgSet <- read.metharray.exp(idat_data_con)
  
} else {
  rgSet <- read.metharray.exp(idat_data_valid)
  
}




##########
# get beta, m, gset, mset - mset and gset might be the same if funnorm
##########
methyl_list <- preprocessMethod(rgSet, preprocess = method)
beta_val <- methyl_list[[1]]
m_val <- methyl_list[[2]]
g_set <- methyl_list[[3]]
m_set <- methyl_list[[4]]


##########
# pheno data
##########
getManifest(rgSet)
getAnnotation(rgSet)

##########
# annotation
##########
anno <- getAnnotation(methyl_set)
names(anno)
anno$DMR
unique(getProbeType(methyl_set))
# only I and II only

# quality control (uses mset mostly)

# qc plot
qc <- getQC(methyl_set)
head(qc)
plotQC(qc, badSampleCutoff = 11.5)


# density plot looks at beta values
densityPlot(methyl_set)
temp_beta <- getBeta(methyl_set)
dim(temp_beta)
temp_beta_2 <- temp_beta[1:10000, 1:50]

densityPlot(temp_beta_2)

# control probes plots

# SNPS

# cross reactive probes

###########
# analysis
###########

# diff methylated positions

# diff methylatied regions




