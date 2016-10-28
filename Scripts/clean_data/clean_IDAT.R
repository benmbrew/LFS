## Script for cleaning IDAT files
# this data could be used to replace methylation from nardine
# this script will be documented heavily as to provide a guide for IDAT files in minfi package to lab.

# initialize libraries
library(minfi)
library(dplyr)

# Use this https://www.bioconductor.org/help/course-materials/2014/BioC2014/minfi_BioC2014.pdf

# initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/classification_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
idat_data <- paste0(methyl_data, '/raw_files')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


####################
# Load id_map to later merge with methylation values. clean it.
####################
id_map <- read.csv(paste0(methyl_data, '/ids_map.csv'), stringsAsFactors = F)

# remove top rows 
id_map <- id_map[-c(1:6),]

# new colnames, lowercase
colnames(id_map) <- tolower(id_map[1,])
id_map <- id_map[-1,]

# combine sentrix_id and sentrix_position 
id_map$identifier <- paste(id_map$sentrix_id, id_map$sentrix_position, sep = '_')
id_map$identifier <- as.factor(id_map$identifier)


####################
# this function will read in each idat file using read.450k.exp from minifi package and store them in a list
####################

rgSetList <- list()
readIDAT <- function(path) {
  
  rgSet <- list()
  directory <- dir(path)
  
  for (folder in directory){
    # read into rgSet
    rgSet[[folder]] <- read.metharray.exp(paste(path, folder, sep = '/'))
  }
  
  return(rgSet)
  
}

rgSetList <- readIDAT(idat_data)

####################
# Loop through list and convert to beta, m, and cn values 
####################
Mset <- list()
ratioSet <- list()
gset <- list()
beta <- list()
m <- list()
cn <- list()

for (dat in 1:length(rgSetList)) {
  Mset[[dat]] <-preprocessRaw(rgSetList[[dat]])
  ratioSet[[dat]] <- ratioConvert(Mset[[dat]], what = 'both', keepCN = T)
  gset[[dat]] <- mapToGenome(ratioSet[[dat]]) 
  beta[[dat]] <- getBeta(gset[[dat]])
  m[[dat]] <- getM(gset[[dat]])
  cn[[dat]] <- getCN(gset[[dat]])
  print(dat)
}

####################
# Function that combines elements of a list into a data frame
####################

combineList <- function(list) {
  
  for (dat in 1:length(list)) {
    list[[dat]] <- t(list[[dat]])
  }
  data <- do.call(rbind, list)
  return(data)
}

beta_values <- combineList(beta)
m_values <- combineList(m)
cn_values <- combineList(cn)


####################
# Function that combines methylation matrices with id_map, to get ids for methylation
####################

findIds <- function(data) {
  data <- as.data.frame(data)
  data$identifier <- rownames(data)
  data$identifier <- as.factor(data$identifier)
  beta_methyl <- inner_join(beta_values, id_map, by = 'identifier')
}

beta_methyl <- findIds(beta_values)

###########################
# Preprocessing and normilization
###########################

# As seen before, it converts a RGChannelSet to a MethylSet by converting the Red and
# Green channels into a matrix of methylated signals and a matrix of unmethylated signals.
# No normalization is performed.

####### PreprocessIllumina
# Convert a RGChannelSet to a MethylSet by implementing the preprocessing choices as
# available in Genome Studio: background subtraction and control normalization. Both of
# them are optional and turning them off is equivalent to raw preprocessing (preprocessRaw):
MSet.illumina <- preprocessIllumina(RGSet, bg.correct = TRUE,
                                    normalize = "controls")


####### PreprocessSWAN

# Perform Subset-quantile within array normalization (SWAN) [6], a within-array normalization
# correction for the technical differences between the Type I and Type II array designs.
# The algorithm matches the Beta-value distributions of the Type I and Type II probes by
# applying a within-array quantile normalization separately for different subsets of probes (divided
#                                                                                            by CpG content). The input of SWAN is a MethylSet, and the function returns
# a MethylSet as well. If an RGChannelSet is provided instead, the function will first call
# preprocessRaw on the RGChannelSet, and then apply the SWAN normalization. We recommend
# setting a seed (using set.seed) before using preprocessSWAN to ensure that the
# normalized intensities will be reproducible.

MSet.swan <- preprocessSWAN(RGSet)


######## PreprocessQuantile

# This function implements stratified quantile normalization preprocessing. The normalization
# procedure is applied to the Meth and Unmeth intensities separately. The distribution of type
# I and type II signals is forced to be the same by first quantile normalizing the type II probes
# across samples and then interpolating a reference distribution to which we normalize the
# type I probes. Since probe types and probe regions are confounded and we know that
# DNAm distributions vary across regions we stratify the probes by region before applying
# this interpolation. Note that this algorithm relies on the assumptions necessary for quantile
# normalization to be applicable and thus is not recommended for cases where global changes
# are expected such as in cancer-normal comparisons. Note that this normalization procedure
# is essentially similar to one previously presented [7]. The different options can be summarized
# into the following list:

# 1) If fixMethOutlier is TRUE, the functions fixes outliers of both the methylated and
# unmethylated channels when small intensities are close to zero.
# 2) If removeBadSamples is TRUE, it removes bad samples using the QC criterion discussed
# previously
# 3) Performs stratified subset quantile normalization if quantileNormalize=TRUE and
# stratified=TRUE
# 4) Predicts the sex (if not provided in the sex argument) using the function getSex and
# normalizes males and females separately for the probes on the X and Y chromosomes

gset.quantile <- preprocessQuantile(RGSet, fixOutliers = TRUE,
                                    removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                    quantileNormalize = TRUE, stratified = TRUE,
                                    mergeManifest = FALSE, sex = NULL)


############## PreprocessFunnorm
# The function preprocessFunnorm implements the functional normalization algorithm developed
# in [8]. Briefly, it uses the internal control probes present on the array to infer
# between-array technical variation. It is particularly useful for studies comparing conditions
# with known large-scale differences, such as cancer/normal studies, or between-tissue studies.
# It has been shown that for such studies, functional normalization outperforms other existing
# approaches [8]. By default, is uses the first two principal components of the control probes
# to infer the unwanted variation.

# gset.funnorm <- preprocessFunnorm(RGSet)

