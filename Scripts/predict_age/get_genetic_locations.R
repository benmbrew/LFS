# this scirpt will retrieve genetic locations using a 450k 
library(minfi)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)



# load data 
getGEOSuppFiles("GSE68777")
untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")
head(list.files("GSE68777/idat", pattern = "idat"))

#idat files
idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

# read into rgSet
rgSet <- read.450k.exp("GSE68777/idat")
rgSet

# apply pData function
pData(rgSet)
## data frame with 0 columns and 40 rows
head(sampleNames(rgSet))

# The rgSet object is a class called RGChannelSet which represents two color data with a green and a red channel, 
# very similar to an ExpressionSet.
# 
# The first step is usually to preprocess the data, using a number of functions including
# 
# preprocessRaw() : do nothing.
# preprocessIllumina() : use Illumina’s standard processing choices.
# preprocessQuantile() : use a version of quantile normalization adapted to methylation arrays.
# preprocessNoob() : use the NOOB background correction method.
# preprocessSWAN() : use the SWAN method.
# preprocessFunnorm() : use functional normalization.
# These functions output different types of objects.

# The class hierarchy in minfi is as follows: data can be stored in an Methylation and Unmethylation 
# channel or in a percent methylation (called Beta) channel. For the first case we have the class MethylSet, 
# for the second case we have the class RatioSet. When you have methylation / unmethylation values you can still 
# compute Beta values on the fly. You convert from a MethylSet to a RatioSet with ratioConvert().

# In addition to these two classes, we have GenomicMethylSet and GenomicRatioSet. The Genomic indicates that the 
# data has been associated with genomic coordinates using the mapToGenome() function.
# 
# The starting point for most analyses ought to be a GenomicRatioSet class. If your preprocessing method of choice 
# does not get you there, use ratioConvert() and mapToGenome() to go the last steps.
# 
# Let us run preprocessQuantile() which arrives at a GenomicRatioSet:
#   
grSet <- preprocessQuantile(rgSet)
grSet

# This is like a SummarizedExperiment; we can get the location of the CpGs by
temp <- granges(grSet)



# The usual methylation measure is called “Beta” values; equal to percent methylation and 
# defined as Meth divided by Meth + Unmeth.

getBeta(grSet)[1:3,1:3]


# CpGs forms clusters known as “CpG Islands”. Areas close to CpG Islands are known as CpG Shores, followed by CpG 
# Shelfs and finally CpG Open Sea probes. An easy way to get at this is to use

head(getIslandStatus(grSet))










