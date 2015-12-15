##### This script looks at different ways in which you can map probe level methylation data to 
# genome level using a variety of R packages..49/34

library(minfi)
library(affy)
library(ArrayExpress)
library(annotate)
library(hgu133a.db)
# Initialize folders
home_folder <- '/home/benbrew/Documents'
project_folder <- paste(home_folder, 'LFS', sep = '/')
data_folder <- paste(project_folder, 'Data', sep = '/')

# set working directory and load data 
setwd(data_folder)
load('cleaned.RData')

# transopse methylation data so that probes are rows and columns are ids
methylation_raw <- as.data.frame(t(methylation_raw), stringsAsFactors = FALSE)

PROBES<- as.character(FCMATRIX$probe)

OUT <- select(rat2302.db, PROBES, c("SYMBOL", "ENTREZID", "GENENAME"))

# Make the first row column names
colnames(methylation_raw) <- methylation_raw[1,]
data <- methylation_raw <- methylation_raw[-1,]
# data <- cbind(id = rownames(methylation_raw), methylation_raw)
# rownames(data) <- NULL
data <- data[, unique(colnames(data))]

# apply the minfi functions to transform data 
# object = Either a MethylSet, a RGChannelSet or a RatioSet.

data <- data[1:100,]
annot_data <- AnnotatedDataFrame(data)

getLocations(annot_data)


