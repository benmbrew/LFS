##### This script maps cgp probe sites to nearest genome using the fdb.infiniumMethylation.hg19 data
# base, and primarily, the getNearestGene function.

library(minfi)
library(affy)
library(FDb.InfiniumMethylation.hg19)
library(dplyr)
# Initialize folders
home_folder <- '/home/benbrew/Documents'
project_folder <- paste(home_folder, 'LFS', sep = '/')
data_folder <- paste(project_folder, 'Data', sep = '/')

# set working directory and load data 
setwd(data_folder)
load('cleaned.RData')

# Load the 450k data from hg19
hm450 <- get450k()

# Get probe names 
probe_names <- as.character(methylation$probe)

# remove probes that have less than 10 characters.
probe_names <- probe_names[nchar(probe_names) ==10]

# subset probes that do not have cg
probes <- hm450[probe_names]

#get the nearest gene to each probe location.
probe_info <- getNearestGene(probes)
probe_info <- cbind(probe = rownames(probe_info), probe_info)
rownames(probe_info) <- NULL

# join probe_info with methylation. 
methyl_gene <- left_join(methylation, probe_info, by = 'probe')


