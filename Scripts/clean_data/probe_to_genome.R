##### This script maps cgp probe sites to nearest genome using the fdb.infiniumMethylation.hg19 data
# base, and primarily, the getNearestGene function.
library(data.table)
library(FDb.InfiniumMethylation.hg19)
library(dplyr)

# Initialize folders
home_folder <- '/home/benbrew/Documents'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')

#################################################################
# Read in cleaned methylation data and match with hm450
#################################################################

# Read in methylation data 
methylation <- read.csv(paste0(methyl_data, '/methylation.csv'), header = TRUE, check.names = FALSE)

# Load the 450k data from hg19 (bioconductor, library is FDb.InfiniumMethylation.hg19)
hm450 <- get450k()

# Get probe names from our methylation data  
probe_names <- as.character(methylation$Probe)

# remove probes that have less than 10 characters.
probe_names <- probe_names[nchar(probe_names) ==10]

# get probes from hm450
probes <- hm450[probe_names]

#get the nearest gene to each probe location.
probe_info <- getNearestGene(probes)
probe_info <- cbind(probe = rownames(probe_info), probe_info)
rownames(probe_info) <- NULL

# join probe_info with methylation. This keeps all of the probes that we could match in hm450 and drops the others.
methyl_gene <- left_join(probe_info, methylation, by = c('probe'= 'Probe'))

# Get rid of extra variables.
methyl_gene$probe <- NULL
methyl_gene$queryHits <- NULL
methyl_gene$subjectHits <- NULL
methyl_gene$distance<- NULL

###################################################################
# Transpose data and put in formate for analysis
###################################################################
# break data into 100 data sets into a list 

n <- methyl_gene$nearestGeneSymbol
temp <- as.data.frame(t(methyl_gene))
names(temp)<- n
temp <-








