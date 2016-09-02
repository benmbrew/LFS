#######################################################################################
# This script will run bumphunter on methylation probes 
library(minfi)
library(bumphunter)
library(dplyr)
library(FDb.InfiniumMethylation.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(impute)


# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/analyze')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

data_name <- '/Chr'



# # load in idat human methylation 450k kanno
# getGEOSuppFiles("GSE68777")
# untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")
# head(list.files("GSE68777/idat", pattern = "idat"))
# 
# #idat files
# idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
# sapply(idatFiles, gunzip, overwrite = TRUE)
# 
# # read into rgSet
# rgSet <- read.450k.exp("GSE68777/idat")
# grSet <- preprocessQuantile(rgSet)
# 
# 
# # get cg site locations 
# cg_locations <- as.data.frame(granges(grSet))

# write.csv(cg_locations, paste0(data_folder, '/cg_locations'))


#
