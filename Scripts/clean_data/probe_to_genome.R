##### This script maps cgp probe sites to nearest genome using the fdb.infiniumMethylation.hg19 data
# base, and primarily, the getNearestGene function.
library(FDb.InfiniumMethylation.hg19)
library(dplyr)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
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
probe_names <- probe_names[nchar(probe_names) == 10]

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

# add _dup for duplicate ids. This is done so dplyr summarise_each will work 
for(i in ncol(methyl_gene):1){
  if(names(methyl_gene)[i] %in% names(methyl_gene)[duplicated(names(methyl_gene), fromLast = FALSE)]){
     names(methyl_gene)[i] <- paste0(names(methyl_gene)[i], '_dup')
    }
}

# get rid of duplicate columns names 
# temp <- methyl_gene[, !duplicated(colnames(methyl_gene))]

# group by gene and sum probe values. 
methyl_summarised <- methyl_gene %>%
  group_by(nearestGeneSymbol) %>%
  summarise_each(funs(mean))

# change _dup back to normal so there are duplicate IDs
for(i in 1:ncol(methyl_summarised)){
  if(grepl( '_dup',names(methyl_summarised)[i])){
    split <- strsplit(names(methyl_summarised)[i], '_dup')
    names(methyl_summarised)[i] <- split
  }
}

###################################################################
# Transpose data and put in formate for analysis
###################################################################
col_names <- methyl_summarised$nearestGeneSymbol
methyl <- as.data.frame(t(methyl_summarised))
names(methyl)<- col_names
methyl <- cbind(x = rownames(methyl), methyl) 
methyl <- methyl[2:nrow(methyl),]
names(methyl)[1] <- 'id'
rownames(methyl) <- NULL
methyl[, 2:ncol(methyl)] <-
  apply(methyl[,2:ncol(methyl)], 2, function(x){as.numeric(as.character(x))})
methyl <- as.data.frame(methyl)

#

write.csv(methyl, paste0(methyl_data, '/methyl.csv'), row.names = FALSE)

# library(IlluminaHumanMethylation27k.db)
# ProbeToSymbol <- IlluminaHumanMethylation27kSYMBOL
# mapped_probes <- mappedkeys(ProbeToSymbol)
# mapped_probes_List <- as.list(ProbeToSymbol[mapped_probes])
# list <- do.call('rbind', mapped_probes_List)
# list <- cbind(rownames(list), list)
# rownames(list) <- NULL
# list <- as.data.frame(list)
# 
# same <- list[list$V2 %in% prprobe_info$nearestGeneSymbol,]
# probes <- probe_info[probe_info$nearestGeneSymbol == 'HOXD3',] 
# probes2 <- list[list$V2 == 'HOXD3',]





