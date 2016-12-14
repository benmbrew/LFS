##### This script maps cgp probe sites to nearest genome using the fdb.infiniumMethylation.hg19 data
# base, and primarily, the getNearestGene function.
# This is the third step in the pipeline
library(FDb.InfiniumMethylation.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)

home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')


#################################################################
# Read in cleaned methylation data and match with hm450
#################################################################

# Read in methylation data 
methylation <- read.csv(paste0(methyl_data, '/methylation.csv'), header = TRUE, check.names = FALSE)
# methylation_tumor <- read.csv(paste0(methyl_data, '/methylation_tumor.csv'), header = TRUE, check.names = FALSE)


# try different method of finding gene
anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other
anno <- as.data.frame(anno)

# make row names that contain probe a column 
anno$Probe <- rownames(anno)

# inner join probe anno and methylation by probe
methyl_gene <- inner_join(anno, methylation, by = 'Probe')

# keep only necessary columns 
methyl_gene$Forward_Sequence <- methyl_gene$SourceSeq <- methyl_gene$Random_Loci <- methyl_gene$Methyl27_Loci <- 
  methyl_gene$UCSC_RefGene_Accession <- methyl_gene$UCSC_RefGene_Group <- methyl_gene$Phantom <-
  methyl_gene$DMR <- methyl_gene$Enhancer <- methyl_gene$HMM_Island <- methyl_gene$Regulatory_Feature_Name <-
  methyl_gene$Regulatory_Feature_Group <- methyl_gene$DHS <- NULL

# lots of probes dont have corresponding genes. for the time being extract first gene and recode blanks
# as "no_nearby_gene"

names(methyl_gene)[1] <- 'gene'

gene_names <- strsplit(methyl_gene$gene, ';')
gene_names <- lapply(gene_names, function(x) x[(length(x) - length(x) +1)])
gene_names <- do.call(rbind, gene_names)

methyl_gene$gene <- as.factor(gene_names)

methyl_gene$Probe <- NULL

# group by gene and take mean 
methyl_summarised <- methyl_gene %>%
  filter(!is.na(gene)) %>%
  group_by(gene) %>%
  summarise_each(funs(mean))

###################################################################
# Transpose data and put in formate for analysis
###################################################################
col_names <- methyl_summarised$gene
methyl <- as.data.frame(t(methyl_summarised))
names(methyl)<- col_names
methyl <- cbind(x = rownames(methyl), methyl) 
methyl <- methyl[2:nrow(methyl),]
names(methyl)[1] <- 'id'
rownames(methyl) <- NULL
methyl[, 2:ncol(methyl)] <-
  apply(methyl[,2:ncol(methyl)], 2, function(x){as.numeric(as.character(x))})
methyl <- as.data.frame(methyl)

# drop duplicates from methylation so LSA work
methyl <- methyl[!duplicated(methyl$id),]
methyl <- methyl[!is.na(methyl$id),]
rownames(methyl) <- methyl[,1]
methyl <- methyl[, -1]


write.csv(methyl, paste0(methyl_data, '/methyl.csv'), row.names = TRUE)

# Load the 450k data from hg19 (bioconductor, library is FDb.InfiniumMethylation.hg19)
hm450 <- get450k()

# Get probe names from our methylation data  
probe_names <- as.character(methylation$Probe)
# probe_names_tumor <- as.character(methylation_tumor$Probe)
#probe_names <- 'cg09087961'


# remove probes that have less than 10 characters.
probe_names <- probe_names[nchar(probe_names) == 10]
# probe_names_tumor <- probe_names_tumor[nchar(probe_names_tumor) == 10]

# get probes from hm450
probes <- hm450[probe_names]
# probes_tumor <- hm450[probe_names_tumor]


#get the nearest gene to each probe location.
probe_info <- getNearestGene(probes)
probe_info <- cbind(probe = rownames(probe_info), probe_info)
rownames(probe_info) <- NULL

# probe_info_tumor <- getNearestGene(probes_tumor)
# probe_info_tumor <- cbind(probe = rownames(probe_info_tumor), probe_info_tumor)
# rownames(probe_info_tumor) <- NULL

# join probe_info with methylation. This keeps all of the probes that we could match in hm450 and drops the others.
methyl_gene <- left_join(probe_info, methylation, by = c('probe'= 'Probe'))
# methyl_gene_tumor <- left_join(probe_info_tumor, methylation_tumor, by = c('probe'= 'Probe'))


# Get rid of extra variables.
methyl_gene$probe <- NULL
methyl_gene$queryHits <- NULL
methyl_gene$subjectHits <- NULL
methyl_gene$distance<- NULL

# Get rid of extra variables.
# methyl_gene_tumor$probe <- NULL
# methyl_gene_tumor$queryHits <- NULL
# methyl_gene_tumor$subjectHits <- NULL
# methyl_gene_tumor$distance<- NULL

# # add _dup for duplicate ids. This is done so dplyr summarise_each will work 
# for(i in ncol(methyl_gene):1){
#   if(names(methyl_gene)[i] %in% names(methyl_gene)[duplicated(names(methyl_gene), fromLast = FALSE)]){
#      names(methyl_gene)[i] <- paste0(names(methyl_gene)[i], '_dup')
#     }
# }

# # add _dup for duplicate ids. This is done so dplyr summarise_each will work 
# for(i in ncol(methyl_gene_tumor):1){
#   if(names(methyl_gene_tumor)[i] %in% names(methyl_gene_tumor)[duplicated(names(methyl_gene_tumor), fromLast = FALSE)]){
#     names(methyl_gene_tumor)[i] <- paste0(names(methyl_gene_tumor)[i], '_dup')
#   }
# }

# get rid of duplicate columns names 
# temp <- methyl_gene[, !duplicated(colnames(methyl_gene))]

# group by gene and sum probe values. 
methyl_summarised <- methyl_gene %>%
  group_by(nearestGeneSymbol) %>%
  summarise_each(funs(mean))

# # group by gene and sum probe values. 
# methyl_summarised_tumor <- methyl_gene_tumor %>%
#   group_by(nearestGeneSymbol) %>%
#   summarise_each(funs(mean))

# # change _dup back to normal so there are duplicate IDs
# for(i in 1:ncol(methyl_summarised)){
#   if(grepl( '_dup',names(methyl_summarised)[i])){
#     split <- strsplit(names(methyl_summarised)[i], '_dup')
#     names(methyl_summarised)[i] <- split
#   }
# }

# # change _dup back to normal so there are duplicate IDs
# for(i in 1:ncol(methyl_summarised_tumor)){
#   if(grepl( '_dup',names(methyl_summarised_tumor)[i])){
#     split <- strsplit(names(methyl_summarised_tumor)[i], '_dup')
#     names(methyl_summarised_tumor)[i] <- split
#   }
# }

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

# drop duplicates from methylation so LSA work
methyl <- methyl[!duplicated(methyl$id),]
methyl <- methyl[!is.na(methyl$id),]
rownames(methyl) <- methyl[,1]
methyl <- methyl[, -1]


write.csv(methyl, paste0(methyl_data, '/methyl.csv'), row.names = TRUE)


# ###################################################################
# # Transpose data and put in formate for analysis
# ###################################################################
# col_names <- methyl_summarised_tumor$nearestGeneSymbol
# methyl_tumor <- as.data.frame(t(methyl_summarised_tumor))
# names(methyl_tumor)<- col_names
# methyl_tumor <- cbind(x = rownames(methyl_tumor), methyl_tumor) 
# methyl_tumor <- methyl_tumor[2:nrow(methyl_tumor),]
# names(methyl_tumor)[1] <- 'id'
# rownames(methyl_tumor) <- NULL
# methyl_tumor[, 2:ncol(methyl_tumor)] <-
#   apply(methyl_tumor[,2:ncol(methyl_tumor)], 2, function(x){as.numeric(as.character(x))})
# methyl_tumor <- as.data.frame(methyl_tumor)
# 
# # drop duplicates from methyl_tumoration so LSA work
# methyl_tumor <- methyl_tumor[!duplicated(methyl_tumor$id),]
# methyl_tumor <- methyl_tumor[!is.na(methyl_tumor$id),]
# rownames(methyl_tumor) <- methyl_tumor[,1]
# methyl_tumor <- methyl_tumor[, -1]
# 

# write.csv(methyl_tumor, paste0(methyl_data, '/methyl_tumor.csv'), row.names = TRUE)

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





