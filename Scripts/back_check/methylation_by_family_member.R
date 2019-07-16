
library(RColorBrewer)
library(dendextend)
# source all_functions.R to load libraries and my functions
source('../predict_age/all_functions.R')

##########
# get base directory for 4 batch
##########
path_to_cases_tor <- '../../Data/methyl_data/cases_toronto'
path_to_cases_mon <- '../../Data/methyl_data/cases_montreal'
path_to_controls <- '../../Data/methyl_data/controls'
path_to_valid <- '../../Data/methyl_data/validation'


##########
# read in meth array - Data/methyl_data/cases_toronto, cases_montreal, controls, validation
##########

# cases 
rgCasesT <- read.metharray.exp(path_to_cases_tor, recursive = T)
rgCasesM <- read.metharray.exp(path_to_cases_mon, recursive = T)

# combine cases arrays 
rgCases <- combineArrays(rgCasesT, rgCasesM)
rm(rgCasesT, rgCasesM)

# controls
rgControls <- read.metharray.exp(path_to_controls, recursive = T)

rgValid <- read.metharray.exp(path_to_valid, recursive = T)

##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
##########
ratio_set <- readRDS('../../Data/model_data/raw_ratio_set.rda')

# get granges object
g_ranges <- as.data.frame(getLocations(ratio_set))

# get probes from rownames
g_ranges$probe <- rownames(g_ranges)

# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]

##########
# read in clinical data
##########
clin <- read.csv('../../Data/clin_data/clinical_two.csv', stringsAsFactors = F)

# clean clinical ids
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

##########
# cases 
##########

# cases batch1
id_map_tor <- read.csv(paste0(path_to_cases_tor, '/SampleSheet.csv'), stringsAsFactors = F)

#cases batch2
id_map_mon <- read.csv(paste0(path_to_cases_mon, '/SampleSheet.csv'), stringsAsFactors = F)
id_map_mon$Project <- NULL

# combine id_map and id_map_other
id_map_cases <- rbind(id_map_tor, id_map_mon)
rm(id_map_tor, id_map_mon)

# clean id map
id_map_cases <- cleanIdMap(id_map_cases)

##########
# Controls batch1
##########
id_map_con <- read.csv(paste0(path_to_controls, '/SampleSheet.csv'), stringsAsFactors = F)

# clean idmap
id_map_con <- cleanIdMap(id_map_con)

##########
# valid
##########
id_map_val <- read.csv(paste0(path_to_valid, '/SampleSheet.csv'))

# homogenize valid map data with cases and controls
id_map_val <- id_map_val[, c('Sample.ID', 'Sample.Group', 'Sentrix.Barcode', 'Sample.Section',
                             'Project', 'Pool_ID', 'Sample_Well')]

# sub out '.' for '_'
colnames(id_map_val) <- gsub('.', '_', colnames(id_map_val), fixed = T)

# change 'Sample_ID' to 'Sample_Name' and 'Sentrix_Barcode' to 'Sentrix_ID'
colnames(id_map_val)[1] <- 'Sample_Name'
colnames(id_map_val)[3] <- 'Sentrix_ID'
colnames(id_map_val)[4] <- 'Sentrix_Position'
colnames(id_map_val)[5] <- 'Sample_Plate'

# clean idmap
id_map_val <- cleanIdMap(id_map_val)

##########
# remove outliers (previously determined) from rgset before normalization
##########
rgControls <- remove_outliers(rgSet = rgControls,
                              id_map = id_map_con,
                              method = 'doesnt_matter',
                              type = 'controls')

rgValid <- remove_outliers(rgSet = rgValid,
                           id_map = id_map_val,
                           method = 'doesnt_matter',
                           type = 'valid')



##########
# subset data - remove controls probes on each data set only if raw preprocessing
##########

# load in gene cpgs
gene_probes <- read_csv('../../Data/all_gene_cpg_loc.csv')

cg_gene_regions <-'Body'

gene_region <- paste(cg_gene_regions, collapse = '|')
# get probe gene region
gene_probes <- gene_probes[grepl(gene_region, gene_probes$focal_gene_regions),]

gene_probes <- as.character(gene_probes$focal_CpGs[!duplicated(gene_probes$focal_CpGs)])

# dont control for gender in model if using funnorm
method = 'noob'
keep_gender <- TRUE
keep_snps <- TRUE
keep_controls <- TRUE

# cases
rg_cases <- subset_rg_set(rg_set = rgCases,
                          keep_gender = keep_gender,
                          keep_controls = keep_controls,
                          keep_snps = keep_snps,
                          get_island = NULL,
                          get_chr = NULL,
                          get_type = NULL,
                          gene_probes = gene_probes)

# controls
rg_controls <- subset_rg_set(rg_set = rgControls,
                             keep_gender = keep_gender,
                             keep_controls = keep_controls,
                             keep_snps = keep_snps,
                             get_island = NULL,
                             get_chr = NULL,
                             get_type = NULL,
                             gene_probes = gene_probes)

# valid
rg_valid <- subset_rg_set(rg_set = rgValid,
                          keep_gender = keep_gender,
                          keep_controls = keep_controls,
                          keep_snps = keep_snps,
                          get_island = NULL,
                          get_chr = NULL,
                          get_type = NULL,
                          gene_probes = gene_probes)



# preprocess controls and valid
beta_cases <- preprocessMethod(rg_cases, preprocess = method)
beta_controls <- preprocessMethod(rg_controls, preprocess = method)
beta_valid <- preprocessMethod(rg_valid, preprocess = method)

# get controls
beta_cases <- process_rg_set_single(beta_data = beta_cases, 
                                    id_map = id_map_cases, 
                                    clin = clin)
# get controls
beta_controls_mod <- process_rg_set_single(beta_data = beta_controls, 
                                           id_map = id_map_con, 
                                           clin = clin)

# get valid
beta_valid_mod <- process_rg_set_single(beta_data = beta_valid, 
                                        id_map = id_map_val, 
                                        clin = clin)

# combine data 
intersecting_names <- Reduce(intersect, list(colnames(beta_cases)[12:ncol(beta_cases)],
                                             colnames(beta_controls_mod)[12:ncol(beta_controls_mod)],
                                             colnames(beta_valid_mod)[12:ncol(beta_valid_mod)]))

# clincal names
clin_names <- colnames(beta_cases)[1:11]

# get right variables 
beta_cases <- beta_cases[,c(clin_names, intersecting_names)]
beta_controls_mod <- beta_controls_mod[,c(clin_names, intersecting_names)]
beta_valid_mod <- beta_valid_mod[,c(clin_names, intersecting_names)]

# get model data in cases for training and test
beta_cases <- getModData(beta_cases)

# get rid of cancer samples in controls 
beta_controls_mod <- beta_controls_mod[grepl('Unaffected', beta_controls_mod$cancer_diagnosis_diagnoses),]

#subset valid - get ids from train and test
case_ids <- beta_cases$ids
beta_valid_mod <- beta_valid_mod[!beta_valid_mod$ids %in% case_ids,]

# remove NAs 
beta_cases <-beta_cases[complete.cases(beta_cases),]

##########
# clustering and heat map
##########

# get cases, controls, and valid
cases <- full_data[grepl('^57|^97', full_data$sentrix_id),]
controls <- full_data[grepl('^2003|2004', full_data$sentrix_id),]
valid <- full_data[grepl('^2009|^201', full_data$sentrix_id),]

# save.image('~/Desktop/temp_methyl.RData ')
# load('~/Desktop/temp_methyl.RData ')

# data should not be biased by any other variable - cancer, age, etc. So just use either cases or contorls. 
# find families that have similar distirbutions among shard variables. Also drop, for now, samples that are 
# single families to compare only samples that have familiy memebrs

# functon that gets a family indicator for each data set 
prepare_data <- function(temp_data) {
  # get family name vector
  family_name <- as.character(temp_data$family_name)
  
  # identify all dups and make indicator for family names if they duplicates and "single_family" if not to 
  # simplify the plotting
  dups_first <- ifelse(duplicated(family_name), 'dup', 'single_family')
  dups_last <- ifelse(duplicated(family_name, fromLast = TRUE), 'dup', 'single_family')
  temp <- as.data.frame(cbind(family_name, dups_first, dups_last))
  temp$new_dup <- ifelse(temp$dups_first == 'dup' | temp$dups_last == 'dup', 'def_dup', 'not_dup')
  temp$new_family_indicator <- ifelse(grepl('def_dup', temp$new_dup), temp$family_name, 'single_family')
  # get data into matrix format, with rownames as family indicator
  rownames(temp_data) <- gsub('[.]', '_', make.names(temp$new_family_indicator, unique=TRUE))
  temp_data <- temp_data[!grepl('single', rownames(temp_data)),]
  temp_data <- as.matrix(temp_data[, 12:ncol(temp_data)])
  
  return(temp_data)
}

cases <- prepare_data(cases)

cols <- rownames(cases)
cols <- unlist(lapply(strsplit(cols, '_'), function(x) x[1]))
colors <- brewer.pal(unique(cols), "Set2")

# prepare hierarchical cluster
dend <- dist(cases, method = 'euclidean') %>%  
  hclust(method = "complete") %>% 
  as.dendrogram %>%
  hang.dendrogram(hang_height=0.1) %>%
  set("labels_cex", 0.5) %>%
  color_branches(k=8, col=colors)

labels(dend) <- paste(as.character(rownames(cases))[order.dendrogram(dend)],"(",labels(dend),")",sep = "")

par(mar = c(3,3,3,7))

#plot clustering 

plot(dend, 
     main = "Clustered Methylation of Cancer Types", 
     horiz =  TRUE,  nodePar = list(cex = .007))
legend("topleft", legend = unique(cancer_type), fill = colors ,cex = 0.5)

# plot heatmap 

hmcols <- colorpanel(2750, "yellow", "black", "blue")

heatmap.2(as.matrix(beta_final),  
          main = "Methylation of LFS Cancer Types",
          srtCol = 60,
          key = FALSE,
          trace="none",
          dendrogram = "row",
          Rowv = dend,
          Colv = "NA", 
          col=hmcols,
          cexCol=0.3,
          cexRow=0.5)     


