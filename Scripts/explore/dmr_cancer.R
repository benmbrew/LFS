
# source all_functions.R to load libraries and my functions
source('../predict_age/all_functions.R')

library(RColorBrewer)
library(ComplexHeatmap)
library(dendextend)
library(colorspace)
library(dendextend)
library(densityP)


##########
# get base directory for 4 batch
##########
path_to_cases_tor <- '../../Data/methyl_data/cases_toronto'
path_to_cases_mon <- '../../Data/methyl_data/cases_montreal'
path_to_controls <- '../../Data/methyl_data/controls'
path_to_valid <- '../../Data/methyl_data/validation'

# set method
method = 'raw'
controls = 'full'
beta_thresh = 0.01

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

# cases
rg_cases <- subset_rg_set(rg_set = rgCases,
                          keep_gender = F,
                          keep_controls = F,
                          keep_snps = F,
                          get_island = NULL,
                          get_chr = NULL,
                          get_type = NULL)

# controls
rg_controls <- subset_rg_set(rg_set = rgControls,
                             keep_gender = F,
                             keep_controls = F,
                             keep_snps = F,
                             get_island = NULL,
                             get_chr = NULL,
                             get_type = NULL)

# valid
rg_valid <- subset_rg_set(rg_set = rgValid,
                          keep_gender = F,
                          keep_controls = F,
                          keep_snps = F,
                          get_island = NULL,
                          get_chr = NULL,
                          get_type = NULL)


rm(rgCases, rgControls, rgValid)


# preprocess controls and valid
beta_cases <- preprocessMethod(rg_cases, preprocess = method)
beta_controls <- preprocessMethod(rg_controls, preprocess = method)
beta_valid <- preprocessMethod(rg_valid, preprocess = method)

rm(rg_cases, rg_controls, rg_valid, ratio_set)

# 
# save.image('~/Desktop/temp_dmr_cancer_raw.RData')
load('~/Desktop/temp_dmr_cancer.RData')


# random selection 


# do cases first (will return list of 2, second element is old controls)
beta_cases <- process_rg_set_single(beta_data = beta_cases[1:200000,], 
                                    id_map = id_map_cases, 
                                    clin = clin)
# get controls
beta_controls_mod <- process_rg_set_single(beta_data = beta_controls[1:200000,],
                                           id_map = id_map_con, 
                                           clin = clin)

# get valid
beta_valid_mod <- process_rg_set_single(beta_data = beta_valid[1:200000,], 
                                        id_map = id_map_val, 
                                        clin = clin)
# 
# if(method == 'raw'){
#   # cases
#   beta_cases <- removeInf(beta_cases, probe_start = 10)
# 
#   # controls 
#   beta_controls_mod <- removeInf(beta_controls_mod, probe_start = 10)
#   
#   # valid 
#   beta_valid_mod <- removeInf(beta_valid_mod, probe_start = 10)
#   
# }

##########
# bumphunter 
##########
intersect_names <- Reduce(intersect, list(colnames(beta_cases)[10:ncol(beta_cases)],
                                          colnames(beta_controls_mod)[10:ncol(beta_controls_mod)],
                                          colnames(beta_valid_mod)[10:ncol(beta_valid_mod)]))

# get clin names
clin_names <- colnames(beta_cases)[1:9]


# train data 
beta_cases <- beta_cases[, c(clin_names,
                             intersect_names)]

# controls data 
beta_controls_mod <- beta_controls_mod[, c(clin_names,
                                           intersect_names)]

# train data 
beta_valid_mod <- beta_valid_mod[, c(clin_names,
                                     intersect_names)]

# get old controls (mut and no cancer from cases)
beta_controls_mod_old <- rbind(subset(beta_cases, p53_germline == 'Mut' & 
                                        cancer_diagnosis_diagnoses == 'Unaffected'))

# get model data in cases for training and test
beta_cases <- getModData(beta_cases)

# full data 
beta_cases <- beta_cases[!duplicated(beta_cases$ids),]

# get rid of cancer samples in controls 
beta_controls_mod <- beta_controls_mod[grepl('Unaffected', beta_controls_mod$cancer_diagnosis_diagnoses),]

#subset valid - get ids from train and test
case_ids <- beta_cases$ids
beta_valid_mod <- beta_valid_mod[!beta_valid_mod$ids %in% case_ids,]

# remove NAs 
beta_cases <-beta_cases[complete.cases(beta_cases),]

##########
# function to plot methylation
##########
beta_thresh = 0.01
cancer_type = 'ERMS'
after_bh = F
diff_meth_vis <- function(beta_thresh, controls, cancer_type, after_bh) {
  
  # determine which controls will be used
  if(controls == 'old') {
    beta_controls_mod <- beta_controls_mod_old
  }
  
  if (controls == 'full') {
    beta_controls_mod <- rbind(beta_controls_mod, beta_controls_mod_old)
  } 
  
  # combine beta_cases and beta_valid_mod
  beta_cases_full <- rbind(beta_cases, beta_valid_mod)
  
  if(remove_cancer_sig) {
    
    bh_feats <- bump_hunter(dat_1 = beta_cases_full, 
                            dat_2 = beta_controls_mod, 
                            bump = 'cancer', 
                            boot_num = 3, 
                            thresh = beta_thresh,
                            g_ranges = g_ranges)
    
    
    # get feature list
    colnames(bh_feats)[1] <- 'chr'
    remove_features <- inner_join(bh_feats, g_ranges)$probe
    
    # subset features by removing these ones
    keep_features <- colnames(beta_cases_full)[!colnames(beta_cases_full) %in% remove_features]
    
    # subset beta_cases_full by removing these features (cancer probes)
    beta_cases_full <- beta_cases_full[, keep_features]
  }
  
  rm(beta_cases,beta_controls, beta_valid)
  # save.image('~/Desktop/temp_dmr_cancer_raw_2.RData')
  load('~/Desktop/temp_dmr_cancer_raw_2.RData')
  cancer_string <- paste('ACC', cancer_type, sep = '|')
  beta_sub <- beta_cases_full[grepl(cancer_string, beta_cases_full$cancer_diagnosis_diagnoses),]
  beta_sub$cancer_diagnosis_diagnoses <- ifelse(grepl(cancer_type, beta_sub$cancer_diagnosis_diagnoses), cancer_type, 'ACC')
  
  # scale data seperately
  beta_sub_acc <- subset(beta_sub, cancer_diagnosis_diagnoses == 'ACC')
  beta_sub_acc[, 10:ncol(beta_sub_acc)] <- scale(beta_sub_acc[, 10:ncol(beta_sub_acc)])
  beta_sub_type <- subset(beta_sub, cancer_diagnosis_diagnoses == cancer_type)
  beta_sub_type[, 10:ncol(beta_sub_type)] <- scale(beta_sub_type[, 10:ncol(beta_sub_type)])
  
  # combine them
  beta_sub_all <- rbind(beta_sub_acc, 
                        beta_sub_type)
  
  # remove 
  rm(beta_sub_acc, beta_sub_type)
  
  # save.image('~/Desktop/temp_dmr_cancer_raw_3.RData')
  load('~/Desktop/temp_dmr_cancer_raw_3.RData')
  
  beta_sub_all <- subset(beta_sub_all, ids != '3301')
  beta_sub_all <- subset(beta_sub_all, ids != '3304')
  beta_sub_all <- subset(beta_sub_all, ids != '3389')
  beta_sub_all <- subset(beta_sub_all, ids != '4329')

  beta_sub_all <- beta_sub_all[, c(1:200, 399, 677, 999, 1221, 2000:3000)]
  
  # make a vector for cancer 
  cancer_type <- beta_sub_all$cancer_diagnosis_diagnoses
  
  # remove unneeded columns 
  beta_sub_all$ids <- beta_sub_all$p53_germline  <- beta_sub_all$age_diagnosis <- 
    beta_sub_all$age_sample_collection<- beta_sub_all$gender<- beta_sub_all$sentrix_id <-
    beta_sub_all$family_name <- beta_sub_all$tm_donor_ <- beta_sub_all$cancer_diagnosis_diagnoses <- NULL
  
  
  beta_mat <- as.matrix(beta_sub_all)
  
  # save.image('~/Desktop/temp_dmr_cancer_raw_group.RData')
  # load('~/Desktop/temp_dmr_cancer_raw_group.RData')
  cancer_type <- ifelse(cancer_type == 'ERMS', 'OS', 'ACC')
  if(density_plot) {
    
    densityPlot(beta_mat, sampGroups = cancer_type, main = "ACC vs OS", xlab = "Beta Value", pal = brewer.pal(8, "Dark2"),
                xlim = c(0,1.0), ylim = c(0,20))
    
  }
  
  # sample 
  col_index <- sample(ncol(beta_mat), 1000, replace = T)
  temp <- as.matrix(beta_mat[, col_index])
  
  if(heat_map) {
    
    d_temp <- dist(beta_mat, method = 'euclidean') # method="man" # is a bit better
    hc_temp <- hclust(d_temp, method = "complete")
    hc_cancer <- cancer_type
    
    dend <- as.dendrogram(hc_temp)
    # order it the closest we can to the order of the observations:
    dend <- rotate(dend, 1:42)
    
    # Color the branches based on the clusters:
    dend <- color_branches(dend, k=2) #, groupLabels=iris_species)
    
    # Manually match the labels, as much as possible, to the real classification of the flowers:
    labels_colors(dend) <-
      rainbow_hcl(2)[sort_levels_values(
        as.numeric(hc_cancer)[order.dendrogram(dend)]
      )]
    
    # We shall add the flower type to the labels:
    labels(dend) <- paste(as.character(cancer_type)[order.dendrogram(dend)],
                          "(",labels(dend),")", 
                          sep = "")
    # We hang the dendrogram a bit:
    dend <- hang.dendrogram(dend,hang_height=0.1)
    # reduce the size of the labels:
    # dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
    dend <- set(dend, "labels_cex", 0.5)
    # And plot:
    par(mar = c(3,3,3,7))
    plot(dend, 
         main = "Clustered methylation
         (the labels give the true cancer type)", 
         horiz =  TRUE,  nodePar = list(cex = .007))
    legend("topleft", legend = unique(cancer_type), fill = rainbow_hcl(3))
    
   
    scaled_temp <- temp %>% as.matrix %>% scale
    # library(gplots)
    gplots::heatmap.2(as.matrix(beta_mat),  
                      main = "Heatmap for LFS",
                      srtCol = 20,
                      dendrogram = "row",
                      Rowv = dend,
                      Colv = "NA", # this to make sure the columns are not ordered
                      trace="none",          
                      denscol = "grey",
                      key = F,
                      density.info = "density") # to add nice colored strips        
    
    
  }
  
  
}


