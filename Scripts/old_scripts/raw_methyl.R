###################################################
# this script is for reading in raw methylation files. 
setwd('/home/benbrew/Documents/LFS/Data/')
library(dplyr)


if("raw_methylation.RData" %in% dir()){
  load("raw_methylation.RData")
}else{
# Load libraries
# Initialize folders
home_folder <- '/home/benbrew/Documents'
project_folder <- paste(home_folder, 'LFS', sep = '/')
data_folder <- paste(project_folder, 'Data', sep = '/')
methyl_data <- paste(data_folder, 'methyl_files/', sep = '/')
data_name <- 'Chr'

# Read in Nardine's methylation data and clinical data 
methyl <- read.table(paste(data_folder,'methyl.txt', sep = '/'))
clin <- read.csv(paste(data_folder, 'clin_all.csv', sep = '/'), header = FALSE)
methyl_17 <- read.csv(paste(data_folder, 'methyl_17.csv', sep = '/'), header = TRUE)

num_sets <- 23 
dat <-  vector('list', num_sets)

# read in raw methylation files. 
for(i in (1:num_sets)){
  dat[[i]] <- read.delim(paste(methyl_data, data_name, i, '.txt', sep = ''))
}

# save Rdata file 
setwd('/home/benbrew/Documents/LFS/Data')
save.image('raw_methylation.RData')
}
### examine elements of the dat list that have ch.i as a column 

### Clean dates in clinical data, column V11

####################################################
# Make a list to store the concatenated clinical ids that have a methylation id within it
# subset clinical data by methylation data from nardin.

if("clinical.RData" %in% dir()){
  load("clinical.RData")
}else{
# helper functions
is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x) # takes all elements that are not null
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)}
#####
clinMatch <- function(clin, methylation_id){
  
  store_list <- list()

  for(i in 1:nrow(clin)){   
    print(i)
    ids <- as.character(clin[,1 ])
    id_temp <- ids[i]  
    
    if(grepl("/", id_temp)){
      sep <- unlist(strsplit(id_temp, "/"))
      
      if(any(sep %in% methylation_id)){
        store_list[[i]] <- id_temp
        store_list <- rmNullObs(store_list)
        print("id found here")
      }
    } 
    
    sep_none <- id_temp
    if(any(sep_none %in% methylation_id)){
      print("already found id")
      
    } else {
      print("no extra ids found")
    }  
  }
  id_ind <- append(unlist(store_list), methylation_id)
  methyl_clin <- clin[clin$V1 %in% id_ind,]
  return(methyl_clin)
  
}


# apply function to match methylation ids from nardin to clinical data
id_methyl <- as.character(methyl[1,])
clin_methyl <- clinMatch(clin, id_methyl)
colnames(clin_methyl) <- c("id", "tp53", "cancer", "cancer_indicator", "age_of_onset",
                        "gdna", "protein", "codon_72", "pin_3", "mdm2","date", "gender")

## 24 of the methylation data Nardin originially sent are in clin


######################################################################
# subset clinical data by raw methylation data ids

# check consistency in list
num_col <- list()
num_row <- list()
column_names <- list()
first_col <- list()

for (i in 1:num_sets){
  num_col[[i]] <- ncol(dat[[i]])
  num_row[[i]] <- nrow(dat[[i]])
  column_names[[i]] <- colnames(dat[[i]])
  first_col[[i]] <- dat[[i]][, 1]
}

# # check if column names and row names have overlap
# full_result <- list()
# for(col1 in 1:num_sets){
#   column1 <- column_names[[col1]]
#   result <- rep.int(0, num_sets)
#   for(col2 in (1:num_sets)[-col1]){
#     column2 <- column_names[[col2]]
#   if(any(column1 %in% column2)){
#     result[col2] <- TRUE
#     }else{
#       result[col2] <- FALSE
#     }
#   }
#   full_result <- rbind(full_result, result)
# }

# the only column name in common is the "X" column indicating the patient IDs. 

# check if first column (IDs) have overlap. there are exactly 139 rows in each set. 
# the 13th set is messed up, so just focuse on the others for now.
####################################################################################
# You want to extract the numbers that follow the last # or "RD-" sign in the first column of every set except
# 13. 
id_raw <- vector('list', num_sets)


for(i in (1:num_sets)[-13]){
  sub_names <- first_col[[i]]
  column_split <- strsplit(as.character(sub_names), '#')
  last_digits <- lapply(column_split, function(x) x[length(x)])
  sub_ids <- unlist(last_digits)
  id_raw[[i]] <- gsub('RD-', '', sub_ids)
}

clin_raw <- clinMatch(clin, id_raw[[1]])

# the ids across all sets are identical. See which one of these matches with malkin ids
colnames(clin_raw) <- c("id", "tp53", "cancer", "cancer_indicator", "age_of_onset",
                           "gdna", "protein", "codon_72", "pin_3", "mdm2","date", "gender")

################################################################################################
# ids raw methylation data with ids from nardine's methylation data by comparing id_methyl and 
# id_raw[[1]]. Because there are repeating ids in id_raw[[1]], check to see how many id_methyl 
# are in id_raw[[1]]

length(which(id_methyl %in% id_raw[[1]]))
# only 18 of Nardine's id_methylation data are in the raw methylation data. 
# What are those IDs that are not in the raw data 
id_overlap <- id_methyl[id_methyl %in% id_raw[[1]]] # these are the 18 in both
id_extra <- id_methyl[!id_methyl %in% id_raw[[1]]] # these are the 7 in Nardine's original 
# methylation data that are not in the raw data she sent. 

# combine clin_methyl and clin_raw using by rbinding the observations not in clin_raw to clin_raw
clin_full <- rbind(clin_raw, clin_methyl[!clin_methyl$id %in% clin_raw$id,])

# Look at 17 patients we already have. 
id_17 <- methyl_17[1,]
id_17 <- as.character(id_17[which(!is.na(id_17))])

# apply clinMatch to 17 methylation data
clin_17 <- clinMatch(clin, id_17)


# the ids across all sets are identical. See which one of these matches with malkin ids
colnames(clin_17) <- c("id", "tp53", "cancer", "cancer_indicator", "age_of_onset",
                        "gdna", "protein", "codon_72", "pin_3", "mdm2","date", "gender")

# how many of clin_17 are in clin full
id_full_overlap <- id_17[id_17 %in% id_raw[[1]]] # 9 of id_17 are in id_raw
id_full_extra <- id_17[!id_17 %in% id_raw[[1]]] # 9 of 1d_17 are not in id_raw
# 18 total because 'patient ID' is first character of the id row. 

id_full_overlap <- id_17[clin_17$id %in% clin_full$id] # 12 of id_17 are in clin_full
id_full_extra <- id_17[!clin_17$id %in% clin_full$id] # 6 of 1d_17 are not in clin_full

# combine clin_17 and clin_full using by rbinding the observations not in clin_raw to clin_raw
clin_full <- rbind(clin_full, clin_17[!clin_17$id %in% clin_full$id,]) # 6 were added

setwd('/home/benbrew/Documents/LFS/Data')
save.image('clinical.RData')
}

# now group by cancer and summarise to get stats 
cancers <- clin_full %>%
  group_by(as.factor(cancer), tp53) %>%
  summarise(n = n())

# Avegage age for each cancer and tp53 status 
age <- clin_full %>%
  group_by(as.factor(cancer), tp53) %>%
  summarise(mean_age = mean(age_of_onset, na.rm = T))

# Cancer by gender 
gender <- clin_full %>% 
  group_by(cancer) %>%
  summarise(female = sum(gender == 1, na.rm = T),
            male = sum(gender == 0, na.rm = T))
            

###########################################################################################
# # read in full clinical data 
# full_clin <- read.csv('full_clin.csv')
# 
# # get ids from full clin and compare to methylation data
# full_clin_ids <- as.character(full_clin[,2])
# 
# # first check if clin ids are in full clin. 
# clin_ids <- as.character(clin[,1])
# clin_ids %in% full_clin_ids
# 
# # not all clin_ids are in full_clin_ids which is odd. Examine those ones closer 
# clin_none <- clin[!clin_ids %in% full_clin_ids,]
# full_clin_none <- full_clin[!full_clin_ids %in% clin_ids,]
# 
# # Use clinMatch to match methylation data ids with full_clin ids 
# # first remove first column of full_clin, so that the new first column is the id, 
# # as the function needs. replace column name with V1 as well. 
# full_clin <- full_clin[,-1]
# colnames(full_clin)[1] <- 'V1'
# 
# full_clin_methyl <- clinMatch(full_clin, id_methyl)
# full_clin_raw <- clinMatch(full_clin, id_raw[[1]])
# full_clin_17 <- clinMatch(full_clin, id_17)
# 
# # combine these three data sets without creating duplicates. 
# full_clin_total <- rbind(full_clin_raw, full_clin_17[!full_clin_17$V1 %in% full_clin_raw$V1,])
# full_clin_total <- rbind(full_clin_total, full_clin_methyl[!full_clin_methyl$V1 %in% full_clin_raw$V1,])
# 
# # full_clin_total has 57 IDs. How many of these are in clin_full
# clin_match <- full_clin_total[!full_clin_total$V1 %in% clin_full$id,]
