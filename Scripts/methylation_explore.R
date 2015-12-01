################################################################
# Script for looking at new methylation data from Nardin and comparing it to methylation data from 
# project with Daniel

# Load in methylation data for patients with ACC (we have 29 WT and 36 Mut), CPC (12 WT, 12 Mut)
# Osteo (4 WT, 4 Mut), Medulloblastoma (? WT, 4 Mut), RMS (1 WT, 13 Mut)
# This is 5 cancers and 115 samples .

library(dplyr)

# set data working directory 
data <- '/home/benbrew/Documents/li_fraumeni/Data'
setwd(data)

############################################################## 
# Read in data
methyl <- read.table('methyl.txt')
clin <- read.csv('clin_all.csv', header = FALSE)
#methyl_short <- read.csv('methyl_17.csv', header = FALSE)

##############################################################
# Subset clinical data by IDs from methylation 
id_methyl <- as.character(methyl[1,])
#id_short_methyl <- as.character(methyl_short[2, ])
#short_methyl_clin <- clin[clin$V1 %in% id_short_methyl,]

# Find IDs separated by "/"

## A helper function that tests whether an object is either NULL _or_ 
## a list of NULLs
is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))

## Recursively step down into list, removing all such objects 
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x) # takes all elements that are not null
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}

# How much of the ids from the methylation data intersect 
# with ids from the clinical data
length(which(id_methyl %in% clin$V1)) # 22 

# but some of clinical data IDs are concatenated with "/". 
# check to see if any of the ids from methyl are in the concatenated ids
# from clinical data

# Make a list to store the concatenated clinical ids that have a methylation id within it
store_list <- list()

for(i in 1:nrow(clin)){   
  print(i)
  ids <- as.character(clin[,1 ])
  id_temp <- ids[i]  

  if(grepl("/", id_temp)){
    sep <- unlist(strsplit(id_temp, "/"))
   
    if(any(sep %in% id_methyl)){
      store_list[[i]] <- id_temp
      #sep <- id_methyl[sep == id_methyl]
      #store_list[[i]] <- sep
      store_list <- rmNullObs(store_list)
      print("id found here")
  }
} 
    
  sep_none <- id_temp
  if(any(sep_none %in% id_methyl)){
    print("already found id")
  
    } else {
    print("no extra ids found")
  }  
}

# add the concatenated ids from clinical to list of methylation ids and subset clinical data by this index
id_ind <- append(unlist(store_list), id_methyl)
methyl_clin <- clin[clin$V1 %in% id_ind,]



###########################################################
# Clean methyl_clin and explore
colnames(methyl_clin) <- c("id", "tp53", "cancer", "cancer_indicator", "age_of_onset",
                           "gdna", "protein", "codon_72", "pin_3", "mdm2","date", "gender")

# What cancers are represented
cancer <- methyl_clin %>% 
  group_by(cancer, tp53) %>% 
  summarise(count = n())

#cancer_short <- short_methyl_clin %>%
  #group_by(id, cancer, tp53) %>%
  #summarise(count = n())

# see if overlapping with 17 methylation patients
# send nardin an email 

