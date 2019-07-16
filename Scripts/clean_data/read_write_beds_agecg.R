
library(sqldf)
# read in .bed file 
bed <- as.data.frame(read.table("../../Data/age_cgs.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
bed <- bed[-1, ]
# write_csv(bed, '../../Data/age_cgs.csv')

# load genomic methyl set (from controls) - you need genetic locations by probe from this object
ratio_set <- readRDS('../../Data/model_data/raw_ratio_set.rda')

# get granges object
g_ranges <- as.data.frame(getLocations(ratio_set))

# get probes from rownames
g_ranges$probe <- rownames(g_ranges)

# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]

#

# use sql data table to merger validators with model_data based on age range
result = sqldf("select * from g_ranges
                 inner join bed
                 on start between bed.V2 and bed.V3")


write_rds(result$probe, '../../Data/age_probes.rda')
