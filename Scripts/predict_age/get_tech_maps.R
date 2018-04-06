

# this script will visualize batch effects with PCA 

library(RColorBrewer)
# set preprocessing method
method <- 'noob'

# old or new data
data_used <- 'new'

# set type of data, beta or m
methyl_type <- 'beta'

# source all_functions.R to load libraries and my functions
source('all_functions.R')

# set data directory
data_dir <- '../../Data/'

# load data 
all_data <- readRDS(paste0(data_dir, 'all_data.rda'))


# subset to remove outlier 
all_data <- all_data[all_data$tm_donor_ != '3955',]


# get patients that have mulitple samples (no dup ids, but dup tm_donor_)
all_data <- all_data[!duplicated(all_data$ids),]
length(which(duplicated(all_data$tm_donor_)))
dup_tms <- all_data$tm_donor_[duplicated(all_data$tm_donor_)]
temp <- all_data[all_data$tm_donor_ %in% dup_tms,]
temp_short <- temp[, 1:40]

# get a cases (450) and valid (850) data set 
# get controls (450) and controls (850)

cases_450 <- all_data[!grepl('Unaffected', all_data$cancer_diagnosis_diagnoses) & grepl('^57|97', all_data$sentrix_id) & all_data$p53_germline == 'Mut',]
controls_450 <- all_data[grepl('Unaffected', all_data$cancer_diagnosis_diagnoses) & grepl('^57|97', all_data$sentrix_id) & all_data$p53_germline == 'Mut',]
cases_850 <- all_data[!grepl('Unaffected', all_data$cancer_diagnosis_diagnoses) & !grepl('^57|97', all_data$sentrix_id) & all_data$p53_germline == 'Mut',]
controls_850 <- all_data[grepl('Unaffected', all_data$cancer_diagnosis_diagnoses) & !grepl('^57|97', all_data$sentrix_id) & all_data$p53_germline == 'Mut',]

# remove duplicates
cases_450 <- cases_450[!duplicated(cases_450$tm_donor_),]
controls_450 <- controls_450[!duplicated(controls_450$tm_donor_),]
cases_850 <- cases_850[!duplicated(cases_850$tm_donor_),]
controls_850 <- controls_850[!duplicated(controls_850$tm_donor_),]


# find ids that are in both cases and controls qu
shared_ids_cases <- unique(cases_450$ids)[unique(cases_450$ids) %in% unique(cases_850$ids)]
shared_ids_controls <- unique(controls_850$ids)[unique(controls_850$ids) %in% unique(controls_450$ids)]

# subset cases 
cases_450_sub <-  cases_450[cases_450$ids %in% shared_ids_cases,]
cases_850_sub <-  cases_850[cases_850$ids %in% shared_ids_cases,]

# controls
controls_450_sub <-  controls_450[controls_450$ids %in% shared_ids_controls,]
controls_850_sub <-  controls_850[controls_850$ids %in% shared_ids_controls,]


##########
# function to plot each id against the other
##########


smoothScatter(x, y = NULL, nbin = 128, bandwidth,
              colramp = colorRampPalette(c("white", blues9)),
              nrpoints = 100, ret.selection = FALSE,
              pch = ".", cex = 1, col = "black",
              transformation = function(x) x^.25,
              postPlotHook = box,
              xlab = NULL, ylab = NULL, xlim, ylim,
              xaxs = par("xaxs"), yaxs = par("yaxs"), ...)


temp_dat <- temp
tm_id <- '2941'

# get_corr <- function(temp_dat, tm_id){
#   sample_tm <- temp[temp$tm_donor_ == tm_id, ]
#   cor_score <- cor(sample_tm[1, 12:ncol(sample_tm)], 
#                    sample_tm[2, 12:ncol(sample_tm)])
#   
# }

plotCaseCon <- 
  function (temp_dat, tm_id) {
    
    sample_tm <- temp_dat[temp_dat$tm_donor_ %in% tm_id,]

    smoothScatter(sample_tm[1,12:ncol(sample_tm) ], 
                  sample_tm[2, 12:ncol(sample_tm)], 
                  main = tm_id,
                  cex= 0.5,
                  cex.axis = 1.5,
                  cex.lab = 1.3,
                  xlab = 'First Sample', 
                  ylab = 'Second Sample',
                  xlim = c(0, 1),
                  ylim = c(0, 1),
                  bty = 'n')
    
  }

# 2941 2942 2921 3872 1169 2828 2188 4613 2973 3704 3714 3615 2850 2942 2188 3714 2941 2921 3872 1169 2828 4613 2973 3704 3615 2850
plotCaseCon(temp, '2941')

str(temp)
# plot data- cases
plotCaseCon(cases_450_sub, cases_850_sub, row_index = 1)
plotCaseCon(cases_450_sub, cases_850_sub, row_index = 2)
plotCaseCon(cases_450_sub, cases_850_sub, row_index = 3)
plotCaseCon(cases_450_sub, cases_850_sub, row_index = 4)
plotCaseCon(cases_450_sub, cases_850_sub, row_index = 5)
plotCaseCon(cases_450_sub, cases_850_sub, row_index = 6)
plotCaseCon(cases_450_sub, cases_850_sub, row_index = 7)
plotCaseCon(cases_450_sub, cases_850_sub, row_index = 8)
plotCaseCon(cases_450_sub, cases_850_sub, row_index = 9)


# plot data- controls
plotCaseCon(controls_450_sub, controls_850_sub, row_index = 1)
plotCaseCon(controls_450_sub, controls_850_sub, row_index = 2)
plotCaseCon(controls_450_sub, controls_850_sub, row_index = 3)
plotCaseCon(controls_450_sub, controls_850_sub, row_index = 4)
plotCaseCon(controls_450_sub, controls_850_sub, row_index = 5)
plotCaseCon(controls_450_sub, controls_850_sub, row_index = 6)
plotCaseCon(controls_450_sub, controls_850_sub, row_index = 7)
plotCaseCon(controls_450_sub, controls_850_sub, row_index = 1)
plotCaseCon(controls_450_sub, controls_850_sub, row_index = 1)

