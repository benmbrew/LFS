
source('helper_functions.R')

##########
# strategy 1
##########
temp_con_norm <- read.csv('~/Desktop/strategy_1/con_table_normal.csv')
temp_val_norm <- read.csv('~/Desktop/strategy_1/valid_table_normal.csv')

temp_con_log <- read.csv('~/Desktop/strategy_1/con_table_log.csv')
temp_val_log <- read.csv('~/Desktop/strategy_1/valid_table_log.csv')

# con normal
plot_pred(temp_con_norm, 
          type = 'con', 
          plot_type = 'age_pred',
          strategy = 'strategy_1', 
          log = FALSE)
plot_pred(temp_con_norm, 
          type = 'con', 
          plot_type = 'ROC',
          strategy = 'strategy_1', 
          log = FALSE)
plot_pred(temp_con_norm, 
          type = 'con', 
          plot_type = 'confusion', 
          strategy = 'strategy_1',
          log = FALSE)

# non log
plot_pred(temp_con_log, 
          type = 'con', 
          plot_type = 'age_pred',
          strategy = 'strategy_1', 
          log = TRUE)
plot_pred(temp_con_log, 
          type = 'con', 
          plot_type = 'ROC',
          strategy = 'strategy_1', 
          log = TRUE)
plot_pred(temp_con_log, 
          type = 'con', 
          plot_type = 'confusion',
          strategy = 'strategy_1', 
          log = TRUE)

# VALID
# val normal
plot_pred(temp_val_norm, 
          type = 'val', 
          plot_type = 'age_pred',
          strategy = 'strategy_1', 
          log = FALSE)
plot_pred(temp_val_norm, 
          type = 'val', 
          plot_type = 'ROC',
          strategy = 'strategy_1', 
          log = FALSE)
plot_pred(temp_val_norm, 
          type = 'val', 
          plot_type = 'valfusion', 
          strategy = 'strategy_1',
          log = FALSE)

# non log
plot_pred(temp_val_log, 
          type = 'val', 
          plot_type = 'age_pred',
          strategy = 'strategy_1', 
          log = TRUE)
plot_pred(temp_val_log, 
          type = 'val', 
          plot_type = 'ROC',
          strategy = 'strategy_1', 
          log = TRUE)
plot_pred(temp_val_log, 
          type = 'val', 
          plot_type = 'valfusion',
          strategy = 'strategy_1', 
          log = TRUE)

############################################################################
##########
# strategy 2
##########
temp_con <- read.csv('~/Desktop/strategy_2/con_table.csv')
temp_val <- read.csv('~/Desktop/strategy_2/valid_table.csv')
# con normal

plot_pred(temp_con, 
          type = 'con', 
          plot_type = 'age_pred',
          strategy = 'strategy_2', 
          log = FALSE)
plot_pred(temp_con, 
          type = 'con', 
          plot_type = 'ROC',
          strategy = 'strategy_2', 
          log = FALSE)
plot_pred(temp_con, 
          type = 'con', 
          plot_type = 'confusion', 
          strategy = 'strategy_2',
          log = FALSE)
# VALID
# val normal
plot_pred(temp_val, 
          type = 'val', 
          plot_type = 'age_pred',
          strategy = 'strategy_2', 
          log = FALSE)
plot_pred(temp_val, 
          type = 'val', 
          plot_type = 'ROC',
          strategy = 'strategy_2', 
          log = FALSE)
plot_pred(temp_val, 
          type = 'val', 
          plot_type = 'confusion', 
          strategy = 'strategy_2',
          log = FALSE)


############################################################################
##########
# strategy 4
##########
temp_con <- read.csv('~/Desktop/strategy_3/null_table_pc.csv')
temp_val <- read.csv('~/Desktop/strategy_3/valid_table_pc_new.csv')

# get controls_age_label
temp_con$controls_age_label <- ifelse(temp_con$age_sample_collection < 72, 'positive', 'negative')
# get valid_age_label
temp_val$valid_age_label <- ifelse(temp_val$age_sample_collection < 72, 'positive', 'negative')

# con normal
plot_pred(temp_con, 
          type = 'con', 
          plot_type = 'age_pred',
          strategy = 'strategy_3', 
          log = FALSE)
plot_pred(temp_con, 
          type = 'con', 
          plot_type = 'ROC',
          strategy = 'strategy_3', 
          log = FALSE)
plot_pred(temp_con, 
          type = 'con', 
          plot_type = 'confusion', 
          strategy = 'strategy_3',
          log = FALSE)
# VALID
# val normal
plot_pred(temp_val, 
          type = 'val', 
          plot_type = 'age_pred',
          strategy = 'strategy_3', 
          log = FALSE)
plot_pred(temp_val, 
          type = 'val', 
          plot_type = 'ROC',
          strategy = 'strategy_3', 
          log = FALSE)
plot_pred(temp_val, 
          type = 'val', 
          plot_type = 'confusion', 
          strategy = 'strategy_3',
          log = FALSE)

temp <- readRDS('data_cv/lfs_bump_probesbeta_no_combat.rda')
temp <- sample(temp, 26)
temp <- as.data.frame(temp)
names(temp) <- 'probe'

temp1 <- readRDS('all_data/all_con_beta.rda')
temp1 <- temp1[temp1$tech == '850k',]

# get intersection
temp_clin <- temp1[, 1:12]
temp_int <- intersect(names(temp1), temp$probe)
temp1 <- temp1[, temp_int]
temp1 <- as.data.frame(cbind(temp_clin, temp1))

# read in all data
temp_con <- read.csv('~/Desktop/read_in.csv')
con_trans <- readRDS(paste0('transform_data_cv/', 'con_transform_','beta', '.rda'))
# get intersection
temp_clin <- con_trans[, 1:12]
temp_int <- intersect(names(con_trans), temp$probe)
con_trans <- con_trans[, temp_int]
con_trans <- as.data.frame(cbind(temp_clin, con_trans))

con_norm <- inner_join(temp1, temp_con, by = 'tm_donor')

tms <- '4749|4520|2942'

con_norm <- con_norm[grepl(tms, con_norm$tm_donor),]
con_trans <- con_trans[grepl(tms, con_trans$tm_donor),]

cpgs <- names(con_trans)[13:ncol(con_trans)]
con_norm$type <- c('false_positive', 'true_positive', 'false_negative')
con_trans$type <- c('false_positive', 'true_positive', 'false_negative')


pdf('~/Desktop/cpgs_plots')
for(i in 1:length(cpgs)){
  i = 3
  this_cpg <- cpgs[i]
  temp_norm <- con_norm[, c('age_sample_collection', 'tm_donor', 'type', this_cpg)]
  temp_trans <- con_trans[, c('age_sample_collection', 'tm_donor', 'type',this_cpg)]
  
  names(temp_norm) <- c('V1', 'V2', 'V3', 'V4')
  names(temp_trans) <- c('V1', 'V2', 'V3', 'V4')
  
  temp_trans$V1 <- temp_norm$V1
  
  ggplot(temp_norm, aes(V1, V4)) + geom_point() +
    geom_label(aes(label = V3)) +
    xlim(c(-5, 70))+
  labs(title = this_cpg,
       x = 'age',
         y = 'cpgs')
  
  ggplot(temp_trans, aes(V1, V4)) + geom_point() +
    geom_label(aes(label = V3)) +
    xlim(c(-5, 70))+
    labs(title = this_cpg,
         x = 'age',
         y = 'cpgs')
}

plot(all_dat$age_sample_collection, all_dat[, this_cpg])


dev.off()
plot(all_dat$age_sample_collection.x[1], all_dat[, 13])










# join temp_con and temp
temp_normal <- inner_join(temp_con, temp, by ='tm_donor')
# join temp_con and con_trans

temp_clin <- temp[, 1:12]
temp_int <- intersect(names(temp), probes)
temp <- temp[, temp_int]
temp <- as.data.frame(cbind(temp_clin, temp))
temp <- read.csv(all_26,'~/Desktop/all_26.csv')

# here
all_dat <- inner_join(temp_con, temp, by = 'tm_donor')
