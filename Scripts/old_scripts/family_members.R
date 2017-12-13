## this script will read in batch data from get_cases, get_controls, or get_valid and explore potential batches and outliers

##########
# initialize libraries
##########
library(tidyverse)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
clin_data <- paste0(data_folder, '/clin_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'noob'
k = 10
combat = F

##########
# load data
##########


if (combat) {
  
  betaFull <-  readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_combat_m.rda')))
  
  
} else {
  betaFull <- readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_m.rda')))
  
}


##########
# get column names
##########
intersect_names <- colnames(betaFull)[9:ncol(betaFull)]

##########
# read in all features from feat data
##########

lfs_feats_m <- readRDS(paste0(feat_data, paste0('/', method, '_', 'lfs_m.rda')))
no_cancer_feats_m <- readRDS(paste0(feat_data, paste0('/', method, '_', 'no_cancer_m.rda')))


# setwd(feat_data)
# file_list_names = list.files()
# 
# 
# # store all raw rda feature lists in feat_list
# feat_list <- lapply(file_list_names, function(x) readRDS(x))
# 
# # order feat list
# feat_list  <- feat_list[order(sapply(feat_list, length), decreasing=F)]
# 
# # select first 10 
# feat_list <- feat_list[5]
# 
# # model_names <- c('enet', 'rf', 'lasso')
# # seeds <- c(1, 2, 3)

model_names <- c('enet')
seeds <- c(1,2,3)
feat_list <- list(no_cancer_feats_m, lfs_feats_m)
file_list_names <- list('no_cancer_m', 'lfs_m')

##########
# read in clinical data
##########
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# clean clinical ids
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

##########
# cases family members
##########
betaCases <- getModData(betaCases)

cases_ids <- as.data.frame(betaCases$ids)
names(cases_ids) <- 'ids'

length(which(cases_ids$ids %in% clin$ids))

temp <- inner_join(cases_ids, clin, by = 'ids')

# remove NA from family 
temp <- temp[!is.na(temp$family_name),]

temp <- temp[!duplicated(temp$tm_donor_),]

summary(as.factor(temp$family_name))
length(unique(temp$family_name))


##########
# controls family member
##########

controls_ids <- as.data.frame(betaControls$ids)
names(controls_ids) <- 'ids'

length(which(controls_ids$ids %in% clin$ids))

temp_c <- inner_join(controls_ids, clin, by = 'ids')

# remove NA from family 
temp_c<- temp_c[!is.na(temp_c$family_name),]

temp_c <- temp_c[!duplicated(temp_c$tm_donor_),]

summary(as.factor(temp_c$family_name))
length(unique(temp_c$family_name))


##########
# valid family member
##########

valid_ids <- as.data.frame(betaValid$ids)
names(valid_ids) <- 'ids'

length(which(valid_ids$ids %in% clin$ids))

temp_v <- inner_join(valid_ids, clin, by = 'ids')

# remove NA from family 
temp_v<- temp_v[!is.na(temp_v$family_name),]

temp_v <- temp_v[!duplicated(temp_v$tm_donor_),]

summary(as.factor(temp_v$family_name))
length(unique(temp_v$family_name))


##########
# overlapping family members
##########

# controls and cases
length(which(temp$family_name %in% temp_c$family_name))
length(which(temp$family_name %in% temp_v$family_name))

fml = as.formula(paste("age", paste0("x", 1:100, collapse=" + "), sep=" ~ "))

mod_dat <- betaCases[, 8:ncol(betaCases)]
mod <- cbind(age = betaCases$age_diagnosis, 
             family = temp$family_name,
             mod_dat)

mod$family <- NULL

mod_temp <- mod[, 1:16660]

fml = as.formula(paste("age", paste0(colnames(mod_temp[, -1]), collapse=" + "), sep=" ~ "))

mod_sub <- mod[, 1:1000]

lm(mod_sub, data = mod)
glmmLasso(mod_sub, rnd = NULL, lambda = 10, data = mod)
lm(mod)



## Not run:
data(knee)
knee[,c(2,4:6)]<-scale(knee[,c(2,4:6)],center=TRUE,scale=TRUE)
knee<-data.frame(knee)

## fit adjacent category model
glm.obj <- glmmLasso(pain ~ time + th + age + sex, rnd = NULL,
                     family = cumulative(), data = knee, lambda=10,
                     switch.NR=TRUE, control=list(print.iter=TRUE))
summary(glm.obj)

## generalized linear mixed model with a
library(glmmLasso)

## Not run:
data("soccer")
soccer[,c(4,5,9:16)]<-scale(soccer[,c(4,5,9:16)],center=TRUE,scale=TRUE)
soccer<-data.frame(soccer)

## linear mixed model
lm1 <- glmmLasso(fmla ,  rnd = list(team=~1),
                 lambda=10, data = soccer, family=binomial(link = "logit"))


soccer$points <-as.numeric(ifelse(soccer$points > 50, 1, 2))


lm1 <- lm(points ~ . - team, data = soccer)


summary(lm1)

fmla <- as.formula(points ~ transfer.spendings + ave.unfair.score
                   + ball.possession + tackles)

## similar linear model without random effects
lm1b <- glmmLasso(fmla, rnd = NULL,
                  lambda=10, data = soccer)
summary(lm1b)
## linear mixed model with slope on ave.attend;
## the coefficient of ave.attend is not penalized;
lm2 <- glmmLasso(points~transfer.spendings + ave.unfair.score
                 + ball.possession + tackles + ave.attend
                 + sold.out, rnd = list(team=~1 + ave.attend), lambda=10,
                 data = soccer, control = list(index=c(1,2,3,4,NA,5),
                                               method="REML",print.iter=TRUE))
summary(lm2)
