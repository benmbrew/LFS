
# source functions script
source('all_functions.R')

# object for beta or m value analysis
methyl_type <- 'm'
combat <- TRUE

# read in all data

if(methyl_type == 'beta'){
  if(combat){
    all_cases_beta_combat <- readRDS('../../Data/all_cases_beta_combat.rda')
    all_con_beta_combat <- readRDS('../../Data/all_con_beta_combat.rda')
  } else {
    all_cases_beta <- readRDS('../../Data/all_cases_beta.rda')
    all_con_beta <- readRDS('../../Data/all_con_beta.rda')
  }
} else {
  if(combat){
    all_cases_m_combat <- readRDS('../../Data/all_cases_m_combat.rda')
    all_con_m_combat <- readRDS('../../Data/all_con_m_combat.rda')
  } else {
    all_cases_m <- readRDS('../../Data/all_cases_m.rda')
    all_con_m <- readRDS('../../Data/all_con_m.rda')
  }
}


# prepare data sets for modelling
