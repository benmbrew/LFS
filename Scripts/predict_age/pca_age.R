# this script will analyze the PCs that are associated with age and remove it.
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
# https://cran.r-project.org/web/packages/superpc/superpc.pdf
# get functions
source('all_functions.R')

# cases
all_cases_beta <- readRDS('../../Data/all_cases_beta.rda')
all_cases_beta_combat <- readRDS('../../Data/all_cases_beta_combat.rda')

# controls
all_con_beta <- readRDS('../../Data/all_con_beta.rda')
all_con_beta_combat <- readRDS('../../Data/all_con_beta_combat.rda')

# read m data
# cases
all_cases_m <- readRDS('../../Data/all_cases_m.rda')
all_cases_m_combat <- readRDS('../../Data/all_cases_m_combat.rda')

# controls
all_con_m <- readRDS('../../Data/all_con_m.rda')
all_con_m_combat <- readRDS('../../Data/all_con_m_combat.rda')

# controls wt
all_con_beta_wt <- readRDS('../../Data/all_con_beta_wt.rda')
all_con_beta_wt_combat <- readRDS('../../Data/all_con_beta_wt_combat.rda')

all_con_m_wt <- readRDS('../../Data/all_con_m_wt.rda')
all_con_m_wt_combat <- readRDS('../../Data/all_con_m_wt_combat.rda')
