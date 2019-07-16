

## get functions
source('all_functions.R')

# controls wt
all_con_beta_wt <- readRDS('../../Data/all_con_beta_wt.rda')
all_con_beta_wt_combat <- readRDS('../../Data/all_con_beta_wt_combat.rda')

# all_con_m_wt <- readRDS('../../Data/all_con_m_wt.rda')
# all_con_m_wt_combat <- readRDS('../../Data/all_con_m_wt_combat.rda')

remove_age_lm <- function(temp_controls) {
  pval_list <- list()
  temp_age <- temp_controls$age_sample_collection
  cg_start <- which(grepl('cg', colnames(temp_controls)))[1]
  temp_controls <- temp_controls[, cg_start:(ncol(temp_controls)) ]
  cg_names <- colnames(temp_controls)
    
  # estimate linear model 
  for (i in 1:ncol(temp_controls)){
    temp_cg <- temp_controls[, i]
    pval_list[[i]] <- tidy(lm(temp_age ~  temp_cg))$p.value[1] 
    print(i)
  }
  p_val <- do.call(rbind, pval_list)
  p_val_probes <-as.data.frame(cbind(p_val = p_val, cg_names = cg_names), stringsAsFactors = FALSE)
  colnames(p_val_probes) <- c('p_val', 'probes')
  p_val_probes$p_val <- as.numeric(p_val_probes$p_val)
  probes <- p_val_probes$probes[p_val_probes$p_val < 0.001]
  return(probes)
}

probes <- remove_age_lm(all_con_beta_wt_combat)

saveRDS(probes, '../../Data/age_cgs_lm.rda')
