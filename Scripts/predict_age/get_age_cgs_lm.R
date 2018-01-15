

# source all_functions.R to load libraries and my functions
source('all_functions.R')

# load data
controls <- readRDS('../../Data/controls_full.rda')

remove_age_lm <- function(temp_controls) {
  pval_list <- list()
  temp_age <- temp_controls$age_sample_collection
  cg_start <- which(grepl('cg', colnames(temp_controls)))[1]
  temp_controls <- temp_controls[, cg_start:(ncol(temp_controls) -1) ]
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
  probes <- p_val_probes$probes[p_val_probes$p_val < 0.05]
  return(probes)
}

probes <- remove_age_lm(controls)

saveRDS(probes, '../../Data/age_cgs_lm.rda')
