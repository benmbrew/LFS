evaluateSurvival <- function(clinicalData, labels) {
  # Retrieve the patient survival times and death status
  survTime <- clinicalData$days_to_death
  deathStatus <- clinicalData$vital_status == "dead"
  
  # Replace missing survival times with days to last follow up
  missingSurvInd <- is.na(survTime)
  lastFollowup <- clinicalData$days_to_last_followup
  survTime[missingSurvInd] <- lastFollowup[missingSurvInd]
  
  # Calculate the p-value of the log-rank test and the concordance
  # index of the Cox proportional hazards model
  survObject <- Surv(survTime, deathStatus)
  survDiff <- survdiff(survObject~labels)
  pval <- 1 - pchisq(survDiff$chisq, length(survDiff$n)-1)
  coxphFit <- coxph(survObject~labels)
  ci <- summary(coxphFit)$concordance
  #pvalcox <- summary(coxphFit)$logtest["pvalue"]
  #   coefcox <- coxphFit$coefficients
  #   
  #   # method 4
  #   library(survcomp)
  #   missInd <- !is.na(survTime)
  #   fit <- coxph(survObject ~ labels)
  #   coxPredict <- predict(fit, type="risk")  
  #   con_index <- concordance.index(x=coxPredict, surv.time=survTime[missInd], 
  #                     surv.event=deathStatus[missInd], method="noether")
  #   con_index_ci <- con_index$c.index
  #   con_index_p <- con_index$p.value
  #   
  #   # method 5
  #   library(rms)
  #   fit.cph <- cph(survObject ~ labels, x=TRUE, y=TRUE, surv=TRUE)
  #   # Get the Dxy
  #   v <- validate(fit.cph, dxy=TRUE, B=1000)
  #   Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
  #   orig_Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.orig"]
  #   
  #   # The c-statistic according to Dxy=2(c-0.5)
  #   bias_corrected_c_index  <- abs(Dxy)/2+0.5
  #   
  results <- c(pval, ci)
  
  return(results)
}

evaluateSurvivalCensor <- function(clinicalData, labels) {
  # Retrieve the patient survival times and death status
  survTime <- clinicalData$days_to_death
  deathStatus <- clinicalData$vital_status == "dead"
  
  # Replace missing survival times with days to last follow up
  missingSurvInd <- is.na(survTime)
  lastFollowup <- clinicalData$days_to_last_followup
  survTime[missingSurvInd] <- lastFollowup[missingSurvInd]
  
  # Calculate the p-value of the log-rank test and the concordance
  # index of the Cox proportional hazards model
  survObject <- Surv(survTime, deathStatus)
  survDiff <- survdiff(survObject~labels)
  pval <- 1 - pchisq(survDiff$chisq, length(survDiff$n)-1)
  coxphFit <- coxph(survObject~labels)
  ci <- summary(coxphFit)$concordance[1]
  results <- c(pval, ci)
  
  return(results)
}
