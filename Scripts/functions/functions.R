# functions to use in LFS 
is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x) # takes all elements that are not null
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)}
#####

clinMatch <- function(clin, methylation_id){
  
  store_list <- list()
  
  for(i in 1:nrow(clin)){   
    print(i)
    ids <- as.character(clin[,1 ])
    id_temp <- ids[i]  
    
    if(grepl("/", id_temp)){
      sep <- unlist(strsplit(id_temp, "/"))
      
      if(any(sep %in% methylation_id)){
        store_list[[i]] <- id_temp
        store_list <- rmNullObs(store_list)
        print("id found here")
      }
    } 
    
    sep_none <- id_temp
    if(any(sep_none %in% methylation_id)){
      print("already found id")
      
    } else {
      print("no extra ids found")
    }  
  }
  id_ind <- append(unlist(store_list), methylation_id)
  methyl_clin <- clin[clin$V1 %in% id_ind,]
  return(methyl_clin)
  
}
