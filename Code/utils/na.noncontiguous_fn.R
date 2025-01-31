#' Find a consecutive sequence of non-NAs
#'
#' Using run length elements (rle), find a sequence greater than or equal to a threshold and set all values outside of this sequence to NA. If multiple contiguous runs are found, selects the last as the data quality is likely more reliable.
#'  
#' @param x Vector of abundance measurements
#' @param consecutive Minimum length of unbroken abundance values.
#'
#' @returns Dataframe
 
na.noncontiguous <- function(x, consecutive) {
  na.rle <- rle(!is.na(x))
  na.rle$values <- na.rle$values & (na.rle$lengths >= consecutive) 
  na.rle$values <- na.rle$values & suppressWarnings(na.rle$lengths == max(na.rle$lengths[na.rle$values],na.rm = TRUE)) #in case of multiple matches longer than duration, extract longest
  csum <- cumsum(na.rle$lengths)
  run.ind <- which(na.rle$values)
  run.ind.len <- length(run.ind)
  if(run.ind.len>1){
    #max.len <- max(na.rle$lengths[run.ind])
    #run.ind <- na.rle$lengths[run.ind][na.rle$lengths[run.ind] == max(na.rle$lengths[run.ind])]
    run.ind <- run.ind[run.ind.len] #take last run as more reliable
  }
  if(run.ind.len == 0){
    x <- NA
  }else if(run.ind !=1){
    x[-((csum[run.ind-1]+1):csum[run.ind])] <- NA
  }else{
    #x[-(csum[1]:csum[run.ind])] <- NA
    x[-(1:csum[1])] <- NA
  }
  
  return(x)
}
