#' Generate data grid for interactive threat models 
#'
#' @param data The data.frame used to fit the bayesian model.
#' @param trend String. Variable name of the trend variable.
#' @param threat String. Desired threat/threat combination to include. E.g. `"pollution"` or `"pollution.invasive"` or `"pollution.climatechange.invasive"`. If an interactive threat, all individual and pairwise combinations also included.
#' @param all_threats Numeric. The time-delay offset to use for time delay embedding. Suggested to be positive here, but if not provided, is set to 10\% the length of the time series.
#' @param nuisance String vector. Additional variables to be excluded from the posterior (i.e. set to NA). Typically the random effects.
#'
#' @returns A data.frame representing the data grid.
#'
prep_data_grid <- function(data,
                           trend = "scaled_year",
                           threat,
                           all_threats = c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                           nuisance = c("series")){
  
  threatcols <- colnames(data)[grepl(paste(all_threats,collapse = "|"),colnames(data))]
  if(!threat %in% threatcols){
    stop(paste(threat,"not in modelled data"))
  }
  
  newdata <- data.frame(tr = seq(min(data[[trend]]),max(data[[trend]]),by=1)) %>%
    select(!!sym(trend) := tr)
  newdata <- cbind(newdata,setNames(lapply(nuisance, function(x){x=NA}), nuisance))
  
  thrts <-unlist(strsplit(threat,"[.]"))
  if(length(thrts) > 1){
    thrts <- c(thrts,combn(thrts,2,FUN = function(x){paste(x,collapse = ".")}),threat)
  }
  thrts <- paste0("^",unique(thrts),"$")
  selcols <- threatcols[grepl(paste0(thrts,collapse = "|"),threatcols)]
  nullcols <- threatcols[!grepl(paste0(thrts,collapse = "|"),threatcols)]
  
  newdata <- cbind(newdata, setNames(lapply(selcols, function(x){x=factor("1",levels = c("0","1"))}), selcols))
  newdata <- cbind(newdata, setNames(lapply(nullcols, function(x){x=factor("0",levels = c("0","1"))}), nullcols)) %>%
    mutate(time = seq_along(scaled_year))
  
  return(newdata)
  
}
