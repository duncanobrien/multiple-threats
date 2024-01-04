#' Generate data grid for interactive threat models 
#'
#' @param data The data.frame used to fit the bayesian model.
#' @param trend String. Variable name of the trend variable.
#' @param threat String. Desired threat/threat combination to include. E.g. `"pollution"` or `"pollution.invasive"` or `"pollution.climatechange.invasive"`. If an interactive threat, all individual and pairwise combinations also included.
#' @param all_threats Numeric. The time-delay offset to use for time delay embedding. Suggested to be positive here, but if not provided, is set to 10\% the length of the time series.
#' @param nuisance String vector. Additional variables to be excluded from the posterior (i.e. set to NA). Typically the random effects.
#' @param nuisance_value What should the `nuisance` columns be set to? `NA` is default but `NULL` may also be used.
#'
#' @returns A data.frame representing the data grid.
#'
prep_data_grid <- function(data,
                           trend = "scaled_year",
                           threat,
                           all_threats = c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                           nuisance = c("series"),
                           nuisance_value = NA){
  
  threatcols <- colnames(data)[grepl(paste(all_threats,collapse = "|"),colnames(data))] #find the columns in the original dataframe that are threat columns
  
  newdata <- data.frame(tr = seq(min(data[[trend]]),max(data[[trend]]),by=1)) %>% #initiate the data grid by setting a default time period to predict over using the original dataset (min to max)
    select(!!sym(trend) := tr) #rename the temporary `$tr` column with the correct name used by the model
  newdata <- cbind(newdata,setNames(lapply(nuisance, function(x){x=nuisance_value}), nuisance)) #add nuisance variables in to the data grid and set to NA. These are typically random effects and must be included so model doesn't break
  ##^^ this may need editing when random slopes are added to model
  
  if(threat == "none"){ #if desired threat is `"none"`, do not select any columns 
    selcols <- ""
    nullcols <- threatcols
  }else{
    if(grepl("\\+",threat)){ #if additive threats are desired
      thrts <- trimws(unlist(strsplit(threat,"\\+"))) #split the `threat` string in to its separate threats and remove whitespace
    }else{ #else if threat is singular or synergistic
      if(!threat %in% threatcols){ #error handling when threat not found
        stop(paste(threat,"not in modelled data"))
      }
      thrts <-unlist(strsplit(threat,"[.]")) #split the `threat` string in to its separate threats
      
      if(length(thrts) > 1){ #if more than one threat in the string
        thrts <- c(thrts,combn(thrts,2,FUN = function(x){paste(x,collapse = ".")}),threat) #find all combinations of those threats
      }
    }

  thrts <- paste0("^",unique(thrts),"$") #paste leading "^" and lagging "$" so the regex will perform perfect matching
  selcols <- threatcols[grepl(paste0(thrts,collapse = "|"),threatcols)] #find the threat columns containing the desired threats 
  nullcols <- threatcols[!grepl(paste0(thrts,collapse = "|"),threatcols)] #and the undesired threat columns
  
  newdata <- cbind(newdata, setNames(lapply(selcols, function(x){x=factor("1",levels = c("0","1"))}), selcols)) #set the desired threat columns to "1"
  }
  newdata <- cbind(newdata, setNames(lapply(nullcols, function(x){x=factor("0",levels = c("0","1"))}), nullcols)) %>% #and the undesired threat columns to "0"
    mutate(time = seq_along(scaled_year)) #before adding a final nuisance variable due to the autocorrelation term
  
  return(newdata)
  
}
