#' Extract posterior predictions for each threat/threat combination
#'
#' @param model The fitted brms model object.
#' @param modelled_data String. Variable name of the trend variable.
#' @param threat_comns String. Vector of threats to extract posterior predictions for. Options include: single string of known threat (e.g. `"climatechange", "pollution.habitatl","climatechange.invasive.exploitation"`) which sets all individual threats and their combinatorial factor to `1` (synergisms), a combination of threats (e.g. `"climatechange + invasive"`, `"invasive + habitatl + disease"`) which only includes the singular threats (additive effect), or `"none"` which sets all threats to `0`, 
#' @param ... Additional arguments to pass to `prep_data_grid()`.
#' @returns A data.frame containing posterior timeseries for the desired threat subsets.
#' 
#' @examples
#' post_draws <- threat_post_draws(model = brms_mod1, mod_data, threat_comns = c("none","climatechange","pollution.habitatl","climatechange + invasive"))

source("Code/prep_data_grid_fn.R")
threat_post_draws <- function(model, modelled_data, threat_comns,...){
  
  out <- do.call("rbind",lapply(threat_comns,function(x){
    
    nw_data <- prep_data_grid(data = modelled_data,
                              threat = x,
                              ...) %>%
      mutate(threat = x)
    
    brms::posterior_epred(model,
                          newdata = nw_data,
                          re.form = NA) %>%
      as.data.frame() %>%
      dplyr::mutate(.draw = 1:NROW(.)) %>%
      tidyr::pivot_longer(-".draw",names_to = ".timepoint",values_to = ".value",
                          names_transform = function(x){as.numeric(gsub("V","",x))}) %>%
      cbind(nw_data)
    
  }))
  
  return(out)
}
