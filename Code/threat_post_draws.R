#' Extract posterior predictions for each threat/threat combination
#'
#' @param model The fitted brms model object.
#' @param threat_comns String. Vector of threats to extract posterior predictions for. Options include: single string of known threat (e.g. `"climatechange", "pollution.habitatl","climatechange.invasive.exploitation"`) which sets all individual threats and their combinatorial factor to `1` (synergisms), a combination of threats (e.g. `"climatechange + invasive"`, `"invasive + habitatl + disease"`) which only includes the singular threats (additive effect), or `"none"` which sets all threats to `0`, 
#' @param ndraws Numeric. Number of draws to take from the posterior.
#' @param n.cores Numeric. Number of cores for parallelisation using parallel
#' @param ... Additional arguments to pass to `prep_data_grid()`.
#' @returns A data.frame containing posterior timeseries for the desired threat subsets.
#' 
#' @examples
#' post_draws <- threat_post_draws(model = brms_mod1, threat_comns = c("none","climatechange","pollution.habitatl","climatechange + invasive"))
require(foreach)
source("Code/prep_data_grid_fn.R")

threat_post_draws <- function(model, threat_comns, ndraws = NULL, n.cores = 4, ...){
  
  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)
  
  out <- foreach::foreach(x = threat_comns, .combine = "rbind",.packages=c("brms","dplyr","tidyr"), .export = "prep_data_grid") %dopar% {
    
    nw_data <- prep_data_grid(data = model$data,
                              threat = x,
                              ...) %>% #generate the data grid
      dplyr::mutate(threat = x) #name the data grid
    
    return(
      brms::posterior_epred(model,
                          newdata = nw_data,
                          re.form = NA,
                          ndraws = ndraws,
                          incl_autocor = FALSE) %>% #extract posterior draws for the above data grid
      as.data.frame() %>%
      dplyr::mutate(.draw = 1:NROW(.)) %>% #label the draw
      tidyr::pivot_longer(-".draw",names_to = ".timepoint",values_to = ".value",
                          names_transform = function(x){as.numeric(gsub("V","",x))}) %>% #pivot in to long format (ggdist friendly) and tidy up colum names
      cbind(nw_data, row.names = NULL) #bind the data grid for completeness
    )
    
  }
  
  parallel::stopCluster(cl)
  return(out)
  
  # out <- do.call("rbind",
  #                lapply(threat_comns,function(x){ #for each threat/combination of threats in `threat_comns`
  #                  
  #                  nw_data <- prep_data_grid(data = model$data,
  #                                            threat = x,
  #                                            ...) %>% #generate the data grid
  #                    mutate(threat = x) #name the data grid
  #                  
  #                  brms::posterior_epred(model,
  #                                        newdata = nw_data,
  #                                        re.form = NA) %>% #extract posterior draws for the above data grid
  #                    as.data.frame() %>%
  #                    dplyr::mutate(.draw = 1:NROW(.)) %>% #label the draw
  #                    tidyr::pivot_longer(-".draw",names_to = ".timepoint",values_to = ".value",
  #                                        names_transform = function(x){as.numeric(gsub("V","",x))}) %>% #pivot in to long format (ggdist friendly) and tidy up colum names
  #                    cbind(nw_data) #bind the data grid for completeness
  #                  
  #                })) #all threat combinations are rbound together
  # 
  # return(out) 
}
