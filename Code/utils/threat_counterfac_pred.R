#' Extract posterior predictions following counterfactual
#'
#' @param model The fitted brms model object.
#' @param threat_comns String. Vector of threats to extract posterior predictions for. Options include: string of single threat(s) (e.g. `c("climatechange","pollution")`).
#' @param ndraws Numeric. Number of draws to take from the posterior.
#' @param n.cores Numeric. Number of cores for parallelisation using parallel.
#' @param ... Additional arguments to pass to `prep_data_grid()`.
#' @returns A data.frame containing posterior timeseries for the desired threat subsets.
#' 
#' @examples
#' post_draws <- threat_post_draws(model = brms_mod1, threat_comns = c("climatechange","pollution","invasive"))
require(foreach)

threat_counterfac_pred <- function(model, threat_comns,
                                    re.form = NULL,
                                    ndraws = 1000, n.cores = 4,...){
  
  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)
  
  threatcols <- colnames(model$data)[grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),collapse = "|"),colnames(model$data))] #find the columns in the original dataframe that are threat columns
  

  out <- foreach::foreach(y = threat_comns, .combine = "rbind",
                          .packages=c("brms","dplyr","tidyr"), 
                          .export = "threatcols") %dopar% {
                            
                            # Split the threats
                            thrts <-unlist(strsplit(y,"[.]"))
                            # Find the threats combinations
                            counterfac_cols <- threatcols[grepl(paste(thrts, collapse = "|"), threatcols)]
                            
                            nw_data <- model$data %>%
                              # Filter out the populations not exposed to threat
                              filter_at(threatcols, any_vars(. != "0")) %>% 
                              # Set the
                              mutate(across(all_of(counterfac_cols),
                                            ~factor("0",levels = c("0","1")))) %>% #generate the data grid
                              dplyr::mutate(counterfac = paste("No", y))  #name the data grid
                            # Make the predictons 
                            tmp <- brms::posterior_epred(model,
                                                         newdata = nw_data,
                                                         re.form = re.form,
                                                         incl_autocor = FALSE,
                                                         sort = TRUE,
                                                         ndraws = ndraws) %>% #extract posterior draws for the above data grid
                              as.data.frame() %>%
                              mutate(.draw = 1:NROW(.)) %>%
                              #extract posterior draws for the data used to create the model
                              pivot_longer(-.draw,names_to = "index",values_to = ".value") %>%
                              cbind(nw_data %>% dplyr::select(series, time), row.names = NULL) %>% 
                              mutate(counterfac=paste("No", y))
                            return(tmp)
                          }
  
  parallel::stopCluster(cl)
  return(out)
}
