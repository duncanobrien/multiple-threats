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

threat_counterfac_draws <- function(model, threat_comns, 
                                    ndraws = 1000, n.cores = 4, ...){
  
   cl <- parallel::makeCluster(n.cores)
   doParallel::registerDoParallel(cl)
  
  # if(any(!threat_comns %in% c("pollution","habitatl","climatechange","invasive", "exploitation","disease"))){
  #   stop("`threat_comns` may only contain: 'pollution','habitatl','climatechange','invasive','exploitation', or 'disease'")
  # }
  threatcols <- colnames(model$data)[grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),collapse = "|"),colnames(model$data))] #find the columns in the original dataframe that are threat columns
  
  # out <- foreach::foreach(x = threat_comns, .combine = "rbind",.packages=c("brms","dplyr","tidyr"), .export = "threatcols") %dopar% {
  #   
  #   counterfac_cols <- threatcols[grepl(x,threatcols)]
  #   
  #   nw_data <- model$data %>%
  #     mutate(across(counterfac_cols,~factor(.x,levels = c("0","1")))) %>% #generate the data grid
  #     dplyr::mutate(counterfac = x) #name the data grid
  #   
    # return(
    #   brms::posterior_epred(model,
    #                         newdata = nw_data,
    #                         re.form = NA,
    #                         ndraws = ndraws) %>% #extract posterior draws for the above data grid
    #     as.data.frame() %>%
    #     dplyr::mutate(.draw = 1:NROW(.)) %>% #label the draw
    #     tidyr::pivot_longer(-".draw",names_to = ".timepoint",values_to = ".value",
    #                         names_transform = function(x){as.numeric(gsub("V","",x))}) %>% #pivot in to long format (ggdist friendly) and tidy up colum names
    #     cbind(nw_data, row.names = NULL) #bind the data grid for completeness
    # )
    
    # return(
    #   brms::posterior_epred(model,
    #                         newdata = nw_data,
    #                         re.form = NA,
    #                         ndraws = ndraws) %>% #extract posterior draws for the above data grid
    #     apply(MARGIN = 1, FUN = function(x){mean(diff(x)/diff(seq_along(x)))}) %>%
    #     as.data.frame() %>%
    #     #cbind(nw_data, row.names = NULL) %>%
    #     mutate(counterfac = x,
    #            .draw = 1:NROW(.)) %>%
    #     rename(.value = ".")
    # )
    
    
  # }
  out <- foreach::foreach(y = threat_comns, .combine = "rbind",
                          .packages=c("brms","dplyr","tidyr"), 
                          .export = "threatcols") %dopar% {
  
  #out <- do.call("rbind",lapply(threat_comns,function(y){
    
    # counterfac_cols <- threatcols[grepl(y,threatcols)]
    
    # thrts <-unlist(strsplit(y,"[.]"))
    # if(length(thrts) > 1){
    #   selcols <- unique(c(thrts,combn(thrts,2,FUN = function(x){paste(x,collapse = ".")}),y)) #find all combinations of those threats
    # }else{
    #   selcols <- thrts
    # }
    # counterfac_cols <- threatcols[grepl(paste(paste0("^",selcols,"$"),collapse = "|"),threatcols)]

    #counterfac_cols <- threatcols[grepl(y,threatcols)]
    
   # mod_data <- model$data %>%
   #   group_by(series) %>%
   #   filter(pollution == "1" |
   #            habitatl == "1" |
   #            climatechange  == "1" |
   #            invasive  == "1" |
   #            exploitation  == "1" |
   #            disease  == "1")
   
    nw_data <- model$data %>%
      mutate(across(all_of(y),
                    ~factor("0",levels = c("0","1")))) %>% #generate the data grid
      dplyr::mutate(counterfac = y)  #name the data grid

    tmp <- brms::posterior_epred(model,
                                 newdata = nw_data,
                                 re.form = NULL,
                                 ndraws = ndraws,
                                 sort = TRUE,
                                 incl_autocor = FALSE) %>% #extract posterior draws for the above data grid
      t()
    tmp <- lapply(split(tmp,nw_data$series), matrix, ncol=NCOL(tmp))  #Split into different time series
    
    tmp <- sapply(tmp, FUN = function(y){ 
      apply(y,MARGIN = 2, FUN = function(x){
        #mean(diff(x)/diff(seq_along(x)))
        
        lx <- length(x)
        diff(c(x[1],x[lx]))/lx
        
        #lm(x~seq_along(x))$coefficients[2]
        })}) %>%
      as.data.frame() %>% 
      mutate(.draw = 1:NROW(.)) %>% 
      pivot_longer(-.draw,
                   names_to = "series",
                   values_to = ".value") %>% 
      merge(y = dplyr::distinct(dplyr::select(nw_data,
                                              -c(y_centered,scaled_year,time))),
            by = "series") 
  
    return(tmp)
  }
  
  parallel::stopCluster(cl)
  return(out)
}
