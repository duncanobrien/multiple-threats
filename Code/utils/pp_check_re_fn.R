#' Posterior predictive check for a brms model that correctly controls for random effect effects
#'  
#' @param object brmsfit
#' @param ndraws number of draws
#'
#' @returns pp_check plot

pp_check_re <- function(object, ndraws = 10,...){
  ranefs <- object$ranef
  re <- unique(ranefs$group)
  re <- na.omit(re[ranefs$cov == ""]) #if there is a specified covariance matrix, can't provide new levels
  resp <- object$formula$formula[[2]]
  
  yrep <- object |> 
    brms::posterior_predict(newdata = object$data |> 
                              dplyr::mutate(dplyr::across(dplyr::all_of(re),~ paste("rep",as.numeric(factor(.x)),sep = "_"))),
                            allow_new_levels = TRUE,
                            sample_new_levels = unique(ranefs$dist),
                            ndraws = ndraws, 
                            ...)
  
  bayesplot::ppc_dens_overlay(object$data[[resp]], yrep)
  
}
