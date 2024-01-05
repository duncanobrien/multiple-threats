# --------------------------------------------------------------------------------------- #
# - FILE NAME:   counterfactual.R         
# - DATE:        03/01/2024
# - DESCRIPTION: Code to test scenarios where the threats have been removed.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila@ub.edu), Duncan O'brien ()
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(brms)
library(tidyverse)
library(tidybayes)
library(parallel)
library(doSNOW)

# Load the model 

m1 <- read_rds("Results/simple-slopes-intvadd2.RDS")

# Load the data 

load("Data/LivingPlanetData2.RData")

# Conterfactual tests ----------------------------------------------------------

# Extract the data used to build the model

pop_data <- m1$data 

# Store the names of threats

threats <- colnames(pop_data)[grepl(paste(c("pollution","habitatl",
                                               "climatechange","invasive", 
                                               "exploitation","disease"),
                                             collapse = "|"),
                                       colnames(pop_data))] 

# Obtain the predicted value for each population 

pop_pred <- posterior_epred(m1,
                            newdata = m1$data,
                            re.form = NULL,
                            ndraws = 1000) %>% #extract posterior draws for the above data grid
  t() %>%
  as.data.frame() %>%
  split(m1$data$series) %>%
  sapply(FUN = function(y){ 
    apply(y,MARGIN = 2, FUN = function(x){mean(diff(x)/diff(seq_along(x)))})}) %>%
  as.data.frame() %>%
  mutate(.draw = 1:n()) %>%
  pivot_longer(-.draw,names_to = "series",values_to = ".value") %>%
  merge(y = dplyr::distinct(dplyr::select(m1$data,-c(y_centered,scaled_year,time))),
        by = "series") %>% 
  # Join with the original data to get a threats column
  left_join(dd_long %>% 
              # Keep only one ID
              distinct(ID,.keep_all = T) %>% 
              # Select only the columns we are interested in
              select(threats, ID) %>% 
              # Transform ID in a character so it can be compared against series
              mutate(ID=as.character(ID)), 
            by=c("series"="ID")) %>% 
  # Count the threats and correct the column threats to match the initial names
  mutate(number=str_count(threats,"[A-Z]"),
         threats=tolower(threats), 
         threats=gsub(":", ".", threats),
         threats=gsub("habitat loss", "habitatl", threats), 
         threats=gsub("invasive spp/genes", "invasive", threats),
         threats=gsub("climate change","climagechange", threats),
         threats=replace_na(threats,"none"))

# Obtain the mean coefficient for each threat 

mean_trends <- m1 %>% 
  gather_draws(`b_.*`, regex = TRUE) %>%
  mean_qi(.width = .95) %>%
  mutate(.variable = gsub("b_", "", .variable), 
         .variable = gsub("threats", "", .variable), 
         .variable = gsub("scaled_year", "year", .variable),
         .variable=gsub("1", "", .variable),
         number=str_count(.variable, '\\.')+1) %>% 
  filter(grepl("year:", .variable))

# Substract the mean coefficient from the predicted values
# We will parallelise the code, so this will take multiple steps

# Define the threats

threats <- colnames(pop_data)[3:39]

# Define the number of cores to use

clus  <- makeCluster(detectCores() - 1)

# register the cluster with doSNOW

registerDoSNOW(clus)

# load in the options using .options.snow

trend_scenarios <- foreach(i = threats, 
                        .combine = 'rbind',
                        .packages = c("tidyverse", "naniar")) %dopar% {
                          
                          # Extract the number of threats
                          
                          number_threats <- length(str_split(i, "\\.")[[1]])
                          
                          # If the number of threats is equal to 1
                          
                          if(number_threats==1){
                            # Create a new dataframe were we substract the mean 
                            # effect of the threats on the population trends 
                            new <- pop_pred %>% 
                              filter(grepl(i,threats)) %>% # Filter by the threat
                              # Substract the mean effect of the single threat
                              mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",i)]) %>% 
                              # Bind to the previous dataset
                              rbind(pop_pred  %>% 
                                      filter(!grepl(i,threats)) %>% 
                                      mutate(new_trend=NA)) %>% 
                              # Add the name of the scenario 
                              mutate(scenario=paste("No", i)) 
                            
                            # If the number of threats is equal to 2
                            
                            } else if(number_threats==2){
                              # Split the names of the interacting threats
                              x2 <- unlist(strsplit(i, 
                                                    ".", 
                                                    fixed=T))
                              # Create a new dataset where we remove the effect
                              # of the two threats, when interacting together
                              new1 <- pop_pred  %>%
                                # Filter for those populations exposed to both
                                # threats
                                filter(grepl(x2[1],threats), 
                                       grepl(x2[2],threats), number!=1) %>% 
                                # Extract the effect of both threats
                                mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",i)])
                              # Create a second dataset where we remove the effect
                              # of the first threat, only 
                              new2 <- pop_pred  %>% 
                                # Filter for the first threat, only those with
                                # more and less than two threats and that do not
                                # overlap with the series in new1
                                filter(grepl(x2[1],threats), 
                                       number!=2, !grepl(i,threats),
                                       !series%in%new1$series) %>%
                                # Substract the effect of the single threat
                                mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",x2[1])]) 
                              # Same as above but with the second threat
                              new <- pop_pred %>% 
                                # Notice we also filer out series in new1 and 
                                # new2
                                filter(grepl(x2[2],threats), number!=2,
                                       !grepl(i,threats), 
                                       !series%in%new1$series,
                                       !series%in%new2$series) %>% 
                                mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",x2[2])])  %>% 
                                # Bind with new1, new2, and pop_pred 
                                rbind(new1) %>% 
                                rbind(new2) 
                              
                              new <- new %>% 
                                rbind(pop_pred  %>% 
                                      filter(!series%in%new$series) %>% 
                                        mutate(new_trend=NA)) %>% 
                              mutate(scenario=paste("No", i))
                              
                              # If the number of threats is equal to 3
                              
                              }else if(number_threats==3){
                                # Separate the threats
                                x3 <- unlist(strsplit(i, 
                                                      ".", 
                                                      fixed=T))
                                # Create new names for the two interactive threats
                                two1 <- paste0(x3[1], ".", x3[2]) 
                                two2 <- paste0(x3[2], ".", x3[3]) 
                                two3 <- paste0(x3[1], ".", x3[3])
                                # Create a data frame where we subtract the 
                                # effect of three threats
                                new1 <- pop_pred  %>% 
                                  filter(grepl(x3[1],threats), 
                                         grepl(x3[2],threats),
                                         grepl(x3[3],threats), 
                                         number==3) %>% 
                                  mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",i)]) 
                                # Subset the population trends with two threats
                                # First and second threat
                                two_threat1 <- mean_trends %>%  
                                  filter(grepl(x3[1], .variable),
                                         grepl(x3[2], .variable),
                                         number==2)
                                # Second and third threat
                                two_threat2 <- mean_trends %>%  
                                  filter(grepl(x3[2], .variable),
                                         grepl(x3[3], .variable),
                                         number==2)
                                # First and third threat
                                two_threat3 <- mean_trends %>%  
                                  filter(grepl(x3[1], .variable),
                                         grepl(x3[3], .variable),
                                         number==2)
                                # If the two threats do not exist 
                                if(dim(two_threat1)[1]==0){
                                  # When there is no value we create a fake one
                                  # with new_trend == NA
                                  new2 <- tryCatch({pop_pred  %>% 
                                      filter(grepl(x3[1], threats),
                                             grepl(x3[2], threats),
                                             !ID%in%new1$ID) %>% 
                                      mutate(new_trend=NA)},
                                      error=function(e){})
                                  # When there are some we do the usual procedure
                                  }else{new2 <- tryCatch({pop_pred  %>% 
                                      filter(grepl(x3[1], threats),
                                             grepl(x3[2], threats),
                                             !series%in%new1$series) %>% 
                                      mutate(new_trend=.value-two_threat1$.value)},
                                      error=function(e){})}
                                # Same here with other two threats
                                if(dim(two_threat2)[1]==0){
                                  new3 <- tryCatch({pop_pred  %>% 
                                      filter(grepl(x3[2], threats),
                                             grepl(x3[3], threats),
                                             !series%in%new1$series) %>% 
                                      mutate(new_trend=NA)},
                                      error=function(e){})
                                  }else{new3 <- tryCatch({pop_pred  %>% 
                                      filter(grepl(x3[2], threats),
                                             grepl(x3[3], threats),
                                             !series%in%new1$series) %>% 
                                      mutate(new_trend=.value-two_threat2$.value)},
                                      error=function(e){})}
                                # Same here with other two threats
                                if(dim(two_threat3)[1]==0){ 
                                  new4 <- tryCatch({pop_pred  %>% 
                                      filter(grepl(x3[1], threats),
                                             grepl(x3[3], threats),
                                             !series%in%new1$series) %>% 
                                      mutate(new_trend=NA)},
                                      error=function(e){})
                                }else{
                                  new4 <- tryCatch({pop_pred  %>% 
                                      filter(grepl(x3[1], threats),
                                             grepl(x3[3], threats),
                                             !series%in%new1$series) %>% 
                                      mutate(new_trend=.value-two_threat3$.value)},
                                      error=function(e){})  
                                }
                                
                                # Now we remove the single effects
                                # For first threat
                                new5 <- pop_pred  %>% 
                                  filter(grepl(x3[1], threats),
                                         !series%in%new1$series, 
                                         !series%in%new2$series,
                                         !series%in%new4$series) %>% 
                                  mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",x3[1])])
                                # For second threat
                                new6 <- pop_pred  %>% 
                                  filter(grepl(x3[2], threats),
                                         !series%in%new1$series, 
                                         !series%in%new2$series,
                                         !series%in%new3$series) %>% 
                                  mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",x3[2])])
                                # Final subset
                                new <- pop_pred  %>% 
                                  filter(grepl(x3[3], threats),
                                         !series%in%new1$series, 
                                         !series%in%new3$series,
                                         !series%in%new4$series) %>% 
                                  mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",x3[3])]) %>% 
                                  # Bind to the previous datasets
                                  rbind(new1) %>% 
                                  rbind(new2) %>% 
                                  rbind(new3) %>% 
                                  rbind(new4) %>% 
                                  rbind(new5) %>% 
                                  rbind(new6) 
                                new <- new %>% 
                                  rbind(pop_pred  %>% 
                                          filter(!series%in%new$series)) %>% 
                                  # Add the scenario name
                                  mutate(scenario=paste("No", i)) 
                                }
                        }

# Stop the parallelisation

stopCluster(clus)

# Save the object 

save(trend_scenarios, "Results/counterfactuals.RData") 

# # Add the control scenario 
# 
# scenarios_data <- 
#   pop_perd %>% mutate(scenario="No Management") %>% 
#   bind_rows(mu_scenarios) %>% 
#   filter(threats!="none") %>% 
#   mutate(prop=((log(lambda)-Estimate)/abs(Estimate))*100,
#          prop=round(prop,2),
#          number_scen= str_count(scenario,"[A-Z]")) %>% 
#   group_by(scenario) %>% 
#   mutate(pos=median(lambda), 
#          pos2=median(prop)) %>% 
#   ungroup()  
# 
# # Add the control scenario 
# 
# scenarios_median <- scenarios_data %>% 
#   group_by(scenario) %>% 
#   summarise(mcontrol=median(Estimate),
#             mcounter=median(log(lambda))) %>% 
#   mutate(prop=((mcounter-(mcontrol))/abs(mcontrol))*100,
#          prop=round(prop,2))

