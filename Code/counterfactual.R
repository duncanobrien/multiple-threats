# --------------------------------------------------------------------------------------- #
# - FILE NAME:   counterfactual.R         
# - DATE:        03/01/2024
# - DESCRIPTION: Code to test scenarios where the threats have been removed.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(brms)
library(tidyverse)
library(tidybayes)

# Conterfactual tests ----------------------------------------------------------
# Load the model 

m1 <- read_rds("Results/simple-slopes_b.RDS")

# Extract the data used to build the model

pop_data <- m1$data %>% 
  mutate(threats = gsub("Invasive spp/genes", "Invasive", threats),
         n.threats = str_count(threats,"[A-Z]"))

# Obtain the predicted value for each population 

pop_perd <- fitted(m1, scale = "response")

# Edit the data to have the necessary variables and join with data values

pop_trend <- pop_perd %>%
  as.data.frame() %>% 
  bind_cols(pop_data[c("threats","time","series")]) %>% # Add the threat affecting the population 
  reframe(mean_trend = mean(diff(Estimate)/diff(time)),
          .by = c(series, threats)) %>%  # Calculate the mean population trend
  mutate(number = str_count(threats,"[A-Z]"),# Get the number of threats
         new_trend=NA) # Create a new column to store the counterfactuals

# Obtain the mean coefficient for each threat 

mean_trends <- m1 %>% 
  gather_draws(`b_.*`, regex = TRUE) %>%
  mean_qi(.width = .95) %>%
  mutate(.variable = gsub("b_", "", .variable), 
         .variable = gsub("threats", "", .variable), 
         .variable = gsub("scaled_time", "time", .variable),
         .variable = gsub("InvasivesppDgenes", "Invasive", .variable),
         .variable = gsub("Habitatloss", "Habitat loss", .variable),
         .variable = gsub("Climatechange", "Climate change", .variable),
         number = str_count(.variable,"[A-Z]")) 

# Substract the mean coefficient from the predicted values
# We will parallelise the code, so this will take multiple steps

# Define the  threats

threats <- unique(pop_data$threats[which(pop_data$n.threat>=1)]) %>% 
  as_tibble() %>% 
  mutate(number=str_count(value,"[A-Z]"))

# Define the number of cores to use

clus  <- makeCluster(detectCores() - 1)

# register the cluster with doParallel

registerDoSNOW(clus)

# load in the options using .options.snow
trend_scenarios <- foreach(i = threats$value, 
                        .combine = 'rbind',
                        .packages = c("tidyverse", "naniar")) %dopar% {
                          
                          # Extract the number of threats
                          
                          number_threats <- threats$number[threats$value==i]
                          
                          # If the number of threats is equal to 1
                          
                          if(number_threats==1){
                            # Create a new dataframe were we substract the mean 
                            # effect of the threats on the population trends 
                            new <- pop_trend %>% 
                              filter(grepl(i,threats)) %>% # Filter by the threat
                              # Substract the mean effect of the single threat
                              mutate(new_trend=mean_trend-mean_trends$.value[mean_trends$.variable==paste0("time:",i)]) %>% 
                              # Bind to the previous dataset
                              rbind(pop_trend %>% 
                                      filter(!grepl(i,threats))) %>% 
                              # Add the name of the scenario 
                              mutate(scenario=paste("No", i)) 
                            
                            # If the number of threats is equal to 2
                            
                            } else if(number_threats==2){
                              # Split the names of the interacting threats
                              x2 <- unlist(strsplit(i, 
                                                    ":", 
                                                    fixed=T))
                              # Create a new dataset where we remove the effect
                              # of the two threats, when interacting together
                              new1 <- pop_trend %>%
                                # Filter for those populations exposed to both
                                # threats
                                filter(grepl(x2[1],threats), 
                                       grepl(x2[2],threats), number!=1) %>% 
                                # Extract the effect of both threats
                                mutate(new_trend=mean_trend-mean_trends$.value[mean_trends$.variable==paste0("time:",i)])
                              # Create a second dataset where we remove the effect
                              # of the first threat, only 
                              new2 <- pop_trend %>% 
                                # Filter for the first threat, only those with
                                # more and less than two threats and that do not
                                # overlap with the series in new1
                                filter(grepl(x2[1],threats), 
                                       number!=2, !grepl(i,threats),
                                       !series%in%new1$series) %>%
                                # Substract the effect of the single threat
                                mutate(new_trend=mean_trend-mean_trends$.value[mean_trends$.variable==paste0("time:",x2[1])]) 
                              # Same as above but with the second threat
                              new <- pop_perd %>% 
                                # Notice we also filer out series in new1 and 
                                # new2
                                filter(grepl(x2[2],threats), number!=2,
                                       !grepl(i,threats), 
                                       !ID%in%new1$ID,
                                       !ID%in%new2$ID) %>% 
                                mutate(new_trend=mean_trend-mean_trends$.value[mean_trends$.variable==paste0("time:",x2[2])])  %>% 
                                # Bind with new1, new2, and pop_trend
                                rbind(new1) %>% 
                                rbind(new2) %>% 
                              rbind(pop_trend %>% 
                                      filter(!series%in%new$series)) %>% 
                              mutate(scenario=paste("No", i))
                              
                              # If the number of threats is equal to 3
                              
                              }else if(number_threats==3){
                                # Separate the threats
                                x3 <- unlist(strsplit(i, 
                                                      ":", 
                                                      fixed=T))
                                # Create new names for the two interactive threats
                                two1 <- paste0(x3[1], ":", x3[2]) 
                                two2 <- paste0(x3[2], ":", x3[3]) 
                                two3 <- paste0(x3[1], ":", x3[3])
                                # Create a data frame where we subtract the 
                                # effect of three threats
                                new1 <- pop_trend %>% 
                                  filter(grepl(x3[1],threats), 
                                         grepl(x3[2],threats),
                                         grepl(x3[3],threats), 
                                         number==3) %>% 
                                  mutate(new_trend=mean_trend-mean_trends$.value[mean_trends$.variable==paste0("time:",i)]) 
                                # Subset the population trends with two threats
                                # First and second threat
                                two_threat1 <- mean_trends %>%  
                                  filter(grepl(x3[1], .variable),
                                         grepl(x3[2], .variable),
                                         grepl("time", .variable),
                                         number==2)
                                # Second and third threat
                                two_threat2 <- mean_trends %>%  
                                  filter(grepl(x3[2], .variable),
                                         grepl(x3[3], .variable),
                                         grepl("time", .variable),
                                         number==2)
                                # First and third threat
                                two_threat3 <- mean_trends %>%  
                                  filter(grepl(x3[1], .variable),
                                         grepl(x3[3], .variable),
                                         grepl("time", .variable),
                                         number==2)
                                # If the two threats do not exist 
                                if(dim(two_threat1)[1]==0){
                                  # When there is no value we create a fake one
                                  # with new_trend == NA
                                  new2 <- tryCatch({pop_trend %>% 
                                      filter(grepl(x3[1], threats),
                                             grepl(x3[2], threats),
                                             !ID%in%new1$ID) %>% 
                                      mutate(new_trend=NA)},
                                      error=function(e){e <- new1$lambda=NA})
                                  # When there are some we do the usual procedure
                                  }else{new2 <- tryCatch({pop_trend %>% 
                                      filter(grepl(x3[1], threats),
                                             grepl(x3[2], threats),
                                             !series%in%new1$series) %>% 
                                      mutate(new_trend=mean_trend-two_threat1$.value)},
                                      error=function(e){e <- new1$new_trend=NA})}
                                # Same here with other two threats
                                if(dim(two_threat2)[1]==0){
                                  new3 <- tryCatch({pop_trend %>% 
                                      filter(grepl(x3[2], threats),
                                             grepl(x3[3], threats),
                                             !series%in%new1$series) %>% 
                                      mutate(new_trend=NA)},
                                      error=function(e){e <- new1$lambda=NA})
                                  }else{new3 <- tryCatch({pop_trend %>% 
                                      filter(grepl(x3[2], threats),
                                             grepl(x3[3], threats),
                                             !series%in%new1$series) %>% 
                                      mutate(new_trend=mean_trend-two_threat2$.value)},
                                      error=function(e){e <- new1$lambda=NA})}
                                # Same here with other two threats
                                if(dim(two_threat3)[1]==0){ 
                                  new4 <- tryCatch({pop_trend %>% 
                                      filter(grepl(x3[1], threats),
                                             grepl(x3[3], threats),
                                             !series%in%new1$series) %>% 
                                      mutate(new_trend=NA)},
                                      error=function(e){e <- new1$lambda=NA})
                                }else{
                                  new4 <- tryCatch({pop_trend %>% 
                                      filter(grepl(x3[1], threats),
                                             grepl(x3[3], threats),
                                             !series%in%new1$series) %>% 
                                      mutate(new_trend=mean_trend-two_threat3$.value)},
                                      error=function(e){e <- new1$lambda=NA})  
                                }
                                
                                # Now we remove the single effects
                                # For first threat
                                new5 <- pop_trend %>% 
                                  filter(grepl(x3[1], threats),
                                         !series%in%new1$series, 
                                         !series%in%new2$series,
                                         !series%in%new4$series) %>% 
                                  mutate(new_trend=mean_trend-mean_trends$.value[mean_trends$.variable==paste0("time:",x3[1])])
                                # For second threat
                                new6 <- pop_trend %>% 
                                  filter(grepl(x3[2], threats),
                                         grepl("time", threats),
                                         !series%in%new1$series, 
                                         !series%in%new2$series,
                                         !series%in%new3$series) %>% 
                                  mutate(new_trend=mean_trend-mean_trends$.value[mean_trends$.variable==paste0("time:",x3[2])])
                                # Final subset
                                new <- pop_trend %>% 
                                  filter(grepl(x3[3], threats),
                                         grepl("time", threats),
                                         !series%in%new1$series, 
                                         !series%in%new3$series,
                                         !series%in%new4$series) %>% 
                                  mutate(new_trend=mean_trend-mean_trends$.value[mean_trends$.variable==paste0("time:",x3[3])]) %>% 
                                  # Bind to the previous datasets
                                  rbind(new1) %>% 
                                  rbind(new2) %>% 
                                  rbind(new3) %>% 
                                  rbind(new4) %>% 
                                  rbind(new5) %>% 
                                  rbind(new6) %>% 
                                  rbind(pop_trend %>% 
                                          filter(!series%in%new$series)) %>% 
                                  # Add the scenario name
                                  mutate(scenario=paste("No", i)) 
                                }
                        }

# Stop the parallelisation

stopCluster(clus)

# Add the control scenario 

scenarios_data <- 
  pop_perd %>% mutate(scenario="No Management") %>% 
  bind_rows(mu_scenarios) %>% 
  filter(threats!="None") %>% 
  mutate(prop=((log(lambda)-Estimate)/abs(Estimate))*100,
         prop=round(prop,2),
         number_scen= str_count(scenario,"[A-Z]")) %>% 
  group_by(scenario) %>% 
  mutate(pos=median(lambda), 
         pos2=median(prop)) %>% 
  ungroup()  

# Add the control scenario 

scenarios_median <- scenarios_data %>% 
  group_by(scenario) %>% 
  summarise(mcontrol=median(Estimate),
            mcounter=median(log(lambda))) %>% 
  mutate(prop=((mcounter-(mcontrol))/abs(mcontrol))*100,
         prop=round(prop,2))



########################
# Duncan counterfactual dydx
########################

source("Code/threat_counterfac_draws.R")
pop_perd <- brms::posterior_epred(mod_intvadd2,
                                  newdata = mod_intvadd2$data,
                                  re.form = NULL,
                                  ndraws = 1000) %>% #extract posterior draws for the above data grid
  t() %>%
  as.data.frame() %>%
  split(mod_intvadd2$data$series) %>%
  sapply(FUN = function(y){ 
    apply(y,MARGIN = 2, FUN = function(x){mean(diff(x)/diff(seq_along(x)))})}) %>%
  as.data.frame() %>%
  mutate(.draw = 1:n()) %>%
  pivot_longer(-.draw,names_to = "series",values_to = ".value") %>%
  merge(y = dplyr::distinct(dplyr::select(mod_intvadd2$data,-c(y_centered,scaled_year,time))),
        by = "series") %>%
  mutate(counterfac = "none")

counter_fac_data <- threat_counterfac_draws(mod_intvadd2,
                                            threat_comns = c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                                            ndraws = 1000) %>%
  rbind(pop_perd) %>%
  mutate(counterfac = fct_relevel(after = "none"))

counterfac_diff <- do.call("rbind",lapply(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),function(x){
  
  counter_fac_data %>%
    subset(counterfac %in% c(x,"none")) %>% #subset to shared additive and interactive threats (e.g. "threat1.threat2" and "threat1 + threat2")
    arrange(series,.draw,counterfac) %>% #as "none" is reference level, sets none first value of group
    reframe(.value = .value[1]-.value[2], .by = c(series,.draw)) %>% #find difference in derivatives between additive and interactive threats
    mutate(counterfac = x) #name threat combination
  
}))
