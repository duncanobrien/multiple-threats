# --------------------------------------------------------------------------------------- #
# - FILE NAME:   counterfactual.R         
# - DATE:        03/01/2024
# - DESCRIPTION: Code to test scenarios where the threats have been removed.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

require(brms)
require(tidyverse)

# Conterfactual tests ----------------------------------------------------------
# Load the model 

m1 <- read_rds("Results/simple-slopes_b.RDS")

# Extract the data used to build the model

pop_data <- m1$data

# Obtain the predicted value for each population 

pop_perd <- fitted(m1, scale = "response")

# Edit the data to have the necessary variables and join with data values

pop_perd2 <- pop_perd %>%
  as.data.frame() %>% 
  bind_cols(pop_data[c("threats", "series")]) %>% # Add the threat affecting the population 
  mutate(number = str_count(threats,"[A-Z]")) # Get the number of threats


# Obtain the mean coefficient for each threat 

mus_mean <- m1 %>% 
  gather_draws(`b_.*`, regex = TRUE) %>%
  mean_qi(.width = .95) %>%
  mutate(.variable = gsub("b_threats", "", .variable), 
         .variable = gsub("InvasivesppDgenes", "Invasive", .variable),
         .variable = gsub("Habitatloss", "Habitat loss", .variable),
         .variable = gsub("Climatechange", "Climate change", .variable),# Create an initial fake population
         number = str_count(.variable,"[A-Z]")) 

# Substract the mean coefficient from the predicted values
# We will parallelise the code, so this will take multiple steps

# Define the  threats

threats <- unique(pops_data$threats[which(pops_data$n.threat>=1)]) %>% 
  as_tibble() %>% 
  mutate(number=str_count(value,"[A-Z]"))

# Define the number of cores to use

clus  <- makeCluster(detectCores() - 1)

# register the cluster with doParallel

registerDoSNOW(clus)

# load in the options using .options.snow
mu_scenarios <- foreach(i = threats$value, 
                        .combine = 'rbind',
                        .packages = c("tidyverse", "naniar")) %dopar% {
                          
                          x <- threats$number[threats$value==i]
                          
                          if(x==1){
                            new <- pop_perd %>% 
                              filter(grepl(i,threats)) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==i])) %>% 
                              rbind(pop_perd %>% 
                                      filter(!grepl(i,threats))) %>% 
                              mutate(scenario=paste("No", i))
                          } else if(x==2){
                            x2 <- unlist(strsplit(i, 
                                                  ":", 
                                                  fixed=T))
                            
                            new1 <- pop_perd %>% 
                              filter(grepl(x2[1],threats), 
                                     grepl(x2[2],threats), number!=1) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==i])) 
                            
                            new2 <- pop_perd %>% 
                              filter(grepl(x2[1],threats), 
                                     number!=2, !grepl(i,threats),
                                     !ID%in%new1$ID) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==x2[1]]))
                            
                            new <- pop_perd %>% 
                              filter(grepl(x2[2],threats), number!=2,
                                     !grepl(i,threats), 
                                     !ID%in%new1$ID,
                                     !ID%in%new2$ID) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==x2[2]])) %>% 
                              rbind(new1) %>% 
                              rbind(new2)  
                            
                            new <- new %>% 
                              rbind(pop_perd %>% 
                                      filter(!ID%in%new$ID)) %>% 
                              mutate(scenario=paste("No", i))
                          }else if(x==3){
                            
                            # Separate the threats 
                            x3 <- unlist(strsplit(i, 
                                                  ":", 
                                                  fixed=T))
                            
                            # Create new names
                            
                            two1 <- paste0(x3[1], ":", x3[2]) 
                            two2 <- paste0(x3[2], ":", x3[3]) 
                            two3 <- paste0(x3[1], ":", x3[3])
                            
                            # Create a data frame with the three threats
                            new1 <- pop_perd %>% 
                              filter(grepl(x3[1],threats), 
                                     grepl(x3[2],threats),
                                     grepl(x3[3],threats), 
                                     number==3) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==i])) 
                            
                            # Subset the mus with two threats
                            two_threat1 <- mus_mean %>%  filter(grepl(x3[1], .variable),
                                                                grepl(x3[2], .variable),
                                                                number==2)
                            two_threat2 <- mus_mean %>%  filter(grepl(x3[2], .variable),
                                                                grepl(x3[3], .variable),
                                                                number==2)
                            two_threat3 <- mus_mean %>%  filter(grepl(x3[1], .variable),
                                                                grepl(x3[3], .variable),
                                                                number==2)
                            
                            # subset two threats and substract
                            if(dim(two_threat1)[1]==0){
                              lam <- mus_mean$.value[mus_mean$.variable==x3[1]]+mus_mean$.value[mus_mean$.variable==x3[2]]
                              new2 <- tryCatch({pop_perd %>% 
                                  filter(grepl(x3[1], threats),
                                         grepl(x3[2], threats),
                                         !ID%in%new1$ID) %>% 
                                  mutate(lambda=exp(Estimate-lam))},
                                  error=function(e){e <- new1$lambda=NA})
                            }else{new2 <- tryCatch({pop_perd %>% 
                                filter(grepl(x3[1], threats),
                                       grepl(x3[2], threats),
                                       !ID%in%new1$ID) %>% 
                                mutate(lambda=exp(Estimate-two_threat1$.value))},
                                error=function(e){e <- new1$lambda=NA})}
                            if(dim(two_threat2)[1]==0){
                              lam <- mus_mean$.value[mus_mean$.variable==x3[2]]+mus_mean$.value[mus_mean$.variable==x3[3]]
                              new3 <- tryCatch({pop_perd %>% 
                                  filter(grepl(x3[2], threats),
                                         grepl(x3[3], threats),
                                         !ID%in%new1$ID) %>% 
                                  mutate(lambda=exp(Estimate-lam))},
                                  error=function(e){e <- new1$lambda=NA})
                            }else{
                              new3 <- tryCatch({pop_perd %>% 
                                  filter(grepl(x3[2], threats),
                                         grepl(x3[3], threats),
                                         !ID%in%new1$ID) %>% 
                                  mutate(lambda=exp(Estimate-two_threat2$.value))},
                                  error=function(e){e <- new1$lambda=NA})  
                            }
                            if(dim(two_threat3)[1]==0){ 
                              lam <- mus_mean$.value[mus_mean$.variable==x3[1]]+mus_mean$.value[mus_mean$.variable==x3[3]]
                              new4 <- tryCatch({pop_perd %>% 
                                  filter(grepl(x3[1], threats),
                                         grepl(x3[3], threats),
                                         !ID%in%new1$ID) %>% 
                                  mutate(lambda=exp(Estimate-lam))},
                                  error=function(e){e <- new1$lambda=NA})
                            }else{
                              new4 <- tryCatch({pop_perd %>% 
                                  filter(grepl(x3[1], threats),
                                         grepl(x3[3], threats),
                                         !ID%in%new1$ID) %>% 
                                  mutate(lambda=exp(Estimate-two_threat3$.value))},
                                  error=function(e){e <- new1$lambda=NA})
                            }
                            new1 %>% replace_with_na()
                            # Remove single effects
                            
                            new5 <- pop_perd %>% 
                              filter(grepl(x3[1], threats), 
                                     !ID%in%new1$ID, 
                                     !ID%in%new2$ID,
                                     !ID%in%new4$ID) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==x3[1]]))
                            
                            
                            new6 <- pop_perd %>% 
                              filter(grepl(x3[2], threats),
                                     !ID%in%new1$ID, 
                                     !ID%in%new2$ID,
                                     !ID%in%new3$ID) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==x3[2]]))
                            
                            # Final subset
                            
                            new <- pop_perd %>% 
                              filter(grepl(x3[3], threats),
                                     !ID%in%new1$ID, 
                                     !ID%in%new3$ID,
                                     !ID%in%new4$ID) %>% 
                              mutate(lambda=exp(Estimate-mus_mean$.value[mus_mean$.variable==x3[3]])) %>% 
                              rbind(new1) %>% 
                              rbind(new2) %>% 
                              rbind(new3) %>% 
                              rbind(new4) %>% 
                              rbind(new5) %>% 
                              rbind(new6)
                            
                            
                            
                            new <- new %>% 
                              rbind(pop_perd %>% 
                                      filter(!ID%in%new$ID)) %>% 
                              mutate(scenario=paste("No", i)) %>% 
                              drop_na(lambda)
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
