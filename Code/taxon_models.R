require(brms)
require(tidyverse)
load("Data/LivingPlanetData2.RData")

dd_long2 <- dd_long %>%
  group_by(ID) %>%  
  # Calculate population change
  mutate(popchange=log(Count+(max(Count, na.rm=T) /100)),
         #ifelse(Count<1, log(Count/(max(Count,na.rm = T))*100+1), 
         #        log(Count+1)),
         Protected_status=gsub(" .*", "", Protected_status)) %>% 
  # Remove any groupings we have created in the pipe
  ungroup()

# Add the biogeographical region

dd_long2 <- dd_long2 %>%
  mutate(T_realm=na_if(T_realm, "NULL"),
         M_realm=na_if(M_realm, "NULL"),
         FW_realm=na_if(FW_realm, "NULL"),
         Realm=coalesce(T_realm, M_realm, FW_realm))

# Find the data ids

dd_final <- dd_long2 %>%
  group_by(ID) %>%
  drop_na(popchange) %>%
  mutate(start=min(Year),
         end=max(Year),
         Duration=end-start) %>% 
  filter(Duration>9)

# Spread the data again

pops <- dd_final %>%
  mutate(year = as.factor(as.character(Year))) %>% 
  dplyr::select(ID, SpeciesName, Class, Order,
                System, Latitude, Longitude,
                Region,
                Protected_status, n.threat, 
                Primary_threat,
                Secondary_threat,
                Tertiary_threat,
                Managed,
                Realm,
                threats, Duration, 
                Year, year, popchange) %>%
  drop_na(popchange) %>%
  group_by(ID) %>%
  dplyr::select(-Year) %>%
  filter(length(unique(year))>=5) %>% 
  pivot_wider(names_from = year, values_from = popchange)

load("Data/SSResults.RData")

pop_data10 <- pop_data %>% 
  mutate(Year=as.numeric(Year)) %>% 
  group_by(ID) %>% 
  filter(length(unique(Year))>=10) %>%
  select(ID,pollution,habitatl,climatechange,
         invasive, exploitation,disease) %>%
  group_by(ID) %>%
  mutate(none = ifelse(all(is.na(c(pollution,habitatl,climatechange,
                                   invasive, exploitation,disease))),
                       "none",NA)) %>%
  filter(ID %in% pops$ID) %>%
  distinct()

mod_dat_full <- pops %>%
  pivot_longer(c("1950":"2019"), names_to = "year", values_to = "y") %>%
  mutate(year = as.numeric(year),
         #time = seq_along(year),
         series = paste(ID), #factor required for autocorrelation estimation
         threats = factor(ifelse(is.na(threats),"none",threats))) %>% #convert stressors to binary
  mutate(threats = fct_relevel(threats, "none")) %>%
  #filter(!all(is.na(value))) %>%
  left_join(pop_data10,by = "ID") %>%
  mutate(across(pollution:none,~ifelse(is.na(.x),"0","1"))) %>% #convert absence of stress to binary
  group_by(ID) %>%
  slice(1:max(which(!is.na(y))))  %>% #remove lagging years containing all NAs (to prevent unbounded interpolation)
  mutate(scaled_time = scales::rescale(year), #rescale start/end of time series between [0-1] for comparability
         time = seq_along(year),
         y_centered = y-(na.omit(y)[1]-0)) %>% #set first value of timeseries to 0 to allow no intercept model
  ungroup(ID) %>%
  as.data.frame()

test_data2 <-  mod_dat_full %>% 
  drop_na(y) %>%
  #filter(series %in% base::sample(unique(series),500)) %>%
  mutate(threats = as.character(threats)) %>%
  group_by(ID) %>%
  mutate(centered_year = year-min(year-1), #set first year with data to year 1
         scaled_time = scales::rescale(year)) %>%
  ungroup() %>%
  mutate(Taxon= ifelse(Class=="Holocephali"|Class=="Elasmobranchii" | 
                                   Class=="Myxini"|Class=="Cephalaspidomorphi"|
                                   Class=="Actinopterygii"|Class=="Sarcopterygii",
                                 "Fish", 
                                 ifelse(Class=="Aves", "Birds",
                                        ifelse(Class=="Mammalia", 
                                               "Mammals",
                                               ifelse(Class=="Amphibia", 
                                                      "Amphibians",
                                                      ifelse(Class=="Reptilia",
                                                             "Reptiles", "NA"))))))

##########################################################################################
## simple slopes ##
##########################################################################################

mod2_am <- brm(bf(y_centered ~  scaled_time*none + 
                     scaled_time*pollution + scaled_time*habitatl+
                     scaled_time*climatechange + scaled_time*invasive+ 
                     scaled_time*exploitation+scaled_time*disease + 
                     (-1 + time|series) + (1|SpeciesName) - 1
                   ,autocor = ~ar(time = time,gr = series,p=1)
                  ),
               data =test_data2 %>% filter(Taxon=="Amphibians"), 
               family = gaussian(),
               iter = 2000,
               refresh=100,
               #backend = "cmdstanr",
               chains = 4,
               control=list(adapt_delta=0.975,max_treedepth = 15),
               cores = 4)
saveRDS(mod2_am,"Results/models/simple-slopes-amph.RDS")

mod2_bir <- brm(bf(y_centered ~  scaled_time*none + 
                    scaled_time*pollution + scaled_time*habitatl+
                    scaled_time*climatechange + scaled_time*invasive+ 
                    scaled_time*exploitation+scaled_time*disease + 
                    (-1 + time|series) + (1|SpeciesName) - 1
                  ,autocor = ~ar(time = time,gr = series,p=1)
                  ),
                data =test_data2 %>% filter(Taxon=="Birds"), 
                family = gaussian(),
                iter = 2000,
                refresh=100,
                #backend = "cmdstanr",
                chains = 4,
                control=list(adapt_delta=0.975,max_treedepth = 15),
                cores = 4)
saveRDS(mod2_bir,"Results/models/simple-slopes-birds.RDS")

mod2_mam <- brm(bf(y_centered ~  scaled_time*none + 
                     scaled_time*pollution + scaled_time*habitatl+
                     scaled_time*climatechange + scaled_time*invasive+ 
                     scaled_time*exploitation+scaled_time*disease + 
                     (-1 + time|series) + (1|SpeciesName) - 1
                   ,autocor = ~ar(time = time,gr = series,p=1)
                   ),
                data =test_data2 %>% filter(Taxon=="Mammals"), 
                family = gaussian(),
                iter = 2000,
                refresh=100,
                #backend = "cmdstanr",
                chains = 4,
                control=list(adapt_delta=0.975,max_treedepth = 15),
                cores = 4)
saveRDS(mod2_mam,"Results/models/simple-slopes-mammals.RDS")

mod2_rep <- brm(bf(y_centered ~  scaled_time*none + 
                     scaled_time*pollution + scaled_time*habitatl+
                     scaled_time*climatechange + scaled_time*invasive+ 
                     scaled_time*exploitation+scaled_time*disease + 
                     (-1 + time|series) + (1|SpeciesName) - 1
                   ,autocor = ~ar(time = time,gr = series,p=1)
                   ),
                data =test_data2 %>% filter(Taxon=="Reptiles"), 
                family = gaussian(),
                iter = 2000,
                refresh=100,
                #backend = "cmdstanr",
                chains = 4,
                control=list(adapt_delta=0.975,max_treedepth = 15),
                cores = 4)
saveRDS(mod2_rep,"Results/models/simple-slopes-rept.RDS")
