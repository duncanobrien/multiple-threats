require(mvgam)
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

pops_data <- pops %>%
  dplyr::select(ID, SpeciesName, Class, Order,
                System, Latitude, Longitude,
                Region,
                Protected_status, n.threat, 
                Primary_threat,
                Secondary_threat,
                Tertiary_threat,
                Managed,
                Realm,
                threats, Duration) 

##########################################################################################
## MVGAM ##
##########################################################################################
sample_ts <- 1:10

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
  filter(ID %in% pops[sample_ts,]$ID) %>%
  distinct()


mod_dat <- pops[sample_ts,] %>%
  pivot_longer(c("1950":"2019"), names_to = "year", values_to = "y") %>%
  mutate(year = as.numeric(year),
         time = seq_along(year),
         series = factor(ID),
         threats = ifelse(is.na(threats),"none",threats)) %>%
  #filter(!all(is.na(value))) %>%
  left_join(pop_data10,by = "ID") %>%
  mutate(across(pollution:none,~factor(ifelse(is.na(.x),0,1)))) %>%
  as.data.frame()

trend_map <- data.frame(series = unique(mod_dat$series),
                        trend = 1:length(unique(mod_dat$series)))

fake_mod <- mvgam(y ~ year*threats ,
                  trend_formula = ~ year,
                  trend_model = 'RW',
                  #trend_map = trend_map,
                  family = gaussian(),
                  data = mod_dat,
                  backend = "cmdstanr",
                  chains = 4,
                  burnin = 500,
                  samples = 500)
unique(mod_dat$series)
plot(x = fake_mod, type = 'trend', series = 2)

