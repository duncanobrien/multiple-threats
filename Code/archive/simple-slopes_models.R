require(brms)
require(tidyverse)
load("Data/LivingPlanetData2.RData")

priors <- c(prior(normal(0, 1), class = b),
            prior(exponential(1), class = sd),
            #prior(exponential(1), class = sigma),
            prior(normal(0,0.25), class = ar))

####################
# raw data taken from PolCap repository
####################
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

pop_data10 <- pop_data %>% # pop10 selected following original repository
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

####################
# further wrangling in to long format and prepare explanatory variables
####################
mod_dat_full <- pops %>%
  pivot_longer(c("1950":"2019"), names_to = "year", values_to = "y") %>%
  mutate(year = as.numeric(year),
         #time = seq_along(year),
         series = paste(ID), #factor required for autocorrelation estimation
         threats = factor(ifelse(is.na(threats),"none",threats))) %>% #convert stressors to binary
  mutate(threats = fct_relevel(threats, "none")) %>%
  left_join(select(pop_data,c(ID,pollution,habitatl,climatechange,
                              invasive, exploitation,disease)),
            multiple = "first",by = "ID") %>% #double check. Currently NOT selected timeseries >= 10 years
  mutate(none = ifelse(all(is.na(c(pollution,habitatl,climatechange,
                                   invasive, exploitation,disease))),
                       "none",NA)) %>%
  mutate(across(pollution:none,~ifelse(is.na(.x),"0","1"))) %>% #convert absence of stress to binary
  group_by(ID) %>%
  slice(1:max(which(!is.na(y))))  %>% #remove lagging years containing all NAs (to prevent unbounded interpolation)
  mutate(scaled_year = c(scale(year,center = TRUE,scale = FALSE)), #center time to 0 for each timeseries
         time = seq_along(year),
         scaled_time = c(scale(time,center = TRUE,scale = FALSE)),
         #y_centered = log(y/max(na.omit(y))), #rescale y by maximum of timeseries and log
         y_centered = y-(na.omit(y)[1]-0) #recenter y so that first value of timeseries is 0 (to allow all intercepts to be removed)
         ) %>% 
  ungroup(ID) %>%
  as.data.frame()

thrts_2 <-combn(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),2) #create all unique two way combinations of threats
thrts_3 <-combn(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),3) 

test_data2 <-  mod_dat_full %>% 
  drop_na(y) #double check how to deal with infrequently sampled timeseries. Current model is able to deal via the ar process
  #filter(series %in% base::sample(unique(series),500)) %>%

test_data3 <- test_data2 %>%
  bind_cols(purrr::pmap_dfc(.l = list(.x = thrts_2[1,], #threat 1
                                      .y = thrts_2[2,]),#threat 2
                            ~ test_data2 %>%
                              select({{.x}},{{.y}},ID) %>% #select just the necessary columns (threat 1, threat 2 and timeseries ID)
                              mutate(!!sym(paste(.x,.y,sep=".")) :=  
                                       ifelse(all(!!sym(.x) == "1") & all(!!sym(.y) == "1"),"1","0"),
                                     .by = ID) %>% #dynamically create new column of whether both threat 1 and threat 2 are 1's
                              select(4))) %>%
  bind_cols(purrr::pmap_dfc(.l = list(thrts_3[1,], #threat 1
                                      thrts_3[2,], #threat 2
                                      thrts_3[3,]),#threat 3
                            function(.x,.y,.z) test_data2 %>%
                              select({{.x}},{{.y}},{{.z}},ID) %>% #select just the necessary columns (threat 1, threat 2 and timeseries ID)
                              mutate(!!sym(paste(.x,.y,.z,sep=".")) :=  
                                       ifelse(all(!!sym(.x) == "1") & all(!!sym(.y) == "1") & all(!!sym(.z) == "1"),"1","0"),
                                     .by = ID) %>% #dynamically create new column of whether both threat 1, threat 2 and threat 3 are 1's
                              select(5)))

ggplot(subset(test_data2,ID %in% sample(ID,5)),aes(x=scaled_time,y=y_centered,col= factor(ID))) + geom_line()
#why negative y values

##########################################################################################
## simple slopes ##
##########################################################################################

mod2b <- brm(bf(y_centered ~  #divide by maximum then log, then center first value on 0
                  scaled_year*threats + #time centered on the mean time 
                 (-1 + scaled_year|series) +
                 (-1 + scaled_year|SpeciesName) + 
                  (-1 + scaled_year|Realm) 
                - 1 #include realm/spp as slopes, x intercepts
               ,autocor = ~ar(time = time,gr = series,p=1) #ou/arima process
               ),
             data =test_data3, 
             family = gaussian(),
             iter = 2000,
             refresh=100,
             #backend = "cmdstanr",
             chains = 4,
             control=list(adapt_delta=0.975,max_treedepth = 15),
             cores = 4)
conditional_effects(mod2b)

saveRDS(mod2b,"Results/models/simple-slopes_b.RDS")

mod2_draws <- emmeans::emtrends(mod2b, ~ threats,var = "scaled_year") |>
  tidybayes::gather_emmeans_draws() |>
  dplyr::rename(.variable = threats) |>
  dplyr::mutate(.variable =  gsub("Invasive spp/genes","Invasive species",.variable),
                .variable =  gsub("none","None",.variable)) |>
  dplyr::mutate(.variable = factor(.variable),
                .variable = fct_relevel(.variable,"None")) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

ggplot(mod2_draws,aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper),interval_size_range = c(0.8, 2)) +
  labs(x="Coefficient", y = "Threat") + 
  scale_y_discrete(limits = rev(levels(mod2_draws$.variable))) + 
   theme_minimal()+
  theme(axis.title.x = element_text(size=12,
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size=12,
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size = 12),
        axis.text.y = element_text(color="black", size = 12),
        strip.text.x = element_text(size = 12),
        axis.ticks = element_line(color="black"),
        plot.title = element_text(hjust = 0.5))
