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

##########################################################################################
## brms ##
##########################################################################################
load("Data/SSResults.RData")

sample_ts <- 1:250
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

mod_dat_full <- pops %>%
  pivot_longer(c("1950":"2019"), names_to = "year", values_to = "y") %>%
  mutate(year = as.numeric(year),
         #time = seq_along(year),
         series = paste(ID), #factor required for autocorrelation estimation
         threats = factor(ifelse(is.na(threats),"none",threats))) %>% #convert stressors to binary
  mutate(threats = fct_relevel(threats, "none")) %>%
  #filter(!all(is.na(value))) %>%
  left_join(pop_data10,by = "ID") %>%
  mutate(across(pollution:none,~ifelse(is.na(.x),0,1))) %>% #convert absence of stress to binary
  group_by(ID) %>%
  slice(1:max(which(!is.na(y))))  %>% #remove lagging years containing all NAs (to prevent unbounded interpolation)
  mutate(scaled_time = scales::rescale(year), #rescale start/end of time series between [0-1] for comparability
         time = seq_along(year),
         y_centered = y-(na.omit(y)[1]-0)) %>% #set first value of timeseries to 0 to allow no intercept model
  ungroup(ID) %>%
  as.data.frame()

mod_dat <- pops[sample_ts,] %>%
  pivot_longer(c("1950":"2019"), names_to = "year", values_to = "y") %>%
  mutate(year = as.numeric(year),
         #time = seq_along(year),
         series = paste(ID), #factor required for autocorrelation estimation
         threats = factor(ifelse(is.na(threats),"none",threats))) %>% #convert stressors to binary
  mutate(threats = fct_relevel(threats, "none")) %>%
  #filter(!all(is.na(value))) %>%
  left_join(pop_data10,by = "ID") %>%
  mutate(across(pollution:none,~factor(ifelse(is.na(.x),0,1)))) %>% #convert absence of stress to binary
  group_by(ID) %>%
  slice(1:max(which(!is.na(y))))  %>% #remove lagging years containing all NAs (to prevent unbounded interpolation)
  mutate(scaled_time = scales::rescale(year), #rescale start/end of time series between [0-1] for comparability
         time = seq_along(year),
         y_centered = y-(na.omit(y)[1]-0)) %>% #set first value of timeseries to 0 to allow no intercept model
  ungroup(ID) %>%
  as.data.frame()

ggplot(mod_dat,aes(x=year,y=y,col=series)) + 
  geom_point() + geom_smooth(method = "lm") + facet_wrap(~series)

test_data2 <-  mod_dat_full %>% 
  drop_na(y) %>%
  #filter(series %in% base::sample(unique(series),500)) %>%
  mutate(threats = as.character(threats)) %>%
  group_by(ID) %>%
  mutate(centered_year = year-min(year),
         scaled_time = scales::rescale(year)) %>%
  ungroup()

mod1 <- brm(bf(y_centered ~ trend  -1,
       trend ~ scaled_time:threats + (-1 + time|series) - 1, 
       nl = TRUE
       , autocor = ~ar(time = time,gr = series,p=1)
       ),
    data = test_data2, family = gaussian(),
    iter = 2000,
    refresh=100,
    #backend = "cmdstanr",
    control=list(adapt_delta=0.985,max_treedepth = 15),
    cores = 4)

mod1_coefs <- broom.mixed::tidy(mod1) |>
  filter(grepl("trend_scaled",term)) |>
  mutate(term = gsub("trend_scaled_time:threats","",term))

ggplot(as.data.frame(mod1_coefs) |>
         arrange(term),
       aes(y = term, x = estimate)) + 
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) + 
  geom_vline(xintercept=0) + 
  ggtitle("non-linear") + 
  theme_bw()

write_csv(mod1_coefs,"/Users/ul20791/Desktop/mod1.csv")

mod1_alt <- brm(bf(y_centered ~ scaled_time*trend  -1,
               trend ~ threats  - 1, nl = TRUE
               #, autocor = ~ar(time = time,gr = series,p=1)
               ),
               data = test_data2, family = gaussian(),
               iter = 2000,
               refresh=100,
               #backend = "cmdstanr",
               control=list(adapt_delta=0.985,max_treedepth = 15),
               cores = 4)

mod3 <- brm(bf(y ~ trend -1,
               trend ~ year:none + 
                 year:pollution + year:habitatl+
                 year:climatechange + year:invasive+ 
                 year:exploitation + year:disease + 
                 (0+year|series) - 1, nl = TRUE
               , autocor = ~ar(time = time,gr = series,p=1)
               ),
            data = mod_dat_full %>%  drop_na(y) %>% filter(System=="Marine") %>% filter(series %in% base::sample(unique(series),50)), 
            family = gaussian(),
            iter = 2000,
            refresh=100,
            #backend = "cmdstanr",
            control=list(adapt_delta=0.975,max_treedepth = 15),
            cores = 4)

ranef(mod1)
tt <- conditional_effects(mod3)


ggplot(test_data2,aes(x=scaled_time,y=y_centered,col=series)) + 
  geom_point() + geom_smooth(method = "lm") + facet_wrap(~series)

mod2 <- brm(bf(y_centered ~  scaled_time*threats + (-1 + time|series) - 1
               ,autocor = ~ar(time = time,gr = series,p=1)
               ),
            data =test_data2, 
            family = gaussian(),
            iter = 2000,
            refresh=100,
            #backend = "cmdstanr",
            chains = 4,
            control=list(adapt_delta=0.975,max_treedepth = 15),
            cores = 4)
simp_slopes <- emmeans::emtrends(mod2, ~ threats,var = "scaled_time")
plot(simp_slopes)
conditional_effects(mod2)
mod2_coefs <- as.data.frame(simp_slopes)

ggplot(as.data.frame(mod2_coefs) |>
         arrange(threats),
       aes(y = threats, x = scaled_time.trend)) + 
  geom_pointrange(aes(xmin = lower.HPD, xmax = upper.HPD)) + 
  geom_vline(xintercept=0) + 
  ggtitle("state space") + 
  theme_bw()

mod2b <- brm(bf(y_centered ~  scaled_time*threats + (-1 + time|series) +
                 (1|SpeciesName) + (1|Realm) - 1
               ,autocor = ~ar(time = time,gr = series,p=1)
               ),
            data =test_data2, 
            family = gaussian(),
            iter = 2000,
            refresh=100,
            #backend = "cmdstanr",
            chains = 4,
            control=list(adapt_delta=0.975,max_treedepth = 15),
            cores = 4)
simp_slopes_b <- emmeans::emtrends(mod2b, ~ threats,var = "scaled_time")
plot(simp_slopes_b)
conditional_effects(mod2b)
mod2b_coefs <- as.data.frame(simp_slopes_b)

ggplot(as.data.frame(mod2b_coefs) |>
         arrange(threats),
       aes(y = threats, x = scaled_time.trend)) + 
  geom_pointrange(aes(xmin = lower.HPD, xmax = upper.HPD)) + 
  geom_vline(xintercept=0) + 
  ggtitle("state space") + 
  theme_bw()


ss_data10 <- pops_data %>% 
  filter(ID%in%test_data2$ID) %>%
  mutate(threats = factor(ifelse(is.na(threats),"none",threats))) %>% #convert stressors to binary
  mutate(threats = fct_relevel(threats, "none"))

mod5 <- brm(mu | se(sigma, sigma = TRUE)  ~ threats-1 , 
            data = ss_data10,
            family = gaussian(), 
            iter = 2000,
            refresh =100,
            chains = 4,
            cores = 4,
            control = list(adapt_delta = .975, max_treedepth = 15))


mod5_coefs <- broom.mixed::tidy(mod5) |>
  filter(grepl("threats",term)) |>
  mutate(term = gsub("threats","",term))

ggplot(as.data.frame(mod5_coefs) |>
         arrange(term),
       aes(y = term, x = estimate)) + 
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) + 
  geom_vline(xintercept=0) + 
  ggtitle("state space") + 
  theme_bw()
