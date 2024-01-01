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
         y_centered = log(y/max(na.omit(y))), #rescale y by maximum of timeseries and log
         y_centered = y-(na.omit(y)[1]-0) #recenter y so that first value of timeseries is 0 (to allow all intercepts to be removed)
  ) %>% 
  ungroup(ID) %>%
  as.data.frame()

test_data3 <-  mod_dat_full %>% 
  drop_na(y) %>%
  #filter(series %in% base::sample(unique(series),500)) %>%
  mutate(threats = as.character(threats)) 

ggplot(subset(test_data3,ID %in% sample(ID,5)),aes(x=scaled_time,y=y_centered,col= factor(ID))) + geom_line()

#why negative y values

##########################################################################################
## Duncan first attempt interactive vs additive ##
##########################################################################################
intvadd_data <- test_data3 %>%
  dplyr::mutate(n_threats = case_when(
                  n.threat == 0 ~ "Additive",
                  n.threat == 1 ~ "Additive",
                  TRUE ~ "Interactive"
                ),
                across(pollution:none,~as.factor(.x)))
                

mod_intvadd <- brm(bf(y_centered ~  #divide by maximum then log, then center first value on 0
                        scaled_year*habitatl*n_threats +
                        scaled_year*climatechange*n_threats +
                        scaled_year*invasive*n_threats +
                        scaled_year*exploitation*n_threats +
                        scaled_year*disease*n_threats + #time centered on the mean time 
                  (-1 + scaled_time|series) +
                  #(-1 + scaled_time|SpeciesName) + 
                    #(-1 + scaled_time|Realm) 
                  - 1 #include realm/spp as slopes, x intercepts
                ,autocor = ~ar(time = time,gr = series,p=1) #ou/arima process
                ),
             data = subset(intvadd_data,ID %in% sample(ID,150)), 
             family = gaussian(),
             iter = 2000,
             refresh=100,
             #backend = "cmdstanr",
             chains = 4,
             control=list(adapt_delta=0.975,max_treedepth = 15),
             cores = 4)

saveRDS(mod_intvadd,"Results/models/simple-slopes-intvadd.RDS")


##########################################################################################
## Tom suggestion interactive vs additive ##
##########################################################################################
thrts_2 <-combn(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),2) #create all unique two way combinations of threats
thrts_3 <-combn(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),3) 

test_data2 <-  mod_dat_full %>% 
  drop_na(y) #double check how to deal with infrequently sampled timeseries. Current model is able to deal via the ar process
#filter(series %in% base::sample(unique(series),500)) %>%

intvadd_dat2 <- test_data2 %>%
  bind_cols(purrr::pmap_dfc(.l = list(.x = thrts_2[1,], #threat 1
                                      .y = thrts_2[2,]),#threat 2
                            ~ test_data2 %>%
                              select({{.x}},{{.y}},ID) %>% #select just the necessary columns (threat 1, threat 2 and timeseries ID)
                              mutate(!!sym(paste(.x,.y,sep=".")) :=  
                                       ifelse(all(!!sym(.x) == "1") & all(!!sym(.y) == "1"),"1","0"),
                                     .by = ID) %>% #dynamically create new column of whether both threat 1 and threat 2 are 1's
                              select(4))%>%
              select_if(function(x)sum(x == "1") > 0) #drop columns where no observed combination of threats
            ) %>%
  bind_cols(purrr::pmap_dfc(.l = list(thrts_3[1,], #threat 1
                                      thrts_3[2,], #threat 2
                                      thrts_3[3,]),#threat 3
                            function(.x,.y,.z) test_data2 %>%
                              select({{.x}},{{.y}},{{.z}},ID) %>% #select just the necessary columns (threat 1, threat 2 and timeseries ID)
                              mutate(!!sym(paste(.x,.y,.z,sep=".")) :=  
                                       ifelse(all(!!sym(.x) == "1") & all(!!sym(.y) == "1") & all(!!sym(.z) == "1"),"1","0"),
                                     .by = ID) %>% #dynamically create new column of whether both threat 1, threat 2 and threat 3 are 1's
                              select(5)) %>%
              select_if(function(x)sum(x == "1") > 0) #drop columns where no observed combination of threats
            ) 

intvadd_dat3 <- subset(intvadd_dat2,ID %in% sample(ID,100)) %>%
  select_if(function(x)!all(x == "0"))

intvadd_dat3 <-intvadd_dat2 %>%
  group_by(ID) %>% 
  filter(length(unique(year))>=10) %>%
  ungroup() %>%
  select_if(function(x)!all(x == "0"))
  
rhs <- paste0(paste("scaled_year*",colnames(intvadd_dat3)[c(21:26,32:NCOL(intvadd_dat3))],sep ="",collapse = " + "),
             " + (-1 + scaled_year|series) + -1")
form <- as.formula(paste("y_centered", "~", rhs))

mod_intvadd2 <- brm(bf(form #include realm/spp as slopes, x intercepts
                      ,autocor = ~ar(time = time,gr = series,p=1) #ou/arima process
                      ),
                   data = intvadd_dat3, 
                   family = gaussian(),
                   iter = 8000,
                   refresh=100,
                   backend = "cmdstanr",
                   prior = priors,
                   chains = 4,
                   control=list(adapt_delta=0.975,max_treedepth = 10),
                   cores = 4)
saveRDS(mod_intvadd2,"Results/models/simple-slopes-intvadd2.RDS")

mod_intvadd2_draws <- emmeans::emtrends(mod_intvadd2,var = "scaled_year"
                                        #,rg.limit = 55000
                                        ) |>
  tidybayes::gather_emmeans_draws() 

threatcols <- colnames(intvadd_dat3)[c(21:26,32:NCOL(intvadd_dat3))]
newdata <- data.frame(scaled_year = seq(min(intvadd_dat3$scaled_year),max(intvadd_dat3$scaled_year),by=1),
                      series = NA)
selcols <- threatcols[grepl("pollution",threatcols)]
nullcols <- threatcols[!grepl("pollution",threatcols)]
newdata <- cbind(newdata, setNames(lapply(selcols, function(x) x=factor("1",levels = c("0","1"))), selcols))
newdata <- cbind(newdata, setNames(lapply(nullcols, function(x) x=factor("0",levels = c("0","1"))), nullcols)) %>%
  mutate(time = seq_along(scaled_year))

postdraws_pollution <- posterior_epred(mod_intvadd2,
                      newdata =  newdata) %>%
  as.data.frame() %>%
  pivot_longer(everything(),names_to = ".draw",values_to = ".value",
               names_transform = function(x){as.numeric(gsub("V","",x))}) %>%
  cbind(newdata)
  


