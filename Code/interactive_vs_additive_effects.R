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
thrts_3 <-combn(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),3) #create all unique three way combinations of threats

test_data2 <-  mod_dat_full %>% 
  drop_na(y) #double check how to deal with infrequently sampled timeseries. Current model is able to deal via the ar process
#filter(series %in% base::sample(unique(series),500)) %>%

intvadd_dat2 <- test_data2 %>% #create a new dataframe with columns containing "0"/"1" for each combination of threats. "0" = combination not present, "1" = combination present
  bind_cols(purrr::pmap_dfc(.l = list(.x = thrts_2[1,], #threat 1
                                      .y = thrts_2[2,]),#threat 2
                            ~ test_data2 %>%
                              select({{.x}},{{.y}},ID) %>% #select just the necessary columns (threat 1, threat 2 and timeseries ID)
                              mutate(!!sym(paste(.x,.y,sep=".")) :=  
                                       ifelse(all(!!sym(.x) == "1") & all(!!sym(.y) == "1"),"1","0"),
                                     .by = ID) %>% #dynamically create new column of whether both threat 1 and threat 2 are 1's
                              select(4))%>%
              select_if(function(x)sum(x == "1") > 0) #drop columns where no observed combination of threats as will prevent model fitting
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

# intvadd_dat3 <- subset(intvadd_dat2,ID %in% sample(ID,100)) %>%
#   select_if(function(x)!all(x == "0"))

intvadd_dat3 <-intvadd_dat2 %>%
  group_by(ID) %>% 
  filter(length(unique(year))>=5) %>% #filter to timeseries > 5 time points
  ungroup() %>%
  select_if(function(x)!all(x == "0"))

rhs <- paste0(paste("scaled_year*",
                    colnames(intvadd_dat3)[
                      grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                                  collapse = "|"),
                            colnames(intvadd_dat3))],
                    sep ="",collapse = " + "),
              " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) + -1") #create the right hand side of the model formula with interactions between scaled_year and all threat combinations
form <- as.formula(paste("y_centered", "~", rhs)) #combine y_centered and rhs into a model formula

mod_intvadd2 <- brm(bf(form #include realm/spp as slopes, x intercepts
                       ,autocor = ~ar(time = time,gr = series,p=1) #ou/arima process
                       ),
                    data = intvadd_dat3, 
                    family = gaussian(),
                    iter = 8000,
                    refresh=100,
                    #backend = "cmdstanr",
                    prior = priors,
                    chains = 4,
                    control=list(adapt_delta=0.975,max_treedepth = 10),
                    cores = 4)
saveRDS(mod_intvadd2,"Results/models/simple-slopes-intvadd3.RDS")
mod_intvadd2 <- readRDS("Results/models/simple-slopes-intvadd2.RDS")
mod_intvadd3 <- readRDS("Results/models/simple-slopes-intvadd3.RDS")

source("Code/prep_data_grid_fn.R")
source("Code/threat_post_draws.R")

all_threats = c("pollution","habitatl","climatechange","invasive", "exploitation","disease")

threatcols <- colnames(mod_intvadd3$data)[grepl(paste(all_threats,collapse = "|"),colnames(mod_intvadd3$data))] 

additive_cols <- do.call("c",lapply(strsplit(threatcols,"[.]"),function(x){
  paste(x,collapse = " + ")
})) #create addtive columns. i.e. "threat1 + threat2"

target_cols <- unique(c(threatcols,additive_cols))

postdraws_intvadd <- threat_post_draws(model = mod_intvadd3,
                                       threat_comns = target_cols,
                                       n.cores = 4) %>%
  mutate(combo_group = case_when(
    grepl("[.]",threat) ~ "interactive",
    grepl("\\+",threat) ~ "additive",
    TRUE ~ "single")) #posterior timeseries estimated for each threat combination: singular (e.g. "threat1"), interactive (e.g. "threat1.threat2") and additive (e.g. "threat1 + threat2")

post_dydx_intvadd <- do.call("rbind",lapply(all_threats,function(x){
  
  out <- postdraws_intvadd %>%
    subset(grepl(x,threat)) %>% #filter to focal threat
    reframe(.value = mean(diff(.value)/diff(time)),.by = c(combo_group,threat,.draw)) %>% #estimate each timeseries' first derivative
    group_by(threat) %>%
    #filter(!any(.value >= abs(0.5))) %>% #drop highly variable threats
    ungroup() %>%
    mutate(threat_group = x) 
  
  return(out)
})) 

dydx_interval_intvadd <- post_dydx_intvadd  %>%
  #mutate(.chain = 1, .iteration = .draw) %>% #add additional columns required by ggdist
  group_by(threat_group,combo_group,threat) %>%
  ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw"))  #extract distribution information

ggplot(data = post_dydx_intvadd , 
       aes(x = .value,y=threat,fill = combo_group)) +
  tidybayes::stat_slab(alpha=0.5) +
  ggdist::geom_pointinterval(data = dydx_interval_intvadd,
                             aes(xmin = .lower, xmax = .upper,col = combo_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  scale_colour_manual(values = c("#F5A433","#8B33F5","#57A06B"), guide = "none") + 
  scale_fill_manual(values = c("#F5A433","#8B33F5","#57A06B"), guide = "none") + 
  #facet_wrap(~threat_group) + 
  coord_cartesian(xlim = c(-0.5,0.5))+
  xlab("Population trend") + 
  ylab("Threat") + 
  theme_minimal()

post_intvadd_diff <- do.call("rbind",lapply(additive_cols[grepl("\\+",additive_cols)],function(x){
  
  post_dydx_intvadd %>%
    subset(threat %in% c(x,gsub(" \\+ ",".",x))) %>% #subset to shared additive and interactive threats (e.g. "threat1.threat2" and "threat1 + threat2")
    #reframe(.value = .value[2]-.value[1], .by = c(threat_group,.draw)) %>% #find difference in derivatives between additive and interactive threats
    reframe(.value = diff(c(.value[2],.value[1])), .by = c(threat_group,.draw)) %>% #find difference in derivatives between additive and interactive threats
    mutate(threats = gsub(" \\+ ",".",x) ) #name threat combination
  
  # post_dydx_intvadd %>%
  #   subset(threat %in% c(x,gsub(" \\+ ",".",x))) %>% #subset to shared additive and interactive threats (e.g. "threat1.threat2" and "threat1 + threat2")
  #   group_by(threat_group,.draw) %>%
  #   summarise(first_threat = threat[1],
  #          second_threat = threat[2],
  #          .value = .value[2]-.value[1]) %>%
  #   ungroup() %>%
  #   mutate(threats = gsub(" \\+ ",".",x) ) #name threat combination
  
}))

post_interval_intvadd_diff <- post_intvadd_diff %>%
  na.omit() %>% #drop missing threats 
  group_by(threat_group,threats) %>%
  ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw")) %>%  #extract distribution information
  mutate(interaction.type = case_when(
    .upper < 0 ~ "synergistic",
    .lower > 0 ~ "antagonistic",
    TRUE ~ "additive"
  ))

ggplot(data = na.omit(post_intvadd_diff), 
       aes(x = .value,y=threats)) +
  tidybayes::stat_slab(data = na.omit(post_intvadd_diff) %>%
                         merge(select(post_interval_intvadd_diff,-.value),
                               by = c("threat_group","threats")) %>%
                         subset(.width == 0.8)
                       ,aes(fill = interaction.type),alpha=0.5) +
  ggdist::geom_pointinterval(data = post_interval_intvadd_diff,
                             aes(xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  coord_cartesian(xlim = c(-0.25,0.25)) + 
  scale_fill_manual(values = c("#694364",
                               "#1E63B3",
                               "#B32315"), name = "") + 
  #facet_wrap(~threat_group) + 
  xlab( expression(paste("Additive ",partialdiff,"y","/",partialdiff,"x"," - interactive ",partialdiff,"y","/",partialdiff,"x"))) + 
  ylab("Threat combination") + 
  theme_minimal()


post_intvadd_diff_syns <- post_intvadd_diff  %>%
  na.omit() %>% #drop missing threats 
  reframe(prop_great_zero = sum(.value>0)/n(),
          prop_less_zero = sum(.value<0)/n(),
          zero_quantile = ecdf(.value)(0),
          .by = c(threat_group,threats)) %>%
  mutate(interaction.type = case_when(
    zero_quantile > 0.9 ~ "synergistic",
    zero_quantile < 0.1 ~ "antagonistic",
    TRUE ~ "additive"
  )) %>%
  group_by(threat_group,interaction.type) %>%
  summarise(n = n()) %>%
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  ungroup()

ggplot(post_intvadd_diff_syns) +
  # Add the stacked bar
  geom_bar(aes(x=as.factor(threat_group), y=freq, fill=interaction.type), 
           stat="identity", alpha=0.8) +
  scale_fill_manual(name = "",
                    values = c("#694364",
                               "#B32315",
                               "#1E63B3"))+
  # Add text showing the freq of each 100/75/50/25 lines
  annotate("text", x = rep(max(post_intvadd_diff_syns$threat_group),5), y = c(20,40, 60, 80,100), 
           label = c("20", "40", "60", "80","100"), 
           color="grey", size=5, angle=0, fontface="bold", hjust=1.5) +
  ylim(-100,220) +
  #theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() 
