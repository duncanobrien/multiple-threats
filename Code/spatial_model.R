require(brms)
require(tidyverse)
require(geosphere)
require(patchwork)
require(ggdist)

source("Code/prep_data_grid_fn.R")
source("Code/threat_post_draws.R")

load("Data/LivingPlanetData2.RData")

norm_range <- function(x){
  (x-min(x))/(max(x)-min(x))
}

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

intvadd_dat3 <-intvadd_dat2 %>%
  group_by(ID) %>% 
  filter(length(unique(year))>=5) %>% #filter to timeseries > 5 time points
  ungroup() %>%
  select_if(function(x)!all(x == "0")) %>%
  mutate(Latitude = round(Latitude)) %>%
  mutate(Site = paste(Latitude,Longitude,sep = "_"))

locs <- distinct(intvadd_dat3[,c("Longitude","Latitude")]) %>%
  mutate(Site = paste(Latitude,Longitude,sep = "_"))

spa_mat_trim = as.matrix(geosphere::distm(locs[,c("Longitude","Latitude")], fun = geosphere::distHaversine))/1000 #km distance between sites
spa_mat_trim = norm_range(spa_mat_trim) #Before this step you could also exponentiate the distances, at the moment we are assuming a linear decay, which is probably fine.
spa_mat_trim = abs(spa_mat_trim - 1)
colnames(spa_mat_trim) = locs$Site
rownames(spa_mat_trim) = locs$Site

rhs <- paste0(paste("scaled_year*",
                    colnames(intvadd_dat3)[
                      grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                                  collapse = "|"),
                            colnames(intvadd_dat3))],
                    sep ="",collapse = " + "),
              " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) +
               (0 + scaled_year|gr(Site, cov = spa_mat)) + -1") #create the right hand side of the model formula with interactions between scaled_year and all threat combinations
form <- as.formula(paste("y_centered", "~", rhs)) #combine y_centered and rhs into a model formula

mod_intvadd_spatial <- brm(bf(form #include realm/spp as slopes, x intercepts
                       ,autocor = ~ar(time = time,gr = series,p=1) #ou/arima process
                       ),
                    data = intvadd_dat3, 
                    data2 = list(spa_mat = spa_mat_trim),
                    family = gaussian(),
                    iter = 4000,
                    refresh=100,
                    #backend = "cmdstanr",
                    prior = priors,
                    chains = 4,
                    control=list(adapt_delta=0.975,max_treedepth = 12),
                    cores = 4)

saveRDS(mod_intvadd_spatial,"Results/models/mod_intvadd_spatial.RDS")

mod_intvadd_spatial <- readRDS("Results/models/mod_intvadd_spatial.RDS")

####################
# combined effects
####################
all_threats = c("none","pollution","habitatl","climatechange","invasive", "exploitation","disease")

threatcols <- colnames(mod_intvadd_spatial$data)[grepl(paste(all_threats,collapse = "|"),colnames(mod_intvadd_spatial$data))] 

postdraws <- threat_post_draws(model = mod_intvadd_spatial,
                               ndraws = 1000,
                               threat_comns = c("none",threatcols)) #estimate posterior draws for all threat singles and combinations
post_dydx <- do.call("rbind",lapply(all_threats,function(x){

  out <- postdraws %>%
    subset(grepl(x,threat)) %>% #subset to focal threat
    reframe(.value = mean(diff(.value)/diff(time)),.by = c(threat,.draw)) %>% #for each posterior timeseries, estimate the first derivative 
    ungroup() %>%
    mutate(int_group = ifelse(grepl("\\.",threat),"combined","single"))%>% #create grouping column for whether the threat is singular or a combination
    mutate(threat_group = x) %>% #overall threat grouping
    group_by(threat) %>%
    #filter(!any(.value >= abs(0.5))) %>% #drop threats with highly uncertain confidence intervals. Typically certain three way interactions
    ungroup() 
  # %>%
  #   subset(!grepl("[[:lower:]]*\\.[[:lower:]]*\\.",threat)) 
  
  return(out)
})) %>%
  subset(threat != "pollution.climatechange")

dydx_interval <- post_dydx  %>%
  #subset(!grepl("[[:lower:]]*\\.[[:lower:]]*\\.",threat)) %>%
  group_by(threat_group,int_group) %>%
  ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw", "threat")) %>% #extract distribution information
  mutate(int_group = ifelse(int_group == "single","Singular","Interactive"))

threat_palette<-c(MetBrewer::met.brewer(name="Hokusai1", n=6, type="continuous"))

palatte <- data.frame(threat_group = unique(subset(post_dydx,threat != "none")$threat_group),
                      fill_col = threat_palette)

plot_dydx_threats <- post_dydx %>%
  subset(threat != "none") %>%
  left_join(palatte,by = "threat_group") %>%
  mutate(fill_col = factor(fill_col),
         int_group = ifelse(int_group == "single","Singular","Interactive"))

p1 <- ggplot(data = plot_dydx_threats, 
             aes(x = .value,y=int_group)) +
  tidybayes::stat_slab(data = subset(plot_dydx_threats,int_group == "Singular"),
                       aes(fill = fill_col,group = threat), alpha=0.3,normalize = "groups") +
  tidybayes::stat_slab(data = subset(plot_dydx_threats,int_group == "Interactive"),
                       aes(fill = fill_col,group = threat,col = fill_col), alpha=0.2,normalize = "panels") +
  ggdist::geom_pointinterval(data = subset(dydx_interval,threat_group != "none"),
                             aes(xmin = .lower, xmax = .upper),position = position_dodge()) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  xlab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
  ylab("Threat combination") + 
  coord_cartesian(xlim = c(-0.2,0.2)) +
  scale_x_continuous(breaks= seq(-0.1,0.1,by=0.1)) + 
  facet_wrap(~threat_group) +
  scale_fill_manual(values = levels(plot_dydx_threats$fill_col),guide = "none") + 
  scale_color_manual(values = levels(plot_dydx_threats$fill_col),guide = "none") + 
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

post_dydx_none <- postdraws %>%
  subset(grepl("none",threat)) %>% #subset to focal threat
  reframe(.value = mean(diff(.value)/diff(time)),.by = c(threat,.draw)) %>% #for each posterior timeseries, estimate the first derivative 
  mutate(int_group = ifelse(grepl("\\.",threat),"combined","single")) #create grouping column for whether the threat is singular or a combination

dydx_none_interval <- post_dydx_none  %>%
  group_by(int_group) %>%
  ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw", "threat")) #extract distribution information

p2 <- ggplot(data = post_dydx_none, 
             aes(x = .value,y=int_group)) +
  tidybayes::stat_slab(aes(), alpha=0.3) +
  ggdist::geom_pointinterval(data = dydx_none_interval,
                             aes(xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  xlab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
  #coord_cartesian(xlim = c(-0.1,0.1)) +
  facet_wrap(~threat) +
  theme_minimal()+
  theme(axis.title.x = element_text(size=12,
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size = 12),
        strip.text.x = element_text(size = 12),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks = element_line(color="black"),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0,0,0,0), "pt"))

p1 + (p2/patchwork::plot_spacer() + plot_layout(heights  = c(1,1))) + plot_layout(widths = c(2,0.8)) 

ggsave("Results/Figure2_full_spatial.pdf",last_plot(),
       height = 8,width = 12,dpi = 300)

####################
# interactive vs additive
####################
all_threats = c("pollution","habitatl","climatechange","invasive", "exploitation","disease")

threatcols <- colnames(mod_intvadd_spatial$data)[grepl(paste(all_threats,collapse = "|"),colnames(mod_intvadd_spatial$data))] 

additive_cols <- do.call("c",lapply(strsplit(threatcols,"[.]"),function(x){
  paste(x,collapse = " + ")
})) #create addtive columns. i.e. "threat1 + threat2"

target_cols <- unique(c(threatcols,additive_cols))

postdraws_intvadd <- threat_post_draws(model = mod_intvadd_spatial,
                                       threat_comns = target_cols,
                                       nuisance = c("series","SpeciesName"),
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
  group_by(threat_group,combo_group,threat) %>%
  ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw"))  #extract distribution information

post_intvadd_diff <- do.call("rbind",lapply(additive_cols[grepl("\\+",additive_cols)],function(x){
  
  post_dydx_intvadd %>%
    subset(threat %in% c(x,gsub(" \\+ ",".",x))) %>% #subset to shared additive and interactive threats (e.g. "threat1.threat2" and "threat1 + threat2")
    #reframe(.value = .value[2]-.value[1], .by = c(threat_group,.draw)) %>% #find difference in derivatives between additive and interactive threats
    reframe(.value = diff(c(.value[2],.value[1])), .by = c(threat_group,.draw)) %>% #find difference in derivatives between additive and interactive threats
    mutate(threats = gsub(" \\+ ",".",x) ) #name threat combination
  
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
                       ,aes(fill = interaction.type),alpha=0.5,normalize = "xy") +
  ggdist::geom_pointinterval(data = post_interval_intvadd_diff,
                             aes(xmin = .lower, xmax = .upper,col = interaction.type),alpha=0.5) +
  geom_point(data = subset(post_interval_intvadd_diff,.width == 0.8), 
             aes(x = .value,col = interaction.type,fill = interaction.type),shape = 21,alpha=0.5,size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  coord_cartesian(xlim = c(-0.25,0.25)) + 
  scale_fill_manual(values = c("#694364",
                               "#1E63B3",
                               "#B32315"), name = "",
                    guide = guide_legend(override.aes = list(color = NA,shape = 2) )) + 
  scale_color_manual(values = c("#694364",
                                "#1E63B3",
                                "#B32315"), guide = "none") + 
  #facet_wrap(~threat_group) + 
  xlab( expression(paste("Additive ",partialdiff,"y","/",partialdiff,"x"," - interactive ",partialdiff,"y","/",partialdiff,"x"))) + 
  ylab("Threat combination") + 
  theme_minimal()

ggsave("Results/Figure3_full_spatial.pdf",last_plot(),
       height = 8,width = 8,dpi = 300)

# Conterfactual tests ----------------------------------------------------------
source("Code/threat_counterfac_draws.R")
# trends for each of the counterfactual scenarios 

threat_cols <- c("pollution","habitatl",
                 "climatechange","invasive", 
                 "exploitation","disease")

threat_cols <-  colnames(mod_intvadd_spatial$data)[grepl(paste(c("pollution","habitatl",
                                                "climatechange","invasive", 
                                                "exploitation","disease"),
                                              collapse = "|"),
                                        colnames(mod_intvadd_spatial$data))]

# Predict the population trends in the "no intervention" scenario

pop_perd <- brms::posterior_epred(mod_intvadd_spatial,
                                  newdata = mod_intvadd_spatial$data %>% 
                                    filter_at(threat_cols, any_vars(. != "0")),
                                  re.form = NA,
                                  incl_autocor = FALSE,
                                  sort = TRUE, 
                                  ndraws = 1000) %>% 
  as.data.frame() %>%
  mutate(.draw = 1:NROW(.)) %>%
  #extract posterior draws for the data used to create the model
  pivot_longer(-.draw,names_to = "index",values_to = ".value") %>%
  cbind(mod_intvadd_spatial$data %>%
          filter_at(threat_cols, any_vars(. != "0")) %>% 
          dplyr::select(series, time), row.names = NULL) %>% 
  reframe(.value = mean(diff(.value)/diff(time)),.by = c(series,.draw)) %>% 
  mutate(counterfac="none")

# Create the different counterfactual scenarios

# We use the function counterfactual draws to estimate the different population 

counter_fac_data <- threat_counterfac_draws(mod_intvadd_spatial,
                                            threat_comns = threat_cols,
                                            re.form = NA,
                                            ndraws = 1000,
                                            center_dydx = "mean",
                                            n.cores = 6) %>%
  # Join with the none counterfactual scenario that we just created
  rbind(pop_perd) %>%
  # Made "none" as the first level of the counterfactual 
  mutate(counterfac = fct_relevel(counterfac, "none"))

threat_palette<-c(MetBrewer::met.brewer(name="Hokusai1", n=6, type="continuous"), "grey60")

# Lets have a look at it 

scenarios_mean <- counter_fac_data %>% 
  group_by(counterfac) %>% 
  summarise(m=median(.value)) %>% 
  arrange(desc(m))

# One threat

(p1 <- counter_fac_data %>% 
    # add a variable counting the number of threats
    mutate(number=str_count(counterfac, '\\.')+1) %>%
    group_by(counterfac) %>% 
    mutate(pos=median(.value)) %>% 
    ungroup() %>%   
    filter(number==1) %>% 
    # Start the plot
    ggplot(aes(x = .value,
               y=reorder(counterfac, pos),
               fill=counterfac, 
               colour=counterfac)) +
    stat_halfeye(adjust = .5, 
                 width = .6, 
                 .width = 0,
                 alpha=.8,
                 justification = -.3, 
                 point_colour = NA,
                 normalize = "groups") + 
    geom_boxplot(width = .2, 
                 outlier.shape = NA, 
                 alpha=.8,
                 colour="black") +
    scale_fill_manual(values=c(threat_palette))+
    scale_colour_manual(values=threat_palette)+
    # geom_boxplot(outlier.shape = NA)+
    # tidybayes::stat_halfeye(aes(group = counterfac)) +
    geom_vline(xintercept = 0,
               linetype = "solid", 
               colour="grey50") +
    geom_vline(xintercept = subset(scenarios_mean,counterfac == "none")$m,
               linetype = "dashed", 
               colour="grey50") +
    xlab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
    ylab("Counterfactual") +
    coord_flip(xlim = c(-0.2,0.2))+
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
          plot.title = element_text(hjust = 0.5),
          legend.position = "none"))

ggsave("Results/counterfactuals_spatial.pdf", p1, 
       width = 8, height = 6)

# Two threats

(p2 <- counter_fac_data %>% 
    # add a variable counting the number of threats
    mutate(number=str_count(counterfac, '\\.')+1) %>%
    group_by(counterfac) %>% 
    mutate(pos=median(.value)) %>% 
    ungroup() %>%   
    filter(number==2) %>% 
    # Start the plot
    ggplot(aes(x = .value,
               y=reorder(counterfac, pos),
               fill=counterfac, 
               colour=counterfac)) +
    stat_halfeye(adjust = .5, 
                 width = .6, 
                 .width = 0,
                 alpha=.8,
                 justification = -.3, 
                 point_colour = NA,
                 normalize = "groups") + 
    geom_boxplot(width = .2, 
                 outlier.shape = NA, 
                 alpha=.8,
                 colour="black") +
    scale_fill_manual(values=MetBrewer::met.brewer(name="Tiepolo", n=15, type="continuous"))+
    scale_colour_manual(values=MetBrewer::met.brewer(name="Tiepolo", n=15, type="continuous"))+
    # geom_boxplot(outlier.shape = NA)+
    # tidybayes::stat_halfeye(aes(group = counterfac)) +
    geom_vline(xintercept = 0,
               linetype = "solid", 
               colour="grey50") +
    geom_vline(xintercept = subset(scenarios_mean,counterfac == "none")$m,
               linetype = "dashed", 
               colour="grey50") +
    xlab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
    ylab("Counterfactual") +
    coord_cartesian(xlim = c(-0.2,0.2))+
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
          plot.title = element_text(hjust = 0.5),
          legend.position = "none"))

ggsave("Results/counterfactualstwothreats_spatial.pdf", p2, 
       width = 6, height = 12)

# Three threats

(p3 <- counter_fac_data %>% 
    # add a variable counting the number of threats
    mutate(number=str_count(counterfac, '\\.')+1) %>%
    group_by(counterfac) %>% 
    mutate(pos=median(.value)) %>% 
    ungroup() %>%   
    filter(number==3) %>% 
    # Start the plot
    ggplot(aes(x = .value,
               y=reorder(counterfac, pos),
               fill=counterfac, 
               colour=counterfac)) +
    stat_halfeye(adjust = .5, 
                 width = .6, 
                 .width = 0,
                 alpha=.8,
                 justification = -.3, 
                 point_colour = NA,
                 normalize = "xy") + 
    geom_boxplot(width = .2, 
                 outlier.shape = NA, 
                 alpha=.8,
                 colour="black") +
    scale_fill_manual(values=MetBrewer::met.brewer(name="Nattier", n=18, type="continuous"))+
    scale_colour_manual(values=MetBrewer::met.brewer(name="Nattier", n=18, type="continuous"))+
    # geom_boxplot(outlier.shape = NA)+
    # tidybayes::stat_halfeye(aes(group = counterfac)) +
    geom_vline(xintercept = 0,
               linetype = "solid", 
               colour="grey50") +
    geom_vline(xintercept = subset(scenarios_mean,counterfac == "none")$m,
               linetype = "dashed", 
               colour="grey50") +
    xlab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
    ylab("Counterfactual") +
    coord_cartesian(xlim = c(-0.2,0.2))+
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
          plot.title = element_text(hjust = 0.5),
          legend.position = "none"))

ggsave("Results/counterfactualsthree_spatial.pdf", p3, 
       width = 6, height = 14)
