# --------------------------------------------------------------------------------------- #
# - FILE NAME:   figure_3.R         
# - DATE:        11/02/2024
# - DESCRIPTION: Code to plot the figures 2 and 3 from the main manuscrpit.
# - AUTHORS: Pol Capdevila Lanzaco (pcapdevila@ub.edu), Duncan O'brien (duncan.a.obrien@gmail.com)
# --------------------------------------------------------------------------------------- #


rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(tidyverse)
library(brms)
library(tidybayes)
library(bayestestR)
library(MetBrewer)
library(wesanderson)
library(cowplot)

# Set default ggplot theme

theme_set(theme_minimal()+
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
                  plot.title = element_text(hjust = 0.5)))

#Create a palette for the threats 

threat_palette<- c(met.brewer(name="Hokusai1", n=6,type="continuous"))

# Color palette for taxons

taxon_pal <- wes_palette("Cavalcanti1", n = 5)

# System palette

system_pal <- c("#A1D6E2","#336B87", "#CC954E")

# Load other customised functions to substract the trends from the plots

source("Code/prep_data_grid_fn.R")
source("Code/threat_post_draws.R")

# Load the models 

mod_amph <- readRDS("Results/models/amph_mod.RDS")
mod_bir <- readRDS("Results/models/bir_mod.RDS")
mod_mam <- readRDS("Results/models/mam_mod.RDS")
mod_rep <- readRDS("Results/models/rep_mod.RDS")
mod_mar <- readRDS("Results/models/mar_mod.RDS")
mod_fre <- readRDS("Results/models/fre_mod.RDS")
mod_ter <- readRDS("Results/models/ter_mod.RDS")

# Define the threats

all_threats <- c("none","pollution","habitatl","climatechange",
                 "invasive", "exploitation","disease")

# Prepare data for system ------------------------------------------------------

system_model_ls <- list("Marine" = mod_mar,
                        "Freshwater" = mod_fre,
                        "Terrestrial" = mod_ter)

# Classify calculate population trends and whether they are additive, synergistic
# or antagonistic

post_mod_system <- lapply(seq_along(system_model_ls), FUN = function(X){
  
  # Extract column names from the data within the list of models
  threatcols <- colnames(system_model_ls[[X]]$data)[grepl(paste(all_threats, collapse = "|"), colnames(system_model_ls[[X]]$data))]
  
  # Create additive columns by combining the column names
  additive_cols <- do.call("c", lapply(strsplit(threatcols, "[.]"), function(x){
    paste(x, collapse = " + ")
  })) 
  
  # Combine unique threat columns and additive columns
  target_cols <- unique(c(threatcols, additive_cols))
  
  # Generate posterior draws for the system model with specified threats and columns
  postdraws_system <- threat_post_draws(model = system_model_ls[[X]],
                                        threat_comns = target_cols,
                                        ndraws = 1000,
                                        nuisance = c("series", "SpeciesName"),
                                        n.cores = 4) %>%
    mutate(combo_group = case_when(
      grepl("[.]", threat) ~ "interactive", # If threat contains a dot, classify as interactive
      grepl("\\+", threat) ~ "additive",    # If threat contains a plus sign, classify as additive
      TRUE ~ "single"))                     # Otherwise, classify as singular
  
})

# Substract the differences

post_mod_system_diff <- do.call("rbind",lapply(post_mod_system , `[[`, 'dydx_diff'))

# Substract the intervals

post_intvadd_interval_realm_diff <- do.call("rbind",lapply(post_mod_system , `[[`, 'interval'))

# Visualise ----------------------------------------------------------

ggplot(data = plot_dydx_realm, 
       aes(x = .value,y=int_group)) +
  tidybayes::stat_slab(data = subset(plot_dydx_realm,int_group == "Singular"),
                       aes(fill = fill_col,group = threat), alpha=0.3,normalize = "groups") +
  tidybayes::stat_slab(data = subset(plot_dydx_realm,int_group == "Interactive"),
                       aes(fill = fill_col,group = threat,col = fill_col), alpha=0.2,normalize = "panels") +
  ggdist::geom_pointinterval(data = subset(post_dydx_interval_realm),
                             aes(xmin = .lower, xmax = .upper),position = position_dodge()) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  xlab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
  ylab("Threat combination") + 
  coord_cartesian(xlim = c(-0.2,0.2)) +
  scale_x_continuous(breaks= seq(-0.1,0.1,by=0.1)) + 
  facet_wrap(~threat_group+System,ncol = 3,scales = "free_y") +
  scale_fill_manual(values = levels(plot_dydx_realm$fill_col),guide = "none") + 
  scale_color_manual(values = levels(plot_dydx_realm$fill_col),guide = "none") + 
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

ggsave("Results/Figure2_system_spatial.pdf",last_plot(),
       height = 10,width = 12,dpi = 300)

ggplot(data = na.omit(post_intvadd_realm_diff), 
       aes(x = .value,y=threats)) +
  ggdist::stat_slab(data = na.omit(post_intvadd_realm_diff) %>%
                      merge(select(post_intvadd_interval_realm_diff,-.value),
                            by = c("System","threat_group","threats")) %>%
                      subset(.width == 0.8)
                    ,aes(fill = interaction.type),alpha=0.5,normalize = "xy") +
  ggdist::geom_pointinterval(data = post_intvadd_interval_realm_diff,
                             aes(xmin = .lower, xmax = .upper,col = interaction.type),alpha=0.5) +
  geom_point(data = subset(post_intvadd_interval_realm_diff,.width == 0.8), 
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
  facet_wrap(~System,ncol = 3) + 
  xlab( expression(paste("Additive ",partialdiff,"y","/",partialdiff,"x"," - interactive ",partialdiff,"y","/",partialdiff,"x"))) + 
  ylab("Threat combination") + 
  theme_minimal()

ggsave("Results/Figure3_system_spatial.pdf",last_plot(),
       height = 10,width = 12,dpi = 300)


# Create a data frame by combining the results for each threat in all_threats

post_dydx_system <- do.call("rbind", lapply(all_threats, function(x){
  
  # Subset the posterior draws to include only the current threat
  out <- postdraws_system %>%
    subset(grepl(x, threat)) %>% 
    # Reframe the data to estimate the first derivative of each timeseries
    reframe(.value = mean(diff(.value) / diff(time)), 
            .by = c(combo_group, threat, .draw)) %>% 
    group_by(threat) %>%
    # Ungroup the data
    ungroup() %>%
    # Add a column indicating the threat group
    mutate(threat_group = x) 
  
  return(out) # Return the modified data
}))

  
  post_system_diff <- do.call("rbind",lapply(additive_cols[grepl("\\+",additive_cols)],function(x){
    
    post_dydx_system %>%
      subset(threat %in% c(x,gsub(" \\+ ",".",x))) %>% #subset to shared additive and interactive threats (e.g. "threat1.threat2" and "threat1 + threat2")
      reframe(.value = diff(c(.value[2],.value[1])), .by = c(threat_group,.draw)) %>% #find difference in derivatives between additive and interactive threats
      mutate(threats = gsub(" \\+ ",".",x) ) #name threat combination
    
  })) %>%
    mutate(System = names(system_model_ls)[X])
  
  post_interval_system_diff <- post_system_diff %>%
    na.omit() %>% #drop missing threats 
    group_by(System,threat_group,threats) %>%
    ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw")) %>%  #extract distribution information
    mutate(interaction.type = case_when(
      .upper < 0 ~ "synergistic",
      .lower > 0 ~ "antagonistic",
      TRUE ~ "additive"
    )) 
  
  return(list("dydx_diff" = post_system_diff,"interval" = post_interval_system_diff))
})

post_mod_system_diff <- do.call("rbind",lapply(post_mod_system, `[[`, 'dydx_diff'))

post_mod_interval_system_diff <- do.call("rbind",lapply(post_mod_system, `[[`, 'interval'))




post_dydx_system <- lapply(seq_along(system_model_ls), FUN = function(X) {
  # Loop through each model in taxon_model_ls
  # Extract column names related to threats in the current model's data
  threatcols <- colnames(system_model_ls[[X]]$data)[grepl(paste(all_threats,
                                                                collapse = "|"),
                                                          colnames(system_model_ls[[X]]$data))] 
  
  # Calculate posterior draws for all combinations of threats
  
  postdraws <- threat_post_draws(model = system_model_ls[[X]],
                                 threat_comns = c("none",threatcols),
                                 ndraws = 1000,
                                 nuisance = c("series","SpeciesName"),
                                 n.cores = 4) 
  
  # Calculate the first derivative of population trends for each threat
  
  post_dydx <- do.call("rbind", lapply(all_threats, function(x) {
    out <- postdraws %>%
      subset(grepl(x, threat)) %>%
      reframe(.value = mean(diff(.value) / diff(time)), .by = c(threat, .draw)) %>%
      ungroup() %>%
      mutate(int_group = ifelse(grepl("\\.", threat), "combined", "single")) %>%
      mutate(threat_group = x) %>%
      group_by(threat) %>%
      ungroup() 
    
    return(out)
  })) %>%
    mutate(System = names(system_model_ls)[X])
  
  # Calculate median and credible intervals of the derivatives
  
  dydx_interval <- post_dydx %>%
    group_by(System, threat_group, int_group) %>%
    ggdist::median_qi(.width = c(.95, .8, .5), .exclude = c(".draw", "threat")) %>%
    mutate(int_group = ifelse(int_group == "single", "Singular", "Interactive"))
  
  # Return a list containing derivatives and their credible intervals
  return(list("dydx" = post_dydx, "interval" = dydx_interval))
})

# Reformat the data to have the trends and intervals in a single data frame

post_dyxd_system <- do.call("rbind",lapply(post_dydx_system, `[[`, 'dydx')) %>%
  mutate(threat_group = fct_relevel(factor(threat_group),"none"), 
         threat_group = fct_recode(threat_group, 
                                   Polluntion="pollution",
                                   `Habitat loss`="habitatl", 
                                   `Climate change`="climatechange", 
                                   `Invasive species`="invasive",
                                   Exploitation="exploitation",
                                   Disease="disease",
                                   None="none"))

post_dydx_interval_system <- do.call("rbind",lapply(post_dydx_system, `[[`, 'interval')) %>%
  mutate(threat_group = fct_relevel(factor(threat_group),"none"),
         threat_group = fct_recode(threat_group, 
                                   Polluntion="pollution",
                                   `Habitat loss`="habitatl", 
                                   `Climate change`="climatechange", 
                                   `Invasive species`="invasive",
                                   Exploitation="exploitation",
                                   Disease="disease",
                                   None="none"))

# Calculate the median population trend for each threat

pos <-  post_dyxd_system %>%
  filter(threat != "none", 
         int_group=="Singular") %>% 
  group_by(threat_group) %>% 
  summarise(pos=median(.value)) %>% 
  arrange(desc(pos)) %>%
  pull(threat_group)

# Classify each threat into single or interactive and associate each threat with
# their median

plot_dydx_system <- post_dyxd_system %>%
  mutate(int_group = ifelse(int_group == "single","Singular","Interactive"),
         threat_group=fct_relevel(threat_group, as.character(pos)),
         threat_group=fct_relevel(factor(threat_group),"None"))


ggplot(data = na.omit(post_mod_taxon_diff), 
       aes(x = .value,y=threats)) +
  ggdist::stat_slab(data = na.omit(post_mod_taxon_diff) %>%
                      merge(select(post_mod_interval_taxon_diff,-.value),
                            by = c("Taxon","threat_group","threats")) %>%
                      subset(.width == 0.8)
                    ,aes(fill = interaction.type),alpha=0.5,normalize = "xy") +
  ggdist::geom_pointinterval(data = post_mod_interval_taxon_diff,
                             aes(xmin = .lower, xmax = .upper,col = interaction.type),alpha=0.5) +
  geom_point(data = subset(post_mod_interval_taxon_diff,.width == 0.8), 
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
  facet_wrap(~Taxon,ncol = 4) + 
  xlab( expression(paste("Additive ",partialdiff,"y","/",partialdiff,"x"," - interactive ",partialdiff,"y","/",partialdiff,"x"))) + 
  ylab("Threat combination") + 
  theme_minimal()

ggsave("Results/Figure3_taxon_spatial.pdf",last_plot(),
       height = 10,width = 12,dpi = 300)

# Extract trends as derivatives for each system ----------------------------------------------------------
source("Code/prep_data_grid_fn.R")
source("Code/threat_post_draws.R")
mod_mar <- readRDS("Results/models/mod_mar_mod_spatial.RDS")

mod_fre <- readRDS("Results/models/mod_fre_mod_spatial.RDS")

mod_ter <- readRDS("Results/models/mod_ter_mod_spatial.RDS")

all_threats = c("none","pollution","habitatl","climatechange","invasive", "exploitation","disease")

realm_model_ls <- list("Marine" = mod_mar,"Freshwater" = mod_fre,"Terrestrial" = mod_ter)

post_dydx_realm <- lapply(seq_along(realm_model_ls), FUN = function(X){
  
  threatcols <- colnames(realm_model_ls[[X]]$data)[grepl(paste(all_threats,collapse = "|"),colnames(realm_model_ls[[X]]$data))] 
  postdraws <- threat_post_draws(model = realm_model_ls[[X]],
                                 threat_comns = c("none",threatcols),
                                 ndraws = 1000,
                                 nuisance = c("series","SpeciesName"),
                                 n.cores = 4) #estimate posterior draws for all threat singles and combinations
  
  post_dydx <- do.call("rbind",lapply(all_threats,function(x){
    out <- postdraws %>%
      subset(grepl(x,threat)) %>% #subset to focal threat
      reframe(.value = mean(diff(.value)/diff(time)),.by = c(threat,.draw)) %>% #for each posterior timeseries, estimate the first derivative 
      ungroup() %>%
      mutate(int_group = ifelse(grepl("\\.",threat),"combined","single"))%>% #create grouping column for whether the threat is singular or a combination
      mutate(threat_group = x) %>% #overall threat grouping
      group_by(threat) %>%
      ungroup() 
    return(out)
  })) %>%
    mutate(System = names(realm_model_ls)[X])
  
  dydx_interval <- post_dydx  %>%
    group_by(System,threat_group,int_group) %>%
    ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw", "threat")) %>% #extract distribution information
    mutate(int_group = ifelse(int_group == "single","Singular","Interactive"))
  
  return(list("dydx" = post_dydx,"interval" = dydx_interval))
})

post_dyxd_realm <- do.call("rbind",lapply(post_dydx_realm, `[[`, 'dydx')) %>%
  mutate(threat_group = fct_relevel(factor(threat_group),"none"))

post_dydx_interval_realm <- do.call("rbind",lapply(post_dydx_realm, `[[`, 'interval')) %>%
  mutate(threat_group = fct_relevel(factor(threat_group),"none"))
threat_palette<-c(MetBrewer::met.brewer(name="Hokusai1", n=6, type="continuous"))

palatte <- data.frame(threat_group = unique(subset(post_dyxd_realm,threat != "none")$threat_group),
                      fill_col = threat_palette)

plot_dydx_realm <- post_dyxd_realm %>%
  #subset(threat != "none") %>%
  left_join(palatte,by = "threat_group") %>%
  mutate(fill_col = factor(fill_col),
         int_group = ifelse(int_group == "single","Singular","Interactive")) 

# Extract additive vs interactive for each realm ----------------------------------------------------------
source("Code/prep_data_grid_fn.R")
source("Code/threat_post_draws.R")

mod_mar <- readRDS("Results/models/mod_mar_mod_spatial.RDS")

mod_fre <- readRDS("Results/models/mod_fre_mod_spatial.RDS")

mod_ter <- readRDS("Results/models/mod_ter_mod_spatial.RDS")

all_threats_intvadd = c("pollution","habitatl","climatechange","invasive", "exploitation","disease")

realm_model_ls <- list("Marine" = mod_mar,"Freshwater" = mod_fre,"Terrestrial" = mod_ter)


# Visualise ----------------------------------------------------------

ggplot(data = plot_dydx_realm, 
       aes(x = .value,y=int_group)) +
  tidybayes::stat_slab(data = subset(plot_dydx_realm,int_group == "Singular"),
                       aes(fill = fill_col,group = threat), alpha=0.3,normalize = "groups") +
  tidybayes::stat_slab(data = subset(plot_dydx_realm,int_group == "Interactive"),
                       aes(fill = fill_col,group = threat,col = fill_col), alpha=0.2,normalize = "panels") +
  ggdist::geom_pointinterval(data = subset(post_dydx_interval_realm),
                             aes(xmin = .lower, xmax = .upper),position = position_dodge()) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  xlab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
  ylab("Threat combination") + 
  coord_cartesian(xlim = c(-0.2,0.2)) +
  scale_x_continuous(breaks= seq(-0.1,0.1,by=0.1)) + 
  facet_wrap(~threat_group+System,ncol = 3,scales = "free_y") +
  scale_fill_manual(values = levels(plot_dydx_realm$fill_col),guide = "none") + 
  scale_color_manual(values = levels(plot_dydx_realm$fill_col),guide = "none") + 
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

ggsave("Results/Figure2_system_spatial.pdf",last_plot(),
       height = 10,width = 12,dpi = 300)

ggplot(data = na.omit(post_mod_realm_diff), 
       aes(x = .value,y=threats)) +
  ggdist::stat_slab(data = na.omit(post_mod_realm_diff) %>%
                      merge(select(post_mod_interval_realm_diff,-.value),
                            by = c("System","threat_group","threats")) %>%
                      subset(.width == 0.8)
                    ,aes(fill = interaction.type),alpha=0.5,normalize = "xy") +
  ggdist::geom_pointinterval(data = post_mod_interval_realm_diff,
                             aes(xmin = .lower, xmax = .upper,col = interaction.type),alpha=0.5) +
  geom_point(data = subset(post_mod_interval_realm_diff,.width == 0.8), 
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
  facet_wrap(~System,ncol = 3) + 
  xlab( expression(paste("Additive ",partialdiff,"y","/",partialdiff,"x"," - interactive ",partialdiff,"y","/",partialdiff,"x"))) + 
  ylab("Threat combination") + 
  theme_minimal()

ggsave("Results/Figure3_system_spatial.pdf",last_plot(),
       height = 10,width = 12,dpi = 300)


mod_mod_spatial <- readRDS("Results/models/mod_mod_spatial.RDS")

####################
# combined effects
####################
all_threats = c("none","pollution","habitatl","climatechange","invasive", "exploitation","disease")

threatcols <- colnames(mod_mod_spatial$data)[grepl(paste(all_threats,collapse = "|"),colnames(mod_mod_spatial$data))] 

postdraws <- threat_post_draws(model = mod_mod_spatial,
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

threatcols <- colnames(mod_mod_spatial$data)[grepl(paste(all_threats,collapse = "|"),colnames(mod_mod_spatial$data))] 

additive_cols <- do.call("c",lapply(strsplit(threatcols,"[.]"),function(x){
  paste(x,collapse = " + ")
})) #create addtive columns. i.e. "threat1 + threat2"

target_cols <- unique(c(threatcols,additive_cols))

postdraws_intvadd <- threat_post_draws(model = mod_mod_spatial,
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

post_mod_diff <- do.call("rbind",lapply(additive_cols[grepl("\\+",additive_cols)],function(x){
  
  post_dydx_intvadd %>%
    subset(threat %in% c(x,gsub(" \\+ ",".",x))) %>% #subset to shared additive and interactive threats (e.g. "threat1.threat2" and "threat1 + threat2")
    #reframe(.value = .value[2]-.value[1], .by = c(threat_group,.draw)) %>% #find difference in derivatives between additive and interactive threats
    reframe(.value = diff(c(.value[2],.value[1])), .by = c(threat_group,.draw)) %>% #find difference in derivatives between additive and interactive threats
    mutate(threats = gsub(" \\+ ",".",x) ) #name threat combination
  
}))

post_interval_mod_diff <- post_mod_diff %>%
  na.omit() %>% #drop missing threats 
  group_by(threat_group,threats) %>%
  ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw")) %>%  #extract distribution information
  mutate(interaction.type = case_when(
    .upper < 0 ~ "synergistic",
    .lower > 0 ~ "antagonistic",
    TRUE ~ "additive"
  ))

ggplot(data = na.omit(post_mod_diff), 
       aes(x = .value,y=threats)) +
  tidybayes::stat_slab(data = na.omit(post_mod_diff) %>%
                         merge(select(post_interval_mod_diff,-.value),
                               by = c("threat_group","threats")) %>%
                         subset(.width == 0.8)
                       ,aes(fill = interaction.type),alpha=0.5,normalize = "xy") +
  ggdist::geom_pointinterval(data = post_interval_mod_diff,
                             aes(xmin = .lower, xmax = .upper,col = interaction.type),alpha=0.5) +
  geom_point(data = subset(post_interval_mod_diff,.width == 0.8), 
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

threat_cols <-  colnames(mod_mod_spatial$data)[grepl(paste(c("pollution","habitatl",
                                                             "climatechange","invasive", 
                                                             "exploitation","disease"),
                                                           collapse = "|"),
                                                     colnames(mod_mod_spatial$data))]

# Predict the population trends in the "no intervention" scenario

pop_perd <- brms::posterior_epred(mod_mod_spatial,
                                  newdata = mod_mod_spatial$data %>% 
                                    filter_at(threat_cols, any_vars(. != "0")),
                                  re.form = NA,
                                  incl_autocor = FALSE,
                                  sort = TRUE, 
                                  ndraws = 1000) %>% 
  as.data.frame() %>%
  mutate(.draw = 1:NROW(.)) %>%
  #extract posterior draws for the data used to create the model
  pivot_longer(-.draw,names_to = "index",values_to = ".value") %>%
  cbind(mod_mod_spatial$data %>%
          filter_at(threat_cols, any_vars(. != "0")) %>% 
          dplyr::select(series, time), row.names = NULL) %>% 
  reframe(.value = mean(diff(.value)/diff(time)),.by = c(series,.draw)) %>% 
  mutate(counterfac="none")

# Create the different counterfactual scenarios

# We use the function counterfactual draws to estimate the different population 

counter_fac_data <- threat_counterfac_draws(mod_mod_spatial,
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


