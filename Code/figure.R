# --------------------------------------------------------------------------------------- #
# - FILE NAME:   figures.R         
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
library(patchwork)

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

taxon_pal <- c(wes_palette("Cavalcanti1", n = 5))

# System palette

system_pal <- c("#A1D6E2","#336B87", "#CC954E")

# Load other customised functions to substract the trends from the plots

source("Code/prep_data_grid_fn.R")
source("Code/threat_post_draws.R")

# Load the models 

mod_amph <- readRDS("Results/models/amph_mod.RDS")
mod_bir <- readRDS("Results/models/bir_mod.RDS")
mod_fish <- readRDS("Results/models/fish_mod.RDS")
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


# Prepare data for taxa --------------------------------------------------------

# Create a list with all the models 

taxon_model_ls <- list("Amphibians" = mod_amph,
                       "Birds" = mod_bir,
                       "Fish" = mod_fish,
                       "Mammals" = mod_mam,
                       "Reptiles" = mod_rep)

# Calculate the population trends for each model

post_dydx_taxon <- lapply(seq_along(taxon_model_ls), FUN = function(X) {
  # Loop through each model in taxon_model_ls
  # Extract column names related to threats in the current model's data
  threatcols <- colnames(taxon_model_ls[[X]]$data)[grepl(paste(all_threats,
                                                               collapse = "|"),
                                                         colnames(taxon_model_ls[[X]]$data))] 
  
  # Calculate posterior draws for all combinations of threats
  
  postdraws <- threat_post_draws(model = taxon_model_ls[[X]],
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
  mutate(Taxon = names(taxon_model_ls)[X])
  
  # Calculate median and credible intervals of the derivatives

  dydx_interval <- post_dydx %>%
  group_by(Taxon, threat_group, int_group) %>%
  ggdist::median_qi(.width = c(.95, .8, .5), .exclude = c(".draw", "threat")) %>%
  mutate(int_group = ifelse(int_group == "single", "Singular", "Interactive"))

# Return a list containing derivatives and their credible intervals
return(list("dydx" = post_dydx, "interval" = dydx_interval))
})

# Reformat the data to have the trends and intervals in a single data frame

post_dyxd_taxon <- do.call("rbind",lapply(post_dydx_taxon, `[[`, 'dydx')) %>%
  mutate(threat_group = fct_relevel(factor(threat_group),"none"), 
         threat_group = fct_recode(threat_group, 
                                   Polluntion="pollution",
                                   `Habitat loss`="habitatl", 
                                   `Climate change`="climatechange", 
                                   `Invasive species`="invasive",
                                   Exploitation="exploitation",
                                   Disease="disease",
                                   None="none"))

post_dydx_interval_taxon <- do.call("rbind",lapply(post_dydx_taxon, `[[`, 'interval')) %>%
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

pos <-  post_dyxd_taxon %>%
  filter(threat != "none", 
         int_group=="Singular") %>% 
  group_by(threat_group) %>% 
  summarise(pos=median(.value)) %>% 
  arrange(desc(pos)) %>%
  pull(threat_group)
  
# Classify each threat into single or interactive and associate each threat with
# their median

plot_dydx_taxon <- post_dyxd_taxon %>%
  mutate(int_group = ifelse(int_group == "single","Singular","Interactive"),
         threat_group=fct_relevel(threat_group, as.character(pos)),
         threat_group=fct_relevel(factor(threat_group),"None"))

# Figure 2: Single and multiple stressors by system and taxa ###################
# Panel a: Single stressors by system ------------------------------------------

(g2a <- plot_dydx_system %>% 
   filter(int_group=="Singular") %>% 
   ggplot(aes(x = threat_group, y=.value)) +
   stat_pointinterval(aes(colour = System, 
                          group = System), 
                      .width = c(.65, .80, .95),
                      position =position_dodge(.7))+
   geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
   ylab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
   xlab("") + 
   coord_cartesian(ylim = c(-0.5,0.5)) +
   scale_color_manual(values = system_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Table S1 ---------------------------------------------------------------------

TableS1 <- as.data.frame(describe_posterior(mod_ter, ci = 0.95, test="none")) %>% 
  mutate(System="Terrestrial") %>% 
  rbind(as.data.frame(describe_posterior(mod_fre, ci = 0.95, test="none")) %>% 
  mutate(System="Freshwater")) %>% 
  rbind(as.data.frame(describe_posterior(mod_mar, ci = 0.95, test="none")) %>% 
  mutate(System="Marine")) %>%  
  mutate(Parameter = gsub("b_", "", Parameter),
         Parameter = gsub("invasive", "Invasive", Parameter),
         Parameter = gsub("habitatl", "Habitat loss", Parameter),
         Parameter = gsub("climatechange", "Climate change", Parameter),
         Parameter = gsub("pollution", "Pollution", Parameter),
         Parameter = gsub("exploitation", "Exploitation", Parameter),
         Parameter = gsub("disease", "Disease", Parameter),
         Parameter = gsub("scaled_year", "Year", Parameter)) 


# Save it 

write.csv2(TableS1, "Results/Table S1.csv")

# Panel b: Single stressors by taxon -------------------------------------------

(g2b <- plot_dydx_taxon %>% 
  filter(int_group=="Singular") %>% 
  ggplot(aes(x = threat_group, y=.value)) +
  tidybayes::stat_pointinterval(aes(colour = Taxon, group = Taxon), 
                                .width = c(.65, .80, .95),position =position_dodge(.7))+
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  ylab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
  xlab("") + 
  coord_cartesian(ylim = c(-0.5,0.5)) +
  scale_color_manual(values = taxon_pal)+
  theme(plot.tag = element_text(size = 14, face = 'bold'),
        plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
        plot.title = element_text(size = 14, hjust = -0.2),
        axis.text.x = element_text(angle=25, hjust=1)))

# Table S2 ---------------------------------------------------------------------

TableS2 <- as.data.frame(describe_posterior(mod_amph, ci = 0.95, test="none")) %>% 
  mutate(Taxon="Amphibians") %>% 
  rbind(as.data.frame(describe_posterior(mod_amph, ci = 0.95, test="none")) %>% 
  mutate(Taxon="Birds")) %>% 
  rbind(as.data.frame(describe_posterior(mod_amph, ci = 0.95, test="none")) %>% 
  mutate(Taxon="Mammals")) %>% 
  rbind(as.data.frame(describe_posterior(mod_amph, ci = 0.95, test="none")) %>% 
  mutate(Taxon="Reptiles")) %>% 
  mutate(Parameter = gsub("b_", "", Parameter),
         Parameter = gsub("invasive", "Invasive", Parameter),
         Parameter = gsub("habitatl", "Habitat loss", Parameter),
         Parameter = gsub("climatechange", "Climate change", Parameter),
         Parameter = gsub("pollution", "Pollution", Parameter),
         Parameter = gsub("exploitation", "Exploitation", Parameter),
         Parameter = gsub("disease", "Disease", Parameter),
         Parameter = gsub("scaled_year", "Year", Parameter)) 


# Save it 

write.csv2(TableS2, "Results/TableS2.csv")

# Panel c: Multiple stressors by system ----------------------------------------

(g2c <- plot_dydx_system %>% 
   filter(int_group=="Interactive"|
            threat_group=="None") %>% 
   ggplot(aes(x = threat_group, y=.value)) +
   stat_pointinterval(aes(colour = System, 
                          group = System), 
                      .width = c(.65, .80, .95),
                      position =position_dodge(.7))+
   geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
   ylab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
   xlab("") + 
   coord_cartesian(ylim = c(-0.5,0.5)) +
   scale_color_manual(values = system_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Panel d: Multiple stressors by taxon -----------------------------------------

(g2d <- plot_dydx_taxon %>% 
   filter(int_group=="Interactive"|
            threat_group=="None") %>% 
   ggplot(aes(x = threat_group, y=.value)) +
   stat_pointinterval(aes(colour = Taxon, 
                          group = Taxon), 
                      .width = c(.65, .80, .95),
                      position =position_dodge(.7))+
   geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
   ylab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
   xlab("") + 
   coord_cartesian(ylim = c(-0.5,0.5)) +
   scale_fill_manual(values = taxon_pal)+
   scale_color_manual(values = taxon_pal)+
   theme(plot.tag = element_text(size = 14, face = 'bold'),
         plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
         plot.title = element_text(size = 14, hjust = -0.2),
         axis.text.x = element_text(angle=25, hjust=1)))

# Combine ----------------------------------------------------------------------

# Create title

title <- ggdraw() + 
  draw_label("Single threats",
             fontface = 'bold',size = 18,
             x = 0,
             hjust = 0) 

# Create subtitles

subtitle1 <- ggdraw() + 
  draw_label("System",
             size = 16,
             x = 0.5,
             vjust = 0) 

subtitle2 <- ggdraw() + 
  draw_label("Taxon",
             x = 0.5,
             size = 16,
             vjust = 0) 
subtitles <- plot_grid(subtitle1, subtitle2, nrow = 1)

# Create row 1

(row1 <- plot_grid(g2a+theme(legend.position = "none",
                             axis.title.x = element_blank()),
                   g2b+theme(legend.position = "none",
                             axis.title.x = element_blank()),
                   labels = "auto"))
# Create title 2

title2 <- ggdraw() + 
  draw_label("Interacting threats",
             fontface = 'bold',
             x = 0, size = 18,
             hjust = 0) 

# Create row 2

(row2 <- plot_grid(g2c+theme(legend.position = "none",
                             axis.title.x = element_blank()),
                   g2d+theme(legend.position = "none",
                             axis.title.x = element_blank()),
                   labels = c("c","d")))

# Greate two legends 

legends <- plot_grid(get_legend(g2a+theme(legend.position = "bottom",
                                          legend.text = element_text(size = 14))),
                     get_legend(g2b+theme(legend.position = "bottom",
                                          legend.text = element_text(size = 14))))

# Put them together 

(figure2<- plot_grid(title, subtitles, row1, title2, row2, legends,
                     rel_heights = c(0.1,0.1,1,.1,1,.1), ncol=1))

# Save it

ggsave("Results/Figure 2.pdf",figure2,
       width = 14, height = 10)

# Figure 3: Interactive threats ################################################
# Classify calculate population trends and whether they are additive, synergistic
# or antagonistic

# Define a function to operate on each element of the system_model_ls list
post_intvadd_system <- lapply(seq_along(system_model_ls), FUN = function(X){
  
  # Extract column names from the data within the current element of the list
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
  
  # Calculate derivatives for each threat
  post_dydx_system <- do.call("rbind", lapply(all_threats, function(x){
    
    out <- postdraws_system %>%
      subset(grepl(x, threat)) %>% # Filter to focal threat
      reframe(.value = mean(diff(.value) / diff(time)), .by = c(combo_group, threat, .draw)) %>% # Estimate each timeseries' first derivative
      group_by(threat) %>%
      ungroup() %>%
      mutate(threat_group = x) 
    
    return(out) # Return the modified data
  })) 
  
  # Calculate differences in derivatives between additive and interactive threats
  post_system_diff <- do.call("rbind", lapply(additive_cols[grepl("\\+", additive_cols)], function(x){
    
    post_dydx_system %>%
      subset(threat %in% c(x, gsub(" \\+ ", ".", x))) %>% # Subset to shared additive and interactive threats (e.g., "threat1.threat2" and "threat1 + threat2")
      reframe(.value = diff(c(.value[2], .value[1])), .by = c(threat_group, .draw)) %>% # Find difference in derivatives between additive and interactive threats
      mutate(threats = gsub(" \\+ ", ".", x) ) # Name threat combination
    
  })) %>%
    mutate(System = names(system_model_ls)[X]) # Add system name to the data
  
  # Calculate interval estimates for the differences in derivatives
  post_interval_system_diff <- post_system_diff %>%
    na.omit() %>% # Drop missing threats 
    group_by(System, threat_group, threats) %>%
    ggdist::median_qi(.width = c(.95, .8, .5), .exclude = c(".draw")) %>%  # Extract distribution information
    mutate(interaction.type = case_when(
      .upper < 0 ~ "synergistic",    # If upper bound of the interval is less than 0, classify as synergistic
      .lower > 0 ~ "antagonistic",   # If lower bound of the interval is greater than 0, classify as antagonistic
      TRUE ~ "additive"              # Otherwise, classify as additive
    )) 
  
  # Return a list containing derivatives and interval estimates
  return(list("dydx_diff" = post_system_diff, "interval" = post_interval_system_diff))
})

# Substract the differences

post_mod_system_diff <- do.call("rbind",lapply(post_intvadd_system, `[[`, 'dydx_diff'))

# Substract the intervals

post_intvadd_interval_system_diff <- do.call("rbind",lapply(post_intvadd_system, `[[`, 'interval'))

# Calculate the position 

pos <- post_mod_system_diff %>% 
  group_by(threats, System) %>% 
  summarise(pos=median(.value))

# Join it

post_mod_system_diff <- post_mod_system_diff %>% 
  mutate(threats = gsub("invasive", "Invasive", threats),
         threats = gsub("habitatl", "Habitat loss", threats),
         threats = gsub("climatechange", "Climate change", threats),
         threats = gsub("pollution", "Pollution", threats),
         threats = gsub("exploitation", "Exploitation", threats),
         threats = gsub("disease", "Disease", threats),
         threat_group = gsub("invasive", "Invasive", threats),
         threat_group = gsub("habitatl", "Habitat loss", threats),
         threat_group = gsub("climatechange", "Climate change", threats),
         threat_group = gsub("pollution", "Pollution", threats),
         threat_group = gsub("exploitation", "Exploitation", threats),
         threat_group = gsub("disease", "Disease", threats)) 

post_intvadd_interval_system_diff <- post_intvadd_interval_system_diff %>% 
  mutate(threats = gsub("invasive", "Invasive", threats),
         threats = gsub("habitatl", "Habitat loss", threats),
         threats = gsub("climatechange", "Climate change", threats),
         threats = gsub("pollution", "Pollution", threats),
         threats = gsub("exploitation", "Exploitation", threats),
         threats = gsub("disease", "Disease", threats),
         threat_group = gsub("invasive", "Invasive", threats),
         threat_group = gsub("habitatl", "Habitat loss", threats),
         threat_group = gsub("climatechange", "Climate change", threats),
         threat_group = gsub("pollution", "Pollution", threats),
         threat_group = gsub("exploitation", "Exploitation", threats),
         threat_group = gsub("disease", "Disease", threats)) 

# Calculate proportions 

data_ad_syst <- post_intvadd_interval_system_diff %>%
  separate_rows(threats, sep = "\\.") %>% 
  group_by(System, threats, interaction.type) %>%
  summarise(n=n()) %>% 
  ungroup() %>% 
  complete(System, threats, interaction.type) %>%
  group_by(System, threats) %>% 
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))) %>% 
  ungroup() 

# Classify calculate population trends and whether they are additive, synergistic
# or antagonistic

# Define a function to operate on each element of the system_model_ls list
post_intvadd_taxon <- lapply(seq_along(taxon_model_ls), FUN = function(X){
  
  # Extract column names from the data within the current element of the list
  threatcols <- colnames(taxon_model_ls[[X]]$data)[grepl(paste(all_threats, collapse = "|"), colnames(taxon_model_ls[[X]]$data))] 
  
  # Create additive columns by combining the column names
  additive_cols <- do.call("c", lapply(strsplit(threatcols, "[.]"), function(x){
    paste(x, collapse = " + ")
  })) 
  
  # Combine unique threat columns and additive columns
  target_cols <- unique(c(threatcols, additive_cols))
  
  # Generate posterior draws for the taxon model with specified threats and columns
  postdraws_taxon <- threat_post_draws(model = taxon_model_ls[[X]],
                                        threat_comns = target_cols,
                                        ndraws = 1000,
                                        nuisance = c("series", "SpeciesName"),
                                        n.cores = 4) %>%
    mutate(combo_group = case_when(
      grepl("[.]", threat) ~ "interactive", # If threat contains a dot, classify as interactive
      grepl("\\+", threat) ~ "additive",    # If threat contains a plus sign, classify as additive
      TRUE ~ "single"))                     # Otherwise, classify as singular
  
  # Calculate derivatives for each threat
  post_dydx_taxon <- do.call("rbind", lapply(all_threats, function(x){
    
    out <- postdraws_taxon %>%
      subset(grepl(x, threat)) %>% # Filter to focal threat
      reframe(.value = mean(diff(.value) / diff(time)), .by = c(combo_group, threat, .draw)) %>% # Estimate each timeseries' first derivative
      group_by(threat) %>%
      ungroup() %>%
      mutate(threat_group = x) 
    
    return(out) # Return the modified data
  })) 
  
  # Calculate differences in derivatives between additive and interactive threats
  post_taxon_diff <- do.call("rbind", lapply(additive_cols[grepl("\\+", additive_cols)], function(x){
    
    post_dydx_taxon %>%
      subset(threat %in% c(x, gsub(" \\+ ", ".", x))) %>% # Subset to shared additive and interactive threats (e.g., "threat1.threat2" and "threat1 + threat2")
      reframe(.value = diff(c(.value[2], .value[1])), .by = c(threat_group, .draw)) %>% # Find difference in derivatives between additive and interactive threats
      mutate(threats = gsub(" \\+ ", ".", x) ) # Name threat combination
    
  })) %>%
    mutate(taxon = names(taxon_model_ls)[X]) # Add taxon name to the data
  
  # Calculate interval estimates for the differences in derivatives
  post_interval_taxon_diff <- post_taxon_diff %>%
    na.omit() %>% # Drop missing threats 
    group_by(taxon, threat_group, threats) %>%
    ggdist::median_qi(.width = c(.95, .8, .5), .exclude = c(".draw")) %>%  # Extract distribution information
    mutate(interaction.type = case_when(
      .upper < 0 ~ "synergistic",    # If upper bound of the interval is less than 0, classify as synergistic
      .lower > 0 ~ "antagonistic",   # If lower bound of the interval is greater than 0, classify as antagonistic
      TRUE ~ "additive"              # Otherwise, classify as additive
    )) 
  
  # Return a list containing derivatives and interval estimates
  return(list("dydx_diff" = post_taxon_diff, "interval" = post_interval_taxon_diff))
})

# Substract the differences

post_mod_taxon_diff <- do.call("rbind",lapply(post_intvadd_taxon, `[[`, 'dydx_diff'))

# Substract the intervals

post_intvadd_interval_taxon_diff <- do.call("rbind",lapply(post_intvadd_taxon, `[[`, 'interval'))

# Correct the data frame

post_mod_taxon_diff <- post_mod_taxon_diff %>% 
  mutate(threats = gsub("invasive", "Invasive", threats),
         threats = gsub("habitatl", "Habitat loss", threats),
         threats = gsub("climatechange", "Climate change", threats),
         threats = gsub("pollution", "Pollution", threats),
         threats = gsub("exploitation", "Exploitation", threats),
         threats = gsub("disease", "Disease", threats),
         threat_group = gsub("invasive", "Invasive", threats),
         threat_group = gsub("habitatl", "Habitat loss", threats),
         threat_group = gsub("climatechange", "Climate change", threats),
         threat_group = gsub("pollution", "Pollution", threats),
         threat_group = gsub("exploitation", "Exploitation", threats),
         threat_group = gsub("disease", "Disease", threats)) 

post_intvadd_interval_taxon_diff <- post_intvadd_interval_taxon_diff %>% 
  mutate(threats = gsub("invasive", "Invasive", threats),
         threats = gsub("habitatl", "Habitat loss", threats),
         threats = gsub("climatechange", "Climate change", threats),
         threats = gsub("pollution", "Pollution", threats),
         threats = gsub("exploitation", "Exploitation", threats),
         threats = gsub("disease", "Disease", threats),
         threat_group = gsub("invasive", "Invasive", threats),
         threat_group = gsub("habitatl", "Habitat loss", threats),
         threat_group = gsub("climatechange", "Climate change", threats),
         threat_group = gsub("pollution", "Pollution", threats),
         threat_group = gsub("exploitation", "Exploitation", threats),
         threat_group = gsub("disease", "Disease", threats)) 

# Calculate the proportions 

data_ad_tax <- post_intvadd_interval_taxon_diff %>% 
  separate_rows(threats, sep = "\\.") %>% 
  group_by(taxon, threats, interaction.type) %>%
  summarise(n=n()) %>% 
  ungroup() %>% 
  complete(taxon, threats, interaction.type) %>%
  group_by(taxon, threats) %>% 
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))) %>% 
  ungroup() 

# Save the tables 

write.csv2(data_ad_syst, "Results/Table S3.csv")
write.csv2(data_ad_tax, "Results/Table S4.csv")

## Panel a: System -------------------------------------------------------------

data <- data_ad_syst %>% 
  mutate(System=as.factor(System),
         freq=freq*100, 
         interaction.type=gsub("synergistic", "Synergy", interaction.type),
         interaction.type=gsub("additive", "Additive", interaction.type),
         interaction.type=gsub("antagonistic", "Antagonistic", interaction.type),
         interaction.type=fct_relevel(interaction.type,
                                      c("Synergy", 
                                        "Antagonistic", 
                                        "Additive"))) 

# Set a number of 'empty bar' to add at the end of each group

empty_bar <- 2
nObsType <- nlevels(as.factor(data$interaction.type))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$System)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$System <- rep(levels(data$System), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(System, threats)
data$id <- rep(seq(1, nrow(data)/nObsType), each=nObsType)

# Get the name and the y position of each label

label_data <- data %>% group_by(id, threats) %>% 
  summarise(tot=sum(freq))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(System) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start 
grid_data <- grid_data[-1,]

# Make the plot

(g3a <- ggplot(data) +
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=freq, 
                 fill=interaction.type), 
             stat="identity", alpha=0.8) +
    # Add a personalised grid
    geom_segment(data=grid_data, aes(x = end-0.5, y = 100, xend = start-0.5, yend = 100),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 80, xend = start-0.5, yend = 80),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 60, xend = start-0.5, yend = 60),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 40, xend = start-0.5, yend = 40),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 20, xend = start-0.5, yend = 20),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    # Colour of the fill
    scale_fill_manual(name = "",
                      values = c("#B32315",
                                 "#1E63B3",
                                 "#694364"))+
    # Add text showing the freq of each 100/75/50/25 lines
    annotate("text", x = rep(max(data$id),5), y = c(20,40, 60, 80,100), 
             label = c("20", "40", "60", "80","100"), 
             color="grey", size=5, angle=0, fontface="bold", hjust=1.5) +
    ylim(-100,220) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    # Add labels on top of each bar
    geom_text(data=label_data, aes(x=id, y=tot+5, label=threats, 
                                   hjust=hjust), 
              color="black", fontface="bold",alpha=0.6, size=5, 
              angle= label_data$angle, inherit.aes = FALSE ) +
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), 
                 colour = system_pal, alpha=0.8, size=2, inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -22, label=System), 
              hjust=c(1,0.5,0), colour = system_pal, alpha=0.8, size=4, 
              fontface="bold", inherit.aes = FALSE))

## Panel b: Taxon -------------------------------------------------------------

data <- data_ad_tax %>% 
  mutate(Taxon=as.factor(taxon),
         freq=freq*100, 
         interaction.type=gsub("synergistic", "Synergy", interaction.type),
         interaction.type=gsub("additive", "Additive", interaction.type),
         interaction.type=gsub("antagonistic", "Antagonistic", interaction.type),
         interaction.type=fct_relevel(interaction.type,
                                      c("Synergy", 
                                        "Antagonistic", 
                                        "Additive"))) 

# Set a number of 'empty bar' to add at the end of each group

empty_bar <- 2
nObsType <- nlevels(as.factor(data$interaction.type))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$Taxon)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$Taxon <- rep(levels(data$Taxon), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(Taxon, threats)
data$id <- rep(seq(1, nrow(data)/nObsType), each=nObsType)

# Get the name and the y position of each label

label_data <- data %>% 
  group_by(id, threats) %>% 
  summarise(tot=sum(freq))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(Taxon) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start 
grid_data <- grid_data[-1,]

# Make the plot

(g3b <- ggplot(data) +
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=freq, fill=interaction.type), 
             stat="identity", alpha=0.8) +
    # Add a personalised grid
    geom_segment(data=grid_data, aes(x = end-0.5, y = 100, xend = start-0.5, yend = 100),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 80, xend = start-0.5, yend = 80),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 60, xend = start-0.5, yend = 60),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 40, xend = start-0.5, yend = 40),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end-0.5, y = 20, xend = start-0.5, yend = 20),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    # Colour of the fill
    scale_fill_manual(name = "",
                      values = c("#B32315",
                                 "#1E63B3",
                                 "#694364"))+
    # Add text showing the freq of each 100/75/50/25 lines
    annotate("text", x = rep(max(data$id),5), y = c(20,40, 60, 80,100), 
             label = c("20", "40", "60", "80","100"), 
             color="grey", size=5, angle=0, fontface="bold", hjust=1) +
    ylim(-100,220) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    # Add labels on top of each bar
    geom_text(data=label_data, aes(x=id, y=tot+5, label=threats, hjust=hjust), 
              color="black", fontface="bold", alpha=0.6, size=5, 
              angle= label_data$angle, inherit.aes = FALSE ) +
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), 
                 colour = taxon_pal, alpha=0.8, size=2 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -22, label=Taxon), 
              hjust=c(0.9,0.8,0.5,0,0), 
              colour = taxon_pal, alpha=0.8, size=4, 
              fontface="bold",inherit.aes = FALSE))

# Combine ----------------------------------------------------------------------

#Create first row 

(row1 <- plot_grid(g3a, g3b, labels="auto"))

# Get the legend

legend <- get_legend(g3a+theme(legend.position = "bottom",
                               legend.text = element_text(size=14)))
# Plot it 

(figure3<- plot_grid(row1, legend, 
                     nrow = 2, rel_heights = c(1,0.05)))


# Save it

ggsave("Results/Figure 3.pdf", figure3, 
       width = 10, height = 8)

# 
# ## Panel a: System -------------------------------------------------------------
# 
# (g3a <- post_mod_system_diff %>%
#   drop_na(.value) %>% 
#   ggplot(aes(x = .value, y=reorder(threats, .value))) +
#   ggdist::stat_slab(data = na.omit(post_mod_system_diff) %>%
#                       merge(select(post_intvadd_interval_system_diff,-.value),
#                             by = c("System","threat_group","threats")) %>%
#                       subset(.width == 0.8)
#                     ,aes(fill = interaction.type),alpha=0.5,normalize = "xy") +
#   ggdist::geom_pointinterval(data = post_intvadd_interval_system_diff,
#                              aes(xmin = .lower, xmax = .upper,col = interaction.type),alpha=0.5) +
#   geom_point(data = subset(post_intvadd_interval_system_diff,.width == 0.8), 
#              aes(x = .value,col = interaction.type,fill = interaction.type),shape = 21,alpha=0.5,size = 3) +
#   geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
#   coord_cartesian(xlim = c(-0.25,0.25)) + 
#   scale_fill_manual(values = c("#694364",
#                                "#1E63B3",
#                                "#B32315"), name = "",
#                     guide = guide_legend(override.aes = list(color = NA,shape = 2) )) + 
#   scale_color_manual(values = c("#694364",
#                                 "#1E63B3",
#                                 "#B32315"), guide = "none") + 
#   facet_wrap(~System,ncol = 3) + 
#   xlab( expression(paste("Additive ",partialdiff,"y","/",partialdiff,"x"," - interactive ",partialdiff,"y","/",partialdiff,"x"))) + 
#   ylab("Threat combination"))
# 
# 
# ## Panel b: Taxon --------------------------------------------------------------
# 
# (g3b <- post_mod_taxon_diff %>%
#    drop_na(.value) %>% 
#    ggplot(aes(x = .value, y=reorder(threats, .value))) +
#    ggdist::stat_slab(data = na.omit(post_mod_taxon_diff) %>%
#                        merge(select(post_intvadd_interval_taxon_diff,-.value),
#                              by = c("taxon","threat_group","threats")) %>%
#                        subset(.width == 0.8)
#                      ,aes(fill = interaction.type),alpha=0.5,normalize = "xy") +
#    ggdist::geom_pointinterval(data = post_intvadd_interval_taxon_diff,
#                               aes(xmin = .lower, xmax = .upper,col = interaction.type),alpha=0.5) +
#    geom_point(data = subset(post_intvadd_interval_taxon_diff,.width == 0.8), 
#               aes(x = .value,col = interaction.type,fill = interaction.type),shape = 21,alpha=0.5,size = 3) +
#    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
#    coord_cartesian(xlim = c(-0.25,0.25)) + 
#    scale_fill_manual(values = c("#694364",
#                                 "#1E63B3",
#                                 "#B32315"), name = "",
#                      guide = guide_legend(override.aes = list(color = NA,shape = 2) )) + 
#    scale_color_manual(values = c("#694364",
#                                  "#1E63B3",
#                                  "#B32315"), guide = "none") + 
#    facet_wrap(~taxon,ncol = 5) + 
#    xlab( expression(paste("Additive ",partialdiff,"y","/",partialdiff,"x"," - interactive ",partialdiff,"y","/",partialdiff,"x"))) + 
#    ylab("Threat combination"))
# 
# # Combine ----------------------------------------------------------------------
# 
# # Create titles
# 
# title1 <- ggdraw() + 
#   draw_label("System",
#              size = 16,
#              x = 0.5,
#              vjust = 0) 
# 
# title2 <- ggdraw() + 
#   draw_label("Taxon",
#              x = 0.5,
#              size = 16,
#              vjust = 0) 
# subtitles <- plot_grid(title1, title2, nrow = 1)
# 
# # Create row 1
# 
# (row1 <- plot_grid(g3a+theme(legend.position = "none",
#                              axis.title.x = element_blank()),
#                    g3b+theme(legend.position = "none",
#                              axis.title.x = element_blank()),
#                    labels = "auto", 
#                    ncol=1))
# 
# # Greate two legends 
# 
# legend <- plot_grid(get_legend(g3a+theme(legend.position = "bottom",
#                                           legend.text = element_text(size = 14))))
# 
# # Put them together 
# 
# (figure3<- plot_grid(row1, legend,
#                      rel_heights = c(1,.1), ncol=1))
# 
# # Save it
# 
# ggsave("Results/Figure 3.pdf",figure3,
#        width = 14, height = 14)
