# --------------------------------------------------------------------------------------- #
# - FILE NAME:   figures_final.R         
# - DATE:        11/02/2024
# - DESCRIPTION: Code to plot the figures 1, 2 and 3 from the main manuscrpit.
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
library(waffle)

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
                  plot.title = element_text(hjust = 0.5),
                  plot.margin = unit(c(0, 0, 0, 0), "cm")))

#Create a palette for the threats 

threat_palette<- c(met.brewer(name="Hokusai1", n=6,type="continuous"))

# Color palette for taxons

taxon_pal <- wes_palette("Cavalcanti1", n = 5)

# Load other customised functions to substract the trends from the plots

source("Code/utils/prep_data_grid_fn.R")
source("Code/utils/threat_post_draws.R")

# Load the models 

mod <- readRDS("Results/models/mod_global_rerun.RDS")

# Load the model data

load("Data/data_models.RData")

# Filter the data to unique time series

data <- mod_dat_full %>% distinct(ID, .keep_all = T)

# Figure 1: Map and trends #####################################################
## Map -------------------------------------------------------------------------

# Set the world map

world <- map_data("world")

# Create map

(g1a <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", linewidth = 0.3) +
    #coord_proj("+proj=wintri") +
    theme_map() + 
    geom_point(data = data, 
               aes(x = Longitude, y = Latitude,
                   fill=Taxon), 
               alpha=0.5, shape=21, size=3) +
    scale_y_continuous(limits = c(-80, 80)) +
    scale_fill_manual("",values = taxon_pal)+
    guides(fill=guide_legend(nrow=2,byrow=TRUE,
                             override.aes = list(alpha = 1)))+
    expand_limits(x = 0, y = 0)+
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(size = 12, hjust =0.5),
          legend.text = element_text(size = 10),
          plot.margin = unit(c(0,-.5,0,0), units = , "cm")))

# Legend

legend1 <- cowplot::get_plot_component(g1a, 'guide-box-top', return_all = TRUE)

# Create distribution 1

anot <- data.frame(x=median(data$Duration), 
                   y=300, 
                   label =paste("Median = ", median(data$Duration)))


(gd1 <- data %>% 
    ggplot() +
    geom_histogram(aes(x = Duration),
                   binwidth = 2,
                   alpha = .6,
                   position = "stack",
                   fill = "grey40", color = "white") +
    geom_vline(aes(xintercept = median(Duration)),
               linetype = "longdash", size=0.5) +
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(breaks = seq(10, 60, 10),
                       expand = c(0,0))+
    ggrepel::geom_text_repel(data = anot, aes(x=x, y=y+100,
                                              label =label),
                             nudge_x = 10, nudge_y = 0)+
    labs(x="Duration (years)", 
         y="Number of datasets")+
    theme(axis.title.x = element_text(margin = margin(t = 0)),
          axis.text.y = element_text(margin = margin(r = 0))))


# Combined map 

(g1a<- ggdraw(g1a+theme(legend.position = c(0.40,0.08))) +
   draw_plot(gd1, x = -0.32, y = -0.3143, scale = 0.35))

## Panel b: threat proportions ------------------------

# Readjust the data frame to contain the threats as individual columns

threat_freq <- mod$data %>% 
  distinct(series, .keep_all=T) %>%
  dplyr::select(!dplyr::contains(".")) %>% 
  dplyr::select(-c(y_centered, scaled_year, SpeciesName, Site, time)) %>% 
  pivot_longer(1:6, names_to="threats") %>% 
  filter(value!=0) %>% 
  group_by(threats) %>% 
  summarise(n=sum(as.numeric(value))) %>% 
  rbind(tibble(threats="None", 
         n=mod$data %>%
           distinct(series, .keep_all = TRUE) %>%
           summarise(count = n()) %>%
           pull(count)-sum(.$n))) %>% 
  mutate(total=mod$data %>%
           distinct(series, .keep_all = TRUE) %>%
           summarise(count = n()) %>%
           pull(count),
         freq=(n/total),
         threats=ifelse(threats=="climatechange", "Climate change",
                        ifelse(threats=="exploitation", 
                               "Exploitation",
                               ifelse(threats=="invasive", 
                                      "Invasive",ifelse(threats=="disease", 
                                                        "Disease",ifelse(threats=="habitatl",
                                                                         "Habitat loss", 
                                                                         ifelse(threats=="pollution",
                                                                                "Pollution", threats))))))) 
# We create the palette

palette <- data.frame(threats = c("None", "Pollution","Habitat loss",
                      "Climate change", "Invasive", 
                      "Exploitation", "Disease"),
                      fill_col = as.factor(c("grey50", threat_palette)))


# Plot

(g1b <- threat_freq %>% 
    left_join(palette, by="threats") %>% 
    mutate(threats=factor(threats,
                          levels=c("None", "Disease","Invasive",
                                   "Climate change", "Pollution",
                                   "Habitat loss", "Exploitation"))) %>% 
    ggplot(aes(fill = fill_col,
               y=freq, x=threats)) + 
    geom_bar(stat = "identity") +
    scale_fill_manual(values = levels(palette$fill_col), 
                      guide=NULL)+
    scale_y_continuous(breaks = seq(0, 1, .1), 
                       label = scales::percent,
                       expand = c(0,0))+ 
    labs(y="Proportion of threats (%)", x="",fill="") +
    # geom_text(aes(label = paste0(round(freq*100, 0), "%")),
    #           size=5,
    #           position = position_stack(vjust = 0.5)) +
    theme(legend.position="bottom",
          legend.text = element_text(size=12),
          strip.text = element_text(hjust = 0),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size = 14), 
          plot.margin = unit(c(0.5,0,0,0), units = "cm")))

## Final figure 1 --------------------------------------------------------------

(fig1 <- plot_grid(g1a, g1b,
                  nrow = 2,rel_heights = c(1,.7),
                  labels = "auto"))

# Save it

ggsave("Results/Figure 1.pdf", fig1, 
       height = 10,width = 10)

# Remove unused things 

rm(fig1, g1a, g1b, gd1, legend1, anot,world)

# Figure 2: Effect of threats to population trends and random effects ##########
## Panel a: Coefficients -------------------------------------------------------
## Prepare data 

# Obtain the coefficients and clean the data

coefs_df <- mod %>% 
  gather_draws(`b_scaled_year.*`,
               regex = TRUE) %>%
  mutate(.variable = gsub("b_", "", .variable),
         .variable = gsub("scaled_year:", "", .variable),
         .variable = gsub("invasive", "Invasive", .variable),
         .variable = gsub("habitatl", "Habitat loss", .variable),
         .variable = gsub("climatechange", "Climate change", .variable),
         .variable = gsub("pollution", "Pollution", .variable),
         .variable = gsub("exploitation", "Exploitation", .variable),
         .variable = gsub("disease", "Disease", .variable),
         .variable = gsub("scaled_year", "None", .variable),
         .variable = gsub(" 1", "", .variable))


# List of threats 

all_threats <- c("None","Pollution","Habitat loss","Climate change",
                 "Invasive", "Exploitation","Disease")

# Obtain the coeficients for each threat group 

coefs_group <- do.call("rbind", lapply(all_threats, function(x){
  out <- coefs_df %>%
    subset(grepl(x, .variable)) %>% # Filter to the current threat
    mutate(int_group = ifelse(grepl("\\.", .variable), "Interactive", "Singular"), # Classify threats as 'singular' or 'interactive'
           threat_group = x) %>% # Assign the current threat to a grouping variable
    ungroup()
  return(out)
})) 

# Calculate median and quantile intervals for the rate of change by threat group and interaction type

coefs_interval <- coefs_group  %>%
  group_by(threat_group, int_group) %>%
  ggdist::median_qi(.width = c(.95, .8, .5), 
                    .exclude = c(".chain", ".iteration", ".draw", ".variable")) %>%
  mutate(int_group = ifelse(int_group == "Singular", "Singular", "Interactive"),
         int = ordered(int_group, c("Singular", "Interactive")))

# Create a data frame matching threat groups to colors

palette <- data.frame(threat_group = unique(coefs_group$threat_group),
                      fill_col = c("grey50", threat_palette))

# Prepare the data for plotting by joining with the palette and adjusting factors

plot_coefs <- coefs_group %>%
  left_join(palette, by = "threat_group") %>%
  mutate(fill_col = factor(fill_col),
         int = ordered(int_group, c("Singular", "Interactive")))

# Single threats

# (g2a <- plot_coefs %>% 
#     ggplot(aes(x = .value,
#                y=reorder(threat_group, .value))) +
#     tidybayes::stat_slab(data = subset(plot_coefs,
#                                        int == "Singular"),
#                          aes(fill = fill_col,
#                              group = .variable), alpha=0.5,
#                          normalize = "groups") +
#     ggdist::geom_pointinterval(data = coefs_interval %>% 
#                                  filter(int=="Singular"),
#                                aes(xmin = .lower, xmax = .upper),
#                                position = position_dodge()) +
#     geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
#     labs(x="Effect of threats on population trends", y=NULL) + 
#     coord_cartesian(xlim = c(-0.2,0.2)) +
#     scale_x_continuous(breaks= seq(-0.2,0.2,by=0.2)) + 
#     ylab("")+
#     scale_fill_manual(values = levels(plot_coefs$fill_col),
#                       guide = "none") + 
#     scale_color_manual(values = levels(plot_coefs$fill_col),
#                        guide = "none"))
# 
# # Create a fake rectangle
# 
# rect <- data.frame(start=0.21, end=0.9)
# 
# # Interactive threats
# 
# (g2 <- plot_coefs %>% 
#     filter(int_group=="Singular" & threat_group=="None"|
#              int_group == "Interactive")%>%
#     mutate(threat_group=ordered(threat_group, levels=c("Invasive",
#                                                        "Disease",
#                                                        "Climate change",
#                                                        "Pollution",
#                                                        "Exploitation",
#                                                        "Habitat loss", 
#                                                        "None")),
#            .value=ifelse(threat_group=="None", NA, .value)) %>% 
#     ggplot(aes(x = .value,
#                y= threat_group)) +
#     tidybayes::stat_slab(aes(fill = fill_col,
#                              group = .variable), alpha=0.5,
#                          normalize = "panels") +
#     ggdist::geom_pointinterval(data = coefs_interval %>% 
#                                  filter(int_group == "Interactive"| 
#                                           int_group=="Single" & threat_group=="None"),
#                                aes(xmin = .lower, xmax = .upper),
#                                position = position_dodge()) +
#     annotate("rect", 
#               xmin=0.21, xmax=2.1,
#              ymin=-Inf, ymax=+Inf, 
#               fill="white")+
#     geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
#     labs(x="Effect of threats on population trends", y=NULL) + 
#     coord_cartesian(xlim = c(-0.2, 1.9)) +
#     scale_x_continuous(breaks= seq(-0.2,0.2,by=0.2)) + 
#     ylab("")+
#     scale_fill_manual(values = levels(plot_coefs$fill_col),
#                       guide = "none") + 
#     scale_color_manual(values = levels(plot_coefs$fill_col),
#                        guide = "none")+
#     theme(axis.text.y=element_blank(),
#           axis.ticks.y = element_blank()))

(g2a <- plot_coefs %>%
    ggplot(aes(x = .value,
               y=reorder(threat_group, .value))) +
    tidybayes::stat_slab(data = plot_coefs %>%
                           filter(int_group == "Singular"),
                         aes(fill = fill_col,
                             group = .variable), alpha=0.5,
                         normalize = "groups") +
    tidybayes::stat_slab(data = subset(plot_coefs,
                                       int_group == "Interactive"),
                         aes(fill = fill_col,
                             group = .variable), alpha=0.5,
                         normalize = "panels") +
    ggdist::geom_pointinterval(data = coefs_interval,
                               aes(xmin = .lower, xmax = .upper),
                               position = position_dodge()) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    labs(x="Effect on population trends", y=NULL) +
    coord_cartesian(xlim = c(-0.2,0.2)) +
    scale_x_continuous(breaks= seq(-0.1,0.1,by=0.1)) +
    facet_wrap(~int) +
    ylab("")+
    scale_fill_manual(values = levels(plot_coefs$fill_col),
                      guide = "none") +
    scale_color_manual(values = levels(plot_coefs$fill_col),
                       guide = "none")+
    theme(axis.title.y=element_blank()))

## Panel b: Population trends explained by interactions vs random effects ------

# Identify the colnames for the threats

threat_cols <-  colnames(mod$data)[grepl(paste(c("pollution","habitatl",
                                                "climatechange","invasive",
                                                "exploitation","disease"),
                                              collapse = "|"),
                                        colnames(mod$data))]

# Get the interaction data

int_dat <- mod$data %>%
  mutate(int_present = if_any(threat_cols,~.x == "1",TRUE),
         y_centered = NA) 

# Generate the null data

null_dat <- int_dat %>%
  mutate(across(all_of(threat_cols),
                ~factor("0",levels = c("0","1"))))

# Base predictions

base_pred <- brms::posterior_epred(mod,
                                   newdata = null_dat,
                                   re.form = NA,
                                   incl_autocor = F,
                                   sort = TRUE,
                                   ndraws = 1000) %>% #extract posterior draws for the above data grid
  as.data.frame() %>%
  mutate(.draw = 1:NROW(.)) %>%
  #extract posterior draws for the data used to create the model
  pivot_longer(-.draw,names_to = "index",values_to = ".value") %>%
  cbind(null_dat %>% dplyr::select(series, time,int_present), row.names = NULL) %>% 
  # reframe(.value = mean(.value/lag(.value)),
  #         .by = c(series,.draw,int_present)) %>% #estimate derivatives per series and .draw
  reframe(.value = mean(diff(.value)/diff(time)),
          .by = c(series,.draw,int_present)) %>% #estimate derivatives per series and .draw
  reframe(.value = mean(.value),.by = c(series,int_present)) #estimate average derivative per series

# Random effects

space_pred <- brms::posterior_epred(mod,
                                   newdata = null_dat,
                                   re.form = NULL,
                                   incl_autocor = F,
                                   sort = TRUE,
                                   ndraws = 1000) %>% #extract posterior draws for the above data grid
  as.data.frame() %>%
  mutate(.draw = 1:NROW(.)) %>%
  #extract posterior draws for the data used to create the model
  pivot_longer(-.draw,names_to = "index",values_to = ".value") %>%
  cbind(null_dat %>% dplyr::select(series, time,int_present), row.names = NULL) %>% 
  # reframe(.value = mean(.value/lag(.value)),
  #         .by = c(series,.draw,int_present)) %>% #estimate derivatives per series and .draw
  reframe(.value = mean(diff(.value)/diff(time)),
          .by = c(series,.draw,int_present)) %>% #estimate derivatives per series and .draw
  reframe(.value = mean(.value),.by = c(series,int_present)) #estimate average derivative per series

# Predictions with interactive effects on

inter_pred <- brms::posterior_epred(mod,
                                   newdata = int_dat,
                                   re.form = NA,
                                   incl_autocor = F,
                                   sort = TRUE,
                                   ndraws = 1000) %>% #extract posterior draws for the above data grid
  as.data.frame() %>%
  mutate(.draw = 1:NROW(.)) %>%
  #extract posterior draws for the data used to create the model
  pivot_longer(-.draw,names_to = "index",values_to = ".value") %>%
  cbind(null_dat %>% dplyr::select(series, time,int_present), row.names = NULL) %>% 
  # reframe(.value = mean(.value/lag(.value)),
  #        .by = c(series,.draw,int_present)) %>% #estimate derivatives per series and .draw
  reframe(.value = mean(diff(.value)/diff(time)),
          .by = c(series,.draw,int_present)) %>% #estimate derivatives per series and .draw
  reframe(.value = mean(.value),.by = c(series,int_present)) #estimate average derivative per series

# Calculate the differences

inter_diff <- rbind(base_pred,inter_pred) %>%
  reframe(beta = diff(.value),.by = c(series,int_present)) %>%
  mutate(type = "Interactive effects")

space_diff <- rbind(base_pred,space_pred) %>%
  reframe(beta = diff(.value),.by = c(series,int_present))%>%
  mutate(type = "Random effects")

# Join them

diff <- rbind(inter_diff,space_diff) %>%
  filter(int_present == TRUE)

# Calculate median and quantile intervals for the rate of change by threat group and interaction type

diff_interval <- diff  %>%
  group_by(type) %>%
  ggdist::median_qi(.width = c(.95, .8, .5), 
                    .exclude = c("series", "int_present", "type")) 

# Plot it

(g2b <- diff %>% 
    ggplot(aes(x = beta, y = type, fill=type)) +
    ggdist::stat_slab(.width = c(0.6,0.8,0.95), alpha=.9) + 
    ggdist::geom_pointinterval(data = diff_interval,
                               aes(xmin = .lower, xmax = .upper),
                               position = position_dodge()) +
    geom_vline(aes(xintercept = 0),  
               linetype = "dashed", 
               colour="grey50") +
    scale_fill_manual("", values=c("#739EF0","#F0C873" ))+
    labs(x = "Trend difference", y="") +
    #xlim(-.3,.3)+
    coord_cartesian(xlim = c(-0.2,0.2)) +
    scale_x_continuous(breaks= seq(-0.2,0.2,by=0.1)) +
    theme(legend.position = "none"))

## Panel c: Population trends projections --------------------------------------

# The threats to project

threats <- c("pollution", "habitatl", 
             "climatechange", "invasive", 
             "exploitation", "disease")

threats_col <-  colnames(mod$data)[grepl(paste(threats,
                                               collapse = "|"),
                                         colnames(mod$data))] 
# Project the population

postdraws <- threat_post_draws(model = mod,
                               threat_comns = c("none",threats_col),
                               ndraws = 1000,
                               nuisance = c("series","SpeciesName"),
                               n.cores = 4)

# Calculate the derivative

post_proj <- do.call("rbind", lapply(c("none",threats), function(x){
  out <- postdraws %>%
    subset(grepl(x, threat)) %>% # Filter to the current threat
    reframe(.value = mean(diff(.value)/diff(time)),
            .by = c(threat, .draw)) %>% # Calculate the first derivative (rate of change) for each time series
    ungroup() %>%
    mutate(int_group = ifelse(grepl("\\.", threat), "combined", "single")) %>% # Classify threats as 'single' or 'combined'
    mutate(threat_group = x) %>% # Assign the current threat to a grouping variable
    # Rename threat groups for clarity
    mutate(threat_group = gsub("invasive", "Invasive", threat_group),
           threat_group = gsub("habitatl", "Habitat\nloss", threat_group),
           threat_group = gsub("climatechange", "Climate\nchange", threat_group),
           threat_group = gsub("pollution", "Pollution", threat_group),
           threat_group = gsub("exploitation", "Exploitation", threat_group),
           threat_group = gsub("disease", "Disease", threat_group),
           threat_group = gsub("none", "None", threat_group)) %>%
    group_by(threat) %>%
    ungroup()
  return(out)
}))

# Calculate median and quantile intervals for the rate of change by threat group and interaction type

dydx_interval <- post_proj  %>%
  group_by(threat_group, int_group) %>%
  ggdist::median_qi(.width = c(.95, .8, .5), 
                    .exclude = c(".draw", "threat", 
                                 "threat_group")) %>% 
  mutate(int_group = ifelse(int_group == "single", "Singular", "Interactive"),
         int = ordered(int_group, c("Singular", "Interactive")))

# Create a data frame matching threat groups to colors

palette <- data.frame(threat_group = unique(post_proj$threat_group),
                      fill_col = c("grey50", threat_palette))

# Prepare the data for plotting by joining with the palette and adjusting factors

plot_dydx_threats <- post_proj %>%
  left_join(palette, by = "threat_group") %>%
  mutate(fill_col = factor(fill_col), 
         int_group = ifelse(int_group == "single", "Singular", "Interactive"),
         int = ordered(int_group, c("Singular", "Interactive")))

# Figure with no threat

# (g2b <- ggplot(data = plot_dydx_threats,
#               aes(x = .value,
#                   y=reorder(threat_group, .value))) +
#     tidybayes::stat_slab(data = plot_dydx_threats %>%
#                            filter(int_group == "Singular"),
#                          aes(fill = fill_col,
#                              group = threat_group), alpha=0.5,
#                          normalize = "groups") +
#     tidybayes::stat_slab(data = subset(plot_dydx_threats,
#                                        int_group == "Interactive"),
#                          aes(fill = fill_col,
#                              group = threat), alpha=0.5,
#                          normalize = "panels") +
#     ggdist::geom_pointinterval(data = dydx_interval,
#                                aes(xmin = .lower, xmax = .upper),
#                                position = position_dodge()) +
#     geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
#     xlab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) +
#     ylab("")+
#     facet_wrap(~int) +
#     scale_fill_manual(values = levels(plot_dydx_threats$fill_col),
#                        guide = "none")+
#     coord_cartesian(xlim = c(-0.2,0.2))+
#     theme(plot.margin = unit(c(0,0.5,0,.5),"cm"),
#           axis.title.y=element_blank()))

(g2c <- ggplot(data = plot_dydx_threats,
               aes(y = .value,
                   x=reorder(threat_group, desc(.value)))) +
    ggdist::stat_halfeye(data = plot_dydx_threats %>%
                           filter(int_group == "Singular"),
                         aes(fill = fill_col,
                             group = threat_group), alpha=0.5,
                         adjust = .5, width = .3, .width = 0, 
                         justification = -.5, 
                         point_colour = NA) +
    geom_boxplot(data = plot_dydx_threats %>%
                   filter(int_group == "Singular"),
                 aes(fill = fill_col,
                     group = threat_group),  
                 width = .2, outliers = F, alpha=.7)+
    ggdist::stat_halfeye(data = plot_dydx_threats %>%
                           filter(int_group == "Interactive"),
                         aes(fill = fill_col,
                             group = threat), alpha=0.5,
                         adjust = .5, width = .3, .width = 0, 
                         justification = -.5, 
                         point_colour = NA) +
    geom_boxplot(data = plot_dydx_threats %>%
                   filter(int_group == "Interactive"),
                 aes(fill = fill_col,
                     group = threat_group),  
                 width = .2, outliers = F, alpha=.7)+
    geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
    facet_wrap(~int) +
    ylab(expression(paste("Population trend (",Delta,"y","/",Delta,"x)"))) +
    xlab("")+
    facet_wrap(~int,scales = "free") +
    scale_fill_manual(values = levels(plot_dydx_threats$fill_col),
                      guide = "none")+
    coord_cartesian(ylim = c(-0.2,0.2))+
    theme(plot.margin = unit(c(0,0.5,0,.5),"cm"),
          axis.title.x=element_blank()))


## Final figure 2 ------------------------------------------------------------

layout <- "
AAABB
CCCCC
"

(fig2 <-g2a+ g2b + free(g2c)+
    plot_layout(design = layout,widths = 1)+
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold")))

# Save it

ggsave("Results/Figure2.pdf", fig2,
       height = 8, width = 14)

# Figure 3: Interactive threats ################################################
# Prepare data -----------------------------------------------------------------
# Define a vector of threat categories

all_threats <- c("pollution", "habitatl", "climatechange", "invasive",
                 "exploitation", "disease")

# Extract column names from the model

threatcols <- colnames(mod$data)[grepl(paste(all_threats, collapse = "|"), 
                                       colnames(mod$data))] 

# Generate additive combinations of threat columns (i.e., "threat1 + threat2")

additive_cols <- do.call("c", lapply(strsplit(threatcols, "[.]"), 
                                     function(x) {
                                       paste(x, collapse = " + ")
                                     }))

# Combine original and additive columns into a unique set

target_cols <- unique(c(threatcols, additive_cols))

# Estimate posterior time series for each threat combination (singular, interactive, additive)

postdraws_intvadd <- threat_post_draws(model = mod,
                                       threat_comns = target_cols,
                                       nuisance = c("series", "SpeciesName"),
                                       n.cores = 4,
                                       ndraws = 1000) %>%
  mutate(combo_group = case_when(
    grepl("[.]", threat) ~ "interactive",
    grepl("\\+", threat) ~ "additive",
    TRUE ~ "single"))

# Calculate the first derivative of the value over time for each threat and categorise them

post_dydx_intvadd <- do.call("rbind",lapply(all_threats,function(x){
  
  out <- postdraws_intvadd %>%
    subset(grepl(x,threat)) %>% #filter to focal threat
    reframe(.value = mean(diff(.value)/diff(time)),
            .by = c(combo_group,threat,.draw)) %>% #estimate each timeseries' first derivative
    group_by(threat) %>%
    #filter(!any(.value >= abs(0.5))) %>% #drop highly variable threats
    ungroup() %>%
    mutate(threat_group = x) 
  
  return(out)
})) 

# Extract distribution information for each threat, categorised by interaction type

dydx_interval_intvadd <- post_dydx_intvadd %>%
  group_by(threat_group, combo_group, threat) %>%
  ggdist::median_qi(.width = c(.95, .8, .5), .exclude = c(".draw"))

# Compare the difference in derivatives between shared additive and interactive threats

post_intvadd_diff <- do.call("rbind", lapply(additive_cols[grepl("\\+", additive_cols)], function(x) {
  ss <- post_dydx_intvadd %>%
    subset(threat %in% c(x, gsub(" \\+ ", ".", x))) %>%
    reframe(.value = diff(c(.value[!grepl("[.]",threat)], .value[grepl("[.]",threat)])), .by = c(threat_group, .draw)) %>%
    mutate(threats = gsub(" \\+ ", ".", x))
}))

# Finalize the comparison, drop missing values, and classify interaction types

post_interval_intvadd_diff <- post_intvadd_diff %>%
  na.omit() %>%
  group_by(threat_group, threats) %>%
  ggdist::median_qi(.width = c(.95, .8, .5), .exclude = c(".draw")) %>%
  mutate(interaction.type = case_when(
    .upper < 0 ~ "synergistic",
    .lower > 0 ~ "antagonistic",
    TRUE ~ "additive"
  ))

# Calculate the proportions 

data_ad <- post_interval_intvadd_diff %>% 
  separate_rows(threats, sep = "\\.") %>% 
  mutate(interaction.type=factor(interaction.type, levels=c("synergy", "additive", "antagonistic"))) %>% 
  group_by(threats, interaction.type) %>%
  summarise(n=n()) %>% 
  ungroup() %>% 
  complete(threats, interaction.type) %>%
  group_by(threats) %>% 
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))) %>% 
  ungroup() 

## Panel a: Waffle plot proportions --------------------------------------------

# We will have to create a separate legend because not all levels are shown given
# that no synergies were found

# Create the base data 

data <- data_ad %>% 
  mutate(freq=freq*100, 
         interaction.type=gsub("synergy", "Synergy", interaction.type),
         interaction.type=gsub("additive", "Additive", interaction.type),
         interaction.type=gsub("antagonistic", "Antagonistic", interaction.type),
         interaction.type=fct_relevel(interaction.type,
                                      c("Synergy", 
                                        "Antagonistic", 
                                        "Additive")),
         threats = gsub("invasive", "Invasive", threats),
         threats = gsub("habitatl", "Habitat loss", threats),
         threats = gsub("climatechange", "Climate change", threats),
         threats = gsub("pollution", "Pollution", threats),
         threats = gsub("exploitation", "Exploitation", threats),
         threats = gsub("disease", "Disease", threats),
         threats=as.factor(threats)) 

# Create a fake data with synergies and plot it

legend <- cowplot::get_plot_component(data %>% 
  mutate(freq=ifelse(freq==0, 1, freq)) %>% 
  ggplot(aes(fill=interaction.type, values=freq))+
  geom_waffle(color = "white", size = .25, 
              n_rows = 10, 
              flip = TRUE, na.rm = F,
              make_proportional = T, 
              show.legend = TRUE) +
  facet_wrap(~threats, nrow = 2, strip.position = "bottom") +
  scale_x_discrete() + 
  scale_y_continuous(labels = function(x) x * 10, # make this multiplyer the same as n_rows
                     expand = c(0,0))+
  scale_fill_manual(name = NULL,
                    values = c("#B32315",
                               "#1E63B3",
                               "#694364"))+
  coord_equal(clip = "off")+
  theme_void()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        strip.text = element_text(size=14),
        plot.margin = margin(-10, 0, -10, 0)),
  'guide-box-top', return_all = TRUE)

# Make the plot 

(g3a <- data %>% 
   ggplot(aes(fill=interaction.type, values=freq))+
  geom_waffle(color = "white", size = .25, 
              n_rows = 10, 
              flip = TRUE, na.rm = F,
              make_proportional = T, 
              show.legend = TRUE) +
  facet_wrap(~threats, nrow = 2, strip.position = "bottom") +
  scale_x_discrete() + 
  scale_y_continuous(labels = function(x) x * 10, # make this multiplyer the same as n_rows
                     expand = c(0,0))+
  scale_fill_manual(name = NULL,
                    values = c("#1E63B3",
                               "#694364"))+
  coord_equal(clip = "off")+
  theme_void()+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        strip.text = element_text(size=14),
        plot.margin = margin(0, 0, 0, 0)))

## Panel b: Summary of the trend difference -------------------------------------

# Correct the data 

diff <- post_intvadd_diff %>%
  mutate(threats = gsub("invasive", "Invasive", threats),
         threats = gsub("habitatl", "Habitat loss", threats),
         threats = gsub("climatechange", "Climate change", threats),
         threats = gsub("pollution", "Pollution", threats),
         threats = gsub("exploitation", "Exploitation", threats),
         threats = gsub("disease", "Disease", threats)) %>% 
  drop_na() %>% 
  arrange(.value)

diff_int <- post_interval_intvadd_diff %>% 
  mutate(threats = gsub("invasive", "Invasive", threats),
         threats = gsub("habitatl", "Habitat loss", threats),
         threats = gsub("climatechange", "Climate change", threats),
         threats = gsub("pollution", "Pollution", threats),
         threats = gsub("exploitation", "Exploitation", threats),
         threats = gsub("disease", "Disease", threats))

slab_data <- na.omit(post_intvadd_diff) %>%
  merge(select(post_interval_intvadd_diff,-.value),
        by = c("threat_group","threats")) %>%
  subset(.width == 0.8) %>% 
  mutate(threats = gsub("invasive", "Invasive", threats),
         threats = gsub("habitatl", "Habitat loss", threats),
         threats = gsub("climatechange", "Climate change", threats),
         threats = gsub("pollution", "Pollution", threats),
         threats = gsub("exploitation", "Exploitation", threats),
         threats = gsub("disease", "Disease", threats))

# Plot it 

(g3b<- ggplot(data = diff, 
       aes(x = .value, y=reorder(threats, -.value))) +
  tidybayes::stat_slab(data = slab_data,
                       aes(fill = interaction.type),alpha=0.5,normalize = "xy") +
  ggdist::geom_pointinterval(data = diff_int,
                             aes(xmin = .lower, xmax = .upper,col = interaction.type),alpha=0.5) +
  geom_point(data = subset(diff_int,.width == 0.8), 
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
  xlab( expression(paste("Additive ",Delta,"y","/",Delta,"x"," - interactive ",Delta,"y","/",Delta,"x"))) + 
  ylab("Threat combination") +
  theme(legend.position = "none"))

# Figure 3 final ---------------------------------------------------------------

# # First column
# 
# layout <- "
# BB
# AA"
# 
# # First column
# 
#  (g3a <- g3+legend & plot_layout(ncol = 1, nrow=2,
#                                 heights = c(0.05,1), 
#                                 design = layout))

# First row

(row1<- plot_grid(g3a,g3b, labels = "auto"))

# Second row 

(figure3 <- plot_grid(row1, legend, ncol = 1, rel_heights = c(1,0.05)))

(figure3 <- g3a/legend + plot_layout(heights=c(1,.1)))

# Final plot 

ggsave("Results/Figure3.pdf", figure3, 
       width = 8, height = 6)
