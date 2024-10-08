
# Make the plot

# (g3a <- data %>% 
#     ggplot(aes(x=threats, y=freq/100, 
#                fill=interaction.type)) +
#     geom_bar(stat="identity", position = "dodge", alpha=0.8) +
#     scale_fill_manual(name = NULL,
#                       values = c("#B32315",
#                                  "#1E63B3",
#                                  "#694364"))+
#     scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
#     labs(y="", x="")+
#     theme(
#       axis.text = element_blank(),
#       axis.title = element_blank(),
#       panel.grid = element_blank()))

## Panel b: System -------------------------------------------------------------

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

(g3b <- ggplot(data) +
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

## Panel c: Taxon -------------------------------------------------------------

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

(g3c <- ggplot(data) +
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

(row1 <- plot_grid(g3a, labels="auto"))

#Create second row 

(row2 <- plot_grid(g3b, g3c, labels=c("b", "c")))

# Plot it 

(figure3<- plot_grid(row1, row2, legend,
                     nrow = 3, rel_heights = c(1,1,0.05)))


# Save it

ggsave("Results/Figure 3.pdf", figure3, 
       width = 10, height = 14)





# Figure 2: Interactions, system, and taxon ####################################
# Load the individual models 

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

# Interval 

dydx_interval_system <- post_dyxd_system  %>%
  group_by(System) %>%
  ggdist::median_qi(.width = c(.95, .8, .5), 
                    .exclude = c(".draw", "threat","threat_group", "int_group"))

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

# Interval 

dydx_interval_taxon <- post_dyxd_taxon  %>%
  group_by(Taxon) %>%
  ggdist::median_qi(.width = c(.95, .8, .5), 
                    .exclude = c(".draw", "threat","threat_group", "int_group"))

# Prepare global ---------------------------------------------------------------

# Create a data frame matching threat groups to colors

palette <- data.frame(threat_group = unique(subset(post_dydx, threat_group != "None")$threat_group),
                      fill_col = threat_palette)

# Prepare the data for plotting by joining with the palette and adjusting factors

plot_dydx_threats <- post_dydx %>%
  left_join(palette, by = "threat_group") %>%
  mutate(fill_col = factor(fill_col),
         int_group = ifelse(int_group == "single", "Singular", "Interactive"),
         int = ordered(int_group, c("Singular", "Interactive")))


# Data frame to plot 

plot_global <- post_dydx %>% 
  filter(int_group=="combined"|
           threat_group=="None") %>% 
  left_join(palette) %>% 
  mutate(fill_col = factor(fill_col),
         threat= gsub("pollution", "Polluntion", threat),
         threat= gsub("habitatl", "Habitat loss", threat),
         threat= gsub("climatechange", "Climate change", threat),
         threat= gsub("invasive", "Invasive species", threat),
         threat= gsub("exploitation", "Exploitation", threat),
         threat= gsub("disease", "Disease", threat),
         threat= gsub("none", "None", threat))

# Interval 

dydx_interval <- post_dydx  %>%
  group_by(threat, int_group) %>%
  ggdist::median_qi(.width = c(.95, .8, .5), 
                    .exclude = c(".draw", "threat_group")) %>% 
  mutate(threat= gsub("pollution", "Polluntion", threat),
         threat= gsub("habitatl", "Habitat loss", threat),
         threat= gsub("climatechange", "Climate change", threat),
         threat= gsub("invasive", "Invasive species", threat),
         threat= gsub("exploitation", "Exploitation", threat),
         threat= gsub("disease", "Disease", threat),
         threat= gsub("none", "None", threat))


# Panel a: Threats -------------------------------------------------------------

(p1 <- plot_global %>% 
   ggplot(aes(x = .value,
              y=reorder(threat, .value))) +
   tidybayes::stat_slab(aes(fill = threat,
                            group = threat), 
                        alpha=0.5,
                        normalize = "groups") +
   ggdist::geom_pointinterval(data =dydx_interval %>% 
                                filter(int_group=="combined"|
                                         threat=="None"),
                              aes(xmin = .lower, xmax = .upper),
                              position = position_dodge()) +
   geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
   xlab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
   coord_cartesian(xlim = c(-0.5,0.5)) +
   scale_x_continuous(breaks= seq(-0.5,0.5,by=0.5)) + 
   ylab("")+
   scale_fill_manual(values=met.brewer(name="Nattier", 
                                       n=length(unique(plot_global$threat)), 
                                       type="continuous"))+
   # scale_color_manual(values = palette$fill_col) +
   theme(axis.title.y=element_blank(),
         legend.position = "none"))

# Panel b: System --------------------------------------------------------------

(p2 <- ggplot(data = plot_dydx_system, 
              aes(x = .value,
                  y=reorder(System, 
                            .value))) +
   tidybayes::stat_slab(aes(fill = System,
                            group = System), 
                        alpha=0.5,
                        normalize = "groups") +
   ggdist::geom_pointinterval(data =dydx_interval_system,
                              aes(xmin = .lower, xmax = .upper),
                              position = position_dodge()) +
   geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
   xlab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
   coord_cartesian(xlim = c(-0.5,0.5)) +
   scale_x_continuous(breaks= seq(-0.5,0.5,by=0.5)) + 
   ylab("")+
   scale_fill_manual(values = system_pal) + 
   # scale_color_manual(values = levels(plot_dydx_threats$fill_col),guide = "none")+
   theme(axis.title.y=element_blank(),
         legend.position = "none"))

# Panel c: Taxon --------------------------------------------------------------

(p3 <- ggplot(data = plot_dydx_taxon, 
              aes(x = .value,
                  y=reorder(Taxon, 
                            .value))) +
   tidybayes::stat_slab(aes(fill = Taxon,
                            group = Taxon), 
                        alpha=0.3,
                        normalize = "groups") +
   ggdist::geom_pointinterval(data =dydx_interval_taxon,
                              aes(xmin = .lower, xmax = .upper),
                              position = position_dodge()) +
   geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
   xlab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
   coord_cartesian(xlim = c(-0.5,0.5)) +
   scale_x_continuous(breaks= seq(-0.5,0.5,by=0.5)) + 
   ylab("")+
   scale_fill_manual(values = taxon_pal) + 
   # scale_color_manual(values = levels(plot_dydx_threats$fill_col),guide = "none")+
   theme(axis.title.y=element_blank(),
         legend.position = "none"))

# Final figure 2 ---------------------------------------------------------------

(fig2 <- (p1+p2/p3)+plot_annotation(tag_levels = "a"))

# Save it

ggsave("Results/Figure 2_test.pdf",fig2,
       width = 14, height = 10)

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
# Threat combinations 

# threatcols <- colnames(mod$data)[grepl(paste(all_threats,collapse = "|"),
#                                        colnames(mod$data))] 
# 
# # Estimate posterior draws for all individual threats and their combinations
# 
# postdraws <- threat_post_draws(model = mod,
#                                threat_comns = c("none", threatcols))
# 
# # Combine posterior draws 
# 
# post_dydx <- do.call("rbind", lapply(all_threats, function(x){
#   out <- postdraws %>%
#     subset(grepl(x, threat)) %>% # Filter to the current threat
#     reframe(.value = mean(diff(.value)/diff(time)), 
#             .by = c(threat, .draw)) %>% # Calculate the first derivative (rate of change) for each time series
#     ungroup() %>%
#     mutate(int_group = ifelse(grepl("\\.", threat), "combined", "single")) %>% # Classify threats as 'single' or 'combined'
#     mutate(threat_group = x) %>% # Assign the current threat to a grouping variable
#     # Rename threat groups for clarity
#     mutate(threat_group = gsub("invasive", "Invasive", threat_group),
#            threat_group = gsub("habitatl", "Habitat loss", threat_group),
#            threat_group = gsub("climatechange", "Climate change", threat_group),
#            threat_group = gsub("pollution", "Pollution", threat_group),
#            threat_group = gsub("exploitation", "Exploitation", threat_group),
#            threat_group = gsub("disease", "Disease", threat_group),
#            threat_group = gsub("none", "None", threat_group)) %>% 
#     group_by(threat) %>%
#     ungroup()
#   return(out)
# })) 
# 
# # Calculate median and quantile intervals for the rate of change by threat group and interaction type
# 
# dydx_interval <- post_dydx  %>%
#   group_by(threat_group, int_group) %>%
#   ggdist::median_qi(.width = c(.95, .8, .5), .exclude = c(".draw", "threat")) %>%
#   mutate(int_group = ifelse(int_group == "single", "Singular", "Interactive"),
#          int = ordered(int_group, c("Singular", "Interactive")))
# 
# 
# # Create a data frame matching threat groups to colors
# 
# palette <- data.frame(threat_group = unique(post_dydx$threat_group),
#                       fill_col = c("grey50", threat_palette))
# 
# # Prepare the data for plotting by joining with the palette and adjusting factors
# 
# plot_dydx_threats <- post_dydx %>%
#   filter(threat_group != "None") %>%
#   left_join(palette, by = "threat_group") %>%
#   mutate(fill_col = factor(fill_col),
#          int_group = ifelse(int_group == "single", "Singular", "Interactive"),
#          int = ordered(int_group, c("Singular", "Interactive")))
# Keep only the none threats
# 
# post_dydx_none <- postdraws %>%
#   filter(grepl("none",threat)) %>% #subset to focal threat
#   reframe(.value = mean(diff(.value)/diff(time)),.by = c(threat,.draw)) %>% #for each posterior timeseries, estimate the first derivative 
#   mutate(int_group = ifelse(grepl("\\.",threat),"combined","Singular"), #create grouping column for whether the threat is singular or a combination
#          threat=gsub("none", "None", threat))
# 
# # Calculate intervals
# 
# dydx_none_interval <- post_dydx_none  %>%
#   group_by(int_group, threat) %>%
#   ggdist::median_qi(.width = c(.95, .8, .5),
#                     .exclude = c(".draw", "threat")) #extract distribution information
# 
# # Figure with no threat 
# 
# (p2 <- ggplot(data = post_dydx_none, 
#               aes(x = .value,
#                   y=threat)) +
#     tidybayes::stat_slab(aes(),fill="grey50", alpha=0.3) +
#     ggdist::geom_pointinterval(data = dydx_none_interval,
#                                aes(xmin = .lower, xmax = .upper)) +
#     geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
#     xlab(expression(paste("Population trend (",partialdiff,"y","/",partialdiff,"x)"))) + 
#     ylab("")+
#     coord_cartesian(xlim = c(-0.2,0.2))+
#     theme(plot.margin = unit(c(0,0.5,0,.5),"cm"),
#           axis.title.y=element_blank()))
# 
# # Combine them
# 
# (g1b <- p1 + (p2/patchwork::plot_spacer() +
#                 plot_layout(heights  = c(1,1))) + 
#     plot_layout(widths = c(1,0.8))) 

# Load the individual models 
# 
# mod_amph <- readRDS("Results/models/amph_mod.RDS")
# mod_bir <- readRDS("Results/models/bir_mod.RDS")
# mod_fish <- readRDS("Results/models/fish_mod.RDS")
# mod_mam <- readRDS("Results/models/mam_mod.RDS")
# mod_rep <- readRDS("Results/models/rep_mod.RDS")
# mod_mar <- readRDS("Results/models/mar_mod.RDS")
# mod_fre <- readRDS("Results/models/fre_mod.RDS")
# mod_ter <- readRDS("Results/models/ter_mod.RDS")
# 
# # Prepare system ---------------------------------------------------------------
# # Store system-specific models in a list
# 
# system_model_ls <- list("Marine" = mod_mar,
#                         "Freshwater" = mod_fre,
#                         "Terrestrial" = mod_ter)
# 
# # Calculate population trends and classify whether they are additive, synergistic
# # or antagonistic
# 
# # Define a function to operate on each element of the system_model_ls list
# post_intvadd_system <- lapply(seq_along(system_model_ls), FUN = function(X){
#   
#   # Extract column names from the data within the current element of the list
#   threatcols <- colnames(system_model_ls[[X]]$data)[grepl(paste(all_threats, collapse = "|"), colnames(system_model_ls[[X]]$data))] 
#   
#   # Create additive columns by combining the column names
#   additive_cols <- do.call("c", lapply(strsplit(threatcols, "[.]"), function(x){
#     paste(x, collapse = " + ")
#   })) 
#   
#   # Combine unique threat columns and additive columns
#   target_cols <- unique(c(threatcols, additive_cols))
#   
#   # Generate posterior draws for the system model with specified threats and columns
#   postdraws_system <- threat_post_draws(model = system_model_ls[[X]],
#                                         threat_comns = target_cols,
#                                         ndraws = 1000,
#                                         nuisance = c("series", "SpeciesName"),
#                                         n.cores = 4) %>%
#     mutate(combo_group = case_when(
#       grepl("[.]", threat) ~ "interactive", # If threat contains a dot, classify as interactive
#       grepl("\\+", threat) ~ "additive",    # If threat contains a plus sign, classify as additive
#       TRUE ~ "single"))                     # Otherwise, classify as singular
#   
#   # Calculate derivatives for each threat
#   post_dydx_system <- do.call("rbind", lapply(all_threats, function(x){
#     
#     out <- postdraws_system %>%
#       subset(grepl(x, threat)) %>% # Filter to focal threat
#       reframe(.value = mean(diff(.value) / diff(time)), .by = c(combo_group, threat, .draw)) %>% # Estimate each timeseries' first derivative
#       group_by(threat) %>%
#       ungroup() %>%
#       mutate(threat_group = x) 
#     
#     return(out) # Return the modified data
#   })) 
#   
#   # Calculate differences in derivatives between additive and interactive threats
#   post_system_diff <- do.call("rbind", lapply(additive_cols[grepl("\\+", additive_cols)], function(x){
#     
#     post_dydx_system %>%
#       subset(threat %in% c(x, gsub(" \\+ ", ".", x))) %>% # Subset to shared additive and interactive threats (e.g., "threat1.threat2" and "threat1 + threat2")
#       reframe(.value = diff(c(.value[2], .value[1])), .by = c(threat_group, .draw)) %>% # Find difference in derivatives between additive and interactive threats
#       mutate(threats = gsub(" \\+ ", ".", x) ) # Name threat combination
#     
#   })) %>%
#     mutate(System = names(system_model_ls)[X]) # Add system name to the data
#   
#   # Calculate interval estimates for the differences in derivatives
#   post_interval_system_diff <- post_system_diff %>%
#     na.omit() %>% # Drop missing threats 
#     group_by(System, threat_group, threats) %>%
#     ggdist::median_qi(.width = c(.95, .8, .5), .exclude = c(".draw")) %>%  # Extract distribution information
#     mutate(interaction.type = case_when(
#       .upper < 0 ~ "synergistic",    # If upper bound of the interval is less than 0, classify as synergistic
#       .lower > 0 ~ "antagonistic",   # If lower bound of the interval is greater than 0, classify as antagonistic
#       TRUE ~ "additive"              # Otherwise, classify as additive
#     )) 
#   
#   # Return a list containing derivatives and interval estimates
#   return(list("dydx_diff" = post_system_diff, "interval" = post_interval_system_diff))
# })
# 
# # Substract the differences
# 
# post_mod_system_diff <- do.call("rbind",lapply(post_intvadd_system, `[[`, 'dydx_diff'))
# 
# # Substract the intervals
# 
# post_intvadd_interval_system_diff <- do.call("rbind",lapply(post_intvadd_system, `[[`, 'interval'))
# 
# # Calculate the position 
# 
# pos <- post_mod_system_diff %>% 
#   group_by(threats, System) %>% 
#   summarise(pos=median(.value))
# 
# # Join it
# 
# post_mod_system_diff <- post_mod_system_diff %>% 
#   mutate(threats = gsub("invasive", "Invasive", threats),
#          threats = gsub("habitatl", "Habitat loss", threats),
#          threats = gsub("climatechange", "Climate change", threats),
#          threats = gsub("pollution", "Pollution", threats),
#          threats = gsub("exploitation", "Exploitation", threats),
#          threats = gsub("disease", "Disease", threats),
#          threat_group = gsub("invasive", "Invasive", threats),
#          threat_group = gsub("habitatl", "Habitat loss", threats),
#          threat_group = gsub("climatechange", "Climate change", threats),
#          threat_group = gsub("pollution", "Pollution", threats),
#          threat_group = gsub("exploitation", "Exploitation", threats),
#          threat_group = gsub("disease", "Disease", threats)) 
# 
# post_intvadd_interval_system_diff <- post_intvadd_interval_system_diff %>% 
#   mutate(threats = gsub("invasive", "Invasive", threats),
#          threats = gsub("habitatl", "Habitat loss", threats),
#          threats = gsub("climatechange", "Climate change", threats),
#          threats = gsub("pollution", "Pollution", threats),
#          threats = gsub("exploitation", "Exploitation", threats),
#          threats = gsub("disease", "Disease", threats),
#          threat_group = gsub("invasive", "Invasive", threats),
#          threat_group = gsub("habitatl", "Habitat loss", threats),
#          threat_group = gsub("climatechange", "Climate change", threats),
#          threat_group = gsub("pollution", "Pollution", threats),
#          threat_group = gsub("exploitation", "Exploitation", threats),
#          threat_group = gsub("disease", "Disease", threats)) 
# 
# # Calculate proportions 
# 
# data_ad_syst <- post_intvadd_interval_system_diff %>%
#   separate_rows(threats, sep = "\\.") %>% 
#   group_by(System, threats, interaction.type) %>%
#   summarise(n=n()) %>% 
#   ungroup() %>% 
#   complete(System, threats, interaction.type) %>%
#   group_by(System, threats) %>% 
#   mutate(n=replace_na(n, 0),
#          freq = (n / sum(n))) %>% 
#   ungroup() 
# 
# # Prepare taxon ----------------------------------------------------------------
# # Create a list with all the models 
# 
# taxon_model_ls <- list("Amphibians" = mod_amph,
#                        "Birds" = mod_bir,
#                        "Fish" = mod_fish,
#                        "Mammals" = mod_mam,
#                        "Reptiles" = mod_rep)
# 
# # Calculate population trends and classify whether they are additive, synergistic
# # or antagonistic
# 
# # Define a function to operate on each element of the system_model_ls list
# post_intvadd_taxon <- lapply(seq_along(taxon_model_ls), FUN = function(X){
#   
#   # Extract column names from the data within the current element of the list
#   threatcols <- colnames(taxon_model_ls[[X]]$data)[grepl(paste(all_threats, collapse = "|"), colnames(taxon_model_ls[[X]]$data))] 
#   
#   # Create additive columns by combining the column names
#   additive_cols <- do.call("c", lapply(strsplit(threatcols, "[.]"), function(x){
#     paste(x, collapse = " + ")
#   })) 
#   
#   # Combine unique threat columns and additive columns
#   target_cols <- unique(c(threatcols, additive_cols))
#   
#   # Generate posterior draws for the taxon model with specified threats and columns
#   postdraws_taxon <- threat_post_draws(model = taxon_model_ls[[X]],
#                                        threat_comns = target_cols,
#                                        ndraws = 1000,
#                                        nuisance = c("series", "SpeciesName"),
#                                        n.cores = 4) %>%
#     mutate(combo_group = case_when(
#       grepl("[.]", threat) ~ "interactive", # If threat contains a dot, classify as interactive
#       grepl("\\+", threat) ~ "additive",    # If threat contains a plus sign, classify as additive
#       TRUE ~ "single"))                     # Otherwise, classify as singular
#   
#   # Calculate derivatives for each threat
#   post_dydx_taxon <- do.call("rbind", lapply(all_threats, function(x){
#     
#     out <- postdraws_taxon %>%
#       subset(grepl(x, threat)) %>% # Filter to focal threat
#       reframe(.value = mean(diff(.value) / diff(time)), .by = c(combo_group, threat, .draw)) %>% # Estimate each timeseries' first derivative
#       group_by(threat) %>%
#       ungroup() %>%
#       mutate(threat_group = x) 
#     
#     return(out) # Return the modified data
#   })) 
#   
#   # Calculate differences in derivatives between additive and interactive threats
#   post_taxon_diff <- do.call("rbind", lapply(additive_cols[grepl("\\+", additive_cols)], function(x){
#     
#     post_dydx_taxon %>%
#       subset(threat %in% c(x, gsub(" \\+ ", ".", x))) %>% # Subset to shared additive and interactive threats (e.g., "threat1.threat2" and "threat1 + threat2")
#       reframe(.value = diff(c(.value[2], .value[1])), .by = c(threat_group, .draw)) %>% # Find difference in derivatives between additive and interactive threats
#       mutate(threats = gsub(" \\+ ", ".", x) ) # Name threat combination
#     
#   })) %>%
#     mutate(taxon = names(taxon_model_ls)[X]) # Add taxon name to the data
#   
#   # Calculate interval estimates for the differences in derivatives
#   post_interval_taxon_diff <- post_taxon_diff %>%
#     na.omit() %>% # Drop missing threats 
#     group_by(taxon, threat_group, threats) %>%
#     ggdist::median_qi(.width = c(.95, .8, .5), .exclude = c(".draw")) %>%  # Extract distribution information
#     mutate(interaction.type = case_when(
#       .upper < 0 ~ "synergistic",    # If upper bound of the interval is less than 0, classify as synergistic
#       .lower > 0 ~ "antagonistic",   # If lower bound of the interval is greater than 0, classify as antagonistic
#       TRUE ~ "additive"              # Otherwise, classify as additive
#     )) 
#   
#   # Return a list containing derivatives and interval estimates
#   return(list("dydx_diff" = post_taxon_diff, "interval" = post_interval_taxon_diff))
# })
# 
# # Substract the differences
# 
# post_mod_taxon_diff <- do.call("rbind",lapply(post_intvadd_taxon, `[[`, 'dydx_diff'))
# 
# # Substract the intervals
# 
# post_intvadd_interval_taxon_diff <- do.call("rbind",lapply(post_intvadd_taxon, `[[`, 'interval'))
# 
# # Correct the data frame
# 
# post_mod_taxon_diff <- post_mod_taxon_diff %>% 
#   mutate(threats = gsub("invasive", "Invasive", threats),
#          threats = gsub("habitatl", "Habitat loss", threats),
#          threats = gsub("climatechange", "Climate change", threats),
#          threats = gsub("pollution", "Pollution", threats),
#          threats = gsub("exploitation", "Exploitation", threats),
#          threats = gsub("disease", "Disease", threats),
#          threat_group = gsub("invasive", "Invasive", threats),
#          threat_group = gsub("habitatl", "Habitat loss", threats),
#          threat_group = gsub("climatechange", "Climate change", threats),
#          threat_group = gsub("pollution", "Pollution", threats),
#          threat_group = gsub("exploitation", "Exploitation", threats),
#          threat_group = gsub("disease", "Disease", threats)) 
# 
# post_intvadd_interval_taxon_diff <- post_intvadd_interval_taxon_diff %>% 
#   mutate(threats = gsub("invasive", "Invasive", threats),
#          threats = gsub("habitatl", "Habitat loss", threats),
#          threats = gsub("climatechange", "Climate change", threats),
#          threats = gsub("pollution", "Pollution", threats),
#          threats = gsub("exploitation", "Exploitation", threats),
#          threats = gsub("disease", "Disease", threats),
#          threat_group = gsub("invasive", "Invasive", threats),
#          threat_group = gsub("habitatl", "Habitat loss", threats),
#          threat_group = gsub("climatechange", "Climate change", threats),
#          threat_group = gsub("pollution", "Pollution", threats),
#          threat_group = gsub("exploitation", "Exploitation", threats),
#          threat_group = gsub("disease", "Disease", threats)) 
# 
# # Calculate the proportions 
# 
# data_ad_tax <- post_intvadd_interval_taxon_diff %>% 
#   separate_rows(threats, sep = "\\.") %>% 
#   group_by(taxon, threats, interaction.type) %>%
#   summarise(n=n()) %>% 
#   ungroup() %>% 
#   complete(taxon, threats, interaction.type) %>%
#   group_by(taxon, threats) %>% 
#   mutate(n=replace_na(n, 0),
#          freq = (n / sum(n))) %>% 
#   ungroup() 
# 
# # Save the tables 
# 
# # write.csv2(data_ad_syst, "Results/Table S3.csv")
# # write.csv2(data_ad_tax, "Results/Table S4.csv")
