# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SupplementaryFigures.R
# - DATE:        07/01/2021
# - DESCRIPTION: Basic diagnostics for the models.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list = ls(all = TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

library(brms)
library(tidyverse)
library(tidybayes)
library(parallel)
library(doSNOW)
library(patchwork)
library(rstan)
library(MetBrewer)

source("Code/utils/prep_data_grid_fn.R")
source("Code/utils/threat_post_draws.R")

# Load the models

m_fr <- read_rds("Results/models/fre_mod.RDS")
m_mar <- read_rds("Results/models/mar_mod.RDS")
m_ter <- read_rds("Results/models/ter_mod.RDS")
m_am <- read_rds("Results/models/amph_mod.RDS")
m_bi <- read_rds("Results/models/bir_mod.RDS")
m_f <- read_rds("Results/models/fish_mod.RDS")
m_ma <- read_rds("Results/models/mam_mod.RDS")
m_rep <- read_rds("Results/models/rep_mod.RDS")

system_ls <- list(m_fr, m_mar, m_ter)
taxon_ls <- list(m_am, m_bi, m_f, m_ma, m_rep)

# Set default ggplot theme

theme_set(
  theme_minimal() +
    theme(
      axis.title.x = element_text(size = 12, margin = margin(
        t = 10,
        r = 0,
        b = 0,
        l = 0
      )),
      axis.title.y = element_text(size = 12, margin = margin(
        t = 0,
        r = 10,
        b = 0,
        l = 0
      )),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 12),
      strip.text.x = element_text(size = 12),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
)

#Create a palette for the threats

threat_palette <- c(met.brewer(name = "Hokusai1", n = 6, type = "continuous"))

# Figure S5: System specific effect of threats to population trends and random effects ##########

# List of threats

all_threats <- c(
  "None",
  "Pollution",
  "Habitat loss",
  "Climate change",
  "Invasive",
  "Exploitation",
  "Disease"
)

# Obtain the coefficients and clean the data

system_coefs <- lapply(system_ls, function(mod) {
  coefs_df <- mod %>%
    gather_draws(`b_scaled_year.*`, regex = TRUE) %>%
    mutate(
      .variable = gsub("b_", "", .variable),
      .variable = gsub("scaled_year:", "", .variable),
      .variable = gsub("invasive", "Invasive", .variable),
      .variable = gsub("habitatl", "Habitat loss", .variable),
      .variable = gsub("climatechange", "Climate change", .variable),
      .variable = gsub("pollution", "Pollution", .variable),
      .variable = gsub("exploitation", "Exploitation", .variable),
      .variable = gsub("disease", "Disease", .variable),
      .variable = gsub("scaled_year", "None", .variable),
      .variable = gsub(" 1", "", .variable)
    )
  
  # Obtain the coeficients for each threat group
  
  coefs_group <- do.call("rbind", lapply(all_threats, function(x) {
    out <- coefs_df %>%
      subset(grepl(x, .variable)) %>% # Filter to the current threat
      mutate(
        int_group = ifelse(grepl("\\.", .variable), "Interactive", "Singular"),
        # Classify threats as 'singular' or 'interactive'
        threat_group = x
      ) %>% # Assign the current threat to a grouping variable
      ungroup()
    return(out)
  }))
  
  # Calculate median and quantile intervals for the rate of change by threat group and interaction type
  
  coefs_interval <- coefs_group  %>%
    group_by(threat_group, int_group) %>%
    ggdist::median_qi(
      .width = c(.95, .8, .5),
      .exclude = c(".chain", ".iteration", ".draw", ".variable")
    ) %>%
    mutate(
      int_group = ifelse(int_group == "Singular", "Singular", "Interactive"),
      int = ordered(int_group, c("Singular", "Interactive"))
    )
  
  # Create a data frame matching threat groups to colors
  
  palette <- data.frame(
    threat_group = unique(coefs_group$threat_group),
    fill_col = c("grey50", threat_palette)
  )
  
  # Prepare the data for plotting by joining with the palette and adjusting factors
  
  plot_coefs <- coefs_group %>%
    left_join(palette, by = "threat_group") %>%
    mutate(fill_col = factor(fill_col),
           int = ordered(int_group, c("Singular", "Interactive")))
  
  (
    g2a <- plot_coefs %>%
      ggplot(aes(
        x = .value, y = reorder(threat_group, .value)
      )) +
      tidybayes::stat_slab(
        data = plot_coefs %>%
          filter(int_group == "Singular"),
        aes(fill = fill_col, group = .variable),
        alpha = 0.5,
        normalize = "groups"
      ) +
      tidybayes::stat_slab(
        data = subset(plot_coefs, int_group == "Interactive"),
        aes(fill = fill_col, group = .variable),
        alpha = 0.5,
        normalize = "panels"
      ) +
      ggdist::geom_pointinterval(
        data = coefs_interval,
        aes(xmin = .lower, xmax = .upper),
        position = position_dodge()
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        colour = "grey50"
      ) +
      labs(x = "Effect on population trends", y = NULL) +
      coord_cartesian(xlim = c(-0.2, 0.2)) +
      scale_x_continuous(breaks = seq(-0.1, 0.1, by = 0.1)) +
      facet_wrap(~ int) +
      ylab("") +
      scale_fill_manual(values = levels(plot_coefs$fill_col), guide = "none") +
      scale_color_manual(values = levels(plot_coefs$fill_col), guide = "none") +
      theme(axis.title.y = element_blank())
  )
  
  return(g2a)
  
})

combinedS5 <-  wrap_plots(system_coefs, ncol = 3) +
  plot_annotation(tag_levels = "a") &
  labs(x = NULL, y = NULL) &
  theme(plot.tag = element_text(face = 'bold'))

# Create a separate plot for the y-axis label

ylabelS5 <- ggplot(data.frame(l = system_coefs[[1]]$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 5, angle = 90) +
  theme_void() +
  coord_cartesian(clip = "off")

# Create a separate plot for the x-axis label

xlabelS5 <- ggplot(data.frame(l = system_coefs[[1]]$labels$x, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 5) +
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the original figure with the y-axis label

(topS5 <- cowplot::plot_grid(ylabelS5, combinedS5, rel_widths = c(1, 25)))

# Combine with the x-axis label

(figureS5 <- cowplot::plot_grid(
  topS5,
  xlabelS5,
  nrow = 2,
  rel_heights = c(25, 1)
))

ggsave("Results/figures/FigureS6.pdf",
       figureS5,
       width = 14,
       height = 4)

# Figure S6: Taxon specific effect of threats to population trends and random effects ##########

# List of threats

all_threats <- c(
  "None",
  "Pollution",
  "Habitat loss",
  "Climate change",
  "Invasive",
  "Exploitation",
  "Disease"
)

# Obtain the coefficients and clean the data

taxon_coefs <- lapply(taxon_ls, function(mod) {
  coefs_df <- mod %>%
    gather_draws(`b_scaled_year.*`, regex = TRUE) %>%
    mutate(
      .variable = gsub("b_", "", .variable),
      .variable = gsub("scaled_year:", "", .variable),
      .variable = gsub("invasive", "Invasive", .variable),
      .variable = gsub("habitatl", "Habitat loss", .variable),
      .variable = gsub("climatechange", "Climate change", .variable),
      .variable = gsub("pollution", "Pollution", .variable),
      .variable = gsub("exploitation", "Exploitation", .variable),
      .variable = gsub("disease", "Disease", .variable),
      .variable = gsub("scaled_year", "None", .variable),
      .variable = gsub(" 1", "", .variable)
    )
  
  # Obtain the coeficients for each threat group
  
  coefs_group <- do.call("rbind", lapply(all_threats, function(x) {
    out <- coefs_df %>%
      subset(grepl(x, .variable)) %>% # Filter to the current threat
      mutate(
        int_group = ifelse(grepl("\\.", .variable), "Interactive", "Singular"),
        # Classify threats as 'singular' or 'interactive'
        threat_group = x
      ) %>% # Assign the current threat to a grouping variable
      ungroup()
    return(out)
  }))
  
  # Calculate median and quantile intervals for the rate of change by threat group and interaction type
  
  coefs_interval <- coefs_group  %>%
    group_by(threat_group, int_group) %>%
    ggdist::median_qi(
      .width = c(.95, .8, .5),
      .exclude = c(".chain", ".iteration", ".draw", ".variable")
    ) %>%
    mutate(
      int_group = ifelse(int_group == "Singular", "Singular", "Interactive"),
      int = ordered(int_group, c("Singular", "Interactive"))
    )
  
  # Create a data frame matching threat groups to colors
  
  palette <- data.frame(
    threat_group = unique(coefs_group$threat_group),
    fill_col = c("grey50", threat_palette)
  )
  
  # Prepare the data for plotting by joining with the palette and adjusting factors
  
  plot_coefs <- coefs_group %>%
    left_join(palette, by = "threat_group") %>%
    mutate(fill_col = factor(fill_col),
           int = ordered(int_group, c("Singular", "Interactive")))
  
  (
    g2a <- plot_coefs %>%
      ggplot(aes(
        x = .value, y = reorder(threat_group, .value)
      )) +
      tidybayes::stat_slab(
        data = plot_coefs %>%
          filter(int_group == "Singular"),
        aes(fill = fill_col, group = .variable),
        alpha = 0.5,
        normalize = "groups"
      ) +
      tidybayes::stat_slab(
        data = subset(plot_coefs, int_group == "Interactive"),
        aes(fill = fill_col, group = .variable),
        alpha = 0.5,
        normalize = "panels"
      ) +
      ggdist::geom_pointinterval(
        data = coefs_interval,
        aes(xmin = .lower, xmax = .upper),
        position = position_dodge()
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        colour = "grey50"
      ) +
      labs(x = "Effect on population trends", y = NULL) +
      coord_cartesian(xlim = c(-0.2, 0.2)) +
      scale_x_continuous(breaks = seq(-0.1, 0.1, by = 0.1)) +
      facet_wrap(~ int) +
      ylab("") +
      scale_fill_manual(values = levels(plot_coefs$fill_col), guide = "none") +
      scale_color_manual(values = levels(plot_coefs$fill_col), guide = "none") +
      theme(axis.title.y = element_blank())
  )
  
  return(g2a)
  
})

combinedS6 <-  wrap_plots(taxon_coefs, ncol = 3) +
  plot_annotation(tag_levels = "a") &
  labs(x = NULL, y = NULL) &
  theme(plot.tag = element_text(face = 'bold'))

# Create a separate plot for the y-axis label

ylabelS6 <- ggplot(data.frame(l = system_coefs[[1]]$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 5, angle = 90) +
  theme_void() +
  coord_cartesian(clip = "off")

# Create a separate plot for the x-axis label

xlabelS6 <- ggplot(data.frame(l = system_coefs[[1]]$labels$x, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 5) +
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the original figure with the y-axis label

(topS6 <- cowplot::plot_grid(ylabelS6, combinedS6, rel_widths = c(1, 25)))

# Combine with the x-axis label

(figureS6 <- cowplot::plot_grid(
  topS6,
  xlabelS6,
  nrow = 2,
  rel_heights = c(25, 1)
))

ggsave("Results/figures/FigureS7.pdf",
       figureS6,
       width = 14,
       height = 5)

# Figure S7: System specific population trends projections --------------------------------------

# The threats to project

threats <- c("pollution",
             "habitatl",
             "climatechange",
             "invasive",
             "exploitation",
             "disease")

system_dydx <- lapply(system_ls, function(mod) {
  threats_col <-  colnames(mod$data)[grepl(paste(threats, collapse = "|"), colnames(mod$data))]
  # Project the population
  
  postdraws <- threat_post_draws(
    model = mod,
    threat_comns = c("none", threats_col),
    ndraws = 1000,
    nuisance = c("series", "SpeciesName"),
    n.cores = 4
  )
  
  # Calculate the derivative
  
  post_proj <- do.call("rbind", lapply(c("none", threats), function(x) {
    out <- postdraws %>%
      subset(grepl(x, threat)) %>% # Filter to the current threat
      reframe(.value = mean(diff(.value) / diff(time)),
              .by = c(threat, .draw)) %>% # Calculate the first derivative (rate of change) for each time series
      ungroup() %>%
      mutate(int_group = ifelse(grepl("\\.", threat), "combined", "single")) %>% # Classify threats as 'single' or 'combined'
      mutate(threat_group = x) %>% # Assign the current threat to a grouping variable
      # Rename threat groups for clarity
      mutate(
        threat_group = gsub("invasive", "Invasive", threat_group),
        threat_group = gsub("habitatl", "Habitat\nloss", threat_group),
        threat_group = gsub("climatechange", "Climate\nchange", threat_group),
        threat_group = gsub("pollution", "Pollution", threat_group),
        threat_group = gsub("exploitation", "Exploitation", threat_group),
        threat_group = gsub("disease", "Disease", threat_group),
        threat_group = gsub("none", "None", threat_group)
      ) %>%
      group_by(threat) %>%
      ungroup()
    return(out)
  }))
  
  # Calculate median and quantile intervals for the rate of change by threat group and interaction type
  
  dydx_interval <- post_proj  %>%
    group_by(threat_group, int_group) %>%
    ggdist::median_qi(
      .width = c(.95, .8, .5),
      .exclude = c(".draw", "threat", "threat_group")
    ) %>%
    mutate(
      int_group = ifelse(int_group == "single", "Singular", "Interactive"),
      int = ordered(int_group, c("Singular", "Interactive"))
    )
  
  # Create a data frame matching threat groups to colors
  
  palette <- data.frame(
    threat_group = unique(post_proj$threat_group),
    fill_col = c("grey50", threat_palette)
  )
  
  # Prepare the data for plotting by joining with the palette and adjusting factors
  
  plot_dydx_threats <- post_proj %>%
    left_join(palette, by = "threat_group") %>%
    mutate(
      fill_col = factor(fill_col),
      int_group = ifelse(int_group == "single", "Singular", "Interactive"),
      int = ordered(int_group, c("Singular", "Interactive"))
    )
  
  (
    g2c <- ggplot(data = plot_dydx_threats, aes(
      y = .value, x = reorder(threat_group, desc(.value))
    )) +
      ggdist::stat_halfeye(
        data = plot_dydx_threats %>%
          filter(int_group == "Singular"),
        aes(fill = fill_col, group = threat_group),
        alpha = 0.5,
        adjust = .5,
        width = .3,
        .width = 0,
        justification = -.5,
        point_colour = NA
      ) +
      geom_boxplot(
        data = plot_dydx_threats %>%
          filter(int_group == "Singular"),
        aes(fill = fill_col, group = threat_group),
        width = .2,
        outliers = F,
        alpha = .7
      ) +
      ggdist::stat_halfeye(
        data = plot_dydx_threats %>%
          filter(int_group == "Interactive"),
        aes(fill = fill_col, group = threat),
        alpha = 0.5,
        adjust = .5,
        width = .3,
        .width = 0,
        justification = -.5,
        point_colour = NA
      ) +
      geom_boxplot(
        data = plot_dydx_threats %>%
          filter(int_group == "Interactive"),
        aes(fill = fill_col, group = threat_group),
        width = .2,
        outliers = F,
        alpha = .7
      ) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "grey50"
      ) +
      facet_wrap(~ int) +
      ylab(expression(
        paste("Population trend (", Delta, "y", "/", Delta, "x)")
      )) +
      xlab("") +
      facet_wrap(~ int, scales = "free") +
      scale_fill_manual(
        values = levels(plot_dydx_threats$fill_col),
        guide = "none"
      ) +
      coord_cartesian(ylim = c(-0.2, 0.2)) +
      theme(
        plot.margin = unit(c(0, 0.5, 0, .5), "cm"),
        axis.title.x = element_blank()
      )
  )
  
})

combinedS7 <-  wrap_plots(system_dydx, nrow = 3) +
  plot_annotation(tag_levels = "a") &
  labs(x = NULL, y = NULL) &
  theme(plot.tag = element_text(face = 'bold'))

# Create a separate plot for the y-axis label

ylabelS7 <- ggplot(data.frame(x = 1, y = 1)) +
  geom_text(
    aes(x, y, label = "Population~trend~(Delta~y/Delta~x)"),
    parse = TRUE,
    size = 5,
    angle = 90
  ) +
  theme_void() +
  coord_cartesian(clip = "off")

# Create a separate plot for the x-axis label

xlabelS7 <- ggplot(data.frame(x = 1, y = 1)) +
  geom_text(aes(x, y, label = system_dydx[[1]]$labels$x), size = 5) +
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the original figure with the y-axis label

(topS7 <- cowplot::plot_grid(ylabelS7, combinedS7, rel_widths = c(1, 25)))

# Combine with the x-axis label

(figureS7 <- cowplot::plot_grid(
  topS7,
  xlabelS7,
  nrow = 2,
  rel_heights = c(25, 1)
))

ggsave("Results/figures/FigureS8.pdf",
       figureS7,
       width = 12,
       height = 8)

# Figure S8: Taxon specific population trends projections --------------------------------------

# The threats to project

threats <- c("pollution",
             "habitatl",
             "climatechange",
             "invasive",
             "exploitation",
             "disease")

taxon_dydx <- lapply(taxon_ls, function(mod) {
  threats_col <-  colnames(mod$data)[grepl(paste(threats, collapse = "|"), colnames(mod$data))]
  # Project the population
  
  postdraws <- threat_post_draws(
    model = mod,
    threat_comns = c("none", threats_col),
    ndraws = 1000,
    nuisance = c("series", "SpeciesName"),
    n.cores = 4
  )
  
  # Calculate the derivative
  
  post_proj <- do.call("rbind", lapply(c("none", threats), function(x) {
    out <- postdraws %>%
      subset(grepl(x, threat)) %>% # Filter to the current threat
      reframe(.value = mean(diff(.value) / diff(time)),
              .by = c(threat, .draw)) %>% # Calculate the first derivative (rate of change) for each time series
      ungroup() %>%
      mutate(int_group = ifelse(grepl("\\.", threat), "combined", "single")) %>% # Classify threats as 'single' or 'combined'
      mutate(threat_group = x) %>% # Assign the current threat to a grouping variable
      # Rename threat groups for clarity
      mutate(
        threat_group = gsub("invasive", "Invasive", threat_group),
        threat_group = gsub("habitatl", "Habitat\nloss", threat_group),
        threat_group = gsub("climatechange", "Climate\nchange", threat_group),
        threat_group = gsub("pollution", "Pollution", threat_group),
        threat_group = gsub("exploitation", "Exploitation", threat_group),
        threat_group = gsub("disease", "Disease", threat_group),
        threat_group = gsub("none", "None", threat_group)
      ) %>%
      group_by(threat) %>%
      ungroup()
    return(out)
  }))
  
  # Calculate median and quantile intervals for the rate of change by threat group and interaction type
  
  dydx_interval <- post_proj  %>%
    group_by(threat_group, int_group) %>%
    ggdist::median_qi(
      .width = c(.95, .8, .5),
      .exclude = c(".draw", "threat", "threat_group")
    ) %>%
    mutate(
      int_group = ifelse(int_group == "single", "Singular", "Interactive"),
      int = ordered(int_group, c("Singular", "Interactive"))
    )
  
  # Create a data frame matching threat groups to colors
  
  palette <- data.frame(
    threat_group = unique(post_proj$threat_group),
    fill_col = c("grey50", threat_palette)
  )
  
  # Prepare the data for plotting by joining with the palette and adjusting factors
  
  plot_dydx_threats <- post_proj %>%
    left_join(palette, by = "threat_group") %>%
    mutate(
      fill_col = factor(fill_col),
      int_group = ifelse(int_group == "single", "Singular", "Interactive"),
      int = ordered(int_group, c("Singular", "Interactive"))
    )
  
  (
    g2c <- ggplot(data = plot_dydx_threats, aes(
      y = .value, x = reorder(threat_group, desc(.value))
    )) +
      ggdist::stat_halfeye(
        data = plot_dydx_threats %>%
          filter(int_group == "Singular"),
        aes(fill = fill_col, group = threat_group),
        alpha = 0.5,
        adjust = .5,
        width = .3,
        .width = 0,
        justification = -.5,
        point_colour = NA
      ) +
      geom_boxplot(
        data = plot_dydx_threats %>%
          filter(int_group == "Singular"),
        aes(fill = fill_col, group = threat_group),
        width = .2,
        outliers = F,
        alpha = .7
      ) +
      ggdist::stat_halfeye(
        data = plot_dydx_threats %>%
          filter(int_group == "Interactive"),
        aes(fill = fill_col, group = threat),
        alpha = 0.5,
        adjust = .5,
        width = .3,
        .width = 0,
        justification = -.5,
        point_colour = NA
      ) +
      geom_boxplot(
        data = plot_dydx_threats %>%
          filter(int_group == "Interactive"),
        aes(fill = fill_col, group = threat_group),
        width = .2,
        outliers = F,
        alpha = .7
      ) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "grey50"
      ) +
      facet_wrap(~ int) +
      ylab(expression(
        paste("Population trend (", Delta, "y", "/", Delta, "x)")
      )) +
      xlab("") +
      facet_wrap(~ int, scales = "free") +
      scale_fill_manual(
        values = levels(plot_dydx_threats$fill_col),
        guide = "none"
      ) +
      coord_cartesian(ylim = c(-0.2, 0.2)) +
      theme(
        plot.margin = unit(c(0, 0.5, 0, .5), "cm"),
        axis.title.x = element_blank()
      )
  )
  
})

combinedS8 <-  wrap_plots(taxon_dydx, nrow = 5) +
  plot_annotation(tag_levels = "a") &
  labs(x = NULL, y = NULL) &
  theme(plot.tag = element_text(face = 'bold'))

# Create a separate plot for the y-axis label

ylabelS8 <- ggplot(data.frame(x = 1, y = 1)) +
  geom_text(
    aes(x, y, label = "Population~trend~(Delta~y/Delta~x)"),
    parse = TRUE,
    size = 5,
    angle = 90
  ) +
  theme_void() +
  coord_cartesian(clip = "off")

# Create a separate plot for the x-axis label

xlabelS8 <- ggplot(data.frame(x = 1, y = 1)) +
  geom_text(aes(x, y, label = system_dydx[[1]]$labels$x), size = 5) +
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the original figure with the y-axis label

(topS8 <- cowplot::plot_grid(ylabelS8, combinedS8, rel_widths = c(1, 25)))

# Combine with the x-axis label

(figureS8 <- cowplot::plot_grid(
  topS8,
  xlabelS8,
  nrow = 2,
  rel_heights = c(25, 1)
))

ggsave("Results/figures/FigureS9.pdf",
       figureS8,
       width = 14,
       height = 12)
