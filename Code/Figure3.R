# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Figure3.R
# - DATE:        07/01/2021
# - DESCRIPTION: Figure 2.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com), Duncan O'Brien (duncan.a.obrien@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list = ls(all = TRUE)) #remove everything

# Libraries

library(brms)
library(tidyverse)
library(tidybayes)
library(parallel)
library(doSNOW)
library(patchwork)
library(rstan)
library(waffle)

source("Code/utils/prep_data_grid_fn.R")
source("Code/utils/threat_post_draws.R")

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

# Load the models

mod <- readRDS("Results/models/mod_global_rerun.RDS")

# Figure 3: Interactive threats ################################################
# Prepare data -----------------------------------------------------------------
# Define a vector of threat categories

all_threats <- c("pollution",
                 "habitatl",
                 "climatechange",
                 "invasive",
                 "exploitation",
                 "disease")

# Extract column names from the model

threatcols <- colnames(mod$data)[grepl(paste(all_threats, collapse = "|"), colnames(mod$data))]

# Generate additive combinations of threat columns (i.e., "threat1 + threat2")

additive_cols <- do.call("c", lapply(strsplit(threatcols, "[.]"), function(x) {
  paste(x, collapse = " + ")
}))

# Combine original and additive columns into a unique set

target_cols <- unique(c(threatcols, additive_cols))

# Estimate posterior time series for each threat combination (singular, interactive, additive)

postdraws_intvadd <- threat_post_draws(
  model = mod,
  threat_comns = target_cols,
  nuisance = c("series", "SpeciesName"),
  n.cores = 4,
  ndraws = 1000
) %>%
  mutate(combo_group = case_when(
    grepl("[.]", threat) ~ "interactive",
    grepl("\\+", threat) ~ "additive",
    TRUE ~ "single"
  ))

# Calculate the first derivative of the value over time for each threat and categorise them

post_dydx_intvadd <- do.call("rbind", lapply(all_threats, function(x) {
  out <- postdraws_intvadd %>%
    subset(grepl(x, threat)) %>% #filter to focal threat
    reframe(.value = mean(diff(.value) / diff(time)),
            .by = c(combo_group, threat, .draw)) %>% #estimate each timeseries' first derivative
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
    reframe(.value = diff(c(.value[2], .value[1])), .by = c(threat_group, .draw)) %>%
    mutate(threats = gsub(" \\+ ", ".", x))
}))

# Finalize the comparison, drop missing values, and classify interaction types

post_interval_intvadd_diff <- post_intvadd_diff %>%
  na.omit() %>%
  group_by(threat_group, threats) %>%
  ggdist::median_qi(.width = c(.95, .8, .5), .exclude = c(".draw")) %>%
  mutate(
    interaction.type = case_when(
      .upper < 0 ~ "synergistic",
      .lower > 0 ~ "antagonistic",
      TRUE ~ "additive"
    )
  )

# Calculate the proportions

data_ad <- post_interval_intvadd_diff %>%
  separate_rows(threats, sep = "\\.") %>%
  mutate(interaction.type = factor(
    interaction.type,
    levels = c("synergy", "additive", "antagonistic")
  )) %>%
  group_by(threats, interaction.type) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  complete(threats, interaction.type) %>%
  group_by(threats) %>%
  mutate(n = replace_na(n, 0), freq = (n / sum(n))) %>%
  ungroup()

## Panel a: Waffle plot proportions --------------------------------------------

# We will have to create a separate legend because not all levels are shown given
# that no synergies were found

# Create the base data

data <- data_ad %>%
  mutate(
    freq = freq * 100,
    interaction.type = gsub("synergy", "Synergy", interaction.type),
    interaction.type = gsub("additive", "Additive", interaction.type),
    interaction.type = gsub("antagonistic", "Antagonistic", interaction.type),
    interaction.type = fct_relevel(interaction.type, c("Synergy", "Antagonistic", "Additive")),
    threats = gsub("invasive", "Invasive", threats),
    threats = gsub("habitatl", "Habitat loss", threats),
    threats = gsub("climatechange", "Climate change", threats),
    threats = gsub("pollution", "Pollution", threats),
    threats = gsub("exploitation", "Exploitation", threats),
    threats = gsub("disease", "Disease", threats),
    threats = as.factor(threats)
  )

# Create a fake data with synergies and plot it

legend <- cowplot::get_plot_component(
  data %>%
    mutate(freq = ifelse(freq == 0, 1, freq)) %>%
    ggplot(aes(fill = interaction.type, values =
                 freq)) +
    geom_waffle(
      color = "white",
      size = .25,
      n_rows = 10,
      flip = TRUE,
      na.rm = F,
      make_proportional = T,
      show.legend = TRUE
    ) +
    facet_wrap( ~ threats, nrow = 2, strip.position = "bottom") +
    scale_x_discrete() +
    scale_y_continuous(
      labels = function(x)
        x * 10,
      # make this multiplyer the same as n_rows
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      name = NULL,
      values = c("#B32315", "#1E63B3", "#694364")
    ) +
    coord_equal(clip = "off") +
    theme_void() +
    theme(
      legend.position = "top",
      legend.text = element_text(size =
                                   12),
      strip.text = element_text(size =
                                  14),
      plot.margin = margin(-10, 0, -10, 0)
    ),
  'guide-box-top',
  return_all = TRUE
)

# Make the plot

(
  g3a <- data %>%
    ggplot(aes(fill = interaction.type, values = freq)) +
    geom_waffle(
      color = "white",
      size = .25,
      n_rows = 10,
      flip = TRUE,
      na.rm = F,
      make_proportional = T,
      show.legend = TRUE
    ) +
    facet_wrap( ~ threats, nrow = 2, strip.position = "bottom") +
    scale_x_discrete() +
    scale_y_continuous(
      labels = function(x)
        x * 10,
      # make this multiplyer the same as n_rows
      expand = c(0, 0)
    ) +
    scale_fill_manual(name = NULL, values = c("#1E63B3", "#694364")) +
    coord_equal(clip = "off") +
    theme_void() +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 12),
      strip.text = element_text(size = 14),
      plot.margin = margin(0, 0, 0, 0)
    )
)

(figure3 <- g3a / legend + plot_layout(heights = c(1, .1)))

# Final plot

ggsave("Results/figures/Figure3.pdf",
       figure3,
       width = 8,
       height = 6)
