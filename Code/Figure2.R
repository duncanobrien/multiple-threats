# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Figure2.R
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
library(MetBrewer)

source("Code/utils/prep_data_grid_fn.R")
source("Code/utils/threat_post_draws.R")

set.seed(43)

# Load the model

mod <- readRDS("Results/models/mod_global_rerun.RDS")

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

# Figure 2: Effect of threats to population trends and random effects ##########
## Panel a: Coefficients -------------------------------------------------------
## Prepare data

# Obtain the coefficients and clean the data

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
  mutate(fill_col = factor(fill_col), int = ordered(int_group, c("Singular", "Interactive")))

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
    facet_wrap( ~ int) +
    ylab("") +
    scale_fill_manual(values = levels(plot_coefs$fill_col), guide = "none") +
    scale_color_manual(values = levels(plot_coefs$fill_col), guide = "none") +
    theme(axis.title.y = element_blank())
)

## Panel b: Population trends explained by interactions vs random effects ------

# Identify the colnames for the threats

threat_cols <-  colnames(mod$data)[grepl(paste(
  c(
    "pollution",
    "habitatl",
    "climatechange",
    "invasive",
    "exploitation",
    "disease"
  ),
  collapse = "|"
), colnames(mod$data))]

# Get the interaction data

int_dat <- mod$data %>%
  mutate(int_present = if_any(threat_cols,  ~ .x == "1", TRUE),
         y_centered = NA)

# Generate the null data

null_dat <- int_dat %>%
  mutate(across(all_of(threat_cols), ~ factor("0", levels = c("0", "1"))))

# Base predictions

base_pred <- brms::posterior_epred(
  mod,
  newdata = null_dat,
  re.form = NA,
  incl_autocor = F,
  sort = TRUE,
  ndraws = 1000
) %>% #extract posterior draws for the above data grid
  as.data.frame() %>%
  mutate(.draw = 1:NROW(.)) %>%
  #extract posterior draws for the data used to create the model
  pivot_longer(-.draw, names_to = "index", values_to = ".value") %>%
  cbind(null_dat %>% dplyr::select(series, time, int_present),
        row.names = NULL) %>%
  # reframe(.value = mean(.value/lag(.value)),
  #         .by = c(series,.draw,int_present)) %>% #estimate derivatives per series and .draw
  reframe(.value = mean(diff(.value) / diff(time)),
          .by = c(series, .draw, int_present)) %>% #estimate derivatives per series and .draw
  reframe(.value = mean(.value),
          .by = c(series, int_present)) #estimate average derivative per series

# Random effects

space_pred <- brms::posterior_epred(
  mod,
  newdata = null_dat,
  re.form = NULL,
  incl_autocor = F,
  sort = TRUE,
  ndraws = 1000
) %>% #extract posterior draws for the above data grid
  as.data.frame() %>%
  mutate(.draw = 1:NROW(.)) %>%
  #extract posterior draws for the data used to create the model
  pivot_longer(-.draw, names_to = "index", values_to = ".value") %>%
  cbind(null_dat %>% dplyr::select(series, time, int_present),
        row.names = NULL) %>%
  # reframe(.value = mean(.value/lag(.value)),
  #         .by = c(series,.draw,int_present)) %>% #estimate derivatives per series and .draw
  reframe(.value = mean(diff(.value) / diff(time)),
          .by = c(series, .draw, int_present)) %>% #estimate derivatives per series and .draw
  reframe(.value = mean(.value),
          .by = c(series, int_present)) #estimate average derivative per series

# Predictions with interactive effects on

inter_pred <- brms::posterior_epred(
  mod,
  newdata = int_dat,
  re.form = NA,
  incl_autocor = F,
  sort = TRUE,
  ndraws = 1000
) %>% #extract posterior draws for the above data grid
  as.data.frame() %>%
  mutate(.draw = 1:NROW(.)) %>%
  #extract posterior draws for the data used to create the model
  pivot_longer(-.draw, names_to = "index", values_to = ".value") %>%
  cbind(null_dat %>% dplyr::select(series, time, int_present),
        row.names = NULL) %>%
  # reframe(.value = mean(.value/lag(.value)),
  #        .by = c(series,.draw,int_present)) %>% #estimate derivatives per series and .draw
  reframe(.value = mean(diff(.value) / diff(time)),
          .by = c(series, .draw, int_present)) %>% #estimate derivatives per series and .draw
  reframe(.value = mean(.value),
          .by = c(series, int_present)) #estimate average derivative per series

# Calculate the differences

inter_diff <- rbind(base_pred, inter_pred) %>%
  reframe(beta = diff(.value), .by = c(series, int_present)) %>%
  mutate(type = "Interactive effects")

space_diff <- rbind(base_pred, space_pred) %>%
  reframe(beta = diff(.value), .by = c(series, int_present)) %>%
  mutate(type = "Random effects")

# Join them

diff <- rbind(inter_diff, space_diff) %>%
  filter(int_present == TRUE)

# Calculate median and quantile intervals for the rate of change by threat group and interaction type

diff_interval <- diff  %>%
  group_by(type) %>%
  ggdist::median_qi(
    .width = c(.95, .8, .5),
    .exclude = c("series", "int_present", "type")
  )

# Plot it

(
  g2b <- diff %>%
    ggplot(aes(
      x = beta, y = type, fill = type
    )) +
    ggdist::stat_slab(.width = c(0.6, 0.8, 0.95), alpha = .9) +
    ggdist::geom_pointinterval(
      data = diff_interval,
      aes(xmin = .lower, xmax = .upper),
      position = position_dodge()
    ) +
    geom_vline(
      aes(xintercept = 0),
      linetype = "dashed",
      colour = "grey50"
    ) +
    scale_fill_manual("", values = c("#739EF0", "#F0C873")) +
    labs(x = "Trend difference", y = "") +
    #xlim(-.3,.3)+
    coord_cartesian(xlim = c(-0.2, 0.2)) +
    scale_x_continuous(breaks = seq(-0.2, 0.2, by = 0.1)) +
    theme(legend.position = "none")
)

## Panel c: Population trends projections --------------------------------------

# The threats to project

threats <- c("pollution",
             "habitatl",
             "climatechange",
             "invasive",
             "exploitation",
             "disease")

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
    facet_wrap( ~ int) +
    ylab(expression(
      paste("Population trend (", Delta, "y", "/", Delta, "x)")
    )) +
    xlab("") +
    facet_wrap( ~ int, scales = "free") +
    scale_fill_manual(
      values = levels(plot_dydx_threats$fill_col),
      guide = "none"
    ) +
    coord_cartesian(ylim = c(-0.2, 0.2)) +
    theme(plot.margin = unit(c(0, 0.5, 0, .5), "cm"), axis.title.x = element_blank())
)


## Final figure 2 ------------------------------------------------------------

layout <- "
AAABB
CCCCC
"

(
  fig2 <- g2a + g2b + free(g2c) +
    plot_layout(design = layout, widths = 1) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))
)

# Save it

ggsave("Results/figures/Figure2.pdf",
       fig2,
       height = 8,
       width = 14)