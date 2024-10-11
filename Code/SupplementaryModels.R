# --------------------------------------------------------------------------------------- #
# - FILE NAME:   ModelDiagnostics.R
# - DATE:        07/01/2021
# - DESCRIPTION: Basic diagnostics for the models.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com), Duncan O'Brien (duncan.a.obrien@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list = ls(all = TRUE)) #remove everything

# Libraries

library(tidyverse)
library(brms)
library(rstan)
library(geosphere)
library(tidybayes)
library(patchwork)
library(MetBrewer)

source("Code/utils/prepare_data_fn.R")
source("Code/utils/prep_data_grid_fn.R")
source("Code/utils/threat_post_draws.R")

# Function to normalise the spatial distance

norm_range <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Load the raw LPI data

load("Data/LivingPlanetData2.RData")

# Set the priors

priors <- c(prior(normal(0, 1), class = b),
            prior(exponential(1), class = sd),
            prior(normal(0, 0.25), class = ar))

# Fit global model upon time series at least 10 years long -------------------------------------------------------
mod_dat_a <- prepare_data(dd_long, duration = 10, ratio = 0.5) |>
  arrange(ID)

mod_dat_10 <- prepare_data(dd_long, duration = 10) |>
  arrange(ID)

# Subset the locations

locs10 <- distinct(mod_dat_10[, c("Longitude", "Latitude")]) %>%
  mutate(Site = paste(Latitude, Longitude, sep = "_"))

# Get the spatial distance matrix

spa_mat_trim10 <- as.matrix(geosphere::distm(locs10[, c("Longitude", "Latitude")], fun = geosphere::distHaversine)) /
  1000 #km distance between sites
# Normalise the distance

spa_mat_trim10 <- norm_range(spa_mat_trim10)

# Get the absolute value

spa_mat_trim10 <- abs(spa_mat_trim10 - 1)

# Set the col and rownames

colnames(spa_mat_trim10) <- locs10$Site
rownames(spa_mat_trim10) <- locs10$Site

# Create the formula

rhs10 <- paste0(
  paste(
    "scaled_year*",
    colnames(mod_dat_10)[grepl(paste(
      c(
        "pollution",
        "habitatl",
        "climatechange",
        "invasive",
        "exploitation",
        "disease"
      ),
      collapse = "|"
    ), colnames(mod_dat_10))],
    sep = "",
    collapse = " + "
  ),
  " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) +
               (0 + scaled_year|gr(Site, cov = spa_mat)) - 1"
)

# Combine y_centered and rhs into a model formula

form10 <- as.formula(paste("y_centered", "~", rhs10))

# Run the model

mod_glob10 <- brm(
  bf(form10 #include realm/spp as slopes, x intercepts
     , autocor = ~ ar(
       time = time, gr = series, p = 1
     )),
  data = mod_dat_10,
  data2 = list(spa_mat = spa_mat_trim10),
  family = gaussian(),
  iter = 5000,
  refresh = 100,
  backend = "cmdstanr",
  silent = 0,
  prior = priors,
  chains = 4,
  control = list(adapt_delta = 0.975, max_treedepth = 12),
  cores = 4
)

# Save the model

saveRDS(mod_glob10, "Results/models/mod_global_10.RDS")

# Fit global model upon time series at least 25 years long -------------------------------------------------------

mod_dat_25 <- prepare_data(dd_long, duration = 25)

# Subset the locations

locs25 <- distinct(mod_dat_25[, c("Longitude", "Latitude")]) %>%
  mutate(Site = paste(Latitude, Longitude, sep = "_"))

# Get the spatial distance matrix

spa_mat_trim25 <- as.matrix(geosphere::distm(locs25[, c("Longitude", "Latitude")], fun = geosphere::distHaversine)) /
  1000 #km distance between sites
# Normalise the distance

spa_mat_trim25 <- norm_range(spa_mat_trim25)

# Get the absolute value

spa_mat_trim25 <- abs(spa_mat_trim25 - 1)

# Set the col and rownames

colnames(spa_mat_trim25) <- locs25$Site
rownames(spa_mat_trim25) <- locs25$Site

# Create the formula

rhs25 <- paste0(
  paste(
    "scaled_year*",
    colnames(mod_dat_25)[grepl(paste(
      c(
        "pollution",
        "habitatl",
        "climatechange",
        "invasive",
        "exploitation",
        "disease"
      ),
      collapse = "|"
    ), colnames(mod_dat_25))],
    sep = "",
    collapse = " + "
  ),
  " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) +
               (0 + scaled_year|gr(Site, cov = spa_mat)) - 1"
)

# Combine y_centered and rhs into a model formula

form25 <- as.formula(paste("y_centered", "~", rhs25))

# Run the model

mod_glob25 <- brm(
  bf(form25 #include realm/spp as slopes, x intercepts
     , autocor = ~ ar(
       time = time, gr = series, p = 1
     )),
  data = mod_dat_25,
  data2 = list(spa_mat = spa_mat_trim25),
  family = gaussian(),
  iter = 5000,
  refresh = 100,
  backend = "cmdstanr",
  silent = 0,
  prior = priors,
  chains = 4,
  control = list(adapt_delta = 0.975, max_treedepth = 12),
  cores = 4
)

# Save the model

saveRDS(mod_glob25, "Results/models/mod_global_25.RDS")

# Create Figure 2 for both models -------------------------------------------------------

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

#Load sensitivity test models
mod_glob10 <- readRDS("Results/models/mod_global_10.RDS")
mod_glob25 <- readRDS("Results/models/mod_global_25.RDS")

model_ls <- list(mod_glob10, mod_glob25)

sensitivity_fig2 <- lapply(model_ls, function(mod) {
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
  
  coefs_group <- do.call("rbind", lapply(c(
    "None",
    "Pollution",
    "Habitat loss",
    "Climate change",
    "Invasive",
    "Exploitation",
    "Disease"
  ), function(x) {
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
  
  # Prepare the data for Panel A by joining with the palette and adjusting factors
  
  plot_coefs <- coefs_group %>%
    left_join(palette, by = "threat_group") %>%
    mutate(fill_col = factor(fill_col),
           int = ordered(int_group, c("Singular", "Interactive")))
  
  g2a <- plot_coefs %>%
    ggplot(aes(x = .value, y = reorder(threat_group, .value))) +
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
    ggdist::geom_pointinterval(data = coefs_interval,
                               aes(xmin = .lower, xmax = .upper),
                               position = position_dodge()) +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               colour = "grey50") +
    labs(x = "Effect on population trends", y = NULL) +
    coord_cartesian(xlim = c(-0.1, 0.1)) +
    scale_x_continuous(breaks = seq(-0.1, 0.1, by = 0.1)) +
    facet_wrap(~ int) +
    ylab("") +
    scale_fill_manual(values = levels(plot_coefs$fill_col), guide = "none") +
    scale_color_manual(values = levels(plot_coefs$fill_col), guide = "none") +
    theme(axis.title.y = element_blank())
  
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
    mutate(int_present = if_any(threat_cols, ~ .x == "1", TRUE),
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
    reframe(beta = diff(.value),
            .by = c(series, int_present)) %>%
    mutate(type = "Interactive effects")
  
  space_diff <- rbind(base_pred, space_pred) %>%
    reframe(beta = diff(.value),
            .by = c(series, int_present)) %>%
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
  
  # Plot panel b
  g2b <- diff %>%
    ggplot(aes(x = beta, y = type, fill = type)) +
    ggdist::stat_slab(.width = c(0.6, 0.8, 0.95), alpha = .9) +
    ggdist::geom_pointinterval(data = diff_interval,
                               aes(xmin = .lower, xmax = .upper),
                               position = position_dodge()) +
    geom_vline(aes(xintercept = 0),
               linetype = "dashed",
               colour = "grey50") +
    scale_fill_manual("", values = c("#739EF0", "#F0C873")) +
    labs(x = "Trend difference", y = "") +
    #xlim(-.3,.3)+
    coord_cartesian(xlim = c(-0.1, 0.1)) +
    scale_x_continuous(breaks = seq(-0.2, 0.2, by = 0.1)) +
    theme(legend.position = "none")
  
  # The threats to project
  
  threats_col <-  colnames(mod$data)[grepl(paste(
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
  # Project the population
  
  postdraws <- threat_post_draws(
    model = mod,
    threat_comns = c("none", threats_col),
    ndraws = 1000,
    nuisance = c("series", "SpeciesName"),
    n.cores = 4
  )
  
  # Calculate the derivative
  
  post_proj <- do.call("rbind", lapply(c(
    "none",
    "pollution",
    "habitatl",
    "climatechange",
    "invasive",
    "exploitation",
    "disease"
  ), function(x) {
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
  
  # Prepare the data for Panel C by joining with the palette and adjusting factors
  
  plot_dydx_threats <- post_proj %>%
    left_join(palette, by = "threat_group") %>%
    mutate(
      fill_col = factor(fill_col),
      int_group = ifelse(int_group == "single", "Singular", "Interactive"),
      int = ordered(int_group, c("Singular", "Interactive"))
    )
  
  g2c <- ggplot(data = plot_dydx_threats, aes(y = .value, x = reorder(threat_group, desc(.value)))) +
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
    geom_hline(yintercept = 0,
               linetype = "dashed",
               colour = "grey50") +
    facet_wrap( ~ int) +
    ylab(expression(paste(
      "Population trend (", Delta, "y", "/", Delta, "x)"
    ))) +
    xlab("") +
    facet_wrap( ~ int, scales = "free") +
    scale_fill_manual(values = levels(plot_dydx_threats$fill_col),
                      guide = "none") +
    coord_cartesian(ylim = c(-0.1, 0.1)) +
    theme(plot.margin = unit(c(0, 0.5, 0, .5), "cm"),
          axis.title.x = element_blank())
  
  layout <- "
  AAABB
  CCCCC
  "
  
  fig2 <- g2a + g2b + free(g2c) +
    plot_layout(design = layout, widths = 1) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))
  
  return(fig2)
})

sensitivity_fig2[[1]] + sensitivity_fig2[[2]]
figureS8 <-  wrap_plots(unlist(sensitivity_fig2), nrow = 2)
