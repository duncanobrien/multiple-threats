# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Figure1.R
# - DATE:        07/01/2021
# - DESCRIPTION: Figure 1.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list = ls(all = TRUE)) #remove everything

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

# Color palette for taxons

taxon_pal <- wes_palette("Cavalcanti1", n = 5)

# Load the model

mod <- readRDS("Results/models/mod_global_rerun.RDS")

# Figure 1: Map and trends #####################################################
## Map -------------------------------------------------------------------------

# Set the world map

world <- map_data("world")

# Create map

(
  g1a <- ggplot() +
    geom_map(
      map = world,
      data = world,
      aes(long, lat, map_id = region),
      color = "gray80",
      fill = "gray80",
      linewidth = 0.3
    ) +
    #coord_proj("+proj=wintri") +
    theme_map() +
    geom_point(
      data = data,
      aes(x = Longitude, y = Latitude, fill = Taxon),
      alpha = 0.5,
      shape = 21,
      size = 3
    ) +
    scale_y_continuous(limits = c(-80, 80)) +
    scale_fill_manual("", values = taxon_pal) +
    guides(fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(alpha = 1)
    )) +
    expand_limits(x = 0, y = 0) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 12, hjust = 0.5),
      legend.text = element_text(size = 10),
      plot.margin = unit(c(0, -.5, 0, 0), units = , "cm")
    )
)

# Legend

legend1 <- cowplot::get_plot_component(g1a, 'guide-box-top', return_all = TRUE)

# Create distribution 1

anot <- data.frame(
  x = median(data$Duration),
  y = 300,
  label = paste("Median = ", median(data$Duration))
)


(
  gd1 <- data %>%
    ggplot() +
    geom_histogram(
      aes(x = Duration),
      binwidth = 2,
      alpha = .6,
      position = "stack",
      fill = "grey40",
      color = "white"
    ) +
    geom_vline(
      aes(xintercept = median(Duration)),
      linetype = "longdash",
      size = 0.5
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(10, 60, 10), expand = c(0, 0)) +
    ggrepel::geom_text_repel(
      data = anot,
      aes(x = x, y = y + 100, label = label),
      nudge_x = 10,
      nudge_y = 0
    ) +
    labs(x = "Duration (years)", y = "Number of datasets") +
    theme(
      axis.title.x = element_text(margin = margin(t = 0)),
      axis.text.y = element_text(margin = margin(r = 0))
    )
)


# Combined map

(g1a <- ggdraw(g1a + theme(legend.position = c(0.40, 0.08))) +
    draw_plot(
      gd1,
      x = -0.32,
      y = -0.3143,
      scale = 0.35
    ))

## Panel b: threat proportions ------------------------

# Readjust the data frame to contain the threats as individual columns

threat_freq <- mod$data %>%
  distinct(series, .keep_all = T) %>%
  dplyr::select(!dplyr::contains(".")) %>%
  dplyr::select(-c(y_centered, scaled_year, SpeciesName, Site, time)) %>%
  pivot_longer(1:6, names_to = "threats") %>%
  filter(value != 0) %>%
  group_by(threats) %>%
  summarise(n = sum(as.numeric(value))) %>%
  rbind(tibble(
    threats = "None",
    n = mod$data %>%
      distinct(series, .keep_all = TRUE) %>%
      summarise(count = n()) %>%
      pull(count) - sum(.$n)
  )) %>%
  mutate(
    total = mod$data %>%
      distinct(series, .keep_all = TRUE) %>%
      summarise(count = n()) %>%
      pull(count),
    freq = (n / total),
    threats = ifelse(
      threats == "climatechange",
      "Climate change",
      ifelse(
        threats == "exploitation",
        "Exploitation",
        ifelse(
          threats == "invasive",
          "Invasive",
          ifelse(
            threats == "disease",
            "Disease",
            ifelse(
              threats == "habitatl",
              "Habitat loss",
              ifelse(threats ==
                       "pollution", "Pollution", threats)
            )
          )
        )
      )
    )
  )
# We create the palette

palette <- data.frame(
  threats = c(
    "None",
    "Pollution",
    "Habitat loss",
    "Climate change",
    "Invasive",
    "Exploitation",
    "Disease"
  ),
  fill_col = as.factor(c("grey50", threat_palette))
)


# Plot

(
  g1b <- threat_freq %>%
    left_join(palette, by = "threats") %>%
    mutate(threats = factor(
      threats,
      levels = c(
        "None",
        "Disease",
        "Invasive",
        "Climate change",
        "Pollution",
        "Habitat loss",
        "Exploitation"
      )
    )) %>%
    ggplot(aes(
      fill = fill_col, y = freq, x = threats
    )) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = levels(palette$fill_col), guide = NULL) +
    scale_y_continuous(
      breaks = seq(0, 1, .1),
      label = scales::percent,
      expand = c(0, 0)
    ) +
    labs(y = "Proportion of threats (%)", x = "", fill = "") +
    # geom_text(aes(label = paste0(round(freq*100, 0), "%")),
    #           size=5,
    #           position = position_stack(vjust = 0.5)) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      strip.text = element_text(hjust = 0),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      plot.margin = unit(c(0.5, 0, 0, 0), units = "cm")
    )
)

## Final figure 1 --------------------------------------------------------------

(fig1 <- plot_grid(
  g1a,
  g1b,
  nrow = 2,
  rel_heights = c(1, .7),
  labels = "auto"
))

# Save it

ggsave("Results/figures/Figure1.pdf",
       fig1,
       height = 10,
       width = 10)
