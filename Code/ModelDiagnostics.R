# --------------------------------------------------------------------------------------- #
# - FILE NAME:   ModelDiagnostics.R
# - DATE:        07/01/2021
# - DESCRIPTION: Basic diagnostics for the models.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com), Duncan O'Brien (duncan.a.obrien@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list = ls(all = TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

library(tidyverse)
library(brms)
library(rstan)
library(bayesplot)
library(patchwork)
library(performance)
library(cowplot)

# Load the models

m1 <- read_rds("Results/models/mod_global_rerun.RDS")
m_fr <- read_rds("Results/models/fre_mod.RDS")
m_mar <- read_rds("Results/models/mar_mod.RDS")
m_ter <- read_rds("Results/models/ter_mod.RDS")
m_am <- read_rds("Results/models/amph_mod.RDS")
m_bi <- read_rds("Results/models/bir_mod.RDS")
m_f <- read_rds("Results/models/fish_mod.RDS")
m_ma <- read_rds("Results/models/mam_mod.RDS")
m_rep <- read_rds("Results/models/rep_mod.RDS")

mod_ls <- list(m1, m_fr, m_mar, m_ter, m_am, m_bi, m_f, m_ma, m_rep)

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
      legend.text = element_text(size = 12)
    )
)

# Normality of residuals -------------------------------------------------------
# Global

p3a <- lapply(mod_ls, function(m) {
  m$data %>%
    mutate(std_resid = residuals(m)[, "Estimate"]) %>%
    ggplot(aes(std_resid)) +
    geom_histogram(aes(y = after_stat(density)),
                   colour = "#9B9B9B",
                   fill = "#BFBFBF") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    labs(x = "Standardised residuals", y = "Density")
})

combined_plotS2 <-  wrap_plots(p3a, nrow = 3) +
  plot_annotation(tag_levels = "a") &
  labs(x = NULL, y = NULL) &
  theme(plot.tag = element_text(face = 'bold'))

# Create a separate plot for the y-axis label

ylabelS2 <- ggplot(data.frame(l = p3a[[1]]$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 5, angle = 90) +
  theme_void() +
  coord_cartesian(clip = "off")

# Create a separate plot for the x-axis label

xlabelS2 <- ggplot(data.frame(l = p3a[[1]]$labels$x, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 5) +
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the original figure with the y-axis label

(topS2 <- cowplot::plot_grid(ylabelS2, combined_plotS2, rel_widths = c(1, 25)))

# Combine with the x-axis label

(figureS2 <- cowplot::plot_grid(
  topS2,
  xlabelS2,
  nrow = 2,
  rel_heights = c(25, 1)
))

ggsave("Results/figures/FigureS2.pdf",
       figureS2,
       width = 12,
       height = 8)

# Homocedasticity --------------------------------------------------------------

# General

p4a <- lapply(mod_ls, function(m) {
  m$data %>%
    mutate(predict_y = predict(m)[, "Estimate"], std_resid = residuals(m)[, "Estimate"]) %>%
    ggplot(aes(predict_y, std_resid)) +
    geom_point(size = 1.5,
               shape = 21,
               fill = "#BFBFBF") +
    stat_smooth(se = FALSE, colour = "#9B9B9B") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    labs(y = "Standardised residuals", x = "Predicted values")
})

combined_plotS3 <-  wrap_plots(p4a, nrow = 3) +
  plot_annotation(tag_levels = "a") &
  labs(x = NULL, y = NULL) &
  theme(plot.tag = element_text(face = 'bold'))

# Create a separate plot for the y-axis label

ylabelS3 <- ggplot(data.frame(l = p4a[[1]]$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 5, angle = 90) +
  theme_void() +
  coord_cartesian(clip = "off")

# Create a separate plot for the x-axis label

xlabelS3 <- ggplot(data.frame(l = p4a[[1]]$labels$x, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 5) +
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the original figure with the y-axis label

(topS3 <- cowplot::plot_grid(ylabelS3, combined_plotS3, rel_widths = c(1, 25)))

# Combine with the x-axis label

(figureS3 <- cowplot::plot_grid(
  topS3,
  xlabelS3,
  nrow = 2,
  rel_heights = c(25, 1)
))

ggsave("Results/figures/FigureS3.pdf",
       figureS3,
       width = 12,
       height = 8)

# Posterior predictive check ---------------------------------------------------

color_scheme_set("darkgray")

p5a <- lapply(mod_ls, function(m) {
  pp_check(m, ndraws = 100) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits = c(-5, 5))
})

figureS4 <-  wrap_plots(p5a, nrow = 3) +
  plot_annotation(tag_levels = "a") &
  #labs(x = NULL, y = NULL) &
  theme(plot.tag = element_text(face = 'bold'))

# Save the combined figure

ggsave("Results/figures/FigureS4.pdf",
       figureS4,
       height = 8,
       width = 10)

# Autocorrelation check ---------------------------------------------------

p6a <- lapply(mod_ls, function(m) {
  # m$data %>%
  #   mutate(std_resid = residuals(m)[, "Estimate"]) %>%
  #   reframe(acf = acf(std_resid, plot = FALSE)$acf, .by = series) %>%
  #   mutate(lag = factor(seq_along(acf) - 1), .by = series) %>%
  #   ggplot(aes(x = lag)) +
  #   geom_point(aes(y = acf),
  #              alpha = 0.5,
  #              colour = "#9B9B9B",
  #              position = position_jitter()) +
  #   geom_boxplot(aes(y = acf), fill = "#BFBFBF", alpha = 0.5) +
  #   coord_cartesian(xlim = c(2, 18)) +
  #   labs(x = "Lag", y = "ACF")
  
  m$data %>%
    mutate(std_resid = residuals(m)[, "Estimate"]) %>%
    reframe(acf = acf(std_resid, plot = FALSE)$acf) %>%
    mutate(lag = seq_along(acf) - 1) %>%
    ggplot(aes(x = lag)) +
    geom_point(aes(y = acf),
               colour = "#9B9B9B") +
    geom_hline(yintercept = c(-0.1,0.1),linetype = "dashed", color = "black") +
    labs(x = "Lag", y = "ACF")
})

figureS5 <-  wrap_plots(p6a, nrow = 3) +
  plot_annotation(tag_levels = "a") &
  #labs(x = NULL, y = NULL) &
  theme(plot.tag = element_text(face = 'bold'))

# Save the combined figure

ggsave("Results/figures/FigureS5.pdf",
       figureS5,
       height = 8,
       width = 10)

# Multicollinearity ------------------------------------------------------------

check_collinearity(m1)
