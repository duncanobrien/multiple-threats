# --------------------------------------------------------------------------------------- #
# - FILE NAME:   ModelDiagnostics.R         
# - DATE:        07/01/2021
# - DESCRIPTION: Basic diagnostics for the models. 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

library(tidyverse)
library(brms)
library(rstan)
library(data.table)
library(dplyr)
library(tidybayes)
library(bayesplot)
library(loo)
library(rstan)
library(patchwork)
library(performance)

# Set working directory

setwd(gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path)))

# Load the model 

m1 <- read_rds("Results/models/mod_global_rerun.RDS")

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
                  legend.text = element_text(size = 12)))

# Model convergence ------------------------------------------------------------

# System rhat

rhats <- rhat(ms1)

mcmc_rhat(rhats)

# Taxa rhat

rhats <- rhat(mt1)
mcmc_rhat(rhats, size = 2)

# Neff for system

ratios <- neff_ratio(ms1)
mcmc_neff(ratios)

# Neff for Taxa

ratios <- neff_ratio(mt1)
mcmc_neff(ratios)

# Normality of residuals -------------------------------------------------------
# Global 

(p1a <- m1$data %>% 
   mutate(std_resid = residuals(m1)[ , "Estimate"]) %>% 
   ggplot(aes(std_resid)) + 
   geom_histogram(aes(y=..density..),
                  colour="#9B9B9B", fill="#BFBFBF")+ 
   scale_y_continuous(expand = c(0,0))+
   scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
   labs(x="Standradised residuals", y="Density"))

# System 

(p2a <- ms1$data %>% 
  mutate(std_resid = residuals(ms1)[ , "Estimate"]) %>% 
  ggplot(aes(std_resid)) + 
  geom_histogram(aes(y=..density..),
                 colour="#9B9B9B", fill="#BFBFBF")+ 
   scale_y_continuous(expand = c(0,0))+
   scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
   labs(x="Standradised residuals", y="Density"))
  
# Taxon

(p3a <- mt1$data %>% 
    mutate(std_resid = residuals(mt1)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# Interactions freshwater

(p4a <- mmu_fr$data %>% 
    mutate(std_resid = residuals(mmu_fr)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# Interactions marine

(p5a <- mmu_mar$data %>% 
    mutate(std_resid = residuals(mmu_mar)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# Interactions terrestrial

(p6a <- mmu_ter$data %>% 
    mutate(std_resid = residuals(mmu_ter)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# Interactions amphibians

(p7a <- mmu_am$data %>% 
    mutate(std_resid = residuals(mmu_am)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# Interactions birds

(p8a <- mmu_bi$data %>% 
    mutate(std_resid = residuals(mmu_bi)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# Interactions fishes

(p9a <- mmu_f$data %>% 
    mutate(std_resid = residuals(mmu_f)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# Interactions mammals

(p10a <- mmu_ma$data %>% 
    mutate(std_resid = residuals(mmu_ma)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# Interactions reptiles

(p11a <- mmu_rep$data %>% 
    mutate(std_resid = residuals(mmu_rep)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# Combine the individual plots 

combined_plot <- (p1a + p2a + p3a + p4a +
                    p5a + p6a + p7a + p8a +
                    p9a + p10a + p11a & 
  labs(x = NULL, y = NULL)) + 
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 3) &
  theme(plot.tag = element_text(face = 'bold'))

# Create a separate plot for the y-axis label 

ylabel <- ggplot(data.frame(l = p1a$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size=5, angle = 90) + 
  theme_void() +
  coord_cartesian(clip = "off")

# Create a separate plot for the x-axis label 

xlabel <- ggplot(data.frame(l = p1a$labels$x, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size=5) + 
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the original figure with the y-axis label 

(top <- cowplot::plot_grid(ylabel, combined_plot,
                                rel_widths = c(1, 25)))

# Combine with the x-axis label 

(figureS3 <- cowplot::plot_grid(top, xlabel, nrow = 2,
                                rel_heights = c(25, 1)))

# Save the final plot 

ggsave("Figure S3.pdf",path = ResultPath, figureS3, width = 12, height = 8)

# Homocedasticity --------------------------------------------------------------

# General

(p1b <- m1$data %>% 
   mutate(predict_y = predict(m1)[ , "Estimate"], 
          std_resid = residuals(m1)[ , "Estimate"]) %>% 
   ggplot(aes(predict_y, std_resid)) + 
   geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
   stat_smooth(se = FALSE,colour="#9B9B9B") +
   scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
   scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
   labs(y="Standardised residuals", 
        x="Predicted values"))

# System 

(p2b <- ms1$data %>% 
   mutate(predict_y = predict(ms1)[ , "Estimate"], 
          std_resid = residuals(ms1)[ , "Estimate"]) %>% 
   ggplot(aes(predict_y, std_resid)) + 
   geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
   stat_smooth(se = FALSE,colour="#9B9B9B") +
   scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
   labs(y="Standardised residuals", 
        x="Predicted values"))

# Taxon

(p3b <- mt1$data %>% 
  mutate(predict_y = predict(mt1)[ , "Estimate"], 
         std_resid = residuals(mt1)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# Interactions freshwater

(p4b <- mmu_fr$data %>% 
    mutate(predict_y = predict(mmu_fr)[ , "Estimate"], 
           std_resid = residuals(mmu_fr)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# Interactions marine

(p5b <- mmu_mar$data %>% 
    mutate(predict_y = predict(mmu_mar)[ , "Estimate"], 
           std_resid = residuals(mmu_mar)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# Interactions terrestrial

(p6b <- mmu_ter$data %>% 
    mutate(predict_y = predict(mmu_ter)[ , "Estimate"], 
           std_resid = residuals(mmu_ter)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# Interactions amphibians

(p7b <- mmu_am$data %>% 
    mutate(predict_y = predict(mmu_am)[ , "Estimate"], 
           std_resid = residuals(mmu_am)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# Interactions birds

(p8b <- mmu_bi$data %>% 
    mutate(predict_y = predict(mmu_bi)[ , "Estimate"], 
           std_resid = residuals(mmu_bi)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# Interactions fishes

(p9b <- mmu_f$data %>% 
    mutate(predict_y = predict(mmu_f)[ , "Estimate"], 
           std_resid = residuals(mmu_f)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# Interactions mammals

(p10b <- mmu_ma$data %>% 
    mutate(predict_y = predict(mmu_ma)[ , "Estimate"], 
           std_resid = residuals(mmu_ma)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# Interactions reptiles

(p11b <- mmu_rep$data %>% 
    mutate(predict_y = predict(mmu_rep)[ , "Estimate"], 
           std_resid = residuals(mmu_rep)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# Combine the plots into a single figure

combined_plot <-  (p1b+p2b+p3b+p4b+
                p5b+p6b+p7b+p8b+
                p9b+p10b+p11b & labs(x = NULL, y = NULL)) +
  plot_annotation(tag_levels = "a")+
  plot_layout(nrow = 3) &
  theme(plot.tag = element_text(face = 'bold'))

# Create a separate plot for the y-axis label 

ylabel <- ggplot(data.frame(l = p1b$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size=5, angle = 90) + 
  theme_void() +
  coord_cartesian(clip = "off")

# Create a separate plot for the x-axis label 

xlabel <- ggplot(data.frame(l = p1b$labels$x, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size=5) + 
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the original figure with the y-axis label 

(top <- cowplot::plot_grid(ylabel, combined_plot,
                           rel_widths = c(1, 25)))

# Combine with the x-axis label 

(figureS4 <- cowplot::plot_grid(top, xlabel, nrow = 2,
                                rel_heights = c(25, 1)))

# Save the combined figure 

ggsave("Figure S4.pdf", figureS4, 
       path = ResultPath,
       height = 8, width = 12)

# Posterior predictive check ---------------------------------------------------

color_scheme_set("darkgray")

# General 

(p1c <- pp_check(m1, ndraws = 100)+ 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

# System 

(p2c <- pp_check(ms1, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

(p2d <- pp_check(ms1, type = "stat_grouped", stat = "mean", group = "System")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)))


# Taxa

(p3c <- pp_check(mt1, ndraws = 100)+ 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))
(p3d <- pp_check(mt1, type = "stat_grouped", 
                 stat = "mean", group = "Taxon")+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)))

# Interactive models freshwater

(p4c <- pp_check(mmu_fr, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

# Interactive models marine

(p5c <- pp_check(mmu_mar, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

# Interactive models terrestrial

(p6c <- pp_check(mmu_ter, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

# Interactive models amphibians

(p7c <- pp_check(mmu_am, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

# Interactive models birds

(p8c <- pp_check(mmu_bi, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

# Interactive models fishes

(p9c <- pp_check(mmu_f, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

# Interactive models mammals

(p10c <- pp_check(mmu_ma, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

# Interactive models reptiles

(p11c <- pp_check(mmu_rep, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

# Combine  plots into a single figure

(combined_figure <-  (p1c+p2c+p3c+p4c+
                p5c+p6c+p7c+p8c+
                p9c+p10c+p11c & 
                labs(x = NULL, y = NULL)) +
  plot_annotation(tag_levels = "a")+
  plot_layout(nrow = 3, guides = "collect") &
    theme(plot.tag = element_text(face = 'bold')))

# Create a separate plot for the y-axis label 

ylabel <- ggplot(data.frame(x = 1, y = 1)) +
  geom_text(aes(x, y, label = "Density"), size=5, angle = 90) + 
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the original figure with the y-axis label and x-axis label

(figureS5 <- cowplot::plot_grid(ylabel, combined_figure,
                                rel_widths = c(1, 25)))

# Save the combined figure 

ggsave("Figure S5.pdf", figureS5, 
       path = ResultPath,
       height = 8, width = 10)

# Multicollinearity ------------------------------------------------------------

check_collinearity(m1)
check_collinearity(ms1)
check_collinearity(mt1)
check_collinearity(mmu_fr)
check_collinearity(mmu_mar)
check_collinearity(mmu_ter)
check_collinearity(mmu_am)
check_collinearity(mmu_bi)
check_collinearity(mmu_f)
check_collinearity(mmu_ma)
check_collinearity(mmu_rep)

