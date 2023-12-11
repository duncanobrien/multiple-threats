require(tidyverse)
require(brms)
require(tidybayes)

mod1 <- readRDS("Results/models/non-linear.RDS")
mod2 <- readRDS("Results/models/simple-slopes.RDS")
mod3 <- readRDS("Results/models/state-space.RDS")

mod1_draws <- mod1 |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  dplyr::mutate(.variable = gsub("b_trend_scaled_time:threats","",.variable),
                .variable =  gsub("InvasivesppDgenes","Invasive species",.variable),
                .variable = gsub("Climatechange","Climate change",.variable),
                .variable =  gsub("Habitatloss","Habitat loss",.variable),
                .variable =  gsub("none","None",.variable)) |>
  dplyr::mutate(.variable = factor(.variable),
                .variable = fct_relevel(.variable,"None")) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

mod2_draws <- emmeans::emtrends(mod2, ~ threats,var = "scaled_time") |>
  tidybayes::gather_emmeans_draws() |>
  dplyr::rename(.variable = threats) |>
  dplyr::mutate(.variable =  gsub("Invasive spp/genes","Invasive species",.variable),
                .variable =  gsub("none","None",.variable)) |>
  dplyr::mutate(.variable = factor(.variable),
                .variable = fct_relevel(.variable,"None")) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

mod3_draws <- mod3 |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  dplyr::mutate(.variable = gsub("b_threats","",.variable),
                .variable =  gsub("InvasivesppDgenes","Invasive species",.variable),
                .variable = gsub("Climatechange","Climate change",.variable),
                .variable =  gsub("Habitatloss","Habitat loss",.variable),
                .variable =  gsub("none","None",.variable)) |>
  dplyr::mutate(.variable = factor(.variable),
                .variable = fct_relevel(.variable,"None")) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

threats_plot_data <- rbind(mod1_draws |>
                             dplyr::mutate(model = "Non-linear"),
                           mod2_draws  |>
                             dplyr::mutate(model = "Simple slopes"),
                           mod3_draws  |>
                             dplyr::mutate(model = "State space"))

ggsave("Results/models/model_comparison.pdf",
ggplot(threats_plot_data,aes(y = .variable, x = .value,col = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper),interval_size_range = c(0.8, 2)) +
  labs(x="Coefficient", y = "Threat") + 
  scale_y_discrete(limits = rev(levels(threats_plot_data$.variable))) + 
  #scale_color_manual(values = c("black","grey80"),guide = "none") + 
  scale_color_manual(values = c("#F5A433","#8B33F5","#57A06B"), guide = "none") +
  facet_wrap(~model,scales = "free_x") + 
  theme_minimal()+
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
        plot.title = element_text(hjust = 0.5)),
width = 10,height = 10)
