require(tidyverse)
require(brms)
require(tidybayes)
require(patchwork)

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

##########################################################################################
## visualise interacting threats ##
##########################################################################################

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

##########################################################################################
## pp_checks ##
##########################################################################################
ggsave("Results/models/model_ppc_comparison.pdf",
       brms::pp_check(mod1) + ggtitle("Non-linear") +
         brms::pp_check(mod2) + ggtitle("Simple slopes") +
         brms::pp_check(mod3) + coord_cartesian(xlim=c(-2.5,2.5)) + ggtitle("State space"),
       width = 10, height = 6)


##########################################################################################
## realm models ##
##########################################################################################
mod2_mar <- readRDS("Results/models/simple-slopes-marine.RDS")
mod2_fre <- readRDS("Results/models/simple-slopes-fresh.RDS")
mod2_ter <- readRDS("Results/models/simple-slopes-terra.RDS")

mod2_mar_draws <- emmeans::emtrends(mod2_mar, ~ none,var = "scaled_time",at = list(none = 1)) |>
  tidybayes::gather_emmeans_draws() |>
  dplyr::rename(".variable" = none) |>
  dplyr::mutate(.variable = "None") |>
  rbind(emmeans::emtrends(mod2_mar, ~ pollution,var = "scaled_time",at = list(pollution = 1)) |>
          tidybayes::gather_emmeans_draws() |>
          dplyr::rename(".variable" = pollution) |>
          dplyr::mutate(.variable = "Pollution"),
        emmeans::emtrends(mod2_mar, ~ habitatl,var = "scaled_time",at = list(habitatl = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = habitatl) |>
          dplyr::mutate(.variable = "Habitat loss"),
        emmeans::emtrends(mod2_mar, ~ climatechange,var = "scaled_time",at = list(climatechange = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = climatechange) |>
          dplyr::mutate(.variable = "Climate change"),
        emmeans::emtrends(mod2_mar, ~ invasive,var = "scaled_time",at = list(invasive = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = invasive) |>
          dplyr::mutate(.variable = "Invasive species"),
        emmeans::emtrends(mod2_mar, ~ exploitation,var = "scaled_time",at = list(exploitation = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = exploitation) |>
          dplyr::mutate(.variable = "Exploitation"),
        emmeans::emtrends(mod2_mar, ~ disease,var = "scaled_time",at = list(disease = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = disease) |>
          dplyr::mutate(.variable = "Disease")) |>
  dplyr::mutate(.variable = factor(.variable),
                .variable = fct_relevel(.variable,"None")) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

mod2_fre_draws <- emmeans::emtrends(mod2_fre, ~ none,var = "scaled_time",at = list(none = 1)) |>
  tidybayes::gather_emmeans_draws() |>
  dplyr::rename(".variable" = none) |>
  dplyr::mutate(.variable = "None") |>
  rbind(emmeans::emtrends(mod2_fre, ~ pollution,var = "scaled_time",at = list(pollution = 1)) |>
          tidybayes::gather_emmeans_draws() |>
          dplyr::rename(".variable" = pollution) |>
          dplyr::mutate(.variable = "Pollution"),
        emmeans::emtrends(mod2_fre, ~ habitatl,var = "scaled_time",at = list(habitatl = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = habitatl) |>
          dplyr::mutate(.variable = "Habitat loss"),
        emmeans::emtrends(mod2_fre, ~ climatechange,var = "scaled_time",at = list(climatechange = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = climatechange) |>
          dplyr::mutate(.variable = "Climate change"),
        emmeans::emtrends(mod2_fre, ~ invasive,var = "scaled_time",at = list(invasive = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = invasive) |>
          dplyr::mutate(.variable = "Invasive species"),
        emmeans::emtrends(mod2_fre, ~ exploitation,var = "scaled_time",at = list(exploitation = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = exploitation) |>
          dplyr::mutate(.variable = "Exploitation"),
        emmeans::emtrends(mod2_fre, ~ disease,var = "scaled_time",at = list(disease = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = disease) |>
          dplyr::mutate(.variable = "Disease")) |>
  dplyr::mutate(.variable = factor(.variable),
                .variable = fct_relevel(.variable,"None")) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

mod2_ter_draws <- emmeans::emtrends(mod2_ter, ~ none,var = "scaled_time",at = list(none = 1)) |>
  tidybayes::gather_emmeans_draws() |>
  dplyr::rename(".variable" = none) |>
  dplyr::mutate(.variable = "None") |>
  rbind(emmeans::emtrends(mod2_ter, ~ pollution,var = "scaled_time",at = list(pollution = 1)) |>
          tidybayes::gather_emmeans_draws() |>
          dplyr::rename(".variable" = pollution) |>
          dplyr::mutate(.variable = "Pollution"),
        emmeans::emtrends(mod2_ter, ~ habitatl,var = "scaled_time",at = list(habitatl = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = habitatl) |>
          dplyr::mutate(.variable = "Habitat loss"),
        emmeans::emtrends(mod2_ter, ~ climatechange,var = "scaled_time",at = list(climatechange = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = climatechange) |>
          dplyr::mutate(.variable = "Climate change"),
        emmeans::emtrends(mod2_ter, ~ invasive,var = "scaled_time",at = list(invasive = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = invasive) |>
          dplyr::mutate(.variable = "Invasive species"),
        emmeans::emtrends(mod2_ter, ~ exploitation,var = "scaled_time",at = list(exploitation = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = exploitation) |>
          dplyr::mutate(.variable = "Exploitation"),
        emmeans::emtrends(mod2_ter, ~ disease,var = "scaled_time",at = list(disease = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = disease) |>
          dplyr::mutate(.variable = "Disease")) |>
  dplyr::mutate(.variable = factor(.variable),
                .variable = fct_relevel(.variable,"None")) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

realms_plot_data <- rbind(mod2_fre_draws |>
                             dplyr::mutate(System = "Freshwater"),
                          mod2_mar_draws  |>
                             dplyr::mutate(System = "Marine"),
                          mod2_ter_draws  |>
                             dplyr::mutate(System = "Terrestrial"))

ggplot(realms_plot_data,aes(y = .value , x = .variable,col = System)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper),interval_size_range = c(0.8, 2),
                                position = position_dodge(0.7)) +
  labs(x="Coefficient", y = "Threat") + 
  scale_color_manual(values = c("#A1D6E2","#336B87", "#CC954E")) +
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
        plot.title = element_text(hjust = 0.5))

##########################################################################################
## taxon models ##
##########################################################################################
mod2_am <- readRDS("Results/models/simple-slopes-amph.RDS")
mod2_bir <- readRDS("Results/models/simple-slopes-amph.RDS")
mod2_mam <- readRDS("Results/models/simple-slopes-mammals.RDS")
mod2_rep <- readRDS("Results/models/simple-slopes-rept.RDS")

mod2_am_draws <- emmeans::emtrends(mod2_am, ~ none,var = "scaled_time",at = list(none = 1)) |>
  tidybayes::gather_emmeans_draws() |>
  dplyr::rename(".variable" = none) |>
  dplyr::mutate(.variable = "None") |>
  rbind(emmeans::emtrends(mod2_am, ~ pollution,var = "scaled_time",at = list(pollution = 1)) |>
          tidybayes::gather_emmeans_draws() |>
          dplyr::rename(".variable" = pollution) |>
          dplyr::mutate(.variable = "Pollution"),
        emmeans::emtrends(mod2_am, ~ habitatl,var = "scaled_time",at = list(habitatl = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = habitatl) |>
          dplyr::mutate(.variable = "Habitat loss"),
        emmeans::emtrends(mod2_am, ~ climatechange,var = "scaled_time",at = list(climatechange = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = climatechange) |>
          dplyr::mutate(.variable = "Climate change"),
        emmeans::emtrends(mod2_am, ~ invasive,var = "scaled_time",at = list(invasive = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = invasive) |>
          dplyr::mutate(.variable = "Invasive species"),
        emmeans::emtrends(mod2_am, ~ exploitation,var = "scaled_time",at = list(exploitation = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = exploitation) |>
          dplyr::mutate(.variable = "Exploitation"),
        emmeans::emtrends(mod2_am, ~ disease,var = "scaled_time",at = list(disease = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = disease) |>
          dplyr::mutate(.variable = "Disease")) |>
  dplyr::mutate(.variable = factor(.variable),
                .variable = fct_relevel(.variable,"None")) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

mod2_bir_draws <- emmeans::emtrends(mod2_bir, ~ none,var = "scaled_time",at = list(none = 1)) |>
  tidybayes::gather_emmeans_draws() |>
  dplyr::rename(".variable" = none) |>
  dplyr::mutate(.variable = "None") |>
  rbind(emmeans::emtrends(mod2_bir, ~ pollution,var = "scaled_time",at = list(pollution = 1)) |>
          tidybayes::gather_emmeans_draws() |>
          dplyr::rename(".variable" = pollution) |>
          dplyr::mutate(.variable = "Pollution"),
        emmeans::emtrends(mod2_bir, ~ habitatl,var = "scaled_time",at = list(habitatl = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = habitatl) |>
          dplyr::mutate(.variable = "Habitat loss"),
        emmeans::emtrends(mod2_bir, ~ climatechange,var = "scaled_time",at = list(climatechange = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = climatechange) |>
          dplyr::mutate(.variable = "Climate change"),
        emmeans::emtrends(mod2_bir, ~ invasive,var = "scaled_time",at = list(invasive = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = invasive) |>
          dplyr::mutate(.variable = "Invasive species"),
        emmeans::emtrends(mod2_bir, ~ exploitation,var = "scaled_time",at = list(exploitation = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = exploitation) |>
          dplyr::mutate(.variable = "Exploitation"),
        emmeans::emtrends(mod2_bir, ~ disease,var = "scaled_time",at = list(disease = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = disease) |>
          dplyr::mutate(.variable = "Disease")) |>
  dplyr::mutate(.variable = factor(.variable),
                .variable = fct_relevel(.variable,"None")) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

mod2_mam_draws <- emmeans::emtrends(mod2_mam, ~ none,var = "scaled_time",at = list(none = 1)) |>
  tidybayes::gather_emmeans_draws() |>
  dplyr::rename(".variable" = none) |>
  dplyr::mutate(.variable = "None") |>
  rbind(emmeans::emtrends(mod2_mam, ~ pollution,var = "scaled_time",at = list(pollution = 1)) |>
          tidybayes::gather_emmeans_draws() |>
          dplyr::rename(".variable" = pollution) |>
          dplyr::mutate(.variable = "Pollution"),
        emmeans::emtrends(mod2_mam, ~ habitatl,var = "scaled_time",at = list(habitatl = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = habitatl) |>
          dplyr::mutate(.variable = "Habitat loss"),
        emmeans::emtrends(mod2_mam, ~ climatechange,var = "scaled_time",at = list(climatechange = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = climatechange) |>
          dplyr::mutate(.variable = "Climate change"),
        emmeans::emtrends(mod2_mam, ~ invasive,var = "scaled_time",at = list(invasive = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = invasive) |>
          dplyr::mutate(.variable = "Invasive species"),
        emmeans::emtrends(mod2_mam, ~ exploitation,var = "scaled_time",at = list(exploitation = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = exploitation) |>
          dplyr::mutate(.variable = "Exploitation"),
        emmeans::emtrends(mod2_mam, ~ disease,var = "scaled_time",at = list(disease = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = disease) |>
          dplyr::mutate(.variable = "Disease")) |>
  dplyr::mutate(.variable = factor(.variable),
                .variable = fct_relevel(.variable,"None")) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

mod2_rep_draws <- emmeans::emtrends(mod2_rep, ~ none,var = "scaled_time",at = list(none = 1)) |>
  tidybayes::gather_emmeans_draws() |>
  dplyr::rename(".variable" = none) |>
  dplyr::mutate(.variable = "None") |>
  rbind(emmeans::emtrends(mod2_rep, ~ pollution,var = "scaled_time",at = list(pollution = 1)) |>
          tidybayes::gather_emmeans_draws() |>
          dplyr::rename(".variable" = pollution) |>
          dplyr::mutate(.variable = "Pollution"),
        emmeans::emtrends(mod2_rep, ~ habitatl,var = "scaled_time",at = list(habitatl = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = habitatl) |>
          dplyr::mutate(.variable = "Habitat loss"),
        emmeans::emtrends(mod2_rep, ~ climatechange,var = "scaled_time",at = list(climatechange = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = climatechange) |>
          dplyr::mutate(.variable = "Climate change"),
        emmeans::emtrends(mod2_rep, ~ invasive,var = "scaled_time",at = list(invasive = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = invasive) |>
          dplyr::mutate(.variable = "Invasive species"),
        emmeans::emtrends(mod2_rep, ~ exploitation,var = "scaled_time",at = list(exploitation = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = exploitation) |>
          dplyr::mutate(.variable = "Exploitation"),
        emmeans::emtrends(mod2_rep, ~ disease,var = "scaled_time",at = list(disease = 1)) |>
          tidybayes::gather_emmeans_draws()|>
          dplyr::rename(".variable" = disease) |>
          dplyr::mutate(.variable = "Disease")) |>
  dplyr::mutate(.variable = factor(.variable),
                .variable = fct_relevel(.variable,"None")) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

taxon_plot_data <- rbind(mod2_am_draws |>
                            dplyr::mutate(Taxon = "Amphibians"),
                         mod2_bir_draws |>
                           dplyr::mutate(Taxon = "Birds"),
                         mod2_mam_draws  |>
                            dplyr::mutate(Taxon = "Mammals"),
                         mod2_rep_draws  |>
                            dplyr::mutate(Taxon = "Reptiles"))

ggplot(taxon_plot_data,aes(y = .value , x = .variable,col = Taxon)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper),interval_size_range = c(0.8, 2),
                                position = position_dodge(0.7)) +
  labs(x="Coefficient", y = "Threat") + 
  scale_color_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5)) +
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
        plot.title = element_text(hjust = 0.5))


