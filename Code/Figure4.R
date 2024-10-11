# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Figure4.R         
# - DATE:        07/01/2021
# - DESCRIPTION: Figure 2. 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com), Duncan O'Brien (duncan.a.obrien@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(brms)
library(tidyverse)
library(tidybayes)
library(parallel)
library(doSNOW)
library(patchwork)
library(rstan)

source("Code/utils/threat_counterfac_draws.R")
source("Code/utils/threat_counterfac_pred.R")

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
                  plot.margin = unit(c(0, 0, 0, 0), "cm")))

# Load the model 

m1 <- read_rds("Results/models/mod_global_rerun.RDS")

# Conterfactual tests ----------------------------------------------------------

# create a vector with all the threats in the dataset

threat_cols <-  colnames(m1$data)[grepl(paste(c("pollution","habitatl",
                                                "climatechange","invasive", 
                                                "exploitation","disease"),
                                              collapse = "|"),
                                        colnames(m1$data))]


# Predict the population trends in the "no intervention" scenario

pop_perd <- brms::posterior_epred(m1,
                                  newdata = m1$data %>% 
                                    filter_at(threat_cols, any_vars(. != "0")),
                                  re.form = NULL,
                                  incl_autocor = FALSE,
                                  sort = TRUE, 
                                  ndraws = 1000) %>% 
  as.data.frame() %>%
  mutate(.draw = 1:NROW(.)) %>%
  #extract posterior draws for the data used to create the model
  pivot_longer(-.draw,names_to = "index",values_to = ".value") %>%
  cbind(m1$data %>%
          filter_at(threat_cols, any_vars(. != "0")) %>% 
          dplyr::select(series, time), row.names = NULL) %>% 
  reframe(.value = mean(diff(.value)/diff(time)),.by = c(series,.draw)) %>% 
  mutate(counterfac="none") %>% 
  reframe(mn = mean(.value), .by=c(series, counterfac))

# Create the counterfactual for all the threats

# We use the function counterfactual draws to estimate the different population 
# trends under the different scenarios removing the threats

counter_all <- threat_counterfac_draws(m1,
                                       threat_comns = paste(c("pollution","habitatl",
                                                              "climatechange","invasive", 
                                                              "exploitation","disease"),
                                                            collapse = "."),
                                       re.form = NULL,
                                       ndraws = 1000,
                                       center_dydx = "mean",
                                       n.cores = 4, trend=T) %>% 
  mutate(counterfac="All")


# Create the different counterfactual scenarios

# We use the function counterfactual draws to estimate the different population 
# trends under the different scenarios removing one threat

counter_fac_data <- threat_counterfac_draws(m1,
                                            threat_comns = threat_cols,
                                            re.form = NULL,
                                            ndraws = 1000,
                                            center_dydx = "mean",
                                            n.cores = 4, trend=T) %>%
  # Join with the none counterfactual scenario that we just created
  rbind(pop_perd, counter_all) %>%
  # Made "none" as the first level of the counterfactual 
  mutate(counterfac = fct_relevel(counterfac, "none"))

# We summarise it 

scenarios_mean <- counter_fac_data %>% 
  group_by(counterfac) %>% 
  summarise(m=median(mn)) %>% 
  arrange(desc(m))

# Plots ------------------------------------------------------------------------

# Create a palette

threat_palette<-c(MetBrewer::met.brewer(name="Hokusai1", n=6, type="continuous"))

palette <- data.frame(counterfac = c("No pollution","No habitat loss", 
                                     "No climate change", "No invasive", 
                                     "No exploitation", "No disease", "All"),
                      fill_col = as.factor(c(threat_palette, "grey50")))

## Panel a: one threat scenarios --------------------------------------------------------

(g4 <- counter_fac_data %>% 
   # add a variable counting the number of threats
   mutate(number=str_count(counterfac, '\\.')+1) %>%
   group_by(counterfac) %>% 
   mutate(pos=median(mn),
          counterfac = gsub("habitatl", "habitat loss", counterfac),
          counterfac = gsub("climatechange", "climate change", counterfac)) %>% 
   ungroup() %>%   
   # Join the palette
   left_join(palette, by = "counterfac") %>%
   # Further modify the levels
   mutate(counterfac = gsub("habitat loss", "habitat\nloss", counterfac),
          counterfac = gsub("climate change", "climate\nchange", counterfac)) %>% 
   # Keep only the single threats
   filter(number==1 | number==6, 
          counterfac!="none") %>% 
   # Start the plot
   ggplot(aes(y = mn,
              x=reorder(counterfac, pos),
              fill=fill_col, 
              colour=fill_col)) +
   # geom_violin(alpha=0.5, colour="white")+
   # geom_boxplot(aes(colour=counterfac),
   #              fill="white", width = .2,
   #              outlier.shape = NA)+
   stat_halfeye(adjust = .5,
                width = .6,
                .width = 0,
                alpha=.7,
                justification = -.3,
                point_colour = NA,
                normalize = "groups") +
   geom_boxplot(width = .2,
                outlier.shape = NA,
                alpha=.7,
                colour="black") +
   gghalves::geom_half_point(show.legend = FALSE,
                             side = "r",
                             range_scale = .2,
                             shape=21, alpha=0.2,
                             colour="grey25") +
   scale_fill_manual(values = levels(palette$fill_col),
                     guide = "none")+
   scale_colour_manual(values = levels(palette$fill_col),
                       guide = "none")+
   geom_hline(yintercept = 0,
              linetype = "solid", 
              colour="grey50") +
   geom_hline(yintercept = subset(scenarios_mean,counterfac == "none")$m,
              linetype = "dashed", 
              colour="grey50") +
   ylab(expression(paste("Mean population trend (",Delta,"y","/",Delta,"x)"))) + 
   xlab("Counterfactual scenarios") +
   ylim(c(-0.2,0.2)))

# Save it

ggsave("Results/figures/Figure4.pdf", g4, 
       width = 10, height = 6)
