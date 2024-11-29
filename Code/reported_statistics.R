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
source("Code/utils/prep_data_grid_fn.R")
source("Code/utils/threat_post_draws.R")
source("Code/utils/threat_counterfac_draws.R")
source("Code/utils/threat_counterfac_pred.R")

#set.seed(43)
set.seed(2024)

# Load the model

mod <- readRDS("Results/models/mod_global_rerun.RDS")

# Load the modelled data

load("Data/data_models.RData")

# The threats to project

threats <- c("pollution",
             "habitatl",
             "climatechange",
             "invasive",
             "exploitation",
             "disease")

threats_col <-  colnames(mod$data)[grepl(paste(threats, collapse = "|"), colnames(mod$data))]

# Number of time series

ts_coverage <- mod_dat_full |> 
  dplyr::mutate(threatened = rowSums(dplyr::across(threats_col) =="1")) |>
  dplyr::mutate(threatened = ifelse(threatened == 0,"0","1")) |>
  dplyr::summarise(n_ts = length(unique(series)),.by = threatened)

# Number of species

spp_coverage <- mod_dat_full |> 
  dplyr::summarise(n_spp = length(unique(series)),.by = Taxon)

# Number of species

sys_coverage <-mod_dat_full |> 
  dplyr::summarise(n_sys = length(unique(series)), .by = System)

## Population trends projections convertion to %--------------------------------------

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
})) |>
  mutate(.value_percent = (exp(.value)-1)*100)

dydx_interval_percent_a <- post_proj  %>%
  select(-.value) %>%
  group_by(threat_group, int_group) %>%
  ggdist::median_qi(
    .width = c(.95, .8, .5),
    .exclude = c(".draw", "threat", "threat_group")
  ) %>%
  mutate(
    int_group = ifelse(int_group == "single", "Singular", "Interactive"),
    int = ordered(int_group, c("Singular", "Interactive"))
  ) |>
  subset(.width == 0.95)

dydx_interval_percent_b <- post_proj  %>%
  select(-.value) %>%
  group_by(threat,threat_group, int_group) %>%
  ggdist::median_qi(
    .width = c(.95, .8, .5),
    .exclude = c(".draw", "threat", "threat_group")
  ) %>%
  mutate(
    int_group = ifelse(int_group == "single", "Singular", "Interactive"),
    int = ordered(int_group, c("Singular", "Interactive"))
  ) |>
  subset(.width == 0.95 & int_group == "Interactive") |>
  arrange(.value_percent)

## Counteractual trends projections convertion to %--------------------------------------

# Predict the population trends in the "no intervention" scenario

pop_perd <- brms::posterior_epred(mod,
                                  newdata = mod$data %>% 
                                    filter_at(threats_col, any_vars(. != "0")),
                                  re.form = NULL,
                                  incl_autocor = FALSE,
                                  sort = TRUE, 
                                  ndraws = 1000) %>% 
  as.data.frame() %>%
  mutate(.draw = 1:NROW(.)) %>%
  #extract posterior draws for the data used to create the model
  pivot_longer(-.draw,names_to = "index",values_to = ".value") %>%
  cbind(mod$data %>%
          filter_at(threats_col, any_vars(. != "0")) %>% 
          dplyr::select(series, time), row.names = NULL) %>% 
  reframe(.value = mean(diff(.value)/diff(time)),.by = c(series,.draw)) %>% 
  mutate(counterfac="none") %>% 
  reframe(mn = mean(.value), .by=c(series, counterfac))

# Create the counterfactual for all the threats

# We use the function counterfactual draws to estimate the different population 
# trends under the different scenarios removing the threats

counter_all <- threat_counterfac_draws(mod,
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

counter_fac_data <- threat_counterfac_draws(mod,
                                            threat_comns = threats_col,
                                            re.form = NULL,
                                            ndraws = 1000,
                                            center_dydx = "mean",
                                            n.cores = 4, trend=T) %>%
  # Join with the none counterfactual scenario that we just created
  rbind(pop_perd, counter_all) %>%
  # Made "none" as the first level of the counterfactual 
  mutate(counterfac = fct_relevel(counterfac, "none"))

# We summarise it 

scenarios_mean_diff <- counter_fac_data %>% 
  group_by(counterfac) %>% 
  summarise(m=median(mn)) %>% 
  arrange(desc(m)) %>%
  mutate(ref_none = m[counterfac=="none"],
         ref_all = m[counterfac=="All"]) |>
  mutate(number=str_count(counterfac, '\\.')+1) |>
  mutate(percent_diff = abs(((m-ref_none)/ref_none)*100), .by = counterfac) |>
  subset(number == 1) |>
  arrange(desc(percent_diff))
