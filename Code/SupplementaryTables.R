# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SupplementaryTables.R
# - DATE:        07/01/2021
# - DESCRIPTION: Basic diagnostics for the models.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list = ls(all = TRUE)) #remove everything

library(tidyverse)
library(bayestestR)
library(rstan)

source("Code/utils/prep_data_grid_fn.R")
source("Code/utils/threat_post_draws.R")

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

sys_model_ls <- list(m_fr, m_mar, m_ter) %>%
  `names<-`(c("Freshwater", "Marine", "Terrestrial"))

taxon_model_ls <- list(m_am, m_bi, m_f, m_ma, m_rep) %>%
  `names<-`(c("Amphibian", "Bird", "Fish", "Mammal", "Reptile"))

all_threats <- c("pollution",
                 "habitatl",
                 "climatechange",
                 "invasive",
                 "exploitation",
                 "disease")

# Table S1: Global model coefficients ##########

tableS1 <- describe_posterior(m1) %>%
  select(-c(pd, CI, ROPE_CI, ROPE_low, ROPE_high, ROPE_Percentage)) %>%
  mutate(Parameter = gsub("b_", "", Parameter)) %>%
  mutate(across(Median:CI_high, ~ round(., digits = 3))) %>%
  mutate(across(Rhat:ESS, ~ round(., digits = 2)))

write.csv(tableS1, "Results/tables/TableS1.csv", row.names = FALSE)

# Table S2: System model coefficients ##########

tableS2 <- lapply(seq_along(sys_model_ls), function(mod) {
  describe_posterior(sys_model_ls[[mod]]) %>%
    select(-c(pd, CI, ROPE_CI, ROPE_low, ROPE_high, ROPE_Percentage)) %>%
    mutate(Parameter = gsub("b_", "", Parameter)) %>%
    mutate(across(Median:CI_high, ~ round(., digits = 3))) %>%
    mutate(across(Rhat:ESS, ~ round(., digits = 2))) %>%
    mutate(system = names(sys_model_ls)[[mod]])
}) %>%
  bind_rows()

write.csv(tableS2, "Results/tables/TableS2.csv", row.names = FALSE)

# Table S3: Taxon model coefficients ##########

tableS3 <- lapply(seq_along(taxon_model_ls), function(mod) {
  describe_posterior(taxon_model_ls[[mod]]) %>%
    select(-c(pd, CI, ROPE_CI, ROPE_low, ROPE_high, ROPE_Percentage)) %>%
    mutate(Parameter = gsub("b_", "", Parameter)) %>%
    mutate(across(Median:CI_high, ~ round(., digits = 3))) %>%
    mutate(across(Rhat:ESS, ~ round(., digits = 2))) %>%
    mutate(system = names(taxon_model_ls)[[mod]])
}) %>%
  bind_rows()

write.csv(tableS3, "Results/tables/TableS3.csv", row.names = FALSE)

# Table S4: Threat interaction types ##########

all_threats <- c("pollution",
                 "habitatl",
                 "climatechange",
                 "invasive",
                 "exploitation",
                 "disease")

# Extract column names from the model

threatcols <- colnames(m1$data)[grepl(paste(all_threats, collapse = "|"), colnames(m1$data))]

# Generate additive combinations of threat columns (i.e., "threat1 + threat2")

additive_cols <- do.call("c", lapply(strsplit(threatcols, "[.]"), function(x) {
  paste(x, collapse = " + ")
}))

# Combine original and additive columns into a unique set

target_cols <- unique(c(threatcols, additive_cols))

# Estimate posterior time series for each threat combination (singular, interactive, additive)

postdraws_intvadd <- threat_post_draws(
  model = m1,
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
    interaction.type = factor(case_when(
      .upper < 0 ~ "synergistic",
      .lower > 0 ~ "antagonistic",
      TRUE ~ "additive"
    ),
    levels = c("synergistic", "antagonistic", "additive"))
  )

# Calculate the proportions

tableS4 <- post_interval_intvadd_diff %>%
  separate_rows(threats, sep = "\\.") %>%
  mutate(interaction.type = factor(
    interaction.type,
    levels = c("synergistic", "antagonistic", "additive"),
    labels = c("Synergistic", "Antagonistic", "Additive")
  )) %>%
  group_by(threats, interaction.type) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  complete(threats, interaction.type) %>%
  group_by(threats) %>%
  mutate(n = replace_na(n, 0), freq = (n / sum(n))) %>%
  ungroup() %>%
  mutate(threats = str_to_title(threats)) %>%
  mutate(
    threats = case_when(
      threats %in% "Climatechange" ~ "Climate change",
      threats %in% "Habitatl" ~ "Habitat loss",
      TRUE ~ threats
    )
  ) %>%
  setNames(c("Threat", "Interaction type", "n", "Frequency"))

write.csv(tableS4, "Results/tables/TableS4.csv", row.names = FALSE)

# Table S5: Threat interaction types by system and taxa ##########

tableS5 <- lapply(seq_along(sys_model_ls), function(mod) {
  # Extract column names from the model
  
  threatcols <- colnames(sys_model_ls[[mod]]$data)[grepl(paste(all_threats, collapse = "|"),
                                                         colnames(sys_model_ls[[mod]]$data))]
  
  # Generate additive combinations of threat columns (i.e., "threat1 + threat2")
  
  additive_cols <- do.call("c", lapply(strsplit(threatcols, "[.]"), function(x) {
    paste(x, collapse = " + ")
  }))
  
  # Combine original and additive columns into a unique set
  
  target_cols <- unique(c(threatcols, additive_cols))
  
  # Estimate posterior time series for each threat combination (singular, interactive, additive)
  
  postdraws_intvadd <- threat_post_draws(
    model = sys_model_ls[[mod]],
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
    ggdist::median_qi(.width = c(.95, .8, .5),
                      .exclude = c(".draw"))
  
  # Compare the difference in derivatives between shared additive and interactive threats
  
  post_intvadd_diff <- do.call("rbind", lapply(additive_cols[grepl("\\+", additive_cols)], function(x) {
    ss <- post_dydx_intvadd %>%
      subset(threat %in% c(x, gsub(" \\+ ", ".", x))) %>%
      reframe(.value = diff(c(.value[2], .value[1])),
              .by = c(threat_group, .draw)) %>%
      mutate(threats = gsub(" \\+ ", ".", x))
  }))
  
  # Finalize the comparison, drop missing values, and classify interaction types
  
  post_interval_intvadd_diff <- post_intvadd_diff %>%
    na.omit() %>%
    group_by(threat_group, threats) %>%
    ggdist::median_qi(.width = c(.95, .8, .5),
                      .exclude = c(".draw")) %>%
    mutate(
      interaction.type = factor(case_when(
        .upper < 0 ~ "synergistic",
        .lower > 0 ~ "antagonistic",
        TRUE ~ "additive"
      ),
      levels = c("synergistic", "antagonistic", "additive"))
    )
  
  # Calculate the proportions
  
  data_ad <- post_interval_intvadd_diff %>%
    separate_rows(threats, sep = "\\.") %>%
    mutate(interaction.type = factor(
      interaction.type,
      levels = c("synergistic", "antagonistic", "additive"),
      labels = c("Synergistic", "Antagonistic", "Additive")
    )) %>%
    group_by(threats, interaction.type) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    complete(threats, interaction.type) %>%
    group_by(threats) %>%
    mutate(n = replace_na(n, 0), freq = (n / sum(n))) %>%
    ungroup() %>%
    mutate(threats = str_to_title(threats)) %>%
    mutate(
      threats = case_when(
        threats %in% "Climatechange" ~ "Climate change",
        threats %in% "Habitatl" ~ "Habitat loss",
        TRUE ~ threats
      )
    ) %>%
    setNames(c("Threat", "Interaction type", "n", "Frequency")) %>%
    add_column(System = names(sys_model_ls)[[mod]], .after = "Threat")
  
  return(data_ad)
  
}) %>%
  bind_rows() %>%
  arrange(Threat, System)

write.csv(tableS5, "Results/tables/TableS5.csv", row.names = FALSE)

# Table S6: Threat interaction types by system and taxa ##########

tableS6 <- lapply(seq_along(taxon_model_ls), function(mod) {
  # Extract column names from the model
  
  threatcols <- colnames(taxon_model_ls[[mod]]$data)[grepl(paste(all_threats, collapse = "|"),
                                                         colnames(taxon_model_ls[[mod]]$data))]
  
  # Generate additive combinations of threat columns (i.e., "threat1 + threat2")
  
  additive_cols <- do.call("c", lapply(strsplit(threatcols, "[.]"), function(x) {
    paste(x, collapse = " + ")
  }))
  
  # Combine original and additive columns into a unique set
  
  target_cols <- unique(c(threatcols, additive_cols))
  
  # Estimate posterior time series for each threat combination (singular, interactive, additive)
  
  postdraws_intvadd <- threat_post_draws(
    model = taxon_model_ls[[mod]],
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
    ggdist::median_qi(.width = c(.95, .8, .5),
                      .exclude = c(".draw"))
  
  # Compare the difference in derivatives between shared additive and interactive threats
  
  post_intvadd_diff <- do.call("rbind", lapply(additive_cols[grepl("\\+", additive_cols)], function(x) {
    ss <- post_dydx_intvadd %>%
      subset(threat %in% c(x, gsub(" \\+ ", ".", x))) %>%
      reframe(.value = diff(c(.value[2], .value[1])),
              .by = c(threat_group, .draw)) %>%
      mutate(threats = gsub(" \\+ ", ".", x))
  }))
  
  # Finalize the comparison, drop missing values, and classify interaction types
  
  post_interval_intvadd_diff <- post_intvadd_diff %>%
    na.omit() %>%
    group_by(threat_group, threats) %>%
    ggdist::median_qi(.width = c(.95, .8, .5),
                      .exclude = c(".draw")) %>%
    mutate(
      interaction.type = factor(case_when(
        .upper < 0 ~ "synergistic",
        .lower > 0 ~ "antagonistic",
        TRUE ~ "additive"
      ),
      levels = c("synergistic", "antagonistic", "additive"))
    )
  
  # Calculate the proportions
  
  data_ad <- post_interval_intvadd_diff %>%
    separate_rows(threats, sep = "\\.") %>%
    mutate(interaction.type = factor(
      interaction.type,
      levels = c("synergistic", "antagonistic", "additive"),
      labels = c("Synergistic", "Antagonistic", "Additive")
    )) %>%
    group_by(threats, interaction.type) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    complete(threats, interaction.type) %>%
    group_by(threats) %>%
    mutate(n = replace_na(n, 0), freq = (n / sum(n))) %>%
    ungroup() %>%
    mutate(threats = str_to_title(threats)) %>%
    mutate(
      threats = case_when(
        threats %in% "Climatechange" ~ "Climate change",
        threats %in% "Habitatl" ~ "Habitat loss",
        TRUE ~ threats
      )
    ) %>%
    setNames(c("Threat", "Interaction type", "n", "Frequency")) %>%
    add_column(Taxon = names(taxon_model_ls)[[mod]], .after = "Threat")
  
  return(data_ad)
}) %>%
  bind_rows() %>%
  arrange(Threat, Taxon)

write.csv(tableS6, "Results/tables/TableS6.csv", row.names = FALSE)
