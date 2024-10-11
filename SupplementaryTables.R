# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SupplementaryTables.R
# - DATE:        07/01/2021
# - DESCRIPTION: Basic diagnostics for the models.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list = ls(all = TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

library(tidyverse)
library(bayestestR)
library(rstan)

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

# Table S1: Global model coefficients ##########

tableS1 <- describe_posterior(m1) %>%
  select(-c(pd, CI, ROPE_CI, ROPE_low, ROPE_high, ROPE_Percentage)) %>%
  mutate(Parameter = gsub("b_", "", Parameter)) %>%
  mutate(across(Median:CI_high,  ~ round(., digits = 3))) %>%
  mutate(across(Rhat:ESS,  ~ round(., digits = 2)))

write.csv(tableS1, "Results/tables/TableS1")

# Table S2: System model coefficients ##########

sys_model_ls <- list(m_fr, m_mar, m_ter) %>%
  `names<-`(c("Freshwater", "Marine", "Terrestrial"))

tableS2 <- lapply(seq_along(sys_model_ls), function(mod) {
  describe_posterior(sys_model_ls[[mod]]) %>%
    select(-c(pd, CI, ROPE_CI, ROPE_low, ROPE_high, ROPE_Percentage)) %>%
    mutate(Parameter = gsub("b_", "", Parameter)) %>%
    mutate(across(Median:CI_high,  ~ round(., digits = 3))) %>%
    mutate(across(Rhat:ESS,  ~ round(., digits = 2))) %>%
    mutate(system = names(sys_model_ls)[[mod]])
  }) %>%
  bind_rows()

write.csv(tableS2, "Results/tables/TableS2")

# Table S3: Taxon model coefficients ##########

taxon_model_ls <- list(m_am, m_bi, m_f, m_ma, m_rep) %>%
  `names<-`(c("Amphibian", "Bird", "Fish", "Mammal", "Reptile"))

tableS3 <- lapply(seq_along(taxon_model_ls), function(mod) {
  describe_posterior(taxon_model_ls[[mod]]) %>%
    select(-c(pd, CI, ROPE_CI, ROPE_low, ROPE_high, ROPE_Percentage)) %>%
    mutate(Parameter = gsub("b_", "", Parameter)) %>%
    mutate(across(Median:CI_high,  ~ round(., digits = 3))) %>%
    mutate(across(Rhat:ESS,  ~ round(., digits = 2))) %>%
    mutate(system = names(taxon_model_ls)[[mod]])
  }) %>%
  bind_rows()

write.csv(tableS3, "Results/tables/TableS3")
