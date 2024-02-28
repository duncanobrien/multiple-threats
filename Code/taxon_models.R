# --------------------------------------------------------------------------------------- #
# - FILE NAME:   taxon_models.R         
# - DATE:        11/02/2024
# - DESCRIPTION: Code to run the individual models for each taxonomic group.
# - AUTHORS: Pol Capdevila Lanzaco (pcapdevila@ub.edu), Duncan O'Brien (duncan.a.obrien@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(tidyverse)
library(brms)

# Function to normalise the spatial distance

norm_range <- function(x){
  (x-min(x))/(max(x)-min(x))
}

# Load the model data

load("Data/data_models.RData")

# Set the priors 

priors <- c(prior(normal(0, 1), class = b),
            prior(exponential(1), class = sd),
            prior(normal(0,0.25), class = ar))

# Taxon specific models ----------------------------------------------------------
## Amphibians ----
# Filter amphibians 

mod_dat_amph <- mod_dat_full %>% 
  filter(Taxon=="Amphibians") %>%
  select_if(function(x)!all(x == "0"))

# Subset the locations 

amph_locs <- distinct(mod_dat_amph[,c("Longitude","Latitude")]) %>%
  mutate(Site = paste(Latitude,Longitude,sep = "_"))

# Get the spatial distance matrix 

spa_mat_trim_amph <- as.matrix(geosphere::distm(amph_locs[,c("Longitude","Latitude")], 
                                                fun = geosphere::distHaversine))/1000 #km distance between sites
# Normalise the distance

spa_mat_trim_amph <- norm_range(spa_mat_trim_amph) 

# Get the absolute value 

spa_mat_trim_amph <- abs(spa_mat_trim_amph - 1)

# Set the col and rownames 

colnames(spa_mat_trim_amph) <- amph_locs$Site
rownames(spa_mat_trim_amph) <- amph_locs$Site

# Create the formula

rhs_amph <- paste0(paste("scaled_year*",
                        colnames(mod_dat_amph)[
                          grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                                      collapse = "|"),
                                colnames(mod_dat_amph))],
                        sep ="",collapse = " + "),
                  " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) +
               (0 + scaled_year|gr(Site, cov = spa_mat)) + -1")
# Combine y_centered and rhs into a model formula

form_amph <- as.formula(paste("y_centered", "~", rhs_amph)) 

# Run the model 

mod_amph <- brm(bf(form_amph #include realm/spp as slopes, x intercepts
                       ,autocor = ~ar(time = time,gr = series,p=1)),
                    data = mod_dat_amph, 
                    data2 = list(spa_mat = spa_mat_trim_amph),
                    family = gaussian(),
                    iter = 5000,
                    refresh=100,
                    backend = "cmdstanr",
                    silent = 0,
                    prior = priors,
                    chains = 4,
                    control=list(adapt_delta=0.975,max_treedepth = 12),
                    cores = 4)
# Save the model 

saveRDS(mod_amph,"Results/models/amph_mod.RDS")

## Birds ----
# Filter birds

mod_dat_bir <- mod_dat_full %>% filter(Taxon=="Birds") %>%
  select_if(function(x)!all(x == "0"))

# Subset the locations 

bir_locs <- distinct(mod_dat_bir[,c("Longitude","Latitude")]) %>%
  mutate(Site = paste(Latitude,Longitude,sep = "_"))

# Get the spatial distance matrix 

spa_mat_trim_bir <- as.matrix(geosphere::distm(bir_locs[,c("Longitude","Latitude")], 
                                                fun = geosphere::distHaversine))/1000 #km distance between sites
# Normalise the distance

spa_mat_trim_bir <- norm_range(spa_mat_trim_bir) 

# Get the absolute value 

spa_mat_trim_bir <- abs(spa_mat_trim_bir - 1)

# Set the col and rownames 

colnames(spa_mat_trim_bir) <- bir_locs$Site
rownames(spa_mat_trim_bir) <- bir_locs$Site

# Create the formula

rhs_bir <- paste0(paste("scaled_year*",
                         colnames(mod_dat_bir)[
                           grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                                       collapse = "|"),
                                 colnames(mod_dat_bir))],
                         sep ="",collapse = " + "),
                   " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) +
               (0 + scaled_year|gr(Site, cov = spa_mat)) + -1")
# Combine y_centered and rhs into a model formula

form_bir <- as.formula(paste("y_centered", "~", rhs_bir)) 

# Run the model 

mod_bir <- brm(bf(form_bir #include realm/spp as slopes, x intercepts
                   ,autocor = ~ar(time = time,gr = series,p=1)),
                data = mod_dat_bir, 
                data2 = list(spa_mat = spa_mat_trim_bir),
                family = gaussian(),
                iter = 5000,
                refresh=100,
                backend = "cmdstanr",
                silent = 0,
                prior = priors,
                chains = 4,
                control=list(adapt_delta=0.975,max_treedepth = 12),
                cores = 4)
# Save the model 

saveRDS(mod_bir,"Results/models/bir_mod.RDS")

## Bony fishes ----
# Filter Bony fishes

mod_dat_fish <- mod_dat_full %>% filter(Taxon=="Fish") %>%
  select_if(function(x)!all(x == "0"))

# Subset the locations 

fish_locs <- distinct(mod_dat_fish[,c("Longitude","Latitude")]) %>%
  mutate(Site = paste(Latitude,Longitude,sep = "_"))

# Get the spatial distance matrix 

spa_mat_trim_fish <- as.matrix(geosphere::distm(fish_locs[,c("Longitude","Latitude")], 
                                               fun = geosphere::distHaversine))/1000 #km distance between sites
# Normalise the distance

spa_mat_trim_fish <- norm_range(spa_mat_trim_fish) 

# Get the absolute value 

spa_mat_trim_fish <- abs(spa_mat_trim_fish - 1)

# Set the col and rownames 

colnames(spa_mat_trim_fish) <- fish_locs$Site
rownames(spa_mat_trim_fish) <- fish_locs$Site

# Create the formula

rhs_fish <- paste0(paste("scaled_year*",
                        colnames(mod_dat_fish)[
                          grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                                      collapse = "|"),
                                colnames(mod_dat_fish))],
                        sep ="",collapse = " + "),
                  " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) +
               (0 + scaled_year|gr(Site, cov = spa_mat)) + -1")

# Combine y_centered and rhs into a model formula

form_fish <- as.formula(paste("y_centered", "~", rhs_fish)) 

# Run the model 

mod_fish <- brm(bf(form_fish #include realm/spp as slopes, x intercepts
                  ,autocor = ~ar(time = time,gr = series,p=1)),
               data = mod_dat_fish, 
               data2 = list(spa_mat = spa_mat_trim_fish),
               family = gaussian(),
               iter = 5000,
               refresh=100,
               backend = "cmdstanr",
               silent = 0,
               prior = priors,
               chains = 4,
               control=list(adapt_delta=0.975,max_treedepth = 12),
               cores = 4)
# Save the model 

saveRDS(mod_fish,"Results/models/fish_mod.RDS")


## Mammals ----
# Filter mammals 

mod_dat_mam <- mod_dat_full %>% filter(Taxon=="Mammals") %>%
  select_if(function(x)!all(x == "0"))

# Subset the locations 

mam_locs <- distinct(mod_dat_mam[,c("Longitude","Latitude")]) %>%
  mutate(Site = paste(Latitude,Longitude,sep = "_"))

# Get the spatial distance matrix 

spa_mat_trim_mam <- as.matrix(geosphere::distm(mam_locs[,c("Longitude","Latitude")], 
                                                fun = geosphere::distHaversine))/1000 #km distance between sites
# Normalise the distance

spa_mat_trim_mam <- norm_range(spa_mat_trim_mam) 

# Get the absolute value 

spa_mat_trim_mam <- abs(spa_mat_trim_mam - 1)

# Set the col and rownames 

colnames(spa_mat_trim_mam) <- mam_locs$Site
rownames(spa_mat_trim_mam) <- mam_locs$Site

# Create the formula

rhs_mam <- paste0(paste("scaled_year*",
                         colnames(mod_dat_mam)[
                           grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                                       collapse = "|"),
                                 colnames(mod_dat_mam))],
                         sep ="",collapse = " + "),
                   " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) +
               (0 + scaled_year|gr(Site, cov = spa_mat)) + -1")

# Combine y_centered and rhs into a model formula

form_mam <- as.formula(paste("y_centered", "~", rhs_mam)) 

# Run the model 

mod_mam <- brm(bf(form_mam #include realm/spp as slopes, x intercepts
                   ,autocor = ~ar(time = time,gr = series,p=1)),
                data = mod_dat_mam, 
                data2 = list(spa_mat = spa_mat_trim_mam),
                family = gaussian(),
                iter = 5000,
                refresh=100,
                backend = "cmdstanr",
                silent = 0,
                prior = priors,
                chains = 4,
                control=list(adapt_delta=0.975,max_treedepth = 12),
                cores = 4)
# Save the model 

saveRDS(mod_mam,"Results/models/mam_mod.RDS")

## Reptiles ----
# Filter reptiles

mod_dat_rep <- mod_dat_full %>% filter(Taxon=="Reptiles") %>%
  select_if(function(x)!all(x == "0"))

# Subset the locations 

rep_locs <- distinct(mod_dat_rep[,c("Longitude","Latitude")]) %>%
  mutate(Site = paste(Latitude,Longitude,sep = "_"))

# Get the spatial distance matrix 

spa_mat_trim_rep <- as.matrix(geosphere::distm(rep_locs[,c("Longitude","Latitude")], 
                                                fun = geosphere::distHaversine))/1000 #km distance between sites
# Normalise the distance

spa_mat_trim_rep <- norm_range(spa_mat_trim_rep) 

# Get the absolute value 

spa_mat_trim_rep <- abs(spa_mat_trim_rep - 1)

# Set the col and rownames 

colnames(spa_mat_trim_rep) <- rep_locs$Site
rownames(spa_mat_trim_rep) <- rep_locs$Site

# Create the formula

rhs_rep <- paste0(paste("scaled_year*",
                         colnames(mod_dat_rep)[
                           grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                                       collapse = "|"),
                                 colnames(mod_dat_rep))],
                         sep ="",collapse = " + "),
                   " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) +
               (0 + scaled_year|gr(Site, cov = spa_mat)) + -1")

# Combine y_centered and rhs into a model formula

form_rep <- as.formula(paste("y_centered", "~", rhs_rep)) 

# Run the model 

mod_rep <- brm(bf(form_rep #include realm/spp as slopes, x intercepts
                   ,autocor = ~ar(time = time,gr = series,p=1)),
                data = mod_dat_rep, 
                data2 = list(spa_mat = spa_mat_trim_rep),
                family = gaussian(),
                iter = 5000,
                refresh=100,
                backend = "cmdstanr",
                silent = 0,
                prior = priors,
                chains = 4,
                control=list(adapt_delta=0.975,max_treedepth = 12),
                cores = 4)
# Save the model 

saveRDS(mod_rep,"Results/models/rep_mod.RDS")

