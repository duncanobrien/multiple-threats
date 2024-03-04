# --------------------------------------------------------------------------------------- #
# - FILE NAME:   realm_models.R         
# - DATE:        11/02/2024
# - DESCRIPTION: Code to run the individual models for each realm.
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

priors <- c(prior(normal(-0.005900686, 0.26), class = b),
            prior(exponential(1), class = sd),
            prior(normal(0,0.25), class = ar))

# Realm specific models --------------------------------------------------------
## Marine ---- 
# Filter for Marine 

mod_dat_mar <- mod_dat_full %>% filter(System=="Marine") %>%
  select_if(function(x)!all(x == "0"))

# Subset the locations 

mar_locs <- distinct(mod_dat_mar[,c("Longitude","Latitude")]) %>%
  mutate(Site = paste(Latitude,Longitude,sep = "_"))

# Get the spatial distance matrix 

spa_mat_trim_mar <- as.matrix(geosphere::distm(mar_locs[,c("Longitude","Latitude")], 
                                                fun = geosphere::distHaversine))/1000 #km distance between sites
# Normalise the distance

spa_mat_trim_mar <- norm_range(spa_mat_trim_mar) 

# Get the absolute value 

spa_mat_trim_mar <- abs(spa_mat_trim_mar - 1)

# Set the col and rownames 

colnames(spa_mat_trim_mar) <- mar_locs$Site
rownames(spa_mat_trim_mar) <- mar_locs$Site

# Create the formula

rhs_mar <- paste0(paste("scaled_year*",
                         colnames(mod_dat_mar)[
                           grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                                       collapse = "|"),
                                 colnames(mod_dat_mar))],
                         sep ="",collapse = " + "),
                   " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) +
               (0 + scaled_year|gr(Site, cov = spa_mat)) + -1")
# Combine y_centered and rhs into a model formula

form_mar <- as.formula(paste("y_centered", "~", rhs_mar)) 

# Run the model 

mod_mar <- brm(bf(form_mar #include realm/spp as slopes, x intercepts
                   ,autocor = ~ar(time = time,gr = series,p=1)),
                data = mod_dat_mar, 
                data2 = list(spa_mat = spa_mat_trim_mar),
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

saveRDS(mod_mar,"Results/models/mar_mod.RDS")

## Freshwater ---- 
# Filter for Freshwater

mod_dat_fre <- mod_dat_full %>% filter(System=="Freshwater") %>%
  select_if(function(x)!all(x == "0"))

# Subset the locations 

fre_locs <- distinct(mod_dat_fre[,c("Longitude","Latitude")]) %>%
  mutate(Site = paste(Latitude,Longitude,sep = "_"))

# Get the spatial distance matrix 

spa_mat_trim_fre <- as.matrix(geosphere::distm(fre_locs[,c("Longitude","Latitude")], 
                                               fun = geosphere::distHaversine))/1000 #km distance between sites
# Normalise the distance

spa_mat_trim_fre <- norm_range(spa_mat_trim_fre) 

# Get the absolute value 

spa_mat_trim_fre <- abs(spa_mat_trim_fre - 1)

# Set the col and rownames 

colnames(spa_mat_trim_fre) <- fre_locs$Site
rownames(spa_mat_trim_fre) <- fre_locs$Site

# Create the formula

rhs_fre <- paste0(paste("scaled_year*",
                        colnames(mod_dat_fre)[
                          grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                                      collapse = "|"),
                                colnames(mod_dat_fre))],
                        sep ="",collapse = " + "),
                  " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) +
               (0 + scaled_year|gr(Site, cov = spa_mat)) + -1")
# Combine y_centered and rhs into a model formula

form_fre <- as.formula(paste("y_centered", "~", rhs_fre)) 

# Run the model 

mod_fre <- brm(bf(form_fre #include realm/spp as slopes, x intercepts
                  ,autocor = ~ar(time = time,gr = series,p=1)),
               data = mod_dat_fre, 
               data2 = list(spa_mat = spa_mat_trim_fre),
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

saveRDS(mod_fre,"Results/models/fre_mod.RDS")

## Terrestrial ---- 
# Filter for Terrestrial

mod_dat_ter <- mod_dat_full %>% filter(System=="Terrestrial") %>%
  select_if(function(x)!all(x == "0"))

# Subset the locations 

ter_locs <- distinct(mod_dat_ter[,c("Longitude","Latitude")]) %>%
  mutate(Site = paste(Latitude,Longitude,sep = "_"))

# Get the spatial distance matrix 

spa_mat_trim_ter <- as.matrix(geosphere::distm(ter_locs[,c("Longitude","Latitude")], 
                                               fun = geosphere::distHaversine))/1000 #km distance between sites
# Normalise the distance

spa_mat_trim_ter <- norm_range(spa_mat_trim_ter) 

# Get the absolute value 

spa_mat_trim_ter <- abs(spa_mat_trim_ter - 1)

# Set the col and rownames 

colnames(spa_mat_trim_ter) <- ter_locs$Site
rownames(spa_mat_trim_ter) <- ter_locs$Site

# Create the formula

rhs_ter <- paste0(paste("scaled_year*",
                        colnames(mod_dat_ter)[
                          grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                                      collapse = "|"),
                                colnames(mod_dat_ter))],
                        sep ="",collapse = " + "),
                  " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) +
               (0 + scaled_year|gr(Site, cov = spa_mat)) + -1")

# Combine y_centered and rhs into a model formula

form_ter <- as.formula(paste("y_centered", "~", rhs_ter)) 

# Run the model 

mod_ter <- brm(bf(form_ter #include realm/spp as slopes, x intercepts
                  ,autocor = ~ar(time = time,gr = series,p=1)),
               data = mod_dat_ter, 
               data2 = list(spa_mat = spa_mat_trim_ter),
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

saveRDS(mod_ter,"Results/models/ter_mod.RDS")
