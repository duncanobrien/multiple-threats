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

# Load the raw LPI data

load("Data/LivingPlanetData2.RData")

mod_dat_full <- prepare_data(dd_long, duration = 10) |>
  arrange(ID)
# Set the priors 

priors <- c(prior(normal(0, 1), class = b),
            prior(exponential(1), class = sd),
            prior(normal(0,0.25), class = ar))

# Subset the locations 

locs <- distinct(mod_dat_full[,c("Longitude","Latitude")]) %>%
  mutate(Site = paste(Latitude,Longitude,sep = "_"))

# Get the spatial distance matrix 

spa_mat_trim <- as.matrix(geosphere::distm(locs[,c("Longitude","Latitude")], 
                                               fun = geosphere::distHaversine))/1000 #km distance between sites
# Normalise the distance

spa_mat_trim <- norm_range(spa_mat_trim) 

# Get the absolute value 

spa_mat_trim <- abs(spa_mat_trim - 1)

# Set the col and rownames 

colnames(spa_mat_trim) <- locs$Site
rownames(spa_mat_trim) <- locs$Site

# Create the formula

rhs <- paste0(paste("scaled_year*",
                        colnames(mod_dat_full)[
                          grepl(paste(c("pollution","habitatl","climatechange","invasive", "exploitation","disease"),
                                      collapse = "|"),
                                colnames(mod_dat_full))],
                        sep ="",collapse = " + "),
                  " + (-1 + scaled_year|SpeciesName) + (-1 + scaled_year|series) +
               (0 + scaled_year|gr(Site, cov = spa_mat)) + 1")

# Combine y_centered and rhs into a model formula

form <- as.formula(paste("y_centered", "~", rhs)) 

# Run the model 

mod_glob <- brm(bf(form #include realm/spp as slopes, x intercepts
                  ,autocor = ~ar(time = time,gr = series,p=1)),
               data = mod_dat_full, 
               data2 = list(spa_mat = spa_mat_trim),
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
saveRDS(mod_glob,"Results/models/mod_global_rerun2.RDS")

saveRDS(mod_glob,"Results/models/mod_global_rerun.RDS")

