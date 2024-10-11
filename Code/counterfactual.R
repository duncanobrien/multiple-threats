# --------------------------------------------------------------------------------------- #
# - FILE NAME:   counterfactual.R         
# - DATE:        03/01/2024
# - DESCRIPTION: Code to test scenarios where the threats have been removed.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila@ub.edu), Duncan O'Brien (duncan.a.obrien@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(brms)
library(tidyverse)
library(ggridges)
library(tidybayes)
library(parallel)
library(doSNOW)
library(patchwork)
library(rstan)

# Load the model 

m1 <- read_rds("Results/models/mod_global_rerun.RDS")

# Load the functions

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

ggsave("Results/Figure 4.pdf", g4, 
       width = 10, height = 6)

## Two threats scenarios -------------------------------------------------------

(p2 <- counter_fac_data %>% 
    # add a variable counting the number of threats
    mutate(number=str_count(counterfac, '\\.')+1,
           counterfac = gsub("habitatl", "habitat loss", counterfac),
           counterfac = gsub("climatechange", "climate change", counterfac)) %>%
    group_by(counterfac) %>% 
    mutate(pos=median(mn)) %>% 
    ungroup() %>%   
    filter(number==2) %>% 
    # Start the plot
    ggplot(aes(x = mn,
               y=reorder(counterfac, pos),
               fill=counterfac, colour=counterfac)) +
   stat_density_ridges(quantile_lines = TRUE, quantiles = 2,
                       scale = 1,vline_size = 1.3,
                       rel_min_height = 0.01,
                       bandwidth = 0.01, alpha=.5) +
    scale_fill_manual(values=MetBrewer::met.brewer(name="Tiepolo", n=15, type="continuous"))+
    scale_colour_manual(values=MetBrewer::met.brewer(name="Tiepolo", n=15, type="continuous"))+
    geom_vline(xintercept = 0,
               linetype = "solid", 
               colour="grey50",linewidth=1) +
    geom_vline(xintercept = subset(scenarios_mean,counterfac == "none")$m,
               linetype = "dashed", 
               colour="grey50",linewidth=1) +
    xlab(expression(paste("Population trend (",Delta,"y","/",Delta,"x)"))) + 
    ylab("Counterfactual") +
    coord_cartesian(xlim = c(-0.2,0.2))+
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
          plot.title = element_text(hjust = 0.5),
          legend.position = "none"))

## Three threats scenarios -----------------------------------------------------

(p3 <- counter_fac_data %>% 
    # add a variable counting the number of threats
    mutate(number=str_count(counterfac, '\\.')+1,
           counterfac = gsub("habitatl", "habitat loss", counterfac),
           counterfac = gsub("climatechange", "climate change", counterfac)) %>%
    group_by(counterfac) %>% 
    mutate(pos=median(mn)) %>% 
    ungroup() %>%   
    filter(number==3) %>% 
    # Start the plot
    ggplot(aes(x = mn,
               y=reorder(counterfac, pos),
               fill=counterfac, 
               colour=counterfac)) +
   stat_density_ridges(quantile_lines = TRUE, quantiles = 2,
                       scale = 1,vline_size = 1.3,
                       rel_min_height = 0.01,
                       bandwidth = 0.01, alpha=.5) +
    scale_fill_manual(values=MetBrewer::met.brewer(name="Nattier", n=18, type="continuous"))+
    scale_colour_manual(values=MetBrewer::met.brewer(name="Nattier", n=18, type="continuous"))+
    # geom_boxplot(outlier.shape = NA)+
    # tidybayes::stat_halfeye(aes(group = counterfac)) +
   geom_vline(xintercept = 0,
              linetype = "solid", 
              colour="grey50",linewidth=1) +
   geom_vline(xintercept = subset(scenarios_mean,counterfac == "none")$m,
              linetype = "dashed", 
              colour="grey50",linewidth=1) +
   xlab(expression(paste("Population trend (",Delta,"y","/",Delta,"x)"))) + 
    ylab("") +
    coord_cartesian(xlim = c(-0.2,0.2))+
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
          plot.title = element_text(hjust = 0.5),
          legend.position = "none"))

# Combine it

(figureS13 <- p2+ p3 +
    plot_annotation(tag_levels = c('a')))

# Save it

ggsave("Results/FigureS13.pdf", figureS13,
       width = 14, height = 12)

# # Estimate the difference in population trend
# 
# counterfac_diff <- do.call("rbind",
#                            lapply(threat_cols,
#                                   function(x){
#                                     counter_fac_data %>%
#                                       # Filter
#                                       subset(counterfac %in% c(x,"none")) %>%
#                                       # as "none" is reference level, sets none 
#                                       # first value of group
#                                       arrange(series,.draw,counterfac) %>%
#                                       # find difference in derivatives between 
#                                       # additive and interactive threats
#                                       reframe(.value = .value[1]-.value[2], 
#                                               .by = c(series,.draw)) %>% 
#                                       #name threat combination
#                                       mutate(counterfac = x) 
#   
# }))
# 
# # Let's plot it 
# 
# x <- counterfac_diff %>% 
#   group_by(counterfac) %>% 
#   summarise(m=mean(.value))
# 
# counterfac_diff %>% 
#   ggplot(aes(x=counterfac, y=.value))+
#   geom_boxplot(outlier.shape = NA)+
#   ylim(-0.06, 0.06)
# 
# # Older version of the code ###############
# # 
# # # Extract the data used to build the model
# # 
# # pop_data <- m1$data 
# # 
# # # Store the names of threats
# # 
# # threats <- colnames(pop_data)[grepl(paste(c("pollution","habitatl",
# #                                                "climatechange","invasive", 
# #                                                "exploitation","disease"),
# #                                              collapse = "|"),
# #                                        colnames(pop_data))] 
# # 
# # # Obtain the predicted value for each population 
# # 
# # pop_pred <- posterior_epred(m1,
# #                             newdata = m1$data,
# #                             re.form = NULL,
# #                             ndraws = 1000) %>% #extract posterior draws for the above data grid
# #   t() %>%
# #   as.data.frame() %>%
# #   split(m1$data$series) %>%
# #   sapply(FUN = function(y){ 
# #     apply(y,MARGIN = 2, FUN = function(x){mean(diff(x)/diff(seq_along(x)))})}) %>%
# #   as.data.frame() %>%
# #   mutate(.draw = 1:n()) %>%
# #   pivot_longer(-.draw,names_to = "series",values_to = ".value") %>%
# #   merge(y = dplyr::distinct(dplyr::select(m1$data,-c(y_centered,scaled_year,time))),
# #         by = "series") %>% 
# #   # Join with the original data to get a threats column
# #   left_join(dd_long %>% 
# #               # Keep only one ID
# #               distinct(ID,.keep_all = T) %>% 
# #               # Select only the columns we are interested in
# #               select(threats, ID) %>% 
# #               # Transform ID in a character so it can be compared against series
# #               mutate(ID=as.character(ID)), 
# #             by=c("series"="ID")) %>% 
# #   # Count the threats and correct the column threats to match the initial names
# #   mutate(number=str_count(threats,"[A-Z]"),
# #          threats=tolower(threats), 
# #          threats=gsub(":", ".", threats),
# #          threats=gsub("habitat loss", "habitatl", threats), 
# #          threats=gsub("invasive spp/genes", "invasive", threats),
# #          threats=gsub("climate change","climagechange", threats),
# #          threats=replace_na(threats,"none"))
# # 
# # # Obtain the mean coefficient for each threat 
# # 
# # mean_trends <- m1 %>% 
# #   gather_draws(`b_.*`, regex = TRUE) %>%
# #   mean_qi(.width = .95) %>%
# #   mutate(.variable = gsub("b_", "", .variable), 
# #          .variable = gsub("threats", "", .variable), 
# #          .variable = gsub("scaled_year", "year", .variable),
# #          .variable=gsub("1", "", .variable),
# #          number=str_count(.variable, '\\.')+1) %>% 
# #   filter(grepl("year:", .variable))
# # 
# # # Substract the mean coefficient from the predicted values
# # # We will parallelise the code, so this will take multiple steps
# # 
# # # Define the threats
# # 
# # threats <- colnames(pop_data)[3:39]
# # 
# # # Define the number of cores to use
# # 
# # clus  <- makeCluster(detectCores() - 1)
# # 
# # # register the cluster with doSNOW
# # 
# # registerDoSNOW(clus)
# # 
# # # load in the options using .options.snow
# # 
# # trend_scenarios <- foreach(i = threats, 
# #                         .combine = 'rbind',
# #                         .packages = c("tidyverse", "naniar")) %dopar% {
# #                           
# #                           # Extract the number of threats
# #                           
# #                           number_threats <- length(str_split(i, "\\.")[[1]])
# #                           
# #                           # If the number of threats is equal to 1
# #                           
# #                           if(number_threats==1){
# #                             # Create a new dataframe were we substract the mean 
# #                             # effect of the threats on the population trends 
# #                             new <- pop_pred %>% 
# #                               filter(grepl(i,threats)) %>% # Filter by the threat
# #                               # Substract the mean effect of the single threat
# #                               mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",i)]) %>% 
# #                               # Bind to the previous dataset
# #                               rbind(pop_pred  %>% 
# #                                       filter(!grepl(i,threats)) %>% 
# #                                       mutate(new_trend=NA)) %>% 
# #                               # Add the name of the scenario 
# #                               mutate(scenario=paste("No", i)) 
# #                             
# #                             # If the number of threats is equal to 2
# #                             
# #                             } else if(number_threats==2){
# #                               # Split the names of the interacting threats
# #                               x2 <- unlist(strsplit(i, 
# #                                                     ".", 
# #                                                     fixed=T))
# #                               # Create a new dataset where we remove the effect
# #                               # of the two threats, when interacting together
# #                               new1 <- pop_pred  %>%
# #                                 # Filter for those populations exposed to both
# #                                 # threats
# #                                 filter(grepl(x2[1],threats), 
# #                                        grepl(x2[2],threats), number!=1) %>% 
# #                                 # Extract the effect of both threats
# #                                 mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",i)])
# #                               # Create a second dataset where we remove the effect
# #                               # of the first threat, only 
# #                               new2 <- pop_pred  %>% 
# #                                 # Filter for the first threat, only those with
# #                                 # more and less than two threats and that do not
# #                                 # overlap with the series in new1
# #                                 filter(grepl(x2[1],threats), 
# #                                        number!=2, !grepl(i,threats),
# #                                        !series%in%new1$series) %>%
# #                                 # Substract the effect of the single threat
# #                                 mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",x2[1])]) 
# #                               # Same as above but with the second threat
# #                               new <- pop_pred %>% 
# #                                 # Notice we also filer out series in new1 and 
# #                                 # new2
# #                                 filter(grepl(x2[2],threats), number!=2,
# #                                        !grepl(i,threats), 
# #                                        !series%in%new1$series,
# #                                        !series%in%new2$series) %>% 
# #                                 mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",x2[2])])  %>% 
# #                                 # Bind with new1, new2, and pop_pred 
# #                                 rbind(new1) %>% 
# #                                 rbind(new2) 
# #                               
# #                               new <- new %>% 
# #                                 rbind(pop_pred  %>% 
# #                                       filter(!series%in%new$series) %>% 
# #                                         mutate(new_trend=NA)) %>% 
# #                               mutate(scenario=paste("No", i))
# #                               
# #                               # If the number of threats is equal to 3
# #                               
# #                               }else if(number_threats==3){
# #                                 # Separate the threats
# #                                 x3 <- unlist(strsplit(i, 
# #                                                       ".", 
# #                                                       fixed=T))
# #                                 # Create new names for the two interactive threats
# #                                 two1 <- paste0(x3[1], ".", x3[2]) 
# #                                 two2 <- paste0(x3[2], ".", x3[3]) 
# #                                 two3 <- paste0(x3[1], ".", x3[3])
# #                                 # Create a data frame where we subtract the 
# #                                 # effect of three threats
# #                                 new1 <- pop_pred  %>% 
# #                                   filter(grepl(x3[1],threats), 
# #                                          grepl(x3[2],threats),
# #                                          grepl(x3[3],threats), 
# #                                          number==3) %>% 
# #                                   mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",i)]) 
# #                                 # Subset the population trends with two threats
# #                                 # First and second threat
# #                                 two_threat1 <- mean_trends %>%  
# #                                   filter(grepl(x3[1], .variable),
# #                                          grepl(x3[2], .variable),
# #                                          number==2)
# #                                 # Second and third threat
# #                                 two_threat2 <- mean_trends %>%  
# #                                   filter(grepl(x3[2], .variable),
# #                                          grepl(x3[3], .variable),
# #                                          number==2)
# #                                 # First and third threat
# #                                 two_threat3 <- mean_trends %>%  
# #                                   filter(grepl(x3[1], .variable),
# #                                          grepl(x3[3], .variable),
# #                                          number==2)
# #                                 # If the two threats do not exist 
# #                                 if(dim(two_threat1)[1]==0){
# #                                   # When there is no value we create a fake one
# #                                   # with new_trend == NA
# #                                   new2 <- tryCatch({pop_pred  %>% 
# #                                       filter(grepl(x3[1], threats),
# #                                              grepl(x3[2], threats),
# #                                              !ID%in%new1$ID) %>% 
# #                                       mutate(new_trend=NA)},
# #                                       error=function(e){})
# #                                   # When there are some we do the usual procedure
# #                                   }else{new2 <- tryCatch({pop_pred  %>% 
# #                                       filter(grepl(x3[1], threats),
# #                                              grepl(x3[2], threats),
# #                                              !series%in%new1$series) %>% 
# #                                       mutate(new_trend=.value-two_threat1$.value)},
# #                                       error=function(e){})}
# #                                 # Same here with other two threats
# #                                 if(dim(two_threat2)[1]==0){
# #                                   new3 <- tryCatch({pop_pred  %>% 
# #                                       filter(grepl(x3[2], threats),
# #                                              grepl(x3[3], threats),
# #                                              !series%in%new1$series) %>% 
# #                                       mutate(new_trend=NA)},
# #                                       error=function(e){})
# #                                   }else{new3 <- tryCatch({pop_pred  %>% 
# #                                       filter(grepl(x3[2], threats),
# #                                              grepl(x3[3], threats),
# #                                              !series%in%new1$series) %>% 
# #                                       mutate(new_trend=.value-two_threat2$.value)},
# #                                       error=function(e){})}
# #                                 # Same here with other two threats
# #                                 if(dim(two_threat3)[1]==0){ 
# #                                   new4 <- tryCatch({pop_pred  %>% 
# #                                       filter(grepl(x3[1], threats),
# #                                              grepl(x3[3], threats),
# #                                              !series%in%new1$series) %>% 
# #                                       mutate(new_trend=NA)},
# #                                       error=function(e){})
# #                                 }else{
# #                                   new4 <- tryCatch({pop_pred  %>% 
# #                                       filter(grepl(x3[1], threats),
# #                                              grepl(x3[3], threats),
# #                                              !series%in%new1$series) %>% 
# #                                       mutate(new_trend=.value-two_threat3$.value)},
# #                                       error=function(e){})  
# #                                 }
# #                                 
# #                                 # Now we remove the single effects
# #                                 # For first threat
# #                                 new5 <- pop_pred  %>% 
# #                                   filter(grepl(x3[1], threats),
# #                                          !series%in%new1$series, 
# #                                          !series%in%new2$series,
# #                                          !series%in%new4$series) %>% 
# #                                   mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",x3[1])])
# #                                 # For second threat
# #                                 new6 <- pop_pred  %>% 
# #                                   filter(grepl(x3[2], threats),
# #                                          !series%in%new1$series, 
# #                                          !series%in%new2$series,
# #                                          !series%in%new3$series) %>% 
# #                                   mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",x3[2])])
# #                                 # Final subset
# #                                 new <- pop_pred  %>% 
# #                                   filter(grepl(x3[3], threats),
# #                                          !series%in%new1$series, 
# #                                          !series%in%new3$series,
# #                                          !series%in%new4$series) %>% 
# #                                   mutate(new_trend=.value-mean_trends$.value[mean_trends$.variable==paste0("year:",x3[3])]) %>% 
# #                                   # Bind to the previous datasets
# #                                   rbind(new1) %>% 
# #                                   rbind(new2) %>% 
# #                                   rbind(new3) %>% 
# #                                   rbind(new4) %>% 
# #                                   rbind(new5) %>% 
# #                                   rbind(new6) 
# #                                 new <- new %>% 
# #                                   rbind(pop_pred  %>% 
# #                                           filter(!series%in%new$series)) %>% 
# #                                   # Add the scenario name
# #                                   mutate(scenario=paste("No", i)) 
# #                                 }
# #                         }
# # 
# # # Stop the parallelisation
# # 
# # stopCluster(clus)
# # 
# # # Save the object 
# # 
# # save(trend_scenarios, "Results/counterfactuals.RData") 
# # 
# # # # Add the control scenario 
# # # 
# # # scenarios_data <- 
# # #   pop_perd %>% mutate(scenario="No Management") %>% 
# # #   bind_rows(mu_scenarios) %>% 
# # #   filter(threats!="none") %>% 
# # #   mutate(prop=((log(lambda)-Estimate)/abs(Estimate))*100,
# # #          prop=round(prop,2),
# # #          number_scen= str_count(scenario,"[A-Z]")) %>% 
# # #   group_by(scenario) %>% 
# # #   mutate(pos=median(lambda), 
# # #          pos2=median(prop)) %>% 
# # #   ungroup()  
# # # 
# # # # Add the control scenario 
# # # 
# # # scenarios_median <- scenarios_data %>% 
# # #   group_by(scenario) %>% 
# # #   summarise(mcontrol=median(Estimate),
# # #             mcounter=median(log(lambda))) %>% 
# # #   mutate(prop=((mcounter-(mcontrol))/abs(mcontrol))*100,
# # #          prop=round(prop,2))
# # 
# # # Panel b: projection ----------------------------------------------------------
# # 
# # # Predict under no scenario
# # 
# # pop_perd <- brms::posterior_epred(m1,
# #                                   newdata = m1$data %>% 
# #                                     filter_at(threat_cols, any_vars(. != "0")),
# #                                   re.form = NULL,
# #                                   incl_autocor = FALSE,
# #                                   sort = TRUE, 
# #                                   ndraws = 250) %>% 
# #   as.data.frame() %>%
# #   mutate(.draw = 1:NROW(.)) %>%
# #   #extract posterior draws for the data used to create the model
# #   pivot_longer(-.draw,names_to = "index",values_to = ".value") %>%
# #   cbind(m1$data %>%
# #           filter_at(threat_cols, any_vars(. != "0")) %>% 
# #           dplyr::select(series, time), row.names = NULL) %>% 
# #   mutate(counterfac="none") 
# # 
# # # Use our function to predict populations under different scenarios
# # 
# # counter_fac_proj <- threat_counterfac_pred(m1,
# #                                             threat_comns = threat_cols,
# #                                             re.form = NULL,
# #                                             ndraws = 250,
# #                                             n.cores = 4) %>%
# #   # Join with the none counterfactual scenario that we just created
# #   rbind(pop_perd) %>%
# #   # Made "none" as the first level of the counterfactual 
# #   mutate(counterfac = fct_relevel(counterfac, "none"),
# #          number=str_count(counterfac, '\\.')+1,
# #          counterfac=gsub("none", "None", counterfac),
# #          counterfac = gsub("habitatl", "habitat loss", counterfac),
# #          counterfac = gsub("climatechange", "climate change", counterfac)) #Create a new variable for number of threats
# # 
# # # Mean trajectory
# # 
# # counter_fac_mean <- counter_fac_proj %>% 
# #   group_by(counterfac, time) %>% 
# #   summarise(m=mean(.value))
# # 
# # counter_interval <-counter_fac_proj %>% 
# #    group_by(counterfac, time) %>% 
# #    ggdist::median_qi(.width = c(.95), 
# #                      .exclude = c(".draw", "index", "series", "number")) #extract distribution information
# # 
# # 
# # # Plot it 
# # 
# # data <- counter_fac_mean %>%
# #   mutate(number=str_count(counterfac, '\\.')+1) %>% 
# #   filter(number==1) 
# # 
# # data %>% 
# #   ggplot(aes(x=time, y=m, group=counterfac)) +
# #   geom_line(linewidth=1, 
# #             aes(colour=counterfac))+
# #   scale_colour_manual(name = NULL, values=threat_palette)+
# #   geom_hline(yintercept=0, linetype="dashed")+
# #   #xlim(0,50)+
# #   #ylim(-0.5,0.5)+
# #   labs(x="Years", y="Relative abundance")+
# #   theme(legend.position = "bottom")
# # 
# # data_int <- counter_interval %>%
# #   mutate(number=str_count(counterfac, '\\.')+1) %>% 
# #   filter(number==1) 
# # 
# # data_int %>%  
# #   mutate(counterfac=gsub("none", "None", counterfac),
# #          counterfac = gsub("habitatl", "habitat loss", counterfac),
# #          counterfac = gsub("climatechange", "climate change", counterfac)) %>% 
# #   filter(counterfac!="None") %>% 
# #   ggplot(aes(x=time, y=.value, group=counterfac)) +
# #   geom_ribbon(aes(ymin = .lower, ymax = .upper, 
# #                   fill=counterfac), 
# #               alpha=0.2, colour=NA)+
# #   geom_line(linewidth=1, 
# #             aes(colour=counterfac))+
# #   scale_colour_manual(name = NULL, values=threat_palette)+
# #   scale_fill_manual(name = NULL, values=threat_palette)+
# #   geom_hline(yintercept=0, linetype="dashed")+
# #   #xlim(0,50)+
# #   #ylim(-0.5,0.5)+
# #   labs(x="Years", y="Relative abundance")+
# #   facet_wrap(~counterfac)+
# #   theme(legend.position = "none")
