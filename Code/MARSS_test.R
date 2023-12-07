require(MARSS)
require(tidyverse)
load("Data/LivingPlanetData2.RData")

dd_long2 <- dd_long %>%
  group_by(ID) %>%  
  # Calculate population change
  mutate(popchange=log(Count+(max(Count, na.rm=T) /100)),
         #ifelse(Count<1, log(Count/(max(Count,na.rm = T))*100+1), 
         #        log(Count+1)),
         Protected_status=gsub(" .*", "", Protected_status)) %>% 
  # Remove any groupings we have created in the pipe
  ungroup()

# Add the biogeographical region

dd_long2 <- dd_long2 %>%
  mutate(T_realm=na_if(T_realm, "NULL"),
         M_realm=na_if(M_realm, "NULL"),
         FW_realm=na_if(FW_realm, "NULL"),
         Realm=coalesce(T_realm, M_realm, FW_realm))

# Find the data ids

dd_final <- dd_long2 %>%
  group_by(ID) %>%
  drop_na(popchange) %>%
  mutate(start=min(Year),
         end=max(Year),
         Duration=end-start) %>% 
  filter(Duration>9)

# Spread the data again

pops <- dd_final %>%
  mutate(year = as.factor(as.character(Year))) %>% 
  dplyr::select(ID, SpeciesName, Class, Order,
                System, Latitude, Longitude,
                Region,
                Protected_status, n.threat, 
                Primary_threat,
                Secondary_threat,
                Tertiary_threat,
                Managed,
                Realm,
                threats, Duration, 
                Year, year, popchange) %>%
  drop_na(popchange) %>%
  group_by(ID) %>%
  dplyr::select(-Year) %>%
  filter(length(unique(year))>=5) %>% 
  pivot_wider(names_from = year, values_from = popchange)

pops_data <- pops %>%
  dplyr::select(ID, SpeciesName, Class, Order,
                System, Latitude, Longitude,
                Region,
                Protected_status, n.threat, 
                Primary_threat,
                Secondary_threat,
                Tertiary_threat,
                Managed,
                Realm,
                threats, Duration) 

# Compile all time-series into a list to estimate state-space population trends##
pop <- lapply(1:NROW(pops),function(x){
  data <- pops[x, 18:85]
  Year <- colnames(data)[is.na(data) == FALSE]
  data <- data[Year]
  N <- as.vector(t(data))
  return(data.frame(pops$ID[i], Year, N))
  }) %>%
  `names<-`(c(pops$ID))

##########################################################################################
## MARSS ##
##########################################################################################
sample_ts <- 1:10
tt_data <- as.matrix(pops[sample_ts,18:85])
years <- as.numeric(colnames(tt_data))
id <- pops$ID[sample_ts]
temp <- matrix(rnorm(50, seq(0,1,1/50), 0.1),nrow=length(sample_ts))

Z.model <- factor(1:length(sample_ts))
R.model <- "diagonal and equal"
C.model <- "unequal"

matrix(c("temp1","temp2"),2,1)

tt <- MARSS(tt_data,model = list(Z = "identity",
                                 R = "diagonal and unequal",
                                 D = "unconstrained",
                                 d = t(as.matrix(years))),
            control = list(maxit=10000))

ggplot2::autoplot(tt, plot.type = "xtT")
trends_tt <- tidy(tt) %>%
  filter(grepl("D.X.",term))


fulldat <- lakeWAplanktonTrans
covariates <- rbind(
  Temp = fulldat[1:length(years), "Temp"],
  TP = fulldat[1:length(years), "TP"]
)
