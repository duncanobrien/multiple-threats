#' Prepare Living Planet data for modelling.
#'
#' Converts 
#' @param data Dataframe of raw Living Planet data
#' @param duration Minimum length of time series to include for analysis (including missing years).
#' @param ratio Numeric. Minimum ratio of total length of non-missing data within a time series relative to the total number of years covered (including missing years)
#' @param latlong_digits Numeric. Number of digits to round Latitude and Longitude coordinates to.
#'
#' @returns Dataframe

prepare_data <- function(data, 
                         duration = 10,
                         #ratio = 0.5,
                         latlong_digits = 0){
  
  source("Code/utils/na.noncontiguous_fn.R")
  dd_long2 <- data %>%
    group_by(ID) %>%  
    # Calculate population change
    mutate(popchange=log(Count+(max(Count, na.rm=T) /100)),
           Protected_status=gsub(" .*", "", Protected_status)) %>% 
    # Remove any groupings we have created in the pipe
    ungroup()
  
  # Filter for ten or more years of data with at least 50% of the years with records
  
  # dd_final <- dd_long2 %>%
  #   group_by(ID) %>%
  #   drop_na(popchange) %>%
  #   mutate(n = n(),
  #          Duration=(max(Year)-min(Year))+1,
  #          Ratio=n/Duration) %>% 
  #   ungroup() %>% 
  #   filter(Duration>=duration,
  #          Ratio>=ratio) 
  
  dd_final <- dd_long2 %>%
    group_by(ID) %>%
    arrange(Year) %>%
    mutate(popchange = na.noncontiguous(popchange,duration)) %>% #find run of unbroken abundances which is at least duration time points long. This function converts all abundances outside of this run to NA
    ungroup() %>%
    drop_na(popchange) %>%
    mutate(Duration = n())
  
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
                  threats, Duration, 
                  Year, year, popchange) %>%
    drop_na(popchange) %>%
    group_by(ID) %>%
    dplyr::select(-Year) %>%
    pivot_wider(names_from = year, values_from = popchange)
  
  # Create threat combinations
  
  thrts_2 <-combn(c("pollution","habitatl","climatechange","invasive", 
                    "exploitation","disease"),2)
  thrts_3 <-combn(c("pollution","habitatl","climatechange",
                    "invasive", "exploitation","disease"),3) 
  
  # Further wrangling in to long format and prepare explanatory variables
  
  mod_dat_full <- pops %>%
    pivot_longer(matches("^[[:digit:]]+$"), #only pivot year columns which are the only digit containing column names
                 names_to = "year", 
                 values_to = "y") %>%
    mutate(year = as.numeric(year),
           series = paste(ID), #factor required for autocorrelation estimation
           threats = factor(ifelse(is.na(threats),"none",threats))) %>% #convert stressors to binary
    mutate(threats = fct_relevel(threats, "none")) %>%
    left_join(dplyr::select(dd_final,c(ID,pollution,habitatl,climatechange,
                                       invasive, exploitation,disease)),
              multiple = "first",by = "ID") %>% 
    mutate(none = ifelse(all(is.na(c(pollution,habitatl,climatechange,
                                     invasive, exploitation,disease))),
                         "none",NA)) %>%
    mutate(across(pollution:none,~ifelse(is.na(.x),"0","1"))) %>% #convert absence of stress to binary
    drop_na(y) %>% 
    group_by(ID) %>%
    mutate(scaled_year = c(scale(year,center = TRUE,scale = FALSE)), #center time to 0 for each timeseries
           time = seq_along(year),
           scaled_time = c(scale(time,center = TRUE,scale = FALSE)),
           #y_centered = y-(na.omit(y)[1]-0) #recenter y so that first value of timeseries is 0 (to allow all intercepts to be removed)
           y_centered = c(scale(y,center = TRUE,scale = FALSE)) #recenter y so that first value of timeseries is 0 (to allow all intercepts to be removed)
    ) %>% 
    ungroup(ID) %>%
    mutate(threats = as.character(threats)) %>%
    mutate(Taxon= ifelse(Class=="Holocephali"|Class=="Elasmobranchii" | 
                           Class=="Myxini"|Class=="Cephalaspidomorphi"|
                           Class=="Actinopterygii"|Class=="Sarcopterygii",
                         "Fish", 
                         ifelse(Class=="Aves", "Birds",
                                ifelse(Class=="Mammalia", 
                                       "Mammals",
                                       ifelse(Class=="Amphibia", 
                                              "Amphibians",
                                              ifelse(Class=="Reptilia",
                                                     "Reptiles", "NA")))))) 
  
  # Create a new dataframe with columns containing "0"/"1" for each combination of 
  # the threats. "0" = combination not present, "1" = combination present.
  
  mod_dat_full <- mod_dat_full %>% 
    bind_cols(purrr::pmap_dfc(.l = list(.x = thrts_2[1,], #threat 1
                                        .y = thrts_2[2,]),#threat 2
                              ~ mod_dat_full %>%
                                select({{.x}},{{.y}},ID) %>% #select just the necessary columns (threat 1, threat 2 and timeseries ID)
                                mutate(!!sym(paste(.x,.y,sep=".")) :=  
                                         ifelse(all(!!sym(.x) == "1") & all(!!sym(.y) == "1"),"1","0"),
                                       .by = ID) %>% #dynamically create new column of whether both threat 1 and threat 2 are 1's
                                select(4))%>%
                select_if(function(x)sum(x == "1") > 0) #drop columns where no observed combination of threats as will prevent model fitting
    ) %>%
    bind_cols(purrr::pmap_dfc(.l = list(thrts_3[1,], #threat 1
                                        thrts_3[2,], #threat 2
                                        thrts_3[3,]),#threat 3
                              function(.x,.y,.z) mod_dat_full %>%
                                select({{.x}},{{.y}},{{.z}},ID) %>% #select just the necessary columns (threat 1, threat 2 and timeseries ID)
                                mutate(!!sym(paste(.x,.y,.z,sep=".")) :=  
                                         ifelse(all(!!sym(.x) == "1") & all(!!sym(.y) == "1") & all(!!sym(.z) == "1"),"1","0"),
                                       .by = ID) %>% #dynamically create new column of whether both threat 1, threat 2 and threat 3 are 1's
                                select(5)) %>%
                select_if(function(x)sum(x == "1") > 0) #drop columns where no observed combination of threats
    ) %>%
    group_by(ID) %>% 
    ungroup() %>%
    select_if(function(x)!all(x == "0")) %>%
    mutate(Latitude = round(Latitude,digits = latlong_digits)) %>%
    mutate(Site = paste(Latitude,Longitude,sep = "_")) %>% 
    as.data.frame()
  
 return(mod_dat_full) 
}