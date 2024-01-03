
all_threats = c("pollution","habitatl","climatechange","invasive", "exploitation","disease")

threatcols <- colnames(intvadd_dat3)[grepl(paste(all_threats,collapse = "|"),colnames(intvadd_dat3))] 

postdraws <- threat_post_draws(model = mod_intvadd2,
                                      modelled_data = intvadd_dat3,
                                      threat_comns = threatcols) #estimate posterior draws for all threat singles and combinations

post_dydx <- do.call("rbind",lapply(all_threats,function(x){

  # out <- postdraws %>%
  #   subset(grepl(x,threat)) %>%
  #   reframe(.value = mean(diff(.value)/diff(time)),.by = c(threat,.draw)) %>%
  #   mutate(int_group = ifelse(grepl("\\.",threat),"combined","single"),
  #          threat_group = x)
  
  out <- postdraws %>%
    subset(grepl(x,threat)) %>% #subset to focal threat
    reframe(.value = mean(diff(.value)/diff(time)),.by = c(threat,.draw)) %>% #for each posterior timeseries, estimate the first derivative 
    mutate(int_group = ifelse(grepl("\\.",threat),"combined","single"))%>% #create grouping column for whether the threat is singular or a combination
    mutate(threat_group = x) %>% #overall threat grouping
    group_by(threat) %>%
    filter(!any(.value >= abs(0.5))) %>% #drop threats with highly uncertain confidence intervals. Typically certain three way interactions
    ungroup() %>%
    subset(!grepl("[[:lower:]]*\\.[[:lower:]]*\\.",threat)) 

  return(out)
  }))

dydx_interval <- post_dydx  %>%
  mutate(.chain = 1, .iteration = .draw) %>% #add additional columns required by ggdist
  #subset(!grepl("[[:lower:]]*\\.[[:lower:]]*\\.",threat)) %>%
  group_by(threat_group,int_group) %>%
  ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".chain", ".iteration", ".draw", "threat")) #extract distribution information


ggplot(data = post_dydx , 
       aes(x = .value,y=int_group)) +
  tidybayes::stat_slab(alpha=0.5) +
  ggdist::geom_pointinterval(data = dydx_interval,
                             aes(xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  xlab("Population trend") + 
  ylab("Threat") + 
  coord_cartesian(xlim = c(-0.5,0.5)) +
  facet_wrap(~threat_group) +
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
        plot.title = element_text(hjust = 0.5))
