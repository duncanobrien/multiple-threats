require(brms)
require(tidyverse)

m1 <- read_rds("Results/models/mod_global.RDS")

source("Code/prep_data_grid_fn.R")
source("Code/threat_post_draws.R")

all_threats = c("pollution","habitatl","climatechange","invasive", "exploitation","disease")

threatcols <- colnames(m1$data)[grepl(paste(all_threats,collapse = "|"),colnames(m1$data))] 

additive_cols <- do.call("c",lapply(strsplit(threatcols,"[.]"),function(x){
  paste(x,collapse = " + ")
})) #create addtive columns. i.e. "threat1 + threat2"

target_cols <- unique(c(threatcols,additive_cols))

postdraws_intvadd <- threat_post_draws(model = m1,
                                       threat_comns = target_cols,
                                       nuisance = c("series","SpeciesName"),
                                       n.cores = 4) %>%
  mutate(combo_group = case_when(
    grepl("[.]",threat) ~ "interactive",
    grepl("\\+",threat) ~ "additive",
    TRUE ~ "single")) #posterior timeseries estimated for each threat combination: singular (e.g. "threat1"), interactive (e.g. "threat1.threat2") and additive (e.g. "threat1 + threat2")

post_dydx_intvadd <- do.call("rbind",lapply(all_threats,function(x){
  
  out <- postdraws_intvadd %>%
    subset(grepl(x,threat)) %>% #filter to focal threat
    reframe(.value = mean(diff(.value)/diff(time)),.by = c(combo_group,threat,.draw)) %>% #estimate each timeseries' first derivative
    group_by(threat) %>%
    #filter(!any(.value >= abs(0.5))) %>% #drop highly variable threats
    ungroup() %>%
    mutate(threat_group = x) 
  
  return(out)
})) 

dydx_interval_intvadd <- post_dydx_intvadd  %>%
  #mutate(.chain = 1, .iteration = .draw) %>% #add additional columns required by ggdist
  group_by(threat_group,combo_group,threat) %>%
  ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw"))  #extract distribution information

ggplot(data = post_dydx_intvadd , 
       aes(x = .value,y=threat,fill = combo_group)) +
  tidybayes::stat_slab(alpha=0.5) +
  ggdist::geom_pointinterval(data = dydx_interval_intvadd,
                             aes(xmin = .lower, xmax = .upper,col = combo_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  scale_colour_manual(values = c("#F5A433","#8B33F5","#57A06B"), guide = "none") + 
  scale_fill_manual(values = c("#F5A433","#8B33F5","#57A06B"), guide = "none") + 
  #facet_wrap(~threat_group) + 
  coord_cartesian(xlim = c(-0.5,0.5))+
  xlab("Population trend") + 
  ylab("Threat") + 
  theme_minimal()

post_intvadd_diff <- do.call("rbind",lapply(additive_cols[grepl("\\+",additive_cols)],function(x){
  
  post_dydx_intvadd %>%
    subset(threat %in% c(x,gsub(" \\+ ",".",x))) %>% #subset to shared additive and interactive threats (e.g. "threat1.threat2" and "threat1 + threat2")
    #reframe(.value = .value[2]-.value[1], .by = c(threat_group,.draw)) %>% #find difference in derivatives between additive and interactive threats
    reframe(.value = diff(c(.value[2],.value[1])), .by = c(threat_group,.draw)) %>% #find difference in derivatives between additive and interactive threats
    mutate(threats = gsub(" \\+ ",".",x) ) #name threat combination
  
  # post_dydx_intvadd %>%
  #   subset(threat %in% c(x,gsub(" \\+ ",".",x))) %>% #subset to shared additive and interactive threats (e.g. "threat1.threat2" and "threat1 + threat2")
  #   group_by(threat_group,.draw) %>%
  #   summarise(first_threat = threat[1],
  #          second_threat = threat[2],
  #          .value = .value[2]-.value[1]) %>%
  #   ungroup() %>%
  #   mutate(threats = gsub(" \\+ ",".",x) ) #name threat combination
  
}))

post_interval_intvadd_diff <- post_intvadd_diff %>%
  na.omit() %>% #drop missing threats 
  group_by(threat_group,threats) %>%
  ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw")) %>%  #extract distribution information
  mutate(interaction.type = case_when(
    .upper < 0 ~ "synergistic",
    .lower > 0 ~ "antagonistic",
    TRUE ~ "additive"
  ))

#######
#Generate Figure 3
#######
ggplot(data = na.omit(post_intvadd_diff), 
       aes(x = .value,y=threats)) +
  tidybayes::stat_slab(data = na.omit(post_intvadd_diff) %>%
                         merge(select(post_interval_intvadd_diff,-.value),
                               by = c("threat_group","threats")) %>%
                         subset(.width == 0.8)
                       ,aes(fill = interaction.type),alpha=0.5,normalize = "xy") +
  ggdist::geom_pointinterval(data = post_interval_intvadd_diff,
                             aes(xmin = .lower, xmax = .upper,col = interaction.type),alpha=0.5) +
  geom_point(data = subset(post_interval_intvadd_diff,.width == 0.8), 
             aes(x = .value,col = interaction.type,fill = interaction.type),shape = 21,alpha=0.5,size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  coord_cartesian(xlim = c(-0.25,0.25)) + 
  scale_fill_manual(values = c("#694364",
                               "#1E63B3",
                               "#B32315"), name = "",
                    guide = guide_legend(override.aes = list(color = NA,shape = 2) )) + 
  scale_color_manual(values = c("#694364",
                               "#1E63B3",
                               "#B32315"), guide = "none") + 
  #facet_wrap(~threat_group) + 
  xlab( expression(paste("Additive ",partialdiff,"y","/",partialdiff,"x"," - interactive ",partialdiff,"y","/",partialdiff,"x"))) + 
  ylab("Threat combination") + 
  theme_minimal()

ggsave("Results/Figure3.pdf",last_plot(),
       height = 8,width = 8,dpi = 300)

post_intvadd_diff_syns <- post_intvadd_diff  %>%
  na.omit() %>% #drop missing threats 
  reframe(prop_great_zero = sum(.value>0)/n(),
          prop_less_zero = sum(.value<0)/n(),
          zero_quantile = ecdf(.value)(0),
          .by = c(threat_group,threats)) %>%
  mutate(interaction.type = case_when(
    zero_quantile > 0.9 ~ "synergistic",
    zero_quantile < 0.1 ~ "antagonistic",
    TRUE ~ "additive"
  )) %>%
  group_by(threat_group,interaction.type) %>%
  summarise(n = n()) %>%
  complete(interaction.type) %>%
  mutate(n=replace_na(n, 0),
         freq = (n / sum(n))*100) %>%
  ungroup()

ggplot(post_intvadd_diff_syns) +
  # Add the stacked bar
  geom_bar(aes(x=as.factor(threat_group), y=freq, fill=interaction.type), 
           stat="identity", alpha=0.8) +
  scale_fill_manual(name = "",
                    values = c("#694364",
                               "#B32315",
                               "#1E63B3"))+
  # Add text showing the freq of each 100/75/50/25 lines
  annotate("text", x = rep(max(post_intvadd_diff_syns$threat_group),5), y = c(20,40, 60, 80,100), 
           label = c("20", "40", "60", "80","100"), 
           color="grey", size=5, angle=0, fontface="bold", hjust=1.5) +
  ylim(-100,220) +
  #theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() 
