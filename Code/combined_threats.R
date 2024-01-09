
all_threats = c("none","pollution","habitatl","climatechange","invasive", "exploitation","disease")

threatcols <- colnames(mod_intvadd3$data)[grepl(paste(all_threats,collapse = "|"),colnames(mod_intvadd3$data))] 

postdraws <- threat_post_draws(model = mod_intvadd3,
                                      threat_comns = c("none",threatcols)) #estimate posterior draws for all threat singles and combinations
post_dydx <- do.call("rbind",lapply(all_threats,function(x){

  # out <- postdraws %>%
  #   subset(grepl(x,threat)) %>%
  #   reframe(.value = mean(diff(.value)/diff(time)),.by = c(threat,.draw)) %>%
  #   mutate(int_group = ifelse(grepl("\\.",threat),"combined","single"),
  #          threat_group = x)
  
  out <- postdraws %>%
    subset(grepl(x,threat)) %>% #subset to focal threat
    reframe(.value = mean(diff(.value)/diff(time)),.by = c(threat,.draw)) %>% #for each posterior timeseries, estimate the first derivative 
    ungroup() %>%
    mutate(int_group = ifelse(grepl("\\.",threat),"combined","single"))%>% #create grouping column for whether the threat is singular or a combination
    mutate(threat_group = x) %>% #overall threat grouping
    group_by(threat) %>%
    #filter(!any(.value >= abs(0.5))) %>% #drop threats with highly uncertain confidence intervals. Typically certain three way interactions
    ungroup() 
  # %>%
  #   subset(!grepl("[[:lower:]]*\\.[[:lower:]]*\\.",threat)) 

  return(out)
  })) %>%
  subset(threat != "pollution.climatechange")

dydx_interval <- post_dydx  %>%
  #subset(!grepl("[[:lower:]]*\\.[[:lower:]]*\\.",threat)) %>%
  group_by(threat_group,int_group) %>%
  ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw", "threat")) #extract distribution information

threat_palette<-c(MetBrewer::met.brewer(name="Hokusai1", n=6, type="continuous"))

palatte <- data.frame(threat_group = unique(subset(post_dydx,threat != "none")$threat_group),
                      fill_col = threat_palette)

plot_dydx_threats <- post_dydx %>%
  subset(threat != "none") %>%
  left_join(palatte,by = "threat_group") %>%
  mutate(fill_col = factor(fill_col))

p1 <- ggplot(data = plot_dydx_threats, 
       aes(x = .value,y=int_group)) +
  tidybayes::stat_slab(data = subset(plot_dydx_threats,int_group == "single"),
                       aes(fill = fill_col,group = threat), alpha=0.3) +
  tidybayes::stat_slab(data = subset(plot_dydx_threats,int_group == "combined"),
                       aes(fill = fill_col,group = threat,col = fill_col), alpha=0.2) +
  ggdist::geom_pointinterval(data = subset(dydx_interval,threat_group != "none"),
                             aes(xmin = .lower, xmax = .upper),position = position_dodge()) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  xlab("Population trend") + 
  ylab("Threat") + 
  coord_cartesian(xlim = c(-0.2,0.2)) +
  facet_wrap(~threat_group) +
  scale_fill_manual(values = levels(plot_dydx_threats$fill_col),guide = "none") + 
  scale_color_manual(values = levels(plot_dydx_threats$fill_col),guide = "none") + 
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

post_dydx_none <- postdraws %>%
    subset(grepl("none",threat)) %>% #subset to focal threat
    reframe(.value = mean(diff(.value)/diff(time)),.by = c(threat,.draw)) %>% #for each posterior timeseries, estimate the first derivative 
    mutate(int_group = ifelse(grepl("\\.",threat),"combined","single")) #create grouping column for whether the threat is singular or a combination

dydx_none_interval <- post_dydx_none  %>%
  group_by(int_group) %>%
  ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw", "threat")) #extract distribution information

p2 <- ggplot(data = post_dydx_none, 
             aes(x = .value,y=int_group)) +
  tidybayes::stat_slab(aes(), alpha=0.3) +
  ggdist::geom_pointinterval(data = dydx_none_interval,
                             aes(xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  xlab("Population trend") + 
  ylab("Threat") + 
  coord_cartesian(xlim = c(-0.2,0.2)) +
  facet_wrap(~threat) +
  theme_minimal()+
  theme(axis.title.x = element_text(size=12,
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size = 12),
        strip.text.x = element_text(size = 12),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks = element_line(color="black"),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0,0,0,0), "pt"))

require(patchwork)

p1 + (p2/plot_spacer() + plot_layout(heights  = c(1,1))) + plot_layout(widths = c(2,1.5)) 

dydx_interval_full <- post_dydx  %>%
  #subset(!grepl("[[:lower:]]*\\.[[:lower:]]*\\.",threat)) %>%
  group_by(threat) %>%
  ggdist::median_qi(.width = c(.95, .8, .5),.exclude = c(".draw", "threat_group","int_group")) %>% #extract distribution information
  rowwise() %>%
  mutate(n.threats =length(unlist(strsplit(threat,"[.]"))))

ggplot(data = post_dydx %>%
         rowwise() %>%
         mutate(n.threats =length(unlist(strsplit(threat,"[.]")))), 
       aes(x = .value,y=threat)) +
  tidybayes::stat_slab(alpha=0.5) +
  ggdist::geom_pointinterval(data = dydx_interval_full,
                             aes(xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  xlab("Population trend") + 
  ylab("Threat") + 
  coord_cartesian(xlim = c(-0.2,0.2)) +
  facet_wrap(~n.threats,scales = "free_y") +
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
