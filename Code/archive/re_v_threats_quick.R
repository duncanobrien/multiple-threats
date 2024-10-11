require(tidyverse)

m1 <- readRDS("Results/models/mod_global_rerun.RDS")

threat_cols <-  colnames(m1$data)[grepl(paste(c("pollution","habitatl",
                                                "climatechange","invasive",
                                                "exploitation","disease"),
                                              collapse = "|"),
                                        colnames(m1$data))]

int_dat = expand.grid(scaled_year = c(0,1),
                      series = unique(m1$data$series)) %>%
  group_by(series) %>%
  mutate(time = seq_along(scaled_year)) %>%
  ungroup() %>%
  left_join(distinct(m1$data %>%
                       select(c(series,threat_cols,Site,SpeciesName))),by = "series") %>%
  mutate(int_present = if_any(threat_cols,~.x == "1",TRUE),
         y_centered = NA) 

null_dat = int_dat %>%
  mutate(across(all_of(threat_cols),
                ~factor("0",levels = c("0","1"))))

base_pred = brms::posterior_epred(m1,
                                  newdata = null_dat,
                                  re.form = NA,
                                  incl_autocor = F,
                                  sort = TRUE,
                                  ndraws = 1000) %>% #extract posterior draws for the above data grid
  colMeans(.) %>%
  cbind(null_dat) %>%
  rename(".value" = ".") %>%
  reframe(.value = diff(.value)/diff(time),.by = c(series)) %>%
  left_join(null_dat %>%
              select(series,int_present) %>%
              distinct(),
            by = "series")

space_pred = brms::posterior_epred(m1,
                                   newdata = null_dat,
                                   re.form = NULL,
                                   incl_autocor = F,
                                   sort = TRUE,
                                   ndraws = 1000) %>% #extract posterior draws for the above data grid
  colMeans(.) %>%
  cbind(null_dat) %>%
  rename(".value" = ".") %>%
  reframe(.value = diff(.value)/diff(time),.by = c(series)) %>%
  left_join(null_dat %>%
              select(series,int_present) %>%
              distinct(),
            by = "series")

inter_pred = brms::posterior_epred(m1,
                                   newdata = int_dat,
                                   re.form = NA,
                                   incl_autocor = F,
                                   sort = TRUE,
                                   ndraws = 1000) %>% #extract posterior draws for the above data grid
  colMeans(.) %>%
  cbind(null_dat) %>%
  rename(".value" = ".") %>%
  reframe(.value = diff(.value)/diff(time),.by = c(series)) %>%
  left_join(null_dat %>%
              select(series,int_present) %>%
              distinct(),
            by = "series")

inter_diff <- rbind(base_pred,inter_pred) %>%
  reframe(beta = diff(.value),.by = c(series,int_present)) %>%
  mutate(type = "Interactions")

space_diff <- rbind(base_pred,space_pred) %>%
  reframe(beta = diff(.value),.by = c(series,int_present))%>%
  mutate(type = "Random effects")

ggplot(rbind(inter_diff,space_diff) %>%
         filter(int_present == TRUE)) +
  geom_jitter(aes(x = beta, y = type), width = 0.01, height = 0.1, alpha = 0.05) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  labs(x = "Difference in trend", y = "") +
  theme_classic()
