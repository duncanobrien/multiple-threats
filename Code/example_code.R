library(brms)
library(ggplot2)

mod = readRDS(file.choose())

threat_cols <-  colnames(mod$data)[grepl(paste(c("pollution","habitatl",
                                                 "climatechange","invasive",
                                                 "exploitation","disease"),
                                               collapse = "|"),
                                         colnames(mod$data))]

dat = mod$data
point1 = subset(dat, scaled_year == 0)
point1$time = 1
point2 = point1
point2$time = 2
point2$scaled_year = 1
base = rbind(point1, point2)
base$y_centered = NA
base$SpeciesName = NA
base$Site = NA
base$series = NA
base[c(9:39)] = 0
base_pred = brms::posterior_epred(mod,
                                  newdata = base,
                                  re.form = NA,
                                  incl_autocor = F,
                                  sort = TRUE,
                                  ndraws = 1000)
base$pred_y = exp(colMeans(base_pred))

space = rbind(point1, point2)
space$y_centered = NA
space[c(9:39)] = 0
space_pred = brms::posterior_epred(mod,
                                   newdata = space,
                                   re.form = NULL,
                                   incl_autocor = F,
                                   sort = TRUE,
                                   ndraws = 1000)
base$pred_yspace = exp(colMeans(space_pred))

inter = rbind(point1, point2)
inter$y_centered = NA
inter$SpeciesName = NA
inter$Site = NA
inter$series = NA
inter_pred = brms::posterior_epred(mod,
                                   newdata = inter,
                                   re.form = NA,
                                   incl_autocor = F,
                                   sort = TRUE,
                                   ndraws = 1000)
base$pred_yinter = exp(colMeans(inter_pred))
r = "1"
base$int_present = apply(inter[,threat_cols], 1, function(r) any(r %in% "1"))



sum_df = data.frame(
  base1 = base$pred_y[c(1:910)],
  base2 = base$pred_y[c(911:1820)],
  space1 = base$pred_yspace[c(1:910)],
  space2 = base$pred_yspace[c(911:1820)],
  inter1 = base$pred_yinter[c(1:910)],
  inter2 = base$pred_yinter[c(911:1820)],
  int_present = base$int_present[c(1:910)]
)
sum_df = subset(sum_df, int_present == TRUE)

sum_df$base = sum_df$base2/sum_df$base1
sum_df$space = sum_df$space2/sum_df$space1
sum_df$inter = sum_df$inter2/sum_df$inter1

sum_df$diff1 = sum_df$space - sum_df$base
sum_df$diff2 = sum_df$inter - sum_df$base

t_df = data.frame(
  type = c(rep("Random effects", nrow(sum_df)), rep("Interactions", nrow(sum_df))),
  beta = c(sum_df$diff1, sum_df$diff2)
)
ggplot(t_df) +
  geom_jitter(aes(x = beta, y = type), width = 0.01, height = 0.1, alpha = 0.05) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  labs(x = "Difference in trend", y = "") +
  theme_classic()