# Calculate model-averaged expected value & create plot
#### Setup ####
suppressPackageStartupMessages({
  library(readr)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(brms)
  library(dplyr)
  # library(stringr)
  library(ggplot2)
  library(cowplot); 
})
  theme_set(theme_cowplot())

order_cop_lake = function(cop_lakes) {
  # Convert cop.lake into a factor
  lake_order = list(gos = "Gossling", ech = "Echo",
                    rob = "Roberts", boo = "Boot", 
                    lau = "Lawier") %>% rev()
  rlang::exec(recode, .x = cop_lakes, !!!lake_order) %>%
    factor(levels = lake_order)
}
# raw_data = read_csv("data/chapter2.copepods.cleaned.csv") %>%
#   mutate(cop.lake = order_cop_lake(cop.lake))

raw_data = read_csv("data/chapter_2_copepod_for_bayes.csv") %>%
  mutate(cop.lake = order_cop_lake(cop.lake))

#### Combine stacking weights, get evs
source("R/weights_stacking_combine.R")
#### Read EVs ####
full_evs = read_csv("out/stacked_evs.csv")

#### Create Plotting Data ####
get_ci = function(x, intervals = c(.5, .95)) {
  # browser()
  tibble(
    interval = rep(intervals, 2),
    quantile = c(0.5 - intervals/2, 0.5 + intervals/2),
    direction = c(".lower",".upper") %>%
      rep(each = length(intervals))
  ) %>% mutate(
    .ci = quantile(x, quantile)
  ) %>% select(-quantile)
}
mean_no_zero = function(x) {
  y = x[x>0]
  if(length(y) == 0) y = NA
  mean(y)
}
se_no_zero = function(x) {
  y = x[x>0]
  n = length(y)
  if(n == 0) y = NA
  sd(y)/sqrt(n)
}

pred_smry_3 = full_evs %>% select(-id) %>% 
  pivot_longer(c(num_worms, worm_freq), 
               names_to = ".variable", values_to = "pred") %>% 
  mutate(.variable = recode(
    .variable, num_worms = "Incidence", worm_freq = "Prevalence")) %>% 
  group_by(.variable, cop.lake, worm.lake) %>%
  summarize(.mean = mean(pred),
            .med = median(pred),
            .sd = sd(pred),
            .ci = get_ci(pred)) %>% 
  ungroup() %>% 
  unpack(c(.ci)) %>%
  mutate(cop.lake = order_cop_lake(cop.lake)) %>%
  pivot_wider(names_from = direction, values_from = .ci) 
  

# Observed mean values
obs_means = raw_data %>% group_by(worm.lake, cop.lake) %>%
  select(numb.worm) %>%
  summarize(Prevalence = mean(numb.worm > 0), prev_se = sd(numb.worm > 0)/sqrt(n()) ,
            Incidence = mean_no_zero(numb.worm),
            inc_se = se_no_zero(numb.worm)) %>%
  ungroup() #%>%
# mutate(cop.lake = order_cop_lake(cop.lake))
obs_means_tidy = obs_means %>%
  gather(key = ".variable", value = ".mean", Prevalence, Incidence)

#### Create the Plots ####
viridis_option = 'D' # C is good so far, though maybe with a bit off the end
# D is better, though the middle color could be slightly shifted?
# So I'm almost surely calculating the CI's improperly
count_plot = raw_data %>% 
    filter(numb.worm>0) %>% 
  ggplot(aes(x = cop.lake, y = numb.worm, group = worm.lake)) + 
  geom_point(aes(group = worm.lake, color = worm.lake, fill = worm.lake),
             # fill aes() is for position
             # color = grey(.3),
             alpha = .4, shape = 1,
             position = position_jitterdodge(
               jitter.width = .05, jitter.height = .3,
               dodge.width = .4, seed = 12345
             )) + 
  geom_linerange(data = pred_smry_3 %>% filter(.variable == "Incidence"),
                 position = position_dodge(width = .4),
                 mapping = aes(y = .med,group = worm.lake,
                               ymin = .lower, 
                               ymax = .upper,
                               size = as.factor(interval)),
                  color = "black") +
  
  geom_point(data = pred_smry_3 %>% filter(.variable == "Incidence"),
             position = position_dodge(width = .4),
             color = "black", size = 3,
             shape = 4,
             mapping = aes(y = .med, group = worm.lake)) + 
  
  geom_point(data = obs_means,
             position = position_dodge(width = .4), 
             mapping = aes(y = Incidence, color = worm.lake),
             shape = "—", size = 5) + 
  theme(legend.position = c(.28, .9), legend.justification = c(1,1))+
  scale_color_viridis_d(name = "Tapeworm\nLake", 
                        label = c("Boot", "Echo", "Gossling"),
                        aesthetics = c("colour", "fill")) + 
  
  scale_x_discrete(name = "Copepod Lake") +
  scale_size_manual(guide = "none", values = c(1.2, .75)) + 
  # label = c("Boot","Echo","Gossling","Lawier", "Roberts"))+
  scale_y_continuous(name = "Number of Tapeworms\nin infected Copepods")
# count_plot

freq_plot = raw_data %>% 
  group_by(cop.lake, worm.lake) %>% 
  summarize(worm_freq = mean(numb.worm > 0)) %>% 
  ggplot(aes(x = cop.lake, y = worm_freq, group = worm.lake)) + 
  ylim(0, 1) + 
  geom_col(aes(fill = worm.lake, color = worm.lake),
           alpha = 0.2, width = 0.35,
           position = position_dodge(width = 0.4)) + 
  scale_color_viridis_d(guide = "none", option = viridis_option,
                        aesthetics = c("colour", "fill")) + 
  # scale_fill_viridis_d(guide = "none", option = viridis_option) + 
  geom_linerange(data = pred_smry_3 %>% 
                   filter(.variable == "Prevalence"),
                 position = position_dodge(width = .4),
                 mapping = aes(y = .med,group = worm.lake,
                               ymin = .lower, 
                               ymax = .upper,
                               size = as.factor(interval)),
                 color = "black") +
  
  geom_point(data = pred_smry_3 %>% 
               filter(.variable == "Prevalence"),
             position = position_dodge(width = .4),
             color = "black", size = 4,
             mapping = aes(y = .med, group = worm.lake)) + 
  ylab("Infection Frequency") + 
  scale_size_manual(guide = "none", values = c(1.2, .75)) + 
  xlab("Copepod Lake")

joint_plot = plot_grid(
  count_plot + theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank()),
  freq_plot, align = "v", axis = "lr", 
  labels = c("A", "B"),
  label_x = .13, label_y = .97,
  ncol = 1, rel_heights = c(1, .6)
  ); ggsave("figures/full_plot.png",joint_plot, 
            width = 6, height = 8, dpi = 300)
cb_plot = colorblindr::cvd_grid(count_plot + scale_y_continuous(name=""))
ggsave("figures/colorblind_plot.png", cb_plot, width = 11, height = 7, dpi = 300)

###### Old; skip  ####


# Calculate expected values for each iteration
# Broken down by overall cop*worm means
# As well as cop*worm by worm_fam and plate means
# mns = list(
#   means = tidy_preds %>%  map(. %>% unnest(cols = c(data)) %>%
#       group_by(cop.lake, worm.lake, .iter) %>%
#         summarize(mean_pred = mean(pred))),
#   fams = tidy_preds %>% map(. %>% unnest(cols = c(data)) %>%
#       group_by(cop.lake, worm.lake, .iter, worm.fam) %>%
#       summarize(fam_pred = mean(pred))),
#   plates = tidy_preds %>% map(. %>% unnest(cols = c(data)) %>%
#       group_by(cop.lake, worm.lake, .iter, plate) %>%
#       summarize(plate_pred = mean(pred)))) %>% transpose()

# Take the previously defined means
# and average over family and plate
# This is still per iteration
# sum_means = mns %>% map(. %>% map(ungroup)) %>%
#   map(~do.call(function(means, fams, plates){
#   # browser()
#   left_join(means,
#       fams, by = c("cop.lake", "worm.lake", ".iter")) %>%
#       left_join(plates, by = c("cop.lake", "worm.lake", ".iter")) %>%
#     group_by(cop.lake, worm.lake, .iter) %>%
#     summarize(mean = unique(mean_pred),
#               sd_fam = sd(fam_pred - mean),
#               sd_plt = sd(plate_pred - mean))
# }, .x )) %>% imap_dfr(~mutate(ungroup(.x), .variable = .y))
# pred_mean_sd = sum_means %>%
#   group_by(.variable, cop.lake, worm.lake) %>%
#   summarize(sd = sd(mean), mean = mean(mean), median = median(mean),
#             fam = median(sd_fam), plt = median(sd_plt))

# Not entirely sure how this is different from sum_means in principle



# Plot mean & standard error
# Should I replace this with 95% Credible Intervals?
# mean_with_se_plots = pred_smry_2  %>%
#   ggplot(aes(x = cop.lake, y = .mean, color = worm.lake)) +
#   facet_wrap(. ~.variable, scales = "free_y", strip.position = "left") +
#   geom_linerange(aes(ymin = .mean - .sd, ymax = .mean+.sd),
#                  position = position_dodge(width = .2)) +
#   geom_point(position = position_dodge(width = .2)) +
#   scale_color_viridis_d(name = "Tapeworm\nLake", #option = "magma", end = .8,
#                         label = c("Boot", "Echo", "Gossling")) +
#   scale_x_discrete(name = "Copepod Lake") +
#                    # label = c("Boot","Echo","Gossling","Lawier", "Roberts"))+
#   # scale_y_continuous(name = "")+
#   theme(strip.placement = "outside",
#         strip.background = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(angle = 30, hjust = 1),
#         legend.position = "top")
# 
# mean_with_se_obs = mean_with_se_plots +
#   geom_point(data = obs_means_tidy, shape = 4, position = position_dodge(width = .2))
# 
# ggsave("figures/mean_with_se.png", mean_with_se_obs,width = 5, height = 8, dpi = 300)
# 
# # TO DO:
# # add model fits to this (with CI?)
# # make raw-data x's
# # Add proportions
# # split raw data into silos?
# 
# rawDatPlot=
# raw_data %>% filter(numb.worm>0) %>%
#   ggplot(aes(x = cop.lake, y = numb.worm)) +
#   geom_jitter(width = .3, height = .3, color = grey(.7), alpha = .5) +
#   geom_linerange(data = obs_means,
#                  position = position_dodge(width = .35),
#                 mapping = aes(y = Incidence, color = worm.lake,
#                               ymin = Incidence - inc_se*1.96,
#                               ymax = Incidence + inc_se*1.96))+
#   geom_point(data = obs_means, position = position_dodge(width = .35),
#              mapping = aes(y = Incidence, color = worm.lake)) +
#   scale_color_viridis_d(name = "Tapeworm\nLake", label = c("Boot", "Echo", "Gossling")) +
#   scale_x_discrete(name = "Copepod Lake") +
#                    # label = c("Boot","Echo","Gossling","Lawier", "Roberts"))+
#   scale_y_continuous(name = "Number of Tapeworms in infected Copepods")
# ggsave("figures/observed_means.png", width = 6, height = 6, dpi = 300)
# # Okay, how to deal with that?
# 
# # Next, average over plate (nobody gives AF)

