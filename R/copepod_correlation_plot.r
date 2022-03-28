# N_DRAWS = draw_details %>% unchop(.draw) %>% nrow()
# Get the credible intervals of effect sizes & the frequency of the factor's inclusion in models

# copepod_genus_effects = local({
  # Separate out genus effect size in presence + absence of copepod lake
  cop_effects = full_effect_sizes %>% 
    filter(factor %in% c("Copepod Lake", "Genus")) %>% 
    ungroup() %>% filter(!absent)
  lake_effects = cop_effects %>% filter(factor == "Copepod Lake")
  genus_effects = cop_effects %>% filter(factor == "Genus")
  
  genus_with_lake = genus_effects %>% 
    semi_join(lake_effects, by = c("component", ".id")) %>% 
    mutate(partition = "Lake")
  genus_without_lake = genus_effects %>% 
    anti_join(lake_effects, by = c("component", ".id")) %>% 
    mutate(partition = "No Lake")
  
  split_genus_effects = bind_rows(genus_with_lake, genus_without_lake) %>% 
    select(-absent, -factor) %>% 
    group_by(component, partition) %>% 
    select(-.id) %>% 
    chop(effect_size) %>% 
    mutate(ci = map(effect_size, get_ci)) %>% 
    select(-effect_size) %>% unnest(ci) %>% 
    ungroup()
  # split_genus_effects %>% 
  #   select(component, partition, median = Q.5, lo = Q.025, hi = Q.975) %>% 
  #   mutate(across(median:hi, ~format(round(.x,3), digits = 3))) %>% 
  #   mutate(CI = glue("[{lo}, {hi}]")) %>% 
  #   select(-lo, -hi) %>% arrange(component, partition)
  
  

  
# })
  
joint_lake_genus = lake_effects %>% 
  select(component, .id, lake = effect_size) %>% 
  left_join(genus_effects %>% 
              select(component, .id, genus = effect_size),
            by = c("component", ".id"))
joint_lake_genus %>% ggplot(aes(x = genus, y = lake)) + 
  geom_density_2d_filled() + 
  scale_x_log10() + scale_y_log10()
# library(ggdist)

cop_eff_cats = cop_effects %>% 
  mutate(class = case_when(
    factor == "Copepod Lake" ~"Copepod Lake",
    .id %in% joint_lake_genus$.id~ "Genus with Lake",
    TRUE ~ "Genus without Lake"
  ))
cop_eff_cats %>% 
  filter(effect_size < 1.5) %>% 
  # filter(effect_size < 2.72) %>% 
  ggplot() + aes(x = class, y = effect_size) + 
  geom_violin() + theme_classic() + 
  facet_wrap(~component, scale = 'free_y')


effect_cutoff = 1.5 #2.714

joint_lake_data_counts = joint_lake_genus %>% 
  bind_rows(genus_without_lake %>% 
              select(component, .id, genus = effect_size) %>% mutate(lake = 0)) %>% 
  filter(component == "count") %>% 
  filter(lake < effect_cutoff, genus < effect_cutoff)

copepod_2d_dens_count =   joint_lake_data_counts%>% 
  # Density estimate
  with(MASS::kde2d(genus, lake, n = 400)) %>% 
  # convert to data frame
  with({
    n = length(x)
    tibble(x = rep(x, times = n),
           y = rep(y, each = n),
           z = as.numeric(z))
    }) %>% 
  mutate(alpha = if_else(z < 1, z^2, 1))


joint_lake_data_counts
library(ggdist)

marginal_copepod_plot_counts_genus = bind_rows(genus_with_lake, genus_without_lake) %>% 
  filter(component == "count") %>% 
  mutate(side = recode(partition, `No Lake` = "top", Lake = "bottom"),
         partition = recode(partition, `No Lake` = "No Lake\nEffect", Lake = "With Lake\nEffect"),) %>% 
  ggplot(aes(x = effect_size, y = partition)) + 
  ggdist::stat_halfeye(aes(side = side), scale = .8, n = 10000) + 
  # ggdist::stat_slabinterval()
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.title.y = element_blank()) +
  # xlim(0, effect_cutoff) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0.05, 0.0)) + 
  annotate("text", x = 1.5, y = 1, label = "Without Lake Effect", hjust = 1, vjust = -0.2, size = 5) +
  annotate("text", x = 1.5, y = 2, label = "With Lake Effect", hjust = 1, vjust = 1.3, size = 5) +
  coord_fixed(xlim = c(0, effect_cutoff), ratio = 1/3) +
  xlab("Copepod Genus (Effect Size)")
marginal_copepod_plot_counts_lake = lake_effects %>%  
  filter(component == "count") %>% 
  ggplot(aes(y = effect_size)) + 
  ggdist::stat_halfeye(side = "left", scale = .9, n = 10000) + 
  # ggdist::stat_slabinterval()
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  # xlim(0, effect_cutoff) + 
  scale_x_discrete(expand = c(0.01, 0.03)) +
  coord_fixed(ylim = c(0, effect_cutoff), ratio = 3) +
  ylab("Copepod Lake (Effect Size)")


counts_2d_plot = copepod_2d_dens_count %>% 
  mutate(alpha = pmin(z, 1)) %>% 
  ggplot(aes(x, y, fill = z,  alpha = alpha)) + 
  geom_raster() + 
  scale_fill_viridis_c("Posterior\nDensity", trans =  "log1p",
                       breaks = c(2, 4, 8, 16, 30),
                       rescaler =~pmax(.x, 1) %>% scales::rescale(),
                       direction = -1, ) + 
  scale_alpha_continuous(range = c(0, 1),
                         # rescaler = ~if_else(.x < 1, .x^2, 1),
                         guide = "none") + 
  cowplot::theme_cowplot() + 
  # scale_fill_viridis_c(direction = -1) + 
  coord_fixed(expand = FALSE, ylim = c(0, effect_cutoff), xlim = c(0, effect_cutoff)) + 
  theme(legend.position = c(1, 1), legend.justification = c(1,1),
        axis.title = element_blank(), 
        axis.text = element_blank())

cop_cor_plot = marginal_copepod_plot_counts_lake + counts_2d_plot +
  plot_spacer() + marginal_copepod_plot_counts_genus + 
  plot_layout(ncol = 2)

ggsave("figures/copepod_effect_posterior_plot.png", cop_cor_plot, width = 8, height = 8, dpi = 300)  

#### Repeat it all, but for prop ####

joint_lake_data_prop = joint_lake_genus %>% 
  bind_rows(genus_without_lake %>% 
              select(component, .id, genus = effect_size) %>% mutate(lake = 0)) %>% 
  filter(component != "count") %>% 
  filter(lake < effect_cutoff, genus < effect_cutoff)

copepod_2d_dens_prop =   joint_lake_data_prop%>% 
  # Density estimate
  with(MASS::kde2d(genus, lake, n = 400)) %>% 
  # convert to data frame
  with({
    n = length(x)
    tibble(x = rep(x, times = n),
           y = rep(y, each = n),
           z = as.numeric(z))
  }) %>% 
  mutate(alpha = if_else(z < 1, z, 1))


joint_lake_data_counts
library(ggdist)
effect_cutoff_prop = 0.6
dens_breaks_prop = c(2, 8, 32, 128)
marginal_copepod_plot_prop_genus = bind_rows(genus_with_lake, genus_without_lake) %>% 
  filter(component != "count") %>% 
  mutate(side = recode(partition, `No Lake` = "top", Lake = "bottom"),
         partition = recode(partition, `No Lake` = "No Lake\nEffect", Lake = "With Lake\nEffect"),) %>% 
  ggplot(aes(x = effect_size, y = partition)) + 
  ggdist::stat_halfeye(aes(side = side), scale = .8, n = 10000) + 
  # ggdist::stat_slabinterval()
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.title.y = element_blank()) +
  # xlim(0, effect_cutoff) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0.05, 0.0)) + 
  annotate("text", x = effect_cutoff_prop, y = 1, label = "Without Lake Effect", hjust = 1, vjust = -0.2, size = 5) +
  annotate("text", x = effect_cutoff_prop, y = 2, label = "With Lake Effect", hjust = 1, vjust = 1.3, size = 5) +
  coord_fixed(xlim = c(0, effect_cutoff_prop), ratio = 1/7.5) +
  xlab("Copepod Genus (Effect Size)")
marginal_copepod_plot_prop_lake = lake_effects %>%  
  filter(component != "count") %>% 
  ggplot(aes(y = effect_size)) + 
  ggdist::stat_halfeye(side = "left", scale = .9, n = 10000) + 
  # ggdist::stat_slabinterval()
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  # xlim(0, effect_cutoff) + 
  scale_x_discrete(expand = c(0.01, 0.03)) +
  coord_fixed(ylim = c(0, effect_cutoff_prop), ratio = 7.5) +
  ylab("Copepod Lake (Effect Size)")


prop_2d_plot = copepod_2d_dens_prop %>% 
  mutate(alpha = pmin(z, 1)) %>% 
  ggplot(aes(x, y, fill = z,  alpha = alpha)) + 
  geom_raster() + 
  scale_fill_viridis_c("Posterior\nDensity", trans =  "log1p",
                       breaks = dens_breaks_prop,
                       rescaler =~pmax(.x, 1) %>% scales::rescale(),
                       direction = -1, ) + 
  scale_alpha_continuous(range = c(0, 1),
                         # rescaler = ~if_else(.x < 1, .x^2, 1),
                         guide = "none") + 
  cowplot::theme_cowplot() + 
  # scale_fill_viridis_c(direction = -1) + 
  coord_fixed(expand = FALSE, ylim = c(0, effect_cutoff_prop), xlim = c(0, effect_cutoff_prop)) + 
  theme(legend.position = c(1, 1), legend.justification = c(1,1),
        axis.title = element_blank(), 
        axis.text = element_blank())

cop_cor_plot_prop = marginal_copepod_plot_prop_lake + prop_2d_plot +
  plot_spacer() + marginal_copepod_plot_prop_genus + 
  plot_layout(ncol = 2)

ggsave("figures/copepod_effect_prop_posterior_plot.png", cop_cor_plot_prop, width = 8, height = 8, dpi = 300)  


combined_cor_fig = 
  marginal_copepod_plot_counts_lake + counts_2d_plot +
  plot_spacer() + marginal_copepod_plot_counts_genus + 
  marginal_copepod_plot_prop_lake + prop_2d_plot +
  plot_spacer() + marginal_copepod_plot_prop_genus + 
  plot_layout(nrow = 4) + 
  plot_annotation(tag_levels = list(c("A", "", "",  "B", "",  "")))



# combined_cor_fig = cop_cor_plot / cop_cor_plot_prop + plot_layout(nrow = 2) + 
#   plot_annotation(tag_levels = list(c("A", "", "", ""), c("B", ""< "")) ) 

ggsave("figures/copepod_effect_posterior_combined.png", combined_cor_fig, width = 8, height = 16, dpi = 300)  
