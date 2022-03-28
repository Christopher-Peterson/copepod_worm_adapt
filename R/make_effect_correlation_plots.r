library(tidyverse)
library(patchwork)
library(ggforce)


theme_set(cowplot::theme_cowplot())

# Load the data ####
full_effect_sizes = read_rds("out/full_effect_sizes.rds")


# set.seed(123)
#### Copepods ####
make_cor_plot = function(parent_factor = "Genus",
                         child_factor = "Copepod Lake", max_value = 3) {S
  parent = rlang::ensym(parent_factor)
  child = rlang::ensym(child_factor)
  # Child is nested within parent
  
  dat_long = full_effect_sizes %>% ungroup() %>% 
    filter(factor %in% c(parent_factor, child_factor)) %>% 
    filter(effect_size <= max_value)
  dat_wide = dat_long %>% filter(!absent) %>% 
    select(-absent) %>% 
    pivot_wider(names_from = factor, values_from = effect_size)
  
  # SPlit wide data into ones with & without the child factor
  dat_wide_no = dat_wide %>% filter(is.na(!!child)) %>% select(-(!!child))
  dat_wide_yes = dat_wide %>% filter(!is.na(!!child)) 
  
  # Create a list of max points to ensure scaling is the same between plots,
  # witin components
  axis_max = dat_long %>% group_by(component) %>% 
    summarise(max = max(effect_size, na.rm = TRUE) ) %>% 
    mutate(!!child := -max/5) 

    # group_split(component) %>% 
    # map(~geom_blank(data = .x))
    # 
  
  # dat_wide %>% mutate() # How to show individual comtributions?
  
  # Make three figures
  # one_factor_plot = dat_wide_no %>% 
  #   ggplot(aes(x = 1, y = !!parent)) + 
  #   geom_sina(alpha = 0.05) +  
  #   # scale_color_brewer(palette = "Set2", guide = 'none') + 
  #   facet_wrap(~component, ncol = 1, strip.position = "left", scales = 'free_y') + 
  #   theme(strip.background = element_blank(),
  #         strip.text = element_blank(),
  #         axis.text.x = element_blank(),
  #         axis.title.x = element_blank(),
  #         axis.ticks.x = element_blank() ) + 
  #   
    # axis_max
  dat_wide_yes %>% 
    ggplot(aes(x = !!child, y = !!parent)) + 
    # geom_point(alpha = 0.05) + 
    stat_density_2d(
      geom = "raster",
      aes(fill = after_stat(ndensity)),
      # contour_var = "ndensity",
      contour = FALSE
    ) + 
    scale_fill_viridis_c() + 
    # scale_color_brewer(palette = "Set2", guide = 'none') + 
    facet_wrap(~component, ncol = 1, strip.position = "right", scale = 'free') + 
    theme(strip.background = element_blank()
          # strip.text = element_blank(),
          # axis.text.y = element_blank(),
          # axis.title.y = element_blank()
          ) + 
      geom_violin(data = dat_wide_no %>% left_join(axis_max)) 
    
    # axis_max
  # browser()
  # two_factor_plot
  # Return:  
  # one_factor_plot + two_factor_plot + plot_layout(widths = c(1, 2))
}
# make_cor_plot()
cop_cor_plot = make_cor_plot(max_value = 1.5)

ggsave("out/copepod_effect_cor.png", cop_cor_plot, dpi = 300, width = 7, height = 5)

worm_cor_plot = make_cor_plot("Worm Lake", "Worm Family (RE)", max_value = 1.5)
ggsave("out/worm_effect_cor.png", worm_cor_plot, dpi = 300, width = 7, height = 5)

## Same, but for worms ####

# worm_cor_plot = local({
#   
#   worm_main = full_effect_sizes %>% ungroup() %>% 
#     filter(factor %in% c("Worm Lake", "Worm Family (RE)"))
#   
#   worm_wide = worm_main %>% filter(!absent) %>% 
#     select(-absent) %>% 
#     pivot_wider(names_from = factor, values_from = effect_size) %>% 
#     rename(Worm_Lake = `Worm Lake`, Worm_Family = `Worm Family (RE)`)
# 
#   worm_wide_no_fam = worm_wide %>% filter(is.na(Worm_Family))
#   worm_wide_yes_fam = worm_wide %>% filter(!is.na(Worm_Family))
#   axis_range = range(worm_main$effect_size, na.rm = TRUE)
#   
#   # Make two figures
#   no_lake_plot = worm_wide_no_fam %>% 
#     ggplot(aes(x = 1, y = Worm_Lake, color = component)) + 
#     geom_sina(height = 0, width = .25, alpha = 0.05) +  
#     scale_color_brewer(palette = "Set2", guide = 'none') + 
#     facet_wrap(~component, ncol = 1, strip.position = "left") + 
#     theme(strip.background = element_blank(), 
#           axis.text.x = element_blank(),
#           axis.title.x = element_blank(),
#           axis.ticks.x = element_blank() ) + 
#     ylim(axis_range[1], axis_range[2])
#   lake_plot = worm_wide_yes_fam %>% 
#     ggplot(aes(x = Worm_Family, y = Worm_Lake, color = component)) + 
#     geom_point(alpha = 0.05) + 
#     scale_color_brewer(palette = "Set2", guide = 'none') + 
#     facet_wrap(~component, ncol = 1, strip.position = "left") + 
#     theme(strip.background = element_blank(), 
#           strip.text = element_blank(),
#           axis.text.y = element_blank(),
#           axis.title.y = element_blank()) + 
#     ylim(axis_range[1], axis_range[2]) + 
#     xlim(axis_range[1], axis_range[2]) + 
#     coord_fixed()
#   # Return:  
#   no_lake_plot + lake_plot + plot_layout(widths = c(1, 2))
# })
