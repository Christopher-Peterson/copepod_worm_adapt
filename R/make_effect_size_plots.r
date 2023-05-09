library(tidyverse)
effects = read_csv("out/marginal_effect_sizes.csv") %>% 
  filter(!factor %in% c("Worm Lake + Family", "Copepod Genus + Lake")) %>% # Something is bugged with that one
  arrange(rev(component), frequency, rev(Q.5)) %>% 
  mutate(factor_fct = fct_inorder(factor)) %>% 
  mutate(factor_fct = recode(factor_fct, Genus = "Copepod Genus")) %>% 
  mutate(component = if_else(component == "count",
                             "Number of Tapeworms", "Infection Frequency")) %>% 
  select(component, factor_fct, frequency, effect_size = Q.5, lower = Q.025, upper = Q.975)

effects_stacked = effects %>% 
  select(factor_fct, component, x = frequency) %>% 
  mutate(measure = "Model Inclusion Frequency") %>% 
  bind_rows(effects %>% select(-frequency) %>% rename(x = effect_size) %>% 
              mutate(measure = "Effect size [with 95% CI]")) %>% 
  mutate(measure = fct_inorder(measure))

theme_set(theme_classic())
library(patchwork)

make_plot = function(Component) {
  effects_stacked %>% filter(component == Component) %>% 
    ggplot(aes(x = x, y = factor_fct)) + 
    geom_linerange(aes(xmin = lower, xmax = upper), color = grey(.4)) +
    geom_point() + 
    facet_wrap(~measure, nrow = 1,  strip.position = "bottom", scales = 'free_x') + 
    theme(strip.placement = 'outside',
          strip.background = element_blank(), 
          axis.title.x = element_blank(),
          strip.text = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 13))+
    ylab(Component)
}

effect_plot = 
  (make_plot("Infection Frequency") + theme(strip.text = element_blank())) /
  make_plot("Number of Tapeworms")

ggsave("figures/effect_size_plot.png", width = 11, height = 8, dpi = 300)
