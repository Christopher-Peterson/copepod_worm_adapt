library(tidyverse)
library(glue)
settings = read_rds("model_settings.rds")

# Break down and figure the following subets out:

model_text_terms = 
  settings %>% select(model_txt, name) %>% #head() %>% 
  unnest_wider(model_txt) %>% 
  mutate(across(pois:hu, ~str_remove(.x, pattern = "^.+ ~ ") %>% 
                      str_split(fixed(" + ")))) 

# So what are the obvious problems

bad_combo = function(data, higher_terms, required_terms) {
  # Returns TRUE for combos where higher terms are present while
  # any of the required terms are missing
  map_lgl(data, ~any(higher_terms %in% .x) && (!all(required_terms %in% .x)) )
}

# First, remove unnested effects
valid_nesting = model_text_terms %>% mutate(
  across(c(pois, hu), list(
  nest_cop = ~bad_combo(.x, "cop.lake", "genus"),
  nest_worm = ~bad_combo(.x, "(1|worm.fam)", "worm.lake")
  ))) %>% 
  filter(!pois_nest_cop, !pois_nest_worm, !hu_nest_cop, !hu_nest_worm) 

# Get rid of the bad interactions
full_subset = valid_nesting %>% mutate(
  across(c(pois, hu), list(
    gw = ~bad_combo(.x, "genus:worm.lake", c("genus", "worm.lake")),
    wc = ~bad_combo(.x,"worm.lake.cop.lake", c("cop.lake", "worm.lake"))
  ))) %>% 
  filter(!pois_gw, !pois_wc, !hu_gw, !hu_wc) %>% 
  select(name, pois, hu)

keep_combos = full_subset$name
filtered_settings = settings %>% filter(name %in% keep_combos)
write_rds(filtered_settings, "model_settings_filtered.rds")
