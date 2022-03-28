## Calculate effect sizes for the Bayesian model
# We'll be defining effect sizes as the standard deviation of the 
  # marginal effects for a factor; it's calculated at each posterior
  # sample, so a posterior distribution can be determined
# This is based on Gelman's suggested effect size for ANOVAs,
  # but we're using marginal effects instead of raw betas
  # because they're more interpretable for non-Gaussian data
library(tidyverse)
library(glue)

## Import data ####

# There's some trouble with readr's automatic column type identifier,
# because there's so many rows of data

# read the first few rows

stacked_draws = read_rds( "out/combined_draws_df.rds")

# stacked_draws_col = read_csv("out/stacked_draws.csv", n_max = 3)
# 
# # Define input column specification, replacing NA's with logicals
# stacked_colspec = map_chr(stacked_draws_col, ~{
#     type = typeof(.x)
#     if(type == "logical") return("d")
#     substr(type, 1, 1)
#   }) %>% 
#   paste(collapse = "")
# # Read the full data
# stacked_draws = read_csv("out/stacked_draws.csv",col_types = stacked_colspec)

# Utility function; Names factors based on parm names
name_factors = function(x) {
  case_when(
    str_detect(x,"r_worm.fam") ~ "Worm Family (RE)",
    str_detect(x,"r_plate") ~  "Plate (RE)",
    str_detect(x, fixed(":worm.lake")) ~ "Copepod x Worm Lake Interaction",
    str_detect(x, fixed(":genus")) ~ "Genus x Worm Lake Interaction",
    str_detect(x, fixed("genus")) ~ "Genus",
    str_detect(x,"cop.lake") ~ "Copepod Lake",
    str_detect(x,"worm.lake") ~ "Worm Lake",
    str_detect(x,"native") ~    "Native")
}
factor_table = tibble(
  parm = names(stacked_draws) %>% str_remove("b_hu_") %>% 
    str_remove("b_") ) %>% 
  mutate( factor = name_factors(parm) ) %>%
  filter(!is.na(factor))

#####################################################
# split data up by model component 
  # (hu = proportion model, no hu = count),
# separate the random effects from the fixed effects;
  # (they need to be processed separately)
# Then re-combine into the full dataset
#####################################################
stacked_draws_relevant = # remove random effects & other unnecessary things
  stacked_draws %>% select(-starts_with("prior_"), 
                           -ends_with("__"), 
                           -contains("r_worm"),
                           -contains("r_plate")) 

stacked_draws_prop = stacked_draws_relevant %>% 
  select(model, contains("_hu"), .chain, .iteration, .draw)
stacked_draws_count = stacked_draws_relevant %>% 
  select(-contains("_hu"))

# These next lines pivot the data into a long form, 
# calculate marginal effects, 
  # which are basically g^-1( beta + Intercept )
# Then standardizes some formatting

long_draws_count = stacked_draws_count %>% 
  select(-model, -starts_with("sd_"), -starts_with(".")) %>% 
  rename(Intercept = b_Intercept) %>% 
  mutate(.id = 1:n()) %>% 
  pivot_longer(-c(.id, Intercept),
               names_to = "parm",
               names_prefix = "b_",
               values_to = "value") %>% 
  mutate(absent = is.na(value)) %>% 
  mutate(linear_value = exp(Intercept + value)) %>% 
  select(-value) %>% 
  left_join(factor_table, by = "parm") %>% 
  mutate(component = "count") 

long_draws_prop = stacked_draws_prop %>% 
  select(-model, -starts_with("sd_"), -starts_with(".")) %>% 
  rename(Intercept = b_hu_Intercept) %>% 
  mutate(.id = 1:n()) %>% 
  pivot_longer(-c(.id, Intercept),
               names_to = "parm",
               names_prefix = "b_hu_",
               values_to = "value") %>% 
  mutate(absent = is.na(value)) %>% 
  mutate(linear_value = gtools::inv.logit(Intercept + value)) %>% 
  select(-value) %>% 
  left_join(factor_table, by = "parm") %>% 
  mutate(component = "prop") 

long_draws_count_r = stacked_draws %>% 
  select(Intercept = b_Intercept, 
         starts_with("r_"), -contains("_hu")) %>% 
  mutate(.id = 1:n()) %>% 
  pivot_longer(-c(.id,Intercept),
               names_to = "parm",
               values_to = "value") %>% 
  mutate(absent = is.na(value)) %>% 
  mutate(linear_value = exp(Intercept + value)) %>% 
  select(-value) %>% 
  left_join(factor_table, by = "parm") %>% 
  mutate(component = "count") 


long_draws_prop_r = stacked_draws %>% 
  select(Intercept = b_hu_Intercept, 
         starts_with("r_") & contains("_hu") ) %>% 
  mutate(.id = 1:n()) %>% 
  pivot_longer(-c(.id,Intercept),
               names_to = "parm",
               values_to = "value") %>% 
  mutate(absent = is.na(value)) %>% 
  mutate(linear_value = gtools::inv.logit(Intercept + value)) %>% 
  select(-value) %>%
  left_join(factor_table, by = "parm") %>% 
  mutate(component = "prop") 
long_draws_list = lst(long_draws_prop_r, long_draws_count_r,
    long_draws_prop, long_draws_count) 
long_draws_list %>% write_rds("out/long_draws_tmp.rds")
# long_draws_list = read_rds("out/long_draws_tmp.rds")
# Pull all of the marginals together
# long_draws_list$long_draws_prop = long_draws_list$long_draws_prop %>%
#   select(-factor) %>%
#   left_join(factor_table, by = "parm")
# long_draws_list$long_draws_count = long_draws_list$long_draws_count %>%
#   select(-factor) %>%
#   left_join(factor_table, by = "parm")

attach(long_draws_list)

nearly_full_draws_marginal =  long_draws_prop %>%
  # Add in the reference class marginals
  mutate(linear_value = gtools::inv.logit(Intercept),
         parm = NA_character_) %>% distinct() %>% 
  bind_rows(long_draws_prop,
            long_draws_prop_r,
            long_draws_count,
            long_draws_count_r) %>% 
  bind_rows(long_draws_count %>% 
              # Add refernce class to these as well
              mutate(linear_value = exp(Intercept),
                     parm = NA_character_) %>% distinct()) %>% 
  filter(!is.na(factor)) # Remove extra parms that got in
nearly_full_draws_marginal %>% write_rds("out/nearly_full_effect_draws_marg.rds")

loo_draw_ids = read_rds("out/loo_draw_ids.rds")
draw_details = loo_draw_ids %>% select(name, model_txt, draw_list) %>%
  unchop(draw_list) %>% mutate(.draw = 1:n()) %>% 
  select(-draw_list) %>% chop(.draw) %>% 
  unnest_wider(model_txt)

### Pool Worm Lake & Worm Family ####
# For count models where both worm lake & family are present, 
# create a combined worm lake + family
get_joint_worm_draw_ids = function(details) {
  details %>% 
    unchop(.draw) %>% group_by(model) %>% chop(.draw) %>% 
    ungroup() %>% 
    filter(str_detect(model, fixed("worm.lake")),
           str_detect(model, fixed("(1|worm.fam)"))) %>% 
    pull(.draw) %>% unlist()
}
joint_worm_ids = list(
  count = draw_details %>% 
    select(model = pois, .draw),
  prop = draw_details %>% 
    select(model = hu, .draw)
) %>% map(get_joint_worm_draw_ids)

# assign a lake to a parameter
assign_wormlake = function(.x) {
  # browser()
  .x = str_remove(.x, "__hu")
  matcher = tibble(.x = unique(.x, )) %>% 
    mutate(out = case_when(
      is.na(.x)            | grepl("fam[boo", .x, fixed = TRUE) ~ "boot",
      .x == "worm.lakeech" | grepl("fam[ech", .x, fixed = TRUE) ~ "echo",
      .x == "worm.lakegos" | grepl("fam[g", .x, fixed = TRUE)   ~ "gos",
      TRUE ~ NA_character_)
    )
  tibble(.x) %>% left_join(matcher, ".x") %>% pull(out)
}

# Not sure if this is working for interactions?

combine_worm_data_levels = function(.data = nearly_full_draws_marginal,
                                    comp = "count", id_list = joint_worm_ids) {
  ids = id_list[[comp]]
  g = switch(comp, count = log, prop = gtools::logit)
  inv_g = switch(comp, count = exp, prop = gtools::inv.logit)
  
  # Subset of the data that has both worm factors
  worm_data = .data %>% 
    filter(component == comp, 
           factor %in% c("Worm Lake", "Worm Family (RE)")) %>% 
    mutate(worm_lake = assign_wormlake(parm)) %>% 
    select(-parm) 
  joint_data = worm_data %>% filter(.id %in% ids)
  
  # Join the worm factors together (by lake) and get the combined linear value
  lake_family_data = joint_data %>% 
    mutate(value = g(linear_value) - Intercept) %>% 
    distinct() %>% 
    select(-linear_value) %>% 
    pivot_wider(names_from = factor, values_from = value,
                values_fn = list) %>% 
    unchop(c(`Worm Lake`, `Worm Family (RE)`)) %>% 
    mutate(linear_value = inv_g(Intercept + `Worm Lake` + `Worm Family (RE)`)) %>% 
    select(-`Worm Lake`, -`Worm Family (RE)`) %>% 
    mutate(factor = "Worm Lake + Family")  
  
  # Now I need to fill in the missing values so that proportions are correct
  # I'm not really sure what this is doing...
  # 
  lake_family_missing = worm_data %>% filter(! .id %in% ids) %>% 
    filter(factor == "Worm Family (RE)") %>% 
    mutate(factor = "Worm Lake + Family") %>% 
    distinct() %>% 
    mutate(absent = TRUE, linear_value = NA_real_)
  
  # browser()
  lake_family_combined = bind_rows(lake_family_missing, lake_family_data) %>% 
    rename(parm = worm_lake) 
  
  # Define worms (pooled), a factor that includes all of the worm stuff 
  # present in that model
  pooled_worm_data = worm_data %>% 
    filter(!absent, ! .id %in% ids) %>%
    bind_rows(lake_family_data %>% filter(!absent)) %>% 
    mutate(factor = "Worm (pooled)") %>% 
    distinct()
    
  no_worm_data = worm_data %>% 
    filter(factor == "Worm Lake", !.id %in% ids) %>% 
    mutate(factor = "Worm (pooled)") %>% distinct()
  
  pooled_worm_combined = bind_rows(pooled_worm_data, no_worm_data) %>% 
    distinct()
  
  bind_rows(.data, pooled_worm_combined,
            lake_family_data)
}

full_draws_with_worm = nearly_full_draws_marginal %>% 
  combine_worm_data_levels("count") %>% 
  combine_worm_data_levels("prop")


combine_worm_data_levels = function(.data = nearly_full_draws_marginal,
                                    comp = "count", id_list = joint_worm_ids) {
  ids = id_list[[comp]]
  g = switch(comp, count = log, prop = gtools::logit)
  inv_g = switch(comp, count = exp, prop = gtools::inv.logit)
  
  # Subset of the data that has both worm factors
  worm_data = .data %>% 
    filter(component == comp, 
           factor %in% c("Worm Lake", "Worm Family (RE)")) %>% 
    mutate(worm_lake = assign_wormlake(parm)) %>% 
    select(-parm) 
  joint_data = worm_data %>% filter(.id %in% ids)
  
  # Join the worm factors together (by lake) and get the combined linear value
  lake_family_data = joint_data %>% 
    mutate(value = g(linear_value) - Intercept) %>% 
    distinct() %>% 
    select(-linear_value) %>% 
    pivot_wider(names_from = factor, values_from = value,
                values_fn = list) %>% 
    unchop(c(`Worm Lake`, `Worm Family (RE)`)) %>% 
    mutate(linear_value = inv_g(Intercept + `Worm Lake` + `Worm Family (RE)`)) %>% 
    select(-`Worm Lake`, -`Worm Family (RE)`) %>% 
    mutate(factor = "Worm Lake + Family")  
  
  # Now I need to fill in the missing values so that proportions are correct
  # I'm not really sure what this is doing...
  # 
  lake_family_missing = worm_data %>% filter(! .id %in% ids) %>% 
    filter(factor == "Worm Family (RE)") %>% 
    mutate(factor = "Worm Lake + Family") %>% 
    distinct() %>% 
    mutate(absent = TRUE, linear_value = NA_real_)
  
  # browser()
  lake_family_combined = bind_rows(lake_family_missing, lake_family_data) %>% 
    rename(parm = worm_lake) 
  
  # Define worms (pooled), a factor that includes all of the worm stuff 
  # present in that model
  pooled_worm_data = worm_data %>% 
    filter(!absent, ! .id %in% ids) %>%
    bind_rows(lake_family_data %>% filter(!absent)) %>% 
    mutate(factor = "Worm (pooled)") %>% 
    distinct()
  
  no_worm_data = worm_data %>% 
    filter(factor == "Worm Lake", !.id %in% ids) %>% 
    mutate(factor = "Worm (pooled)") %>% distinct()
  
  pooled_worm_combined = bind_rows(pooled_worm_data, no_worm_data) %>% 
    distinct()
  
  bind_rows(.data, pooled_worm_combined,
            lake_family_data)
}
### Pool Copepods ####
# full_draws_with_worm = read_rds("out/full_effect_size_draws_marginal.rds")

get_joint_cop_draw_ids = function(details) {
  # browser()
  details %>% 
    unchop(.draw) %>% group_by(model) %>% chop(.draw) %>% 
    ungroup() %>% 
    filter(str_detect(model, fixed("cop.lake")),
           str_detect(model, fixed("genus"))) %>% 
    pull(.draw) %>% unlist()
}
joint_cop_ids = list(
  count = draw_details %>% 
    select(model = pois, .draw),
  prop = draw_details %>% 
    select(model = hu, .draw)
) %>% map(get_joint_cop_draw_ids)


# full_draws_with_worm %>% filter(!absent, factor == "Copepod Lake") %>% arrange(.id, parm)
# 
# full_draws_with_worm %>% filter(!absent, factor == "Genus") %>% arrange(.id, parm)
assign_copgenus = function(parm) {
  # browser()
  genus_key = tribble(
    ~"parm",       ~"out",
    "genusM"         , "M", 
    "Genus",           "A",
    "cop.lakeech"    , "A",
    "cop.lakegos"    , "A",
    "cop.lakelau"    , "M",
    "cop.lakerob"    , "M",
    "Copepod Lake"  , "A"
  )
  tibble(parm = parm) %>% 
    left_join(genus_key, "parm") %>% pull(out)
}
combine_cop_data_levels = function(.data = nearly_full_draws_marginal,
                                    comp = "count", id_list = joint_cop_ids) {
  ids = id_list[[comp]]
  g = switch(comp, count = log, prop = gtools::logit)
  inv_g = switch(comp, count = exp, prop = gtools::inv.logit)
  
  # if(comp == "prop") browser()
  # NA is either boom (M) or (A)
  cop_data = .data %>% 
    filter(component == comp, 
           factor %in% c("Copepod Lake", "Genus")) %>% 
    # For NA's in parm, key it by the factor name
    # to distinguish when mergin w/ the genus key
    mutate(parm = coalesce(parm, factor),
           Genus_ID = assign_copgenus(parm)) %>% 
    select(-parm) 
  # Subset of the data that has both copepod factors
  joint_data = cop_data %>% filter(.id %in% ids)
  
  # Join the copepod factors together (by lake) and get the combined linear value
  lake_genus_data = joint_data %>% 
    mutate(value = g(linear_value) - Intercept) %>% 
    distinct() %>% 
    select(-linear_value) %>% 
    pivot_wider(names_from = factor, values_from = value,
                values_fn = list) %>% 
    unchop(c(Genus, `Copepod Lake`)) %>% 
    mutate(linear_value = inv_g(Intercept + `Copepod Lake` + `Genus`)) %>% 
    select(-Genus, -`Copepod Lake`) %>% 
    mutate(factor = "Copepod Genus + Lake")  
  
  # Now I need to fill in the missing values so that proportions are correct
  # I'm not really sure what this is doing...
   
  lake_genus_missing = cop_data %>% filter(! .id %in% ids) %>% 
    filter(factor == "Copepod Lake") %>% 
    mutate(factor = "Copepod Genus + Lake")  %>% 
    distinct() %>% 
    mutate(absent = TRUE, linear_value = NA_real_)
  
  # browser()
  lake_genus_combined = bind_rows(lake_genus_missing, lake_genus_data) %>% 
    rename(parm = Genus_ID)  # Not sure How I feel about this one...
  
  # Define worms (pooled), a factor that includes all of the worm stuff 
  # present in that model
  pooled_cop_data = cop_data %>% 
    filter(!absent, ! .id %in% ids) %>%
    bind_rows(lake_genus_data %>% filter(!absent)) %>% 
    mutate(factor = "Copepod (pooled)") %>% 
    distinct()
  
  no_cop_data = cop_data %>% 
    filter(factor == "Genus", absent) %>% 
           # !.id %in% ids) %>% 
    mutate(factor = "Copepod (pooled)") %>% distinct()
  
  pooled_cop_combined = bind_rows(pooled_cop_data, no_cop_data) %>% 
    distinct() %>% rename(parm = Genus_ID)
  # Counts are 5 for lake_genus, 2 for the others, which is what they should be
  bind_rows(.data, pooled_cop_combined,
            lake_genus_data)
}
full_draws_marginal = full_draws_with_worm %>% 
  combine_cop_data_levels("count") %>% 
  combine_cop_data_levels("prop") %>% 
  distinct()
# full_draws_marginal = full_draws_with_worm # consider doing the same for genus + copepod lake

full_draws_marginal %>% write_rds("out/full_effect_size_draws_marginal.rds")  

if(exists("full_draw_marginal")) rep("======= DONE ========\n", 10)
# full_draws_marginal = read_rds("out/full_effect_size_draws_marginal.rds")  


############################################
# Calculate effect size CI's & inclusion frequencies, 
# then save the output
############################################

# Calculate effect sizes for each factor & sample
full_effect_sizes = full_draws_marginal %>% 
  group_by(component, factor, .id) %>% 
  # distinct() %>% 
  summarise(effect_size = sd(linear_value, na.rm = TRUE),
            absent = any(is.na(linear_value)))

write_rds(full_effect_sizes, "out/full_effect_sizes.rds")
# full_effect_sizes = read_rds("out/full_effect_sizes.rds")

get_ci = function(.x, probs = c(.025, .05, .25, .5, .75, .95, .975)){
  quants = quantile(.x, probs, na.rm = TRUE)
  names(quants) = paste0("Q", substring(probs, 2))
  tibble(!!!quants)
}
# N_DRAWS = draw_details %>% unchop(.draw) %>% nrow()
# Get the credible intervals of effect sizes & the frequency of the factor's inclusion in models
effect_size_ci = full_effect_sizes %>% 
  select(-.id) %>% 
  chop(c(effect_size, absent)) %>% 
  mutate(
    frequency = 1 - map_dbl(absent, mean),
    # frequency2 = 1 - map_dbl(absent, sum)/N_DRAWS,
    ci = map(effect_size, get_ci)
  ) %>% 
  select(-effect_size, -absent) %>% unnest(ci) %>% 
  ungroup() %>% 
  arrange(component, desc(frequency), factor) 


# Save it
write_csv(effect_size_ci, "out/marginal_effect_sizes.csv")  

effect_size_ci = read_csv("out/marginal_effect_sizes.csv")  

total_n = n_distinct(full_effect_sizes$.id)

# Format for manuscript
pretty_ci = effect_size_ci %>% 
  filter(!factor %in% c("Worm Lake + Family", "Copepod Genus + Lake")) %>% # Something is bugged with that one
  select(component:frequency, median = Q.5, lo = Q.025, hi = Q.975) %>% 
  mutate(across(frequency:hi, ~format(round(.x,3), digits = 3))) %>% 
  mutate(CI = glue("[{lo}, {hi}]")) %>% 
  select(-lo, -hi) %>% 
  mutate(factor = recode(factor, Genus = "Copepod Genus")) %>% 
  as.data.frame()

pretty_ci
pretty_ci %>% write_csv("writing/effect_sizes.csv")

### Results Table of effect sizes ####
# component                          factor frequency median             CI
#     count                Copepod (pooled)     0.849  0.640 [0.421, 1.289]
#     count                   Copepod Genus     0.849  0.621 [0.019, 2.201]
#     count                   Worm (pooled)     0.817  0.187 [0.027, 1.297]
#     count                       Worm Lake     0.817  0.182 [0.023, 1.374]
#     count                    Copepod Lake     0.522  0.503 [0.190, 2.464]
#     count                      Plate (RE)     0.495  0.095 [0.004, 0.421]
#     count                          Native     0.473  0.126 [0.005, 0.798]
#     count                Worm Family (RE)     0.420  0.113 [0.006, 0.575]
#     count   Genus x Worm Lake Interaction     0.354  0.616 [0.052, 4.309]
#     count Copepod x Worm Lake Interaction     0.208  0.721 [0.079, 7.482]
#      prop                Copepod (pooled)     0.871  0.255 [0.213, 0.364]
#      prop                   Copepod Genus     0.871  0.274 [0.043, 0.385]
#      prop                   Worm (pooled)     0.840  0.058 [0.020, 0.118]
#      prop                       Worm Lake     0.840  0.060 [0.018, 0.128]
#      prop                    Copepod Lake     0.593  0.143 [0.072, 0.242]
#      prop                      Plate (RE)     0.497  0.032 [0.002, 0.078]
#      prop                Worm Family (RE)     0.409  0.018 [0.001, 0.060]
#      prop                          Native     0.372  0.039 [0.002, 0.142]
#      prop Copepod x Worm Lake Interaction     0.235  0.095 [0.045, 0.164]
#      prop   Genus x Worm Lake Interaction     0.218  0.053 [0.009, 0.175]

#### Copepod genus effects ####
copepod_genus_effects = local({
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
  split_genus_effects %>% 
    select(component, partition, median = Q.5, lo = Q.025, hi = Q.975) %>% 
    mutate(across(median:hi, ~format(round(.x,3), digits = 3))) %>% 
    mutate(CI = glue("[{lo}, {hi}]")) %>% 
    select(-lo, -hi) %>% arrange(component, partition)
  
  
})
# These are all for genus w/ or w/o cop lake
# component partition median CI            
# count     Lake      0.342  [0.011, 2.714]
# count     No Lake   0.712  [0.518, 0.958]
# prop      Lake      0.210  [0.031, 0.395]
# prop      No Lake   0.338  [0.303, 0.375]

# Break down incidence variables by co-occurrence in ensemble
# full_effect_sizes %>% 
#   ungroup() %>% 
#   filter(component == "count", factor != "Copepod Lake") %>% 
#   mutate(included = !absent) %>%
#   select(factor, included, .id) %>% 
#   distinct() %>% 
#   pivot_wider(names_from = factor, values_from = included) %>% 
#   select(-.id) %>% 
#   group_by(across(everything())) %>% 
#   summarise(count = n()) %>% 
#   ungroup() %>% 
#   mutate(freq = count / total_n ) %>% 
#   arrange(desc(freq)) %>% 
#   nest(fcts = 1:4) %>% 
#   mutate(terms = map_chr(fcts, 
#     ~names(.x)[as.logical(.x)] %>% paste(collapse = " + "))) %>% 
#   select(-fcts, -count)

### Include this in the manuscript ####

# worm_inclusion  freq
# both            0.209
# family          0.342
# lake            0.191
# none            0.258



# Table X: Importance of each factor for the prevalence (rate of infection) and 
  # intensity (number of worms, given infection) model components.
  # Inclusion is the weighted frequency of models that included the factor.
  # Effect size was calculated as the standard deviation of the factor's marginal effects 
  # for each posterior sample in which the factor was included; 
  # effect sizes are presented with medians and 95% credible intervals 
  # and are in units of worms (Intensity) or infection frequency (prevalence).
  # Random effects are indicated with "(RE)"

# Component     Factor           Inclusion    Effect Size [95% CI]            
# Prevalence    Copepod Lake        1         0.244  [0.220, 0.267]
# Prevalence    Worm Lake           1         0.036  [0.017, 0.058]
# Prevalence    Native              0.252     0.013  [0.001, 0.036]
# Prevalence    Plate (RE)          0.0699    0.019  [0.001, 0.048]
# Prevalence    Worm Family (RE)    0.0283    0.022  [0.002, 0.050]
# Intensity     Copepod Lake        1         0.510  [0.384, 0.664]
# Intensity     Worm Family (RE)    0.551     0.041  [0.004, 0.116]
# Intensity     Worm Lake           0.399     0.051  [0.013, 0.128]
# Intensity     Native              0.391     0.044  [0.005, 0.112]
# Intensity     Plate (RE)          0.223     0.030  [0.002, 0.105]

# Updated w/ pooled worm intensity estimates:
# Component     Factor           Inclusion    Effect Size [95% CI]            
# Prevalence    Copepod Lake          1.000   0.244  [0.220, 0.267]
# Prevalence    Worm Lake             1.000   0.036  [0.017, 0.058]
# Prevalence    Native                0.252   0.013  [0.001, 0.036]
# Prevalence    Plate (RE)            0.070   0.019  [0.001, 0.048]
# Prevalence    Worm Family (RE)      0.028   0.022  [0.002, 0.050]
# Intensity     Copepod Lake          1.000   0.510  [0.384, 0.664]
# Intensity     Worm (pooled)         0.742   0.047  [0.010, 0.109]
# Intensity     Native                0.391   0.044  [0.005, 0.112]
# Intensity     Worm Family (RE)      0.342   0.042  [0.006, 0.103]
# Intensity     Plate (RE)            0.223   0.030  [0.002, 0.105]
# Intensity     Worm Lake + Family    0.209   0.054  [0.020, 0.118]
# Intensity     Worm Lake             0.191   0.050  [0.016, 0.111]



# Filter out the damn genus:worm.lake interactions
# 


# This tries to parse apart the non-copepod terms for the intensity model:
# I'm not sure if this belongs in the manuscript

# Frequency Terms
# 0.190     Worm Lake                                           
# 0.175     Worm Family (RE)                                    
# 0.171     Native                                              
# 0.112     Native + Worm Family (RE)                           
# 0.110     Worm Family (RE) + Worm Lake                        
# 0.0872    Native + Plate (RE)                                 
# 0.0800    Plate (RE) + Worm Family (RE) + Worm Lake           
# 0.0537    Plate (RE) + Worm Family (RE)                       
# 0.0190    Native + Worm Family (RE) + Worm Lake               
# 0.00178   Native + Plate (RE) + Worm Family (RE)              
# 0.000381  None                                                  
# 0.000212  Plate (RE) + Worm Lake                              
# 0.000127  Native + Worm Lake                                  
# 0.0000424 Native + Plate (RE) + Worm Family (RE) + Worm Lake  
# 
# 
# 
# Be clear that this is mostly confounded with copepod species. 
  # The exception is Boot and Lawier. This is the one case where 
  # you have the same species in two lakes. So it might be worth 
  # a focused analysis paragraph specifically on these two. 
  # Do the worms from Boot Lake do better on Boot copepods, or 
  # Lawier copepods, or similarly on the two? This is your best 
  # shot at really separating out the LA issue, from the HS issue.

# 
# 