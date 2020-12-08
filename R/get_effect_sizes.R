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
stacked_draws_col = read_csv("out/stacked_draws.csv", n_max = 3)

# Define input column specification, replacing NA's with logicals
stacked_colspec = map_chr(stacked_draws_col, ~{
    type = typeof(.x)
    if(type == "logical") return("d")
    substr(type, 1, 1)
  }) %>% 
  paste(collapse = "")
# Read the full data
stacked_draws = read_csv("out/stacked_draws.csv",col_types = stacked_colspec)

# Utility function; Names factors based on parm names
name_factors = function(x) {
  case_when(
    str_detect(x,"r_worm.fam") ~ "Worm Family (RE)",
    str_detect(x,"r_plate") ~  "Plate (RE)",
    str_detect(x,"cop.lake") ~ "Copepod Lake",
    str_detect(x,"worm.lake") ~ "Worm Lake",
    str_detect(x,"native") ~    "Native")
}
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
  select(1:3, model, contains("_hu"))
stacked_draws_count = stacked_draws_relevant %>% 
  select(-contains("_hu"))

# These next lines pivot the data into along form, 
# calculate marginal effects, 
  # which are basically g^-1( beta + Intercept )
# Then standardizes some formatting

long_draws_count = stacked_draws_count %>% 
  select(-model, -starts_with("sd_"), -(1:3)) %>% 
  rename(Intercept = b_Intercept) %>% 
  mutate(.id = 1:n()) %>% 
  pivot_longer(-c(.id, Intercept),
               names_to = "parm",
               names_prefix = "b_",
               values_to = "value") %>% 
  mutate(absent = is.na(value)) %>% 
  mutate(linear_value = exp(Intercept + value)) %>% 
  select(-value) %>% 
  mutate(factor = name_factors(parm),
         component = "count") 

long_draws_prop = stacked_draws_prop %>% 
  select(-model, -starts_with("sd_"), -(1:3)) %>% 
  rename(Intercept = b_hu_Intercept) %>% 
  mutate(.id = 1:n()) %>% 
  pivot_longer(-c(.id, Intercept),
               names_to = "parm",
               names_prefix = "b_hu_",
               values_to = "value") %>% 
  mutate(absent = is.na(value)) %>% 
  mutate(linear_value = gtools::inv.logit(Intercept + value)) %>% 
  select(-value) %>% 
  mutate(factor = name_factors(parm),
         component = "prop") 

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
  mutate(factor = name_factors(parm),
         component = "count") 


long_draws_prop_r = stacked_draws %>% 
  select(Intercept = b_hu_Intercept, 
         starts_with("r_"), contains("_hu")) %>% 
  mutate(.id = 1:n()) %>% 
  pivot_longer(-c(.id,Intercept),
               names_to = "parm",
               values_to = "value") %>% 
  mutate(absent = is.na(value)) %>% 
  mutate(linear_value = gtools::inv.logit(Intercept + value)) %>% 
  select(-value) %>% 
  mutate(factor = name_factors(parm),
         component = "prop") 


# Pull all of the marginals together
nearly_full_draws_marginal = long_draws_prop %>%
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

# For count models where both worm lake & family are present, 
# create a combined worm lake + family
worm_fam_count_ids = nearly_full_draws_marginal %>% 
  filter(factor %in% c("Worm Lake", "Worm Family (RE)"),
         component == "count") %>% 
  select(-linear_value, - parm) %>% 
  distinct() %>% 
  pivot_wider(names_from = factor, values_from = absent) %>% 
  filter(!(`Worm Lake` | `Worm Family (RE)`)) %>% 
  pull(.id)

# assign a lake to a parameter
assign_wormlake = function(.x) {
  case_when(
    is.na(.x)            | grepl("fam[boo", .x, fixed = TRUE) ~ "boot",
    .x == "worm.lakeech" | grepl("fam[ech", .x, fixed = TRUE) ~ "echo",
    .x == "worm.lakegos" | grepl("fam[g", .x, fixed = TRUE)   ~ "gos",
    TRUE ~ NA_character_)
}

worm_lake_fam_dat = nearly_full_draws_marginal %>% 
  filter(.id %in% worm_fam_count_ids, component == "count",
         factor %in% c("Worm Lake", "Worm Family (RE)")) %>% 
  mutate(worm_lake = assign_wormlake(parm),
         value = log(linear_value) - Intercept) %>% 
  select(-parm, -linear_value) %>% 
  pivot_wider(names_from = factor, values_from = value,
              values_fn = list) %>% 
  unchop(c(`Worm Lake`, `Worm Family (RE)`)) %>% 
  mutate(linear_value = exp(Intercept + `Worm Lake` + `Worm Family (RE)`)) %>% 
  select(-`Worm Lake`, -`Worm Family (RE)`) %>% 
  rename(parm = worm_lake) %>% 
  mutate(factor = "Worm Lake + Family")
  
# Now I need to fill in the missing values so that proportions are correct
worm_lake_fam_missing = 
  nearly_full_draws_marginal %>% 
  filter(.id %in% worm_fam_count_ids, component == "count",
         factor %in% c("Worm Lake", "Worm Family (RE)")) %>% 
  bind_rows(
    nearly_full_draws_marginal %>% 
      filter(! (.id %in% worm_fam_count_ids), 
             component == "count",
             factor == c("Worm Family (RE)")) %>% 
      mutate(factor = "Worm Lake + Family") %>% 
      distinct()
  ) %>% mutate(absent = TRUE, linear_value = NA_real_)
  

# also add a worm (pooled) option
pool_worms = function(data) {
  worm_dat = data %>%  
    filter(component == "count",
      factor %in% c("Worm Lake", "Worm Lake + Family", "Worm Family (RE)"))
  pooled_worm = worm_dat %>% filter(!absent) %>% 
    mutate(factor = "Worm (pooled)")
  new_dat = worm_dat %>% filter(! .id %in% pooled_worm$.id) %>% 
    mutate(factor = "Worm (pooled)") %>% distinct() %>% 
    bind_rows(pooled_worm)
  bind_rows(data, new_dat)
}

full_draws_marginal = nearly_full_draws_marginal %>% 
  # Remove lake_fam ids
  filter(!(component == "count" &
           .id %in% worm_fam_count_ids &
           factor %in% c("Worm Lake", "Worm Family (RE)"))) %>% 
  bind_rows(worm_lake_fam_missing, worm_lake_fam_dat) %>% 
  pool_worms()
  

############################################
# Calculate effect size CI's & inclusion frequencies, 
# then save the output
############################################

# Calculate effect sizes for each factor & sample
full_effect_sizes = full_draws_marginal %>% 
  group_by(component, factor, .id) %>% 
  summarise(effect_size = sd(linear_value, na.rm = TRUE),
            absent = any(is.na(linear_value)))


get_ci = function(.x, probs = c(.025, .05, .25, .5, .75, .95, .975)){
  quants = quantile(.x, probs, na.rm = TRUE)
  names(quants) = paste0("Q", substring(probs, 2))
  tibble(!!!quants)
}

# Get the credible intervals of effect sizes & the frequency of the factor's inclusion in models
effect_size_ci = full_effect_sizes %>% 
  select(-.id) %>% 
  chop(c(effect_size, absent)) %>% 
  mutate(
    frequency = 1 - map_dbl(absent, mean),
    ci = map(effect_size, get_ci)
  ) %>% 
  select(-effect_size, -absent) %>% unnest(ci) %>% 
  ungroup() %>% 
  arrange(component, desc(frequency), factor) 

# Save it
write_csv(effect_size_ci, "out/marginal_effect_sizes.csv")  

total_n = n_distinct(full_effect_sizes$.id)

# Format for manuscript
effect_size_ci %>% 
  select(component:frequency, median = Q.5, lo = Q.025, hi = Q.975) %>% 
  mutate(across(median:hi, ~format(round(.x,3), digits = 3))) %>% 
  mutate(CI = glue("[{lo}, {hi}]")) %>% 
  select(-lo, -hi)

# Parse out whether worm lake, worm family, or both were in the intensity model
stacked_draws %>% 
  select(lake = b_worm.lakeech, fam = sd_worm.fam__Intercept) %>% 
  mutate_all(~!is.na(.x)) %>% 
  mutate(worm_inclusion = case_when(
    lake & fam ~ "both",
    lake ~ "lake",
    fam ~ "family",
    TRUE ~ "none")) %>% 
  group_by(worm_inclusion) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(freq = count/ sum(count))



# Break down incidence variables by co-occurrence in ensemble
full_effect_sizes %>% 
  ungroup() %>% 
  filter(component == "count", factor != "Copepod Lake") %>% 
  mutate(included = !absent) %>%
  select(factor, included, .id) %>% 
  distinct() %>% 
  pivot_wider(names_from = factor, values_from = included) %>% 
  select(-.id) %>% 
  group_by(across(everything())) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(freq = count / total_n ) %>% 
  arrange(desc(freq)) %>% 
  nest(fcts = 1:4) %>% 
  mutate(terms = map_chr(fcts, 
    ~names(.x)[as.logical(.x)] %>% paste(collapse = " + "))) %>% 
  select(-fcts, -count)

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