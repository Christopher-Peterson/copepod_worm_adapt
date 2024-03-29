suppressPackageStartupMessages({
  library(readr)
  library(gtools)
  library(readr)
  library(purrr)
  library(dplyr)
  library(tidyr)
  library(brms)
  library(rlang)
  library(stringr)
  library(glue)
})

if(!exists("argv")) argv = commandArgs(TRUE)

chunk = argv[1] %>% as.integer() %|% 1L
chunk_settings = read_rds(glue("out/stacking_chunks_settings/chunk_{chunk}.rds"))
out_file = glue("out/stacking_chunks/chunk_{chunk}.rds")


# DO NOT RUN DIRECTLY; This is called by get_ev.r
#### Setup ####
order_cop_lake = function(cop_lakes) {
  # Convert cop.lake into a factor
  lake_order = list(gos = "Gossling", ech = "Echo",
                    rob = "Roberts", boo = "Boot", 
                    lau = "Lawier") %>% rev()
  rlang::exec(recode, .x = cop_lakes, !!!lake_order) %>%
    factor(levels = lake_order)
}
unorder_cop_lake = function(cop_lakes) {
  # Convert cop.lake into a factor
  lake_order = c(Gossling = "gos", Echo = "ech",
                 Roberts = "rob", Boot = "boo", 
                 Lawier = "lau") %>% sort()
  rlang::exec(recode, .x = as.character(cop_lakes),
              !!!lake_order) #%>%
  # factor(levels = lake_order)
}

raw_data = read_csv("data/chapter_2_copepod_for_bayes.csv") %>%
  mutate(cop.lake = order_cop_lake(cop.lake))



unique_data = raw_data %>% 
  mutate(id = 1:n()) %>%
  group_by(cop.lake, worm.fam, worm.lake, plate, genus) %>%
  filter(id == min(id)) %>% select(id) %>% ungroup() %>% 
  mutate(cop.lake = unorder_cop_lake(cop.lake),
         native = cop.lake == worm.lake)

## Define Functions ####

# Convert linpreds to expected values for both number & frequency of worms
ev_freq = function(.x, draws, par_ids) {
  # x = .x[draws, par_ids , drop = FALSE]
  # BRMS has hurdle model linear predictor on reverse scale
  1 - gtools::inv.logit(.x)
}
ev_num = function(.x, draws, par_ids) {
  # x = .x[draws, par_ids , drop = FALSE]
  # Truncated poisson expected value (see Wikipedia entry for formula)
  # .x is log_lambda
  lambda = exp(.x)
  lambda / (1 - exp(-lambda))
}

#' Read in a model, & get the expected values for selected draws
#' Meant to be used with pmap_dfr() and the output of add_stacking_draws()
#' @param name model name
#' @param draws list of indices to keep from model
#' @param key_data data frame of unique observation combinations
#' @param model_dir location of model fits
#' @param model_suffix filename suffix of model fit rds files
#' @return a data frame of expected values for both parameters
get_expected_value_draws = function(name, draws, ..., 
                                    key_data = unique_data,
                                    model_dir =  "out/model_runs/", 
                                    model_suffix = "_model.rds")  {
  model_name = paste0(model_dir, name, model_suffix)
  model_fit = read_rds(model_name)
  key_ids = key_data$id
  n_keys = length(key_ids)
  if(length(draws) == 0) {
    return(NULL)
  }
  # These are the linear predictors
  num_worms = brms::posterior_linpred(
    model_fit, scale = "linear", summary = FALSE,
    subset = draws, newdata = key_data) %>%
    ev_num(draws, key_ids)
  worm_freq = brms::posterior_linpred(
    model_fit, scale = "linear", dpar="hu", summary = FALSE,
    subset = draws, newdata = key_data) %>%
    ev_freq(draws, key_ids)
  
  # browser()
  tibble(num_worms = as.numeric(num_worms), 
         worm_freq = as.numeric(worm_freq)) %>% 
    mutate(.draw = seq_along(draws) %>% rep(n_keys),
           id = rep(key_ids, each = length(draws)),
           model = name) %>% 
    left_join(key_data, by = "id")
}

#### Run it ####

full_draws = pmap_dfr(chunk_settings, get_expected_value_draws)
write_rds(full_draws, out_file)