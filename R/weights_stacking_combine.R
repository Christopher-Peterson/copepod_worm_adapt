suppressPackageStartupMessages({
  library(readr)
  library(gtools)
  library(readr)
  library(purrr)
  library(dplyr)
  library(tidyverse)
  library(tidyr)
  library(glue)
  library(brms)
  library(stringr)
})
N_CORES = 48L
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

stacked_draws = read_rds( "out/combined_draws_df.rds")

# Replace NA's with 0 in stacked_draws
stacked_draws_nz = stacked_draws %>% 
  mutate(across(where(is.double), rlang::`%|%`, 0))

# loo_stacking = read_csv("out/loo_stacking.csv")

unique_data = raw_data %>% 
  mutate(id = 1:n()) %>%
  group_by(cop.lake, worm.fam, worm.lake, plate, genus) %>%
  filter(id == min(id)) %>% select(id) %>% ungroup() %>% 
  mutate(cop.lake = unorder_cop_lake(cop.lake),
         native = cop.lake == worm.lake)

# Dummy code the fixed effects
# model_matrix = 
#   unique_data %>% select(id, worm.fam, worm.lake) %>% 
#   bind_cols(
#     as_tibble(model.matrix(~cop.lake + worm.lake + genus + native, data = unique_data) ))



# Helper Functions ####
ev_hu = function(.x) {
  # x = .x[draws, par_ids , drop = FALSE]
  # BRMS has hurdle model linear predictor on reverse scale
  1 - gtools::inv.logit(.x)
}
ev_count = function(.x) {
  # x = .x[draws, par_ids , drop = FALSE]
  # Truncated poisson expected value (see Wikipedia entry for formula)
  # .x is log_lambda
  lambda = exp(.x)
  lambda / (1 - exp(-lambda))
}

# Convert a vector of parameters names into an expression that
# calculates expected value
#' @param parms character vector of parameter names
#' @param ev_func either("ev_hu" or "ev_count")
#' @return an unevaluated expression to be run in the context of stacked_draws
parse_ev_expr = function(parms, ev_func) {
  lin_pred = glue("`{parms}`") %>% # wrap in backticks to protect non-syntactic names
    glue_collapse(" + ") # linear predictor is sum of parms
  glue("{ev_func}({lin_pred})") %>%  # Apply ev_func()
  rlang::parse_expr() # Convert to an unevaluated expression
}
# Designed to be pmapped to unique data;
# returns a list of expected value expressions

get_expected_value_exprs = function(cop.lake, worm.fam, worm.lake, plate, 
                               genus, native, id, draws = stacked_draws_nz) {
    # browser()
  native = as.character(native)
  parms = tribble(~"parm", ~"var", ~"post",
         "b_Intercept",     "",          "",   
         "b_hu_Intercept",  "",          "",                 
         "b_cop.lake",      cop.lake,    "",   
         "b_worm.lake",     worm.lake,   "",    
         "b_genus",         genus,       "",
         "b_native",        native,      "", 
         "b_hu_cop.lake",   cop.lake,    "",      
         "b_hu_worm.lake",  worm.lake,   "",       
         "b_hu_genus",      genus,       "",   
         "b_hu_native",     native,      "",    
         "r_plate[",        plate,       ",Intercept]",
         "r_worm.fam[",     worm.fam,    ",Intercept]",
         "r_plate__hu[",    plate,       ",Intercept]",
         "r_worm.fam__hu[", worm.fam,    ",Intercept]") %>% 
    glue_data("{parm}{var}{post}") %>% 
    as.character()
  # Some of these are invalid options; remove them
  keep_parms = parms[parms %in% posterior::variables(draws)]
  # separate out the two parameter types
  hu_parms_lgl = keep_parms %>% str_detect("_hu")
  # Create a set of unevaluated expressions that will find the linear preds
  list(
    id = id,
    hu_expr = parse_ev_expr(keep_parms[hu_parms_lgl], "ev_hu"),
    count_expr = parse_ev_expr(keep_parms[!hu_parms_lgl], "ev_count")
  )
}
exec_ev_exprs = function(expr_list, draws = stacked_draws_nz, ...) {
  # browser()
  with(expr_list, draws %>%
      transmute(count_ev = !!count_expr,
                hu_ev = !!hu_expr,
                id = id, 
                .draw = 1:n())
    )
}



### TEST BLOCK ####
# 
# test_mdls = c(huh_mdl = "combo_run_4766", bad_mdl = "combo_run_8326")
# 
# draw_details = loo_draw_ids %>% select(name, model_txt, draw_list) %>%
#   unchop(draw_list) %>% mutate(.draw = 1:n())
# test_draws = draw_details %>% filter(name %in% test_mdls) %>% pull(.draw)
# stacked_draws_test = stacked_draws %>% 
#   mutate(.draw = 1:n()) %>% 
#   filter(.draw %in% test_draws)
# 
# test_ev_exprs = pmap(unique_data, get_expected_value_exprs, draws = stacked_draws_test)
# test_ev_list = lapply(test_ev_exprs, exec_ev_exprs, draws = stacked_draws_test)


# Run Functions ####
# Returns a list that
ev_exprs = pmap(unique_data, get_expected_value_exprs)
ev_list = lapply(ev_exprs, exec_ev_exprs, mc.cores = N_CORES)

expected_values = bind_rows(ev_list) %>% 
  left_join(unique_data, by = "id")
write_rds(expected_values, "out/expected_values_long.rds")
expected_values = read_rds("out/expected_values_long.rds")
# Next, I need to pool worm fam within worm lake,
# # genus within cop lake & drop plate

# Create summaries
summarise_one_ev = function(ev_vec, .resp = c("count", "hu")) {
  # ev_vec = ev_list[[1]]$count_ev
  quantiles = c(.99, .95, .9, .8, .5)
  full_quantiles = c((1-quantiles)/2, rev(quantiles + (1-quantiles)/2))
  out_vec = c(
    quantile(ev_vec, full_quantiles), 
    mean = mean(ev_vec), median = median(ev_vec))
  tibble(!!!out_vec, .resp = .resp)
}

summarize_ev = function(ev_df) {
  bind_rows(
    summarise_one_ev(ev_df$count_ev, "count"),
    summarise_one_ev(ev_df$hu_ev, "hu")
  ) %>% mutate(id = ev_df$id[[1]])
}

expected_value_smry = map_dfr(ev_list, summarize_ev)
expected_value_full_smry = expected_value_smry %>% left_join(unique_data, by = "id")
 
expected_value_full_smry %>% write_rds("out/expected_values_smry.rds")


# not_na_props = expected_values %>% 
#   group_by(cop.lake, worm.fam, worm.lake, plate, genus, native) %>% 
#   summarise(prop_count = sum(!is.na(count_ev)) / n(),
#             prop_hu = sum(!is.na(hu_ev)) / n())
# 
# expected_values$.draw %>% n_distinct()
# 
# 
# draw_details = loo_draw_ids %>% select(name, model_txt, draw_list) %>%
  # unchop(draw_list) %>% mutate(.draw = 1:n())
# 
# draw_na = expected_values %>% group_by(.draw) %>% 
#   summarise(prop_count = sum(!is.na(count_ev)) / n(),
#             prop_hu = sum(!is.na(hu_ev)) / n())
# na_info = left_join(draw_details, draw_na, by = ".draw")
# tmp = draw_na %>%
#   select(-.draw) %>% 
#   mutate_all(~round(.x, 3)) %>% 
#   table()
# tmp / sum(tmp)
# 
# good_models = draw_na %>% filter( prop_hu > 0.8, prop_count > 0.8) %>% 
#   left_join(draw_details, by = ".draw") %>%
#   pull(model_txt) %>% unique() %>% 
#   unlist() %>%
#   unique() %>% 
#   str_remove(".+ ~ ") %>% 
#   unique() %>% 
#   sort() 
# 
# tst = draw_na %>% filter(prop_hu < 0.2, prop_hu > 0) %>% 
#   left_join(expected_values, by = ".draw") 
# tst %>%   group_by(cop.lake, native) %>% 
#   summarise(prop_hu = sum(!is.na(hu_ev)))
# ALL OF THEM ARE IN THE REFERENCE CATEGORIES! 
  # YES!!!





# 
# bad_models = draw_na %>% filter( prop_hu == 0, prop_count == 0) %>% 
#   left_join(draw_details, by = ".draw") %>%
#   pull(model_txt) %>% unique() %>% 
#   unlist() %>%
#   unique() %>% 
#   str_remove(".+ ~ ") %>% 
#   unique() %>% 
#   sort() 
# 
# good_models %in% bad_models
# 
# 
# draw_na %>% 
#   filter( prop_hu < 1, prop_hu > 0) %>% 
#   left_join(draw_details, by = ".draw") %>%
#   pull(model_txt) %>% 
#   map_chr("hu") %>% 
#   unique() %>% 
#   str_remove(".+ ~ ") %>% 
#   unique() %>% 
#   sort() 
# huh_models = draw_na %>% 
#   filter( prop_count < 1, prop_count > 0) %>% 
#   left_join(draw_details, by = ".draw") %>%
#   pull(name) %>% 
#   pull(model_txt) %>% 
#   map_chr("pois") %>% 
#   unique() %>% 
#   str_remove(".+ ~ ") %>% 
#   unique() %>% 
#   sort() 
# 
# essential_terms = good_models %>% str_split(fixed(" + ")) %>% 
#   reduce(~.x[.x %in% .y])

# Essential Terms:
# "cop.lake"  "worm.lake"  "genus"  "native" "(1|plate)" "(1|worm.fam)"
# All Present in all of the good models;
# none of the bad or huh models had them all present
# 
# 
# huh_terms = huh_models %>% str_split(fixed(" + ")) %>% 
#   reduce(~.x[.x %in% .y]) # Missing cop.lake and native
# essential_terms[! essential_terms %in% huh_terms]
# 
# 
# huh_models %>% str_split(fixed(" + ")) %>% 
#   map_lgl(~all(essential_terms %in% .x))



# SO I need to get 
  # A) Model components
  # B) Summarized expected values
  # C) Effect Size Table
  # D) Make the figure w/ expected values
  # 



## Define Functions ####
#' 
#' #' Determine which draws to use for stacking weights
#' #' @param stacking_df data farme containing column `weight`
#' #' @param n_samples number of samples per model
#' #' @param n_draws total number of draws for the full model; if NA, goes for the max possible
#' #' @return `stacking_df` with an extra `draws` column, wihich is a list of ints
#' add_stacking_draws <- function(stacking_df, n_samples, n_draws = NA) {
#'   K <- nrow(stacking_df)
#'   if(is.na(n_draws)) n_draws = n_samples / max(stacking_df$weight)
#'   w <- round(n_draws * stacking_df$weight) # expected number of draws from each model
#'   
#'   draws = purrr::map(w, ~sort(sample(seq_len(n_samples), 
#'                                      size = .x, replace = FALSE)))
#'   stacking_df$draws = draws
#'   stacking_df
#' }
#' 
#' # Convert linpreds to expected values for both number & frequency of worms
#' 
#' #' Read in a model, & get the expected values for selected draws
#' #' Meant to be used with pmap_dfr() and the output of add_stacking_draws()
#' #' @param name model name
#' #' @param draws list of indices to keep from model
#' #' @param key_data data frame of unique observation combinations
#' #' @param model_dir location of model fits
#' #' @param model_suffix filename suffix of model fit rds files
#' #' @return a data frame of expected values for both parameters
#' get_expected_value_draws = function(name, draws, ..., 
#'                                     key_data = unique_data,
#'                                     model_dir =  "out/model_runs/", 
#'                                     model_suffix = "_model.rds")  {
#'   model_name = paste0(model_dir, name, model_suffix)
#'   model_fit = read_rds(model_name)
#'   key_ids = key_data$id
#'   n_keys = length(key_ids)
#'   if(length(draws) == 0) {
#'     return(NULL)
#'   }
#'     # These are the linear predictors
#'   num_worms = brms::posterior_linpred(
#'       model_fit, scale = "linear", summary = FALSE,
#'       subset = draws, newdata = key_data) %>%
#'     ev_num(draws, key_ids)
#'   worm_freq = brms::posterior_linpred(
#'       model_fit, scale = "linear", dpar="hu", summary = FALSE,
#'       subset = draws, newdata = key_data) %>%
#'     ev_freq(draws, key_ids)
#'   
#'   # browser()
#'   tibble(num_worms = as.numeric(num_worms), 
#'          worm_freq = as.numeric(worm_freq)) %>% 
#'     mutate(.draw = seq_along(draws) %>% rep(n_keys),
#'            id = rep(key_ids, each = length(draws)),
#'            model = name) %>% 
#'     left_join(key_data, by = "id")
#' }
#' 
#' # pmap_dfr in parallel
#' par_pmap_dfr = function(.l, .f,  .cores = 48) {
#'   output = rlang::exec(parallel::mcmapply, FUN = .f, !!!.l, mc.cores = .cores)
#'   bind_rows(output)
#' }
#' #### assemble effect table ####
#' 
#' 
#' sum_weight_by_term = function(data, term) {
#'   term = enexpr(term)
#'   term_type = as.character(as.symbol(term))
#'   data %>% 
#'     select(weight, !!term) %>% 
#'     unchop(!!term) %>% 
#'     mutate(term := trimws(!!term)) %>% 
#'     group_by(term) %>% 
#'     summarise("{term_type}_weight" := sum(weight)) %>% 
#'     ungroup()
#' }
#' term_weights = local({
#'   split_names = loo_stacking %>% 
#'     mutate(across(c(hu, numb.worm),
#'                   str_split, pattern = fixed("+")))
#'   left_join(
#'     split_names %>% sum_weight_by_term(numb.worm),
#'     split_names %>% sum_weight_by_term(hu),
#'     by = "term")
#' })
#' write_csv(term_weights, "out/stacking_term_weights.csv")
#' 
#' 
#' #### Run it ####
#' max_draws = 5000000L
#' set.seed(109342); loo_with_draws = 
#'   add_stacking_draws(loo_stacking,n_samples = 4000L, n_draws = max_draws) %>%
#'   mutate(n = map_int(draws, length)) %>% filter(n > 0)
#' write_rds(loo_with_draws, "out/loo_stacking_draws.rds")
#' 
#' dir.create("out/stacking_chunks")
#' dir.create("out/stacking_chunks_settings")
#' # Split into chunks
#' N_chunks = 500
#' chunked_draws = loo_with_draws %>% mutate(chunk = (1:n()) %% N_chunks) %>% group_split(chunk)
#' walk(chunked_draws, ~write_rds(.x, glue::glue("out/stacking_chunks_settings/chunk_{chunk}.rds", chunk = .x$chunk[1])))
#' 
#' jobfile = paste("R-stan R/weights_stacking_draw_chunks.R", 1:500)
#' writeLines(jobfile, "jobs/weights_stacking_draw_chunks.job")
# full_evs = pmap_dfr(loo_with_draws, get_expected_value_draws)

# write_csv(full_evs, "out/stacked_evs.csv")


# full_draws = pmap_dfr(loo_with_draws, pick_draws)
