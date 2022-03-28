# Runs full bayesian stacking on the subset of models with decent PBMA+ weights

N_CORES = 48
suppressPackageStartupMessages({
  library(readr)
  library(purrr)
  library(dplyr)
  library(glue)
  library(rlang)
  library(tidyr)
  library(posterior)
})

if(!exists("argv")) argv <- commandArgs(TRUE)
threshold = argv[1] %>% as.numeric() %|% 0
seed = argv[2] %>% as.integer() %|% 572806L

options(mc.cores = N_CORES)
source("R/weights_functions.r")

run_settings = read_rds("model_settings_filtered.rds")
# Read results from previous version
first_pass = fs::dir_ls("out/stacking_weight_chunks", glob = "*rds") %>% 
  map_dfr(read_rds) 



filtered_pass = run_settings %>% 
  mutate(file_name = glue::glue("out/loo/{name}.rds") %>% as.character()) %>% 
  left_join(first_pass)
  
# Filter them down
  # filter(weight >= threshold)

loo_lpd_matrix = read_loo_matrix(filtered_pass$file_name, cores = N_CORES)

loo_weights = stacking_weights_fast(loo_lpd_matrix, optim_method = "BFGS", optim_control = list()) 

write_rds(loo_weights, "out/loo_weights_raw.rds")


loo_weights = read_rds("out/loo_weights_raw.rds")

loo_weight_df =  tibble(weight = as.numeric(loo_weights),
                        name = basename(filtered_pass$file_name) %>% stringr::str_remove(".rds") ) %>% 
  left_join(run_settings, by = "name") %>%  
  arrange(desc(weight))
write_rds(loo_weight_df, "out/loo_weights.rds")
loo_weight_df = read_rds("out/loo_weights.rds")
# Problem w/ that is that the things aren't in the correct order


loo_text_tbl = loo_weight_df %>% 
  select(weight, name, model_txt) %>%
  # double unnest
  unnest(cols = c(model_txt)) %>% 
  unnest(cols = c(model_txt)) %>%
  separate(model_txt, into = c("y", "x"), sep = " ~ ") %>%
  spread(key = y, value = x) %>%
  arrange(desc(weight)) 
write_csv(loo_text_tbl, "out/loo_stacking.csv")

# Select draw ids
set.seed(52399832); loo_draw_ids = loo_weight_df %>% 
  mutate(n_draws = round(weight/weight[1] * 4000),
         draw_list = map(n_draws, ~sample(1:4000, size = .x)))
write_rds(loo_draw_ids, "out/loo_draw_ids.rds")
loo_draw_ids = read_rds("out/loo_draw_ids.rds")

read_and_draw = function(file_name, draw_list, ...) {
  fit = read_rds(file_name)
  brms::posterior_samples(fit, subset = draw_list) %>% 
    posterior::as_draws()
}
    # Probably need to change sim$[samples, chains, iter, warmup, n_save, warmup2, and permutation]
    

combined_data = loo_draw_ids %>% 
  mutate(file_name = glue("out/model_runs/{name}_model.rds")) %>% 
  select(file_name, draw_list) %>% 
  with({
    parallel::mcmapply(read_and_draw, file_name, draw_list, mc.cores = N_CORES)
  })
write_rds(combined_data, "out/combined_draws_raw.rds")


# Fill in zeroes for missing parameters before combining
add_null_parms = function(fit, model, ref_names = variables(combined_data[[1]])) {
  missing_vars = ref_names[! ref_names %in% variables(fit)]
  fit[missing_vars] = NA_real_
  fit %>% mutate(model = model)
}

padded_data = combined_data %>% imap(add_null_parms)
bound_draws = bind_rows(padded_data)
write_rds(bound_draws, "out/combined_draws_df.rds")

