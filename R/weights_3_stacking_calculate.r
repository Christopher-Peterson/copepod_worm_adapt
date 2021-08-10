# Runs full bayesian stacking on the subset of models with decent PBMA+ weights

N_CORES = 48
suppressPackageStartupMessages({
  library(readr)
  library(purrr)
  library(dplyr)
  library(glue)
  library(rlang)
})

if(!exists("argv")) argv <- commandArgs(TRUE)
seed = argv[1] %>% as.integer() %|% 572806L
threshold = argv[2] %>% as.numeric() %|% 1e-5

options(mc.cores = N_CORES)
source("R/weights_functions.r")

# Read results from previous version
first_pass = fs::dir_ls("out/stacking_weight_chunks", glob = "*rds") %>% 
  map_dfr(read_rds) %>% 
# Filter them down
  filter(weight >= threshold)

loo_lpd_matrix = read_loo_matrix(first_pass$file_name, cores = N_CORES)

loo_weights = stacking_weights_fast(loo_lpd_matrix, optim_method = "BFGS", optim_control = list()) 


loo_weight_df =  tibble(weight = as.numeric(loo_weights), name = nms) %>% 
  left_join(run_settings, by = "name") %>%  
  arrange(desc(weights))
write_rds(loo_weight_df, "out/loo_weights.rds")


loo_text_tbl = loo_weight_df %>% 
  select(weight, name, model_txt) %>%
  unnest() %>% unnest() %>%
  separate(model_txt, into = c("y", "x"), sep = " ~ ") %>%
  spread(key = y, value = x) %>%
  arrange(desc(weight)) 
write_csv(loo_text_tbl, "out/loo_stacking.csv")

# If this doesn't save names, then just map the out/loo dir for them.