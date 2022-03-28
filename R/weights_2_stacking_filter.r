# This script runs stacking on a subset of the overall data
# The next one will load all the subsets, filter out the irrelevant runs, then run it again on the whole set

N_CORES = 48
suppressPackageStartupMessages({
  library(readr)
  library(purrr)
  library(dplyr)
  library(glue)
  library(tidyr)
  library(rlang)
})

if(!exists("argv")) argv <- commandArgs(TRUE)
chunk = argv[1] %>% as.integer() %|% 1L
seed = argv[2] %>% as.integer() %|% 56789896L

output_file =glue("out/stacking_weight_chunks/{chunk}.rds")
set.seed(seed)
input_file = glue("out/pbma_weights/stacking_files_{chunk}.rds")
loo_file_names = read_rds(input_file)
source("R/weights_functions.r")

options(mc.cores = N_CORES)
loo_lpd_matrix = read_loo_matrix(loo_file_names, cores = N_CORES)

loo_weights = stacking_weights_fast(loo_lpd_matrix, optim_method = "BFGS", optim_control = list()) 


loo_weight_df =  tibble(weight = as.numeric(loo_weights), file_name = loo_file_names) 
write_rds(loo_weight_df, output_file)


# loo_weights = loo_model_weights(loo_files, cores = N_CORES)



run_settings = read_rds("model_settings_filtered.rds")
loo_text_tbl = loo_weight_df %>% 
  select(weight, name, model_txt) %>%
  unnest() %>% unnest() %>%
  separate(model_txt, into = c("y", "x"), sep = " ~ ") %>%
  spread(key = y, value = x) %>%
  arrange(desc(weight)) 
write_csv(loo_text_tbl, "out/loo_stacking.csv")

# If this doesn't save names, then just map the out/loo dir for them.