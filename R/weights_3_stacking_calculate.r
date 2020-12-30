# Runs full bayesian stacking on the subset of models with decent PBMA+ weights

N_CORES = 48
library(readr)
library(purrr)
library(dplyr)
library(glue)
library(loo)

loo_file_names = read_lines("out/loo_files.txt")

loo_files = parallel::mclapply(loo_file_names, read_rds, mc.cores = N_CORES)

options(mc.cores = N_CORES)
loo_weights = loo_model_weights(loo_files2)

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