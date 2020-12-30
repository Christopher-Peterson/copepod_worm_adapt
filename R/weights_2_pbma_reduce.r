# Pass2 combines the output of pass 1 into a single PBMA+ model weight set, then filters
# all of the candidates down to a subset that have non-negligible weights

library(readr)
library(purrr)
library(dplyr)

PBMA_CUTOFF = 0.005
# Need to see if there are names still
rough_weights = dir("out/pbma_weights", full.names = TRUE) %>% 
  lapply(read_rds) %>% do.call(cbind, .) %>% rowMeans()
# Not sure if there are names here
write_rds(rough_weights, "out/pbma_weights.rds")

# rw2 = rough_weights %>% as.matrix %>% as.numeric
# First run pbma to filter down the options
stacking_names = read_rds("model_settings.rds") %>% mutate(pbma = rough_weights)  %>%
  arrange(desc(pbma)) %>% 
  # select(id, model_txt, pbma, name) %>% 
  filter(pbma > PBMA_CUTOFF) %>% 
  pull(name)
stacking_files = glue("out/loo/{stacking_names}.rds")
write_lines(stacking_files, "out/loo_files.txt")

# pbma_options = read_rds("out/pbma_filter.rds") %>% select(-model_txt) %>% 
#   left_join(run_settings, by = "id")
# 
# nms = pbma_options %>% 
#   filter(pbma > .005) %>% pull(name) 
# loo_file_names2 = glue("out/loo/{nms}.rds")
# 
# loo_files2 = parallel::mclapply(loo_file_names2, read_rds, mc.cores = N_CORES)
# 
# options(mc.cores = N_CORES)
# loo_weights = loo_model_weights(loo_files2, cores = N_CORES)
# loo_weights %>% as.numeric %>% tibble(weight = ., name = nms) -> sv
# 
# left_join(sv, run_settings, by = "name") -> dn
# dn %>% arrange(desc(weights))
# write_rds(dn, "out/loo_weights.rds")
# 
# loo_weights = read_rds("out/loo_weights.rds")
# 
# loo_weights %>% select(weight, name, model_txt) %>% 
#   unnest %>% unnest %>% 
#   separate(model_txt, into = c("y", "x"), sep = " ~ ") %>% 
#   spread(key = y, value = x) %>% 
#   arrange(desc(weight)) -> loo_tbl
# write_csv(loo_tbl, "out/loo_stacking.csv")
# 
# # If this doesn't save names, then just map the out/loo dir for them.