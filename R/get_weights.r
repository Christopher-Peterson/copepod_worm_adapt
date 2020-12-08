# Once all runs are complete, calculate model weights

N_CORES = 48
library(readr)
library(loo)
library(brms)
library(purrr)
library(parallel)
library(glue)
library(dplyr)
run_settings = read_rds("model_settings.rds")


loo_file_names = glue("out/loo/{run_settings$name}.rds")


loo_files = parallel::mclapply(loo_file_names, read_rds, mc.cores = N_CORES)


# First run pseudo-Bayes model averaging to get a rough approximation
# of the weights; this is used as a first-step filter, to 
# filter the model set to the potentially useful ones 
# before running the more expensive stacking weights
options(mc.cores = N_CORES)
rough_weights = loo_model_weights(loo_files, cores = N_CORES, method = "pseudobma")

write_rds(rough_weights, "out/pbma_weights.rds")
rough_weights = read_rds("out/pbma_weights.rds")
rw2 = rough_weights %>% as.matrix %>% as.numeric
# First run pbma to filter down the options
run_settings %>% mutate(pbma = rw2) -> rough_set
rough_set %>% arrange(desc(pbma)) %>% select(id, model_txt, pbma )-> rs2
rs2 %>% filter(pbma > .001) %>% write_rds("out/pbma_filter.rds")

pbma_options = read_rds("out/pbma_filter.rds") %>% select(-model_txt) %>% 
  left_join(run_settings, by = "id")

nms = pbma_options %>% 
  filter(pbma > .005) %>% pull(name) 
loo_file_names2 = glue("out/loo/{nms}.rds")

loo_files2 = parallel::mclapply(loo_file_names2, read_rds, mc.cores = N_CORES)

options(mc.cores = N_CORES)
loo_weights = loo_model_weights(loo_files2, cores = N_CORES)
loo_weights %>% as.numeric %>% tibble(weight = ., name = nms) -> sv

left_join(sv, run_settings, by = "name") -> dn
dn %>% arrange(desc(weights))
write_rds(dn, "out/loo_weights.rds")

loo_weights = read_rds("out/loo_weights.rds")

loo_weights %>% select(weight, name, model_txt) %>% 
  unnest %>% unnest %>% 
  separate(model_txt, into = c("y", "x"), sep = " ~ ") %>% 
  spread(key = y, value = x) %>% 
  arrange(desc(weight)) -> loo_tbl
write_csv(loo_tbl, "out/loo_stacking.csv")

# If this doesn't save names, then just map the out/loo dir for them.