# Calculating model weights is 2 part job
# Part 1 is to use PBMA model weights, which are fast & parallelizable but not the best
# From these, the top portion of models will be selected to do full bayesian stacking

suppressPackageStartupMessages({
library(readr)
library(loo)
library(brms)
library(purrr)
library(parallel)
library(glue)
library(rlang)
library(dplyr)
})
argv = commandArgs(TRUE)
# Run arguments
SEED = argv[1] %>% as.numeric() %|% 1234
stacking_job = argv[2] %|% "jobs/stacking_p1.job"
N_CORES = argv[3] %>% as.numeric() %|%  48
N_REPS = argv[4] %>% as.numeric() %|%  60
STACKING_CHUNKS = argv[5] %>% as.numeric() %|% 4L
PBMA_CUTOFF_LOG = argv[6] %>% as.numeric() %|% -100 # CHANGE THIS


run_settings = read_rds("model_settings_filtered.rds")



dir.create("out/pbma_weights")
out_file = file.path("out/pbma_weights", paste0("pbma_part_seed_", SEED, ".rds"))

loo_file_names = glue("out/loo/{run_settings$name}.rds")


loo_lpd_list = parallel::mclapply(loo_file_names, 
                               function(.x) {
                                 read_rds(.x)$pointwise[, "elpd_loo"]
                               }, mc.cores = N_CORES)
loo_lpd_matrix = do.call(cbind, loo_lpd_list)
# First run pseudo-Bayes model averaging to get a rough approximation
# of the weights; this is used as a first-step filter, to 
# filter the model set to the potentially useful ones 
# before running the more expensive stacking weights
options(mc.cores = N_CORES)


# By default the pseudobma_weights runs a for loop over BB_n,
# then eventually col sums the results together. This can be map-reduced

pbma_func = function(iter, lpd_mat =loo_lpd_matrix,  ...) {
  pseudobma_weights(lpd_mat, BB = TRUE, ...)
}
set.seed(SEED)


rough_weights_list = mclapply(seq_len(N_CORES), 
                              pbma_func, mc.cores = N_CORES,
                              lpd_mat = loo_lpd_matrix, BB_n = N_REPS)
pbma_weights = rowMeans(do.call(cbind, unclass(rough_weights_list)))
write_rds(pbma_weights, out_file)

#### Determine which models to move onto stacking ####


str(loo_lpd_matrix)


stacking_names = run_settings %>% 
  mutate(pbma = pbma_weights)  %>%
  arrange(desc(pbma)) %>% 
  # select(id, model_txt, pbma, name) %>% 
  filter(log10(pbma) > PBMA_CUTOFF_LOG) %>% 
  mutate(chunk = (((1:n()) %% STACKING_CHUNKS) + 1),
         loo_files = glue("out/loo/{name}.rds") %>% as.character()) %>% 
  select(chunk, loo_files) %>% 
  tidyr::chop(loo_files) %>% 
  mutate(file_name = glue("out/pbma_weights/stacking_files_{chunk}.rds"),
         seed = rdunif(n(), 65536))
# write output
with(stacking_names, walk2(loo_files, file_name, write_rds))

glue_data(stacking_names, "R-stan R/weights_2_stacking_filter.r {chunk} {seed}") %>% write_lines(stacking_job)

dir.create("out/stacking_weight_chunks")





