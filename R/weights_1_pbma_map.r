# Calculating model weights is 2 part job
# Part 1 is to use PBMA model weights, which are fast & parallelizable but not the best
# From these, the top portion of models will be selected to do full bayesian stacking

N_CORES = 48
N_REPS = 6
library(readr)
library(loo)
library(brms)
library(purrr)
library(parallel)
library(glue)
library(dplyr)
run_settings = read_rds("model_settings.rds")
argv = commandArgs(TRUE)

SEED = argv[1] %>% as.numeric()


dir.create("out/pbma_weights")
out_file = file.path("out/pbma_weights", paste0("pbma_part_seed_", "SEED", ".rds"))

loo_file_names = glue("out/loo/{run_settings$name}.rds")


loo_lpd_matrix = parallel::mclapply(loo_file_names, 
                               function(.x) {
                                 read_rds(.x)$pointwise[, "elpd_loo"]
                               }, mc.cores = N_CORES)


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

rough_weights_list = parallel::mclapply(seq_len(N_CORES), 
                                           pseudo_bma_weights, mc.cores = N_CORES,
                                           lpd_mat = loo_lpd_matrix, BB_n = N_REPS)
pbma_weights = rowMeans(do.call(cbind, unclass(rough_weights_list)))
write_rds(pbma_weights, out_file)
