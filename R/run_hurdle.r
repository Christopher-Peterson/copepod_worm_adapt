# This takes a combination of model predictors, fits the model, and calculates LOO-IC
# It should be run from the command line once for each model combo.
suppressPackageStartupMessages({
  library(dplyr);library(tidyr); library(purrr); library(tibble); 
  library(readr); library(brms); library(loo)
})
N_CORES = 4
args = commandArgs(TRUE)
model_num = as.integer(args)
run_settings = read_rds("model_settings.rds")[model_num,]

copepods = read_csv("data/chapter2.copepods.cleaned.csv") %>% 
  mutate(native = cop.lake == worm.lake)

## Function to run models ====================
run_hurdle_model = function(formula, priors, name,
                            data = copepods, cores = N_CORES,
                            save_dir = "out/model_runs", ...) {
  # This function takes a formula and runs the hurdle model on it
  model_fit = brm(formula,
                  family = hurdle_poisson(), data = data, cores = cores,
                  prior = priors, sample_prior = "yes", ...)
  # It saves the output here
  write_rds(model_fit, file.path(save_dir, paste0(name, "_model.rds")))
  invisible(model_fit)
}

set.seed(run_settings$rng_seed) # make random numbers reproduceable
ad = .99 #run_settings$adapt_delta
fit = run_hurdle_model(run_settings$formula[[1]], run_settings$priors[[1]], 
                       run_settings$name, 
                       control = list(adapt_delta = ad))

# out = capture.output(try(run_hurdle_model(run_settings$formula[[1]], run_settings$priors[[1]], 
#                                       run_settings$name, 
#                                       control = list(adapt_delta = ad)), "message"))

# Determine how many divergences the model has
divergences =  nuts_params(fit, pars = "divergent__")$Value
# if(sum(divergences) > 0) {
#   # re run model if divergent
#   new_ad = ifelse(ad == .95, .99, .9999)
#   fit = run_hurdle_model(run_settings$formula[[1]], run_settings$priors[[1]], 
#                          run_settings$name, 
#                          control = list(adapt_delta = new_ad))
# # are there still divergences?
#   new_divergences =  nuts_params(fit, pars = "divergent__")$Value
#   if(sum(new_divergences)>0) stop("DIVERGENT TRANSITIONS in model ", model_num, " persist after re-run")
# }

# Get LOO info
looic = brms::loo(fit, cores = N_CORES, reloo = TRUE)
write_rds(looic, file.path("out/loo", paste0(run_settings$name, ".rds")))



# I apologize for my general lack of contact; it's asdf 
# 